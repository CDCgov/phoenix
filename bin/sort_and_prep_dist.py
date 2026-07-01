#!/usr/bin/env python3

"""
Description: Grabs the best species match based on %/read hits from the kraken tool run.
             Simplified for nextflow inclusion.

Usage: ./sort_and_prep_dists.py -x dist_file [-a assembly] [-o outdir] [-t terra] [-V]

Output location: same directory as input dist_file

Modules required: None

Created by Nick Vlachos (nvx4@cdc.gov)
Converted to Python: (05/06/2025)
"""

import argparse
import re
import shutil
import ssl
import sys
import urllib.request
from pathlib import Path

__version__ = "2.1"
MAX_HITS = 40       # Target is 20, allow twice as many for poor-quality isolates
MIN_KMERS = 5       # Minimum kmer threshold to filter low-confidence hits
CUTOFF_RANK = 20    # Row index (1-based) used to determine the distance cutoff
NCBI_FTP_BASE = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF"

# Regex patterns
GCF_FULL_RE  = re.compile(r"(GCF_\d{9}\.\d+[^/]*)")
GCF_ACCN_RE  = re.compile(r"(GCF_\d{9}\.\d+)")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser( description="Sort MASH distances and download best-matching reference genomes.")
    parser.add_argument("-x", "--dist-file",    required=True, help="Path to MASH distance .txt file")
    parser.add_argument("-a", "--assembly-file", default="",   help="Assembly file path")
    parser.add_argument("-o", "--outdir",        default=".",  help="Output directory for downloaded genomes")
    parser.add_argument("-V", "--version", action="version", version=f"%(prog)s: {__version__}")
    return parser.parse_args()


def build_ssl_context() -> ssl.SSLContext:
    """
    Return an SSL context. Certificate verification is intentionally
    disabled to match the original --no-check-certificate behaviour.
    """
    ctx = ssl.create_default_context()
    ctx.check_hostname = False
    ctx.verify_mode    = ssl.CERT_NONE
    return ctx


def parse_dist_line(line: str) -> tuple[str, float, int] | None:
    """
    Parse one line from the MASH dist file.

    Expected whitespace-delimited columns:
        0: source/reference name
        1: query name
        2: distance  (may be scientific notation, e.g. 1.5e-04)
        3: p-value
        4: kmer_matches/total  (e.g. 450/1000)

    Returns (source, dist, kmers) or None if the line is malformed.
    """
    parts = line.split()
    if len(parts) < 5:
        return None
    try:
        source = parts[0]
        dist   = float(parts[2])                    # handles scientific notation natively
        kmers  = int(parts[4].split("/")[0])
        return source, dist, kmers
    except (ValueError, IndexError):
        return None


def sort_dist_file(dist_path: Path) -> tuple[Path, list[tuple[str, float, int]]]:
    """
    Read, parse, and sort the dist file by distance (ascending).
    Writes a companion _sorted.txt file.

    Returns (sorted_path, sorted_rows).
    """
    rows = []
    with dist_path.open() as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            parsed = parse_dist_line(line)
            if parsed:
                rows.append(parsed)

    rows.sort(key=lambda r: r[1])   # sort by distance, ascending

    sorted_path = dist_path.with_name(dist_path.stem + "_sorted.txt")
    with sorted_path.open("w") as fh:
        for source, dist, kmers in rows:
            fh.write(f"{source}\t\t{dist}\t\t{kmers}\n")

    return sorted_path, rows


def build_ncbi_url(source: str) -> tuple[str, str] | None:
    """
    Construct the NCBI FTP URL for a GCF accession found in *source*.

    Returns (url, filename) or None if no GCF accession is found.
    """
    m_full  = GCF_FULL_RE.search(source)
    m_accn  = GCF_ACCN_RE.search(source)
    if not m_accn:
        return None

    accession = m_accn.group(1)                     # e.g. GCF_000005845.2
    filename  = m_full.group(1) if m_full else accession
    filename  = filename.removesuffix(".fna")

    # Strip GCF_ prefix and version to get the 9-digit number
    num_part = accession.removeprefix("GCF_").split(".")[0]   # e.g. "000005845"
    alpha    = num_part[0:3]
    bravo    = num_part[3:6]
    charlie  = num_part[6:9]

    url = f"{NCBI_FTP_BASE}/{alpha}/{bravo}/{charlie}/{filename}/{filename}_genomic.fna.gz"
    return url, filename


def download_genome(url: str, dest: Path, ssl_ctx: ssl.SSLContext) -> bool:
    """
    Download *url* to *dest*. Returns True on success (file exists and is non-empty).
    """
    print(f"Trying - {url} -> {dest}")
    try:
        with urllib.request.urlopen(url, context=ssl_ctx) as resp, dest.open("wb") as out:
            shutil.copyfileobj(resp, out, length=1 << 20)  # 1 MB chunks
    except Exception as exc:
        print(f"  Download error: {exc}")
        dest.unlink(missing_ok=True)
        return False

    if dest.stat().st_size == 0:
        dest.unlink()
        return False

    return True


def main() -> None:
    args = parse_args()
    dist_path = Path(args.dist_file)
    if not dist_path.is_file():
        print(f"dist file not found: {dist_path}", file=sys.stderr)
        sys.exit(1)

    outdir     = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    ssl_ctx    = build_ssl_context()
    sample_name = dist_path.stem                    # basename without extension

    hits_file = Path(f"{sample_name}_best_MASH_hits.txt")
    log_file  = Path(f"{sample_name}_genome_download_log.txt")
    log_file.write_text("")                         # reset log

    # Sort
    _sorted_path, rows = sort_dist_file(dist_path)

    if len(rows) < CUTOFF_RANK:
        print("Warning: fewer rows than cutoff rank; using last available row.")
    cutoff = rows[min(CUTOFF_RANK, len(rows)) - 1][1]
    print(f"Cutoff IS: {cutoff}")

    matches = 0
    with hits_file.open("w") as hits_fh, log_file.open("a") as log_fh:
        for source, dist, kmers in rows:
            print(f"{source}")
            print(f"dist-{dist} - {source}")

            # Apply filters
            if dist > cutoff or kmers <= MIN_KMERS or matches >= MAX_HITS:
                break

            dest = outdir / f"{source}.gz"

            # Already cached locally
            if dest.is_file():
                hits_fh.write(f"{dest}\n")
                matches += 1
                continue

            # Attempt NCBI download
            result = build_ncbi_url(source)
            if result is None:
                msg = f"⚠ No valid GCF accession found in: {source}"
                print(msg)
                log_fh.write(msg + "\n")
                continue

            url, filename = result
            print(f"Copying - {filename}")

            if download_genome(url, dest, ssl_ctx):
                hits_fh.write(f"{dest}\n")
                matches += 1
                log_fh.write(f"✓ Successfully downloaded: {source}\n")
            else:
                log_fh.write(f"✗ {source} did not download correctly\n")

    if hits_file.stat().st_size == 0:
        hits_file.unlink()

    print(f"Done. {matches} genome(s) written to {hits_file}")


if __name__ == "__main__":
    main()