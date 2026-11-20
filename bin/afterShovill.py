#!/usr/bin/env python3

"""
Prepare Shovill assembly outputs for PHoeNIx.

Converts Shovill assembly files into the file names
expected by the existing PHoeNIx assembly workflow.
"""

import argparse
import gzip
import shutil
from pathlib import Path

__version__ = "1.0"


def parse_args():
    parser = argparse.ArgumentParser(
        description="Prepare Shovill outputs for PHoeNIx."
    )
    parser.add_argument(
        "--prefix",
        required=True,
    )
    parser.add_argument(
        "--outdir",
        default="shovill_output",
        type=Path,
    )
    parser.add_argument(
        "-V",
        "--version",
        action="version",
        version=f"%(prog)s: {__version__}",
    )
    return parser.parse_args()


def gzip_file(source, destination):
    """Compress a file reproducibly without a gzip timestamp."""
    with source.open("rb") as f_in:
        with destination.open("wb") as f_out:
            with gzip.GzipFile(
                filename="",
                mode="wb",
                fileobj=f_out,
                mtime=0,
            ) as gz:
                shutil.copyfileobj(
                    f_in,
                    gz,
                )


def main():
    args = parse_args()

    prefix = args.prefix
    outdir = args.outdir

    outcome_file = Path(
        f"{prefix}_spades_outcome.csv"
    )

    final_assembly = (
        outdir
        / "contigs.fa"
    )

    raw_assembly = (
        outdir
        / "spades.fasta"
    )

    # Preserve the status field written by the
    # Nextflow process before this script runs.
    status = outcome_file.read_text().strip()

    # Shovill's final polished assembly becomes the
    # PHoeNIx scaffold output.
    scaffolds_created = False

    if (
        final_assembly.is_file()
        and final_assembly.stat().st_size > 0
    ):
        gzip_file(
            final_assembly,
            Path(
                f"{prefix}.scaffolds.fa.gz"
            ),
        )

        scaffolds_created = True

    # Preserve the underlying SPAdes assembly as the
    # PHoeNIx contig output.
    contigs_created = False

    if (
        raw_assembly.is_file()
        and raw_assembly.stat().st_size > 0
    ):
        gzip_file(
            raw_assembly,
            Path(
                f"{prefix}.contigs.fa.gz"
            ),
        )

        contigs_created = True

    status += (
        ",scaffolds_created"
        if scaffolds_created
        else ",no_scaffolds"
    )

    status += (
        ",contigs_created"
        if contigs_created
        else ",no_contigs"
    )

    # Preserve the FAIry summary if the final
    # assembly was not produced.
    if not scaffolds_created:
        old = Path(
            f"{prefix}_summary_old_3.txt"
        )

        if old.is_file():
            old.replace(
                Path(
                    f"{prefix}_trimstats_summary.txt"
                )
            )

    outcome_file.write_text(
        status
    )

    # Shovill may name the assembly graph after
    # the underlying assembler.
    graph = (
        outdir
        / "spades.gfa"
    )

    if (
        graph.is_file()
        and graph.stat().st_size > 0
    ):
        gzip_file(
            graph,
            Path(
                f"{prefix}.assembly.gfa.gz"
            ),
        )


if __name__ == "__main__":
    main()
