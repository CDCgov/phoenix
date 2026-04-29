#!/usr/bin/env python

import os
import sys
import errno
import argparse
from pathlib import Path
import pandas as pd

def get_version():
    return "1.0.0"

def parse_args(args=None):
    parser = argparse.ArgumentParser(
    description="Validate a file-of-filepaths (one filepath per line).",  epilog="Example usage: python check_samplesheet.py <FILE_IN> <FILE_OUT>", )
    parser.add_argument("FILE_IN", help="Input file (one filepath per line).")
    parser.add_argument("FILE_OUT", help="Output file.")
    parser.add_argument("--version", action="version", version=get_version())
    return parser.parse_args(args)

def make_dir(path):
    if path:
        try:
            os.makedirs(path)
        except OSError as e:
            if e.errno != errno.EEXIST:
                raise

def read_paths(file_in: str) -> pd.Series:
    # Read one path per line; keep as strings; drop empty/whitespace-only lines
    df = pd.read_csv(file_in, header=None, dtype=str, keep_default_na=False)
    paths = df[0].astype(str).str.strip()
    paths = paths[paths != ""]
    return paths

def check_for_duplicates(paths: pd.Series):
    dup_mask = paths.duplicated(keep=False)
    if dup_mask.any():
        dups = paths[dup_mask].tolist()
        raise ValueError(
            "Duplicate filepath(s) found (file locations must be unique):\n"
            + "\n".join(sorted(set(dups)))
        )

def check_files_exist(paths: pd.Series):
    missing = [p for p in paths if not Path(p).exists()]
    if missing:
        raise ValueError(
            "The following file(s) do not exist but are required:\n"
            + "\n".join(missing)
        )

def check_file_types(paths: pd.Series):
    suffixes = paths.apply(lambda p: Path(p).suffix.lower())
    allowed = {".tsv", ".xlsx"}

    invalid = sorted(set(suffixes) - allowed)
    if invalid:
        raise ValueError(f"Invalid file type(s) found: {invalid}. Allowed: {sorted(allowed)}")

    unique = sorted(set(suffixes))
    if len(unique) != 1:
        raise ValueError(f"Mixed file types detected: {unique}. Files must all be the same type.")

    return unique[0]  # '.tsv' or '.xlsx'

def check_samplesheet(file_in, file_out):
    paths = read_paths(file_in)

    if paths.empty:
        raise ValueError("Input file contains no file paths.")

    check_for_duplicates(paths)
    check_files_exist(paths)
    file_type = check_file_types(paths)

    out_dir = os.path.dirname(file_out)
    make_dir(out_dir)

    # Write cleaned/normalized paths back out (one per line)
    with open(file_out, "w", encoding="utf-8") as fout:
        for p in paths:
            fout.write(p + "\n")

    return file_type

def main(args=None):
    args = parse_args(args)
    file_type = check_samplesheet(args.FILE_IN, args.FILE_OUT)
    print(f"OK: validated {args.FILE_IN} (all files are {file_type})")

if __name__ == "__main__":
    sys.exit(main())
