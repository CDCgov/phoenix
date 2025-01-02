#!/usr/bin/env python


"""Provide a command line tool to validate and transform tabular samplesheets."""


import argparse
import csv
import glob
import logging
import sys
from collections import Counter
from pathlib import Path

# Function to get the script version
def get_version():
    return "1.0.0"

def parse_args(argv=None):
    """Define and immediately parse command line arguments."""
    parser = argparse.ArgumentParser( description="Validate and transform a tabular samplesheet.",
        epilog="Example: python check_samplesheet.py samplesheet.csv samplesheet.valid.csv", )
    parser.add_argument("file_in", metavar="FILE_IN", type=Path, help="Tabular input samplesheet in CSV or TSV format.", )
    parser.add_argument("file_out", metavar="FILE_OUT", type=Path, help="Transformed output samplesheet in CSV format.", )
    parser.add_argument("-l", "--log-level", help="The desired log level (default WARNING).", choices=("CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG"), default="WARNING",)
    parser.add_argument('--version', action='version', version=get_version())# Add an argument to display the version
    return parser.parse_args(argv)

logger = logging.getLogger()


class RowChecker:
    """
    Define a service that can validate and transform each given row.

    Attributes:
        modified (list): A list of dicts, where each dict corresponds to a previously
            validated and transformed row. The order of rows is maintained.

    """
    def __init__(
        self,
        sample_col="sample",
        first_col="directory",
        **kwargs,
    ):
        """
        Initialize the row checker with the expected column names.

        Args:
            sample_col (str): The name of the column that contains the sample name (default "sample").
            first_col (str): The name of the column that contains the only the directory of a sample to update (default "directory").
        """
        super().__init__(**kwargs)
        self._sample_col = sample_col
        self._first_col = first_col
        self._seen = set()
        self.modified = []

    def validate_and_transform(self, row):
        """
        Perform all validations on the given row and insert the read pairing status.

        Args:
            row (dict): A mapping from column headers (keys) to elements of that row
                (values).

        """
        self._validate_sample(row)
        self._validate_first(row)
        self._seen.add((row[self._sample_col], row[self._first_col]))
        self.modified.append(row)

    def _validate_sample(self, row):
        """Assert that the sample name exists and convert spaces to underscores."""
        if len(row[self._sample_col]) <= 0:
            raise AssertionError("Sample input is required.")
        # Sanitize samples slightly.
        row[self._sample_col] = row[self._sample_col].replace(" ", "_")

    def _validate_first(self, row):
        """Assert that the directory entry is non-empty and has the right format."""
        if len(row[self._first_col]) <= 0:
            raise AssertionError("At least a valid directory path is required.")

    def validate_unique_samples(self):
        """
        Assert that the combination of sample name and directory name is unique.
        """
        if len(self._seen) != len(self.modified):
            raise AssertionError("The pair of sample name and directory must be unique.")
        seen = Counter()
        samples_list = []
        for row in self.modified:
            sample = row[self._sample_col]
            if sample not in samples_list: # check that sample names are unique.
                samples_list.append(sample)
            else:
                raise AssertionError("ERROR: {} is used as a sample name more than once. Samples IDs should be unique.".format(sample))
            seen[sample] += 1
            #row[self._sample_col] = f"{sample}_T{seen[sample]}"
            row[self._sample_col] = f"{sample}"

    def validate_needed_files(self):
        """
        Checking for necesary files to run the pipeline
        """
        if len(self._seen) != len(self.modified):
            raise AssertionError("The pair of sample name and directory must be unique.")
        
        files = []
        # Define the file path
        sample_name = self._sample_col
        if str(self._first_col).endswith('/'):
            path = str(self._first_col)[:-1]
        else:
            path = str(self._first_col)
        files.append(path + "/file_integrity/" + sample_name + "_scaffolds_summary.txt")
        files.append(path + "/fastp_trimd/" + sample_name + "_1.trim.fastq.gz")
        files.append(path + "/fastp_trimd/" + sample_name + "_2.trim.fastq.gz")
        files.append(path + "/assembly/" + sample_name + ".filtered.scaffolds.fa.gz")
        files.append(path + "/annotation/" + sample_name + ".faa")
        files.append(path + "/annotation/" + sample_name + ".gff")
        files.append(path + "/" + sample_name + ".tax")
        files.append(path + "/" + sample_name + "_summaryline.tsv")
        files.append(path + "/" + sample_name + ".synopsis")
        files.append(glob.glob(path + "../" + "*_GRiPHin_Summary.xlsx")[0])
        files.append(glob.glob(path + "../" + "*_GRiPHin_Summary.tsv")[0])
        files.append(glob.glob(path + "../" + "*_Phoenix_Summary.tsv")[0])
        for file_path in files:
            # Check if the file exists
            if Path(file_path).exists():
                pass
            else:
                logger.critical("The file {} does not exist.".format(file_path))
                sys.exit(1)


def read_head(handle, num_lines=10):
    """Read the specified number of lines from the current position in the file."""
    lines = []
    for idx, line in enumerate(handle):
        if idx == num_lines:
            break
        lines.append(line)
    return "".join(lines)


def sniff_format(handle, file_in):
    """
    Detect the tabular format.

    Args:
        handle (text file): A handle to a `text file`_ object. The read position is
        expected to be at the beginning (index 0).

    Returns:
        csv.Dialect: The detected tabular format.

    .. _text file:
        https://docs.python.org/3/glossary.html#term-text-file

    """
    peek = read_head(handle)
    handle.seek(0)
    sniffer = csv.Sniffer()
    dialect = sniffer.sniff(peek)
    return dialect


def check_samplesheet(file_in, file_out):
    """
    Check that the tabular samplesheet has the structure expected by nf-core pipelines.

    Validate the general shape of the table, expected columns, and each row. Also add
    an additional column which records whether one or two FASTQ reads were found.

    Args:
        file_in (pathlib.Path): The given tabular samplesheet. The format can be either
            CSV, TSV, or any other format automatically recognized by ``csv.Sniffer``.
        file_out (pathlib.Path): Where the validated and transformed samplesheet should
            be created; always in CSV format.

    Example:
        This function checks that the samplesheet follows the following structure,
        see also the `viral recon samplesheet`_::

            sample,assembly
            SAMPLE_1,Directory
            SAMPLE_2,Directory
    """
    required_columns = {"sample", "directory"}
    # See https://docs.python.org/3.9/library/csv.html#id3 to read up on `newline=""`.
    with file_in.open(newline="") as in_handle:
        reader = csv.DictReader(in_handle, dialect=sniff_format(in_handle, file_in))
        # Validate the existence of the expected header columns.
        if not required_columns.issubset(reader.fieldnames):
            req_cols = ", ".join(required_columns)
            logger.critical(f"The sample sheet **must** contain these column headers: {req_cols}.")
            sys.exit(1)
        # Validate each row.
        checker = RowChecker()
        for i, row in enumerate(reader):
            try:
                checker.validate_and_transform(row)
            except AssertionError as error:
                logger.critical(f"{str(error)} On line {i + 2}.")
                sys.exit(1)
        checker.validate_unique_samples()
    header = list(reader.fieldnames)
    #header.insert(1, "single_end")
    # See https://docs.python.org/3.9/library/csv.html#id3 to read up on `newline=""`.
    with file_out.open(mode="w", newline="") as out_handle:
        writer = csv.DictWriter(out_handle, header, delimiter=",")
        writer.writeheader()
        for row in checker.modified:
            writer.writerow(row)


def main(argv=None):
    """Coordinate argument parsing and program execution."""
    args = parse_args(argv)
    logging.basicConfig(level=args.log_level, format="[%(levelname)s] %(message)s")
    if not args.file_in.is_file():
        logger.error(f"The given input file {args.file_in} was not found!")
        sys.exit(2)
    args.file_out.parent.mkdir(parents=True, exist_ok=True)
    check_samplesheet(args.file_in, args.file_out)


if __name__ == "__main__":
    sys.exit(main())
