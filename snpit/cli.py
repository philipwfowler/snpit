"""Whole genome SNP based identification of members of the Mycobacterium tuberculosis
complex.
SNP-IT allows rapid Mycobacterial speciation of VCF files aligned to NC000962 (H37Rv).
"""
import argparse
import sys
from pathlib import Path
from .version import __version__

DEFAULT_THRESHOLD = 10.0


def cli():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="""Path to the VCF or FAST(A/Q) file to read and classify. File can be 
        multi-sample and/or compressed.""",
        type=Path,
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        help="""Path to output results to. Default is STDOUT (-).""",
        default="-",
    )
    parser.add_argument(
        "--threshold",
        type=float,
        help="""The percentage of snps above which a sample is considered to belong to 
        a lineage. [{}]""".format(
            DEFAULT_THRESHOLD
        ),
        default=DEFAULT_THRESHOLD,
    )
    parser.add_argument(
        "--filter", help="Whether to adhere to the FILTER column.", action="store_true"
    )
    parser.add_argument(
        "--status",
        help="""Whether to adhere to the STATUS column. This is a custom 
        field that gives more fine-grained control over whether a sample passes a 
        user-defined filtering criterion, even if the record has PASS in FILTER.""",
        action="store_true",
    )
    parser.add_argument(
        "-v",
        "--version",
        help="""Show the program's version number and exit.""",
        action="version",
        version=__version__,
    )
    args = parser.parse_args()

    if args.output == "-":
        args.output = sys.stdout
    else:
        p = Path(args.output)
        if not p.parent.is_dir():
            raise NotADirectoryError(
                "Directory specified for output file does not exist: {}".format(
                    p.parent
                )
            )
        args.output = p.open("w")

    if not args.input.is_file():
        raise FileNotFoundError("Input file {} does not exist!".format(str(args.input)))

    return args
