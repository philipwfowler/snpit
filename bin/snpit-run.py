"""Whole genome SNP based identification of members of the Mycobacterium tuberculosis
complex.
SNP-IT allows rapid Mycobacterial speciation of VCF files aligned to NC000962 (H37Rv).
"""

import argparse
import sys
from pathlib import Path

from snpit.core import SnpIt, output_results

DEFAULT_THRESHOLD = 10.0


def cli():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input",
        required=True,
        help="""Path to the VCF or FASTA file to read and classify 
            (can be bzip2ed/gzipped)""",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        help="""Path to output results to.
                            Default is STDOUT (-).""",
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

    if not Path(args.input).is_file():
        raise FileNotFoundError("Input file {} does not exist!".format(args.input))

    return args


if __name__ == "__main__":
    options = cli()

    # create an instance (this loads all the lineages)
    snpit = SnpIt(
        threshold=options.threshold,
        ignore_filter=not options.filter,
        ignore_status=not options.status,
    )

    compressed = True if options.input.endswith("gz") else False
    suffixes = Path(options.input).suffixes

    if ".vcf" in suffixes:
        results = snpit.classify_vcf(Path(options.input))
    elif any(ext in suffixes for ext in [".fa", ".fasta"]):
        results = snpit.classify_fasta_sequence(options.input, compression=compressed)
    else:
        raise Exception(
            "Only VCF and FASTA files are allowed as inputs (may be compressed with gzip,bzip2)"
        )

    output_results(options.output, results)
    options.output.close()
