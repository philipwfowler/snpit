#!/usr/bin/env python

"""Whole genome SNP based identification of members of the Mycobacterium tuberculosis
complex.
SNP-IT allows rapid Mycobacterial speciation of VCF files aligned to NC000962 (H37Rv).
"""

from pathlib import Path
import sys
from snpit import SnpIt
import argparse

THRESHOLD = 10.0


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
            THRESHOLD
        ),
        default=THRESHOLD,
    )
    parser.add_argument(
        "--filter",
        help="Whether to adhere to the FILTER column.",
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
    tb = SnpIt(
        threshold=options.threshold,
        ignore_filter=not options.filter,
    )

    compressed = True if options.input.endswith("gz") else False
    suffixes = Path(options.input).suffixes

    if ".vcf" in suffixes:
        results = tb.classify_vcf(options.input)
    elif any(ext in suffixes for ext in [".fa", ".fasta"]):
        results = tb.load_fasta(options.input, compression=compressed)
    else:
        raise Exception(
            "Only VCF and FASTA files are allowed as inputs (may be compressed with gzip,bzip2)"
        )

    print("Species\tLineage\tSublineage\tPercentage\tInput", file=options.output)
    output = ""
    for sample_name, classification in results.items():
        output += f"{classification}\n"
    # if tb.percentage is not None:
    #     lineage = tb.lineage or "N/A"
    #     sublineage = tb.sublineage or "N/A"
    #     percentage = round(tb.percentage, 2)
    #     output = "{species}\t{lineage}\t{sublineage}\t{percentage}\t{input_file}".format(
    #         species=tb.species,
    #         lineage=lineage,
    #         sublineage=sublineage,
    #         percentage=percentage,
    #         input_file=options.input,
    #     )
    # else:
    #     output = "{0}\t{0}\t{0}\t{0}\t{1}".format("N/A", options.input)
    print(output, file=options.output)
    options.output.close()
