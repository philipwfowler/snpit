from snpit.cli import cli
from snpit.snpit import SnpIt, output_results


def main(args=None):
    if args is None:
        args = cli()

    snpit = SnpIt(
        threshold=args.threshold,
        ignore_filter=not args.filter,
        ignore_status=not args.status,
    )

    suffixes = args.input.suffixes

    if ".vcf" in suffixes:
        results = snpit.classify_vcf(args.input)
    elif any(ext in suffixes for ext in [".fa", ".fasta"]):
        results = snpit.classify_fasta(args.input)
    else:
        raise Exception(
            "Only VCF or FASTA files are allowed as input (may be compressed)"
        )

    output_results(args.output, results)
    args.output.close()


if __name__ == "__main__":
    main()
