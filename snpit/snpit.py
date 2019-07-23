import csv
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Tuple, List, Dict, IO

import pysam

from .genotype import Genotype, UnexpectedGenotypeError
from .lineage import Lineage

LIBRARY_DIR = Path(__file__).parent.parent / "lib/"


class SnpIt(object):
    """
    The snpit class is designed to take a VCF file and return the most likely lineage
    based on Sam Lipworth's SNP-IT.

    The methods have been separated so it can be incorporated into single Python
    scripts that processes multiple VCF files.
    """

    def __init__(self, threshold: float, ignore_filter=True, ignore_status=True):

        """
        Args:
            threshold: The percentage of snps above which a sample is
            considered to belong to a lineage.
            ignore_filter: Whether to ignore the FILTER column in VCF input.
            ignore_status: Whether to ignore the STATUS field (if present).
        """
        self.threshold = threshold
        self.ignore_filter = ignore_filter
        self.ignore_status = ignore_status

        # library file which contains a list of all the lineages and sub-lineages
        library_csv = LIBRARY_DIR / "library.csv"

        # self.lineages = load_lineages_from_csv(library_csv)
        self.lineages, self.lineage_positions = load_lineages_from_csv(library_csv)

    def classify_vcf(self, vcf_path: Path) -> Dict[str, Tuple[float, Lineage]]:
        """Loads the vcf file and then, for each lineage, identify the base at each of
        the identifying positions in the genome.
        """
        sample_lineage_counts = self.count_lineage_classifications_for_samples_in_vcf(
            vcf_path
        )

        results = {}
        for sample_name, lineage_counts in sample_lineage_counts.items():
            results[sample_name] = self.determine_lineage(lineage_counts)
        return results

    def count_lineage_classifications_for_samples_in_vcf(
        self, vcf_path: Path
    ) -> defaultdict:
        """For each sample in the given VCF, count the number of records that match
        each lineage.
        """
        vcf = pysam.VariantFile(str(vcf_path))

        sample_lineage_counts = defaultdict(Counter)

        for record in vcf:
            if self.is_record_invalid(record):
                continue

            for sample_idx, (sample_name, sample_info) in enumerate(
                record.samples.items()
            ):
                if not self.ignore_status and sample_info["STATUS"] == "FAIL":
                    continue

                try:
                    genotype = Genotype(*sample_info["GT"])
                except TypeError as err:
                    genotype = minos_gt_in_wrong_position_fix(record, sample_idx)
                    if genotype is None:
                        raise err

                variant = self.get_variant_for_genotype_in_vcf_record(genotype, record)

                if not variant:
                    continue

                lineages_sharing_variant_with_sample = self.lineage_positions.get(
                    record.pos, {}
                ).get(variant, [])

                sample_lineage_counts[sample_name].update(
                    lineages_sharing_variant_with_sample
                )

        # need to do this so that sample appears in the results even if
        # there are no counts for any lineages
        for sample_name in vcf.header.samples:
            if sample_name not in sample_lineage_counts:
                sample_lineage_counts[sample_name].update([])

        return sample_lineage_counts

    def is_record_invalid(self, record: pysam.VariantRecord) -> bool:
        return record.pos not in self.lineage_positions or (
            not self.ignore_filter and "PASS" not in record.filter.keys()
        )

    @staticmethod
    def get_variant_for_genotype_in_vcf_record(
        genotype: Genotype, record: pysam.VariantRecord
    ) -> str:
        """Retrieves the variant a genotype maps to for a given record.

        Args:
            genotype: The genotype call for the sample.
            record: A VCF record object.
        Returns:
            str: A hyphen if the call is null (ie ./.) or the alt variant if
            the call is alt. Returns an empty string if the call is ref or heterozygous.
        """
        if genotype.is_reference() or genotype.is_heterozygous():
            variant = ""
        elif genotype.is_alt():
            alt_call = max(genotype.call())
            variant = record.alleles[alt_call]
        elif genotype.is_null():
            variant = "-"
        else:
            raise UnexpectedGenotypeError(
                """Got a genotype for which a Ref/Alt/Null call could not be 
                    determined: {}.\nPlease raise this with the developers.""".format(
                    genotype.call()
                )
            )
        return variant

    def classify_fasta(self, fasta_path: Path) -> Dict[str, Tuple[float, Lineage]]:
        fasta = pysam.FastaFile(str(fasta_path))
        sample_lineage_counts = dict()

        for sample_name in fasta.references:
            sequence = fasta.fetch(sample_name)
            sample_lineage_counts[sample_name] = self.classify_fasta_sequence(sequence)

        results = {}
        for sample_name, lineage_counts in sample_lineage_counts.items():
            results[sample_name] = self.determine_lineage(lineage_counts)
        return results

    def classify_fasta_sequence(self, sequence: str) -> Counter:
        """For each lineage, identify the base at each of the identifying positions
        """
        counts = Counter()
        seq_len = len(sequence)

        for pos in self.lineage_positions:
            if pos >= seq_len:
                continue

            base_in_sequence = sequence[pos - 1]  # all lineages positions are 1-based
            base_in_known_variants = base_in_sequence in self.lineage_positions[pos]
            if not base_in_known_variants:
                continue

            counts.update(self.lineage_positions[pos][base_in_sequence])

        return counts

    def determine_lineage(self, lineage_counts: Counter) -> Tuple[float, Lineage]:
        """
        Having read the VCF file, for each lineage, calculate the percentage of SNP
        present in the sample.
        Note that this means the percentages will not add up to 100%.

        Returns:
            tuple of (lineage,percentage)
        """
        lineage_shared_percentage = self.shared_variants_percentage_for_each_lineage_from_counts(
            lineage_counts
        )
        lineage_shared_percentage.sort(reverse=True)

        if not lineage_shared_percentage:
            identified_lineage = Lineage()
            identified_lineage_percentage = 0
        else:
            identified_lineage = lineage_shared_percentage[0][1]
            identified_lineage_percentage = lineage_shared_percentage[0][0]

        if identified_lineage_percentage > self.threshold:

            if (
                identified_lineage.lineage == "Lineage 4"
                and not identified_lineage.sublineage
                and len(lineage_shared_percentage) > 1
            ):
                if self.is_next_best_lineage4(*lineage_shared_percentage[1]):
                    identified_lineage = lineage_shared_percentage[1][1]

        return identified_lineage_percentage, identified_lineage

    def shared_variants_percentage_for_each_lineage_from_counts(
        self, counts: Counter
    ) -> List[Tuple[float, Lineage]]:
        """Calculates the percentage of shared variants across all lineages.

        Args:
            counts: A Counter containing the number of times a lineage had a shared
            variant with a sample.
        Returns:
            A sorted list of tuples where in the first element in a tuple is the number
            of shared SNPs divided by the total number of SNPs that define a lineage
            and the second element is the lineage for the percentage.
        """
        percentage_shared_snps_for_lineages = []

        for lineage_name, shared_calls in counts.items():
            lineage = self.lineages[lineage_name]
            ref_calls = len(lineage.snps)

            percentage_shared_snps_for_lineages.append(
                ((shared_calls / ref_calls) * 100, lineage)
            )

        return percentage_shared_snps_for_lineages

    def is_next_best_lineage4(
        self, next_best_lineage_percentage: float, next_best_lineage: Lineage
    ) -> bool:
        return (
            next_best_lineage.lineage == "Lineage 4"
            and next_best_lineage.sublineage != ""
            and next_best_lineage_percentage > self.threshold
        )


def load_lineages_from_csv(filepath: Path) -> Tuple[dict, dict]:
    """Load lineage metadata from a CVS file and return as a list of Lineages."""
    library = csv.DictReader(filepath.open())
    lineages = dict()
    position_map = dict()

    for entry in library:
        lineage = Lineage.from_csv_entry(entry)
        lineages[lineage.name] = lineage
        lineage_variants_file = LIBRARY_DIR / lineage.name
        if not lineage_variants_file.exists():
            print(
                "Lineage file {} does not exist for lineage {}".format(
                    str(lineage_variants_file), lineage.name
                ),
                file=sys.stderr,
            )
            continue

        lineage.add_snps(lineage_variants_file)

        for position, variant in lineage.snps.items():
            if position not in position_map:
                position_map[position] = {variant: [lineage.name]}
            elif variant in position_map[position]:
                position_map[position][variant].append(lineage.name)
            else:
                position_map[position][variant] = [lineage.name]

    return lineages, position_map


def output_results(outfile: IO[str], results: Dict[str, Tuple[float, Lineage]]):
    print("Sample\tSpecies\tLineage\tSublineage\tName\tPercentage", file=outfile)

    for sample_name, (percentage, lineage) in results.items():
        output = format_output_string(sample_name, round(percentage, 2), lineage)
        print(output, file=outfile)


def format_output_string(sample_name: str, percentage: float, lineage: Lineage) -> str:
    if not percentage:  # either there was no classification or was below threshold
        return "{1}\t{0}\t{0}\t{0}\t{0}\t0".format("N/A", sample_name)

    species = lineage.species or "N/A"
    lineage_name = lineage.lineage or "N/A"
    sublineage = lineage.sublineage or "N/A"
    name = lineage.name or "N/A"

    return "{sample_name}\t{species}\t{lineage_name}\t{sublineage}\t{name}\t{percentage}".format(
        sample_name=sample_name,
        species=species,
        lineage_name=lineage_name,
        sublineage=sublineage,
        name=name,
        percentage=percentage,
    )


def minos_gt_in_wrong_position_fix(record, sample_idx):
    """A version of minos had GT in the second column instead of the first"""
    info = str(record).strip().split("\t")[9 + sample_idx]
    for field in info.split(":"):
        if "/" in field:
            return Genotype.from_string(field)
