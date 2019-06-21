#! /usr/bin/env python
import csv
import gzip
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Tuple, List, Iterator
import numpy as np
import pysam
from Bio import SeqIO

from .genotype import Genotype, UnexpectedGenotypeError
from .lineage import Lineage

LIBRARY_DIR = Path(__file__).parent.parent / "lib"


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

    def add_snps_to_all_lineages(self):
        """Add all of the variants that make up a lineage to the Lineage objects."""
        for lineage_name in self.lineages:
            lineage_variants_file = LIBRARY_DIR / lineage_name

            if not lineage_variants_file.exists():
                print(
                    "Lineage file {} does not exist for lineage {}".format(
                        lineage_variants_file, lineage_name
                    ),
                    file=sys.stderr,
                )
                continue

            self.lineages[lineage_name].add_snps(lineage_variants_file)

    def classify_vcf(self, vcf_path: Path) -> dict:
        """Loads the vcf file and then, for each lineage, identify the base at each of
        the identifying positions in the genome.

        Args:
            vcf_file: Path to the VCF file to be read
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
        vcf = pysam.VariantFile(vcf_path)

        sample_lineage_counts = defaultdict(Counter)

        for record in vcf:
            if record.pos not in self.lineage_positions:
                continue
            if not self.ignore_filter and "PASS" not in record.filter.keys():
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

                if genotype.is_reference() or genotype.is_heterozygous():
                    continue
                elif genotype.is_alt():
                    alt_call = max(genotype.call())
                    variant = record.alleles[alt_call]
                elif genotype.is_null():
                    variant = "-"
                else:
                    raise UnexpectedGenotypeError(
                        f"""Got a genotype for which a Ref/Alt/Null call could not be 
                            determined: {genotype.call()}.\nPlease raise this with the 
                            developers."""
                    )

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

    def classify_record(self, record):
        """Determine what lineage(s) (if any) a VCF record belongs to.

        Args:
            record (vcf.model._Record): A VCF record object.
        """

        classifications = defaultdict(Counter)
        for lineage in self.lineages.values():
            if record.POS not in lineage.snps:
                continue

            for call in record.samples:
                genotype = Genotype.from_string(call["GT"])
                variant = self.get_variant_for_genotype_in_vcf_record(genotype, record)
                is_variant_shared_with_lineage = (
                    variant is not None and variant == lineage.snps[record.POS]
                )

                if is_variant_shared_with_lineage:
                    classifications[call.sample].update([lineage.name])

        return classifications

    @staticmethod
    def get_variant_for_genotype_in_vcf_record(genotype, record):
        """Retrieves the variant a genotype maps to for a given record.

        Args:
            genotype (Genotype): The genotype call for the sample.
            record (vcf.model._Record): A VCF record object.
        Returns:
            str or None: A hyphen if the call is null (ie ./.) or the alt variant if
            the call is alt. Returns None is the call is ref or heterozygous.
        """
        if genotype.is_reference() or genotype.is_heterozygous():
            return None
        elif genotype.is_null():
            # record a hyphen which won't match, regardless of the reference
            return "-"
        elif genotype.is_alt():
            # replace the H37Rv base with the actual base from the VCF file
            alt_variant = record.ALT[int(genotype.call()[0]) - 1]
            return alt_variant

    def load_fasta(self, fasta_file, compression=False):
        """
        Loads a supplied fasta file and then, for each lineage, identify the base at
        each of the identifying positions

        Args:
         fasta_file (str): Path to the fasta file to be read
         compression (bool): whether the fasta file is compressed by gz or bzip2
        """

        # setup the dictionaries of expected SNPs for each lineage
        # self._reset_lineage_snps()
        self.sample_snps = {}

        # open the fasta file for reading
        open_fn = gzip.open if compression else open

        with open_fn(fasta_file, "rt") as fasta_file:
            fasta_reader = SeqIO.read(fasta_file, "fasta")

        #  iterate through the lineages
        for lineage_name in self.lineages:

            self.sample_snps[lineage_name] = {}

            # iterate over the positions in the reference set of snps for that lineage
            for pos in self.reference_snps[lineage_name]:
                if pos in self.reference_snps[lineage_name].keys():

                    # CAUTION the GenBank File is 1-based, but the lineage files are 0-based
                    # Remember the nucleotide at the defining position
                    self.sample_snps[lineage_name][int(pos)] = fasta_reader.seq[
                        int(pos) - 1
                    ]

    # todo: clean up this function
    def determine_lineage(self, lineage_counts: Counter) -> Tuple[float, Lineage]:
        """
        Having read the VCF file, for each lineage, calculate the percentage of SNP
        present in the sample.
        Note that this means the percentages will not add up to 100%.

        Returns:
            tuple of (lineage,percentage)
        """
        percentage_shared_snps_for_lineages = []

        for name, lineage in self.lineages.items():
            shared_calls = lineage_counts[name]
            ref_calls = len(lineage.snps)

            # thereby calculate the percentage of SNPs in this sample that match the lineage
            percentage_shared_snps_for_lineages.append(
                ((shared_calls / ref_calls) * 100, lineage)
            )

        percentage_shared_snps_for_lineages.sort(reverse=True)

        identified_lineage = percentage_shared_snps_for_lineages[0][1]
        identified_lineage_percentage = percentage_shared_snps_for_lineages[0][0]

        # if the top lineage is above the specified threshold, return the classification
        if identified_lineage_percentage > self.threshold:

            # look at the next-highest lineage if the top one is Lineage 4 but with no sublineage
            if (
                identified_lineage.lineage == "Lineage 4"
                and identified_lineage.sublineage == ""
            ):

                next_best_lineage = percentage_shared_snps_for_lineages[1][1]
                next_best_lineage_percentage = percentage_shared_snps_for_lineages[1][0]

                # if the next best lineage is ALSO lineage 4, but this one has a
                # sublineage and is above the threshold, report that one instead
                if (
                    next_best_lineage.lineage == "Lineage 4"
                    and next_best_lineage.sublineage != ""
                    and next_best_lineage_percentage > self.threshold
                ):

                    identified_lineage = next_best_lineage

            return identified_lineage_percentage, identified_lineage
        else:
            return 0, Lineage()


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
                    lineage_variants_file, lineage.name
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


def output_results(outfile, results):
    print("Sample\tSpecies\tLineage\tSublineage\tName\tPercentage", file=outfile)

    for sample_name, (percentage, lineage) in results.items():
        output = format_output_string(sample_name, percentage, lineage)
        print(output, file=outfile)


def format_output_string(sample_name, percentage, lineage):
    if not percentage:  # either there was no classification or was below threshold
        return "{1}\t{0}\t{0}\t{0}\t{0}\t{0}".format("N/A", sample_name)

    return "{sample}\t{species}\t{lineage}\t{sublineage}\t{name}\t{percentage}".format(
        species=lineage.species,
        lineage=(lineage.lineage or "N/A"),
        sublineage=(lineage.sublineage or "N/A"),
        percentage=round(percentage, 2),
        name=lineage.name,
        sample=sample_name,
    )


def minos_gt_in_wrong_position_fix(record, sample_idx):
    """A version of minos had GT in the second column instead of the first"""
    info = str(record).strip().split("\t")[9 + sample_idx]
    for field in info.split(":"):
        if "/" in field:
            return Genotype.from_string(field)
