from pathlib import Path


class Lineage:
    def __init__(self, name="", species="", lineage="", sublineage=""):
        self.name = name
        self.species = species
        self.lineage = lineage
        self.sublineage = sublineage
        self.snps = dict()

    def __eq__(self, other):
        return self.name == other.name

    def __lt__(self, other):
        return self.name < other.name

    def __gt__(self, other):
        return self.name > other.name

    def __hash__(self):
        return hash(self.name)

    def add_snps(self, lineage_variants_file: Path):
        """Adds SNP information to lineage from a variant file. The variant file must
        be in the format position\tbase
        """
        with lineage_variants_file.open() as lineage_file:
            for line in lineage_file:
                lineage_variant = line.rstrip().split("\t")
                position = int(lineage_variant[0])
                base = lineage_variant[1]

                self.snps[position] = base

    @staticmethod
    def from_csv_entry(entry: dict) -> "Lineage":
        """Create a Lineage object from a CSV DictReader entry. Designed to work with
        the package lib/library.csv file.
        """
        return Lineage(
            species=entry.get("species", ""),
            lineage=entry.get("lineage", ""),
            sublineage=entry.get("sublineage", ""),
            name=entry.get("id", ""),
        )
