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
