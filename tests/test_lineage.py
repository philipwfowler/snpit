from snpit.lineage import Lineage
import csv
from pathlib import Path


def test_equalityOperator_twoEqualReturnsTrue():
    lineage1 = Lineage(name="test", lineage="Lineage 5")
    lineage2 = Lineage(name="test", lineage="Lineage 6")

    assert lineage1 == lineage2


def test_equalityOperator_twoNonEqualReturnsFalse():
    lineage1 = Lineage(name="foo", lineage="Lineage 5")
    lineage2 = Lineage(name="test", lineage="Lineage 6")

    assert not (lineage1 == lineage2)


def test_inequalityOperator_twoEqualReturnsFalse():
    lineage1 = Lineage(name="test", lineage="Lineage 5")
    lineage2 = Lineage(name="test", lineage="Lineage 6")

    assert not (lineage1 != lineage2)


def test_inequalityOperator_twoNonEqualReturnsTrue():
    lineage1 = Lineage(name="foo", lineage="Lineage 5")
    lineage2 = Lineage(name="test", lineage="Lineage 6")

    assert lineage1 != lineage2


def test_fromCsvEntry_emptyEntryReturnsEmptyLineage():
    entry = dict()

    actual = Lineage.from_csv_entry(entry)
    expected = Lineage()

    assert actual == expected


def test_fromCsvEntry_realEntryEntryReturnsLineage():
    library = csv.DictReader(Path("test_cases/test_library.csv").open())
    entry = next(library)

    actual = Lineage.from_csv_entry(entry)
    expected = Lineage(
        species="M. tuberculosis",
        lineage="Lineage 1",
        sublineage="Sublineage 7",
        name="indo-oceanic",
    )

    assert actual == expected
