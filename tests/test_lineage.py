from snpit.lineage import Lineage
import csv
from pathlib import Path

TEST_CASE_DIR = Path("tests/test_cases")


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


def test_lessThanOperator_xLessThanYReturnsTrue():
    x = Lineage(name="bar")
    y = Lineage(name="foo")

    assert x < y


def test_lessThanOperator_xGreaterThanYReturnsFalse():
    x = Lineage(name="foo")
    y = Lineage(name="bar")

    assert not (x < y)


def test_lessThanOperator_xEqualsYReturnsFalse():
    x = Lineage(name="bar")
    y = Lineage(name="bar")

    assert not (x < y)


def test_greaterThanOperator_xLessThanYReturnsFalse():
    x = Lineage(name="bar")
    y = Lineage(name="foo")

    assert not (x > y)


def test_greaterThanOperator_xGreaterThanYReturnsTrue():
    x = Lineage(name="foo")
    y = Lineage(name="bar")

    assert x > y


def test_greaterThanOperator_xEqualsYReturnsFalse():
    x = Lineage(name="bar")
    y = Lineage(name="bar")

    assert not (x > y)


def test_addSnps_emptyFileSnpsRemainEmpty():
    lineage_variant_file = TEST_CASE_DIR / "empty.tsv"
    lineage = Lineage()
    lineage.add_snps(lineage_variant_file)

    actual = lineage.snps
    expected = dict()

    assert actual == expected


def test_addSnps_realFileSnpsContainAllEntries():
    lineage_variant_file = TEST_CASE_DIR / "test_lineage.tsv"
    lineage = Lineage()
    lineage.add_snps(lineage_variant_file)

    actual = lineage.snps
    expected = {1011511: "C", 1022003: "C", 1028217: "A", 1034758: "T", 1071966: "G"}

    assert actual == expected


def test_fromCsvEntry_emptyEntryReturnsEmptyLineage():
    entry = dict()

    actual = Lineage.from_csv_entry(entry)
    expected = Lineage()

    assert actual == expected


def test_fromCsvEntry_realEntryEntryReturnsLineage():
    library = csv.DictReader(TEST_CASE_DIR.joinpath("test_library.csv").open())
    entry = next(library)

    actual = Lineage.from_csv_entry(entry)
    expected = Lineage(
        species="M. tuberculosis",
        lineage="Lineage 1",
        sublineage="Sublineage 7",
        name="Indo_Oceanic",
    )

    assert actual == expected
