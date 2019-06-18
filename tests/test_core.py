import vcf
from snpit.genotype import Genotype
from snpit.core import *


def get_record(type):
    v = vcf.Reader(open("test_cases/test.vcf"))
    records = [record for record in v]

    if type == "ref":
        idx = 1
    elif type == "alt":
        idx = 0
    elif type == "null":
        idx = 3
    elif type == "het":
        idx = 2
    else:
        return None
    return records[idx]


def test_getSampleGenotypedVariant_refCallReturnsNone():
    record = get_record("ref")
    genotype = Genotype("0")

    actual = snpit.get_sample_genotyped_variant(genotype, record)
    expected = None

    assert actual == expected


def test_getSampleGenotypedVariant_altCallReturnsVariant():
    record = get_record("alt")
    genotype = Genotype("1")

    actual = snpit.get_sample_genotyped_variant(genotype, record)
    expected = "A"

    assert actual == expected


def test_getSampleGenotypedVariant_nullCallReturnsHyphen():
    record = get_record("null")
    genotype = Genotype(".")

    actual = snpit.get_sample_genotyped_variant(genotype, record)
    expected = "-"

    assert actual == expected


def test_getSampleGenotypedVariant_hetCallReturnsNone():
    record = get_record("het")
    genotype = Genotype("0", "1")

    actual = snpit.get_sample_genotyped_variant(genotype, record)
    expected = None

    assert actual == expected


def test_loadLineageMetadataFromFile_emptyFileReturnsEmptyDict():
    filepath = Path("test_cases/empty.csv")

    actual = load_lineage_metadata_from_file(filepath)
    expected = dict()

    assert actual == expected


def test_loadLineageMetadataFromFile_fileWithTwoEntriesReturnsDictWithTwoLineages():
    filepath = Path("test_cases/test_library.csv")

    actual = load_lineage_metadata_from_file(filepath)
    expected = {
        "indo-oceanic": dict(
            species="M. tuberculosis", lineage="Lineage 1", sublineage="Sublineage 7"
        ),
        "beijing": dict(species="M. tuberculosis", lineage="Lineage 2", sublineage=""),
    }

    assert actual == expected
