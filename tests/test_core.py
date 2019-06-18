import vcf
from pathlib import Path
from snpit.genotype import Genotype
from snpit.lineage import Lineage
from snpit.core import load_lineages_from_csv, SnpIt


def get_record(record_type):
    v = vcf.Reader(open("test_cases/test.vcf"))
    records = [record for record in v]

    if record_type == "ref":
        idx = 1
    elif record_type == "alt":
        idx = 0
    elif record_type == "null":
        idx = 3
    elif record_type == "het":
        idx = 2
    else:
        return None
    return records[idx]


def test_getSampleGenotypedVariant_refCallReturnsNone():
    record = get_record("ref")
    genotype = Genotype("0")

    actual = SnpIt.get_variant_for_genotype_in_vcf_record(genotype, record)
    expected = None

    assert actual == expected


def test_getSampleGenotypedVariant_altCallReturnsVariant():
    record = get_record("alt")
    genotype = Genotype("1")

    actual = SnpIt.get_variant_for_genotype_in_vcf_record(genotype, record)
    expected = "A"

    assert actual == expected


def test_getSampleGenotypedVariant_nullCallReturnsHyphen():
    record = get_record("null")
    genotype = Genotype(".")

    actual = SnpIt.get_variant_for_genotype_in_vcf_record(genotype, record)
    expected = "-"

    assert actual == expected


def test_getSampleGenotypedVariant_hetCallReturnsNone():
    record = get_record("het")
    genotype = Genotype("0", "1")

    actual = SnpIt.get_variant_for_genotype_in_vcf_record(genotype, record)
    expected = None

    assert actual == expected


def test_loadLineagesFromCsv_emptyFileReturnsEmptyDict():
    filepath = Path("test_cases/empty.csv")

    actual = load_lineages_from_csv(filepath)
    expected = []

    assert actual == expected


def test_loadLineagesFromCsv_fileWithTwoEntriesReturnsDictWithTwoLineages():
    filepath = Path("test_cases/test_library.csv")

    actual = load_lineages_from_csv(filepath)
    expected = [
        Lineage(
            name="indo-oceanic",
            species="M. tuberculosis",
            lineage="Lineage 1",
            sublineage="Sublineage 7",
        ),
        Lineage(
            name="beijing",
            species="M. tuberculosis",
            lineage="Lineage 2",
            sublineage="",
        ),
    ]

    assert actual == expected


def test_classifyVcf_exampleVcfReturnsCorrectClassification():
    snpit = SnpIt(10, True)

    actual = snpit.classify_vcf("../example/example.vcf")
    expected = {"example": ("M. tuberculosis", "Lineage 2", "", 97.368_421_052_631_58)}

    assert actual == expected
