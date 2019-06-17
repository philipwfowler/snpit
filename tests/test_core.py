import vcf
from snpit.genotype import Genotype
from snpit.core import snpit


def get_record(type):
    v = vcf.Reader(open("test.vcf"))
    records = [record for record in v]
    idx = 0
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
