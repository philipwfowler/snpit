import pytest

from snpit.genotype import Genotype, InvalidGenotypeString


def test_equalityOperator_returnTrueForTwoIdenticalGenotypes():
    assert Genotype(0) == Genotype(0)


def test_equalityOperator_returnTrueForTwoIdenticalGenotypesInDifferentOrder():
    assert Genotype(1, 0) == Genotype(0, 1)


def test_equalityOperator_returnTrueForHaploidRefAndDiploidRef():
    assert Genotype(0) == Genotype(0, 0)


def test_equalityOperator_returnTrueForNullCalls():
    assert Genotype(-1, -1) == Genotype(-1)


def test_equalityOperator_returnFalseForRefVsAltCalls():
    assert Genotype(0) != Genotype(1)


def test_equalityOperator_returnFalseForRefVsHetCalls():
    assert Genotype(0) != Genotype(1, 0)


def test_callReturnsTupleWithNone():
    g = Genotype(0)

    actual = g.call()
    expected = [0, 0]

    assert actual == expected


def test_callReturnsTupleWithBothCalls():
    g = Genotype(0, 1)

    actual = g.call()
    expected = [0, 1]

    assert actual == expected


def test_isNull_returnsTrueForFullStopCall():
    g = Genotype(-1)

    assert g.is_null()


def test_isNull_returnsFalseForNumberedCall():
    g = Genotype(0)

    assert not g.is_null()


def test_isNull_returnsTrueForNoneCall():
    g = Genotype(-1)

    assert g.is_null()


def test_isHeterozygous_returnsTrueForRefAltHetCall():
    g = Genotype(0, 1)

    assert g.is_heterozygous()


def test_isHeterozygous_returnsTrueForAltRefCall():
    g = Genotype(4, 0)

    assert g.is_heterozygous()


def test_isHeterozygous_returnsFalseForRefCall():
    g = Genotype(0, 0)

    assert not g.is_heterozygous()


def test_isHeterozygous_returnsFalseForAltCall():
    g = Genotype(1, 1)

    assert not g.is_heterozygous()


def test_isHeterozygous_returnsFalseForSingleCall():
    g = Genotype(1)

    assert not g.is_heterozygous()


def test_isReference_returnsTrueForRefNoneCall():
    g = Genotype(0)

    assert g.is_reference()


def test_isReference_returnsTrueForRefRefCall():
    g = Genotype(0, 0)

    assert g.is_reference()


def test_isReference_returnsFalseForAltNoneCall():
    g = Genotype(1)

    assert not g.is_reference()


def test_isReference_returnsFalseForHetCall():
    g = Genotype(0, 1)

    assert not g.is_reference()


def test_isReference_returnsFalseForNullCall():
    g = Genotype(-1, -1)

    assert not g.is_reference()


def test_isReference_returnsFalseForAltAltCall():
    g = Genotype(1, 1)

    assert not g.is_reference()


def test_isAlt_returnsFalseForRefNoneCall():
    g = Genotype(0)

    assert not g.is_alt()


def test_isAlt_returnsFalseForRefRefCall():
    g = Genotype(0, 0)

    assert not g.is_alt()


def test_isAlt_returnsTrueForAltNoneCall():
    g = Genotype(1)

    assert g.is_alt()


def test_isAlt_returnsFalseForHetCall():
    g = Genotype(0, 1)

    assert not g.is_alt()


def test_isAlt_returnsFalseForNullCall():
    g = Genotype(-1, -1)

    assert not g.is_alt()


def test_isAlt_returnsTrueForAltAltCall():
    g = Genotype(1, 1)

    assert g.is_alt()


def test_fromString_validStringReturnsDiploidRefGenotype():
    actual = Genotype.from_string("0/0")
    expected = Genotype(0, 0)

    assert actual == expected


def test_fromString_validStringReturnsDiploidHetGenotype():
    actual = Genotype.from_string("1/0")
    expected = Genotype(0, 1)

    assert actual == expected


def test_fromString_validStringReturnsDiploidAltGenotype():
    actual = Genotype.from_string("1/1 ")
    expected = Genotype(1, 1)

    assert actual == expected


def test_fromString_validStringReturnsDiploidNullGenotype():
    actual = Genotype.from_string("./.")
    expected = Genotype(-1, -1)

    assert actual == expected


def test_fromString_validStringReturnsHaploidGenotype():
    actual = Genotype.from_string("0/")
    expected = Genotype(0)

    assert actual == expected


def test_fromString_invalidStringNoSlashRaisesError():
    with pytest.raises(InvalidGenotypeString):
        g = Genotype.from_string("0")


def test_fromString_invalidStringTooManySlashesRaisesError():
    with pytest.raises(InvalidGenotypeString):
        g = Genotype.from_string("0/0/0")
