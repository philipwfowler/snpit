from collections import Counter, defaultdict
from io import StringIO
from pathlib import Path

import pysam
import pytest

from snpit.snpit import load_lineages_from_csv, SnpIt, output_results
from snpit.genotype import Genotype
from snpit.lineage import Lineage

TEST_CASE_DIR = Path("tests/test_cases")


def get_record(record_type):
    test_vcf = TEST_CASE_DIR / "test.vcf"
    vcf = pysam.VariantFile(str(test_vcf))
    records = [record for record in vcf]

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


def create_test_snpit():
    positions = [2, 5, 8, 11, 14, 16, 20]
    names = [
        [("L1", "A")],
        [("L2", "C")],
        [("L1", "A"), ("L4", "C")],
        [("L2", "T"), ("L3", "G")],
        [("L3", "A")],
        [("L4", "G")],
        [("L1", "C"), ("L3", "T"), ("L4", "G")],
    ]
    lineages = {
        "L1": Lineage(name="L1", species="S1", lineage="Lin1", sublineage="SL1"),
        "L2": Lineage(name="L2", species="S2", lineage="Lin2", sublineage="SL2"),
        "L3": Lineage(name="L3", species="S3", lineage="Lin3", sublineage="SL3"),
        "L4": Lineage(name="L4", species="S4", lineage="Lineage 4", sublineage="SL4"),
    }
    position_map = dict()

    for i in range(len(positions)):
        position_map[positions[i]] = {}
        for name, variant in names[i]:
            if variant not in position_map[positions[i]]:
                position_map[positions[i]][variant] = [name]
            else:
                position_map[positions[i]][variant].append(name)

    snpit = SnpIt(threshold=10)
    snpit.lineages = lineages
    snpit.lineage_positions = position_map
    return snpit


def test_getSampleGenotypedVariant_refCallReturnsEmpty():
    record = get_record("ref")
    genotype = Genotype(0)

    actual = SnpIt.get_variant_for_genotype_in_vcf_record(genotype, record)
    expected = ""

    assert actual == expected


def test_getSampleGenotypedVariant_altCallReturnsVariant():
    record = get_record("alt")
    genotype = Genotype(1)

    actual = SnpIt.get_variant_for_genotype_in_vcf_record(genotype, record)
    expected = "A"

    assert actual == expected


def test_getSampleGenotypedVariant_nullCallReturnsHyphen():
    record = get_record("null")
    genotype = Genotype(None)

    actual = SnpIt.get_variant_for_genotype_in_vcf_record(genotype, record)
    expected = "-"

    assert actual == expected


def test_getSampleGenotypedVariant_hetCallReturnsNone():
    record = get_record("het")
    genotype = Genotype(0, 1)

    actual = SnpIt.get_variant_for_genotype_in_vcf_record(genotype, record)
    expected = ""

    assert actual == expected


def test_loadLineagesFromCsv_emptyFileReturnsEmptyDict():
    filepath = TEST_CASE_DIR / "empty.csv"

    actual = load_lineages_from_csv(filepath)
    expected = (dict(), dict())

    assert actual == expected


def test_loadLineagesFromCsv_fileWithTwoEntriesReturnsDictWithTwoLineages():
    filepath = TEST_CASE_DIR / "test_library.csv"
    beijing_lineage = Lineage(
        name="beijing", species="M. tuberculosis", lineage="Lineage 2", sublineage=""
    )
    indo_lineage = Lineage(
        name="Indo_Oceanic",
        species="M. tuberculosis",
        lineage="Lineage 1",
        sublineage="Sublineage 7",
    )

    actual_lineages, actual_position_map = load_lineages_from_csv(filepath)

    assert len(actual_position_map) == 714
    assert actual_position_map[1083755] == {"G": [indo_lineage.name]}
    assert actual_position_map[1288698] == {"A": [beijing_lineage.name]}

    assert actual_lineages["beijing"] == beijing_lineage
    assert actual_lineages["Indo_Oceanic"] == indo_lineage


def test_classifyVcf_exampleVcfReturnsCorrectClassification():
    snpit = SnpIt(10, True)
    vcf_path = TEST_CASE_DIR / "example.vcf"

    actual = snpit.classify_vcf(vcf_path)
    expected = {
        "example": (
            97.2972972972973,
            Lineage(lineage="Lineage 2", species="M. tuberculosis", name="beijing"),
        )
    }

    assert actual == expected


def test_countLineageClassificationsForSamplesInVcfWithNoSharedPositions_ReturnEmpty():
    snpit = create_test_snpit()
    snpit.ignore_filter = True
    snpit.ignore_status = True
    vcf_path = TEST_CASE_DIR / "empty_multisample.vcf"

    actual = snpit.count_lineage_classifications_for_samples_in_vcf(vcf_path)
    expected = defaultdict()
    expected["Sample1"] = Counter([])
    expected["Sample2"] = Counter([])

    assert actual == expected


def test_countLineageClassificationsForSamplesInVcf_ignoreFilterIgnoreStatus():
    snpit = create_test_snpit()
    snpit.ignore_filter = True
    snpit.ignore_status = True
    vcf_path = TEST_CASE_DIR / "multisample_test.vcf"

    actual = snpit.count_lineage_classifications_for_samples_in_vcf(vcf_path)
    expected = defaultdict()
    expected["Sample1"] = Counter(["L4", "L3", "L3"])
    expected["Sample2"] = Counter(["L1", "L1", "L1", "L2", "L3", "L3"])

    assert actual == expected


def test_countLineageClassificationsForSamplesInVcf_notIgnoreFilterIgnoreStatus():
    snpit = create_test_snpit()
    snpit.ignore_filter = False
    snpit.ignore_status = True
    vcf_path = TEST_CASE_DIR / "multisample_test.vcf"

    actual = snpit.count_lineage_classifications_for_samples_in_vcf(vcf_path)
    expected = defaultdict()
    expected["Sample1"] = Counter(["L3", "L3"])
    expected["Sample2"] = Counter(["L1", "L2", "L3"])

    assert actual == expected


def test_countLineageClassificationsForSamplesInVcf_ignoreFilterNotIgnoreStatus():
    snpit = create_test_snpit()
    snpit.ignore_filter = True
    snpit.ignore_status = False
    vcf_path = TEST_CASE_DIR / "multisample_test.vcf"

    actual = snpit.count_lineage_classifications_for_samples_in_vcf(vcf_path)
    expected = defaultdict()
    expected["Sample1"] = Counter(["L4", "L3", "L3"])
    expected["Sample2"] = Counter(["L1", "L2", "L3", "L3"])

    assert actual == expected


def test_countLineageClassificationsForSamplesInVcf_notIgnoreFilterNotIgnoreStatus():
    snpit = create_test_snpit()
    snpit.ignore_filter = False
    snpit.ignore_status = False
    vcf_path = TEST_CASE_DIR / "multisample_test.vcf"

    actual = snpit.count_lineage_classifications_for_samples_in_vcf(vcf_path)
    expected = defaultdict()
    expected["Sample1"] = Counter(["L3", "L3"])
    expected["Sample2"] = Counter(["L2", "L3"])

    assert actual == expected


def test_determineLineage_emptyCountsInReturnsEmptyResults():
    snpit = create_test_snpit()
    snpit.lineages = {}

    actual_percentage, actual_lineage = snpit.determine_lineage(Counter())
    expected_percentage, expected_lineage = (0, Lineage())

    assert actual_percentage == expected_percentage
    assert actual_lineage == actual_lineage


def test_determineLineage_L3MostCountsReturnsL3():
    snpit = create_test_snpit()
    dummy_snps = {x: "A" for x in range(10)}
    snpit.lineages["L3"].snps = dummy_snps
    snpit.lineages["L1"].snps = dummy_snps

    actual_percentage, actual_lineage = snpit.determine_lineage(
        Counter(["L3", "L1", "L3"])
    )
    expected_percentage, expected_lineage = (
        20.0,
        Lineage(name="L3", species="S3", lineage="Lin3", sublineage="SL3"),
    )

    assert actual_percentage == expected_percentage
    assert actual_lineage == actual_lineage


def test_determineLineage_L4MostCountsReturnsL4():
    snpit = create_test_snpit()
    dummy_snps = {x: "A" for x in range(10)}
    snpit.lineages["L4"].snps = dummy_snps
    snpit.lineages["L1"].snps = dummy_snps

    actual_percentage, actual_lineage = snpit.determine_lineage(
        Counter(["L4", "L1", "L4"])
    )
    expected_percentage, expected_lineage = (
        20.0,
        Lineage(name="L4", species="S3", lineage="Lin3", sublineage="SL3"),
    )

    assert actual_percentage == expected_percentage
    assert actual_lineage == actual_lineage


def test_determineLineage_L4NoSublineageMostCountsReturnsL4():
    snpit = create_test_snpit()
    dummy_snps = {x: "A" for x in range(10)}
    snpit.lineages["L4"].snps = dummy_snps
    snpit.lineages["L1"].snps = dummy_snps
    snpit.lineages["L4"].sublineage = ""

    actual_percentage, actual_lineage = snpit.determine_lineage(
        Counter(["L4", "L1", "L4"])
    )
    expected_percentage, expected_lineage = (
        20.0,
        Lineage(name="L4", species="S3", lineage="Lineage 4"),
    )

    assert actual_percentage == expected_percentage
    assert actual_lineage == actual_lineage


def test_determineLineage_L4NoSublineageMostCountsNextBestL4WithSublineageReturnsNextBest():
    snpit = create_test_snpit()
    dummy_snps = {x: "A" for x in range(10)}
    snpit.lineages["L4"].snps = dummy_snps
    snpit.lineages["L1"].snps = dummy_snps
    snpit.lineages["L4"].sublineage = ""
    snpit.lineages["L1"].lineage = "Lineage 4"
    snpit.lineages["L1"].sublineage = "corner case"

    actual_percentage, actual_lineage = snpit.determine_lineage(
        Counter(["L4", "L1", "L4"])
    )
    expected_percentage, expected_lineage = (
        20.0,
        Lineage(name="L1", species="S1", lineage="Lineage 4", sublineage="corner case"),
    )

    assert actual_percentage == expected_percentage
    assert actual_lineage == actual_lineage


def test_outputResults_emptyResultswritesJustTheHeader():
    outfile = StringIO()
    results = dict()
    output_results(outfile, results)
    outfile.seek(0)

    actual = outfile.read()
    expected = "Sample\tSpecies\tLineage\tSublineage\tName\tPercentage\n"

    assert actual == expected


def test_outputResults_emptySampleResultsWritesSampleWithNAs():
    outfile = StringIO()
    results = dict(Sample1=(0, Lineage()))
    output_results(outfile, results)
    outfile.seek(0)

    actual = outfile.read()
    expected = (
        "Sample\tSpecies\tLineage\tSublineage\tName\tPercentage\n"
        "Sample1\tN/A\tN/A\tN/A\tN/A\t0\n"
    )

    assert actual == expected


def test_outputResults_twoSamplesResultsWritesTwoSamples():
    outfile = StringIO()
    results = {
        "Sample1": (25.2533333, Lineage(name="L1", sublineage="SL1")),
        "Sample2": (75.75899, Lineage(name="L2", species="S2", lineage="Lin2")),
    }
    output_results(outfile, results)
    outfile.seek(0)

    actual = outfile.read()
    expected = (
        "Sample\tSpecies\tLineage\tSublineage\tName\tPercentage\n"
        "Sample1\tN/A\tN/A\tSL1\tL1\t25.25\n"
        "Sample2\tS2\tLin2\tN/A\tL2\t75.76\n"
    )

    assert actual == expected


def test_classifyFasta_emptyFastaRaisesOSError():
    empty_fasta = TEST_CASE_DIR / "empty.fa"
    snpit = create_test_snpit()

    with pytest.raises(OSError):
        snpit.classify_fasta(empty_fasta)


def test_classifyFasta_singleSampleFastaReturnSingleSample():
    fasta_path = TEST_CASE_DIR / "single_sample.fa"
    snpit = create_test_snpit()
    snpit.lineages["L1"].snps = {x: "A" for x in range(10)}
    snpit.lineages["L3"].snps = {x: "A" for x in range(20)}
    snpit.lineages["L4"].snps = {x: "A" for x in range(2)}

    actual = snpit.classify_fasta(fasta_path)
    expected = {"Sample1": (50.0, Lineage(name="L4"))}

    assert actual == expected


def test_classifyFasta_multiSampleFastaReturnMultiSample():
    fasta_path = TEST_CASE_DIR / "multi_sample.fa"
    snpit = create_test_snpit()
    snpit.lineages["L1"].snps = {x: "A" for x in range(10)}
    snpit.lineages["L2"].snps = {x: "A" for x in range(20)}
    snpit.lineages["L3"].snps = {x: "A" for x in range(20)}
    snpit.lineages["L4"].snps = {x: "A" for x in range(2)}

    actual = snpit.classify_fasta(fasta_path)
    expected = {
        "Sample1": (50.0, Lineage(name="L4")),
        "Sample2": (20.0, Lineage(name="L1")),
    }

    assert actual == expected
