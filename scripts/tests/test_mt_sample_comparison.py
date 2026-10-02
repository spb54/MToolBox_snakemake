import numpy as np
import pytest

from scripts.mt_sample_comparison import (
    benjamini_hochberg, genome_wide_error, min_detectable_hf, read_sample_sheet,
    sam_differences, wilson_interval
)


def test_wilson_interval():
    low, high = wilson_interval(20, 100)
    assert low == pytest.approx(0.1334, abs=1e-3)
    assert high == pytest.approx(0.2888, abs=1e-3)
    # the interval narrows with depth at the same HF
    assert wilson_interval(200, 1000)[0] > low
    assert wilson_interval(0, 0) == (0.0, 1.0)


def test_min_detectable_hf_falls_with_depth():
    values = [min_detectable_hf(d, 0.03) for d in (50, 200, 1000, 5000)]
    assert values == sorted(values, reverse=True)
    assert values[-1] < 0.05 < values[0]
    assert min_detectable_hf(0, 0.03) == 1.0


def test_benjamini_hochberg():
    q = benjamini_hochberg([0.01, 0.04, 0.03, 0.5])
    assert list(np.round(q, 4)) == [0.04, 0.0533, 0.0533, 0.5]


def test_default_sample_groups():
    groups = read_sample_sheet(None, ["Keck_14175_tumour", "Keck_14175_cfDNA", "P2_normal"])
    assert groups["Keck_14175_tumour"] == ("Keck_14175", "tumour")
    assert groups["Keck_14175_cfDNA"] == ("Keck_14175", "cfDNA")
    assert groups["P2_normal"] == ("P2", "normal")


def test_genome_wide_error():
    # 3 positions x 1000 reads: 2 errors at the first two, a real 50%
    # variant at the third (excluded)
    plus = np.array([[998, 0, 500], [2, 0, 0], [0, 998, 500], [0, 2, 0]])
    counts = {"+": plus, "-": np.zeros_like(plus)}
    assert genome_wide_error(counts, 0.03) == pytest.approx(4 / 2000 / 3)


def test_sam_differences():
    mt = "ACGTACGTACGTACGTACGT"
    # aligned at 3 with a 1-base deletion; differences at mt 5 (A>T) and 13 (A>C)
    sam = ["numt1\t0\tchrM\t3\t60\t4M1D6M\t*\t0\t0\tGTTCTACGTC\t*"]
    diffs = sam_differences(sam, mt)
    assert set(diffs) == {(5, "A", "T"), (13, "A", "C")}
    assert diffs[(5, "A", "T")] == {"numt1"}
