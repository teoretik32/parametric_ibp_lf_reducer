"""Tests for the label lattice: enumerate_box (arbitrary N/M), id map, complexity."""

from __future__ import annotations

import pytest

from parametric_ibp_lf_reducer import (
    LabelIndex,
    enumerate_box,
    label_complexity,
    make_label,
    split_label,
    zero_label,
)


def test_make_split_zero_roundtrip():
    lab = make_label([0, 1, 2], [-1, 3])
    assert lab == (0, 1, 2, -1, 3)
    n, m = split_label(lab, 3, 2)
    assert n == (0, 1, 2) and m == (-1, 3)
    assert zero_label(3, 2) == (0, 0, 0, 0, 0)


def test_split_length_mismatch_raises():
    with pytest.raises(ValueError):
        split_label((0, 1, 2), 2, 2)


def test_enumerate_box_broadcast_range_counts():
    # N=2, M=1: n in [0,1]^2 (2*2), m in [-1,0] (2) -> 8 labels of dimension 3.
    labels = list(enumerate_box(2, 1, (0, 1), (-1, 0)))
    assert len(labels) == 8
    assert all(len(lab) == 3 for lab in labels)
    assert (0, 0, -1) in labels and (1, 1, 0) in labels


def test_enumerate_box_arbitrary_dimensions_no_hardcode():
    # N=5, M=3 with per-axis and broadcast mixes; count must be the product of axis sizes.
    n_ranges = [(0, 1), (0, 0), (-1, 1), (0, 2), (0, 0)]  # sizes 2,1,3,3,1
    m_range = (-1, 0)  # size 2 per m-axis, 3 axes
    labels = list(enumerate_box(5, 3, n_ranges, m_range))
    expected = (2 * 1 * 3 * 3 * 1) * (2**3)
    assert len(labels) == expected
    assert all(len(lab) == 8 for lab in labels)


def test_enumerate_box_zero_polynomials():
    labels = list(enumerate_box(2, 0, (0, 1), (0, 0)))
    assert len(labels) == 4
    assert all(len(lab) == 2 for lab in labels)


def test_enumerate_box_bad_range_raises():
    with pytest.raises(ValueError):
        list(enumerate_box(2, 1, (2, 0), (0, 1)))  # lo>hi
    with pytest.raises(ValueError):
        list(enumerate_box(2, 1, [(0, 1)], (0, 1)))  # wrong per-axis count


def test_label_index_bijection():
    idx = LabelIndex.from_box(2, 1, (0, 1), (0, 1))
    assert len(idx) == 8
    for i, lab in enumerate(idx):
        assert idx.id(lab) == i
        assert idx.label(i) == lab
        assert lab in idx
    assert not idx.has((9, 9, 9))


def test_label_index_rejects_duplicates():
    with pytest.raises(ValueError):
        LabelIndex([(0, 0), (0, 0)])


def test_label_complexity_monotone():
    base = zero_label(2, 1)
    assert label_complexity(base, 2, 1) == 0
    # more positive n-shift -> larger complexity
    assert label_complexity((3, 0, 0), 2, 1) > label_complexity((1, 0, 0), 2, 1)
    # deeper negative m-depth is penalised more heavily than a positive m-shift
    assert label_complexity((0, 0, -2), 2, 1) > label_complexity((0, 0, 2), 2, 1)
