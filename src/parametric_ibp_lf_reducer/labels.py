"""Label lattice utilities.

A *label* is an integer tuple of length ``N + M`` encoding the discrete shift of a family
member relative to the base integrand (spec §2)::

    label = (n_1, ..., n_N, m_1, ..., m_M)
    F_label = F_base * prod_i x_i^(n_i) * prod_l G_l^(m_l)

Nothing here is hardcoded to a particular ``N`` or ``M``; all sizes come from arguments. This
module provides label construction/splitting, a search-box enumeration, an id map, and a
structural complexity metric (a utility — it does NOT decide masters; ranking is a later pass).
"""

from __future__ import annotations

import itertools
from collections.abc import Iterable, Iterator
from dataclasses import dataclass

Label = tuple[int, ...]
Range = tuple[int, int]


def make_label(n: Iterable[int], m: Iterable[int]) -> Label:
    return (*(int(x) for x in n), *(int(x) for x in m))


def split_label(label: Label, nvars: int, npolys: int) -> tuple[Label, Label]:
    if len(label) != nvars + npolys:
        raise ValueError(f"label length {len(label)} != nvars+npolys = {nvars + npolys}")
    return tuple(label[:nvars]), tuple(label[nvars:])


def zero_label(nvars: int, npolys: int) -> Label:
    return (0,) * (nvars + npolys)


def _normalize_ranges(spec, count: int, what: str) -> list[Range]:
    if count == 0:
        return []
    # A single (lo, hi) pair broadcast to every axis...
    if (
        isinstance(spec, tuple)
        and len(spec) == 2
        and all(isinstance(x, int) for x in spec)
    ):
        pairs = [(int(spec[0]), int(spec[1]))] * count
    else:  # ...or an explicit per-axis sequence of (lo, hi) pairs.
        seq = list(spec)
        if len(seq) != count:
            raise ValueError(f"{what} ranges: expected {count} pairs, got {len(seq)}")
        pairs = [(int(lo), int(hi)) for lo, hi in seq]
    for lo, hi in pairs:
        if lo > hi:
            raise ValueError(f"{what} range has lo>hi: ({lo},{hi})")
    return pairs


def enumerate_box(nvars: int, npolys: int, n_range, m_range) -> Iterator[Label]:
    """Enumerate labels over a rectangular box of ``n``- and ``m``-shifts.

    ``n_range``/``m_range`` are each either a single ``(lo, hi)`` pair (broadcast to every axis)
    or a per-axis sequence of ``(lo, hi)`` pairs. Works for any ``nvars >= 1`` and ``npolys >= 0``.
    """
    if nvars < 1:
        raise ValueError("nvars must be >= 1")
    if npolys < 0:
        raise ValueError("npolys must be >= 0")
    n_ranges = _normalize_ranges(n_range, nvars, "n")
    m_ranges = _normalize_ranges(m_range, npolys, "m")
    axes = [range(lo, hi + 1) for lo, hi in n_ranges] + [range(lo, hi + 1) for lo, hi in m_ranges]
    for combo in itertools.product(*axes):
        yield tuple(combo)


@dataclass
class LabelIndex:
    """A stable bijection between labels and contiguous integer ids."""

    labels: list[Label]

    def __post_init__(self) -> None:
        self.labels = list(self.labels)
        self._id: dict[Label, int] = {}
        for i, lab in enumerate(self.labels):
            if lab in self._id:
                raise ValueError(f"duplicate label {lab}")
            self._id[lab] = i

    @classmethod
    def from_box(cls, nvars: int, npolys: int, n_range, m_range) -> "LabelIndex":
        return cls(list(enumerate_box(nvars, npolys, n_range, m_range)))

    def id(self, label: Label) -> int:
        return self._id[label]

    def has(self, label: Label) -> bool:
        return label in self._id

    def label(self, idx: int) -> Label:
        return self.labels[idx]

    def __len__(self) -> int:
        return len(self.labels)

    def __contains__(self, label: Label) -> bool:
        return label in self._id

    def __iter__(self) -> Iterator[Label]:
        return iter(self.labels)


# Default structural weights for label complexity (utility only, not a ranking decision).
DEFAULT_COMPLEXITY_WEIGHTS = {
    "pos_n": 1.0,  # total positive n-shift (raising integration-variable powers)
    "neg_m": 2.0,  # total negative m-depth (deeper polynomial denominators)
    "pos_m": 1.0,  # total positive m-shift
    "abs_n": 0.5,  # total |n| spread
}


def label_complexity(
    label: Label, nvars: int, npolys: int, weights: dict[str, float] | None = None
) -> float:
    """A structural size proxy for a label. Higher = more complex.

    This is a plain metric over the shift lattice; it deliberately does not look at local
    finiteness or master status (those belong to valuations/ranking in later passes).
    """
    w = weights or DEFAULT_COMPLEXITY_WEIGHTS
    n, m = split_label(label, nvars, npolys)
    pos_n = sum(x for x in n if x > 0)
    abs_n = sum(abs(x) for x in n)
    neg_m = sum(-x for x in m if x < 0)
    pos_m = sum(x for x in m if x > 0)
    return (
        w["pos_n"] * pos_n
        + w["abs_n"] * abs_n
        + w["neg_m"] * neg_m
        + w["pos_m"] * pos_m
    )
