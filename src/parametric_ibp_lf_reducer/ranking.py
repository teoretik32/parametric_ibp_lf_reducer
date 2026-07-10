"""Label ranking: the column elimination order for the reduction (spec §5.8, method review §6).

Ranking assigns each label an *elimination priority* used later to choose pivot columns in the
modular RREF. It does NOT solve anything and it NEVER removes a label — it only orders them.

Ordering is ``(tier, -complexity, label)`` ascending, so index 0 is the highest-priority pivot
(eliminate first) and the tail is what we prefer to leave *free* (the masters):

    tier 0 : the target label            -> eliminated first (we want its normal form)
    tier 1 : NON locally finite / Unknown -> eliminated before any LF label (never left free)
    tier 2 : generic locally finite       -> complex ones eliminated, simplest kept free
    tier 3 : user-preferred LF masters     -> kept free with highest preference (ordered last)

Within a tier, higher structural complexity is eliminated earlier, so the simplest locally
finite integrands drift to the free tail. A *simple* non-LF label still ranks in tier 1, so it
is eliminated before any LF master — simplicity never rescues a divergent integral.
"""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass

from .family import ParametricFamily
from .labels import Label, label_complexity
from .valuations import is_locally_finite

# tier constants
TIER_TARGET = 0
TIER_NON_LF = 1
TIER_LF = 2
TIER_PREFERRED = 3


@dataclass
class RankedLabels:
    """Result of ranking: the elimination order plus the per-label metadata used to build it."""

    ordered: list[Label]  # index 0 = eliminate first; tail = prefer free (masters)
    tiers: dict[Label, int]
    complexity: dict[Label, float]
    lf: dict[Label, object]  # True | False | "Unknown"

    def pivot_priority(self, label: Label) -> int:
        """Position in the elimination order (0 = eliminated first)."""
        return self.ordered.index(label)

    def free_preference(self) -> list[Label]:
        """Labels ordered from most-preferred-free (masters) to least."""
        return list(reversed(self.ordered))


def rank_labels(
    family: ParametricFamily,
    labels: Iterable[Label],
    target: Label | None = None,
    preferred_masters: Iterable[Label] = (),
    lf_map: dict[Label, object] | None = None,
    weights: dict[str, float] | None = None,
    targets: Iterable[Label] = (),
) -> RankedLabels:
    """Order ``labels`` by elimination priority. Returns all input labels (deduplicated), reordered.

    ``lf_map`` may supply precomputed local-finiteness verdicts to avoid recomputation; any label
    not present is evaluated via :func:`is_locally_finite`.

    Perf.5: ``targets`` optionally names *several* target labels; every one of them lands in
    tier 0 (eliminated before everything else), ordered among themselves by the same
    ``(-complexity, label)`` rule as within any tier. ``targets=(t,)`` is bit-identical to
    ``target=t``; ``target`` and ``targets`` may be combined (the union is used).
    """
    unique = list(dict.fromkeys(labels))  # dedup, preserve first-seen order
    preferred = set(preferred_masters)
    lf_map = dict(lf_map or {})
    target_set = set(targets)
    if target is not None:
        target_set.add(target)

    lf: dict[Label, object] = {}
    complexity: dict[Label, float] = {}
    tiers: dict[Label, int] = {}
    for label in unique:
        verdict = lf_map[label] if label in lf_map else is_locally_finite(family, label)
        lf[label] = verdict
        complexity[label] = label_complexity(label, family.nvars, family.npolys, weights)
        tiers[label] = _tier(label, verdict, target_set, preferred)

    ordered = sorted(unique, key=lambda lab: (tiers[lab], -complexity[lab], lab))
    return RankedLabels(ordered=ordered, tiers=tiers, complexity=complexity, lf=lf)


def _tier(label: Label, verdict, targets: set[Label], preferred: set[Label]) -> int:
    if label in targets:
        return TIER_TARGET
    if verdict is not True:  # False or "Unknown" -> must not be left free
        return TIER_NON_LF
    if label in preferred:
        return TIER_PREFERRED
    return TIER_LF


def ordered_labels(
    family: ParametricFamily,
    labels: Iterable[Label],
    target: Label | None = None,
    preferred_masters: Iterable[Label] = (),
    lf_map: dict[Label, object] | None = None,
    weights: dict[str, float] | None = None,
) -> list[Label]:
    """Convenience wrapper returning just the elimination-ordered label list."""
    return rank_labels(
        family, labels, target, preferred_masters, lf_map, weights
    ).ordered
