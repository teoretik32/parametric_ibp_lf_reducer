"""Lightweight stage timing diagnostics (Perf.0).

Pure observability: a :class:`StageTimings` dict accumulates wall-clock seconds per named
stage via the :meth:`StageTimings.stage` context manager (``time.perf_counter``). Timings
never influence control flow, statuses, certificate/LF gates, or mathematical results —
they are attached to ``diagnostics.extra["timings"]`` and exported in the CLI JSON only.

All stage keys are pre-seeded to ``0.0`` by :func:`new_stage_timings` so consumers see a
stable schema even when a stage never ran (e.g. row generation under ``reduce_rows_once``).
"""

from __future__ import annotations

import time
from contextlib import contextmanager

#: Stable schema of timed stages (seconds, accumulated).
STAGE_KEYS = (
    # row generation (reduce_family_once only; 0.0 for ready-made rows)
    "row_generation_total",
    "algebraic_rows",
    "coordinate_rows",
    "tangent_fields",
    "tangent_rows",
    # label local-finiteness flags
    "lf_flags",
    # modular record collection (records_total) and its inner per-point stages
    "records_total",
    "ranking_once",  # Perf.1: ranking hoisted out of the per-record loop (built once per run)
    "assemble_rows_mod_p",
    "ranking",  # legacy per-record ranking (0.0 when the hoisted path is used)
    "rref_mod_p",
    "extract_normal_form",
    # downstream stages
    "reconstruction",
    "certificate_total",
    "certificate_points_total",
)


class StageTimings(dict):
    """``stage name -> accumulated seconds`` map with a timing context manager."""

    @contextmanager
    def stage(self, name: str):
        """Accumulate the elapsed wall-clock time of the ``with`` body into ``self[name]``."""
        t0 = time.perf_counter()
        try:
            yield self
        finally:
            self[name] = self.get(name, 0.0) + (time.perf_counter() - t0)


def new_stage_timings() -> StageTimings:
    """A fresh accumulator with every key in :data:`STAGE_KEYS` pre-seeded to ``0.0``."""
    timings = StageTimings()
    for key in STAGE_KEYS:
        timings[key] = 0.0
    return timings
