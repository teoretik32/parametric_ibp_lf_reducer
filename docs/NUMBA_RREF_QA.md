# Numba RREF backend — QA record (v0.1.4)

Release theme: optional, auto-selectable Numba mod-p RREF backend
(`numba_int_array_experimental`) with certified full-pipeline validation.
This document is the QA snapshot backing the v0.1.4 release; details and
per-step history live in `notes/PERF_STATUS.md` (Perf.7–Perf.13) and
`notes/RREF_BACKEND_DESIGN.md`.

## What ships

- `parametric_ibp_lf_reducer.sparse_rref` — backend registry with three
  choices: `dict` (default, historical reference implementation),
  `numba_int_array_experimental` (opt-in native kernel), and `auto`
  (per-matrix heuristic selection).
- Selection surfaced through `ReducerConfig(rref_backend=...)`, the Python
  API, and the CLI flag `--rref-backend`.
- `numba` is an **optional** dependency via the `[speed]` extra
  (`pip install -e ".[speed]"`); base installs are unaffected.

## Behavior contract

| Requested | numba installed | numba missing |
|---|---|---|
| `dict` | dict | dict |
| `numba_int_array_experimental` | numba | **hard error** (`BackendUnavailable` / CLI `EXIT_USAGE`) — no silent substitution |
| `auto` | numba iff matrix clears gates AND `prime < 2^31`, else dict | dict (silent, by design) |
| unknown name | hard error | hard error |

Auto gates (conservative, unchanged in this release):
`min_rows=500`, `min_cols=400`, `min_nnz=3000`. The decision and reason are
recorded in the run diagnostics (`backend_selection_reason`), so every run
is auditable after the fact.

## Equivalence evidence

1. **Unit/property parity** — `tests/test_rref_numba_backend.py`: identical
   pivot columns, pivot order, inversion counts, and reduced rows vs the
   `dict` backend on random dense/sparse/degenerate systems, duplicate and
   zero rows, non-trivial rank deficiency, and edge primes; results are
   bit-identical, not approximate.
2. **Auto-selection suite** — `tests/test_rref_backend_auto.py`: gate
   thresholds, numba-absent fallback, explicit-request failure semantics,
   diagnostics recording; runs in both CI jobs (with and without numba).
3. **Tiny full-pipeline parity** — `reduce` on
   `examples/tiny_success_input.wl.txt` under all three backends produces
   byte-identical Wolfram output (also enforced in CI).
4. **Certified full-box validation (Perf.13)** — corrected Example 4\*,
   972 labels / 12360 rows / rank 9924 / 36 of 36 records valid, run
   end-to-end under `dict`, explicit numba, and `auto`:
   - identical combined result: `Success`, `AllLocallyFinite=True`,
     2 certified coefficients, certificate `Passed`;
   - all cross-backend record statuses, selected ranks, and certificate
     points identical;
   - wall time 3963.4s (dict) → 803.8s (numba, 4.93×) → 766.5s (auto,
     5.17×); `rref_mod_p` 3124.1s → 656.1s.

## Known limitations

- Numba backend requires `prime < 2^31` (int64 arithmetic headroom); larger
  primes always use `dict`.
- First call in a process pays JIT compile cost (~seconds); mitigated by the
  per-process kernel cache and irrelevant for large boxes.
- `auto` never selects numba below the gate thresholds even when it might
  win; thresholds favor predictability over peak speedup on mid-size
  matrices.
- Backend is marked `_experimental` in its name deliberately: the public
  stability promise is only that `dict` remains the reference and that
  explicit requests are honored-or-fail, never silently substituted.
