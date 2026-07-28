# External Int2 HyperInt runtime status

Chronological record. It merges the Method.12B Phase 0 runtime audit (when no
Maple was present) with the later successful runtime installation, and states
the safety status that Method.12R/13 imposes on the External Int2 masters.

Machine-readable runtime summary: `validation/hyperint/hyperint_smoke.json`.

## Historical Method.12B Phase 0

Date: 2026-07-27. Host: Windows 10 Pro 10.0.19044, working copy
`B:\Soft\math_scratch`, branch `feature/external-int2-hyperint-integration`
(created from `main` @ `4f84c6a`).

**Verdict at that time: BLOCKED — Maple was not installed on this machine and
the HyperInt package was not present.** Per the Method.12B gate no integration
was attempted, and no result was faked or substituted.

Probes performed, all negative:

1. Executables on PATH (`maple`, `cmaple`, `maple.exe`, `cmaple.exe`, `wmaple`,
   `maple2023/2024/2025`) — none found.
2. Install directories (`C:\Program Files\Maple*`, `C:\Program Files (x86)\Maple*`,
   `B:\soft\Maple*`, root-level `Maple*` on C:, D:, E:, B:) — none found.
3. Registry `HKLM:\SOFTWARE\Maplesoft` — key absent.
4. Environment variables matching `maple|hyperint` — none set.
5. WSL — no Linux distributions installed, so no Linux-side Maple either.
6. `HyperInt*.mpl` under `B:\soft`, `C:\Users\Teoretik` and the repo tree —
   none found; only the four prepared *input* files existed.

**This remains a true historical record** of the host at that date. It was
superseded by the installation below, not retracted.

The Method.12B static audit of the four prepared `.mpl` inputs (integrand text,
`intOrder := [x2, x5, x7]`, epsilon orders L1→ep^2, L2→ep^3, L3/L4→ep^4, `r`
kept symbolic) also stands as a description of those files — but see the safety
status below: the basis they encode has since been revoked.

## Runtime installation completed later

| Item | Value |
|------|-------|
| Maple version | Maple 2025.1, X86 64 WINDOWS, Jun 12 2025, Build ID 1932578 |
| `MAPLE_HOME` | `K:\_TOOLS\Maple2025` |
| cmaple | `K:\_TOOLS\Maple2025\bin.X86_64_WINDOWS\cmaple.exe` |
| HyperInt path | `K:\_TOOLS\HyperInt` (outside the repo, **not committed**) |
| HyperInt origin | `https://bitbucket.org/PanzerErik/hyperint.git` (Erik Panzer, official) |
| HyperInt revision | `ce15b287022e698d3e3b884d5c827620d20499bd` (2023-07-31) |
| Reported version | HyperInt, version 1.0 |
| `HYPERINT_HOME` | `K:\_TOOLS\HyperInt` |

Smoke results:

- **Maple smoke = 2** — `kernelopts(version); printf("MAPLE_SMOKE=%d\n",1+1)`.
- **HyperInt smoke integral = 1** — `hyperInt(1/(1+x)^2, x=0..infinity)`,
  after `fibrationBasis`, exactly the expected value. The documented API needed
  no adaptation for this revision.

Located files: `HyperInt.mpl`, `periodLookups.m` (4.8 MB MZV/alternating-sum
lookups), `HyperTests.mpl`, documentation in `doc/` plus `Manual.mw`, `README.md`.

Loading through an absolute path — `HyperInt.mpl` guards its period table with
`if not assigned(...)`, so the table can be redirected before the read:

```maple
_hyper_autoload_periods := ["K:/_TOOLS/HyperInt/periodLookups.m"]:
read "K:/_TOOLS/HyperInt/HyperInt.mpl":
```

Environment variables are set **process-local** per command; nothing is written
to the user profile, and no license file, license-server address or activation
data is copied into the repository.

A driver probe (no integration) confirmed the runner's absolute-path loading
against the L1 input with the `hyperInt` call omitted: `epsOrder = 2`,
`intOrder = [x2, x5, x7]`, `fser` 21 terms, `hyperInt`/`fibrationBasis` defined.
That validates the runtime plumbing only — it says nothing about the validity of
the L1 integrand.

## Current safety status after Method.12R/13

- The **Method.11c four-master LF basis was revoked**. The four-term relation
  `J[Target] == C1 J[L1] + ... + C4 J[L4]` is not a valid integral identity:
  the RHS carries an `ep^-3` pole `-1/(6 r^2)` the target cannot have, and in
  the convergence chamber `ep=-3/5, r=1` target `3.9267` vs RHS `-0.2891`.
- **L4 `[0,0,1,-1,0,0,-1]` fails on the mixed infinity ray `(0,-1,-1)`**
  (`base_score == 0`; `x5, x7 -> Infinity` at fixed `x2`) — it is not locally
  finite. Method.13 then showed the defect was systemic in the surface filters
  and corrected them to the complete toric criteria.
- **Prepared L1–L4 inputs are historical/setup fixtures only.** They are not a
  runnable workload, and the Method.12A linear-reducibility audit attached to
  them describes now-revoked integrands.
- **No real master integration may run** unless a new basis artifact explicitly
  has all four of:

  | Field | Required |
  |-------|----------|
  | `Status` | `Success` |
  | `AllLocallyFinite` | `True` |
  | `IntegralIdentityStatus` | valid / not revoked |
  | `SurfaceValidationStatus` | `Passed` |

What must happen before any real master integration:

1. a NEW External Int2 LF basis derived under the corrected complete-ray LF
   gate (`newton_wall_normals` + `complete_polyhedral_rays`);
2. that basis passing surface validation (Method.13 complete toric criteria);
3. an integral-identity check of the new relation — the check Method.12R shows
   the old one fails;
4. a fresh linear-reducibility audit for the new integrands, replacing the
   Method.12A audit;
5. a basis-status artifact recording all four fields above, so the runtime gate
   opens.

The repository's current `validation/external_int2_lf_result.m` records
`Status -> "Revoked(Method.12R)"`, `AllLocallyFinite -> False`,
`IntegralIdentityStatus -> "Revoked"`, `SurfaceValidationStatus -> "Failed"`,
with `FormalRowSpanCertificate -> "Passed"` retained (the modular row-span
certificate remains true as a formal statement).

## Available setup infrastructure

- **The generic Maple/HyperInt runtime smoke remains valid** — it is
  independent of which basis is integrated, and is unaffected by Method.12R/13.
- **The runner and dry-run remain reusable.** `scripts/run_external_int2_hyperint.py`
  is stdlib-only and never imports RREF, modular-record or normal-form code.
  Resolution: cmaple via `MAPLE_CLI` → `MAPLE_HOME/bin.*/cmaple[.exe]` → `PATH`;
  HyperInt via `HYPERINT_HOME`. Deterministic paths per master `<M>`:

  | Role | Path |
  |------|------|
  | Input (historical fixture) | `validation/hyperint/external_int2_lf_<M>.mpl` |
  | Generated driver | `outputs/hyperint/<M>/driver_<M>.mpl` |
  | Result | `validation/hyperint/results/external_int2_lf_<M>_result.mpl` |
  | Meta | `validation/hyperint/results/external_int2_lf_<M>_meta.json` |
  | Log | `outputs/hyperint/<M>/<M>_run.log` |

  `outputs/` is gitignored, so logs are never committed. The driver sets
  `currentdir` to `HYPERINT_HOME` and reads the input by absolute path, so the
  input's own relative `read "HyperInt.mpl"` and the `periodLookups.m` autoload
  resolve without editing the input file. The `hyperInt`/`fibrationBasis` calls
  are appended by the driver using exactly the `fser` and `intOrder` the input
  defines; no mathematics is changed.

- **Actual External Int2 L1 execution is refused by the provenance gate.**
  Every non-dry-run invocation calls `check_basis_status()` first; a revoked,
  legacy-unchecked, missing or ambiguous artifact prints a refusal naming each
  failed field, exits non-zero (code `4`), and never invokes cmaple. Any
  revocation marker (`Revoked`, `Invalidated`, `Method.12R`) fails the gate
  outright.
- **No override is provided.** `--force` governs overwriting a completed result
  only. `--basis-status PATH` selects *which* artifact is checked, for a future
  corrected basis, but that artifact must still satisfy all four conditions.

Dry run remains the only supported mode for the L1–L4 fixtures:

```
python scripts/run_external_int2_hyperint.py --master L1 --dry-run
```
