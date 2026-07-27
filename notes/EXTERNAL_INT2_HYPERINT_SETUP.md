# External Int2 Method.12B — Phase 0 runtime audit (HyperInt/Maple setup report)

Date: 2026-07-27. Host: Windows 10 Pro 10.0.19044, working copy `B:\Soft\math_scratch`,
branch `feature/external-int2-hyperint-integration` (created from `main` @ `4f84c6a`).

## Verdict

**BLOCKED: Maple is not installed on this machine; the HyperInt package is not present.**
Per the Method.12B gate, no integration was attempted, no results were faked or
substituted. Phases 1-5 have not been started.

## Probes performed (all negative unless stated)

1. Executables on PATH: `maple`, `cmaple`, `maple.exe`, `cmaple.exe`, `wmaple`,
   `maple2023`, `maple2024`, `maple2025` — none found (`command -v`, PowerShell
   `Get-Command '*maple*'`).
2. Install directories: `C:\Program Files\Maple*`, `C:\Program Files (x86)\Maple*`,
   `B:\soft\Maple*`, root-level `Maple*` on drives C:, D:, E:, B: — none found.
3. Registry: `HKLM:\SOFTWARE\Maplesoft` — key absent.
4. Environment variables matching `maple|hyperint` — none set.
5. WSL: no Linux distributions installed (`wsl -l -q` returns usage text only),
   so no Linux-side Maple either.
6. HyperInt package files (`HyperInt*.mpl`) under `B:\soft`, `C:\Users\Teoretik`,
   repo tree — none found (only our four prepared *input* files exist, see below).

## What is required to unblock

1. **Maple** (command-line `cmaple`/`maple` is sufficient; GUI not needed).
   HyperInt is pure Maple code; any reasonably recent Maple (>= 2016) works.
2. **HyperInt** by E. Panzer: `HyperInt.mpl` (+ optional `periodLookups.m` tables)
   from the author's public repository (bitbucket: PanzerErik/hyperint).
   Place `HyperInt.mpl` where the inputs can `read` it (repo root or
   `validation/hyperint/`, or patch the `read` path in the inputs).
3. Smoke test (item 1.4 of Phase 0), to be run before the L1 pilot:
   `read "HyperInt.mpl": hyperInt(1/(1+x)^2, [x]);` must return `1`, and
   `fibrationBasis(hyperInt(1/((1+x)*(x+y)), [x]), [y]);` must return a
   `-ln(y)/(y-1)`-equivalent expression, both without errors.
4. Invocation pattern planned for the pilot (deterministic, logged):
   `cmaple -q scripts/hyperint/run_L1.mpl > outputs/external_int2_hyperint_L1.log`
   where the runner `read`s the certified input, uncomments nothing by hand
   (activation happens in the runner, not by editing the certified input),
   `save`s the exact result to `validation/hyperint/external_int2_lf_L1_result.m`
   and writes `..._L1_meta.json` (wall time via `time[real]()`, `kernelopts(bytesalloc)`
   peak, integration order, epsOrder, alphabet, output size, completion flag).

## Phase 0 item 3 — static audit of the four prepared inputs (PASS with noted gaps)

Files: `validation/hyperint/external_int2_lf_L{1,2,3,4}.mpl` (12 lines each).

Checks that PASS for all four:
- full master integrand is defined inline, exactly matching the certified basis
  (G0=x2+1, G1=x5+1, G2=x7+1, G3=r*x2*x5+x2*x7+x7+1):
  - L1 `[-1,0,0,-1,-1,0,0]`: `x2^ep*(x2+1)^(ep-1)*(x5+1)^(ep-1)*(x7+1)^(-ep-1)*G3^(ep-1)`
  - L2 `[-1,0,0,0,-1,0,-1]`: `x2^ep*(x2+1)^ep*(x5+1)^(ep-1)*(x7+1)^(-ep-1)*G3^(ep-2)`
  - L3 `[0,0,0,-1,0,0,-1]`: `x2^(ep+1)*(x2+1)^(ep-1)*(x5+1)^ep*(x7+1)^(-ep-1)*G3^(ep-2)`
  - L4 `[0,0,1,-1,0,0,-1]`: as L3 with extra `x7` factor;
- intended variable order `intOrder := [x2, x5, x7]` (from the certified
  linear-reducibility audit, 6/6 orders valid — fallback orders available);
- correct epsilon orders: L1 -> ep^2, L2 -> ep^3, L3 -> ep^4, L4 -> ep^4
  (`series(f, ep=0, epsOrder+1)` keeps through ep^epsOrder inclusive);
- `r` stays symbolic; no numeric substitution anywhere (regex-checked).

Gaps to close in the runner (NOT in the certified inputs) once Maple exists:
- `hyperInt`/`fibrationBasis` calls are commented out by design ("audit artifact
  only; review before running") — activation belongs to the Phase 1 runner;
- no deterministic output file is written yet (no `save`/`writeto`) — the runner
  must add `..._Lk_result.m` + `..._Lk_meta.json`;
- no logging/checkpointing yet — the runner must tee to `outputs/` (not committed).

## Cost projection

Not possible honestly before the L1 pilot. Qualitative expectation only: all
letters are linear with roots in {0, -1, r}; L1 (weight <= 3, 5 factors) should
be the cheapest; L3/L4 (through ep^4, weight <= 5) dominate.

## Oracle isolation statement

The Laurent oracle (`validation/external_int2_full_laurent_audit.json`, AnsvInt2)
was NOT read during this phase and remains comparison-only for Phase 4.
