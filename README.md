# parametric_ibp_lf_reducer

*Читайте это на русском: [README.ru.md](README.ru.md).*

Pure-Python reducer for **parametric IBP (integration-by-parts) systems** that
produces reductions whose every term is **locally finite at `epsilon = 0`**, with an
independent exact-modular **row-span certificate** gating each `Success`.

Wolfram/Mathematica syntax appears **only as a text exchange format** for
input/output documents. No Mathematica/Wolfram runtime is required or ever invoked.

## What "locally finite at epsilon = 0" means

A reduction term `C(params, ep) * Integral[F]` is *locally finite* when the
integrand `F` produces an integral whose parametric representation stays
integrable term-by-term as `ep -> 0` (no `1/ep` blow-up hidden inside a single
term). The reducer refuses to report `Success` if it only reached a *formal*
normal form whose terms are not all locally finite — that case is an honest
typed `Failure` (`NormalFormNotLocallyFinite`), never a silent success.

## Input requirement: explicit family

Input documents must carry an **explicit parametric family** (propagators /
family definition). A document with only a raw `Integrand` is *not* guessed at:
the run returns a typed `Failure` with reason `ParserNeedsExplicitFamily`.

## Install

```bash
python -m venv .venv
source .venv/bin/activate            # PowerShell: .venv\Scripts\Activate.ps1
pip install -e .                     # base install (dict RREF backend)
pip install -e ".[speed]"            # + optional Numba RREF backend
python -m pip install -e '.[dev,speed]'   # contributors: dev tools + speed
python -m pytest                     # fast suite
ruff check .
```

Requires Python >= 3.11. `speed` (numba) is optional; without it the package
runs fully on the default `dict` backend.

## CLI quick start

Tiny example (~1–2 s, full pipeline including the certificate):

```bash
python -m parametric_ibp_lf_reducer reduce examples/tiny_success_input.wl.txt \
    --out result.m --diagnostics-json diagnostics.json
```

Expected: exit code 0, `"Status" -> "Success"`, one locally finite term.

D4 heavy example (**slow, ~10–15 minutes**, the certified release baseline):

```bash
python -m parametric_ibp_lf_reducer reduce examples/d4_cli_example_input.wl.txt \
    --out d4_result.m --diagnostics-json d4_diagnostics.json
```

Expected: `Success`, `AllLocallyFinite -> True`, certificate `Passed`, and a
**3-term** LF basis {M1, M2, M3} (see Limitations).

## Python API quick start

```python
from parametric_ibp_lf_reducer import api

result = api.reduce_wolfram_style_input(
    open("examples/tiny_success_input.wl.txt").read()
)
print(result.status)              # "Success"
print(result.wolfram_style_text) # Wolfram-like result document
```

## Optional Numba RREF backend (v0.1.4)

The mod-p RREF kernel dominates heavy runs. v0.1.4 adds an optional,
auto-selectable Numba backend:

```bash
python -m parametric_ibp_lf_reducer reduce input.wl.txt --rref-backend auto
python -m parametric_ibp_lf_reducer reduce input.wl.txt --rref-backend numba_int_array_experimental
```

- **Default is still `dict`** — nothing changes unless you opt in.
- **`auto` is the recommended opt-in for large systems** when the `[speed]`
  extra is installed; it picks Numba per matrix only when the system is large
  enough (thresholds: `min_rows=500`, `min_cols=400`, `min_nnz=3000`) and
  `prime < 2^31`. Small systems normally stay on `dict`.
- `auto` **falls back to `dict`** if Numba is unavailable; an *explicit*
  Numba request instead **fails with a clear error**.
- Certified benchmark (corrected Example 4\* full box, 972 labels,
  12360 rows): wall **3963.4s (dict) → 766.5s (auto → Numba), ~5.17×**, with
  bit-identical mathematical results and certificate `Passed`. See
  [docs/PERFORMANCE.md](docs/PERFORMANCE.md) and
  [docs/NUMBA_RREF_QA.md](docs/NUMBA_RREF_QA.md).

## Statuses

| status | meaning |
|--------|---------|
| `Success` | all terms locally finite **and** the certificate gate passed |
| `Failure` (exported coarse status; concrete reason in `Error`/JSON `status`) | honest typed failure — e.g. `TargetNotReducible`, `InterpolationFailed`, `NormalFormNotLocallyFinite`, `ResourceLimitReached` |
| `VerificationFailed` | reduction was found formally, but independent verification (reconstruction check / row-span certificate) rejected it |
| `ParserNeedsExplicitFamily` | input document lacks an explicit parametric family; nothing was reduced |

Exit codes: `0` = Success, `1` = typed failure (result + JSON still written),
`2` = usage / malformed document. Diagnostics JSON fields are documented in
[docs/USAGE.md](docs/USAGE.md).

## Limitations (Release.1)

- **Explicit-family input only** — no integrand auto-factorization.
- **Adaptive search is a bounded schedule, not a prover** — opt-in `--adaptive`
  runs a deterministic escalation of fixed passes
  ([docs/ADAPTIVE_SEARCH.md](docs/ADAPTIVE_SEARCH.md)); exhausting it never
  proves that no reduction exists. Without the flag: one fixed pass per
  invocation, unchanged.
- **Dense multivariate reconstruction limits** — high coefficient degree or many
  parameters needs more scattered samples than the defaults.
- **D4 reduces to a 3-term LF basis** {M1, M2, M3}, equivalent to (and certified
  against the row span containing) the 5-term reference basis M1..M5; the
  reference basis is deliberately not forced.
- **No Mathematica runtime dependency** — hence no symbolic Wolfram cross-check;
  verification is modular row-span certification at rational points.

## Documentation

- [docs/USAGE.md](docs/USAGE.md) — full user path (English); [docs/USAGE.ru.md](docs/USAGE.ru.md) — по-русски.
- [docs/FINAL_QA.md](docs/FINAL_QA.md) — release QA record; [CHANGELOG.md](CHANGELOG.md).
- `docs/0*_ru.md` — original specification / method / validation notes (Russian).

## License

No license has been selected yet — see
[docs/LICENSE_NOTE.md](docs/LICENSE_NOTE.md) before reusing this code.
