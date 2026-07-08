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
python -m pip install -e '.[dev,speed]'
python -m pytest                     # fast suite
ruff check .
```

Requires Python >= 3.11. `speed` (numba) is optional.

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
- **No adaptive search** — one fixed pass per invocation; nothing enlarges the
  label box / degrees / samples on failure.
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
