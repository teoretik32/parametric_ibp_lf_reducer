from fractions import Fraction
from parametric_ibp_lf_reducer import ReducerConfig, parse_family_text, reduce_family_once

PRIMES = [2_147_483_647, 2_147_483_629, 2_147_483_587]

TEXT = """
IBPInput = <|
  "Variables" -> {u, v},
  "Parameters" -> {ep},
  "Regulators" -> {ep},
  "Domain" -> "PositiveOrthant",
  "Polynomials" -> <| "P0" -> 1 + u, "P1" -> 1 + v |>,
  "MonomialExponents" -> <| u -> -1 - ep, v -> ep |>,
  "PolynomialExponents" -> <| "P0" -> -1 + ep, "P1" -> -2 - ep |>,
  "TargetMultiplier" -> 1
|>;
"""

fam = parse_family_text(TEXT)
T = (0, 0, 0, 0)

# deterministic scattered samples (candidate default generator)
DENOMS = (7, 11, 13, 17, 19, 23)
def scattered(params, n=12):
    out = []
    for k in range(n):
        pt = {}
        for j, name in enumerate(params):
            d = DENOMS[j % len(DENOMS)]
            if j == 0:
                num = 3 * k + 1
            else:
                num = (11 * k + 5 * (j + 1)) % 37
            pt[name] = Fraction(2) + Fraction(num, d)
        out.append(pt)
    return out

samples = scattered(fam.parameters)
print("samples:", [str(s["ep"]) for s in samples])

for box in [((0, 0), (-1, 0)), ((-1, 0), (-1, 0)), ((0, 1), (-1, 0)), ((-1, 1), (-2, 0))]:
    cfg = ReducerConfig(primes=PRIMES, samples=samples, label_box=box, max_ibp_degree=1)
    res = reduce_family_once(fam, T, cfg)
    ex = res.diagnostics.extra
    print(f"box={box}: status={res.status} all_lf={res.all_locally_finite} "
          f"terms={[t.label for t in res.terms]} err={res.error!r}")
    print("   cert:", ex.get("certificate", {}).get("certificate_status"),
          "msgs:", res.diagnostics.messages)
