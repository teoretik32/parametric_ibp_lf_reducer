"""Row generation: algebraic and coordinate-IBP relations over integer labels (spec §5.6).

A *row* is one linear relation ``sum_a coeff_a(params) * J[label_a] = 0`` between family members,
stored sparsely as ``{label: ParamExpr}``. Two kinds are produced here (tangent/syzygy rows are
a later pass):

A. **Algebraic rows** — from ``G_l = sum_a c_{l,a} x^a``::

       J(n,m) - sum_a c_{l,a} J(n+a, m-e_l) = 0

   These are exact identities between integrals (no surface term), always kept.

B. **Coordinate IBP rows** — from ``0 = integral of d/dx_i (P F_label)`` with ``P = x^p`` a
   monomial ansatz. Expanding ``d/dx_i(P F) = (p_i + e_i) x^(p-1_i) F
   + sum_l f_l * (P d_i G_l / G_l) F`` gives::

       0 = (p_i + e_i) J[n + p - 1_i, m]
           + sum_l sum_{b in supp G_l, b_i>0} f_l b_i c_{l,b} J[n + p + b - 1_i, m - e_l]

   A coordinate row is kept ONLY if ``P F_label`` is surface-free at ``x_i = 0`` and
   ``x_i = infinity`` (component-local check, spec §7.1). Anything else — including an
   ``"Unknown"`` surface verdict — is conservatively rejected and recorded with a reason.
"""

from __future__ import annotations

from collections.abc import Iterable, Iterator
from dataclasses import dataclass, field

from .coefficients import ParamExpr
from .family import ParametricFamily
from .labels import Label
from .surface import coordinate_primitive_surface_free, vector_field_surface_free
from .tangent_fields import TangentField
from .wolfram_text_export import coeff_to_wolfram_text


@dataclass
class Row:
    """One sparse linear relation ``sum coeff * J[label] = 0``."""

    kind: str
    provenance: dict
    terms: dict[Label, ParamExpr] = field(default_factory=dict)

    def add_term(self, label: Label, coeff: ParamExpr) -> None:
        self.terms[label] = self.terms[label] + coeff if label in self.terms else coeff

    def normalized(self) -> "Row":
        self.terms = {lab: c for lab, c in self.terms.items() if not c.is_zero}
        return self

    def is_trivial(self) -> bool:
        return len(self.terms) == 0

    def labels(self) -> list[Label]:
        return list(self.terms.keys())

    def dedup_key(self) -> frozenset:
        return frozenset(self.terms.items())


@dataclass
class RejectedRow:
    kind: str
    provenance: dict
    reason: str  # "surface_not_free" | "surface_unknown" | "trivial_row"


@dataclass
class RowGenerationResult:
    rows: list[Row] = field(default_factory=list)
    rejected: list[RejectedRow] = field(default_factory=list)

    def __len__(self) -> int:
        return len(self.rows)

    def __iter__(self) -> Iterator[Row]:
        return iter(self.rows)

    def all_labels(self) -> set[Label]:
        out: set[Label] = set()
        for row in self.rows:
            out.update(row.terms.keys())
        return out


# ---- label arithmetic -------------------------------------------------------
def _shift(label: Label, nvars: int, dn, dm) -> Label:
    n = [label[i] + dn[i] for i in range(nvars)]
    m = [label[nvars + j] + dm[j] for j in range(len(dm))]
    return (*n, *m)


def monomials_up_to(nvars: int, max_degree: int) -> Iterator[tuple[int, ...]]:
    """All exponent vectors in ``N^nvars`` with total degree in ``[0, max_degree]``."""

    def rec(pos: int, remaining: int) -> Iterator[tuple[int, ...]]:
        if pos == nvars - 1:
            for e in range(remaining + 1):
                yield (e,)
            return
        for e in range(remaining + 1):
            for rest in rec(pos + 1, remaining - e):
                yield (e, *rest)

    if nvars < 1:
        raise ValueError("nvars must be >= 1")
    yield from rec(0, max_degree)


# ---- algebraic rows ---------------------------------------------------------
def algebraic_row(family: ParametricFamily, label: Label, poly_index: int) -> Row:
    """``J(n,m) - sum_a c_{l,a} J(n+a, m-e_l) = 0`` for polynomial ``l = poly_index``."""
    nvars, npolys = family.nvars, family.npolys
    params = family.parameters
    gl = family.polynomials[family.poly_names[poly_index]]
    row = Row("algebraic", {"seed": label, "poly": poly_index})
    row.add_term(label, ParamExpr.one(params))
    dm = tuple(-1 if j == poly_index else 0 for j in range(npolys))
    for a, ca in gl.terms.items():
        row.add_term(_shift(label, nvars, a, dm), ca.scale_int(-1))
    return row.normalized()


def generate_algebraic_rows(
    family: ParametricFamily, seed_labels: Iterable[Label], dedup: bool = True
) -> RowGenerationResult:
    result = RowGenerationResult()
    seen: set[frozenset] = set()
    for label in seed_labels:
        for poly_index in range(family.npolys):
            row = algebraic_row(family, label, poly_index)
            if row.is_trivial():
                result.rejected.append(RejectedRow("algebraic", row.provenance, "trivial_row"))
                continue
            key = row.dedup_key()
            if dedup and key in seen:
                continue
            seen.add(key)
            result.rows.append(row)
    return result


# ---- coordinate IBP rows ----------------------------------------------------
def coordinate_ibp_primitive_row(
    family: ParametricFamily, label: Label, var_index: int, multiplier_exps
) -> Row:
    """Raw expansion of ``0 = integral d/dx_i(P F_label)`` (no surface check applied)."""
    nvars, npolys = family.nvars, family.npolys
    params = family.parameters
    p = tuple(int(e) for e in multiplier_exps)
    e, f = family.exponent_at_label(label)
    unit = tuple(1 if k == var_index else 0 for k in range(nvars))
    row = Row("coordinate_ibp", {"seed": label, "var": var_index, "P": p})

    # (p_i + e_i) J[n + p - 1_i, m]
    dn0 = tuple(p[k] - unit[k] for k in range(nvars))
    coeff0 = e[var_index] + ParamExpr.from_int(p[var_index], params)
    row.add_term(_shift(label, nvars, dn0, (0,) * npolys), coeff0)

    # sum_l sum_b f_l b_i c_{l,b} J[n + p + b - 1_i, m - e_l]
    for pidx, name in enumerate(family.poly_names):
        fl = f[pidx]
        dm = tuple(-1 if j == pidx else 0 for j in range(npolys))
        for b, cb in family.polynomials[name].terms.items():
            bi = b[var_index]
            if bi == 0:
                continue
            dn = tuple(p[k] + b[k] - unit[k] for k in range(nvars))
            row.add_term(_shift(label, nvars, dn, dm), fl * cb.scale_int(bi))
    return row.normalized()


def generate_coordinate_ibp_rows(
    family: ParametricFamily,
    seed_labels: Iterable[Label],
    max_degree: int,
    eps_direction: str = "minus",
    dedup: bool = True,
) -> RowGenerationResult:
    """Generate surface-filtered coordinate IBP rows for monomial multipliers up to ``max_degree``."""
    result = RowGenerationResult()
    seen: set[frozenset] = set()
    for label in seed_labels:
        for var_index in range(family.nvars):
            for p in monomials_up_to(family.nvars, max_degree):
                prov = {"seed": label, "var": var_index, "P": p}
                verdict = coordinate_primitive_surface_free(
                    family, label, var_index, multiplier_exps=p, eps_direction=eps_direction
                )
                if verdict is not True:
                    reason = "surface_unknown" if verdict == "Unknown" else "surface_not_free"
                    result.rejected.append(RejectedRow("coordinate_ibp", prov, reason))
                    continue
                row = coordinate_ibp_primitive_row(family, label, var_index, p)
                if row.is_trivial():
                    result.rejected.append(RejectedRow("coordinate_ibp", prov, "trivial_row"))
                    continue
                key = row.dedup_key()
                if dedup and key in seen:
                    continue
                seen.add(key)
                result.rows.append(row)
    return result


# ---- tangent (syzygy) IBP rows ---------------------------------------------
def tangent_ibp_primitive_row(
    family: ParametricFamily, label: Label, field: TangentField
) -> Row:
    """Raw expansion of ``0 = integral div(Q F_label)`` for a tangent field ``Q`` (no m-shift).

    Uses ``div(Q F) = F [ div Q + sum_i e_i Q_i / x_i + sum_l f_l H_l ]``, where ``H_l`` are the
    field's stored multipliers (``Q.grad G_l = H_l G_l``) — the divisibility is NOT recomputed
    here. All terms keep the source label's ``m`` tuple.
    """
    nvars, npolys = family.nvars, family.npolys
    e, f = family.exponent_at_label(label)
    dm0 = (0,) * npolys
    row = Row("tangent_ibp", {"seed": label, "field": field})

    # div Q = sum_i d_i Q_i
    for i in range(nvars):
        for c, cc in field.components[i].derivative(i).terms.items():
            row.add_term(_shift(label, nvars, c, dm0), cc)
    # sum_i e_i Q_i / x_i
    for i in range(nvars):
        ei = e[i]
        for c, cc in field.components[i].terms.items():
            dn = tuple(c[k] - (1 if k == i else 0) for k in range(nvars))
            row.add_term(_shift(label, nvars, dn, dm0), ei * cc)
    # sum_l f_l H_l
    for li in range(npolys):
        fl = f[li]
        for d, dc in field.multipliers[li].terms.items():
            row.add_term(_shift(label, nvars, d, dm0), fl * dc)
    return row.normalized()


def generate_tangent_ibp_rows(
    family: ParametricFamily,
    seed_labels: Iterable[Label],
    fields: Iterable[TangentField],
    eps_direction: str = "minus",
    dedup: bool = True,
) -> RowGenerationResult:
    """Generate surface-filtered tangent IBP rows (``div(Q F)``) for verified tangent fields.

    A field that is not tangent to all ``G_l`` is rejected (never turned into a row). A row is
    kept only if its toric flux is surface-free (``vector_field_surface_free`` is ``True``);
    otherwise it is rejected with a reason. No local-finiteness filtering of intermediate labels
    is applied.
    """
    result = RowGenerationResult()
    seen: set[frozenset] = set()
    seed_labels = list(seed_labels)
    for tf in fields:
        if not tf.is_tangent(family):
            result.rejected.append(
                RejectedRow("tangent_ibp", {"field": tf}, "field_not_tangent")
            )
            continue
        for label in seed_labels:
            prov = {"seed": label, "field": tf}
            verdict = vector_field_surface_free(
                family, label, list(tf.components), eps_direction=eps_direction
            )
            if verdict is not True:
                reason = "surface_unknown" if verdict == "Unknown" else "surface_not_free"
                result.rejected.append(RejectedRow("tangent_ibp", prov, reason))
                continue
            row = tangent_ibp_primitive_row(family, label, tf)
            if row.is_trivial():
                result.rejected.append(RejectedRow("tangent_ibp", prov, "trivial_row"))
                continue
            key = row.dedup_key()
            if dedup and key in seen:
                continue
            seen.add(key)
            result.rows.append(row)
    return result


# ---- rendering (debug / diagnostics) ---------------------------------------
def render_row(family: ParametricFamily, row: Row) -> str:
    """Human-readable Wolfram-like rendering of a row as ``(coeff) * J[label] + ... = 0``."""
    parts = []
    for label, coeff in sorted(row.terms.items(), key=lambda kv: kv[0]):
        parts.append(f"({coeff_to_wolfram_text(coeff)})*J[{family.label_to_wolfram_text(label)}]")
    return " + ".join(parts) + " = 0"
