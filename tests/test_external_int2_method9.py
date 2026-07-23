"""Tests for External Int2 Method.9 (surface-policy audit) and the tangent-module export."""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pytest
import sympy as sp

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(REPO_ROOT / "src"))

from parametric_ibp_lf_reducer import surface  # noqa: E402
from parametric_ibp_lf_reducer.row_generation import (  # noqa: E402
    generate_coordinate_ibp_rows,
    monomials_up_to,
)
from parametric_ibp_lf_reducer.tangent_fields import generate_tangent_fields  # noqa: E402


def _load_script(name: str):
    if name in sys.modules:
        return sys.modules[name]
    spec = importlib.util.spec_from_file_location(name, REPO_ROOT / "scripts" / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def m9():
    return _load_script("run_external_int2_method9")


@pytest.fixture(scope="module")
def family(m9):
    return m9.T2.parse_family_text(m9.INPUT_FILE.read_text(encoding="utf-8"))


@pytest.fixture(scope="module")
def fields(m9, family):
    return generate_tangent_fields(family, m9.TANGENT_BLOCKS)


EP = sp.Symbol("ep")


class TestChamberSign:
    def test_signs_at_chamber_point(self, m9):
        sign = m9.chamber_sign_factory(sp.Rational(-3, 5))
        assert sign(EP + 1, ["ep"]) == "pos"
        assert sign(-3 * EP - 1, ["ep"]) == "pos"  # the Method.9 flip score
        assert sign(EP, ["ep"]) == "neg"
        assert sign(5 * EP + 3, ["ep"]) == "zero"

    def test_foreign_symbol_is_unknown(self, m9):
        # lowercase "unknown": chamber_sign is a drop-in for surface.regulated_sign
        sign = m9.chamber_sign_factory(sp.Rational(-3, 5))
        assert sign(EP + sp.Symbol("mu"), ["ep"]) == "unknown"

    def test_limit_vs_chamber_disagree_on_flip_score(self, m9):
        # production limit policy: value at ep=0 is -1 -> "neg"; chamber: 4/5 -> "pos"
        assert surface.regulated_sign(-3 * EP - 1, ["ep"], "minus") == "neg"
        assert m9.chamber_sign_factory(sp.Rational(-3, 5))(-3 * EP - 1, ["ep"]) == "pos"


class TestSignPolicy:
    def test_restores_original(self, m9):
        orig = surface.regulated_sign
        with m9.sign_policy(lambda *a, **k: "pos"):
            assert surface.regulated_sign is not orig
        assert surface.regulated_sign is orig

    def test_restores_on_exception(self, m9):
        orig = surface.regulated_sign
        with pytest.raises(RuntimeError):
            with m9.sign_policy(lambda *a, **k: "pos"):
                raise RuntimeError("boom")
        assert surface.regulated_sign is orig


class TestConfirmedFlip:
    SEED = (-1, 0, 1, 0, -1, -1, 0)

    def test_field2_chamber_only(self, m9, family, fields):
        tf = fields[2]
        v_limit = surface.vector_field_surface_free(
            family, self.SEED, list(tf.components), eps_direction="minus"
        )
        with m9.sign_policy(m9.chamber_sign_factory(sp.Rational(-3, 5))):
            v_chamber = surface.vector_field_surface_free(
                family, self.SEED, list(tf.components), eps_direction="minus"
            )
        assert v_limit is False
        assert v_chamber is True


class TestLimitModeMatchesGenerators:
    """classify() under the production policy must reproduce the library accept sets."""

    TINY_SEEDS = [
        (0, 0, 0, 0, 0, 0, 0),
        (-1, 0, 1, 0, -1, -1, 0),
        (1, 1, 0, -1, 0, 0, -1),
        (0, 1, 1, -1, -1, 0, 0),
    ]

    def test_coordinate_accept_set(self, m9, family):
        accepted = set()
        for label in self.TINY_SEEDS:
            for var_index in range(family.nvars):
                for p in monomials_up_to(family.nvars, m9.MAX_IBP_DEGREE):
                    v = surface.coordinate_primitive_surface_free(
                        family, label, var_index, multiplier_exps=p, eps_direction="minus"
                    )
                    if v is True:
                        accepted.add((label, var_index, p))
        gen = generate_coordinate_ibp_rows(family, self.TINY_SEEDS, m9.MAX_IBP_DEGREE, dedup=False)
        gen_accept = {(r.provenance["seed"], r.provenance["var"], r.provenance["P"]) for r in gen.rows}
        gen_accept |= {
            (rej.provenance["seed"], rej.provenance["var"], rej.provenance["P"])
            for rej in gen.rejected
            if rej.reason == "trivial_row"
        }
        assert accepted == gen_accept


class TestTangentFields:
    def test_three_fields_all_tangent(self, family, fields):
        assert len(fields) == 3
        assert all(tf.is_tangent(family) for tf in fields)

    def test_independent_sympy_residuals(self, family, fields):
        xs = [sp.Symbol(v) for v in family.variables]
        for tf in fields:
            q = [c.to_sympy(family.variables) for c in tf.components]
            for li, name in enumerate(family.poly_names):
                g = family.polynomials[name].to_sympy(family.variables)
                h = tf.multipliers[li].to_sympy(family.variables)
                resid = sp.expand(
                    sum(q[i] * sp.diff(g, xs[i]) for i in range(family.nvars)) - h * g
                )
                assert resid == 0, (name, tf.degree_block)


class TestSingularExport:
    def test_render_matches_artifact(self, family):
        exporter = _load_script("export_external_int2_tangent_module")
        rendered = exporter.render_sing(family)
        assert "ring R = (0,r),(x2,x5,x7),dp;" in rendered
        assert "module K = syz(module(A));" in rendered
        artifact = (REPO_ROOT / "scripts" / "external_int2_tangent_module.sing").read_text(
            encoding="utf-8"
        )
        assert artifact.replace("\r\n", "\n") == rendered.replace("\r\n", "\n")
