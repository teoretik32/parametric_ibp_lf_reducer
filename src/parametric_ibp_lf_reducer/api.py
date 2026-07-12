"""Text-in/text-out public API (Pass 2I.3, CLAUDE.md "Public API target").

``reduce_wolfram_style_input`` takes an *explicit-family* Wolfram-like document (spec §3.1),
builds a :class:`ReducerConfig` from the document's ``"Options"`` association plus Python
``overrides``, runs one fixed :func:`reduce_family_once` pass, and returns the typed
:class:`ReductionResult`; ``reduce_wolfram_style_input_to_text`` renders it as Wolfram-like
association text (``^`` never ``**``).

Honesty contract:

* an input carrying only a whole ``Integrand`` (no explicit ``Polynomials`` /
  ``MonomialExponents`` / ``PolynomialExponents``) maps to a typed
  ``Failure/ParserNeedsExplicitFamily`` result — no factorization is guessed (spec §3.2);
* no ``Success`` is stamped here — everything goes through ``reduce_family_once`` and the
  strict gate in :mod:`result`; the row-span certificate gate stays default-on;
* default samples come from a deterministic *scattered* (non-lattice) generator — a product
  grid is not an independent validation basis (assumption A30);
* unknown Python override keys are rejected (``ValueError``); unknown ``"Options"`` keys in the
  input document are ignored (the document may carry hints for other consumers).

The default entry point stays the non-adaptive single pass. Opt-in adaptive search (a
deterministic schedule of fixed passes, same gates) lives in :mod:`adaptive` and is exposed
here as :func:`reduce_wolfram_style_input_adaptive`.
"""

from __future__ import annotations

from collections.abc import Mapping
from fractions import Fraction

from .adaptive import AdaptiveSearchConfig, reduce_family_adaptive
from .family import ParametricFamily
from .input_parser import (
    ParserNeedsExplicitFamily,
    parse_explicit_family,
    parse_mathematica_association,
)
from .labels import Label, zero_label
from .reducer import ReducerConfig, reduce_family_once
from .sparse_rref import RREF_BACKEND_CHOICES
from .result import (
    FAILURE_PARSER_NEEDS_EXPLICIT_FAMILY,
    ReductionDiagnostics,
    ReductionResult,
)

# --- small, deterministic defaults (test-friendly; override for real problems) ----------------
_DEFAULT_PRIMES = (2_147_483_647, 2_147_483_629, 2_147_483_587)
_DEFAULT_N_SAMPLES = 12
_DEFAULT_LABEL_BOX = ((0, 0), (-1, 0))  # n-shifts fixed at 0, each m-shift in {-1, 0}
_SCATTER_DENOMS = (7, 11, 13, 17, 19, 23)


def default_scattered_samples(parameters, n_samples: int = _DEFAULT_N_SAMPLES) -> list[dict]:
    """Deterministic scattered rational sample points for the given parameter names.

    The first parameter walks a strictly increasing sequence ``2 + (3k+1)/7``; further
    parameters jump via ``2 + ((11k + 5(j+1)) mod 37)/d_j`` with distinct odd denominators, so
    the point set has no product-lattice structure and no low-degree curve through it
    (assumption A30 — a degenerate grid can validate a wrong interpolant). Values within each
    coordinate are pairwise distinct for ``n_samples <= 37``.
    """
    parameters = list(parameters)
    if not parameters:
        return [{}]
    points: list[dict] = []
    for k in range(n_samples):
        pt: dict = {}
        for j, name in enumerate(parameters):
            denom = _SCATTER_DENOMS[j % len(_SCATTER_DENOMS)]
            num = 3 * k + 1 if j == 0 else (11 * k + 5 * (j + 1)) % 37
            pt[name] = Fraction(2) + Fraction(num, denom)
        points.append(pt)
    return points


# --- value coercion (Options arrive as loose text; overrides as Python objects) ---------------
def _is_automatic(value) -> bool:
    return value is None or (isinstance(value, str) and value.strip() == "Automatic")


def _as_int(value) -> int:
    if isinstance(value, bool):
        raise ValueError(f"expected an integer, got boolean {value!r}")
    if isinstance(value, int):
        return value
    return int(str(value).strip())


def _as_ints(value) -> tuple[int, ...]:
    if isinstance(value, (str, bytes, Mapping)):
        raise ValueError(f"expected a sequence of integers, got {value!r}")
    return tuple(_as_int(x) for x in value)


def _as_bool(value) -> bool:
    if isinstance(value, bool):
        return value
    text = str(value).strip()
    if text == "True":
        return True
    if text == "False":
        return False
    raise ValueError(f"expected True/False, got {value!r}")


def _as_fraction(value) -> Fraction:
    if isinstance(value, (Fraction, int)):
        return Fraction(value)
    return Fraction(str(value).strip())  # exact rationals only ("3/2", "5"); no floats


def _as_label(value) -> Label:
    if isinstance(value, str):
        inner = value.strip().strip("{}")
        return tuple(int(part) for part in inner.split(",") if part.strip())
    return tuple(_as_int(x) for x in value)


def _as_labels(value) -> tuple[Label, ...]:
    if isinstance(value, (str, bytes, Mapping)):
        raise ValueError(f"expected a sequence of labels, got {value!r}")
    return tuple(_as_label(x) for x in value)


def _intify(node):
    if isinstance(node, (list, tuple)):
        return tuple(_intify(x) for x in node)
    return _as_int(node)


def _as_label_box(value):
    box = _intify(value)
    if not isinstance(box, tuple) or len(box) != 2:
        raise ValueError("label_box must be a pair (n_range, m_range)")
    return box


def _as_degree_blocks(value) -> list[tuple[int, int]]:
    blocks = _intify(value)
    if not isinstance(blocks, tuple) or any(
        not isinstance(b, tuple) or len(b) != 2 for b in blocks
    ):
        raise ValueError("tangent degree blocks must be a sequence of (deg_num, deg_den) pairs")
    return [tuple(b) for b in blocks]


def _as_samples(value) -> list[dict]:
    if isinstance(value, Mapping):
        raise ValueError("samples must be a sequence of {parameter -> rational} points")
    out: list[dict] = []
    for pt in value:
        if not isinstance(pt, Mapping):
            raise ValueError(f"each sample point must map parameter -> rational, got {pt!r}")
        out.append({str(name): _as_fraction(v) for name, v in pt.items()})
    return out


# Wolfram-style "Options" keys -> canonical snake_case setting names.
_OPTION_KEYS = {
    "TargetLabel": "target_label",
    "Labels": "labels",
    "LabelBox": "label_box",
    "MaxIBPDegree": "max_ibp_degree",
    "TangentDegreeBlocks": "tangent_degree_blocks",
    "TangentDegrees": "tangent_degree_blocks",  # alias used by the example documents
    "PreferredMasters": "preferred_masters",
    "Primes": "primes",
    "Samples": "samples",
    "CertificatePoints": "certificate_points",
    "MinValidRecords": "min_valid_records",
    "RREFBackend": "rref_backend",  # Perf.11: opt-in RREF implementation selector
}


def _as_rref_backend(value) -> str:
    """Validate an RREF backend name against :data:`sparse_rref.RREF_BACKEND_CHOICES`.

    Perf.11 added the concrete backends; Perf.12 adds ``"auto"`` (experimental heuristic
    that resolves to dict or numba per matrix — selection only, identical results).
    """
    name = str(value).strip().strip('"')
    if name not in RREF_BACKEND_CHOICES:
        raise ValueError(f"unknown RREF backend {name!r}; expected one of {RREF_BACKEND_CHOICES}")
    return name


# Canonical setting -> coercer. Everything here (minus ``target_label``) is a ReducerConfig field.
_COERCERS = {
    "target_label": _as_label,
    "labels": _as_labels,
    "label_box": _as_label_box,
    "max_ibp_degree": _as_int,
    "tangent_degree_blocks": _as_degree_blocks,
    "preferred_masters": _as_labels,
    "primes": _as_ints,
    "samples": _as_samples,
    "certificate_points": _as_samples,
    "certificate_primes": _as_ints,
    "min_valid_records": _as_int,
    "min_certificate_points": _as_int,
    "require_certificate_for_success": _as_bool,  # explicit opt-out only; default stays ON
    "eps_direction": str,
    "jobs": _as_int,  # Perf.3: worker processes for record collection (1 = serial)
    "rref_backend": _as_rref_backend,  # Perf.11: backend selection only, identical results
}

# Python-side override aliases (CLAUDE.md's preferred API spells ``tangent_degrees``).
_OVERRIDE_ALIASES = {"tangent_degrees": "tangent_degree_blocks"}


def build_reducer_config(
    family: ParametricFamily, overrides: Mapping | None = None
) -> tuple[Label, ReducerConfig]:
    """Merge defaults < document ``Options`` < ``overrides`` into ``(target_label, config)``.

    Unknown document option keys are ignored; unknown override keys raise ``ValueError``;
    ``Automatic``/``None`` values mean "use the default". The certificate gate is default-on.
    """
    settings: dict = {}

    for wl_key, node in (family.options or {}).items():
        snake = _OPTION_KEYS.get(wl_key)
        if snake is None or _is_automatic(node):
            continue  # foreign/informational option keys are not ours to reject
        settings[snake] = _COERCERS[snake](node)

    for key, value in (overrides or {}).items():
        snake = _OVERRIDE_ALIASES.get(key, key)
        if snake not in _COERCERS:
            raise ValueError(
                f"unknown override {key!r}; supported: {sorted(_COERCERS)} "
                f"(aliases: {sorted(_OVERRIDE_ALIASES)})"
            )
        if _is_automatic(value):
            continue
        settings[snake] = _COERCERS[snake](value)

    width = family.nvars + family.npolys
    target = settings.pop("target_label", None)
    if target is None:
        target = zero_label(family.nvars, family.npolys)
    elif len(target) != width:
        raise ValueError(f"target_label has length {len(target)}, expected nvars+npolys = {width}")

    settings.setdefault("primes", list(_DEFAULT_PRIMES))
    settings.setdefault("samples", default_scattered_samples(family.parameters))
    if "labels" not in settings and "label_box" not in settings:
        settings["label_box"] = _DEFAULT_LABEL_BOX

    return target, ReducerConfig(**settings)


# --- typed parser failure (never raises through the text API for this expected case) ----------
def _parser_needs_explicit_family_result(detail: str) -> ReductionResult:
    return ReductionResult(
        status=FAILURE_PARSER_NEEDS_EXPLICIT_FAMILY,
        target_label=(),
        all_locally_finite="Unknown",
        terms=(),
        formal_success=False,
        error=detail,
        diagnostics=ReductionDiagnostics(messages=(detail,)),
    )


# --- public API --------------------------------------------------------------------------------
def reduce_wolfram_style_input(
    input_text: str, overrides: Mapping | None = None, **keyword_overrides
) -> ReductionResult:
    """Parse an explicit-family document, run one reduction pass, return the typed result.

    ``overrides`` (a dict) and/or keyword arguments supply Python-side config values that win
    over the document's ``"Options"``; both accept the keys listed in ``_COERCERS`` (plus the
    ``tangent_degrees`` alias). Malformed documents raise :class:`ParserError`; a document that
    needs integrand auto-factorization returns an honest ``ParserNeedsExplicitFamily`` failure.
    """
    merged = dict(overrides or {})
    merged.update(keyword_overrides)

    raw = parse_mathematica_association(input_text)
    try:
        family = parse_explicit_family(raw)
    except ParserNeedsExplicitFamily as exc:
        return _parser_needs_explicit_family_result(str(exc))

    target, config = build_reducer_config(family, merged)
    return reduce_family_once(family, target, config)


def reduce_wolfram_style_input_to_text(
    input_text: str, overrides: Mapping | None = None, **keyword_overrides
) -> str:
    """Text-in/text-out convenience: the Wolfram-like rendering of the reduction result."""
    result = reduce_wolfram_style_input(input_text, overrides, **keyword_overrides)
    return result.wolfram_style_text


def reduce_wolfram_style_input_adaptive(
    input_text: str,
    overrides: Mapping | None = None,
    *,
    search: AdaptiveSearchConfig | None = None,
    **keyword_overrides,
) -> ReductionResult:
    """Opt-in adaptive variant of :func:`reduce_wolfram_style_input` (Pass Adaptive.1).

    Same parsing, same config merge, same strict Success gates — but instead of one fixed pass
    it runs the deterministic :func:`adaptive.reduce_family_adaptive` schedule (default:
    :func:`adaptive.default_search_levels` derived from the document/override config).
    ``search`` customizes the schedule and the resource limits. The per-level history is
    attached at ``result.diagnostics.extra["adaptive"]``.
    """
    merged = dict(overrides or {})
    merged.update(keyword_overrides)

    raw = parse_mathematica_association(input_text)
    try:
        family = parse_explicit_family(raw)
    except ParserNeedsExplicitFamily as exc:
        return _parser_needs_explicit_family_result(str(exc))

    target, config = build_reducer_config(family, merged)
    return reduce_family_adaptive(family, target, config, search=search)
