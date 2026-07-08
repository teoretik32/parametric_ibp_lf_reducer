"""Parser for the Wolfram-like/Mathematica-style textual input format.

Wolfram/Mathematica syntax here is *only a text exchange format* — there is no Mathematica
backend. This module performs a small, self-contained recursive-descent parse of the subset
used by ``IBPInput`` associations::

    Association :  <| (key -> value (, key -> value)*)? |>
    List        :  { (value (, value)*)? }
    String      :  "..."
    key         :  String | bare-symbol
    value       :  Association | List | String | raw-math-expression

Raw math expressions (polynomials, exponents, coefficients) are captured as text, then
converted to SymPy exactly once (``^`` -> ``**``), during this setup phase only. Polynomials
become :class:`SparsePoly`; exponents/coefficients become :class:`ParamExpr`.
"""

from __future__ import annotations

from dataclasses import dataclass

import sympy as sp

from .coefficients import ParamExpr
from .family import ParametricFamily
from .sparse_poly import SparsePoly


class ParserError(ValueError):
    """Raised on malformed Wolfram-like input."""


class ParserNeedsExplicitFamily(ParserError):
    """Raised when a whole integrand cannot be unambiguously factored (spec §3.2)."""


# ---------------------------------------------------------------------------
# Generic AST nodes
# ---------------------------------------------------------------------------
@dataclass
class MString:
    value: str


@dataclass
class MRaw:
    """A raw math-expression substring (not yet converted to SymPy)."""

    text: str


@dataclass
class MList:
    items: list


@dataclass
class MAssoc:
    entries: list[tuple[str, object]]

    def __post_init__(self) -> None:
        self._map = {k: v for k, v in self.entries}

    def __contains__(self, key: str) -> bool:
        return key in self._map

    def __getitem__(self, key: str):
        if key not in self._map:
            raise ParserError(f"missing required key {key!r}")
        return self._map[key]

    def get(self, key: str, default=None):
        return self._map.get(key, default)

    def keys(self) -> list[str]:
        return [k for k, _ in self.entries]


# ---------------------------------------------------------------------------
# Structural parser
# ---------------------------------------------------------------------------
_OPENERS = "([{"
_CLOSERS = ")]}"


class _Scanner:
    def __init__(self, text: str) -> None:
        self.s = text
        self.i = 0
        self.n = len(text)

    def _skip_ws(self) -> None:
        while self.i < self.n and self.s[self.i] in " \t\r\n":
            self.i += 1

    def _starts(self, tok: str) -> bool:
        return self.s.startswith(tok, self.i)

    def parse_document(self) -> MAssoc:
        """Parse ``[name =] <| ... |> [;]`` and return the top-level association."""
        self._skip_ws()
        # Optional ``identifier =`` assignment prefix.
        j = self.i
        while j < self.n and (self.s[j].isalnum() or self.s[j] in "_$"):
            j += 1
        k = j
        while k < self.n and self.s[k] in " \t\r\n":
            k += 1
        if j > self.i and k < self.n and self.s[k] == "=" and not self.s.startswith("==", k):
            self.i = k + 1
        self._skip_ws()
        if not self._starts("<|"):
            raise ParserError("expected top-level Association '<| ... |>'")
        value = self.parse_assoc()
        if not isinstance(value, MAssoc):  # pragma: no cover
            raise ParserError("top-level value is not an Association")
        return value

    def parse_value(self, stops: set[str]):
        self._skip_ws()
        if self._starts("<|"):
            return self.parse_assoc()
        if self._starts("{"):
            return self.parse_list()
        if self._starts('"'):
            return self.parse_string()
        return self.read_raw(stops)

    def parse_assoc(self) -> MAssoc:
        assert self._starts("<|")
        self.i += 2
        entries: list[tuple[str, object]] = []
        self._skip_ws()
        if self._starts("|>"):
            self.i += 2
            return MAssoc(entries)
        while True:
            key = self.parse_key()
            self._skip_ws()
            if not self._starts("->"):
                raise ParserError(f"expected '->' after key {key!r} near index {self.i}")
            self.i += 2
            value = self.parse_value({",", "|>"})
            entries.append((key, value))
            self._skip_ws()
            if self._starts("|>"):
                self.i += 2
                return MAssoc(entries)
            if self.i < self.n and self.s[self.i] == ",":
                self.i += 1
                continue
            raise ParserError(f"expected ',' or '|>' in Association near index {self.i}")

    def parse_list(self) -> MList:
        assert self._starts("{")
        self.i += 1
        items: list = []
        self._skip_ws()
        if self.i < self.n and self.s[self.i] == "}":
            self.i += 1
            return MList(items)
        while True:
            items.append(self.parse_value({",", "}"}))
            self._skip_ws()
            if self.i < self.n and self.s[self.i] == "}":
                self.i += 1
                return MList(items)
            if self.i < self.n and self.s[self.i] == ",":
                self.i += 1
                continue
            raise ParserError(f"expected ',' or '}}' in List near index {self.i}")

    def parse_string(self) -> MString:
        assert self._starts('"')
        self.i += 1
        start = self.i
        while self.i < self.n and self.s[self.i] != '"':
            if self.s[self.i] == "\\":
                self.i += 1
            self.i += 1
        if self.i >= self.n:
            raise ParserError("unterminated string literal")
        value = self.s[start : self.i]
        self.i += 1
        return MString(value)

    def parse_key(self):
        self._skip_ws()
        if self._starts('"'):
            return self.parse_string().value
        raw = self.read_raw({"->"})
        return raw.text

    def read_raw(self, stops: set[str]) -> MRaw:
        start = self.i
        depth = 0
        while self.i < self.n:
            if depth == 0:
                for st in stops:
                    if self._starts(st):
                        return MRaw(self.s[start : self.i].strip())
            if self._starts("<|"):
                depth += 1
                self.i += 2
                continue
            if self._starts("|>"):
                depth -= 1
                self.i += 2
                continue
            c = self.s[self.i]
            if c == '"':
                self.parse_string()
                continue
            if c in _OPENERS:
                depth += 1
            elif c in _CLOSERS:
                depth -= 1
            self.i += 1
        return MRaw(self.s[start:].strip())


def parse_mathematica_association(text: str) -> MAssoc:
    """Parse a Wolfram-like association string into a generic :class:`MAssoc` AST."""
    return _Scanner(text).parse_document()


# ---------------------------------------------------------------------------
# Family builder (explicit format, spec §3.1)
# ---------------------------------------------------------------------------
def _symbol_list(node) -> list[str]:
    if not isinstance(node, MList):
        raise ParserError("expected a List of symbols")
    out: list[str] = []
    for item in node.items:
        if isinstance(item, MRaw):
            out.append(item.text.strip())
        elif isinstance(item, MString):
            out.append(item.value)
        else:
            raise ParserError(f"expected a bare symbol in list, got {type(item).__name__}")
    return out


def _sympify_raw(node, local_dict):
    if isinstance(node, MRaw):
        text = node.text
    elif isinstance(node, MString):
        text = node.value
    else:
        raise ParserError(f"expected a math expression, got {type(node).__name__}")
    try:
        return sp.sympify(text.replace("^", "**"), locals=local_dict)
    except (sp.SympifyError, SyntaxError, TypeError) as exc:
        raise ParserError(f"cannot parse expression {text!r}: {exc}") from exc


def _node_to_py(node):
    """Loose conversion of Options / Assumptions nodes to plain Python."""
    if isinstance(node, MString):
        return node.value
    if isinstance(node, MRaw):
        return node.text
    if isinstance(node, MList):
        return [_node_to_py(x) for x in node.items]
    if isinstance(node, MAssoc):
        return {k: _node_to_py(v) for k, v in node.entries}
    return node


def parse_explicit_family(raw: MAssoc) -> ParametricFamily:
    """Build a :class:`ParametricFamily` from an explicit-format association (spec §3.1)."""
    if not isinstance(raw, MAssoc):
        raise ParserError("parse_explicit_family expects an MAssoc")
    if "Polynomials" not in raw or "MonomialExponents" not in raw or "PolynomialExponents" not in raw:
        raise ParserNeedsExplicitFamily(
            "explicit family requires 'Polynomials', 'MonomialExponents' and "
            "'PolynomialExponents' keys"
        )

    variables = _symbol_list(raw["Variables"])
    parameters = _symbol_list(raw["Parameters"])
    regulators = _symbol_list(raw["Regulators"]) if "Regulators" in raw else []
    if not variables:
        raise ParserError("at least one integration variable is required")

    local_dict = {name: sp.Symbol(name) for name in (*variables, *parameters, *regulators)}

    domain_node = raw.get("Domain")
    domain = domain_node.value if isinstance(domain_node, MString) else "PositiveOrthant"

    polys_node = raw["Polynomials"]
    if not isinstance(polys_node, MAssoc):
        raise ParserError("'Polynomials' must be an Association")
    poly_names: list[str] = []
    polynomials: dict[str, SparsePoly] = {}
    for name, node in polys_node.entries:
        expr = _sympify_raw(node, local_dict)
        polynomials[name] = SparsePoly.from_sympy(expr, variables, parameters)
        poly_names.append(name)

    mono_node = raw["MonomialExponents"]
    if not isinstance(mono_node, MAssoc):
        raise ParserError("'MonomialExponents' must be an Association")
    monomial_exponents: list[ParamExpr] = []
    for v in variables:
        if v not in mono_node:
            raise ParserError(f"missing MonomialExponents entry for variable {v!r}")
        monomial_exponents.append(
            ParamExpr.from_sympy(_sympify_raw(mono_node[v], local_dict), parameters)
        )

    pexp_node = raw["PolynomialExponents"]
    if not isinstance(pexp_node, MAssoc):
        raise ParserError("'PolynomialExponents' must be an Association")
    polynomial_exponents: list[ParamExpr] = []
    for name in poly_names:
        if name not in pexp_node:
            raise ParserError(f"missing PolynomialExponents entry for polynomial {name!r}")
        polynomial_exponents.append(
            ParamExpr.from_sympy(_sympify_raw(pexp_node[name], local_dict), parameters)
        )

    target_node = raw.get("TargetMultiplier")
    if target_node is None:
        target_multiplier = ParamExpr.one(parameters)
    else:
        target_multiplier = ParamExpr.from_sympy(_sympify_raw(target_node, local_dict), parameters)

    assumptions_node = raw.get("Assumptions")
    assumptions = tuple(_node_to_py(assumptions_node)) if isinstance(assumptions_node, MList) else ()

    options_node = raw.get("Options")
    options = _node_to_py(options_node) if isinstance(options_node, MAssoc) else {}

    return ParametricFamily(
        variables=tuple(variables),
        parameters=tuple(parameters),
        regulators=tuple(regulators),
        domain=domain,
        poly_names=tuple(poly_names),
        polynomials=polynomials,
        monomial_exponents=tuple(monomial_exponents),
        polynomial_exponents=tuple(polynomial_exponents),
        target_multiplier=target_multiplier,
        assumptions=tuple(str(a) for a in assumptions),
        options=options,
    )


def try_factor_integrand(raw: MAssoc) -> ParametricFamily:
    """Best-effort auto-factoring of a whole ``Integrand`` (spec §3.2).

    Not implemented in Pass 1A: conservatively refuse instead of guessing an ambiguous
    factorization. Callers should fall back to the explicit family format (§3.1).
    """
    raise ParserNeedsExplicitFamily(
        "automatic integrand factorization is not available yet; provide an explicit family "
        "with 'Polynomials', 'MonomialExponents' and 'PolynomialExponents'"
    )


def parse_family_text(text: str) -> ParametricFamily:
    """Convenience: parse a full explicit-family document string into a family."""
    return parse_explicit_family(parse_mathematica_association(text))
