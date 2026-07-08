"""Reduction result / diagnostics model + the strict ``Success`` gate (spec §6, §7, §10).

This module holds the *output* contract of the reducer and the single place that is allowed to
stamp ``Status="Success"``. It deliberately does NOT orchestrate a search or run reconstruction —
it consumes already-produced reconstruction + local-finiteness evidence and decides, strictly,
whether that evidence is good enough to be a physical success.

The gate is intentionally paranoid (CLAUDE.md §5-6): a formal normal form that merely *reduces*
the target is not a success. ``Status="Success"`` requires all of:

* reconstruction/interpolation was verified;
* an independent-sample validation passed;
* every right-hand-side term is provably locally finite (``lf_flags[label] is True``);
* there are no ``"Unknown"`` local-finiteness terms and no non-locally-finite terms;
* the term list is non-empty, or an explicit valid zero-reduction is represented.

Anything short of that is a typed :class:`ReductionResult` failure that still records
``FormalSuccess`` and the formal terms for diagnostics. No ``Success`` is produced anywhere else.
"""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass, field

from .labels import Label
from .wolfram_text_export import coeff_to_wolfram_text, label_to_wolfram_text

# --- status / failure-reason constants -------------------------------------------------------
STATUS_SUCCESS = "Success"
STATUS_FAILURE = "Failure"  # coarse exported status; the concrete reason goes in ``Error``

FAILURE_TARGET_NOT_REDUCIBLE = "TargetNotReducible"
FAILURE_INTERPOLATION_FAILED = "InterpolationFailed"
FAILURE_NORMAL_FORM_NOT_LOCALLY_FINITE = "NormalFormNotLocallyFinite"
FAILURE_VERIFICATION_FAILED = "VerificationFailed"
FAILURE_RESOURCE_LIMIT_REACHED = "ResourceLimitReached"
# Input-level failure (Pass 2I.3): the text API refuses to guess an integrand factorization
# (spec §3.2) — an input without an explicit family is an honest typed Failure, never a stub run.
FAILURE_PARSER_NEEDS_EXPLICIT_FAMILY = "ParserNeedsExplicitFamily"


class FailureReason:
    """Namespace of failure statuses (kept as strings so they render straight to Wolfram text)."""

    TARGET_NOT_REDUCIBLE = FAILURE_TARGET_NOT_REDUCIBLE
    INTERPOLATION_FAILED = FAILURE_INTERPOLATION_FAILED
    NORMAL_FORM_NOT_LOCALLY_FINITE = FAILURE_NORMAL_FORM_NOT_LOCALLY_FINITE
    VERIFICATION_FAILED = FAILURE_VERIFICATION_FAILED
    RESOURCE_LIMIT_REACHED = FAILURE_RESOURCE_LIMIT_REACHED
    PARSER_NEEDS_EXPLICIT_FAMILY = FAILURE_PARSER_NEEDS_EXPLICIT_FAMILY


ALL_FAILURE_REASONS = frozenset(
    {
        FAILURE_TARGET_NOT_REDUCIBLE,
        FAILURE_INTERPOLATION_FAILED,
        FAILURE_NORMAL_FORM_NOT_LOCALLY_FINITE,
        FAILURE_VERIFICATION_FAILED,
        FAILURE_RESOURCE_LIMIT_REACHED,
        FAILURE_PARSER_NEEDS_EXPLICIT_FAMILY,
    }
)


# --- dataclasses -----------------------------------------------------------------------------
@dataclass(frozen=True)
class ReductionTerm:
    """One right-hand-side master term ``C_a(params) * Integral[F_a]`` of the reduction."""

    label: Label
    coefficient_text: str  # Wolfram-like text (``^`` not ``**``)
    integrand_text: str  # relative factor, e.g. ``x2*x3/(G0^2*G1)``
    locally_finite: object  # True | False | "Unknown"
    coefficient: object = None  # raw SymPy/ParamExpr kept for programmatic use


@dataclass
class ReductionDiagnostics:
    """Non-authoritative counters/evidence attached to a result (never decides success)."""

    formal_success: bool = False
    reconstruction_verified: bool = False
    independent_validation_passed: bool = False
    n_terms: int = 0
    non_lf_terms: tuple = ()
    unknown_lf_terms: tuple = ()
    n_records: int = 0
    n_skipped_records: int = 0
    zero_reduction: bool = False
    messages: tuple = ()
    extra: dict = field(default_factory=dict)


@dataclass
class ReductionResult:
    """Typed reduction outcome. ``status`` is ``"Success"`` or a :class:`FailureReason` value."""

    status: str
    target_label: Label
    all_locally_finite: object  # True | False | "Unknown"
    terms: tuple[ReductionTerm, ...] = ()
    formal_success: bool = False
    error: str | None = None
    diagnostics: ReductionDiagnostics = field(default_factory=ReductionDiagnostics)

    @property
    def success(self) -> bool:
        return self.status == STATUS_SUCCESS

    @property
    def exported_status(self) -> str:
        """Coarse status for the Wolfram text: ``"Success"`` or ``"Failure"``.

        The precise failure reason stays in :attr:`status` (and is exported as ``Error``).
        """
        return STATUS_SUCCESS if self.success else STATUS_FAILURE

    @property
    def failure_reason(self) -> str | None:
        """The concrete failure reason (a :class:`FailureReason` value), or ``None`` on success."""
        return None if self.success else self.status

    @property
    def wolfram_style_text(self) -> str:
        return result_to_wolfram_text(self)


# --- the strict Success gate -----------------------------------------------------------------
def build_reduction_result_from_reconstruction(
    family,
    target_label: Label,
    reconstructed_coefficients: dict | None,
    lf_flags: dict,
    *,
    reconstruction_verified: bool = False,
    independent_validation_passed: bool = False,
    formal_success: bool = False,
    interpolation_failed: bool = False,
    target_reducible: bool = True,
    verification_failed: bool = False,
    resource_limit_reached: bool = False,
    allow_zero_reduction: bool = False,
    n_records: int = 0,
    n_skipped_records: int = 0,
    messages: Iterable[str] = (),
) -> ReductionResult:
    """Assemble a :class:`ReductionResult`, applying the strict ``Success`` gate.

    ``reconstructed_coefficients`` maps ``label -> C_label`` (SymPy/ParamExpr) or is ``None`` when
    reconstruction produced nothing. ``lf_flags`` maps ``label -> True|False|"Unknown"``. The gate
    order is: resource limit -> target not reducible -> interpolation failed -> verification not
    passed -> local-finiteness -> (non-empty or explicit zero) -> ``Success``.
    """
    messages = tuple(messages)

    # Build the formal terms (if any) so failures still expose the formal normal form + LF flags.
    terms: list[ReductionTerm] = []
    non_lf: list[Label] = []
    unknown: list[Label] = []
    if reconstructed_coefficients is not None:
        for label in sorted(reconstructed_coefficients):
            coeff = reconstructed_coefficients[label]
            flag = lf_flags.get(label, "Unknown")
            terms.append(
                ReductionTerm(
                    label=label,
                    coefficient_text=coeff_to_wolfram_text(coeff),
                    integrand_text=family.label_to_wolfram_text(label),
                    locally_finite=flag,
                    coefficient=coeff,
                )
            )
            if flag is False:
                non_lf.append(label)
            elif flag is not True:
                unknown.append(label)
    terms = tuple(terms)

    if reconstructed_coefficients is None:
        all_lf: object = "Unknown"
    elif non_lf:
        all_lf = False
    elif unknown:
        all_lf = "Unknown"
    else:
        all_lf = True

    zero_reduction = reconstructed_coefficients is not None and not terms

    diag = ReductionDiagnostics(
        formal_success=formal_success,
        reconstruction_verified=reconstruction_verified,
        independent_validation_passed=independent_validation_passed,
        n_terms=len(terms),
        non_lf_terms=tuple(non_lf),
        unknown_lf_terms=tuple(unknown),
        n_records=n_records,
        n_skipped_records=n_skipped_records,
        zero_reduction=zero_reduction,
        messages=messages,
    )

    def _result(status: str, effective_lf: object, error: str | None) -> ReductionResult:
        # A failure must never advertise AllLocallyFinite=True: without a verified reduction the
        # local finiteness of the (unconfirmed) terms is not certified -> downgrade to "Unknown".
        if status != STATUS_SUCCESS and effective_lf is True:
            effective_lf = "Unknown"
        return ReductionResult(
            status=status,
            target_label=target_label,
            all_locally_finite=effective_lf,
            terms=terms,
            formal_success=formal_success,
            error=error,
            diagnostics=diag,
        )

    # --- gate (fail fast, most-fundamental problem first) ---
    if resource_limit_reached:
        return _result(
            FAILURE_RESOURCE_LIMIT_REACHED,
            all_lf,
            "resource limit reached before the reduction could be completed",
        )
    if not target_reducible:
        return _result(
            FAILURE_TARGET_NOT_REDUCIBLE,
            all_lf,
            "target integral is not reducible (not a pivot in the generated system)",
        )
    if interpolation_failed or reconstructed_coefficients is None:
        return _result(
            FAILURE_INTERPOLATION_FAILED,
            all_lf,
            "coefficient reconstruction/interpolation did not validate",
        )
    if verification_failed or not (reconstruction_verified and independent_validation_passed):
        return _result(
            FAILURE_VERIFICATION_FAILED,
            all_lf,
            "reconstruction was not verified on independent samples",
        )
    if non_lf:
        return _result(
            FAILURE_NORMAL_FORM_NOT_LOCALLY_FINITE,
            False,
            f"normal form contains non-locally-finite terms: "
            f"{[label_to_wolfram_text(x) for x in non_lf]}",
        )
    if unknown:
        return _result(
            FAILURE_NORMAL_FORM_NOT_LOCALLY_FINITE,
            "Unknown",
            f"normal form contains terms whose local finiteness is unknown: "
            f"{[label_to_wolfram_text(x) for x in unknown]}",
        )
    if not terms and not allow_zero_reduction:
        return _result(
            FAILURE_VERIFICATION_FAILED,
            all_lf,
            "empty normal form without an explicit zero-reduction",
        )

    # Every precondition met: this is the only path that stamps Success.
    return _result(STATUS_SUCCESS, True, None)


# --- Wolfram-like text export ----------------------------------------------------------------
def _wl_tribool(value: object) -> str:
    if value is True:
        return "True"
    if value is False:
        return "False"
    return "Unknown"


def _wl_label_list(labels: Iterable[Label]) -> str:
    return "{" + ",".join(label_to_wolfram_text(x) for x in labels) + "}"


def result_to_wolfram_text(result: ReductionResult) -> str:
    """Render a :class:`ReductionResult` as Wolfram-like association text (``^`` never ``**``)."""
    d = result.diagnostics
    lines: list[str] = ["<|"]
    # Coarse, unambiguous status; the precise reason (when failing) is carried by "Error".
    lines.append(f'  "Status" -> "{result.exported_status}",')
    lines.append(f'  "TargetLabel" -> {label_to_wolfram_text(result.target_label)},')
    lines.append(f'  "AllLocallyFinite" -> {_wl_tribool(result.all_locally_finite)},')

    if result.terms:
        term_lines = [
            '    <| "Integrand" -> {i}, "Coefficient" -> {c}, "LocallyFinite" -> {lf} |>'.format(
                i=t.integrand_text, c=t.coefficient_text, lf=_wl_tribool(t.locally_finite)
            )
            for t in result.terms
        ]
        lines.append('  "Terms" -> {\n' + ",\n".join(term_lines) + "\n  },")
    else:
        lines.append('  "Terms" -> {},')

    if not result.success:
        lines.append(f'  "Error" -> "{result.failure_reason}",')  # concrete FailureReason value
        if result.error:
            lines.append(f'  "ErrorDetail" -> "{result.error}",')

    diag_lines = [
        f'    "FormalSuccess" -> {_wl_tribool(bool(result.formal_success))},',
        f'    "ReconstructionVerified" -> {_wl_tribool(bool(d.reconstruction_verified))},',
        f'    "IndependentValidationPassed" -> {_wl_tribool(bool(d.independent_validation_passed))},',
        f'    "NumTerms" -> {d.n_terms},',
        f'    "NonLFTerms" -> {_wl_label_list(d.non_lf_terms)},',
        f'    "UnknownLFTerms" -> {_wl_label_list(d.unknown_lf_terms)}',
    ]
    lines.append('  "Diagnostics" -> <|\n' + "\n".join(diag_lines) + "\n  |>")
    lines.append("|>")
    return "\n".join(lines)
