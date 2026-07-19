"""Parametric IBP locally finite reducer.

Pure-Python package. Wolfram/Mathematica syntax is used only as a textual I/O format; there is
no Mathematica runtime dependency.

The public entry point is ``reduce_wolfram_style_input`` (text in, typed result out) /
``reduce_wolfram_style_input_to_text`` (text in, Wolfram-like text out), a thin layer over the
non-adaptive ``reduce_family_once`` pass. All ``Success`` decisions go through the strict gate
in :mod:`result`; the row-span certificate gate is default-on.
"""

from __future__ import annotations

from .adaptive import (
    AdaptiveLevelReport,
    AdaptiveSearchConfig,
    AdaptiveSearchDiagnostics,
    SearchLevel,
    default_search_levels,
    reduce_family_adaptive,
)
from .api import (
    build_reducer_config,
    default_scattered_samples,
    reduce_wolfram_style_input,
    reduce_wolfram_style_input_adaptive,
    reduce_wolfram_style_input_to_text,
)
from .certificate import (
    STATUS_IN_SPAN,
    STATUS_NOT_IN_SPAN,
    CertificateResult,
    verify_reduction_relation_mod_p,
)
from .coefficients import ParamExpr
from .family import IntegrandFactor, ParametricFamily
from .finite_field import (
    batch_inverse,
    generate_primes,
    inv_mod,
    is_probable_prime,
    powmod,
)
from .input_parser import (
    ParserError,
    ParserNeedsExplicitFamily,
    parse_explicit_family,
    parse_family_text,
    parse_mathematica_association,
    try_factor_integrand,
)
from .labels import (
    LabelIndex,
    enumerate_box,
    label_complexity,
    make_label,
    split_label,
    zero_label,
)
from .lf_feasibility import (
    LFFeasibilityResult,
    feasibility_to_payload,
    lf_reduction_coefficients_mod_p,
    lf_reduction_feasible_mod_p,
)
from .modular_normal_form import (
    BadSpecialization,
    NormalFormResult,
    assemble_rows_mod_p,
    modular_normal_form,
)
from .ranking import RankedLabels, ordered_labels, rank_labels
from .records import (
    NormalFormRecord,
    collect_normal_form_records,
    collect_normal_form_records_multi,
    record_from_result,
    summarize_records,
)
from .reconstruction import (
    RANK_POLICY_ALL,
    RANK_POLICY_MAX_RANK,
    InterpolationFailed,
    collect_value_table,
    interpolate_multivariate,
    interpolate_univariate,
    rational_reconstruction,
    reconstruct_coefficients,
    reconstruct_rational,
    select_records_for_reconstruction,
)
from .result import (
    ALL_FAILURE_REASONS,
    FAILURE_INTERPOLATION_FAILED,
    FAILURE_NORMAL_FORM_NOT_LOCALLY_FINITE,
    FAILURE_PARSER_NEEDS_EXPLICIT_FAMILY,
    FAILURE_RESOURCE_LIMIT_REACHED,
    FAILURE_TARGET_NOT_REDUCIBLE,
    FAILURE_VERIFICATION_FAILED,
    STATUS_FAILURE,
    STATUS_SUCCESS,
    FailureReason,
    ReductionDiagnostics,
    ReductionResult,
    ReductionTerm,
    build_reduction_result_from_reconstruction,
    result_to_wolfram_text,
)
from .reducer import (
    CERTIFICATE_FAILED,
    CERTIFICATE_INSUFFICIENT,
    CERTIFICATE_NOT_RUN,
    CERTIFICATE_PASSED,
    ReducerConfig,
    ReducerRunDiagnostics,
    reduce_family_once,
    reduce_rows_multi,
    reduce_rows_once,
)
from .row_generation import (
    RejectedRow,
    Row,
    RowGenerationResult,
    algebraic_row,
    coordinate_ibp_primitive_row,
    generate_algebraic_rows,
    generate_coordinate_ibp_rows,
    generate_tangent_ibp_rows,
    render_row,
    tangent_ibp_primitive_row,
)
from .sparse_poly import SparsePoly
from .sparse_rref import RREFResult, rref_mod_p
from .surface import (
    coordinate_primitive_surface_free,
    regulated_sign,
    vector_field_surface_free,
)
from .tangent_fields import TangentField, generate_tangent_fields, verify_tangent
from .valuations import (
    LocalFinitenessReport,
    Ray,
    RayVerdict,
    ShiftRecommendation,
    base_score,
    compute_candidate_rays,
    explain_local_finiteness,
    is_locally_finite,
    report_to_payload,
    valuation_poly,
)
from .wolfram_text_export import (
    coeff_to_wolfram_text,
    integrand_to_wolfram_text,
    label_to_wolfram_text,
    sympy_to_wolfram_text,
)

__all__ = [
    # foundation (Pass 1A)
    "ParamExpr",
    "SparsePoly",
    "ParametricFamily",
    "parse_mathematica_association",
    "parse_explicit_family",
    "parse_family_text",
    "try_factor_integrand",
    "ParserError",
    "ParserNeedsExplicitFamily",
    # labels + projections + export (Pass 1B)
    "make_label",
    "split_label",
    "zero_label",
    "enumerate_box",
    "LabelIndex",
    "label_complexity",
    "rank_labels",
    "ordered_labels",
    "RankedLabels",
    "IntegrandFactor",
    "coeff_to_wolfram_text",
    "integrand_to_wolfram_text",
    "label_to_wolfram_text",
    "sympy_to_wolfram_text",
    # valuations / local finiteness + surface (Pass 1C)
    "Ray",
    "compute_candidate_rays",
    "valuation_poly",
    "base_score",
    "is_locally_finite",
    # Method.1 directional LF audit (External Int2)
    "RayVerdict",
    "ShiftRecommendation",
    "LocalFinitenessReport",
    "explain_local_finiteness",
    "report_to_payload",
    "coordinate_primitive_surface_free",
    "vector_field_surface_free",
    "regulated_sign",
    # row generation (Pass 2A)
    "Row",
    "RejectedRow",
    "RowGenerationResult",
    "algebraic_row",
    "generate_algebraic_rows",
    "coordinate_ibp_primitive_row",
    "generate_coordinate_ibp_rows",
    "tangent_ibp_primitive_row",
    "generate_tangent_ibp_rows",
    "render_row",
    # tangent fields (Pass 2B)
    "TangentField",
    "generate_tangent_fields",
    "verify_tangent",
    # finite field + sparse RREF (Pass 2E)
    "inv_mod",
    "powmod",
    "batch_inverse",
    "is_probable_prime",
    "generate_primes",
    "rref_mod_p",
    "RREFResult",
    # matrix assembly + modular normal form (Pass 2F)
    "assemble_rows_mod_p",
    "modular_normal_form",
    "NormalFormResult",
    "BadSpecialization",
    # LF-constrained feasibility (Method.1, External Int2)
    "LFFeasibilityResult",
    "lf_reduction_feasible_mod_p",
    "lf_reduction_coefficients_mod_p",
    "feasibility_to_payload",
    # multi-sample record collector (Pass 2G.1)
    "NormalFormRecord",
    "collect_normal_form_records",
    "collect_normal_form_records_multi",
    "record_from_result",
    "summarize_records",
    # coefficient reconstruction (Pass 2G)
    "rational_reconstruction",
    "reconstruct_rational",
    "collect_value_table",
    "interpolate_univariate",
    "interpolate_multivariate",
    "reconstruct_coefficients",
    "InterpolationFailed",
    # rank-consistency record selection (Pass D4.3)
    "select_records_for_reconstruction",
    "RANK_POLICY_MAX_RANK",
    "RANK_POLICY_ALL",
    # modular row-span certificate (Pass D4.4)
    "verify_reduction_relation_mod_p",
    "CertificateResult",
    "STATUS_IN_SPAN",
    "STATUS_NOT_IN_SPAN",
    # reduction result / diagnostics + strict Success gate (Pass 2I.1)
    "ReductionTerm",
    "ReductionDiagnostics",
    "ReductionResult",
    "FailureReason",
    "ALL_FAILURE_REASONS",
    "STATUS_SUCCESS",
    "STATUS_FAILURE",
    "FAILURE_TARGET_NOT_REDUCIBLE",
    "FAILURE_INTERPOLATION_FAILED",
    "FAILURE_NORMAL_FORM_NOT_LOCALLY_FINITE",
    "FAILURE_VERIFICATION_FAILED",
    "FAILURE_RESOURCE_LIMIT_REACHED",
    "FAILURE_PARSER_NEEDS_EXPLICIT_FAMILY",
    "build_reduction_result_from_reconstruction",
    "result_to_wolfram_text",
    # reducer orchestration MVP (Pass 2I.2)
    "ReducerConfig",
    "ReducerRunDiagnostics",
    "reduce_family_once",
    "reduce_rows_multi",
    "reduce_rows_once",
    # reducer certificate gate (Pass D4.5)
    "CERTIFICATE_PASSED",
    "CERTIFICATE_FAILED",
    "CERTIFICATE_INSUFFICIENT",
    "CERTIFICATE_NOT_RUN",
    # text-in/text-out public API (Pass 2I.3)
    "reduce_wolfram_style_input",
    "reduce_wolfram_style_input_to_text",
    "build_reducer_config",
    "default_scattered_samples",
    # opt-in adaptive search (Pass Adaptive.1)
    "SearchLevel",
    "AdaptiveSearchConfig",
    "AdaptiveLevelReport",
    "AdaptiveSearchDiagnostics",
    "default_search_levels",
    "reduce_family_adaptive",
    "reduce_wolfram_style_input_adaptive",
]
