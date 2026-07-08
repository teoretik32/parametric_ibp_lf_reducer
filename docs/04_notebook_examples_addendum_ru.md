# Addendum: examples extracted from `Examples_for_IBP_parametric.nb`

Источник: Mathematica notebook `Examples_for_IBP_parametric.nb`.

Этот файл добавляет notebook-примеры как fixtures, но **не все notebook-примеры являются IBP-разложениями в локально конечные интегралы**. Mathematica/Wolfram здесь остается только текстовым форматом и математической нотацией; backend Mathematica не должен становиться зависимостью Python-пакета.

## Исправленная классификация

### Group A: examples with explicit IBP/LF decomposition in the notebook

1. `Example 1` / `I4exampl1`

   Notebook содержит:
   - initial integrand;
   - `JIntegrandList`;
   - `JCoefficientListt`;
   - правило вида `Sum[coeff_i*JIntegrand_i]`.

   Это можно использовать как notebook-defined known decomposition fixture. Он имеет ту же base-family, что canonical D=4 case в main spec, но basis в notebook другой: 6 terms, а не canonical 5-term target from spec. Поэтому **не заменять canonical D=4 regression 11.3**. В Pass 1A использовать только parser/coefficient fixtures; e2e позже как separate optional known-decomposition regression.

2. `Example 2` / `I3exampl2`

   Notebook содержит:
   - initial integrand;
   - `JIntegrandList`;
   - `JCoefficientListt`;
   - правило вида `Sum[coeff_i*JIntegrand_i]`.

   Это можно использовать как notebook-defined N=3 known-decomposition fixture. Не утверждать автоматически, что это именно старый regression 11.2 из main spec, пока mapping к старым обозначениям `A,B,C` не задокументирован явно. В Pass 1A использовать только parser/coefficient fixtures; e2e позже как separate optional known-decomposition regression или как 11.2 only after audit.

### Group B: starred examples with known integral value/epsilon expansion, not known LF decomposition

3. `Example 3*` / `I4exampl3`

   Notebook содержит initial integrand и аналитический ответ/epsilon expansion для интеграла. Notebook **не содержит** `JIntegrandList`/`JCoefficientListt` и **не задаёт разложение по locally finite integrands**.

   Следовательно:
   - не считать это missing family for main-spec 11.4;
   - не подставлять expected one-term LF factor/coefficient from 11.4;
   - не делать e2e reducer regression для AllLocallyFinite на основе этого примера;
   - можно использовать как parser/sparse-poly fixture и как future value-level/numerical cross-check only.

4. `Example 4*` / `tm1Int` / `I3exampl4`

   Notebook содержит initial integrand, prefactor/gamma expression и аналитический ответ/epsilon expansion для интеграла. Notebook **не содержит** разложение по locally finite integrands.

   Следовательно:
   - не считать это missing family for main-spec 11.5;
   - не использовать main-spec expected basis inclusion как expected для этого notebook example;
   - не использовать `tm1Coeff` как reducer coefficient: это integral-value prefactor/gamma object, not an IBP LF coefficient;
   - можно использовать `tm1Int` в Pass 1A как parser/sparse-poly fixture, особенно для проверки monomial `x4^2`.

## Pass 1A policy

In Pass 1A these examples must not expand the scope into reduction. Add only parser/sparse-poly/coefficient fixtures:

- parse Example 1 and Example 2 explicit-family files;
- parse starred Example 3* and 4* input-only files as explicit families if present;
- check variables/parameters/regulators are preserved;
- check polynomial support, including `x4^2` in starred Example 4 input;
- check rational coefficient parsing/eval_mod_p only for Group A coefficients.

No end-to-end `Success` tests should be added in Pass 1A. Full e2e tests remain xfail/integration until row_generation, RREF, reconstruction and reducer orchestration exist.

## Deprecated files from an older patch

If these files exist from a previous patch, treat them as deprecated/misleading:

- `examples/id4example2_candidate_explicit_family.wl.txt`
- `examples/id3example3_x4_squared_candidate_family.wl.txt`
- `validation/id4example2_expected_one_term.json`
- `validation/id3example3_expected_basis_inclusion.json`

They must not be used as proof that starred notebook examples are main-spec 11.4/11.5. Prefer the corrected names:

- `examples/notebook_star_example3_known_value_input.wl.txt`
- `examples/notebook_star_example4_known_value_input.wl.txt`
- `validation/notebook_star_example3_known_value_expansion.txt`
- `validation/notebook_star_example4_known_value_expansion.txt`
