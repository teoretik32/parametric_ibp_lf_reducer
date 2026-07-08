Скопируй этот prompt в Claude Code из корня scratch-проекта после распаковки correction patch.

```text
Важная поправка к notebook examples patch.

Прочитай:
- docs/04_notebook_examples_addendum_ru.md
- docs/05_star_examples_policy_correction_ru.md
- validation/README_STAR_EXAMPLES_POLICY_RU.md
- notes/notebook_examples_audit_template.md

Исправь предыдущую интерпретацию, если она уже попала в notes/tests:

1. Starred examples `Example 3*` и `Example 4*` — это отдельный тип примеров: известен аналитический ответ/epsilon expansion для интеграла, но НЕ известно IBP-разложение по locally finite integrands.
2. Они никак не связаны с Example 1 и Example 2 как LF-decomposition examples.
3. Не считать `Example 3*` missing family for main-spec 11.4.
4. Не считать `Example 4*` missing family for main-spec 11.5.
5. Не использовать `validation/id4example2_expected_one_term.json` и `validation/id3example3_expected_basis_inclusion.json` как expected для starred examples; если эти файлы есть, считать их deprecated.
6. Main-spec 11.4 и 11.5 остаются pending/xfail до получения настоящих input families и LF decomposition expectations.
7. Starred examples можно использовать только как parser/sparse-poly fixtures в Pass 1A, особенно Example 4* для проверки monomial `x4^2`.
8. `tm1Coeff`/Gamma-prefactor из Example 4* не является reducer coefficient.

Обнови минимально:
- notes/assumptions.md;
- notes/test_strategy.md;
- любые уже созданные тесты/xfail markers, если они ошибочно связывают starred examples с 11.4/11.5.

Затем продолжай только Pass 1A:
- package skeleton;
- input_parser.py explicit family;
- coefficients.py rational expr + eval_mod_p;
- sparse_poly.py canonical sparse dict, arbitrary N, derivative, valuation, eval_mod_p;
- pytest на parser/coefficient/sparse_poly.

Не делай e2e reduction для starred examples. После изменений запусти:

python -m pytest

В отчёте явно напиши:
- что starred examples classified as known-value-only;
- какие tests активные;
- какие tests xfail/pending;
- что будет Pass 1B.
```
