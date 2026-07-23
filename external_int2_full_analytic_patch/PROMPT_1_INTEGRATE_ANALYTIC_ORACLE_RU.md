# Prompt 1 для Claude Code / Fable: интеграция аналитического oracle External Int2

```text
Pure symbolic mathematics / parametric Feynman-integral task.
No security/network code.

Работай в существующем репозитории parametric_ibp_lf_reducer.
Не начинай проект заново.
Не меняй reducer core, LF semantics, certificate gates или RREF.
Не запускай новые многосчасовые label-box расчёты.

В корень проекта распакован patch со следующими файлами:
- scripts/audit_external_int2_full_laurent.py
- tests/test_external_int2_full_laurent.py
- notes/EXTERNAL_INT2_FULL_LAURENT_DERIVATION.md
- validation/external_int2_full_laurent_audit.json
- validation/external_int2_full_laurent_result.m
- examples/external_int2_source_reference.wl.txt
- examples/source/ParametricInt_examples_4_ChatGPT_v2.nb

Цель:
интегрировать независимый аналитический oracle для External Int2, который воспроизводит Laurent-ряд через ep^0, но НЕ заявляет LF-basis decomposition.

1. Сначала прочитай:
- notes/HANDOFF.md
- notes/EXTERNAL_INT2_AUDIT.md
- notes/EXTERNAL_INT2_FULL_LAURENT_DERIVATION.md
- scripts/audit_external_int2_full_laurent.py
- source notebook из examples/source/

2. Проверь patch:
- exact x7 preintegration;
- ODE recurrence for s0..s3;
- assembly of coefficients ep^-4..ep^0;
- no import of parametric_ibp_lf_reducer in the audit script;
- machine-readable JSON output;
- Wolfram result syntax.

3. Сверь compact formula с AnsvInt2 из notebook независимо:
- r=s/t;
- GPL scaling G(-t,...;s) -> H_{-1,...}(r);
- G(0,...,0;z)=Log[z]^n/n!;
- compare every Laurent order ep^-4, ep^-3, ep^-2, ep^-1, ep^0;
- add an explicit comparison result to JSON and notes;
- do not silently trust the patch formulas.

4. Адаптируй test loading/imports к фактической структуре репозитория.
Не переносить example-specific formulas в reducer core.

5. Add/confirm tests:
- exact x7 identity;
- ODE residual zero through s3;
- Laurent assembly residual zero through ep^0;
- source-notebook comparison at all five orders;
- numerical fingerprints at two positive (s,t) points;
- audit script has no reducer-core imports.

6. Update:
- CHANGELOG.md Unreleased;
- notes/HANDOFF.md;
- notes/EXTERNAL_INT2_AUDIT.md with a short Method.Analytic section.

7. Run:
python scripts/audit_external_int2_full_laurent.py --numeric
python -m pytest
ruff check .

8. Commit only after green:
git add scripts/audit_external_int2_full_laurent.py \
  tests/test_external_int2_full_laurent.py \
  notes/EXTERNAL_INT2_FULL_LAURENT_DERIVATION.md \
  validation/external_int2_full_laurent_audit.json \
  validation/external_int2_full_laurent_result.m \
  examples/external_int2_source_reference.wl.txt \
  examples/source/ParametricInt_examples_4_ChatGPT_v2.nb \
  CHANGELOG.md notes/HANDOFF.md notes/EXTERNAL_INT2_AUDIT.md

git commit -m "test: add full External Int2 Laurent oracle"

Не tag/release и не push без отдельного разрешения.

Отчёт до 25 строк:
- files changed;
- tests/ruff;
- notebook comparison per Laurent order;
- numeric residuals;
- confirmation that this solves the value/series but not yet the LF-basis task;
- commit hash;
- working tree status.
```
