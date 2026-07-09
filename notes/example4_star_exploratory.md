# Example 4* — exploratory certified LF reduction (Success)

Дата: 2026-07-08. Статус: **exploratory** — НЕ часть certified baseline (единственная
curated end-to-end конфигурация по-прежнему D4). Записано по итогам фонового прогона
(exit 0, ~45–50 мин wall на этой машине).

## Вход и команда

- Вход: `examples/example4_star_input.wl.txt` (explicit family; Options в документе:
  `TargetLabel -> {0,0,0, 0,0,0,0}`, `MaxIBPDegree -> 2`, `TangentDegrees -> {{1,1},{2,2}}`,
  `LabelBox -> {{{0,1},{0,1},{0,1}}, {{-2,0},{-2,0},{-2,0},{-2,0}}}`).
- Команда:

```bash
python -m parametric_ibp_lf_reducer reduce examples/example4_star_input.wl.txt \
    --out outputs/example4_star_result.m \
    --diagnostics-json outputs/example4_star_diagnostics.json
```

## Результат

`status = Success`, `certificate_status = Passed` (3/3, rank-filtered 0,
rank-exceeded 0, bad 0), `selected_rank = 7257`, `n_records = 36`,
`n_skipped_records = 0`, `zero_reduction = false`,
`reconstruction_verified = true`, `independent_validation_passed = true`.

Target `{0,0,0,0,0,0,0}` → 2 терма, `all_locally_finite = true`:

| label | integrand | coefficient |
|---|---|---|
| `{1,1,0,-1,0,0,0}` | `x4*x7/G0` | `-(3703*ep^3 - 521*ep^2 - 57*ep - 1)/(5500*ep^3)` |
| `{1,1,0,0,0,-1,0}` | `x4*x7/G2` | `-(ep + 1)*(17*ep + 1)*(48*ep + 1)/(5500*ep^3)` |

Артефакты: `validation/example4_star_result.m`, `validation/example4_star_diagnostics.json`.

## Ограничения интерпретации (политика docs/05)

- Example 4* — **known-value-only** (`KnownIntegralValueNotLFDecomposition`):
  в ноутбуке известно только значение интеграла
  (`validation/notebook_star_example4_known_value_expansion.txt`, ε-разложение в `ee`),
  а НЕ референсная LF-декомпозиция.
- Полный численный кросс-чек редукции против known value невозможен без значений
  master-интегралов `{1,1,0,-1,0,0,0}` и `{1,1,0,0,0,-1,0}`; ε-разложение **нельзя**
  использовать как reducer coefficient (явное предупреждение в validation-файле).
- Сравнение с deprecated артефактами
  (`validation/id3example3_expected_basis_inclusion.json` и т.п.) не выполнялось и
  выполняться не должно.
- Прогон тяжёлый (~45–50 мин) — в regression-набор не включён; повтор вручную по
  команде выше.

## Corrected Example 4* (v0.1.1)

Исправленная версия примера: `examples/example4_star_corrected_input.wl.txt`
(множитель `15*ep + 24*ep*x7`, по линейности `15*ep*J[{0,0,0,0,0,0,0}] +
24*ep*J[{0,1,0,0,0,0,0}]`, оркестрация `scripts/run_example4_star_corrected.py`).
Результат: `Success`, certificate `Passed`, `selected_rank=9924`; артефакты
`validation/example4_star_corrected_result.m` /
`validation/example4_star_corrected_diagnostics.json`; тесты
`tests/test_example4_star_corrected.py`. Caveats этого документа (known-value-only,
не baseline, без численного кросс-чека) остаются в силе.
