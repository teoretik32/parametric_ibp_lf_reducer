# USAGE — parametric_ibp_lf_reducer (по-русски)

*English version: [USAGE.md](USAGE.md).*

Чистый Python-редьюсер локально-конечных параметрических IBP. Wolfram/Mathematica
синтаксис встречается ТОЛЬКО как текстовый формат обмена (входные/выходные
документы); Wolfram runtime никогда не вызывается.

## Установка и тесты

```bash
python -m venv .venv
source .venv/bin/activate          # PowerShell: .venv\Scripts\Activate.ps1
python -m pip install -e '.[dev,speed]'
python -m pytest                   # полный быстрый набор
ruff check .
```

## Tiny-пример (~1–2 с, полный конвейер с сертификатом)

```bash
python -m parametric_ibp_lf_reducer reduce examples/tiny_success_input.wl.txt \
    --out result.m --diagnostics-json diagnostics.json
```

Ожидается: exit code 0, `"Status" -> "Success"`, один локально конечный терм.
Детали: `examples/tiny_success_expected_notes.md`.

Эквивалентный вызов API:

```python
from parametric_ibp_lf_reducer import api
result = api.reduce_wolfram_style_input(open("examples/tiny_success_input.wl.txt").read())
print(result.status, result.wolfram_style_text)
```

## Тяжёлый пример D4 (МЕДЛЕННО: ~10–15 минут)

`examples/d4_cli_example_input.wl.txt` — валидационная семья D4 с полной проверенной
конфигурацией, встроенной в её `"Options"`-ассоциацию (36 scattered non-lattice
sample-точек, label box, preferred masters M1..M5, off-sample rank-generic
certificate-точки). Файл самодостаточен:

```bash
python -m parametric_ibp_lf_reducer reduce examples/d4_cli_example_input.wl.txt \
    --out d4_result.m --diagnostics-json d4_diagnostics.json
```

Ожидается (зафиксированный, сертифицированный исход): exit 0,
`"Status" -> "Success"`, `"AllLocallyFinite" -> True`, `certificate_status "Passed"`
и **3-термный** LF-базис с лейблами {M1, M2, M3} =
{(0,1,1,0,-2,-1,0), (1,1,0,0,-2,-1,0), (0,1,1,0,-3,-1,0)}.
Почему это отличается от 5-термного референса — см. «Ограничения».

Соответствующий opt-in тест (то же предупреждение — ~10–15 мин на модуль):

```bash
RUN_D4_FULL=1 python -m pytest tests/test_d4_vertical.py tests/test_d4_cli_e2e.py
```

(В PowerShell сначала: `$env:RUN_D4_FULL = "1"`.)

## Exit-коды

| код | смысл |
|-----|-------|
| 0 | редукция достигла `Status -> "Success"` (пройдены LF- и certificate-гейты) |
| 1 | редукция отработала, но `Success` не достигнут — честный типизированный отказ; текст результата и JSON всё равно записываются; причина в stderr и в `"error"` |
| 2 | usage / I/O / битый документ — ничего не редуцировалось |

## Статусы

| статус | смысл |
|--------|-------|
| `Success` | все термы локально конечны и сертификат пройден |
| `Failure` (грубый экспортируемый статус; причина — в `Error` / JSON `status`) | `TargetNotReducible`, `InterpolationFailed`, `NormalFormNotLocallyFinite`, `ResourceLimitReached` |
| `VerificationFailed` | формальная редукция найдена, но независимая верификация её отвергла |
| `ParserNeedsExplicitFamily` | нет явной параметрической семьи во входе; ничего не редуцировалось |

## Поля diagnostics JSON

- `status` / `exported_status` / `success` / `error` — типизированный исход;
  `status` — либо `"Success"`, либо одна из констант `ALL_FAILURE_REASONS`.
- `certificate_status` — `"Passed"` / `"Failed"` / `"Insufficient"` / `"NotRun"`
  (row-span certificate gate включён по умолчанию; `Success` подразумевает `"Passed"`).
- `certificate` — скалярные счётчики сертификата (`n_certificate_points`,
  `n_certificate_points_passed/failed`, `selected_rank`, ...).
- `target_label`, `all_locally_finite`
- `terms[]` — `label`, `coefficient` (текст, степени через `^`), `integrand` (текст),
  `locally_finite` (True/False/"Unknown").
- `diagnostics` — `formal_success`, `reconstruction_verified`,
  `independent_validation_passed`, `n_terms`, `non_lf_terms`, `unknown_lf_terms`,
  `n_records`, `n_skipped_records`, `zero_reduction`, `messages`.

## Ограничения (Release.1)

- **Только explicit-family вход.** Документ с одним лишь `Integrand` возвращает
  типизированный `Failure/ParserNeedsExplicitFamily` — авто-факторизация integrand
  не угадывается.
- **Нет adaptive search.** Один фиксированный проход на вызов: label box / степени /
  сэмплы берутся из `Options` документа (или дефолтов); при отказе ничего не
  расширяется автоматически.
- **Нет зависимости от Mathematica runtime** — а значит, нет и символьного кросс-чека
  через Wolfram; верификация — модулярная row-span сертификация в рациональных точках.
- **Пределы плотной многомерной интерполяции.** Реконструкция коэффициентов — плотный
  перебор степеней по scattered рациональным сэмплам; высокая степень или много
  параметров требуют больше сэмплов, чем дефолты, а product-lattice сетки могут
  «подтвердить» неверные интерполянты (используйте scattered-точки, как в примерах).
- **D4 редуцируется к 3-термному LF-базису** {M1,M2,M3}, эквивалентному
  (и сертифицированному против row span, содержащего) 5-термному референсному базису
  M1..M5; референсный базис сознательно НЕ форсируется.

## Example 4* (exploratory, known-value-only)

`examples/example4_star_input.wl.txt` — пример «только известное значение»:
для него известно ε-разложение самого интеграла
(`validation/notebook_star_example4_known_value_expansion.txt`), но НЕ референсная
LF-декомпозиция. Exploratory-прогон редьюсера вернул сертифицированный `Success`
(certificate `Passed` 3/3) с 2-термной полностью локально-конечной редукцией —
артефакты: `outputs/example4_star_result.m`,
`outputs/example4_star_diagnostics.json`; детали:
`notes/example4_star_exploratory.md`. Это **не** часть certified baseline
(единственная curated end-to-end конфигурация — D4), и численный кросс-чек против
известного значения невозможен без значений master-интегралов.

## Релизная sanity-проверка

```bash
python -m pytest
ruff check .
# опциональный тяжёлый acceptance (~25–30 мин суммарно):
RUN_D4_FULL=1 python -m pytest tests/test_d4_vertical.py tests/test_d4_cli_e2e.py
```

или `scripts/final_check.ps1` / `scripts/final_check.sh` (флаг `-Heavy` / `--heavy`
добавляет D4-прогоны).

Примечание: `ruff format --check` всё ещё хочет переформатировать 22 файла,
предшествующих текущему конфигу форматирования; это pre-existing и сознательно не
смешивается с логическими изменениями. `ruff check .` — чистый.
