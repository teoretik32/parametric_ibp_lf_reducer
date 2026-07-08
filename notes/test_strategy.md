# Test Strategy — `parametric_ibp_lf_reducer`

Подход: **test-first по слоям**. Каждый pass добавляет unit-тесты до перехода к следующему;
`python -m pytest` держим зелёным. Тяжёлые прогоны — под маркером `integration`, дублируются
лёгкими unit-тестами.

## 0. Общие правила

- Никакой тест не проходит за счёт hardcode имён `x4/x7/H/A/B/C` или `N=3/4` в ядре.
  Тесты параметризуются по `Variables`/`Parameters` из входа.
- Политика статусов: если слой ещё не умеет — честный `Failure`/partial-статус проверяется
  тестом, а фейковый `Success` при `AllLocallyFinite=False` считается провалом (есть
  анти-тест на это).
- Коэффициенты проверяются **и** символьно (SymPy `simplify(expected-got)==0`), **и**
  независимыми finite-field/rational samples. «По одной точке» — запрещено (анти-паттерн).
- LF проверяется строго при `eps=0`.
- Маркеры: `unit` (по умолчанию, быстрые), `integration` (тяжёлые: D=4 end-to-end),
  `xfail`/`skip` (pending input families).

## 1. Матрица тестов по модулям

| Модуль | Unit-тесты |
|---|---|
| `input_parser` | Парс `<|...|>`, `{}`, `->`, `^ * + - /`; explicit family из `examples/d4`; список `Variables/Parameters/Regulators`; ошибочный/неоднозначный ввод → `ParserNeedsExplicitFamily`; идемпотентность имён |
| `coefficients` | `from_sympy` round-trip; `eval_mod_p` совпадает с прямым вычислением рационального выражения mod p; деление на 0 в знаменателе mod p → корректный skip-сигнал |
| `sparse_poly` | `add/mul/pow_small` на произвольном `N`; мономы степени >1 (`x4^2`, `x^3`); `derivative` (в т.ч. по переменной, которой нет → 0); `valuation(ray)` = min по support; `eval_mod_p` |
| `labels` | `enumerate_box` размеры; `label_id` биекция; `complexity` монотонна по n-shift/m-depth/support |
| `family` | `label_to_factor` для сдвигов n и m; `exponent_at_label`; `label_to_wolfram` синтаксис; `specialize(sample,prime)` |
| `valuations` | `compute_candidate_rays` включает `±e_i` и фасетные нормали; `base_score`; `is_locally_finite` True/False/Unknown (см. §2) |
| `surface` | coordinate IBP проверяет только `x_i=0,∞` (не отбрасывает валидные строки); vector IBP — toric flux; кеш по `(field,ray)` |
| `tangent_fields` | На `G=1+x+y` degree-блок находит поле, эквивалентное `(xy,−xy)`; отбрасывание нулевых/эквивалентных полей |
| `row_generation` | algebraic row из `G_l=Σ c x^a`; coordinate IBP раскрытие; tangent IBP без m-shift; dedup by hash; surface-filter реально отсекает |
| `ranking` | non-LF выдавливаются раньше LF; target исключаем; простые LF остаются свободными |
| `finite_field` | `inv*a==1`; `batch_inverse`; `powmod` vs Python `pow` |
| `sparse_rref` | RREF vs плотный эталон (numpy mod p); normal form target; устойчивость к порядку строк; переиспользование структуры |
| `interpolation` | Восстановление известных rational functions (uni + multi); union-support; degree search |
| `reconstruction` | CRT + rational reconstruction точных рациональных; независимая валидация ловит overfit |
| `wolfram_text_export` | `3/7` не `0.428…`; `^` не `**`; factorized форма; round-trip парс→export→парс |
| `reducer`/API/CLI | Статусы; exit code 0 только при Success; `result.m` валиден; `diagnostics.json` заполнен |

## 2. Критические LF-тесты (анти-ошибка «eps как LF»)

- Кейс, где `base_score(rho,eps=0)=0`, но `d(base_score)/d(eps)<0` (сходится только при
  `Re eps<0`): `is_locally_finite → False`. Это прямой тест на Ошибку 2 из method review.
- Кейс с внутренним нулём/неполнотой набора лучей: random safety net находит дивергенцию →
  `"Unknown"`, никогда `True`; `Success` при таком term невозможен.
- Положительный контроль: явно LF интегранд (все `base_score>0`) → `True`.

## 2b. Notebook example fixtures (addendum `docs/04` + поправка 2026-07-06)

**Два типа notebook examples** (ИСПРАВЛЕНО):
- **known-decomposition** (Example 1/2): есть `JIntegrandList`/`JCoefficientListt` → явное
  LF-разложение известно;
- **known-value-only** (Example 3*/4*): известен только аналитический ответ / ε-expansion,
  LF-разложения **нет**; **не** являются families для main-spec 11.4/11.5.

В **Pass 1A** все четыре используются **только** как parser / sparse-poly / coefficient
фикстуры — никакого e2e/reduction:

| Файл | Тип | Роль в Pass 1A | e2e-статус |
|---|---|---|---|
| `examples/notebook_example1_d4_alt_explicit_family.wl.txt` | known-decomp | parser/coeff fixture; **не заменяет** canonical 5-term D=4 (alt 6-term) | позже: optional preferred-master |
| `examples/notebook_example2_n3_five_term_explicit_family.wl.txt` | known-decomp | parser/coeff fixture; tentative 11.2 (решать позже) | `pending`/`xfail` |
| `examples/id4example2_candidate_explicit_family.wl.txt` | known-value-only (Ex.3*) | parser/sparse-poly fixture; **не** family для 11.4 | не привязан к 11.4 |
| `examples/id3example3_x4_squared_candidate_family.wl.txt` | known-value-only (Ex.4*) | parser/sparse-poly fixture; **основная `x4^2`** проверка; **не** family для 11.5 | не привязан к 11.5 |

**Deprecated/misleading (от старого patch, НЕ использовать как expected output):**
`validation/id4example2_expected_one_term.json`, `validation/id3example3_expected_basis_inclusion.json`.
Имена файлов starred-примеров (`id4example2_candidate…`, `id3example3_…candidate…`) —
исторические и вводят в заблуждение; связи с 11.4/11.5 нет.

Coefficient-файлы **known-decomposition** (`expected_d4_coefficients.json`,
`notebook_example1/2_*_expected.json`) в Pass 1A проверяются только на рациональный парсинг +
`eval_mod_p` (кросс-сверка с sympy), не как reduction output.

**Main-spec 11.4 и 11.5** остаются `pending`/`xfail` до отдельной поставки настоящих input
families + LF-decomposition ожиданий (starred examples их не закрывают).

## 3. Regression cases (spec §11)

| Кейс | Статус в MVP | Что проверяем |
|---|---|---|
| 11.1 N=2 tangent | **активен** | `Q=(xy,−xy)` тангенциально к `G=1+x+y`; сгенерированная строка без m-shift; коэффициенты n-сдвигов `(nu1+n1)`, `(nu2+n2)` |
| 11.3 D=4 | **активен (integration)** | LF-базис `M1..M5`; `C1..C5` символьно + независимые modular samples + точка `ep=-3/4,r=7/5`; `AllLocallyFinite=True` |
| 11.2 N=3 five-master | `xfail` (pending family) | ожидаемый базис `1/A,1/(AB),1/(AC),1/(BC),1/(BC^2)` — активировать при наличии A,B,C |
| 11.4 ID4example2 | `xfail` (pending family) | one-term `M=x5x7/(G1^3 G2)`, `C(ep)` — coordinate-surface не слишком строгий |
| 11.5 ID3example3 | `xfail` (pending family) | `x4^2` в полиноме; базис с `x4x7/((1+x4)(1+x7))` и т.д. |

Примечание: 11.2/11.4/11.5 требуют входных families, которых нет в docs. Как только они
предоставлены — снимаем `xfail`.

## 4. Данные и фикстуры

- `examples/d4_explicit_family.wl.txt` — вход для 11.3.
- `validation/expected_d4_coefficients.json` — ожидаемые `C1..C5`, LF-базис, точка проверки.
- Мелкие семейства для unit-тестов генерируются в коде тестов (не хардкодятся в ядре).

## 5. Гейты качества

- `python -m pytest -q` зелёный после каждого pass.
- `ruff` без ошибок.
- Ни одного `Success` при `AllLocallyFinite=False` (анти-тест).
- End-to-end 11.3 в pass 5 — обязательный гейт MVP.
