# HANDOFF — `parametric_ibp_lf_reducer`

Живой handoff-документ фактического состояния (обновляется в конце каждого pass). Полный
инженерный контекст — в `notes/assumptions.md` (A1–A26). План — `notes/implementation_plan.md`.

## Текущий статус (2026-07-11)

- **Последний завершённый pass: Adaptive.1 — opt-in adaptive search** (ветка
  `feature/adaptive-search-mvp`): `adaptive.py` — детерминированное расписание эскалации
  обычных fixed-проходов (`SearchLevel`: label box m-deepening / `max_ibp_degree` /
  tangent blocks / `extra_samples` / `extra_primes`), стоп на первом *сертифицированном*
  `Success`, иначе — детерминированный best-partial отказ + полная история и
  рекомендации в `diagnostics.extra["adaptive"]`. Лимиты: `max_labels` (pre-flight),
  `max_rows` (post-level), `timeout_sec` (только между уровнями) → типизированный
  `ResourceLimitReached`, никогда не фабрикуют успех. API:
  `reduce_family_adaptive` / `reduce_wolfram_style_input_adaptive` /
  `AdaptiveSearchConfig` / `default_search_levels`; CLI: `--adaptive`,
  `--adaptive-max-levels` (без флага путь байт-в-байт прежний). Каждый уровень идёт
  через прежний строгий gate (certificate + reconstruction verification + LF);
  исчерпание расписания ≠ доказательство нередуцируемости. Тесты —
  `tests/test_adaptive_search.py` (loop-policy на стабах, math-тесты на tiny-семье
  без моков); docs — `docs/ADAPTIVE_SEARCH.md` / `.ru.md`.
  **Adaptive.1a (доводка, merge-readiness):** opt-in маска `expand_n` в
  `default_search_levels` (отмеченные n-оси расширяются симметрично на дельту уровня;
  **требует** build-time guard `max_labels` — каждый планируемый уровень обязан
  укладываться, иначе `ValueError`; guard отличается от runtime-предпроверки, которая
  пропускает уровень); в per-level отчёт добавлено поле `error` (детерминированная
  строка ≤500 символов из diagnostic messages попытки, `None` при успехе; полные
  неуспешные результаты намеренно не сохраняются); в docs/докстринге явно зафиксировано,
  что ни один лимит не жёстко-превентивен (уровни атомарны).
- **Предыдущий pass:** **Pass Verify.1 (= D4.5) — certificate gate в редьюсере**
  (поверх D4.4; см. `notes/D4_STATUS.md` §D4.5, assumptions **A29–A31**). `Success` теперь по
  умолчанию требует row-span сертификации реконструированного соотношения в независимых
  off-sample точках (`require_certificate_for_success=True`); провал/неинформативность →
  `Error="VerificationFailed"`. Verify.1-доводка: `certificate_primes`,
  `certificate_rank_policy="selected_rank"` (только она), **rank > selected_rank в
  certificate-точке = жёсткий провал** (`n_certificate_rank_exceeded`: selected rank был не
  generic). Regression product-grid false-success + corrupt-coefficient + rank-exceeded +
  mixed-points — в `tests/test_certificate_gate.py`.
- **Tests:** **203 passed, 3 skipped** (default; skipped = heavy `RUN_D4_FULL`-трио),
  `ruff check .` — clean. Heavy opt-in при `RUN_D4_FULL=1`: **13 passed** (~16 мин, один общий
  full-config прогон через module fixture; D4-Success теперь certified by construction —
  explicit rank-generic `certificate_points` в конфиге фикстуры).
- **D4 итог: сертифицированный `Success`.** Full-config (35 scattered rational сэмплов + 1
  планово rank-deficient `(2,3)`, 3 primes, deg2, tangent, preferred_masters=M1..M5):
  terms = **`{M1,M2,M3}`**, все LF; **reducer-output относение row-span-сертифицировано в 3
  off-sample rank-generic точках**, и там же сертифицирован reference C1..C5 ⟹ **эквивалентность
  modulo generated rows доказана**. M4/M5 диагностически редуцируются к {M1,M2,M3} (почему
  базис 3-term). Forced/protected-basis **не потребовался**.
- **Важный урок D4.4 (A30):** сертификат ПОЙМАЛ реальный дефект D4.3-прогона — на 6×6 integer
  lattice интерполяция прошла on-lattice holdout, но давала **неверные коэффициенты вне решётки**
  (`Π(ep−k)` deg 6 = max_deg исчезает на всей решётке). Фикс: scattered non-lattice sample grid
  (test-config; core не менялся). Отсюда правило: holdout на решётке ≠ независимая валидация;
  приёмка — row-span certificate в off-sample rank-generic точках
  (`certificate.verify_reduction_relation_mod_p`, без Success-стампа).
- **D4.3 механика:** `select_records_for_reconstruction(records, rank_policy="max_rank")` в
  `reconstruction.py`; только `Reduced`+formal_success records с максимальным наблюдённым RREF
  rank кормят reconstruction (missing label в max-rank record = честный 0; low-rank records
  нулей не вносят, skip+count). `rank_policy="all"` = старое поведение (tests/debug).
  `reducer` проверяет `min_valid_records` по post-filter count и экспортирует
  `record_selection`-диагностику в `result.diagnostics.extra` + message.
- **Success pipeline / текущие условия Success:** единственное место — строгий gate
  `result.build_reduction_result_from_reconstruction`; для reducer-путей `Success` требует ВСЕ из:
  target reducible; reconstruction (по max-rank records, A29) валидирован на holdout;
  **row-span certificate `Passed` в off-sample rank-generic точках (A31, default-on)**;
  все термы `LocallyFinite=True` (без False/Unknown); terms непусты или явный zero-reduction.
  `reduce_wolfram_style_input`/CLI-success/adaptive — ещё нет (Pass 2I.3+).

## Result / strict Success gate (Pass 2I.1)

### result.py — новый output-слой
- Dataclasses: `ReductionTerm` (label, coefficient_text, integrand_text, locally_finite, raw coeff),
  `ReductionDiagnostics` (formal_success / reconstruction_verified / independent_validation_passed /
  non_lf_terms / unknown_lf_terms / counts / zero_reduction / messages),
  `ReductionResult` (status, target_label, all_locally_finite, terms, formal_success, error, diagnostics;
  `.success` и `.wolfram_style_text`).
- Status/failure constants: `STATUS_SUCCESS` + `FailureReason` (`TargetNotReducible`,
  `InterpolationFailed`, `NormalFormNotLocallyFinite`, `VerificationFailed`, `ResourceLimitReached`).
- `build_reduction_result_from_reconstruction(family, target, coeffs|None, lf_flags, *, verified,
  validated, formal_success, interpolation_failed, target_reducible, verification_failed,
  resource_limit_reached, allow_zero_reduction, …)` — **единственное** место, стампящее `Success`.
  Gate (fail-fast): resource → target-not-reducible → interpolation-failed → not-verified →
  non-LF → Unknown-LF → (empty без zero-reduction) → **Success**. `Success` требует: reconstruction
  verified **и** independent validation passed **и** все RHS labels `lf_flags[label] is True`
  (никаких `False`/`"Unknown"`), и terms непусты **или** явный zero-reduction.
- `result_to_wolfram_text` — Wolfram-like association (`^`, не `**`). **Contract (Pass 2I.1a):**
  `Status` = coarse `"Success"`/`"Failure"`; при failure — `Error -> "<FailureReason>"`
  (`NormalFormNotLocallyFinite`/`TargetNotReducible`/`InterpolationFailed`/`VerificationFailed`/
  `ResourceLimitReached`) + `ErrorDetail` (человекочитаемо); `AllLocallyFinite` для failure
  **никогда** не `True` (coerced в `"Unknown"`); `FormalSuccess` всегда в `Diagnostics`.
  Внутренний `ReductionResult.status` хранит точный reason; `.exported_status`/`.failure_reason` —
  для маппинга. **formal_success ≠ success** (formal normal form с non-LF термами →
  `FormalSuccess=True`, `Status="Failure"`, `Error="NormalFormNotLocallyFinite"`).
- Экспортировано из пакета; `reduce_wolfram_style_input`/CLI success **не** добавлены (Pass 2I.2+).

## Что реально работает (реализовано и покрыто тестами)

| Слой | Модуль | Статус |
|------|--------|--------|
| parser (explicit family) | `input_parser.py` | ✅ |
| rational coeffs + `eval_mod_p` | `coefficients.py` (`ParamExpr`) | ✅ |
| sparse polynomials (arbitrary N) | `sparse_poly.py` | ✅ |
| family container / projections / labels | `family.py`, `labels.py` | ✅ |
| valuations / local-finiteness (True/False/Unknown) | `valuations.py` | ✅ |
| surface (coord + toric-flux, conservative) | `surface.py` | ✅ |
| row generation (algebraic + coord IBP + tangent IBP) | `row_generation.py` | ✅ |
| tangent fields (SymPy syzygy setup) | `tangent_fields.py` | ✅ |
| column ranking (pivot order) | `ranking.py` | ✅ |
| finite field + sparse RREF (GF(p)) | `finite_field.py`, `sparse_rref.py` | ✅ |
| single-sample modular normal form | `modular_normal_form.py` | ✅ |
| **multi-sample record collector** | `records.py` | ✅ |
| coeff reconstruction (**uni + multivariate**) | `reconstruction.py` | ✅ |
| **rank-consistency record selection (max-rank)** | `reconstruction.py` | ✅ (Pass D4.3) |
| **modular row-span certificate (relation verifier)** | `certificate.py` | ✅ (Pass D4.4) |
| **certificate gate в редьюсере (default-on)** | `reducer.py` | ✅ (Pass D4.5) |
| Wolfram-like text export | `wolfram_text_export.py` | ✅ (coeff/label/integrand) |
| **result/diagnostics + strict Success gate** | `result.py` | ✅ (Pass 2I.1/2I.1a) |
| **reducer orchestration (single fixed pass)** | `reducer.py` | ✅ (Pass 2I.2) |

## Записи и reconstruction (Pass 2G.1 + 2G)

### records.py (Pass 2G.1) — collector есть
- `NormalFormRecord` — плоская сериализуемая точка: `prime, sample, target_label, status,
  formal_success, coeffs: dict[label,int mod p], support: tuple[label], all_terms_lf,
  non_lf_terms, unknown_lf_terms, rank, diagnostics`.
- `collect_normal_form_records(family, rows, target, primes, samples, ...)` — оркестратор:
  прогоняет `modular_normal_form` по `samples × primes` (samples внешний цикл, детерминированный
  порядок), **каждую** точку записывает честно: BadSpecialization/TargetNotReducible/EmptySystem
  сохраняют реальный `status` (coeffs пустой), **не выбрасываются**. Reconstruction решает, что
  потреблять.
- `record_from_result`, `summarize_records` (диагностика по статусам).

### reconstruction.py (Pass 2G + 2H) — **uni + multivariate**
- ✅ `rational_reconstruction`, `reconstruct_rational` (multi-prime CRT + стабильность).
- ✅ `collect_value_table(results_or_records)` — **consumer**: принимает и `NormalFormResult`
  (`.terms`), и `NormalFormRecord` (`.coeffs`) через `_record_coeffs`. Union support, 0-fill,
  non-Reduced → skip+count (`n_skipped`).
- ✅ `interpolate_univariate` — degree search + holdout-валидация.
- ✅ `interpolate_multivariate(values, params, max_deg=6)` (**Pass 2H**) — dense linear-algebra
  ансатц `N/D`: моном-базисы степеней `num_deg/den_deg`, нуль-пространство рациональной матрицы
  (SymPy, exact `Fraction`), **degree search от простого к сложному**, требует **1-мерного**
  nullspace + **независимую holdout-валидацию** (иначе следующая степень / `InterpolationFailed`).
- ✅ `reconstruct_coefficients(records, params)` — диспетчер: 1 параметр → univariate,
  ≥2 → multivariate. Потребляет `NormalFormRecord` из collector.

Явно **отсутствует / неполно** (фиксируем, не откатываем):

- ⚠️ Multivariate — **dense**, не sparse/Zippel: число неизвестных ~ `binom(k+d, d)` по num и den,
  растёт экспоненциально по числу параметров/степени; практично для малых `k` (2–3) и умеренных
  степеней. Нужно достаточно sample-точек (`≥ #num_mon + #den_mon − 1` для fit + holdout).
- ❌ **Reducer Success pipeline / `AllLocallyFinite` gate / D4 e2e** — впереди (нужны orchestration
  + LF-gate + export; reconstruction-сторона для D4 `(ep,r)` уже готова).

## Активные stubs / honest failures (Success нигде)

- `family.specialize(sample, prime)` → `NotImplementedError` (частичный кусок —
  `specialize_polynomials`, только `G_l` mod p).
- `input_parser.try_factor_integrand` → `ParserNeedsExplicitFamily` (conservative refuse).
- `__main__.main()` (CLI) → `SystemExit` exit 1, TODO (CLI-success — Pass 2I.2+).
- `reduce_wolfram_style_input` — **не существует / не экспортирован** (Pass 2I.3). Оркестратор уже
  есть (`reducer.reduce_family_once`/`reduce_rows_once`) и вызывает strict gate, но text-in/text-out
  обёртки и CLI-success ещё нет.
- `family.specialize(...)` по-прежнему `NotImplementedError` (reducer использует
  `modular_normal_form`/`collect_normal_form_records`, а не `specialize`).
- Validation-кейсы 11.2 / 11.4 / 11.5 — `xfail`/`pending` (нет входных families); см. A7, A14, A15.

## Starred examples (важная поправка — не терять)

- **Example 1 / Example 2** — known-decomposition fixtures (есть `JIntegrandList`/`JCoefficientListt`);
  только parser/coeff фикстуры. Example 1 = alt 6-term, **не заменяет** canonical D=4 5-term.
- **Example 3\* / Example 4\*** — **known-value-only** (известен лишь аналитический ответ / eps
  expansion; LF-разложения нет). **Не** являются missing families для 11.4 / 11.5. Файлы
  `id4example2_expected_one_term.json` и `id3example3_expected_basis_inclusion.json` —
  **deprecated/misleading**, не использовать как expected. Example 4\* — только `x4^2` fixture,
  `tm1Coeff` (Gamma-prefactor) **не** reducer-коэффициент. См. A12–A15.

## Дорожная карта дальше

1. ✅ **Pass 2G.1 — multi-sample normal-form records collector** (ЗАВЕРШЁН). `records.py`.
2. ✅ **Pass 2H — multivariate reconstruction / interpolation** (ЗАВЕРШЁН). `interpolate_multivariate`
   + диспетчер `reconstruct_coefficients`.
2b. ✅ **Pass 2H.1 — audit + hardening multivariate reconstruction** (ЗАВЕРШЁН). +8 тестов
   (bivariate poly/rational, denom без const-члена, corrupted holdout, insufficient points,
   numeric prefactor, union-support special-zero, bad/non-pivot records skip).
3. ✅ **Pass 2I.1 — result/diagnostics skeleton + strict Success gate** (ЗАВЕРШЁН). `result.py`
   (dataclasses + `build_reduction_result_from_reconstruction` + `result_to_wolfram_text`).
4. ✅ **Pass 2I.2 — reducer orchestration MVP** (ЗАВЕРШЁН). `reducer.py`:
   - `ReducerConfig` (primes/samples/labels|label_box/max_ibp_degree/tangent_degree_blocks/
     min_valid_records/preferred_masters/eps_direction), `ReducerRunDiagnostics` (n_labels/n_rows/
     n_records/n_reduced/n_bad_specializations/n_target_not_pivot/row+recon diag).
   - `reduce_family_once(family, target, config)` — один фикс. pass: enumerate → LF-flags
     (`is_locally_finite`) → rows (algebraic + coord surface-filtered + tangent если настроено) →
     `collect_normal_form_records` → `reconstruct_coefficients` → **strict gate из result.py**.
   - `reduce_rows_once(...)` — orchestration на готовых synthetic rows (для тестов без row-gen).
   - Failure mapping: нет reduced/target-not-pivot → `TargetNotReducible`; recon не валидируется →
     `InterpolationFailed`; non-LF терм → `NormalFormNotLocallyFinite`; Unknown-LF → failure;
     bad specialization → skip+count (не патчится). Run-diagnostics в `result.diagnostics.extra`.
   - **Не** реализовано (намеренно): `reduce_wolfram_style_input`, CLI-success, D4 e2e,
     adaptive expansion, perf-tuning. `Success` — только через gate 2I.1.
4b. ✅ **Pass D4.1/D4.2/D4.3 — D4 vertical validation** (ЗАВЕРШЕНЫ). Row-span certificate (3
   точки) → диагноз rank-poisoning → **max-rank record selection** (`reconstruction.py`, A29)
   → первый честный D4 full-config `Success` с 3-term LF basis `{M1,M2,M3}` (см.
   `notes/D4_STATUS.md`).
4c. ✅ **Pass D4.4 — D4 acceptance + equivalence certificate** (ЗАВЕРШЁН). Generic
   `certificate.verify_reduction_relation_mod_p` (без Success-стампа); reference-сертификат
   переведён на helper; reducer-output сертификат + эквивалентность reference↔reducer modulo
   row span (opt-in); диагностика M4/M5 → {M1,M2,M3}. `protected_masters`/forced basis —
   **не потребовался** для математического Success. См. `notes/D4_STATUS.md` §D4.4.
4d. ✅ **Pass Verify.1 (= D4.5) — certificate gate в редьюсере** (ЗАВЕРШЁН). Default-on
   `require_certificate_for_success`; auto/explicit `certificate_points` + `certificate_primes`;
   `certificate_rank_policy="selected_rank"`; классификация точек (rank-filtered /
   **rank-exceeded = жёсткий провал** / bad / pass / fail), `certificate_status`
   Passed/Failed/Insufficient/NotRun; маппинг на `VerificationFailed`; regression
   product-grid false-success + corrupt-coefficient + rank-exceeded + mixed-points. A31.
5. **Pass 2I.3 (СЛЕДУЮЩИЙ, рекомендован) — `reduce_wolfram_style_input` API + Wolfram-in/out
   склейка** поверх `reduce_family_once` (parse text → config defaults → result →
   `wolfram_style_text`), всё ещё без adaptive/CLI-success.
6. **Pass 2I.4+ — adaptive schedule, CLI success; опционально forced/protected-basis mode**,
   если когда-нибудь потребуется именно reference-базис (для Success не нужен).

## Критические инварианты (неизменны)

- `Success` **только** при `AllLocallyFinite=True` + reconstruction провалидирован на независимых
  samples; `formal_success ≠ success`.
- **`Success` требует row-span certificate `Passed` в off-sample rank-generic точках**
  (Verify.1, default-on); on-grid holdout сам по себе — НЕ независимая валидация (A30/A31).
- Нелокально финитные термы в normal form → диагностика/Failure, **не** Success.
- SymPy разрешён только в parse / tangent-setup / factorization / reconstruction; **не** в
  row-gen / RREF hot loop (там int mod p).
- Нет hardcode `N/M/x4/x7/H/A/B/C` и validation-кейсов; при неоднозначности — conservative
  `Unknown`/`Failure`.
- Bad specialization отвергается (skip+count), **не** патчится; нет reconstruction по одной точке.

## Pass 2I.3 — Public API layer (CLOSED, 2026)

Status: **complete, verified**.

- `parametric_ibp_lf_reducer/api.py` — public entry points:
  - `reduce_wolfram_style_input(...)` — parses Wolfram-like text, runs `reduce_family_once`, returns `ReducerRunResult`;
  - `reduce_wolfram_style_input_to_text(...)` — same pipeline, renders Wolfram-like text output;
  - explicit-family variants honor `ParserNeedsExplicitFamily` behavior (raise without family, succeed with one).
- API is re-exported from the package root (`__init__.py`); no CLI/Mathematica/Wolfram runtime dependency.
- Tests: `tests/test_api.py` — 5 tests covering callable surface, text→text round-trip, explicit-family requirement, and error propagation.

### Final verification (Pass 2I.3 gate)
- `python -m pytest` → **208 passed, 3 skipped** (skips: pre-existing, unrelated to API layer).
- `ruff check .` → **All checks passed!**
- API tests present and green: `python -m pytest tests/test_api.py -q` → 5 passed.

## Example 4* exploratory run (2026-07-08)

- `examples/example4_star_input.wl.txt` → **Success**, certificate **Passed 3/3**
  (selected_rank=7257, n_records=36, 0 skipped), 2-term all-LF reduction:
  `x4*x7/G0` и `x4*x7/G2` (коэффициенты в `validation/example4_star_result.m`).
- Детали и caveats: `notes/example4_star_exploratory.md`. Exploratory, не baseline;
  known-value-only пример (политика docs/05 соблюдена: без сравнения с deprecated
  артефактами, ε-разложение не использовалось как коэффициент).

## Corrected Example 4* (v0.1.1 candidate)

- `examples/example4_star_corrected_input.wl.txt`: исправленный множитель
  `15*ep + 24*ep*x7`, по линейности `15*ep*J[{0,0,0,0,0,0,0}] +
  24*ep*J[{0,1,0,0,0,0,0}]` (оркестрация: `scripts/run_example4_star_corrected.py`;
  ядро `src/` не менялось — обработка `lhs_terms` полностью generic).
- Результат: **Success**, certificate **Passed**, `selected_rank=9924`;
  артефакты: `validation/example4_star_corrected_result.m`,
  `validation/example4_star_corrected_diagnostics.json`.
- Тесты: `tests/test_example4_star_corrected.py`. Статус: v0.1.1 candidate,
  тег локальный, **не запушено**.

## Perf status (Perf.0–Perf.6, 2026-07; детали — `notes/PERF_STATUS.md`)

- **Perf.0 — stage timing diagnostics** (`timing.py`, `StageTimings`): чистая наблюдаемость,
  wall-clock по стадиям в `diagnostics.extra["timings"]` + CLI JSON. Математика/гейты не тронуты.
- **Perf.1 — ranking hoist** (принят): ranking вынесен из per-record цикла, строится один раз
  на прогон (`ranking_once`); verdict-equality harness подтвердил идентичность рангов/лейблов.
- **Perf.2 — ranking profiling** (принят, без дальнейшего копания): ranking — не доминанта
  на representative configs.
- **Perf.3 — parallel normal-form records (ОТРИЦАТЕЛЬНЫЙ результат, принят):**
  `collect_normal_form_records(..., jobs=N)` / `ReducerConfig.jobs` / CLI `--jobs` —
  корректно (равенство serial↔parallel: `tests/test_perf3_jobs_equality.py`, 12 тестов),
  но на Windows (`spawn`) старт пула + импорты (0.68 s @1 worker … 1.56 s @8) больше всей
  параллелизуемой record-работы (~0.96 s на heavy 2×49). Замеры: fast 2×9 jobs=1 0.814 s
  vs jobs=4 3.929 s; heavy 2×49 jobs=1 1.636 s vs jobs=2 2.693 s / jobs=4 4.107 s / jobs=8
  6.329 s. **Default `jobs=1` остаётся; `--jobs` — experimental**, пересматривать только
  когда record-работа на прогон достигнет десятков секунд.
- Tests после Perf.3: **245 passed, 7 skipped**; `ruff check .` — clean.
- **Perf.4 — heavy profile Corrected Example 4*** (принят): script-level stage timings в
  `scripts/run_example4_star_corrected.py` (`perf4_timings` в diagnostics JSON). Heavy run
  ~2h24m; доминанты: `rref_mod_p` (5631.8 s суммарно по двум таргетам, ~77% редукций) и
  certificate-работа (~40% wall). Вывод: цель Perf.5 — multi-target reuse, НЕ multiprocessing.
- **Perf.5 — multi-target / linear-LHS normal-form reuse (принят):**
  `collect_normal_form_records_multi` / `reduce_rows_multi` — ОДИН shared pipeline
  (ranking, assemble, `rref_mod_p`, record collection, certificate points) на общей
  row system для нескольких таргетов; per-target только selection / reconstruction /
  Success gate / certificate verdicts. Heavy run corrected Example 4*: **wall ~1h22m**
  (было ~2h24m, ~1.75x); `rref_mod_p` 2715.1 s ОДИН раз (было 5631.8 s);
  результаты идентичны certified baseline (2 terms, rank 9924, certificate Passed 5/5).
  Тесты: `tests/test_perf5_multi_target.py` — 15 passed (multi↔serial equality);
  full suite + ruff clean. Остаточные hotspots: один большой mod-p RREF (~2715 s)
  и certificate-работа (~2070 s) — дальше нужны новые kernel-дизайны, не reshuffling.
- **Perf.6 — certificate-point RREF reuse (принят):** combined certificate
  переиспользует RREF-ы, уже посчитанные для multi-pass certificate points
  (`rref_cache` в `_run_certificate_step` / cache-aware путь в `certificate.py`;
  cache miss = прежнее поведение). Heavy run: лог
  `combined certificate: 5 points, reusing 3 RREF(s) ... (Perf.6)`;
  `combined_certificate` **518.7 s** (было 1293.3 s, ~2.5x на стадии),
  wall ~1h15m; результаты идентичны certified baseline (2 terms, rank 9924
  на всех 5 точках, certificate Passed 5/5). Тесты: certificate-gate +
  perf-сьюты + full suite + ruff — clean. Остаточный hotspot: сам shared
  `rref_mod_p` (~2900 s, ~2/3 wall) — дальше только быстрый mod-p RREF kernel.

## Perf.7/Perf.8 merged to main (2026-07-11, merge 3a70cef, no tag)

- Branch `perf/rref-backend-prototype` merged --no-ff into `main`; pushed.
- NO release: no tag, no v0.1.4, version stays as-is.
- RREF backend default is still `"dict"` (`DEFAULT_RREF_BACKEND` in
  `sparse_rref.py`); `int_sparse_experimental` is opt-in only (0.69x on
  the real medium matrix, ~0.90x synthetic). Dict backend NOT removed.
- `collect_stats=True` counters + `scripts/bench_rref_backends.py` +
  `scripts/profile_rref_real_matrix.py` now on main; real profile at
  `validation/rref_real_matrix_profile.json` (fill-in 10.5x, rows stay
  sparse). Candidate B (pure-Python int-array rows) measured 1.5-1.9x
  SLOWER and rejected — see `docs/RREF_BACKEND_PLAN.md`.
- Next RREF work is DESIGN ONLY: candidate C/D triggers and a Numba/E
  merge-kernel sketch; no implementation under the pure-Python
  constraint.

## Perf.10 — Numba RREF backend prototype (branch `perf/numba-rref-backend`, НЕ merged)

- Pure-Python constraint снят для *опционального* бэкенда: новый модуль
  `src/parametric_ibp_lf_reducer/sparse_rref_numba.py`, opt-in backend
  `numba_int_array_experimental` (lazy import — пакет без numba работает
  как раньше; numba уже объявлен в extras `speed = ["numba>=0.59"]`).
- Default остаётся `"dict"`; LF/certificate gates НЕ тронуты; математика
  бит-в-бит повторяет dict pivot loop (тот же выбор пивотов, порядок
  элиминации, 1 инверсия на пивот; guard `p < 2**31` для int64).
- Замеры (editable install, эта машина): synthetic fast bench —
  tiny 10.07x, medium 1000x800 **25.34x** vs dict; реальная ранкинг-матрица
  512x917 — dict 0.828s → numba **0.083s** (~10x); JSON обновлён в
  `validation/rref_real_matrix_profile.json`. Подробности:
  `notes/PERF_STATUS.md` (секция Perf.10).
- Тесты: `tests/test_rref_numba_backend.py` (skip без numba) — parity vs
  dict на random/dense/rank-deficient/partial-order/малых простых +
  граница `2**31 - 1`, reject `>= 2**31`, stats/plain-int, end-to-end
  `modular_normal_form` parity. Full suite green, ruff clean.
- Коммит `ab861e1`; ветка слита `--no-ff` в `main` (`81ec173`),
  запушена в origin, ветка удалена.

## Perf.11 — выбор RREF-бэкенда через ReducerConfig/API/CLI (ветка `perf/rref-backend-plumbing`)

- `ReducerConfig.rref_backend` (default `"dict"`) прокинут через
  `records.py` (serial-путь и worker-процессы через per-point
  context), `modular_normal_form.py`, `reducer.py` и `certificate.py`
  — ranking, normal-form и certificate RREF уважают выбор. Только
  селекция бэкенда: выбор пивотов, порядок элиминации и LF/certificate
  gates НЕ тронуты; результаты идентичны по построению (parity-тесты
  Perf.10 это дополнительно фиксируют).
- Имя бэкенда валидируется по `sparse_rref.RREF_BACKENDS` до начала
  работы; неизвестное имя — быстрая понятная ошибка.
- API принимает `rref_backend` и передаёт в `ReducerConfig` как есть.
- CLI: флаг `--rref-backend` (choices = `RREF_BACKENDS`, default
  `dict`); неизвестное значение — usage error (`EXIT_USAGE`).
- Тесты: `tests/test_cli.py` (+2: accepted/rejected значение), numba
  parity suite без изменений. Full suite green, ruff clean.
- Ветка слита `--no-ff` в `main`, запушена в origin, ветка удалена.

## Perf.12 — `rref_backend="auto"`: эвристический выбор dict vs numba (ветка `perf/rref-auto-backend`)

- Новое имя `"auto"` (`AUTO_RREF_BACKEND`); валидация теперь по
  `RREF_BACKEND_CHOICES = (*RREF_BACKENDS, "auto")`. Default остаётся
  `"dict"` — auto строго opt-in.
- `select_rref_backend()` резолвит `"auto"` per-matrix, до начала
  работы (serial и worker-процессы ведут себя одинаково): numba только
  если доступна, `prime < 2**31` и матрица проходит все гейты
  `AUTO_RREF_THRESHOLDS` (`min_rows: 500`, `min_cols: 400`,
  `min_nnz: 3000`; консервативно, может меняться); иначе `"dict"`.
  Отсутствие numba — тихий fallback ТОЛЬКО для `"auto"`; явный запрос
  бэкенда по-прежнему падает быстро и явно.
- Результаты backend-идентичны (parity Perf.10); в stats пишется
  решение: `requested_rref_backend`, `selected_rref_backend`,
  `backend_selection_reason`, `numba_available`,
  `auto_thresholds_used`.
- CLI: `--rref-backend auto` принимается (в help помечен
  experimental, Perf.12).
- Тесты: новый `tests/test_rref_backend_auto.py` (гейты
  размера/prime, fallback без numba, явный запрос не подменяется,
  end-to-end) + `tests/test_cli.py` для `auto`. Full suite green,
  ruff clean.

## Perf.13 — полная валидация auto/Numba RREF-бэкенда на полном боксе

- `scripts/validate_rref_backend_full.py`, полный бокс 972 лейбла
  (12360 строк), оба таргета, бэкенды dict / numba / auto.
- Итог: dict 3963с → numba 804с (4.93×) → auto 766с (5.17×, auto
  выбирает numba по порогам). `identical_across_backends: true`,
  mismatches пусто; сертификат Passed 3/3, rank 9924, 36/36 записей,
  коэффициенты побитово совпадают с сертифицированным Example 4★.
- Вывод: auto безопасен на реальной нагрузке, ~5× ускорение
  end-to-end; пороги оставлены консервативными. Perf.13 закрыт.
