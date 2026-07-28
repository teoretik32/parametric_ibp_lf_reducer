# External Int2 Method.13 — complete mixed-boundary surface validation

Дата: 2026-07-28. Ветка `feature/external-int2-hyperint-integration`.
НЕ закоммичено — ожидает review (по условию задачи).

Артефакты: `docs/TORIC_SURFACE_VALIDATION.md` (вывод критериев),
`scripts/audit_surface_validation.py`,
`validation/external_int2_method13_surface_audit.json`,
`validation/surface_reaudit_manifest.json`, тесты
`tests/test_method13_surface_validation.py` (+ правки существующих).

## Phase 0 — provenance

- `validation/external_int2_lf_result.m`: добавлены
  `FormalRowSpanCertificate -> "Passed"`, `IntegralIdentityStatus -> "Revoked"`,
  `SurfaceValidationStatus -> "Failed"` (различение: модульный row-span
  сертификат жив, интегральное тождество отозвано, surface-обоснование строк
  провалено).
- Контракт результата (`result.py`): новое поле
  `SurfaceValidationStatus ∈ {Passed, Failed, Unknown, LegacyUnchecked}`,
  default **LegacyUnchecked**. `Success` теперь требует: reconstruction
  verified + независимый row-span certificate + все RHS-термы LF=True +
  `SurfaceValidationStatus=Passed`. Новый typed отказ
  `SurfaceValidationNotPassed`. Живой редьюсер (`reducer._finalize_target`)
  ставит `Passed` — его строки прошли исправленные фильтры этого же прогона;
  legacy-сборки (например `build_reduction_result_from_reconstruction` из
  записанных коэффициентов) получают default и НЕ проходят гейт молча.

## Phase A — полный набор лучей (детерминированный)

`newton_wall_normals` + `complete_polyhedral_rays` (введены в 12R,
финализированы здесь): координатные ±e_i в порядке переменных, затем
полиэдральные лучи в отсортированном лексикографическом порядке; примитивная
целочисленная нормализация; дедупликация; бюджет читается в момент вызова
(патчится в тестах) — переполнение ⇒ `complete=False` ⇒ `is_locally_finite`
и surface-фильтры возвращают `Unknown`, никогда `True`; случайная сетка может
дать только дополнительный `False`. Для External Int2 набор — ровно 18
ожидаемых лучей, включая `(0,-1,-1)`; тест пиннит полный список.

## Phase B — координатный поверхностный критерий (вывод в docs/)

Для `∫ d/dx_i (P F)`: поверхностный член на фасете x_i — это
`(N-1)`-мерный поперечный интеграл. В лог-координатах поток поля
`V = (P F / x_i) Π x_j e_i` через сферу на бесконечности в направлении луча
`d` масштабируется как `λ^{score(P F, d) − d_i}`, откуда точный критерий:

```
для каждого луча d полного набора с d_i ≠ 0:
    surface_score(d, i) = score(P F, d) − d_i > 0   (строго)
```

На ±e_i это в точности старый компонентный тест (`exp_zero > 0`,
`exp_inf < 0`) — production-фильтр `coordinate_primitive_surface_free`
заменён на полный фасетный критерий (компонент-локальный допуск удалён).
Точная арифметика: ноль ⇒ Failed; неразрешимый знак ⇒ Unknown; бюджет ⇒
Unknown. Обе политики (limit/chamber) используют ОДИН полный набор лучей.
Regression: строка 12R (`[-1,0,0,-3,0,-3,0]`, x5, P=x5², луч `(1,-1,0)`,
score −1) отвергается при обеих политиках.

## Phase C — нормальный поток для div(Q F)

В лог-координатах поток через грань луча `d` управляется **нормальной
компонентой** `N_d(x) = Σ_i d_i Q_i/x_i` (компоненты с d_i = 0 тангенциальны
и не дают условий). `N_d` собирается точной арифметикой коэффициентов
(`ParamExpr`): мономы разных компонент могут сокращаться ТОЛЬКО точно —
не по модулю и не численно. Каждый выживший моном обязан иметь
`score(x^m F, d) > 0`. Синтетика в тестах: `Q=(x,−y)` — компоненты по
отдельности маргинальны, поток сокращается ⇒ валидно; `Q=(x,y)` —
покомпонентно ок, смешанный поток выживает со score 0 ⇒ невалидно;
коэффициенты `(1, 2147483646)` сокращались бы mod 2147483647, но не точно ⇒
невалидно; `(ep, −ep)` сокращается точно в кольце коэффициентов ⇒ валидно.

## Phase D — пересбор строк Int2 Level-0 (без решения)

Старый (реплицированный) набор воспроизводит Method.10 **точно**:
limit 46737 (12288 alg / 25232 coord / 9217 tangent), chamber 49439
(12288 / 26256 / 10895); движок сверен с production-фильтрами на 2×300
случайных парах — 0 расхождений.

| policy | old | new | invalidated (coord/tangent) | gained (tangent) |
|---|---|---|---|---|
| limit | 46737 | **37365** | 9523 (6841/2682) = 20.4% | 151 |
| chamber | 49439 | **40034** | 9532 (6229/3303) = 19.3% | 127 |

- «Gained» — поля, которые старый по-термный тест отвергал, а точный
  нормальный поток (с точными сокращениями и без тангенциальных компонент)
  легитимно принимает.
- Chamber-only (2702 = 1024 coord + 1678 tangent): **1554 остаются валидны**
  (589 coord + 965 tangent), 1148 отвергнуты.
- Witness-скрининг (Method.7, 6 точек): breaks **1685 → 1383** среди
  выживших chamber-only строк — обоснование Method.10 rerun сохраняется в
  ослабленном виде.
- Сравнение с предварительными оценками 12R: chamber-only 1148 vs 1170,
  limit 9523 vs 9865. Разница объяснена: 12R-оценка была строже по tangent
  (по-термно, без точных сокращений, с тангенциальными d_i=0 компонентами);
  координатная часть критериев идентична.
- LF-вердикты в артефакте: Target False, L1/L2/L3 True, L4 False.

Preflight (единственный RREF аудита, 1 простое × 3 сэмпла на исправленном
chamber-наборе, 40034 строки): **Feasible(modular) на всех трёх точках**,
ранг 20690 (стабилен), LF-True меток в боксе при строгом гейте — 1741/3072.
То есть исправленная система по-прежнему генерически LF-feasible; Phase F
item 20 (bounded reconstruction + новый сертификат) **разрешён по гейту, но
не запускался** — ждёт review этого аудита (масштаб: ~40 сэмплов × ~4 мин
RREF). Positive-control 2mh-семейство от Леонида (item 21) не требуется:
preflight не Obstructed.

## Phase E — project-wide re-audit

`validation/surface_reaudit_manifest.json` (все мастера всех кейсов
LF=True под полным набором лучей; все таргеты честно False):

| кейс | rows old→new (invalid.) | статус |
|---|---|---|
| D4 CLI | 2092→1176 (918; +2 gained) | **UnchangedValid** — пере-сертифицирован живьём в этой сессии (in-span сертификаты на 3 точках + reducer/CLI Success с `SurfaceValidationStatus=Passed` против исправленной системы строк) |
| notebook Example 2 | 1116→954 (162) | **UnchangedValid** — живой adaptive Success + certificate Passed + коэффициенты symbolически равны notebook-fixture |
| notebook Example 1 | — (fixture-only) | **LegacyInsufficientData** — закоммиченного Success-артефакта нет; мастера LF=True подтверждены |
| Example 4* | 12360→8140 (**4220 = 34%**) | **NeedsRecompute** — закоммиченный Success legacy; e2e re-run env-gated; независимое known-value свидетельство записано |
| External Int1 | 818→686 (132 = 16%) | **NeedsRecompute** — legacy сертификат; независимое численное/Laurent свидетельство (item 19) сохранено в манифесте |

## Ограничения и запреты

- Реконструкция/новый LF-поиск НЕ запускались (Phase F — по итогам preflight
  и review).
- RREF в аудите — только preflight.
- Никакие обязательства случайных лучей не сертифицируют True.
- Ничего не закоммичено до review.
