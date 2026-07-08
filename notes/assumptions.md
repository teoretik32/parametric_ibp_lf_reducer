# Assumptions — `parametric_ibp_lf_reducer`

Инженерные допущения, принятые без уточняющих вопросов (по требованию задачи). Каждое
задокументировано с обоснованием; при поступлении иных данных — пересматриваем.

## A1. Target label = базовый интегранд
Целевой интеграл = `(n=0, m=0)` × `TargetMultiplier`, т.е. сам `F_base` из входа.
*Обоснование:* spec §3 определяет исходный интегранд как
`Product[x_i^MonomialExponents] * Product[G_l^PolynomialExponents] * TargetMultiplier`;
редукция ищет его представление через LF-мастера.

## A2. Один основной регулятор `epsilon`
MVP оптимизирован под один регулятор `ep`, но архитектура (list `Regulators`, параметрические
коэффициенты) допускает несколько.
*Обоснование:* spec §2 явно разрешает «MVP может быть оптимизирован под один epsilon с
архитектурной возможностью расширения».

## A3. LF-тест — ray-scaling при `eps=0` + random safety net
`is_locally_finite(label)` возвращает `True` только если `base_score(rho, eps=0) > 0` для всех
детерминированных лучей-кандидатов **и** случайная выборка лучей не нашла дивергенцию.
Иначе `False` (детерминированная дивергенция) или `"Unknown"` (найдена только random-net'ом
или тест неполон). `"Unknown"` **никогда** не считается `True`; `Success` требует `True`.
*Обоснование:* spec §4.4, §8; method review §5. Полный нормальный фан в общем N дорог —
safety net обеспечивает soundness (см. `risk_register.md` R1).

## A4. Reconstruction: multi-prime CRT → точные рациональные → функциональная интерполяция
Схема: для набора рациональных param-точек получаем точное рациональное значение каждого
коэффициента через multi-prime CRT + rational reconstruction; затем интерполируем rational
function (dense-rational ansatz с degree search); проверяем на независимых samples.
*Обоснование:* spec §5.10, method review §7. Разделяет modular-integer recovery и
функциональную интерполяцию — проще и строже для MVP, чем полный Zippel (тот — full target).

## A5. Границы использования SymPy
SymPy разрешён в: `input_parser` (парсинг), `tangent_fields` (разовый малый syzygy nullspace
в setup), `wolfram_text_export` (финальная факторизация), `reconstruction`/`interpolation`
(финал). **Запрещён** в `row_generation`/`sparse_rref` hot loop (там только int mod p).
*Обоснование:* CLAUDE.md п.7, spec §9.3.

## A6. scipy — опциональная зависимость
Используется только для фасетных нормалей Newton-политопов (`ConvexHull.equations`). Есть
чистый fallback (координатные + edge-нормали + random net) без scipy.
*Обоснование:* обязательные deps — только `sympy`, `numpy` (`pyproject.toml`). Не тянем scipy
в обязательные.

## A7. Кейсы 11.2 / 11.4 / 11.5 — без входных families
В docs даны только ожидаемые LF-базисы/коэффициенты, но не входные полиномы (A/B/C,
ID4example2, ID3example3). Оформляем как `xfail`/`skip`-фикстуры «pending input family».
MVP-приёмка на них не завязана (якорь — 11.1 + 11.3).
*Обоснование:* `docs/03_validation_cases_ru.md` §11.2/11.4/11.5 не содержат входных семейств.

## A8. `epsilon-direction` — для surface-регуляризации, не для LF
`epsilon-direction` (`minus`/`plus` из Options/CLI) используется только в surface-free тестах
регулируемой области. Финальная LF-проверка мастеров всегда строго при `eps=0`.
*Обоснование:* spec §7.3, §8; method review §4.3.

## A9. Область интегрирования — `R_+^N`
По умолчанию и в MVP. Bulk singularities вне boundary не обрабатываются автоматически без
assumptions пользователя, гарантирующих отсутствие внутренних нулей (тогда LF → `Unknown`).
*Обоснование:* spec §2 (ограничения MVP), §4.4.

## A10. Коэффициенты `G_l` — из `Q(params)`
Коэффициенты полиномов `G_l` — рациональные функции внешних параметров; степени `x_i`, `G_l`
— аффинные/рационально-простые по регуляторам/параметрам.
*Обоснование:* spec §2.

## Классификация notebook examples (ИСПРАВЛЕНО)
Два разных типа:
- **Example 1, Example 2 — known-decomposition fixtures.** Есть initial integrand +
  `JIntegrandList` + `JCoefficientListt` (явное разложение по LF-интеграндам). В Pass 1A —
  только parser/coefficient тесты, **не** full e2e reduction.
- **Example 3*, Example 4* — known-value-only.** Известен только аналитический ответ /
  epsilon expansion самого интеграла; LF-разложения **нет** (нет `JIntegrandList` /
  `JCoefficientListt`). Они **не связаны** с Example 1/2 и **не являются** missing families
  для main-spec 11.4/11.5. В Pass 1A — только parser/sparse-poly фикстуры.

## A12. Example 1 (I4exampl1) — known-decomposition, alt 6-term, не canonical
Та же base family, что и canonical D=4 (11.3), notebook даёт явное **6-term** разложение
(`JIntegrandList`/`JCoefficientListt`). Не заменяет canonical `validation/expected_d4_coefficients.json`
(5-term). Pass 1A — parser/coeff фикстура; позже — опциональный preferred-master вариант.
*Обоснование:* `docs/04_notebook_examples_addendum_ru.md` п.1 + уточнение классификации.

## A13. Example 2 (I3exampl2) — known-decomposition, tentative 11.2
Concrete N=3 family с явным 5-term разложением. Возможный кандидат на regression 11.2, но это
решение принимается **позже** (не в Pass 1A); e2e остаётся `pending`/`xfail`. Mapping к старым
`A,B,C` документируется в самом тесте при активации.
*Обоснование:* addendum п.2 + уточнение классификации.

## A14. Example 3* (I4exampl3) — known-value-only, НЕ family для 11.4
Известен только аналитический ответ / epsilon expansion интеграла; LF-decomposition нет.
**Не** трактуется как missing input family для main-spec 11.4 (ID4example2). В Pass 1A — только
parser/sparse-poly фикстура. Файл `validation/id4example2_expected_one_term.json` (от старого
patch) **deprecated/misleading** и **не** используется как expected output для Example 3*.
Main-spec 11.4 остаётся `pending`/`xfail` до отдельной поставки настоящей input family +
LF-decomposition ожиданий.
*Обоснование:* поправка пользователя от 2026-07-06 (перекрывает addendum п.3).

## A15. Example 4* (tm1Int) — known-value-only, НЕ family для 11.5; `x4^2` фикстура
Известен только аналитический ответ; LF-decomposition нет. **Не** трактуется как missing input
family для main-spec 11.5 (ID3example3). В Pass 1A используется **только** как parser/sparse-poly
фикстура — в частности для проверки монома `x4^2` в `G3`. `tm1Coeff` (Gamma-prefactor) **не**
используется как reducer-коэффициент. Файл `validation/id3example3_expected_basis_inclusion.json`
(от старого patch) **deprecated/misleading** и **не** используется как expected output для
Example 4*. Main-spec 11.5 остаётся `pending`/`xfail` до отдельной поставки настоящей input
family + LF-decomposition ожиданий.
*Обоснование:* поправка пользователя от 2026-07-06 (перекрывает addendum п.4).

## A16. Границы Pass 1A
Pass 1A реализует только: package skeleton, `input_parser` (explicit family), `coefficients`
(`ParamExpr`+`eval_mod_p`), `sparse_poly` (canonical sparse dict, произвольный N, derivative,
valuation, `eval_mod_p`) и **минимальный** контейнер `family`. НЕ реализуются (Pass 1B+):
`valuations`, `surface`, `row_generation`, `tangent_fields`, `ranking`, `sparse_rref`,
`reconstruction`, reducer success-pipeline, полноценные `labels`/`family`-методы редукции.
CLI/reducer-заглушки **не** возвращают `Success`. При неоднозначности — conservative
`Failure`/`Unknown`.
*Обоснование:* явные инструкции текущего задания.

## A17. Парсинг Wolfram-like выражений
Математические подвыражения (полиномы, экспоненты, коэффициенты) парсятся структурным
recursive-descent парсером (association/list/string + raw-math), затем raw-math один раз
конвертируется в SymPy (`^`→`**`, символы из объявленных Variables/Parameters). Полиномы →
`SparsePoly` через `Poly(expr, *vars)`; коэффициенты/экспоненты → `ParamExpr` через
`fraction(together(...))`. Это setup-фаза (SymPy разрешён); hot loop не затрагивается.
Проверка: free-symbols полинома ⊆ Variables∪Parameters, экспонент/коэффициентов ⊆ Parameters;
иначе `ParserError`.
*Обоснование:* CLAUDE.md, spec §5.1; conservative parsing.

## A18. Границы и решения Pass 1B
Pass 1B добавляет только label/projection/export слой:
- `labels.py`: `Label=tuple[int,...]` длины `N+M`; `enumerate_box` (spec диапазонов — единая
  пара `(lo,hi)` broadcast на все оси **или** per-axis последовательность; произвольные N/M);
  `LabelIndex` (биекция label↔id); `label_complexity` — **структурная утилита**, не решение о
  masters (ranking — позже).
- `family`: `label_to_factor` возвращает **относительный** фактор (целые сдвиги `n_i,m_l`
  относительно base), совпадающий с нотацией `lf_basis` в validation (`x2*x3/(G0^2*G1)`);
  `exponent_at_label` → `(e_i,f_l)` как `ParamExpr` (`base + сдвиг`).
- `wolfram_text_export`: коэффициенты через `sympy.factor` + `sstr().replace("**","^")`
  (rational как `p/q`, `^` не `**`); знаменатель интегранда в скобках при >1 факторе.
- `specialize(sample,prime)` — честный `NotImplementedError` (полная FamilyModP нужна row-gen
  слою). Реальный частичный кусок — `specialize_polynomials` (только `G_l` mod p через
  `SparsePoly.eval_mod_p`, без экспонент).
НЕ реализовано в Pass 1B: `valuations`, `surface`, local-finiteness decision, `tangent_fields`,
`row_generation`, `ranking`, `finite_field`, `sparse_rref`, `reconstruction`, reducer success.
*Обоснование:* явные инструкции текущего задания; conservative stubs вместо fake behavior.

## A19. Решения Pass 1C (valuations + surface)
- **base_score** вдоль луча `rho`: `Σ_i rho_i(e_i+1) + Σ_l f_l·val_rho(G_l)`, экспоненты при
  `eps=0`. `val_rho` = min по support (min-конвенция, без учёта cancellation).
- **Candidate rays (MVP)**: `±e_i` (coord0/coordInf) + support-мономы каждого `G_l` (и их
  отрицания) + диагональ по переменным каждого `G_l` (и отрицание); primitive-reduced, dedup.
  Полнота не гарантирована в общем N → random safety-net + Unknown (см. R1).
- **`is_locally_finite` → True/False/"Unknown"** (строго при `eps=0`):
  `True` только если по всем candidate- и random-лучам `base_score>0`, **все** экспоненты
  numeric при `eps=0`, и все знаменатели provably positive (`_bulk_safe`). Любой `base_score≤0`
  ⇒ `False` (строгое правило: ровно 0 — **не** LF). Symbolic/неполнота/возможная bulk-ноль ⇒
  `"Unknown"` (никогда True).
- **`_bulk_safe`**: учитываются только полиномы-знаменатели (`f0<0` или возможно-отрицательный);
  positivity через SymPy с symbol'ами, помеченными positive из assumptions вида `X > 0`.
- **surface.regulated_sign**: для surface-тестов используется предел регулятора `eps→0^∓`
  (leading value при `eps=0`, затем first-order коэффициент `eps` + направление). Это **только**
  для surface (spec §7.3); финальная LF-проверка всё равно строго при `eps=0`. `eps_direction`
  по умолчанию `"minus"`.
- **coordinate_primitive_surface_free**: component-local — только `x_i=0` (`p_i+e_i+Σf_l·minpow_i>0`)
  и `x_i=∞` (`p_i+e_i+Σf_l·maxpow_i<0`); **не** требует зануления по всем toric-лучам (spec §7.1).
- **vector_field_surface_free**: нормальный flux по toric candidate-лучам; `neg`/`zero`/`unknown`
  обрабатываются консервативно (`zero` ⇒ строка не surface-free ⇒ False; `unknown` ⇒ Unknown).
- Row generation в Pass 1C **не** реализовано; tangency 11.1 проверяется прямой sparse-poly
  алгеброй (`Σ_i Q_i ∂_iG = 0`), без отдельного tangent_fields-модуля.
*Обоснование:* spec §4.4/§5.4/§7/§8, method review §4–§5; conservative-by-default.

## A20. Решения Pass 2A (row_generation: algebraic + coordinate IBP)
- **Row** — разреженная связь `{label: ParamExpr}`, `Σ coeff·J[label]=0`, над целочисленными
  labels; множество переменных системы = объединение labels во всех строках (не фикс. box).
- **Algebraic rows** `J(n,m) − Σ_a c_{l,a} J(n+a, m−e_l) = 0` — точные тождества интегралов,
  **всегда** принимаются (без surface-check).
- **Coordinate IBP rows** `0=∫∂_{x_i}(P F)`, ансатц `P=x^p`, degree `0..max_degree`
  (`MaxIBPDegree`). Раскрытие: `(p_i+e_i)J[n+p−1_i,m] + Σ_l Σ_{b:b_i>0} f_l b_i c_{l,b}
  J[n+p+b−1_i, m−e_l]`. Принимается **только** если
  `coordinate_primitive_surface_free(...)==True`; иначе reject с причиной `surface_not_free`
  (False) или **консервативно** `surface_unknown` (Unknown). Пустые строки → `trivial_row`.
- **eps_direction** по умолчанию `"minus"` (регулируемая область surface); прокидывается в
  surface-check. Финальная LF-проверка мастеров всё равно строго при `eps=0` (Pass 1C).
- **Dedup** — по точному множеству термов (`frozenset((label, ParamExpr))`). Строки,
  эквивалентные с точностью до общего скаляра, в MVP **не** дедуплицируются (важно для
  производительности, не для корректности).
- Tangent/syzygy rows **не** реализованы в Pass 2A (Pass 2B). Rejected-строки сохраняются с
  провенансом `{seed,var,P}`/`{seed,poly}` и причиной — для диагностики.
*Обоснование:* spec §5.6/§7, method review §2/§4; conservative surface-фильтр.

## A21. Решения Pass 2B (tangent_fields, SymPy MVP backend)
- **TangentField** = `(components Q_i, multipliers H_l, degree_block)`; `Q·∇G_l = H_l G_l` для
  **всех** l. Поле принимается только если `is_tangent(family)` (все defects `Σ_i Q_i∂_iG_l −
  H_l G_l` нулевые) — проверяется по построению и дополнительно guard'ом.
- **Solver**: degree-блоки `(d_Q,d_H)`; ансатц `Q_i` (deg≤d_Q), `H_l` (deg≤d_H) с неизвестными
  коэффициентами; тождество `Σ_i Q_i∂_iG_l − H_l G_l ≡ 0` в `x` → однородная линейная система;
  **nullspace через SymPy** над полем `Q(params)` (setup-фаза, не hot loop).
- **Dedup**: нулевое поле отбрасывается; поля, пропорциональные с точностью до **параметрического
  скаляра** (не зависящего от `x`), схлопываются (`_proportional` через 2×2 миноры коэфф.-векторов).
- **verify_tangent(family, Q)**: `Q·∇G_l` делится на `G_l` (через `sympy.cancel`; `x` не должен
  оставаться в знаменателе) → возвращает `H_l`; standalone-проверка (для 11.1 и как guard).
- **Ограничения SymPy MVP backend** (документированы):
  1. nullspace вычисляется для **generic** значений параметров; спец-точки, где зануляется
     ведущий коэффициент, отдельно не разбираются.
  2. Только dense single-block ансатц; нет Singular/Sage syzygy-backend'а.
  3. Пропорциональность/делимость через `sympy.simplify`/`cancel` — может быть дорогой для
     больших степеней (приемлемо: setup-фаза, малые степени).
  4. Финитная специализация полей mod p (для hot loop row-gen) — **не** здесь (Pass 2C+).
- **Не** реализовано в Pass 2B: `generate_tangent_ibp_rows`, интеграция
  `vector_field_surface_free` в row_generation, ranking/finite_field/RREF/reconstruction/reducer.
*Обоснование:* spec §5.7/§3, method review §3; SymPy разрешён в setup (CLAUDE.md п.7).

## A22. Решения Pass 2C (tangent IBP rows)
- **Раскрытие** `0=∫div(Q F_label)` для tangent-поля: `div(QF)=F[divQ + Σ_i e_i Q_i/x_i +
  Σ_l f_l H_l]`, где `e_i=a_i+n_i`, `f_l=b_l+m_l` (из `exponent_at_label`), а `H_l` — **сохранённые**
  `TangentField.multipliers`. **Никакого m-сдвига**: все термы сохраняют m-tuple source label.
  Делимость `Q·∂G/G` внутри строки **не** пересчитывается (используется stored H).
- **Q_i/x_i** для монома с `x_i`-степенью 0 → сдвиг `n_i−1` (валиден, не пропускается).
- **generate_tangent_ibp_rows**: поле проходит guard `is_tangent(family)` (по построению +
  проверка defect через SparsePoly-умножение, один раз на поле); non-tangent → reject
  `field_not_tangent`, строка **не** добавляется (controlled, без crash). Surface-gate —
  `vector_field_surface_free` (toric flux); verdict≠True → reject `surface_not_free`/`surface_unknown`.
  **LF промежуточных labels не фильтруется** (non-LF термы сохраняются).
- **Interop**: компоненты `Q_i` и `H_l` должны иметь те же `parameters`, что и family
  (иначе `SparsePoly` mul бросает mismatch). `generate_tangent_fields` уже это обеспечивает.
- Radial-подобные поля (напр. `Q=(x+y+1,0)`, `H=(1)`) алгебраически дают строку, но обычно
  **отклоняются surface-gate** (ненулевой boundary flux) — что и защищает от ложного тождества.
*Обоснование:* method review §3/§4, spec §5.6/§7.2.

## A23. Решения Pass 2D (ranking)
- **ranking = только порядок исключения колонок**, ничего не решает и **не удаляет** labels
  (`ordered` — перестановка входа). Индекс 0 = eliminate-first (pivot), хвост = prefer free (masters).
- **Tiers**: `0` target (исключается первым — нужна его normal form); `1` non-LF **или**
  `"Unknown"` (исключаются раньше любых LF, никогда не free); `2` generic LF; `3` preferred
  masters (пользовательские, самые «свободные», в хвосте).
- **Sort key** = `(tier, −complexity, label)`: внутри tier более сложные исключаются раньше,
  простейшие LF-интегранды дрейфуют в свободный хвост. `complexity` из `labels.label_complexity`.
- **Simple non-LF не остаётся free**: простой (низкая complexity) non-LF всё равно в tier 1 →
  раньше любого LF-мастера. Простота не «спасает» расходящийся интеграл.
- `"Unknown"` трактуется как non-LF (tier 1) — консервативно, не оставляем free (success
  требует True).
- `formal_success`/`success` здесь **не** вводятся; ranking — предпосылка для pivot-порядка RREF.
*Обоснование:* spec §5.8, method review §6 (Ошибка 4: не оставлять простые, но расходящиеся
интегралы свободными).

## A24. Решения Pass 2E (finite_field + minimal sparse_rref)
- **finite_field.py** — чистая int-арифметика mod p (без SymPy): `inv_mod` (Ферма, `pow(a,p-2,p)`),
  `powmod` (в т.ч. отрицательная степень через инверсию), `batch_inverse` (Montgomery trick,
  1 инверсия на пакет; ноль → `ZeroDivisionError`), `add/sub/mul_mod`, `is_probable_prime`
  (детерминированный Miller-Rabin для `n<3.3e24`), `generate_primes(count, upper)` (по убыванию).
- **sparse_rref.py** — минимальный sparse RREF над GF(p):
  - строки — sparse `dict[column,int]`; **column = любой hashable+orderable** ключ (int-индекс
    или сам label-tuple), чтобы matrix-assembly (2F) мог работать прямо над labels.
  - `rref_mod_p(rows, prime, column_order)` — pivot выбирается по `column_order` (порядок из
    ranking: eliminate-first); колонки вне списка пивотятся последними в sorted-порядке. Так
    low-priority колонки (мастера) остаются free.
  - поддерживает full RREF (каждая pivot-колонка только в своей строке; forward+backward);
    зависимые строки схлопываются (rank корректен). Система однородна (`A·J=0`).
  - `RREFResult`: `pivots{col→row}`, `pivot_order`, `free_cols`, `all_cols`, `rank`.
- **Не** реализовано в 2E: matrix assembly из `Row` (через `ParamExpr.eval_mod_p` на sample),
  извлечение normal form target, reconstruction, reducer Success. numba/parallel — full target.
*Обоснование:* spec §5.9; conservative minimal RREF, проверенный против dense-эталона в тестах.

## A25. Решения Pass 2F (matrix assembly + modular normal form)
- **assemble_rows_mod_p**: `Row` → sparse `dict[label,int]` через `ParamExpr.eval_mod_p(sample,p)`.
  Нулевой числитель → терм отбрасывается; **нулевой знаменатель → `BadSpecialization`** (sample
  отвергается целиком, **не** патчится, не пропускается тихо). Только int mod p (no floats).
- **modular_normal_form(family, rows, target, sample, prime, preferred_masters, lf_map)** →
  `NormalFormResult`: assemble → `rank_labels(target=target)` → `rref_mod_p(column_order=ranked.ordered)`;
  если `target` — pivot, normal form = `target = Σ (−v) J[free_col]` (RREF-строка target).
- **Статусы**: `Reduced` (target выражен через free), `TargetNotReducible` (target отсутствует/
  не pivot), `BadSpecialization`, `EmptySystem`.
- **formal_success** (target редуцирован в одной точке) — **отдельно** от success; LF термов
  только **диагностируется**: `all_terms_lf ∈ {True, False, "Unknown"}`, плюс `non_lf_terms`,
  `unknown_lf_terms`. **Success нигде не возвращается.** formal_success=True при non-LF/Unknown
  термах — это диагностика, не физический успех.
- **Детерминизм**: ranking + RREF детерминированы; `terms` строятся из sorted-строки.
- Reconstruction по многим (prime,sample) и Success pipeline — Pass 2G+.
*Обоснование:* spec §6/§7/§10; строгое разделение formal_success vs success (CLAUDE.md п.5–6).

## A26. Решения Pass 2G (reconstruction, 151 строка)
- **rational_reconstruction(a,m)**: Wang-стиль, `|p|,q ≤ sqrt(m/2)`; `None` если модуль мал
  (недостаточно primes) — не выдаёт неверное значение.
- **reconstruct_rational({prime:res})**: инкрементальный CRT + rational reconstruction; требует
  **стабильности** (совпадение на двух подряд primes) → защита от ложного восстановления по
  одному пункту.
- **collect_value_table(results)**: группировка Reduced-результатов по sample, восстановление
  точного рационального значения по many primes. **Union support с 0-fill** (отсутствующий терм
  = 0 в поле, а не «тихо выброшен»). Non-Reduced (в т.ч. `BadSpecialization`) **пропускаются и
  считаются** (`n_skipped`), не искажают результат.
- **interpolate_univariate**: dense rational ансатц через `sympy.rational_interpolate` + degree
  search + **независимая holdout-валидация** (≥`min_validation` точек); не прошло → `InterpolationFailed`.
  Требует ≥ `min_validation+2` точек → **никакой реконструкции по одной точке**.
- **reconstruct_coefficients(results, param_names)**: univariate (1 параметр) поддержан;
  **multivariate — консервативный refuse** (`InterpolationFailed`), не угадываем (Pass 2H).
- Работает **только над записями** `NormalFormResult` (family не нужна); нет hardcode имён/N/M.
  **Success не собирается** — это отдельный reducer-гейт (Pass 2H+).
*Обоснование:* spec §5.10/§7; method review §7 (Ошибка 5: не восстанавливать по одной точке).

## A27. Решения Pass 2G.1 (multi-sample normal-form records collector)
- **`records.py`** — новый producer-слой (архитектурно чище отдельного модуля, чем расширение
  `reconstruction.py`). `NormalFormRecord` — плоская сериализуемая точка `(prime, sample,
  target_label, status, formal_success, coeffs: dict[label,int mod p], support, all_terms_lf,
  non_lf_terms, unknown_lf_terms, rank, diagnostics)`. `coeffs` = переименованный `terms`
  (label → coeff mod prime); отсутствующий label = **точный ноль** в этой точке (не «выброшен»).
- **`collect_normal_form_records(family, rows, target, primes, samples, ...)`** — прогон
  `modular_normal_form` по `samples × primes` (samples внешний цикл, детерминированный порядок
  входа). **Каждая** точка записывается честно: `BadSpecialization`/`TargetNotReducible`/
  `EmptySystem` сохраняют реальный `status` (coeffs пустой), **не выбрасываются** — фильтрацию
  делает reconstruction (`collect_value_table` берёт только Reduced+formal_success, остальное —
  `n_skipped`). Никакого reconstruction/Success внутри collector.
- **Минимальная адаптация `reconstruction.py`** (не переписан): `collect_value_table` через
  `_record_coeffs(r)` принимает и `NormalFormResult` (`.terms`), и `NormalFormRecord` (`.coeffs`);
  union-support/0-fill/skip-логика без изменений. Univariate-only ограничение сохранено.
- **Не** реализовано в 2G.1: multivariate reconstruction (Pass 2H), reducer Success pipeline,
  D4 e2e. Recos остаются known-value-only (не LF-decomposition expectations); starred examples —
  см. A12–A15.
*Обоснование:* spec §5.10/§6/§7; строгое разделение producer(records)/consumer(reconstruction);
честная запись всех точек (bad не патчится, не теряется).

## A28. Решения Pass 2H (multivariate reconstruction/interpolation)
- **`interpolate_multivariate(values, params, max_deg=6, min_validation=2)`** — dense
  linear-algebra ансатц `C = N/D`: моном-базисы total-degree `≤ num_deg` (числитель) и
  `≤ den_deg` (знаменатель); в каждой sample-точке `s` уравнение `N(s) − v_s·D(s) = 0` →
  однородная **рациональная** матрица; решение = базис нуль-пространства (`sympy.Matrix.nullspace`,
  точные `Fraction`, setup-фаза — SymPy разрешён A5). `N/D` через `sympy.cancel`.
- **Degree search от простого к сложному** (по `num_deg+den_deg`, затем `num_deg`): возвращаем
  первую пару степеней, у которой nullspace **ровно 1-мерный** и кандидат **проходит независимую
  holdout-валидацию** (≥`min_validation` точек, `D≠0` на них). Требование dim==1 отсекает
  недоопределённые/неоднозначные степени; недостаточно точек (`< #num_mon+#den_mon−1`) → эта
  степень пропускается. Ничего не провалидировало в пределах `max_deg` → `InterpolationFailed`
  (**не угадываем**).
- **`reconstruct_coefficients` — диспетчер**: 1 параметр → `interpolate_univariate` (без изменений);
  ≥2 → `interpolate_multivariate`. Прежнее «multivariate → refuse» (A26) **снято**. Потребляет
  `NormalFormRecord` (collector, A27) через `collect_value_table`.
- **Ограничения (документированы):** dense, **не** sparse/Zippel — число неизвестных `~binom(k+d,d)`
  по num и den, экспоненциально по `k`/степени; практично для малых `k` (2–3) и умеренных степеней.
  Нужно достаточно различных sample-точек. Спец-точки, где зануляется знаменатель, приходят как
  `BadSpecialization` из collector и не попадают в values (skip+count, A25/A27).
- **Не** реализовано в 2H: reducer Success pipeline, `AllLocallyFinite` gate, D4 e2e (Pass 2I).
*Обоснование:* spec §5.10/§7; method review §7; строгая независимая валидация вместо угадывания.

## A29. Решения Pass D4.3 (rank-consistency record selection)
- **Только `Reduced`+`formal_success` records — coefficient records**; всё остальное никогда не
  попадает в value table (было и раньше, теперь оформлено в
  `reconstruction.select_records_for_reconstruction`).
- **Default `rank_policy="max_rank"`**: reconstruction потребляет только valid records с
  **максимальным наблюдённым RREF rank**. Математическое основание: ранг специализации может
  только *падать* ниже generic rank (миноры, тождественно равные нулю, нулю и остаются), поэтому
  max-rank records = generic points. Rank-deficient точка решает *другую* (меньшую) систему — её
  сжатый/сдвинутый support **нельзя** union-0-fill'ить (ложные нули ломают интерполяцию; симптом
  D4.2). Такие records skip+count, **не** патчатся.
- **Отсутствующий label в max-rank record — по-прежнему честный нуль** (special zero, union
  support); low-rank records нулей **не** вносят.
- **`rank_policy="all"`** — старое поведение (для тестов/отладки), явно, не по умолчанию.
- **Диагностика селекции обязательна**: `n_valid_records_before_rank_filter`, `selected_rank`,
  `n_selected_records`, `n_rank_filtered_records`, `rank_histogram`,
  `support_after_rank_filter` — в `reducer` попадает в
  `result.diagnostics.extra["record_selection"]`; `min_valid_records` проверяется по
  **post-filter** количеству.
*Обоснование:* D4.2 diagnosis (`notes/D4_STATUS.md`): generic rank 2041 vs rank-deficient сэмплы
(1995/2011) → несовместимые supports → InterpolationFailed. Generic-политика вместо
D4-хардкода; LF gate не ослаблен; `Success` — только через прежний строгий gate.

## A30. Решения Pass D4.4 (row-span certificate + non-degenerate sampling)
- **`certificate.verify_reduction_relation_mod_p(family, rows, target, terms, sample, prime)`** —
  generic модульный сертификат: собрать rows mod p → точно вычислить заявленные коэффициенты
  (SymPy/`ParamExpr`/`int`/`Fraction`) → редуцировать вектор `J[target] − Σ C_i·J[label_i]`
  пивотами RREF. `InSpan`/`NotInSpan` (+честный residual) / `BadSpecialization` (полюс в
  коэффициенте или в строке — точка отвергается, не патчится) / `EmptySystem`. **Никогда не
  стампит Success**; SymPy — только разовая оценка коэффициентов (не hot loop).
- **Урок (интерполяционная деградация на решётке):** product-grid из ≤ `max_deg+1` значений на
  ось **вырожден** для dense degree search — `Π(ep−a_i)` исчезает на всей решётке, включая
  holdout (он лежит на той же решётке) → «валидированный», но неверный вне решётки кандидат.
  Гейтовое `independent validation` ровно настолько независимо, насколько независимы точки
  holdout. Поэтому: (a) sample grids для multivariate reconstruction — **scattered/без
  product-структуры и без низкостепенной кривой через точки**; (b) финальная приёмка —
  **row-span certificate в off-sample rank-generic точках** (rank == generic; в rank-deficient
  точках сертификат reducer-output неинформативен).
- Проверка эквивалентности двух разложений разного базиса: оба соотношения in-span в одних и
  тех же точках ⟹ их разность in-span ⟹ эквивалентность **modulo generated rows**;
  покоэффициентное сравнение не требуется.
*Обоснование:* D4.4 (`notes/D4_STATUS.md`): сертификат поймал неверные off-grid коэффициенты
формально «успешного» D4.3-прогона на 6×6 integer lattice; honest-Failure/verification-first
принципы CLAUDE.md §5–6.

## A31. Решения Pass D4.5 (certificate gate в редьюсере)
- **`Success` по умолчанию требует row-span сертификации**: `ReducerConfig`
  (`certificate_points`, `require_certificate_for_success=True`, `min_certificate_points=1`);
  reducer сертифицирует реконструированное соотношение в независимых off-sample точках ДО gate.
  `certificate_status != "Passed"` при включённом флаге → `verification_failed=True` →
  `Error="VerificationFailed"` через существующий строгий gate (сам gate не менялся;
  `FormalSuccess` в Diagnostics остаётся честным).
- **Классификация точек** (Verify.1): rank **<** `selected_rank` реконструкции → rank-filtered
  (неинформативна — сертификат на другой страте ничего не доказывает); rank **>**
  `selected_rank` → **жёсткий провал** (`n_certificate_rank_exceeded`): ранг специализации не
  может превышать generic → значит `selected_rank` реконструкции был НЕ generic и её
  коэффициентам верить нельзя; полюс в строке/коэффициенте → bad (точка отвергается); иначе
  pass/fail по residual. Провал любой информативной точки или rank-exceeded → `Failed`; ноль
  информативных → `Insufficient` (никогда не «молчаливый pass»).
- **Конфиг (Verify.1)**: `certificate_primes` (по умолчанию — primes редукции),
  `certificate_rank_policy="selected_rank"` (единственная поддерживаемая политика, иное →
  `ValueError`).
- **Auto-точки** (когда `certificate_points` пуст и gate включён): детерминированные,
  строго за пер-параметровым максимумом сэмплов (гарантированно off-sample), с разными
  нечётными знаменателями (13/17/19), primes ротируются. Явные точки предпочтительны для
  дорогих кейсов (D4: probe-verified rank-generic).
- **Opt-out только явный**: `require_certificate_for_success=False` — старое поведение
  (для отладки/тестов); при переданных явных точках сертификат всё равно считается и
  репортится, просто не гейтит.
*Обоснование:* урок A30 должен жить в pipeline, а не только в тестах; regression
`tests/test_certificate_gate.py` воспроизводит product-grid false-success синтетически и
доказывает, что gate его ловит.

## A11. Формат вывода
Основной вывод — Wolfram-like текст (`^`, `/`, `<||>`, `->`, рациональные как `3/7`,
факторизованные коэффициенты). Это текстовый формат обмена, не исполняемый Wolfram-код.
*Обоснование:* CLAUDE.md, docs/00, spec §4, §5.11.
