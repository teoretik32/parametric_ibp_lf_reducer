# Risk Register — `parametric_ibp_lf_reducer`

Оценка: L=likelihood, I=impact (1..3). Живой документ.

## R1. Полнота набора лучей для LF в общем N — L2 I3
Ложный `True` для LF, если пропущен релевантный луч (внутренние нули, невыпуклые эффекты
при смешанных знаках коэффициентов).
**Митигация:** детерминированный набор (координатные `±e_i` + фасетные нормали Newton-политопов
`G_l` + нормали Минковски-суммы) **плюс** adaptive random-ray safety net. Любая найденная
дивергенция вне детерминированного набора → `LocallyFinite="Unknown"` (никогда `True`).
`Success` требует `True` у всех terms. Документировано в `assumptions.md` (§3).

## R2. Newton-полигон facets без scipy — L2 I2
Точные фасетные нормали в общем N требуют convex hull; scipy не в обязательных deps.
**Митигация:** scipy — опциональная зависимость (`ConvexHull.equations`). Fallback без scipy:
координатные лучи + нормали к рёбрам (разности exponent-векторов) + random net. Для
validation-кейсов (симплициальные политопы `1+x1+…`) этого достаточно.

## R3. D=4 end-to-end не сходится при узком label-box — L2 I3
Редукция может требовать более широкого диапазона n/m или большей степени IBP/tangent.
**Митигация:** adaptive search schedule (расширение n-range, m-range, `MaxIBPDegree`,
`TangentDegrees`, числа samples/primes). При невозможности — честный `Failure` с
`Recommendations` (какие параметры увеличить). Никакого фейкового Success.

## R4. Overfit в reconstruction — L2 I3
Rational-function ansatz может «подогнаться» под конечный набор точек.
**Митигация:** обязательная валидация на held-out independent samples (другие простые/точки).
Не прошло → `InterpolationFailed`. Degree search с запасом; проверка на нескольких primes.

## R5. Слишком строгий surface-check выбрасывает валидные coordinate rows — L2 I3
Прямо описанный в spec §7.1 / method review Ошибка (толкает к неверной или пустой системе).
**Митигация:** coordinate IBP `∂_{x_i}(P F)` проверяет **только** границы `x_i=0` и `x_i=∞`,
а не все toric-лучи. Toric flux — только для vector/tangent IBP. Есть регрессионный тест
на «не слишком строгий» (аналог 11.4, активируется при наличии family).

## R6. Ложные IBP-строки (игнор boundary terms) — L2 I3
Формальный total derivative без surface-проверки усиливает систему и даёт красивую, но
неверную редукцию (method review Ошибка 1).
**Митигация:** строка попадает в систему только после surface-free фильтра. Строки с
ненулевым boundary term не используются в MVP (можно логировать отдельно).

## R7. Производительность / взрыв выражений — L2 I2
Прямая работа над `Q(ep,r,…)` не масштабируется.
**Митигация:** modular/finite-field контур; mod-p int в hot loop; row templates строятся
один раз и переиспользуются между samples; кеш derivative/valuation/label-id; sparse rows;
dedup by hash. numba/parallel/streaming — full target.

## R8. Отсутствуют input families 11.2/11.4/11.5 — L3 I2
В docs только ожидаемые результаты, без входных полиномов.
**Митигация:** `xfail`-фикстуры «pending input family». MVP-приёмка на 11.1 + 11.3.
Активация при получении families от пользователя или корректной реконструкции.

## R9. Tangent syzygy solve тяжёлый символьно — L2 I2
Nullspace над `Q(params)` может расти.
**Митигация:** setup-only (вне hot loop), ограничен degree-блоками `TangentDegrees`.
При росте — mod-p ansatz на random param-sample + lift/переинстанцирование. Отбрасывание
нулевых/эквивалентных полей.

## R10. Смешение pointwise identity и integral identity — L1 I3
RHS — сумма **интегралов**, не поточечное равенство интеграндов (method review Ошибка 6).
**Митигация:** архитектура оперирует labels/интегралами; export формирует `Int[...] -> Σ C Int[...]`;
опциональный divergence certificate (full target) для строгой проверки.

## R11. Reconstruction по одной точке — L1 I3
Даёт лишь числовую редукцию (method review Ошибка 5).
**Митигация:** обязательный слой rational-function reconstruction; тест запрещает
single-point acceptance.

## R12. Неверный target label — L1 I2
Неоднозначность, что считать целевым интегралом.
**Митигация:** допущение — target = базовый интегранд `(n=0,m=0)` × `TargetMultiplier`
(`assumptions.md` §1); задокументировано и покрыто тестом.
