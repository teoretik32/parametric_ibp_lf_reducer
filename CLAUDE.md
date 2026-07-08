# CLAUDE.md

## Основная роль

Ты Claude Code в этом репозитории. Твоя задача — реализовать Python-пакет `parametric_ibp_lf_reducer`, который делает parametric IBP/algebraic reduction до локально конечных интеграндов.

## Самое важное

1. Писать только Python-код.
2. Wolfram/Mathematica-style syntax в документах — это только текстовый формат ввода/вывода и нотация примеров.
3. Не вызывать Mathematica/Wolfram как обязательный backend.
4. Не хардкодить `N=3`, `N=4`, `x4`, `x7`, `H`, `A`, `B`, `C` или любой validation case.
5. `Success` разрешён только если `AllLocallyFinite=True`, все финальные terms локально конечны при `epsilon=0`, surface checks пройдены, reconstruction проверен на независимых samples.
6. Формальная normal form с нелокально конечными terms — это Failure/partial diagnostics, не success.
7. SymPy можно использовать для парсинга, финальной факторизации и reconstruction, но не в массовом row-generation/RREF hot loop.
8. Основной вычислительный контур должен быть sparse/modular/finite-field.

## Порядок чтения

1. `TASK.md`
2. `docs/00_READ_ME_FIRST_RU.md`
3. `docs/01_clean_spec_ru.md`
4. `docs/02_method_review_ru.md`
5. `docs/03_validation_cases_ru.md`

## Рабочий протокол

- Сначала создай/обнови `notes/implementation_plan.md` с коротким планом модулей и порядком тестов.
- Делай small commits/steps логически: parser → sparse polynomials → valuations/surface → rows → modular linear algebra → reconstruction → export/CLI.
- После каждого слоя добавляй тесты.
- Запускай `python -m pytest`; если тест тяжёлый, добавь lightweight unit tests и отдельный integration marker.
- Если capability ещё не реализован, возвращай честный Failure/status, не фальсифицируй Success.
- Для сложных/больших cases используй adaptive diagnostics и recommendations.

## Public API target

Preferred API:

```python
from parametric_ibp_lf_reducer import reduce_wolfram_style_input

result = reduce_wolfram_style_input(input_text, max_ibp_degree=2, tangent_degrees=[(1, 1), (2, 2)])
print(result.wolfram_style_text)
```


## CLI target

```bash
python -m parametric_ibp_lf_reducer reduce input.m --out result.m --log result.log --diagnostics-json diagnostics.json
```

`result.m` — это текстовый Wolfram-like output, а не кодовая база на Mathematica.


## Сценарий этого архива

Пиши реализацию с нуля. В этой папке намеренно нет старого нерабочего Python-кода. Опираться нужно на cleaned spec, method review и validation cases.
