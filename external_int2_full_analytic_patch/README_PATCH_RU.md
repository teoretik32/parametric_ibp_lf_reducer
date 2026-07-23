# External Int2 full analytic oracle patch

Этот patch добавляет в проект независимый аналитический oracle для External Int2:

- точное интегрирование по `x7`;
- ODE/epsilon-рекурсию для `S=rQ`;
- компактный Laurent-ряд через `ep^0`;
- быстрые символические тесты;
- machine-readable validation JSON;
- сохранённый source notebook.

Oracle **не является LF-разложением**. Он нужен как точный acceptance reference для дальнейшего исправления/расширения parametric-IBP реализации.

## Установка

Распаковать ZIP поверх корня репозитория, затем:

```powershell
python scripts/audit_external_int2_full_laurent.py --numeric
python -m pytest tests/test_external_int2_full_laurent.py
ruff check .
```

Далее выполнить последовательно:

1. `PROMPT_1_INTEGRATE_ANALYTIC_ORACLE_RU.md`
2. после green/commit — `PROMPT_2_FIND_TRUE_LF_BASIS_RU.md`

## Важное ограничение

Текущий snapshot проекта, доступный при подготовке patch, может быть старее вашей рабочей ветки Method.8+. Patch специально состоит в основном из новых файлов. Fable должен адаптировать импорты/документацию к фактическому HEAD и не откатывать более новые изменения.
