# Scenario 1: from scratch

Готовый репозиторий-заготовка для запуска Claude Code над задачей `parametric_ibp_lf_reducer`.

Главное: это **Python-only** проект. Wolfram/Mathematica syntax встречается только как текстовый формат входа/выхода и как нотация в математических примерах.

## С чего начать Claude Code

1. Прочитать `CLAUDE.md`.
2. Прочитать `TASK.md`.
3. Прочитать `docs/00_READ_ME_FIRST_RU.md`, затем cleaned spec/method/validation.
4. Реализовать Python-пакет в `src/parametric_ibp_lf_reducer/`.
5. Добавить pytest tests в `tests/` и прогонять их регулярно.

## Рекомендуемые команды

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -e '.[dev,speed]'
python -m pytest
```

## Важное ограничение

Если алгоритм пока нашёл только formal normal form, но финальные terms не локально конечны при `epsilon=0`, он обязан вернуть Failure/partial diagnostics, а не Success.

## Usage (Release.1)

Полный пользовательский путь (CLI/API, tiny + D4 heavy examples, exit codes,
diagnostics JSON, limitations): **`docs/USAGE.md`**.

Быстрый старт:

```bash
python -m parametric_ibp_lf_reducer reduce examples/tiny_success_input.wl.txt \
    --out result.m --diagnostics-json diagnostics.json
```

Release sanity: `scripts/final_check.sh` (или `scripts/final_check.ps1`;
флаг `--heavy` / `-Heavy` добавляет ~25-30 мин D4 acceptance-прогоны).
