# parametric_ibp_lf_reducer

*Read this in English: [README.md](README.md).*

Чистый Python-редьюсер для **параметрических IBP-систем** (integration by parts),
который строит редукции, где **каждый терм локально конечен при `epsilon = 0`**,
а каждый `Success` дополнительно проходит независимый точный модулярный
**row-span сертификат**.

Wolfram/Mathematica синтаксис используется **только как текстовый формат обмена**
(входные/выходные документы). Mathematica/Wolfram runtime не требуется и никогда
не вызывается.

## Что значит «локально конечен при epsilon = 0»

Терм редукции `C(params, ep) * Integral[F]` *локально конечен*, если подынтегральное
выражение `F` даёт интеграл, чьё параметрическое представление остаётся интегрируемым
потермово при `ep -> 0` (внутри отдельного терма не прячется расходимость `1/ep`).
Если алгоритм нашёл только *формальную* нормальную форму, чьи термы не все локально
конечны, он обязан вернуть типизированный `Failure`
(`NormalFormNotLocallyFinite`), а не «тихий» Success.

## Требование к входу: explicit family

Входной документ обязан содержать **явную параметрическую семью** (пропагаторы /
определение семьи). Документ, где есть только сырой `Integrand`, не «угадывается»:
запуск возвращает типизированный `Failure` с причиной `ParserNeedsExplicitFamily`.

## Установка

```bash
python -m venv .venv
source .venv/bin/activate            # PowerShell: .venv\Scripts\Activate.ps1
pip install -e .                     # базовая установка (dict RREF-бэкенд)
pip install -e ".[speed]"            # + опциональный Numba RREF-бэкенд
python -m pip install -e '.[dev,speed]'   # для разработки: dev-инструменты + speed
python -m pytest                     # быстрый набор тестов
ruff check .
```

Требуется Python >= 3.11. `speed` (numba) — опционально; без него пакет
полностью работает на дефолтном `dict`-бэкенде.

## Быстрый старт: CLI

Tiny-пример (~1–2 с, полный конвейер, включая сертификат):

```bash
python -m parametric_ibp_lf_reducer reduce examples/tiny_success_input.wl.txt \
    --out result.m --diagnostics-json diagnostics.json
```

Ожидается: exit code 0, `"Status" -> "Success"`, один локально конечный терм.

Тяжёлый пример D4 (**медленно, ~10–15 минут**, сертифицированный релизный baseline):

```bash
python -m parametric_ibp_lf_reducer reduce examples/d4_cli_example_input.wl.txt \
    --out d4_result.m --diagnostics-json d4_diagnostics.json
```

Ожидается: `Success`, `AllLocallyFinite -> True`, сертификат `Passed`,
**3-термный** LF-базис {M1, M2, M3} (см. «Ограничения»).

## Быстрый старт: Python API

```python
from parametric_ibp_lf_reducer import api

result = api.reduce_wolfram_style_input(
    open("examples/tiny_success_input.wl.txt").read()
)
print(result.status)              # "Success"
print(result.wolfram_style_text) # Wolfram-подобный документ результата
```

## Опциональный Numba RREF-бэкенд (v0.1.4)

Ядро модулярного RREF доминирует в тяжёлых прогонах. v0.1.4 добавляет
опциональный, автоматически выбираемый Numba-бэкенд:

```bash
python -m parametric_ibp_lf_reducer reduce input.wl.txt --rref-backend auto
python -m parametric_ibp_lf_reducer reduce input.wl.txt --rref-backend numba_int_array_experimental
```

- **Дефолт по-прежнему `dict`** — без явного opt-in ничего не меняется.
- **`auto` — рекомендуемый opt-in для больших систем**, когда установлен
  extra `[speed]`; Numba выбирается по-матрично, только если система
  достаточно велика (пороги: `min_rows=500`, `min_cols=400`, `min_nnz=3000`)
  и `prime < 2^31`. Малые системы обычно остаются на `dict`.
- `auto` **откатывается на `dict`**, если Numba недоступна; *явный* запрос
  Numba, наоборот, **падает с понятной ошибкой**.
- Сертифицированный бенчмарк (corrected Example 4\*, полный бокс, 972 лейбла,
  12360 строк): wall **3963.4с (dict) → 766.5с (auto → Numba), ~5.17×**,
  математические результаты побитово идентичны, сертификат `Passed`. См.
  [docs/PERFORMANCE.md](docs/PERFORMANCE.md) и
  [docs/NUMBA_RREF_QA.md](docs/NUMBA_RREF_QA.md).

## Статусы

| статус | смысл |
|--------|-------|
| `Success` | все термы локально конечны **и** пройден certificate gate |
| `Failure` (грубый экспортируемый статус; конкретная причина — в `Error` / JSON `status`) | честный типизированный отказ: `TargetNotReducible`, `InterpolationFailed`, `NormalFormNotLocallyFinite`, `ResourceLimitReached` |
| `VerificationFailed` | редукция найдена формально, но независимая проверка (reconstruction check / row-span сертификат) её отвергла |
| `ParserNeedsExplicitFamily` | во входном документе нет явной параметрической семьи; ничего не редуцировалось |

Exit-коды: `0` = Success, `1` = типизированный отказ (result и JSON всё равно
записываются), `2` = usage / битый документ. Поля diagnostics JSON описаны в
[docs/USAGE.ru.md](docs/USAGE.ru.md).

## Ограничения (Release.1)

- **Только explicit-family вход** — авто-факторизация integrand не выполняется.
- **Нет adaptive search** — один фиксированный проход на вызов; при отказе ничего
  не расширяет label box / степени / сэмплы автоматически.
- **Пределы плотной многомерной реконструкции** — высокая степень коэффициентов или
  много параметров требуют больше scattered-сэмплов, чем в дефолтах.
- **D4 редуцируется к 3-термному LF-базису** {M1, M2, M3}, эквивалентному
  (и сертифицированному против row span, содержащего) 5-термному референсному
  базису M1..M5; референсный базис сознательно не форсируется.
- **Нет зависимости от Mathematica runtime** — и потому нет символьного кросс-чека
  через Wolfram; верификация — модулярная row-span сертификация в рациональных точках.

## Документация

- [docs/USAGE.ru.md](docs/USAGE.ru.md) — полный пользовательский путь (рус.);
  [docs/USAGE.md](docs/USAGE.md) — in English.
- [docs/FINAL_QA.md](docs/FINAL_QA.md) — релизный QA; [CHANGELOG.md](CHANGELOG.md).
- `docs/0*_ru.md` — исходная спецификация / метод / валидационные заметки.

## Лицензия

Лицензия пока не выбрана — см. [docs/LICENSE_NOTE.md](docs/LICENSE_NOTE.md)
перед любым переиспользованием кода.
