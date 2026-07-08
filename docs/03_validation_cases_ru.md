# Validation cases and acceptance checks

# Уточнение для Claude Code: Python-only, Mathematica/Wolfram только как текстовая нотация

В этой папке нужно писать **только Python-пакет**. Все упоминания Mathematica, Mathematica-style, Wolfram-style, `^`, `{}`, `<| |>`, `->`, `Int[...]` и файлов вроде `result.m` означают **только текстовый формат обмена и нотацию примеров**.

Запрещено:

- писать реализацию на Wolfram Language / Mathematica;
- создавать notebook как основной результат;
- требовать установленную Mathematica как runtime;
- подменять Python-алгоритм вызовами Mathematica;
- воспринимать математические формулы ниже как исполняемый код.

Разрешено и нужно:

- реализовать парсер/экспорт такого синтаксиса на Python;
- использовать SymPy только для парсинга, финальной факторизации и reconstruction, но не в hot loop;
- экспортировать результат строкой в Wolfram-like/Mathematica-style синтаксисе, если это нужно для совместимости downstream-notebook.



Эти кейсы вынесены отдельно, чтобы Claude Code не искал их по длинному ТЗ. Они должны стать `pytest`/fixture tests и независимыми validation samples.

## Input contract excerpt

3. Формат входа
---------------

Код должен поддерживать два режима входа.

3.1. Явный Wolfram-like/Mathematica-style текстовый формат семейства — обязательный режим

Пример:

    IBPInput = <|
      "Variables" -> {x1,x2,x3,x4},
      "Parameters" -> {ep,r},
      "Regulators" -> {ep},
      "Domain" -> "PositiveOrthant",

      "Polynomials" -> <|
        "G0" -> 1 + x1 + x2 + x3,
        "G1" -> 1 + x4,
        "G2" -> r*x1*x2 + r*x2*x3 + r*x1*x2*x4 + x3*x4
      |>,

      "MonomialExponents" -> <|
        x1 -> 0,
        x2 -> 0,
        x3 -> 0,
        x4 -> -ep
      |>,

      "PolynomialExponents" -> <|
        "G0" -> 2*ep,
        "G1" -> 1 + 3*ep,
        "G2" -> -2 - ep
      |>,

      "TargetMultiplier" -> 1,
      "Assumptions" -> {r > 0},

      "Options" -> <|
        "MaxIBPDegree" -> 2,
        "TangentDegrees" -> {{1,1},{2,2}},
        "LabelBox" -> Automatic,
        "SurfaceMode" -> "CoordinateForCoordinateIBP_ToricForVectorIBP",
        "OutputFormat" -> "WolframStyleText"
      |>
    |>

Здесь исходный интегранд понимается как

    Product[x_i^MonomialExponents[x_i]] * Product[G_l^PolynomialExponents[G_l]] * TargetMultiplier.

3.2. Авторазбор цельного Wolfram-like/Mathematica-style текстового интегранда — желательный frontend

Пример:

    IBPInput = <|
      "Variables" -> {x1,x2,x3,x4},
      "Parameters" -> {ep,r},
      "Regulators" -> {ep},
      "Integrand" -> x4^(-ep)*(1+x1+x2+x3)^(2*ep)*(1+x4)^(1+3*ep)
                     *(r*x1*x2+r*x2*x3+r*x1*x2*x4+x3*x4)^(-2-ep),
      "Domain" -> "PositiveOrthant",
      "Assumptions" -> {r > 0}
    |>

Frontend должен попытаться разложить Integrand на:

    - мономиальные степени x_i;
    - степени явно встречающихся полиномиальных факторов;
    - общий рациональный множитель.

Если факторизация неоднозначна, frontend должен вернуть ParserNeedsExplicitFamily и предложить пользователю явный формат 3.1. Ядро редуктора не должно зависеть от эвристического авторазбора.

## Regression cases

11. Регрессионные тесты
-----------------------

Тесты должны быть оформлены как pytest + сохранённые expected-output файлы.

11.1. Минимальный N=2 tangent-field sanity test

G = 1+x+y.
Проверить, что поле Q=(xy,-xy) даёт строку без m-shift.

11.2. Старый N=3 пример с пятёркой локально конечных интегралов

Ожидаемый локально конечный basis в старых обозначениях:

    1/A,
    1/(A*B),
    1/(A*C),
    1/(B*C),
    1/(B*C^2).

Тест должен проверять:

    - target reducible;
    - all_lf=True;
    - coefficients match expected symbolic or modular samples.

11.3. D=4 пример с рациональными коэффициентами от epsilon и r

Семейство:

    F = x4^(-ep) * G0^(2 ep) * G1^(1+3 ep) * H^(-2-ep),

    G0 = 1+x1+x2+x3,
    G1 = 1+x4,
    H  = r*x1*x2 + r*x2*x3 + r*x1*x2*x4 + x3*x4.

Ожидаемый LF basis:

    M1 factor = x2*x3/(G0^2*G1),
    M2 factor = x1*x2/(G0^2*G1),
    M3 factor = x2*x3/(G0^3*G1),
    M4 factor = x1*x2/(G0^3*G1),
    M5 factor = x2*x3/(G0^4*G1).

Ожидаемые коэффициенты:

    C1 = (2*ep - 1)*(6*ep^2*r + 3*ep^2 - 12*ep*r - ep + 4*r)/ep^3

    C2 = 2*(2*ep - 1)*(3*ep*r + 2*ep - r)/ep^2

    C3 = 6*(ep - 1)^2*(2*ep - 1)*(2*ep*r - ep - 4*r)/ep^3

    C4 = 4*(ep - 1)*(2*ep - 1)*(3*ep*r - r + 3)/ep^2

    C5 = -2*(ep - 1)*(2*ep - 3)*(2*ep - 1)*(3*ep - 4)*(r + 1)/ep^3

Проверить независимые finite-field samples и одну рациональную точку, например ep=-3/4, r=7/5.

11.4. ID4example2 coordinate-surface regression

Проверяет, что coordinate IBP surface-check не является чрезмерно строгим.
Ожидаемый one-term LF result:

    I = C(ep) * M,

где

    M factor = x5*x7/(G1^3*G2),

    C(ep) = -3*(6*ep-1)*(6*ep+1)*(9*ep-1)*(34*ep-1)/(112000*ep^4*(8*ep+1)).

11.5. ID3example3 с x4^2 в полиноме

Проверяет поддержку мономов степени выше 1.
Базис в проверенных точках должен включать:

    x4*x7/((1+x4)*(1+x7)),
    x4*x7/((1+x7)*G2),
    x4*x7*x8/G3.

## CLI/API contract excerpt

12. Минимальная CLI-форма
-------------------------

Команда:

    python -m parametric_ibp_lf_reducer reduce input.m --out result.m --log result.log

Опции:

    --max-ibp-degree 2
    --n-range auto
    --m-range auto
    --tangent-degrees 1,1 2,2
    --prime-count 8
    --sample-count 64
    --epsilon-direction minus
    --jobs 16
    --memory-limit-gb 64
    --timeout-sec 0
    --preferred-masters preferred.m
    --diagnostics-json diagnostics.json

Поведение:

    - result.m содержит только Wolfram-like текстовую result association;
    - result.log содержит читаемый лог;
    - diagnostics.json содержит машинно-читаемую диагностику;
    - exit code 0 только при Status -> Success.


13. Python API
--------------

Минимальный API:

    from parametric_ibp_lf_reducer import reduce_wolfram_style_input

    result = reduce_wolfram_style_input(
        input_text,
        max_ibp_degree=2,
        tangent_degrees=[(1,1),(2,2)],
        jobs=8,
        output_format="mathematica"
    )

    print(result.wolfram_style_text)

Объекты:

    ReductionResult:
      status: str
      terms: list[ReductionTerm]
      diagnostics: Diagnostics
      wolfram_style_text: str

    ReductionTerm:
      coeff: RationalFunction
      label: tuple[int,...]
      integrand: Expr
      locally_finite: bool

## Дополнительные правила для тестов

- Тесты не должны проходить за счёт hardcode конкретных имён `x4`, `x7`, `H`, `A`, `B`, `C` в ядре.
- Если реализация пока не умеет получить редукцию для тяжёлого кейса, она должна честно вернуть Failure-статус, а не Success с `AllLocallyFinite -> False`.
- D=4 coefficients должны проверяться как символьные выражения или через независимые finite-field / rational samples.
- Локальная конечность проверяется при `epsilon = 0`, а не за счёт знака регулятора.
