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

Техническое задание на Python-код: универсальный редуктор IBP в параметрическом представлении до локально конечных интеграндов
=============================================================================================================================

Версия: 2026-07-06
Рабочее название: parametric_ibp_lf_reducer
Целевой пользователь: исследователь, работающий с параметрическими представлениями многократных интегралов Фейнмана.


1. Цель
-------

Нужно разработать Python-код, который принимает на вход параметрический интеграл в Wolfram-like/Mathematica-style текстовом формате, строит IBP/algebraic редукцию в параметрическом представлении и выдаёт представление исходного интеграла через набор локально конечных интегралов.

Вход:

    - Wolfram-like/Mathematica-style текстовое описание переменных интегрирования x_i;
    - список внешних параметров: epsilon, кинематические переменные, массы, ratios и т.п.;
    - исходный интегранд или явно разобранное семейство полиномов G_l и степеней;
    - область интегрирования: по умолчанию R_+^N;
    - целевая метка интеграла;
    - настройки поиска редукции.

Выход:

    - Wolfram-like текстовый экспорт с коэффициентами C_a(parameters);
    - Wolfram-like текстовый список локально конечных подынтегральных выражений F_a(x,parameters);
    - правило уровня интегралов

          Integral[OriginalIntegrand] -> Sum_a C_a Integral[F_a];

    - диагностический блок: labels, local-finiteness flags, surface checks, samples, primes, статус, ошибки.

Принципиальные требования:

    1. Количество переменных интегрирования N не должно быть зашито.
    2. Количество полиномов G_l не должно быть зашито.
    3. Количество внешних параметров не должно быть зашито.
    4. Коэффициенты должны восстанавливаться как рациональные функции внешних параметров, а не вычисляться только в одной числовой точке.
    5. При success=True все интегранды в RHS обязаны быть локально конечными при epsilon=0.
    6. Код не должен выдавать formal normal form как успешную физическую редукцию, если all_lf=False.
    7. Все алгоритмы должны быть рассчитаны на большие промежуточные выражения и большие разреженные системы.


2. Математическая форма семейства
---------------------------------

Базовая форма семейства:

    J(n,m) = ∫_{R_+^N} d^N x
             prod_i x_i^(a_i(parameters) + n_i)
             prod_l G_l(x,parameters)^(b_l(parameters) + m_l).

Или, эквивалентно,

    F_label = F_base * prod_i x_i^n_i * prod_l G_l^m_l,

где

    label = (n_1,...,n_N,m_1,...,m_M).

Внутреннее представление должно поддерживать:

    - N ≥ 1 произвольное;
    - M ≥ 0 произвольное;
    - sparse-полиномы G_l с мономами любой степени, например x4^2, x1*x2*x4, x7^3 и т.п.;
    - коэффициенты полиномов G_l как рациональные функции внешних параметров;
    - степени x_i и G_l как рациональные/аффинные выражения от регуляторов и внешних параметров;
    - несколько регуляторов в перспективе, но MVP может быть оптимизирован под один epsilon с архитектурной возможностью расширения.

Ограничение MVP, которое допустимо явно зафиксировать:

    - область интегрирования R_+^N;
    - степени зависят от epsilon линейно или рационально-просто;
    - полиномы G_l имеют коэффициенты из Q(parameters);
    - bulk singularities вне boundary не обрабатываются автоматически, если пользователь не задал assumptions, гарантирующие отсутствие внутренних нулей.


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


4. Формат выхода
----------------

Выход должен быть обычным текстом Wolfram-like export-представления. Основной объект:

    IBPReductionResult = <|
      "Status" -> "Success",
      "AllLocallyFinite" -> True,
      "Variables" -> {x1,x2,x3,x4},
      "Parameters" -> {ep,r},
      "Assumptions" -> {r > 0},

      "OriginalIntegrand" -> (...),

      "Terms" -> {
        <|
          "Coeff" -> C1[ep,r] in explicit Wolfram-like text syntax,
          "Integrand" -> localIntegrand1,
          "Label" -> {n1,n2,n3,n4,m1,m2,m3},
          "LocallyFinite" -> True
        |>,
        ...
      },

      "IntegralRule" -> HoldForm[
        Int[OriginalIntegrand,{x1,0,Infinity},{x2,0,Infinity},...] ->
          C1 Int[localIntegrand1,{x1,0,Infinity},...] +
          C2 Int[localIntegrand2,{x1,0,Infinity},...] + ...
      ],

      "Diagnostics" -> <|
        "NVariables" -> N,
        "NPolynomials" -> M,
        "NLabels" -> ...,
        "NRows" -> ...,
        "Rank" -> ...,
        "Primes" -> {...},
        "ParameterSamples" -> {...},
        "SurfaceRowsRejected" -> ...,
        "NonLFLabelsEliminated" -> ...,
        "InterpolationVerified" -> True,
        "IndependentChecks" -> {...}
      |>
    |>;

Для failure:

    IBPReductionResult = <|
      "Status" -> "Failure",
      "Error" -> "NormalFormNotLocallyFinite" | "TargetNotReducible" | "InterpolationFailed" | ...,
      "FormalSuccess" -> True|False,
      "AllLocallyFinite" -> False,
      "PartialNormalForm" -> {...},
      "Recommendations" -> {...}
    |>;

Запрещается возвращать Status -> "Success", если AllLocallyFinite -> False.


5. Основные модули Python-пакета
--------------------------------

Рекомендуемая структура:

    parametric_ibp_lf_reducer/
      __init__.py
      input_parser.py
      wolfram_text_export.py
      sparse_poly.py
      family.py
      valuations.py
      surface.py
      row_generation.py
      tangent_fields.py
      ranking.py
      finite_field.py
      sparse_rref.py
      interpolation.py
      reconstruction.py
      reducer.py
      diagnostics.py
      tests/

5.1. input_parser.py

Функции:

    parse_mathematica_association(text: str) -> RawInput
    parse_explicit_family(raw: RawInput) -> ParametricFamily
    try_factor_integrand(raw: RawInput) -> ParametricFamily | ParserFailure

Требования:

    - поддержка операторов Wolfram-like text syntax ^, *, +, -, /;
    - поддержка списков { ... } и Association <| ... |>;
    - сохранение имён переменных и параметров;
    - отсутствие вычислительной алгебры в hot path;
    - все выражения коэффициентов переводятся в sympy AST или внутренний рациональный AST.

5.2. sparse_poly.py

Внутренний формат полинома:

    SparsePoly = dict[ExponentTuple, CoeffExpr]

где

    ExponentTuple = tuple[int, ...] length N.

Функции:

    add, mul, pow_by_small_int, derivative(var_index), monomial_mul,
    evaluate_mod_prime(sample, prime),
    support(), degree(), valuation(ray).

Требования:

    - никакого expanded SymPy в hot loop;
    - мономы хранятся отсортированно/канонически;
    - коэффициенты параметрические, но в finite-field фазе быстро специализируются в int mod p;
    - поддержка степеней x_i^2, x_i^k без специальных случаев.

5.3. family.py

Класс:

    ParametricFamily:
      variables: tuple[str]
      parameters: tuple[str]
      regulators: tuple[str]
      polynomials: tuple[SparsePoly]
      monomial_exponents: tuple[CoeffExpr]
      polynomial_exponents: tuple[CoeffExpr]
      base_integrand_expr: original symbolic expression
      assumptions: list

Методы:

    label_to_factor(label) -> sparse/rational factor
    label_to_mathematica_integrand(label) -> str
    exponent_at_label(label) -> exponents of x_i and G_l
    specialize(sample, prime) -> FamilyModP

5.4. valuations.py

Назначение: проверка локальной конечности и surface behavior.

Функции:

    compute_candidate_rays(family) -> list[Ray]
    valuation_poly(poly, ray) -> int
    scaling_score_integrand(label, ray, eps_value=0) -> Score
    is_locally_finite(label) -> bool

Определение локальной конечности:

    is_locally_finite(label) == True

только если для всех rays выполнено

    base_score(label, ray, epsilon=0) > 0.

Нельзя учитывать знак epsilon-регулятора как локальную конечность.

Для вычисления rays:

    - coordinate zero rays;
    - infinity/compactification rays;
    - нормали граней Newton polytopes для G_l;
    - смешанные rays из common refinement нормальных фанов;
    - fallback: adaptive random/integer rays для поиска пропущенных дивергенций.

5.5. surface.py

Функции:

    coordinate_primitive_surface_free(family, label, var_index, multiplier) -> bool
    vector_field_surface_free(family, label, vector_field) -> bool

Правила:

    - для coordinate IBP ∂_{x_i}(Q_i F) проверять только границы компоненты x_i=0 и x_i=∞;
    - для vector/tangent IBP div(QF) проверять нормальный flux по toric rays;
    - surface-free тест может использовать регулируемую область epsilon -> 0^- или epsilon -> 0^+, заданную в Options;
    - локальная конечность финальных masters проверяется отдельно и строго при epsilon=0.

5.6. row_generation.py

Генерировать строки трёх типов.

A. Coordinate IBP rows:

    0 = ∫ ∂_{x_i}( P(x) F_label ).

P(x) — моном или малый sparse-полином из заданного ansatz. После раскрытия:

    ∂_{x_i}(P F) = (∂P) F + P F [a_i/x_i + sum_l b_l ∂_iG_l/G_l].

Термы переводятся в labels:

    monomial shift -> n-shift,
    division by G_l -> m_l - 1,
    multiplication by ∂G_l monomials -> n-shifts.

Строка добавляется только если primitive прошёл surface-free фильтр.

B. Algebraic rows:

    J(n,m) - sum_a c_{l,a} J(n+a,m-e_l) = 0.

C. Tangent/syzygy IBP rows:

    0 = ∫ div(Q F_label),
    Q · ∂G_l = H_l G_l.

Если Q касателен ко всем G_l, строка не сдвигает m. Если касателен к подмножеству, строка может использоваться в hybrid mode с явным контролем m-shifts.

5.7. tangent_fields.py

Два режима:

    1. External CAS mode:
       вызвать Singular/Sage для syzygy/intersection modules.

    2. Finite-field ansatz mode:
       для degree bounds d_Q, d_H построить линейную систему по коэффициентам Q_i и H_l:

           sum_i Q_i ∂_iG_l - H_l G_l = 0.

       Решить над GF(p), получить basis vector fields, поднять или использовать как modular row generators.

Требования:

    - не зашивать число переменных;
    - кешировать derivative tables;
    - поддерживать несколько degree-блоков;
    - отбрасывать нулевые/эквивалентные поля;
    - проверять surface-free flux перед добавлением строк.

5.8. ranking.py

Ранжирование переменных для исключения.

Priority tiers:

    Tier 0: target label, если нужно получить его normal form.
    Tier 1: labels, которые не локально конечны при epsilon=0.
    Tier 2: локально конечные, но сложные labels.
    Tier 3: preferred local finite masters, если пользователь их задал.
    Tier 4: самые простые локально конечные labels, которые можно оставить свободными.

Сложность label:

    complexity(label) = weighted sum of
      total positive n-shift,
      total negative m-depth,
      numerator total degree,
      denominator polynomial count,
      support size of factor,
      expected HyperInt difficulty proxy,
      user penalties.

Нельзя оставлять нелокально конечный label свободным только потому, что у него маленькая сумма индексов.

5.9. finite_field.py и sparse_rref.py

Требования:

    - арифметика GF(p) без SymPy в hot loop;
    - строки в sparse dict или CSR-like формате;
    - row normalization по pivot;
    - modular inverse через pow(a, p-2, p) или batch inversion;
    - возможность forward elimination и backward substitution;
    - сохранение pivot map для normal form target;
    - повторное использование структуры матрицы между parameter samples;
    - потоковая генерация строк для больших систем;
    - row hashing и удаление дублей.

Рекомендуемые оптимизации:

    - numba для hot loops;
    - хранение exponent tuples как compact int arrays;
    - block/sector ordering;
    - parallel row generation;
    - parallel finite-field samples;
    - optional C++/Rust backend для RREF при больших системах.

5.10. interpolation.py и reconstruction.py

Задача: восстановить C_a(parameters) из finite-field normal forms.

Алгоритмы:

    - CRT по нескольким простым p;
    - rational reconstruction для Q;
    - univariate Padé/Thiele для одного параметра;
    - multivariate Zippel sparse interpolation для многих параметров;
    - adaptive numerator/denominator degree search;
    - validation on independent samples;
    - factorized output.

Интерфейс:

    interpolate_normal_form_coefficients(records, param_symbols) -> dict[label, RationalFunction]

Где records содержат:

    prime,
    sample point,
    label -> coeff_mod_prime.

Требования:

    - если support normal form меняется между sample points из-за специальных нулей коэффициентов, нужно брать union support и считать отсутствующий коэффициент равным нулю в данной точке;
    - специальные точки, зануляющие знаменатели, отбрасывать;
    - если interpolation не прошла независимую проверку, вернуть InterpolationFailed.

5.11. wolfram_text_export.py

Функции:

    coeff_to_wolfram_text(expr) -> str
    integrand_to_wolfram_text(label) -> str
    result_to_wolfram_text(result) -> str

Требования:

    - использовать Wolfram-like текстовый синтаксис: ^, *, /, { }, <| |>, ->;
    - в Wolfram-like текстовом выходе не использовать Python-синтаксис **;
    - рациональные числа как 3/7, а не 0.428571;
    - факторизованные коэффициенты, например

          -3*(6*ep-1)*(6*ep+1)*(9*ep-1)*(34*ep-1)/(112000*ep^4*(8*ep+1))

    - интегранды выводить в форме исходный F_base * multiplier или полностью раскрытой factorized-формой;
    - labels выводить как списки целых чисел.


6. Главный алгоритм reducer.py
------------------------------

Псевдокод:

    def reduce_to_locally_finite(input_text):
        raw = parse_mathematica_association(input_text)
        family = parse_explicit_family(raw) or try_factor_integrand(raw)

        prepare_sparse_polynomials(family)
        rays = compute_candidate_rays(family)
        derivative_cache = build_derivative_cache(family)
        valuation_cache = build_valuation_cache(family, rays)

        for search_level in adaptive_search_schedule(raw.Options):
            labels = enumerate_labels(family, target, search_level)
            lf_flags = {label: is_locally_finite(label) for label in labels}

            row_templates = []
            row_templates += generate_algebraic_rows(labels, family)
            row_templates += generate_coordinate_ibp_rows(labels, family, search_level, surface_filter=True)

            if search_level.use_tangent_fields:
                tangent_fields = generate_tangent_fields(family, search_level.tangent_degrees)
                row_templates += generate_tangent_ibp_rows(labels, family, tangent_fields, surface_filter=True)

            ranking = build_ranking(labels, target, lf_flags, preferred_masters)

            records = []
            for prime in primes:
                for sample in parameter_samples:
                    if bad_denominator(sample, prime):
                        continue
                    matrix = instantiate_rows_mod_p(row_templates, sample, prime)
                    nf = sparse_rref_normal_form(matrix, ranking, target)
                    if nf.formal_success:
                        records.append((prime, sample, nf))

            if not enough_records(records):
                continue

            coeffs = interpolate_normal_form_coefficients(records, family.parameters)
            if not coeffs.verified:
                continue

            terms = build_terms(coeffs, family)
            all_lf = all(is_locally_finite(term.label) for term in terms)
            if all_lf:
                return export_success(family, terms, diagnostics)
            else:
                save_partial_failure(...)
                continue

        return export_failure(best_partial_result)

Адаптивный schedule должен расширять:

    - диапазоны n_i;
    - диапазоны m_l;
    - степень monomial multipliers P(x);
    - степень tangent-field ansatz;
    - число parameter samples;
    - число primes;
    - density of preferred LF candidate basis.


7. Surface-free проверки: точные требования
-------------------------------------------

7.1. Coordinate IBP

Для строки

    ∂_{x_i}(P F_label)

primitive равен

    P F_label.

Нужно проверить, что вклад на гранях x_i=0 и x_i=∞ исчезает в выбранной регулируемой области. Проверка должна учитывать:

    - степень x_i в P;
    - степень x_i в monomial part F_label;
    - valuation всех G_l при x_i=0 и x_i=∞;
    - знак epsilon в выбранном пределе;
    - возможные нули leading coefficient при assumptions.

7.2. Vector/tangent IBP

Для строки

    div(Q F_label)

проверяется нормальный поток через boundary ray rho:

    Flux_rho ~ (rho · Q/x) F_label.

Требуется положительный boundary exponent в регулируемой области. Все результаты кэшируются по паре (vector_field, ray).

7.3. Запрет ложных строк

Если surface_check=False, строка не попадает в систему как 0=... . Допустимо сохранить её отдельно как строку с boundary term, но в MVP такие строки не используются.


8. Локальная конечность финальных интегралов
--------------------------------------------

Функция:

    is_locally_finite(label) -> bool

должна проверять только epsilon=0.

Пример неправильной логики:

    base_score == 0 and eps_coeff < 0  -> считать локально конечным.

Это запрещено. Такой интеграл может быть регулируемо сходящимся, но не locally finite.

Правильная логика:

    locally finite iff base_score > 0 for every relevant ray.

Если тест неполон из-за неизвестной bulk singularity, результат должен быть:

    LocallyFinite -> "Unknown"

а не True. Для success все terms должны иметь True.


9. Производительность и работа с большими выражениями
-----------------------------------------------------

Общие правила:

    1. Не раскрывать глобально большие рациональные выражения.
    2. В hot loop работать только с целыми числами mod p и exponent tuples.
    3. SymPy использовать для парсинга, финальной факторизации и interpolation, но не для массового row generation/RREF.
    4. Все rows хранить разреженно.
    5. Row templates строить один раз и переиспользовать между parameter samples.
    6. Derivatives, valuations, monomial shifts, label ids кэшировать.
    7. Использовать streaming row generation, чтобы не держать лишние строки в памяти.
    8. Дедуплицировать строки по hash нормализованного sparse pattern.
    9. Поддерживать parallel execution по primes/samples.
    10. Не использовать Mathematica как обязательный backend. Wolfram-like текстовый output — только формат ввода/вывода.

Целевые ориентиры для первой рабочей реализации:

    - N до 6-8 на реальных задачах с тысячами labels;
    - N до 20 на синтетических sparse-тестах surface/valuation слоя;
    - 10^3-10^5 labels в зависимости от плотности rows;
    - 10^4-10^6 nonzero matrix entries;
    - finite-field solve должен быть главным вычислительным режимом.

Эти числа являются ориентирами производительности, а не математическим ограничением архитектуры.


10. Диагностика и статусы
-------------------------

Обязательные статусы:

    Success
    ParserFailure
    ParserNeedsExplicitFamily
    SurfaceRowsInsufficient
    TargetNotReducible
    NormalFormNotLocallyFinite
    InterpolationFailed
    VerificationFailed
    ResourceLimitReached

Для каждого failure нужно сохранять:

    - лучший найденный partial normal form;
    - список нелокально конечных terms, которые остались;
    - сколько rows было построено и сколько rejected;
    - какие search levels пробовались;
    - рекомендации: увеличить m-range, n-range, maxdeg, tangent degree, добавить preferred masters, добавить assumptions.

Пример failure-output:

    <|
      "Status" -> "Failure",
      "Error" -> "NormalFormNotLocallyFinite",
      "FormalSuccess" -> True,
      "AllLocallyFinite" -> False,
      "NonLocallyFiniteTerms" -> {...},
      "SuggestedNextOptions" -> <|
        "IncreaseMaxIBPDegree" -> 3,
        "IncreaseTangentDegree" -> {2,3},
        "ExpandMRange" -> {-4,1}
      |>
    |>


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


14. Критерии приёмки
--------------------

Код считается соответствующим ТЗ, если выполнены следующие условия.

A. Универсальность формы:

    - ни в одном ядровом модуле нет hardcoded N=3, N=4 или конкретных имён x4,x7,H,A,B,C;
    - arbitrary N проходит через tuple/list variables;
    - arbitrary parameters проходят через tuple/list parameters;
    - sparse-polynomial layer поддерживает произвольные мономиальные степени.

B. Корректность редукции:

    - все IBP rows проходят surface-free проверку;
    - all_lf=True для всех terms при success;
    - formal_success отделён от success;
    - локальная конечность проверяется при epsilon=0;
    - коэффициенты являются рациональными функциями parameters;
    - есть независимая finite-field валидация коэффициентов.

C. Производительность:

    - row generation и RREF не используют SymPy в hot loop;
    - repeated parameter samples переиспользуют row templates;
    - есть parallel execution по samples/primes;
    - есть memory-aware режим и streaming rows.

D. Вывод:

    - result.m может быть вставлен в Wolfram/Mathematica notebook без ручной правки как текстовый output;
    - coefficients и integrands выводятся в Wolfram-like syntax как текстовый формат;
    - labels и diagnostics сохраняются.

E. Регрессионные тесты:

    - N=2 sanity test проходит;
    - N=3 five-master test проходит;
    - D=4 symbolic coefficients test проходит;
    - ID4example2 coordinate-surface test проходит;
    - ID3example3 x4^2 sparse-polynomial test проходит.


15. Нежелательные реализации
----------------------------

Запрещённые или нежелательные подходы:

    - решать систему сразу над Q(ep,s,t,...) без modular layer;
    - использовать полное Expand больших выражений перед каждой строкой;
    - считать локально конечными интегралы, которые сходятся только благодаря epsilon;
    - считать редукцию успешной при all_lf=False;
    - хардкодить конкретные полиномы из тестовых примеров;
    - завязываться на установленную Mathematica как обязательный runtime;
    - восстанавливать коэффициенты по одной точке;
    - выбрасывать допустимые coordinate IBP из-за слишком строгой toric surface-check;
    - игнорировать surface terms;
    - не сохранять диагностику failure.


16. Итоговая спецификация результата
------------------------------------

Для успешного запуска пользователь должен получить файл result.m примерно такого вида:

    ClearAll[IBPReductionResult];

    IBPReductionResult = <|
      "Status" -> "Success",
      "AllLocallyFinite" -> True,
      "Variables" -> {x1,x2,x3,x4},
      "Parameters" -> {ep,r},
      "OriginalIntegrand" -> x4^(-ep)*(1+x1+x2+x3)^(2*ep)*(1+x4)^(1+3*ep)*
                             (r*x1*x2+r*x2*x3+r*x1*x2*x4+x3*x4)^(-2-ep),
      "Terms" -> {
        <|"Coeff" -> C1, "Integrand" -> F1, "Label" -> {...}, "LocallyFinite" -> True|>,
        <|"Coeff" -> C2, "Integrand" -> F2, "Label" -> {...}, "LocallyFinite" -> True|>
      },
      "IntegralRule" -> HoldForm[
        Int[OriginalIntegrand,{x1,0,Infinity},{x2,0,Infinity},{x3,0,Infinity},{x4,0,Infinity}]
        -> C1*Int[F1,{x1,0,Infinity},{x2,0,Infinity},{x3,0,Infinity},{x4,0,Infinity}]
         + C2*Int[F2,{x1,0,Infinity},{x2,0,Infinity},{x3,0,Infinity},{x4,0,Infinity}]
      ],
      "Diagnostics" -> <|...|>
    |>;

Это и есть конечный продукт: не просто список master labels, а готовое Wolfram-like текстовое представление исходного интеграла через локально конечные подынтегральные выражения и рациональные коэффициенты.
