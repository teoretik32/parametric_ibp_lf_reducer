# Implementation Plan — `parametric_ibp_lf_reducer`

Статус: черновик к первому coding pass. Живой документ — обновляется после каждого слоя.

## 0. Цель и жёсткие инварианты

Пакет принимает параметрический интеграл Фейнмана в Wolfram-like тексте, строит
IBP/algebraic редукцию в параметрическом представлении и выдаёт
`Integral[F0] -> Sum_a C_a(params) Integral[F_a]`, где все `F_a` **локально конечны при
`epsilon=0`**.

Инварианты (нарушение = баг, не «фича»):
1. `Status=Success` ⟺ `AllLocallyFinite=True` + surface-checks пройдены + reconstruction
   проверен на независимых samples.
2. Formal normal form с не-LF terms → `Failure`/partial, **не** Success.
3. LF проверяется строго при `eps=0`. Знак регулятора ≠ LF.
4. SymPy запрещён в hot loop (row-generation/RREF). Разрешён: parsing, tangent-syzygy setup
   (малый, разовый), финальная факторизация, reconstruction.
5. Основной контур — sparse / modular / finite-field.
6. Ничего не хардкодить: `N`, `M`, имена переменных, конкретные полиномы кейсов.

## 1. Module map (`src/parametric_ibp_lf_reducer/`)

Структура spec §5 + два добавления: `coefficients.py`, `labels.py`.

| Модуль | Ответственность | Ключевые функции/классы |
|---|---|---|
| `input_parser.py` | Wolfram-like association parser; explicit family (3.1); best-effort integrand factor (3.2) | `parse_mathematica_association`, `parse_explicit_family`, `try_factor_integrand` |
| `coefficients.py` | Параметрический rational expr; конверсия из SymPy 1 раз; быстрый mod-p eval | `ParamExpr`, `from_sympy`, `eval_mod_p(sample, prime)` |
| `sparse_poly.py` | `SparsePoly=dict[ExpTuple, CoeffId]`; арифметика без SymPy | `add/mul/pow_small/derivative/monomial_mul/support/degree/valuation/eval_mod_p` |
| `labels.py` | Label=`tuple[int,...]` длины `N+M`; box enumeration; id-map; complexity | `LabelSpace`, `enumerate_box`, `label_id`, `complexity` |
| `family.py` | `ParametricFamily`; переход label→factor/text | `label_to_factor`, `label_to_wolfram`, `exponent_at_label`, `specialize` |
| `valuations.py` | Rays + LF при `eps=0` | `compute_candidate_rays`, `valuation_poly`, `base_score`, `is_locally_finite` |
| `surface.py` | Surface-free фильтры | `coordinate_primitive_surface_free`, `vector_field_surface_free` |
| `tangent_fields.py` | Syzygy ansatz `ΣQ_i∂_iG_l−H_lG_l=0` по degree-блокам | `generate_tangent_fields` |
| `row_generation.py` | Row templates: algebraic / coordinate IBP / tangent IBP; dedup; surface-filter | `generate_algebraic_rows`, `generate_coordinate_ibp_rows`, `generate_tangent_ibp_rows` |
| `ranking.py` | Tiers + `complexity`; порядок исключения колонок | `build_ranking` |
| `finite_field.py` | GF(p) арифметика | `inv`, `batch_inverse`, `powmod` |
| `sparse_rref.py` | Sparse modular RREF с заданным порядком pivot; normal form target | `sparse_rref_normal_form` |
| `interpolation.py` | Univariate (Thiele/Padé) + multivariate dense-rational ansatz | `interpolate_rational` |
| `reconstruction.py` | CRT + rational reconstruction; сбор records → coeff functions; валидация | `interpolate_normal_form_coefficients` |
| `reducer.py` | Оркестратор spec §6 | `reduce_to_locally_finite` |
| `diagnostics.py` | Статусы, счётчики, recommendations, partial | `Diagnostics`, статус-энум |
| `wolfram_text_export.py` | Wolfram-like текст (`^`, `/`, `<||>`, `->`; `3/7`; factorized) | `coeff_to_wolfram_text`, `integrand_to_wolfram_text`, `result_to_wolfram_text` |
| `__init__.py` | Публичный API | `reduce_wolfram_style_input`, `ReductionResult`, `ReductionTerm`, `Diagnostics` |
| `__main__.py` | CLI | `main`, `reduce` subcommand + опции spec §12 |

## 2. Внутренние представления

- **ExponentTuple**: `tuple[int,...]` длины `N` (степени `x_i` в мономе). Канонически
  отсортированы в `SparsePoly` по ключам dict.
- **CoeffId**: индекс в таблице `ParamExpr` (интернирование коэффициентов полиномов `G_l`,
  чтобы hot loop работал с int-id, а не с выражениями).
- **Label**: `tuple[int,...]` длины `N+M` = `(n_1..n_N, m_1..m_M)`.
- **Row template**: список `(label_id, ParamExpr_coeff)` — параметрический, инстанцируется
  mod p на каждом sample. Хранится sparse; ключ дедупа — hash нормализованного паттерна.
- **Sample**: dict `param_name -> rational/int`. **Prime**: int.

## 3. Математические контракты (реализационные)

- **Coordinate IBP** `0=∫∂_{x_i}(P·F)`:
  `∂_i(PF)=(∂_iP)F + P·F·[(a_i+n_i)/x_i + Σ_l(b_l+m_l)∂_iG_l/G_l]`.
  `/G_l` → `m_l→m_l−1`; мономы `∂_iG_l`, `∂_iP`, `1/x_i` → n-сдвиги. Строка входит только
  если `P·F` surface-free на `x_i=0,∞` (см. `surface.coordinate_primitive_surface_free`).
- **Algebraic** `J(n,m)−Σ_a c_{l,a}J(n+a,m−e_l)=0` из `G_l=Σ_a c_{l,a}x^a`. Точное, всегда генерим.
- **Tangent IBP** `0=∫div(QF)`, `Q·∂G_l=H_lG_l` ⇒ без m-сдвига:
  `∫F[divQ + Σ_iQ_i(a_i+n_i)/x_i + Σ_l(b_l+m_l)H_l]`.
- **Local finiteness (eps=0)**: `base_score(label,rho)=Σ_i rho_i(e_i+1)+Σ_l f_l·val_rho(G_l)`,
  `e_i=a_i+n_i`, `f_l=b_l+m_l` при `eps=0`, `val_rho(G_l)=min_{a∈supp G_l}(a·rho)`.
  LF ⟺ `base_score>0` для всех образующих лучей общего измельчения нормальных фанов.
  Кандидаты лучей: `±e_i`; фасетные нормали Newton-политопов `G_l` и Минковски-суммы;
  adaptive random rays как safety net → `Unknown` (не True) при находке дивергенции вне набора.

## 4. Порядок реализации (passes)

**Pass 1 (foundation, первый coding pass):** `input_parser`, `coefficients`, `sparse_poly`,
`labels`, `family`, `valuations`, `surface` + тесты. Линейной алгебры нет.

**Pass 2:** `tangent_fields`, `row_generation`, `ranking` + row-тесты и полный 11.1.

**Pass 3:** `finite_field`, `sparse_rref` + RREF/normal-form тесты.

**Pass 4:** `reconstruction`, `interpolation` + тесты на известных rational functions.

**Pass 5:** `reducer`, `diagnostics`, `wolfram_text_export`, API, CLI + end-to-end 11.3.

## 5. Test queue (укрупнённо; детали в `test_strategy.md`)

1. Parser: `examples/d4_explicit_family.wl.txt` → корректная family; ошибочный ввод → диагностика.
2. Sparse-poly: произвольный `N`, мономы степени >1 (`x4^2`), derivative, valuation, mod-p eval.
3. LF при `eps=0`: положительные/отрицательные кейсы; «сходится из-за eps» ⇒ не LF; Unknown-путь.
4. Surface: coordinate IBP проверяет только `x_i=0,∞`; toric flux для vector IBP.
5. 11.1 tangency: `G=1+x+y`, `Q=(xy,−xy)`, `Q·∇G=0`, row без m-shift (полностью — pass 2).
6. RREF (pass 3), reconstruction (pass 4).
7. End-to-end 11.3 (pass 5): LF-базис `M1..M5`, коэффициенты `C1..C5`, `ep=-3/4,r=7/5`.

## 6. MVP vs Full target

**MVP:** всё выше до зелёного end-to-end 11.3 + 11.1; честные статусы; export/API/CLI.
**Full:** robust integrand auto-factor; Zippel/Ben-Or–Tiwari sparse interpolation;
несколько регуляторов; external CAS syzygy; numba/parallel/streaming; divergence certificate;
numerical euclidean check; фикстуры 11.2/11.4/11.5 при наличии их input families.

## 7. Известные пробелы во входных данных

11.2 (N=3 five-master), 11.4 (ID4example2), 11.5 (ID3example3): в docs есть только ожидаемые
LF-базисы/коэффициенты, **нет входных families**. Оформляем как `xfail`-фикстуры «pending
input family». MVP-приёмка на них не завязана (якорь — 11.1 + 11.3).
