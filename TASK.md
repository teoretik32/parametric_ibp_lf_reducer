# TASK: сценарий 1 — пиши код с нуля и проверяйся на примерах

Сделать Python-пакет `parametric_ibp_lf_reducer` с архитектурой из cleaned spec.

## Цель

Принять параметрический интеграл в Wolfram-like/Mathematica-style текстовом формате, построить IBP/algebraic редукцию в параметрическом представлении и выдать представление исходного интеграла через локально конечные интегралы.

## Что делать

1. Начать с архитектурного плана в `notes/implementation_plan.md`.
2. Реализовать Python-модули в `src/parametric_ibp_lf_reducer/`.
3. Добавить parser explicit family для Wolfram-like text association.
4. Реализовать sparse-poly core без SymPy в hot loop.
5. Реализовать local-finiteness при `epsilon=0` и surface-free checks.
6. Реализовать row generators: algebraic, coordinate IBP, tangent/syzygy IBP.
7. Реализовать ranking, sparse finite-field RREF, modular samples.
8. Реализовать rational-function reconstruction и independent checks.
9. Реализовать Wolfram-like text export и CLI/API.
10. Сделать pytest regression tests из `docs/03_validation_cases_ru.md`.

## Чего не делать

- Не писать Mathematica/Wolfram implementation.
- Не требовать Mathematica runtime.
- Не хардкодить validation cases.
- Не считать epsilon-regulated convergence локальной конечностью.
- Не возвращать Success при `AllLocallyFinite=False`.
