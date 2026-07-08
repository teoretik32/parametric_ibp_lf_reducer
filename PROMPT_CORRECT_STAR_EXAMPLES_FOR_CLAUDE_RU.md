Вставь в Claude Code немедленно, если предыдущий notebook patch уже был применён:

```text
Стоп и поправка по notebook examples.

Я получил уточнение: examples со звездой (`Example 3*`, `Example 4*`) были поняты неправильно. Это отдельный тип примеров: там известен answer/epsilon expansion для самого интеграла, но НЕ известно разложение по locally finite integrands. Они никак не связаны с Example 1/2 как LF-decomposition examples.

Сделай correction pass перед продолжением Pass 1A:

1. Прочитай docs/05_star_examples_policy_correction_ru.md, если файл есть. Если файла нет, всё равно следуй этой инструкции.
2. Обнови notes/assumptions.md: starred examples are known-value-only, not LF decomposition fixtures.
3. Обнови notes/test_strategy.md:
   - Example 3* и 4* только parser/sparse-poly fixtures сейчас;
   - no e2e reducer expected outputs for them;
   - main-spec 11.4/11.5 remain xfail/pending input family unless separately supplied.
4. Удали или пометь xfail любые tests, которые пытаются редуцировать Example 3* к ID4example2 expected one-term или Example 4* к ID3example3 basis inclusion.
5. Не используй validation/id4example2_expected_one_term.json и validation/id3example3_expected_basis_inclusion.json как expected для starred examples; если они есть, считать deprecated.
6. Не используй tm1Coeff/Gamma-prefactor как reducer coefficient.
7. Сохрани полезное: Example 4* можно использовать для parser/sparse-poly test на monomial x4^2.

После correction pass продолжай только Pass 1A и запусти python -m pytest.
```
