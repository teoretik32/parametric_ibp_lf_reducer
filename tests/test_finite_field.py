"""Tests for GF(p) arithmetic helpers (Pass 2E)."""

from __future__ import annotations

import random

import pytest

from parametric_ibp_lf_reducer import (
    batch_inverse,
    generate_primes,
    inv_mod,
    is_probable_prime,
    powmod,
)
from parametric_ibp_lf_reducer.finite_field import add_mod, mul_mod, sub_mod

P = 2_147_483_647  # 2**31 - 1


def test_inv_mod_roundtrip():
    rng = random.Random(1)
    for _ in range(200):
        a = rng.randrange(1, P)
        assert a * inv_mod(a, P) % P == 1


def test_inv_mod_zero_raises():
    with pytest.raises(ZeroDivisionError):
        inv_mod(0, P)
    with pytest.raises(ZeroDivisionError):
        inv_mod(P, P)  # 0 mod P


def test_powmod_matches_builtin_and_handles_negative():
    rng = random.Random(2)
    for _ in range(100):
        a = rng.randrange(1, P)
        e = rng.randrange(0, 50)
        assert powmod(a, e, P) == pow(a, e, P)
    a = 12345
    assert powmod(a, -1, P) == inv_mod(a, P)
    assert powmod(a, -3, P) == pow(inv_mod(a, P), 3, P)


def test_basic_mod_ops():
    assert add_mod(P - 1, 5, P) == 4
    assert sub_mod(3, 10, P) == (P - 7)
    assert mul_mod(P - 1, P - 1, P) == 1


def test_batch_inverse_matches_individual():
    rng = random.Random(3)
    vals = [rng.randrange(1, P) for _ in range(50)]
    got = batch_inverse(vals, P)
    assert got == [inv_mod(v, P) for v in vals]
    assert batch_inverse([], P) == []


def test_batch_inverse_rejects_zero():
    with pytest.raises(ZeroDivisionError):
        batch_inverse([3, 0, 5], P)


def test_is_probable_prime():
    assert is_probable_prime(P)  # 2**31 - 1 is prime
    assert is_probable_prime(2) and is_probable_prime(97)
    assert not is_probable_prime(1)
    assert not is_probable_prime(P - 2)  # composite
    assert not is_probable_prime(561)  # Carmichael number


def test_generate_primes():
    ps = generate_primes(5, upper=1000)
    assert len(ps) == 5
    assert all(is_probable_prime(p) for p in ps)
    assert len(set(ps)) == 5
    assert all(p < 1000 for p in ps)
    assert ps == sorted(ps, reverse=True)  # descending
