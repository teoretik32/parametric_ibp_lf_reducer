"""GF(p) arithmetic helpers for the modular reduction contour (spec §5.9).

Pure integer arithmetic only — no SymPy. These are the primitives used by the sparse modular
RREF and (later) by matrix assembly and reconstruction.
"""

from __future__ import annotations

from collections.abc import Sequence


def add_mod(a: int, b: int, p: int) -> int:
    return (a + b) % p


def sub_mod(a: int, b: int, p: int) -> int:
    return (a - b) % p


def mul_mod(a: int, b: int, p: int) -> int:
    return (a * b) % p


def inv_mod(a: int, p: int) -> int:
    """Modular inverse via Fermat's little theorem (``p`` must be prime)."""
    a %= p
    if a == 0:
        raise ZeroDivisionError("inverse of 0 modulo p")
    return pow(a, p - 2, p)


def powmod(a: int, e: int, p: int) -> int:
    """``a**e mod p`` supporting negative exponents (via the modular inverse)."""
    if e < 0:
        return pow(inv_mod(a, p), -e, p)
    return pow(a % p, e, p)


def batch_inverse(values: Sequence[int], p: int) -> list[int]:
    """Invert many nonzero residues with a single modular inversion (Montgomery's trick)."""
    n = len(values)
    if n == 0:
        return []
    prefix = [1] * (n + 1)
    for i, a in enumerate(values):
        prefix[i + 1] = prefix[i] * (a % p) % p
    if prefix[n] == 0:
        raise ZeroDivisionError("batch_inverse requires all values nonzero modulo p")
    running = inv_mod(prefix[n], p)
    out = [0] * n
    for i in range(n - 1, -1, -1):
        out[i] = running * prefix[i] % p
        running = running * (values[i] % p) % p
    return out


# --- prime utilities (small deterministic Miller-Rabin) ----------------------
_MR_WITNESSES = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37)


def is_probable_prime(n: int) -> bool:
    """Deterministic Miller-Rabin for all ``n < 3.3e24`` (covers 31/63-bit primes)."""
    if n < 2:
        return False
    for w in _MR_WITNESSES:
        if n % w == 0:
            return n == w
    d = n - 1
    r = 0
    while d % 2 == 0:
        d //= 2
        r += 1
    for a in _MR_WITNESSES:
        x = pow(a, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(r - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def generate_primes(count: int, upper: int = 2**31 - 1) -> list[int]:
    """Return ``count`` distinct primes just below ``upper``, descending."""
    primes: list[int] = []
    n = upper if upper % 2 else upper - 1
    while len(primes) < count and n > 2:
        if is_probable_prime(n):
            primes.append(n)
        n -= 2
    return primes
