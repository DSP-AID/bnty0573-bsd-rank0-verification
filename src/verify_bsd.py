"""High-precision verification of the Birch and Swinnerton-Dyer leading-term
formula for elliptic curves over Q of analytic rank 0.

For every curve in data/curves.json (input: a minimal Weierstrass model, i.e.
the coefficients [a1,a2,a3,a4,a6], plus the stored torsion order and the stored
conjectural order of Sha from LMFDB), this script computes from FIRST
PRINCIPLES, at working precision ~45 decimal digits:

  1. a_p by direct point counting over F_p (no tables, no external L-series);
  2. the local data at every bad prime p >= 5 by an explicit Tate-algorithm
     node test on the reduced cubic: split/nonsplit multiplicative type, the
     Hecke eigenvalue a_p = +/-1 and the Tamagawa number c_p;
  3. the real period Omega_E = integral over E(R) of the invariant
     differential dx/(2y + a1 x + a3), in closed form via the Legendre reduction
     2*K(m)/sqrt(rho) evaluated by AGM (all curves in
     the dataset have Delta < 0, so E(R) is connected and Omega_E is the full
     real period, with no component-group ambiguity);
  4. L(E,1) through the rapidly convergent series obtained from the completed
     L-function. For a weight-2 newform f of squarefree-ish level N with
     Fricke eigenvalue f|_2 W_N = -f (i.e. global root number w_E = +1, which
     holds automatically for analytic rank 0), one has

         L(E,1) = 2*pi * ( I_slow + 2*J ),
         I_slow = sum_{n>=1} a_n (e^{-2 pi n/N} - e^{-2 pi n}) / (2 pi n),
         J      = sum_{n>=1} a_n e^{-2 pi n} / (2 pi n),

     with a rigorous Hecke-type tail bound |a_n| <= d(n) sqrt(n).

and then checks the BSD prediction

  L(E,1) * #E(Q)_tor^2 / ( Omega_E * prod_p c_p ) == | Sha(E/Q) |

against the stored order of Sha. Agreement at >= 20 significant digits for
every curve is the verification result.

Usage:
    python src/verify_bsd.py [--data data/curves.json] [--out results/summary.json]
"""

from __future__ import annotations

import argparse
import json
import math
import os
import sys

import mpmath as mp

WORK_DP = 45  # working decimal digits (target agreement: >= 20 digits)


# ---------------------------------------------------------------- basic data

def b_invariants(ainvs):
    a1, a2, a3, a4, a6 = ainvs
    b2 = a1 * a1 + 4 * a2
    b4 = 2 * a4 + a1 * a3
    b6 = a3 * a3 + 4 * a6
    b8 = a1 * a1 * a6 + 4 * a2 * a6 - a4 * a4 + a1 * a3 * a4
    # Delta = disc(F)/16 for F(x) = 4x^3 + b2 x^2 + 2 b4 x + b6 (the model
    # (2y + a1 x + a3)^2 = F(x)); the cubic-discriminant formula is exact.
    A, B, C, D = 4, b2, 2 * b4, b6
    disc_f = (
        18 * A * B * C * D
        - 4 * B**3 * D
        + B * B * C * C
        - 4 * A * C**3
        - 27 * A * A * D * D
    )
    assert disc_f % 16 == 0
    delta = disc_f // 16
    c4 = b2 * b2 - 24 * b4
    return b2, b4, b6, b8, delta, c4


def legendre(v, p):
    r = pow(v % p, (p - 1) // 2, p)
    return 0 if r == p - 1 else r


# ------------------------------------------------------------- local primes

def count_points(ainvs, p):
    """#E(F_p) set-theoretic (includes the singular point when p is bad)."""
    a1, a2, a3, a4, a6 = ainvs
    if p == 2:  # direct enumeration; Legendre symbols degenerate mod 2
        cnt = 1
        for x in range(2):
            for y in range(2):
                lhs = y * y + a1 * x * y + a3 * y
                rhs = x**3 + a2 * x * x + a4 * x + a6
                if (lhs - rhs) % 2 == 0:
                    cnt += 1
        return cnt
    cnt = 1
    for x in range(p):
        d = ((a1 * x + a3) ** 2 + 4 * (x**3 + a2 * x**2 + a4 * x + a6)) % p
        cnt += 1 + (0 if d == 0 else 1 if legendre(d, p) == 1 else -1)
    return cnt


def local_bad_data(ainvs, N, p):
    """Tate node test for semistable bad primes p >= 5.

    Returns (a_p, c_p, v_p(Delta)). Semistable + p>=5 forces multiplicative
    reduction of type I_n with n = v_p(Delta) >= 1 and v_p(N) = 1.
    split <=> the two tangent lines at the node are defined over F_p
          <=> F''(x0)/2 is a square mod p, where x0 is the node abscissa of
            Y^2 = F(x) = 4x^3 + b2 x^2 + 2 b4 x + b6.
    """
    b2, b4, b6, _b8, delta, c4 = b_invariants(ainvs)
    assert delta % p == 0 and c4 % p != 0, "not multiplicative at p"
    v = 0
    d = abs(delta)
    while d % p == 0:
        d //= p
        v += 1
    assert N % p == 0 and (N // p) % p != 0, "bad prime must divide N exactly once"

    def F(x):
        return (4 * x**3 + b2 * x * x + 2 * b4 * x + b6) % p

    def Fp(x):
        return (12 * x * x + 2 * b2 * x + 2 * b4) % p

    node = [x for x in range(p) if F(x) == 0 and Fp(x) == 0]
    assert len(node) == 1, "node not unique mod p"
    x0 = node[0]
    half_fpp = (12 * x0 + b2) % p  # F''(x0)/2
    split = legendre(half_fpp, p) == 1
    a_p = 1 if split else -1
    c_p = v if split else math.gcd(2, v)

    # independent cross-check: #E_ns(F_p) = p - a_p (split: F_p^*, p-1 pts;
    # nonsplit: kernel of norm, p+1 pts)
    total = count_points(ainvs, p)
    ns = total - 1  # remove the singular node
    assert ns == p - a_p, "node test vs point count mismatch"
    return a_p, c_p, v


# -------------------------------------------------------------- coefficients

def hecke_coefficients(ainvs, N, bad_ap, nmax):
    """a_n for 1 <= n <= nmax from the a_p (p <= nmax) by multiplicativity.

    Good p: a_{p^{k+1}} = a_p a_{p^k} - p a_{p^{k-1}}; bad semistable p:
    a_{p^k} = a_p^k with a_p = +/-1 supplied in bad_ap.
    """
    primes = [
        q
        for q in range(2, nmax + 1)
        if all(q % r for r in range(2, math.isqrt(q) + 1))
    ]
    ap = {
        p: (bad_ap[p] if p in bad_ap else p + 1 - count_points(ainvs, p))
        for p in primes
    }
    a = [0] * (nmax + 1)
    a[1] = 1
    for p in primes:
        pk, e = p, 1
        while pk <= nmax:
            if p in bad_ap:
                val = ap[p] ** e
            elif e == 1:
                val = ap[p]
            else:
                val = ap[p] * a[p ** (e - 1)] - p * a[p ** (e - 2)]
            for m in range(1, nmax // pk + 1):
                if math.gcd(m, p) == 1:
                    a[m * pk] = a[m] * val
            pk *= p
            e += 1
    return a


# --------------------------------------------------------------- real period

def real_period(ainvs):
    """Omega_E = 2 * int_{e1}^inf dx / sqrt(F(x)), F = 4x^3+b2x^2+2b4x+b6.

    Requires Delta < 0 (F has one real root e1; E(R) connected). Closed form:
    deflate F = (x - e1)(4x^2 + c x + d); the quadratic factor has a
    conjugate pair alpha, conj(alpha). With rho = |alpha - e1| and
    m = (rho + Re(alpha) - e1)/(2 rho), Omega_E = 2 K(m)/sqrt(rho) (Legendre
    reduction, evaluated by AGM). Adaptive quadrature on the u-substituted
    integrand fails for near-singular models (77.b1: |Im alpha| ~ 0.01 gives
    an integrand spike of width ~1e-3); the closed form is immune.
    """
    b2, b4, b6, _b8, delta, _c4 = b_invariants(ainvs)
    assert delta < 0, "real period routine requires Delta < 0 (one component)"
    with mp.workdps(WORK_DP + 15):
        # Delta < 0 => F has exactly one real root; find it by bisection
        # (robust where polyroots iteration fails to converge).
        B = 2 + (abs(b2) + abs(2 * b4) + abs(b6)) / 4

        def f(x):
            return 4 * x**3 + b2 * x * x + 2 * b4 * x + b6

        lo, hi = mp.mpf(-B), mp.mpf(B)
        assert f(lo) < 0 < f(hi)
        for _ in range(350):
            mid = (lo + hi) / 2
            if f(mid) > 0:
                hi = mid
            else:
                lo = mid
        e1 = (lo + hi) / 2

        # exact deflation F = (x - e1)(4x^2 + c x + d): match x^2 and x coeffs
        c = b2 + 4 * e1
        d = 2 * b4 + e1 * c
        # constant term consistency: -e1*d == b6 up to the root residual
        assert abs(-e1 * d - b6) < mp.mpf(10) ** (-WORK_DP + 5), "deflation mismatch"
        p = -c / 8  # Re(alpha)
        q = mp.sqrt(16 * d - c * c) / 8  # |Im(alpha)| > 0 since Delta < 0
        rho = mp.sqrt((p - e1) ** 2 + q * q)
        m_par = (rho + p - e1) / (2 * rho)  # Legendre parameter, 0 < m < 1
        val = 2 * mp.ellipk(m_par) / mp.sqrt(rho)
    return mp.re(val)


# --------------------------------------------------------------------- L(1)

def l_value_at_1(ainvs, N, a_coeffs, nmax):
    """L(E,1) = 2*pi*(I_slow + 2J); returns (value, tail_bound)."""
    with mp.workdps(WORK_DP + 10):
        two_pi = 2 * mp.pi
        i_slow = mp.fsum(
            a_coeffs[n]
            * (mp.exp(-two_pi * n / N) - mp.exp(-two_pi * n))
            / (two_pi * n)
            for n in range(1, nmax + 1)
        )
        j_sum = mp.fsum(
            a_coeffs[n] * mp.exp(-two_pi * n) / (two_pi * n)
            for n in range(1, nmax + 1)
        )
        l1 = two_pi * (i_slow + 2 * j_sum)

        # rigorous tail bound: |a_n| <= d(n)*sqrt(n) <= 2n, so the omitted
        # tail of both sums is bounded by sum_{n>nmax} 2*e^{-c n}/pi with
        # c = 2*pi/N (the e^{-2 pi n} part is smaller).
        c = two_pi / N
        tail = 2 * mp.e ** (-c * (nmax + 1)) / (mp.pi * (1 - mp.e ** (-c)))
    return l1, tail


def choose_nmax(N, digits=WORK_DP - 5):
    """smallest nmax with e^{-2 pi nmax/N}/(pi(1-e^{-c})) < 10^-digits"""
    c = 2 * math.pi / N
    target = math.log(10) * (digits + 2)
    return math.ceil(target / c) + 8


# ------------------------------------------------------------------ verify

def verify_curve(curve, dps=WORK_DP):
    with mp.workdps(dps):
        ainvs = curve["ainvs"]
        n = curve["conductor"]
        t = curve["torsion"]
        s_stored = curve["sha"]
        bad_ap, tamagawa = {}, 1
        for p in curve["bad_primes"]:
            a_p, c_p, _v = local_bad_data(ainvs, n, p)
            bad_ap[p] = a_p
            tamagawa *= c_p
        nmax = choose_nmax(n)
        a_n = hecke_coefficients(ainvs, n, bad_ap, nmax)
        omega = real_period(ainvs)
        l1, tail = l_value_at_1(ainvs, n, a_n, nmax)
        ratio = (l1 * t * t) / (omega * tamagawa)  # predicted |Sha|
        pred = int(mp.nint(ratio))
        ok_square = pred >= 1 and math.isqrt(pred) ** 2 == pred
        rel = abs(ratio - s_stored) / s_stored
        passed = (
            pred == s_stored
            and ok_square
            and rel < mp.mpf("1e-20")
            and tail < mp.mpf("1e-30")
        )
        return {
            "label": curve["label"],
            "cremona": curve.get("cremona"),
            "conductor": n,
            "torsion": t,
            "tamagawa_product": tamagawa,
            "nmax_used": nmax,
            "sha_stored": s_stored,
            "sha_predicted": pred,
            "l1_digits": mp.nstr(l1, 30),
            "omega_digits": mp.nstr(omega, 30),
            "relative_error": mp.nstr(rel, 6),
            "tail_bound": mp.nstr(tail, 4),
            "pass": bool(passed),
        }


def main(argv=None):
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    pa = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    pa.add_argument("--data", default=os.path.join(root, "data", "curves.json"))
    pa.add_argument("--out", default=None)
    args = pa.parse_args(argv)

    with open(args.data, encoding="utf-8") as f:
        curves = json.load(f)

    print(f"BSD leading-term verification on {len(curves)} rank-0 curves " f"(dps={WORK_DP})")
    print(
        "label      N     T  prod c   Sha_pred Sha_stored  rel.err    tail     verdict"
    )
    results, nfail = [], 0
    for cv in curves:
        r = verify_curve(cv)
        results.append(r)
        if not r["pass"]:
            nfail += 1
        print(
            f"{r['label']:<10} {r['conductor']:<5} {r['torsion']:<3} "
            f"{r['tamagawa_product']:<6} {r['sha_predicted']:<8} "
            f"{r['sha_stored']:<10} {r['relative_error']:<9} "
            f"{r['tail_bound']:<9} {'PASS' if r['pass'] else 'FAIL'}"
        )
    print(f"\n{len(results) - nfail}/{len(results)} curves verified; failures: {nfail}")
    if args.out:
        os.makedirs(os.path.dirname(args.out), exist_ok=True)
        with open(args.out, "w", encoding="utf-8") as f:
            json.dump(results, f, indent=1)
    return 0 if nfail == 0 and results else 1


if __name__ == "__main__":
    sys.exit(main())
