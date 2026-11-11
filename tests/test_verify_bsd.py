"""Self-checks for the BSD verifier on curves 11.a2 (Cremona 11a1) and 77.b1.

All reference values are hardcoded; no network access is required.
"""

import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

import mpmath as mp
from verify_bsd import (
    b_invariants,
    count_points,
    local_bad_data,
    real_period,
    verify_curve,
)

AINVS = [0, -1, 1, -10, -20]  # y^2 + y = x^3 - x^2 - 10x - 20


def test_delta_and_c4():
    _b2, _b4, _b6, _b8, delta, c4 = b_invariants(AINVS)
    assert delta == -161051  # -11^5
    assert c4 == 496


def test_point_counts_11a1():
    # a_p = p + 1 - #E(F_p); matches the LMFDB q-expansion of 11.a2
    expected = {2: -2, 3: -1, 5: 1, 7: -2, 13: 4}
    for p, a_p in expected.items():
        assert p + 1 - count_points(AINVS, p) == a_p


def test_local_data_at_11():
    a_p, c_p, v = local_bad_data(AINVS, 11, 11)
    assert (a_p, c_p, v) == (1, 5, 5)  # split multiplicative type I_5


def test_real_period_11a1():
    with mp.workdps(35):
        omega = real_period(AINVS)
        ref = mp.mpf("1.2692093042795534216887946168")
        assert abs(omega - ref) < mp.mpf("1e-24")


def test_real_period_77b1_spike():
    """Regression: 77.b1 (y^2+y=x^3+x^2-89x+295) is near-singular, |Im alpha|
    ~ 0.01; naive quadrature of the period integral is off by ~7%. Reference
    value from LMFDB (30 digits)."""
    with mp.workdps(35):
        omega = real_period([0, 1, 1, -89, 295])
        ref = mp.mpf("4.6467301969045698455666134453")
        assert abs(omega - ref) < mp.mpf("1e-28")


def test_full_bsd_check_11a1():
    curve = {
        "label": "11.a2",
        "cremona": "11a1",
        "conductor": 11,
        "ainvs": AINVS,
        "torsion": 5,
        "sha": 1,
        "bad_primes": [11],
    }
    res = verify_curve(curve)
    assert res["pass"] is True
    assert res["sha_predicted"] == 1
    assert float(res["relative_error"]) < 1e-20
