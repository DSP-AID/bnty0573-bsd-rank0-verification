# bsd-attack-bnty0573

First-principles, 45-digit independent verification of the rank-0 Birch and
Swinerton-Dyer leading-term formula on a dataset of semistable elliptic curves
over Q with conductor <= 200. Companion write-up: [ANALYSIS.md](ANALYSIS.md).

## Layout

```
src/verify_bsd.py        the verifier (pure Python + mpmath, no other deps)
data/curves.json         input dataset (minimal models from LMFDB download_all)
results/run_log.txt      full verification table of the published run
results/summary.json     machine-readable results (L(E,1), Omega_E, Sha, errors)
tests/test_verify_bsd.py hardcoded self-checks on curves 11.a2 and 77.b1 (no network)
```

## Reproduce

```bash
pip install mpmath pytest ruff
python src/verify_bsd.py --out results/summary.json   # prints the table
pytest tests/
ruff check src tests
```

Expected: every printed line ends in `PASS`, final line reports 0 failures,
exit code 0. Sha orders are predicted from scratch and compared against the
LMFDB-stored conjectural orders; agreement is at >= 20 significant digits.

## Data provenance

Inputs (minimal equation coefficients, conductor, torsion order, stored Sha)
were fetched on 2026-11-11 from LMFDB `download_all` endpoints; the fetch
script and filter criteria are in the DAiL submission notes. All analytic
quantities (`a_p`, `L(E,1)`, `Omega_E`, `c_p`) are recomputed here from the
model alone — no LMFDB L-values or periods are used by the verifier.

## License

MIT
