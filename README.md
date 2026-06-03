# 6N Cousin & Sexy Prime-Pair Distortion (Part IX)

Does the twin's conditional omega-enrichment carry over to cousin (difference 4)
and sexy (difference 6) prime pairs? It does not — it distorts into three
distinct, geometry-determined shapes.

**Geometry on the 6N skeleton.** All primes >3 lie on 6N±1.
- **twin** (diff 2): (6N−1, 6N+1) — single centre N, both wings.
- **cousin** (diff 4): (6N+1, 6(N+1)−1) — straddles N, N+1. (The other option
  6N+3 is always divisible by 3.)
- **sexy-A** (diff 6): (6N−1, 6(N+1)−1) — straddles N, N+1, both left wings.
- **sexy-B** (diff 6): (6N+1, 6(N+1)+1) — straddles N, N+1, both right wings.

Pairs are binned by the **left centre's** omega₍>3₎(N), the analogue of the twin
attribution.

**Result (S₁₀, 1.5×10⁹ centres, >5×10⁸ pairs), normalised to omega=1:**

| omega | twin | cousin | sexy-A | sexy-B |
|------:|-----:|-------:|-------:|-------:|
| 1 | 1.000 | 1.000 | 1.000 | 1.000 |
| 4 | 1.936 | 0.804 | 0.804 | 1.198 |
| 6 | 3.524 | 0.420 | 0.414 | 0.881 |
| 7 | 4.808 | 0.127 | 0.150 | 0.402 |

Three distortions:
- **twin rises** monotonically (single centre — a factor-rich N protects both wings).
- **cousin and sexy-A fall** monotonically, and are **numerically identical** at
  every omega. Cause: they share the right member 6(N+1)−1; after normalisation
  the left member's omega-response divides out, leaving the shared right member's
  response on the *adjacent* centre N+1 — common to both, hence coincident. The
  fall (vs the twin's rise) reflects that the partner lives on N+1, not N.
- **sexy-B is non-monotone**, peaking near omega=4 then falling — two competing
  effects, not separated here.

The distortion's sign and shape are set by the **pairing geometry** (which
centres, which wings), not by the difference alone: diff-4 and diff-6 can coincide
(cousin ≡ sexy-A), while two diff-6 forms differ (sexy-A vs sexy-B).

> **Scope.** Experimental / computational number theory; phenomenon + geometric
> account. A quantitative two-centre mechanism is left as the open problem. No
> claim about the infinitude of twin, cousin, or sexy primes, or any k-tuple
> conjecture.

**Open problem.** Predict the three curves quantitatively from a two-centre
conditional model over (N, N+1) using the modular-shift/lockdown mechanism of
Part V: twin's rise, the fall and exact coincidence of cousin and sexy-A, and the
non-monotone peak of sexy-B.

Part I: doi:10.5281/zenodo.20470367 · II: doi:10.5281/zenodo.20477664 ·
V: doi:10.5281/zenodo.20510700 · VIII: doi:10.5281/zenodo.20519998

---

## Layout

```
.
├── README.md
├── LICENSE                 (MIT)
├── CITATION.cff
├── data/
│   └── cousin_sexy_S10_data.csv   omega, Ncenters, rates + normalised, S10
├── code/
│   ├── cousin_sexy.py      streaming scan; counts twin/cousin/sexy-A/sexy-B by
│   │                       left-centre omega; emits cousin_sexy_S{K}_data.csv
│   └── make_cs_fig.py      builds the 2-panel figure from ../data
├── figures/                fig_paper9_cousin_sexy.{pdf,png}
└── paper/                  Chen_6N_Paper9.{tex,pdf} + figure
```

## Reproducing

Requirements: Python 3.8+, `numpy`, `matplotlib`.

```bash
pip install numpy matplotlib

# 1. Scan and count. Default S10 (~16 min). Memory-light streaming (does NOT
#    store all centres; carries one-centre overlap across segments for the
#    straddling pairs). Emits cousin_sexy_S{K}_data.csv.
python code/cousin_sexy.py            # S10
MAXK=9 python code/cousin_sexy.py     # S9 (faster, for validation)

# 2. Figure (reads ../data/cousin_sexy_S10_data.csv).
cd code && python make_cs_fig.py
```

### Conventions (same as Parts I–VIII)

- Twin centre N: 6N−1, 6N+1 both prime. omega₍>3₎(N) = count of distinct prime
  factors >3 of N.
- Straddling pairs (cousin, sexy) are attributed to the LEFT centre's omega.
- Rates are pair counts per left-centre, normalised to omega=1 for shape comparison.
- Engine: complete segmented-sieve factorisation + deterministic interval-sieve
  primality; S₁₀ twin count 23,988,173 matches Part I.

## License

MIT — see `LICENSE`.
