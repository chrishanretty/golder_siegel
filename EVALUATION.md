# Evaluation: how well does the R replication match GGS (2012)?

The full numbers are in `outputs/results.md`, and the figures are in `outputs/fig1.png`,
`outputs/fig2.png` and `outputs/fig3.png`.

The main batch is 50,000 party systems × 2,500 formations each. It took 22.8 minutes
on 8 cores.

## How the match was evaluated

There are two benchmarks:

1. **The published tables.**
2. **The authors' own simulated data** from Dataverse. Our analysis code reproduces
   *every* published number from these data exactly, including Tables 1–4, the mean
   delays and the ENPP weights. Any remaining difference is therefore in the model,
   not in the analysis.

Sampling noise is small at n = 50,000. The SE of a proportion is about 0.0015 and
the SE of mean delay is about 0.16. Unlike the paper's statistics, the
replication's are therefore precise enough to show small systematic gaps.

## Results (main specification: Euclidean voting, as stated in the paper)

| | Published | Replication |
|---|---:|---:|
| **Table 1, unweighted:** SPMaj / MWC / MinCoal / SPMin / Surplus | .086 / .196 / .461 / .150 / .107 | .074 / .198 / .464 / .152 / .112 |
| **Table 1, weighted** | .157 / .257 / .254 / .230 / .102 | .135 / .264 / .262 / .230 / .109 |
| **Table 2:** slope / intercept, unweighted | .744 / .102 | .742 / .102 |
| **Table 2:** slope / intercept, weighted | .771 / .117 | .767 / .117 |
| **Table 3:** formateur bonus | .076 | .078 |
| **Table 3:** seat share × m / × ENPP / × diversity / × β | −.098 / .046 / −.130 / .040 | −.089 / .044 / −.136 / .030 |
| **Table 4, OLS:** β / m / ENPP / diversity | 19.5 / 8.82 / −20.3 / −40.7 | 18.7 / 7.92 / −17.8 / −43.7 |
| **Table 4, Cox:** β / m / ENPP / diversity | −.82 / −.27 / .53 / 1.42 | −.82 / −.25 / .48 / 1.57 |
| **Mean delay:** unweighted / weighted | 29.13 / 29.84 | 27.04 / 27.32 |

Further checks:

- **Figure 1** (the 5-party system in footnote 16): slope and intercept at
  formation 2,500 are 0.86 / 0.06 in both. Mean delay, the formateur distribution
  and the government-size distribution match GGS's run.
- **Figure 2:** the GAM curves for government type against β overlap GGS's within
  confidence bands.
- **Figure 3:** the delay density is nearly identical to GGS's, with a slightly
  smaller spike at 100 periods.
- **Conditional on m:** every pattern matches, including the U-shape of delay in m,
  the type mix and the decline of the Gamson slope with m.

Every sign and every qualitative conclusion in the paper replicates.

## Residual differences and diagnosis

1. **Single-party majority is 1.2 points too low** (2.2 points weighted). This is
   entirely an election-stage effect. Whenever a party holds a seat majority, both
   GGS and we always end with a single-party majority government. GGS's party
   systems simply contain a majority party more often (47% vs 43% at m=3), so their
   seat shares are more dispersed than Voronoi cells of a uniform voter grid on
   [−3,3]² give. Mean and sd of ENPP by m still match.
   - **Other voter layouts:** grids of other shapes and sizes, random voters, a
     disc, a torus and L1 voting all fail to match the seat distribution.
   - **Chebyshev (L∞) voting** reproduces the frequency of majority parties.
2. **Mean delay is about 2 periods (7%) too short.** The gap persists within
   m × ENPP strata and after regression adjustment (−1.7), and is largest at small m.
   - **Ruled out by experiment** (`R/experiments.R`): the cabinet-draw method
     (Dirichlet vs normalised uniforms), a wider policy-proposal range, and
     rescaling β (both give the wrong shape across m).
   - **L∞ voting closes about two-thirds of it**, which points to the same
     election-geometry detail as gap 1.
3. **Sensitivity run with L∞ voting** (20,000 runs):
   - **Better:** single-party majority (.087 vs .086), the m=3–4 delays and the
     Table 4 coefficients all match.
   - **Worse:** some Table 3 interaction terms (ENPP, diversity).
   - **Conclusion:** GGS's exact voter geometry is unknown, since the deposit holds
     no model code. We keep Euclidean voting, as stated in the paper, as the main
     specification.

## Calibrated or inferred details

These are not stated in the paper and were inferred from the authors' deposited data:

- **1,000 voters:** every vote share is an exact multiple of 0.001 and the shares
  sum exactly to 1.
- **σ_c = 2.5** in the diversity measure: the ratio to GGS's variable is constant
  across m.
- **The previous formateur is kept** when the caretaker is reinstalled.

The ENPP weights replicate GGS's `np::npudens` likelihood-CV kernel density.
