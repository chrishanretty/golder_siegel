# Summary of the Golder, Golder & Siegel (2012) ABM

Golder, M., Golder, S. N., & Siegel, D. A. (2012). "Modeling the Institutional
Foundation of Parliamentary Government Formation." *Journal of Politics* 74(2): 427–445.

## Idea

This is a "zero-intelligence" null model of government formation. Parties make
**random** government proposals. The only constraints are the two constitutional
constraints common to every parliamentary system:

1. an incumbent (caretaker) government always exists and is the status quo (SQ);
2. a new government needs the support of its members and of a legislative majority.

## Party system (the election stage, which involves no choices)

| Parameter | Distribution |
|---|---|
| number of parties *m* | uniform on {3,…,10} |
| policy weight β (common to all parties) | U[0,1] |
| party ideal points *z_i* | U([-3,3]²) |

- Voters are evenly spaced over [-3,3]², so the median voter is at (0,0). They are
  proximity voters, with ties broken by a fair die.
- Seat share equals vote share (pure PR).

## Utility

u_i = cs_i − β‖x − z_i‖²

- cs_i is party i's cabinet share (0 if it is not in government).
- x is the government's policy, a point in R².

## Government formation (a Markov chain over status-quo governments)

Each period:

1. **Coalition.** Every party proposes a coalition that contains itself, drawn
   uniformly from all 2^(m−1) subsets of the other parties.
2. **Portfolios and policy.** It also draws a cabinet distribution uniformly on the
   simplex over the coalition members, and a policy uniformly on [-3,3]².
3. **Qualifying test.** A proposal qualifies if every coalition member strictly
   prefers it to the SQ, and the parties that strictly prefer it hold a simple
   majority of seats.
4. **Installation.** If one or more proposals qualify, the head of state installs
   one at random. Its proposer is the *formateur*, and it becomes the new SQ.
   Otherwise every proposal is discarded and the parties go back to step 1.
5. **Cut-off.** After 100 failed periods the caretaker (SQ) government is
   reinstalled.

**Delay** is the number of periods (1–100) until a government is installed.

- The initial SQ is policy (0,0) with zero cabinet shares for every party.
- A **run** is one party system followed by 2,500 successive formations. Outcomes
  are taken at the 2,500th formation, which acts as burn-in to the ergodic limit.
- The paper simulates 50,000 runs, one per random party system.

## Outcome measures

- **Government type:**
  - single-party majority;
  - single-party minority;
  - minority coalition;
  - minimal winning coalition (every member pivotal);
  - surplus majority (including "unity" governments of all parties).
- **Gamson regression:** cabinet share on the share of the government's seats a
  party contributes. Unit = governing party; SEs clustered by government.
- **Effective number of parties:** ENPP = 1/Σs².
- **Ideological diversity:** Σ s_j ‖z_j − z̄‖ / σ_c, where z̄ is the seat-weighted
  mean of party positions.
- **Weights:** each run is weighted by a kernel density estimate (likelihood-CV
  bandwidth) of the ENPP of 406 West European governments, 1945–98, evaluated at
  the run's ENPP.

## Ambiguities and how they were resolved

| Issue | Choice | Basis |
|---|---|---|
| Number and layout of voters | 1,000 voters on a 40×25 cell-centred grid | Every vote share in the GGS data is a multiple of 0.001. Mean and sd of ENPP by *m* match the GGS data. |
| σ_c in the diversity measure | 2.5 | The ratio of our raw measure to GGS's variable is 2.50 for every *m*. |
| "Simple majority" | Seat share strictly greater than 0.5 | Standard definition. |
| Formateur after the caretaker is reinstalled | The previous formateur is kept | Not specified in the paper. |
| Unity governments | Merged into surplus | As in the GGS replication script. |
| Weights | Kernel density of real-world ENPP, likelihood-CV bandwidth (h = 0.162) | GGS replication script (np::npudens). |

The authors' replication deposit is Harvard Dataverse, doi:10.7910/DVN/KIFRVV. It
holds their simulated output, the real-world data and the analysis script, but no
model code. It is stored in `data/original/`.
