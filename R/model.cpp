// Core government-formation loop of the Golder, Golder & Siegel (2012)
// zero-intelligence model. Uses R's RNG (unif_rand) so that results are
// reproducible with set.seed() and identical to the pure-R reference
// implementation in ggs_model.R (same draws in the same order).
#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

// Government type codes: 0 = none (initial status quo), 1 = single-party
// majority, 2 = single-party minority, 3 = minority coalition,
// 4 = minimal winning coalition, 5 = surplus majority, 6 = unity (all parties,
// surplus).
static int classify(const std::vector<int>& in_gov, const NumericVector& seats,
                    double quota) {
  int m = seats.size(), k = 0;
  double tot = 0.0;
  for (int j = 0; j < m; ++j) if (in_gov[j]) { ++k; tot += seats[j]; }
  if (k == 0) return 0;
  bool maj = tot > quota;
  if (k == 1) return maj ? 1 : 2;
  if (!maj) return 3;
  for (int j = 0; j < m; ++j)                 // any non-pivotal member?
    if (in_gov[j] && tot - seats[j] > quota) return (k == m) ? 6 : 5;
  return 4;
}

// Arguments: z (m x 2 ideal points), seats (seat shares), beta (policy weight),
// n_form formations, max_periods before caretaker reinstalled, half_width of
// the policy space, quota for a majority, initial status-quo policy, and
// history = TRUE to return per-formation running Gamson coefficients and
// cabinet shares.
// [[Rcpp::export]]
List ggs_run_cpp(NumericMatrix z, NumericVector seats, double beta,
                 int n_form = 2500, int max_periods = 100,
                 double half_width = 3.0, double quota = 0.5,
                 NumericVector sq_policy = NumericVector::create(0.0, 0.0),
                 bool history = false,
                 int cab_method = 0, double policy_half_width = -1.0) {
  // cab_method: 0 = uniform on the simplex (default, as in the paper);
  //             1 = normalised iid uniforms (sensitivity variant).
  // policy_half_width: range of proposed policies; < 0 means same as half_width.
  const double phw = policy_half_width < 0 ? half_width : policy_half_width;
  const int m = seats.size();
  // status quo
  std::vector<double> sq_cab(m, 0.0), u_sq(m);
  double sq_x = sq_policy[0], sq_y = sq_policy[1];
  int sq_form = -1;                          // formateur of SQ (0-based), -1 none
  // proposals
  std::vector<int> mem(m * m);
  std::vector<double> cab(m * m), px(m), py(m), u(m);
  std::vector<int> qual(m), sq_in(m, 0);

  IntegerVector delay(n_form), type(n_form), formateur(n_form), n_gov(n_form);
  // running Gamson regression sums (Figure 1)
  double Sx = 0, Sy = 0, Sxx = 0, Sxy = 0, Sn = 0;
  NumericVector run_slope(history ? n_form : 0), run_int(history ? n_form : 0);
  NumericMatrix cab_hist(history ? n_form : 0, m);

  for (int f = 0; f < n_form; ++f) {
    for (int j = 0; j < m; ++j) {
      double dx = sq_x - z(j, 0), dy = sq_y - z(j, 1);
      u_sq[j] = sq_cab[j] - beta * (dx * dx + dy * dy);
    }
    int chosen = -1, t = 0;
    for (t = 1; t <= max_periods; ++t) {
      // 1-2. every party draws a proposal
      for (int i = 0; i < m; ++i) {
        int* M = &mem[i * m];
        double* C = &cab[i * m];
        for (int j = 0; j < m; ++j) {
          if (j == i) M[j] = 1;
          else M[j] = unif_rand() < 0.5;
        }
        double s = 0.0;
        for (int j = 0; j < m; ++j) {
          if (M[j]) { C[j] = cab_method == 0 ? -std::log(unif_rand()) : unif_rand(); s += C[j]; }
          else C[j] = 0.0;
        }
        for (int j = 0; j < m; ++j) C[j] /= s;
        px[i] = -phw + 2.0 * phw * unif_rand();
        py[i] = -phw + 2.0 * phw * unif_rand();
      }
      // 3. compare to status quo
      int nq = 0;
      for (int i = 0; i < m; ++i) {
        const int* M = &mem[i * m];
        const double* C = &cab[i * m];
        bool ok = true;
        double support = 0.0;
        for (int j = 0; j < m; ++j) {
          double dx = px[i] - z(j, 0), dy = py[i] - z(j, 1);
          double uj = C[j] - beta * (dx * dx + dy * dy);
          bool pref = uj > u_sq[j];
          if (M[j] && !pref) { ok = false; break; }
          if (pref) support += seats[j];
        }
        if (ok && support > quota) qual[nq++] = i;
      }
      // 4. head of state installs one at random
      if (nq > 0) {
        chosen = (nq == 1) ? qual[0] : qual[(int)(unif_rand() * nq)];
        break;
      }
    }
    if (t > max_periods) t = max_periods;    // 5. caretaker reinstalled
    if (chosen >= 0) {
      for (int j = 0; j < m; ++j) {
        sq_cab[j] = cab[chosen * m + j];
        sq_in[j] = mem[chosen * m + j];
      }
      sq_x = px[chosen]; sq_y = py[chosen]; sq_form = chosen;
    }
    delay[f] = t;
    formateur[f] = sq_form + 1;
    type[f] = classify(sq_in, seats, quota);
    int k = 0; double gs = 0.0;
    for (int j = 0; j < m; ++j) if (sq_in[j]) { ++k; gs += seats[j]; }
    n_gov[f] = k;
    if (history) {
      for (int j = 0; j < m; ++j) cab_hist(f, j) = sq_cab[j];
      for (int j = 0; j < m; ++j) if (sq_in[j]) {
        double x = seats[j] / gs, y = sq_cab[j];
        Sx += x; Sy += y; Sxx += x * x; Sxy += x * y; Sn += 1;
      }
      double den = Sn * Sxx - Sx * Sx;
      run_slope[f] = den > 1e-12 ? (Sn * Sxy - Sx * Sy) / den : NA_REAL;
      run_int[f] = den > 1e-12 ? (Sy - run_slope[f] * Sx) / Sn : NA_REAL;
    }
  }
  NumericVector fcab(m);
  IntegerVector fin(m);
  for (int j = 0; j < m; ++j) { fcab[j] = sq_cab[j]; fin[j] = sq_in[j]; }
  List out = List::create(_["delay"] = delay, _["type"] = type,
                          _["formateur"] = formateur, _["n_gov"] = n_gov,
                          _["final_cab"] = fcab, _["final_in"] = fin,
                          _["final_policy"] = NumericVector::create(sq_x, sq_y));
  if (history) { out["run_slope"] = run_slope; out["run_int"] = run_int; out["cab_hist"] = cab_hist; }
  return out;
}
