// =============================================================================
// src/gfe_exact_bb.cpp  --  compiled inner loop of the exact GFE algorithm
//
// C++ (Rcpp) version of one stage of the repetitive branch and bound in
// R/gfe_exact.R, which itself translates Bonhomme and Manresa's
// Exact algorithm/Brusco_cov.f90 (their extension of Brusco, 2006, to
// covariates). The authors ran this algorithm as compiled Fortran; R is too
// slow for N much above 20 because the number of nodes grows quickly, so this
// file plays the role of their Fortran. The logic, step numbers, and variable
// names follow the R function line by line.
//
// Compiled on first use by Rcpp::sourceCpp() (Windows users need Rtools).
//
// R translation Copyright (C) 2026 Eliot R. Weinstein. Licensed under the GNU
// GPL, version 2 or (at your option) any later version; see LICENSE. The
// license covers the translation only, not Bonhomme and Manresa's method or
// original code.
// =============================================================================
#include <Rcpp.h>
#include <vector>
#include <cmath>
using namespace Rcpp;

// Concentrated partial objective a - b' C^{-1} b from v = (a, b, vec(C)).
// Returns NAN if C is (numerically) singular.
static double partial_obj(const std::vector<double>& v, int K) {
  if (K == 0) return v[0];
  if (K == 1) return (v[2] > 1e-12) ? v[0] - v[1] * v[1] / v[2] : NAN;
  if (K == 2) {
    double c11 = v[3], c21 = v[4], c22 = v[6];
    double dt = c11 * c22 - c21 * c21;
    if (!(dt > 1e-12 * c11 * c22)) return NAN;
    return v[0] - (c22 * v[1] * v[1] - 2.0 * c21 * v[1] * v[2] + c11 * v[2] * v[2]) / dt;
  }
  // General K: Gaussian elimination with partial pivoting on [C | b].
  std::vector<double> C(K * K), b(K), sol(K);
  for (int k = 0; k < K; k++) {
    b[k] = v[1 + k];
    for (int l = 0; l < K; l++) C[k + K * l] = v[1 + K + l * K + k];
  }
  std::vector<double> Cw = C, bw = b;
  double scale = 0; for (int k = 0; k < K; k++) scale = std::max(scale, std::fabs(C[k + K * k]));
  for (int c = 0; c < K; c++) {
    int piv = c; double best = std::fabs(Cw[c + K * c]);
    for (int r = c + 1; r < K; r++) if (std::fabs(Cw[r + K * c]) > best) { best = std::fabs(Cw[r + K * c]); piv = r; }
    if (best <= 1e-12 * scale) return NAN;
    if (piv != c) { for (int l = 0; l < K; l++) std::swap(Cw[c + K * l], Cw[piv + K * l]); std::swap(bw[c], bw[piv]); }
    for (int r = c + 1; r < K; r++) {
      double f = Cw[r + K * c] / Cw[c + K * c];
      for (int l = c; l < K; l++) Cw[r + K * l] -= f * Cw[c + K * l];
      bw[r] -= f * bw[c];
    }
  }
  for (int r = K - 1; r >= 0; r--) {
    double s = bw[r];
    for (int l = r + 1; l < K; l++) s -= Cw[r + K * l] * sol[l];
    sol[r] = s / Cw[r + K * r];
  }
  double q = 0; for (int k = 0; k < K; k++) q += b[k] * sol[k];
  return v[0] - q;
}

// [[Rcpp::export]]
List bb_stage_cpp(NumericVector PWsub, int E, int G, int K, double Z,
                  NumericVector opt) {
  const int Qd = 1 + K + K * K;
  // PWsub is an E x E x Qd array (column-major): PW(i, j, q).
  auto PW = [&](int i, int j, int q) { return PWsub[i + E * j + E * E * q]; };

  int P = 0, Q = G, M = 1;               // 1-based object depth and cluster
  std::vector<int> n(G + 1, 0), x(E + 1, 0), XB;
  std::vector<double> sums((G + 1) * Qd, 0.0), v(Qd);
  double nodes = 0;
  bool have = false;

  auto add_obj = [&](int P, int M, double sign) {
    // add/remove pairs (P, j), j < P, j in cluster M  (SUMA, SUMB, SUMC)
    for (int j = 1; j < P; j++) if (x[j] == M)
      for (int q = 0; q < Qd; q++) sums[M * Qd + q] += sign * PW(P - 1, j - 1, q);
  };

  int step = 1;
  double ZB = NAN;
  while (true) {
    if (step == 1) {                       // STEP 1: increment search depth
      P++; M = 1;
      n[M]++; if (n[M] == 1) Q--;
      x[P] = M;
      add_obj(P, M, +1.0);
      step = 2;
    } else if (step == 2) {                // STEP 2: feasibility
      step = (E - P < Q) ? 5 : 3;
    } else if (step == 3) {                // STEP 3: bound test
      nodes++;
      int used = G - Q;
      for (int q = 0; q < Qd; q++) {
        double s = 0;
        for (int m = 1; m <= used; m++) s += sums[m * Qd + q] / n[m];
        v[q] = s;
      }
      ZB = partial_obj(v, K);
      if (!std::isnan(ZB) && Z <= ZB + opt[E - P]) { step = 5; continue; }
      if (P != E) { step = 1; continue; }  // STEP 4: update incumbent
      if (!std::isnan(ZB)) { Z = ZB; XB.assign(x.begin() + 1, x.end()); have = true; }
      step = 5;
    } else if (step == 5) {                // STEP 5: action after fathoming
      step = (M == G || n[M] == 1) ? 7 : 6;
    } else if (step == 6) {                // STEP 6: branch right
      add_obj(P, M, -1.0);
      n[M]--;
      M++;
      n[M]++; if (n[M] == 1) Q--;
      x[P] = M;
      add_obj(P, M, +1.0);
      step = 2;
    } else {                               // STEP 7: retract depth
      x[P] = 0;
      add_obj(P, M, -1.0);
      n[M]--;
      P--;
      if (n[M] == 0) Q++;
      if (P > 0) { M = x[P]; step = 5; } else break;
    }
    if (((long long) nodes & 0xFFFFF) == 0) Rcpp::checkUserInterrupt();
  }
  return List::create(_["Z"] = Z, _["found"] = have,
                      _["XB"] = have ? IntegerVector(XB.begin(), XB.end()) : IntegerVector(0),
                      _["nodes"] = nodes);
}
