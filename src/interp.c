/*
 * interp.c — faithful port of FLOWMath's Akima spline (with its endpoint
 * slope estimation and clamped-index extrapolation), recursive 2D/3D
 * interpolation, and trapezoidal integration.
 */
#include <math.h>
#include <stdlib.h>

#include "ccblade.h"

/* index i such that xvec[i] <= x < xvec[i+1] (0-based), clamped to
 * [0, n-2] to allow extrapolation via the end segments (FLOWMath.findindex) */
static int findindex(const double* xvec, int n, double x) {
  int lo = 0, hi = n - 1;
  /* largest i with xvec[i] <= x (searchsortedlast) */
  if (x < xvec[0]) return 0;
  while (hi - lo > 1) {
    int mid = (lo + hi) / 2;
    if (xvec[mid] <= x)
      lo = mid;
    else
      hi = mid;
  }
  if (lo >= n - 1) lo = n - 2;
  return lo;
}

CCAkima* cc_akima(const double* x, const double* y, int n) {
  if (n < 2) return NULL;

  CCAkima* s = (CCAkima*)malloc(sizeof(CCAkima));
  if (!s) return NULL;
  s->n = n;
  s->xdata = (double*)malloc(sizeof(double) * n);
  s->p0 = (double*)malloc(sizeof(double) * (n - 1));
  s->p1 = (double*)malloc(sizeof(double) * (n - 1));
  s->p2 = (double*)malloc(sizeof(double) * (n - 1));
  s->p3 = (double*)malloc(sizeof(double) * (n - 1));
  if (!s->xdata || !s->p0 || !s->p1 || !s->p2 || !s->p3) {
    cc_akima_free(s);
    return NULL;
  }
  for (int i = 0; i < n; i++) s->xdata[i] = x[i];

  /* segment slopes with extended endpoints; m_ext[k] corresponds to
   * Julia's m[k-2] (Julia indices run -1..n+1). Zero-initialized like
   * FLOWMath's OffsetVector(zeros(...)) so the n=2 endpoint estimates
   * match Julia exactly. */
  double* m = (double*)calloc(n + 4, sizeof(double));
  if (!m) {
    cc_akima_free(s);
    return NULL;
  }
  for (int i = 1; i <= n - 1; i++)
    m[i + 2] = (y[i] - y[i - 1]) / (x[i] - x[i - 1]);
  /* Julia: m[0] = 2*m[1] - m[2]; m[-1] = 2*m[0] - m[1];
   *        m[n] = 2*m[n-1] - m[n-2]; m[n+1] = 2*m[n] - m[n-1] */
  m[2] = 2.0 * m[3] - m[4];             /* Julia m[0]  */
  m[1] = 2.0 * m[2] - m[3];             /* Julia m[-1] */
  m[n + 2] = 2.0 * m[n + 1] - m[n];     /* Julia m[n]  */
  m[n + 3] = 2.0 * m[n + 2] - m[n + 1]; /* Julia m[n+1] */

  const double eps = 1e-30;
  double* t = (double*)malloc(sizeof(double) * (n + 1));
  if (!t) {
    free(m);
    cc_akima_free(s);
    return NULL;
  }
  for (int i = 1; i <= n; i++) {
    /* Julia m[i-2], m[i-1], m[i], m[i+1] -> here m[i], m[i+1], m[i+2], m[i+3]
     */
    double m1 = m[i], m2 = m[i + 1], m3 = m[i + 2], m4 = m[i + 3];
    double w1 = fabs(m4 - m3);
    double w2 = fabs(m2 - m1);
    if (w1 < eps && w2 < eps)
      t[i] = 0.5 * (m2 + m3);
    else
      t[i] = (w1 * m2 + w2 * m3) / (w1 + w2);
  }

  for (int i = 1; i <= n - 1; i++) {
    double dx = x[i] - x[i - 1];
    double t1 = t[i], t2 = t[i + 1];
    double mi = m[i + 2];
    s->p0[i - 1] = y[i - 1];
    s->p1[i - 1] = t1;
    s->p2[i - 1] = (3.0 * mi - 2.0 * t1 - t2) / dx;
    s->p3[i - 1] = (t1 + t2 - 2.0 * mi) / (dx * dx);
  }

  free(m);
  free(t);
  return s;
}

double cc_akima_eval(const CCAkima* s, double x) {
  int j = findindex(s->xdata, s->n, x);
  double dx = x - s->xdata[j];
  return s->p0[j] + s->p1[j] * dx + s->p2[j] * dx * dx +
         s->p3[j] * dx * dx * dx;
}

double cc_akima_deriv(const CCAkima* s, double x) {
  int j = findindex(s->xdata, s->n, x);
  double dx = x - s->xdata[j];
  return s->p1[j] + 2.0 * s->p2[j] * dx + 3.0 * s->p3[j] * dx * dx;
}

void cc_akima_free(CCAkima* s) {
  if (!s) return;
  free(s->xdata);
  free(s->p0);
  free(s->p1);
  free(s->p2);
  free(s->p3);
  free(s);
}

/* one-shot akima: build, evaluate, free */
static double akima1(const double* x, const double* y, int n, double xq) {
  CCAkima* s = cc_akima(x, y, n);
  if (!s) return NAN;
  double v = cc_akima_eval(s, xq);
  cc_akima_free(s);
  return v;
}

double cc_interp2d_akima(const double* x, const double* y, const double* f2,
                         int nx, int ny, double xq, double yq) {
  /* FLOWMath.interp2d: for each j over y, interp over x; then interp the
   * results over y. */
  double* tmp = (double*)malloc(sizeof(double) * nx);
  double* yinterp = (double*)malloc(sizeof(double) * ny);
  if (!tmp || !yinterp) {
    free(tmp);
    free(yinterp);
    return NAN;
  }
  for (int j = 0; j < ny; j++) {
    for (int i = 0; i < nx; i++) tmp[i] = f2[i + nx * j];
    yinterp[j] = akima1(x, tmp, nx, xq);
  }
  double out = akima1(y, yinterp, ny, yq);
  free(tmp);
  free(yinterp);
  return out;
}

double cc_interp3d_akima(const double* x, const double* y, const double* z,
                         const double* f3, int nx, int ny, int nz, double xq,
                         double yq, double zq) {
  /* FLOWMath.interp3d: for each k over z, interp2d over (x, y); then
   * interp the results over z. */
  double* f2 = (double*)malloc(sizeof(double) * nx * ny);
  double* zinterp = (double*)malloc(sizeof(double) * nz);
  if (!f2 || !zinterp) {
    free(f2);
    free(zinterp);
    return NAN;
  }
  for (int k = 0; k < nz; k++) {
    for (int j = 0; j < ny; j++)
      for (int i = 0; i < nx; i++) f2[i + nx * j] = f3[i + nx * (j + ny * k)];
    zinterp[k] = cc_interp2d_akima(x, y, f2, nx, ny, xq, yq);
  }
  double out = akima1(z, zinterp, nz, zq);
  free(f2);
  free(zinterp);
  return out;
}

double cc_trapz(const double* x, const double* y, int n) {
  double acc = 0.0;
  for (int i = 0; i < n - 1; i++) acc += (x[i + 1] - x[i]) * (y[i + 1] + y[i]);
  return 0.5 * acc;
}
