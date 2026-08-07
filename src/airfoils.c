/*
 * airfoils.c — port of CCBlade.jl's airfoil data handling: file parsing,
 * AFType constructors/evaluation, Mach/Reynolds/rotation corrections,
 * tip-loss factors, linear lift-curve extraction, and Viterna extrapolation.
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ccblade.h"

/* ---------------- file parsing ---------------- */

/* Port of CCBlade.parsefile. Format:
 *   line 1: info (skipped)
 *   line 2: Re
 *   line 3: Mach
 *   remaining lines: alpha cl cd [ignored...]
 * If radians==0, alpha is converted from degrees.
 * On success returns 0 and fills malloc'd arrays (caller frees). */
static int parsefile(const char* filename, int radians, char** info, double* Re,
                     double* Mach, double** alpha, double** cl, double** cd,
                     int* n) {
  FILE* f = fopen(filename, "r");
  if (!f) return -1;

  char line[4096];
  if (!fgets(line, sizeof line, f)) {
    fclose(f);
    return -1;
  }
  line[strcspn(line, "\r\n")] = '\0';
  char* info_copy = (char*)malloc(strlen(line) + 1);
  if (!info_copy) {
    fclose(f);
    return -1;
  }
  strcpy(info_copy, line);

  if (!fgets(line, sizeof line, f)) {
    free(info_copy);
    fclose(f);
    return -1;
  }
  *Re = strtod(line, NULL);
  if (!fgets(line, sizeof line, f)) {
    free(info_copy);
    fclose(f);
    return -1;
  }
  *Mach = strtod(line, NULL);

  int cap = 64, cnt = 0;
  double* a = (double*)malloc(sizeof(double) * cap);
  double* c1 = (double*)malloc(sizeof(double) * cap);
  double* c2 = (double*)malloc(sizeof(double) * cap);
  if (!a || !c1 || !c2) {
    free(info_copy);
    free(a);
    free(c1);
    free(c2);
    fclose(f);
    return -1;
  }

  while (fgets(line, sizeof line, f)) {
    char* p = line;
    char* end;
    double va = strtod(p, &end);
    if (end == p) continue; /* blank line */
    p = end;
    double vc1 = strtod(p, &end);
    if (end == p) continue;
    p = end;
    double vc2 = strtod(p, &end);
    if (end == p) continue;
    if (cnt == cap) {
      cap *= 2;
      double* na = (double*)realloc(a, sizeof(double) * cap);
      double* nc1 = (double*)realloc(c1, sizeof(double) * cap);
      double* nc2 = (double*)realloc(c2, sizeof(double) * cap);
      if (!na || !nc1 || !nc2) {
        free(na ? na : a);
        free(nc1 ? nc1 : c1);
        free(nc2 ? nc2 : c2);
        free(info_copy);
        fclose(f);
        return -1;
      }
      a = na;
      c1 = nc1;
      c2 = nc2;
    }
    a[cnt] = va;
    c1[cnt] = vc1;
    c2[cnt] = vc2;
    cnt++;
  }
  fclose(f);

  if (!radians) {
    for (int i = 0; i < cnt; i++) a[i] *= CC_PI / 180.0;
  }

  *info = info_copy;
  *alpha = a;
  *cl = c1;
  *cd = c2;
  *n = cnt;
  return 0;
}

/* Port of CCBlade.writefile. Returns 0 on success. */
static int writefile(const char* filename, const char* info, double Re,
                     double Mach, const double* alpha, const double* cl,
                     const double* cd, int n, int radians) {
  FILE* f = fopen(filename, "w");
  if (!f) return -1;
  fprintf(f, "%s\n", info);
  fprintf(f, "%.17g\n", Re);
  fprintf(f, "%.17g\n", Mach);
  double factor = radians ? 1.0 : 180.0 / CC_PI;
  for (int i = 0; i < n; i++)
    fprintf(f, "%.17g\t%.17g\t%.17g\n", alpha[i] * factor, cl[i], cd[i]);
  fclose(f);
  return 0;
}

/* ---------------- helpers ---------------- */

static char* cc_strdup(const char* s) {
  size_t n = strlen(s) + 1;
  char* p = (char*)malloc(n);
  if (p) memcpy(p, s, n);
  return p;
}

static double* cc_copy(const double* src, int n) {
  double* p = (double*)malloc(sizeof(double) * n);
  if (p) memcpy(p, src, sizeof(double) * n);
  return p;
}

static CCAirfoil* af_alloc(CCAFKind kind) {
  CCAirfoil* af = (CCAirfoil*)calloc(1, sizeof(CCAirfoil));
  if (af) af->kind = kind;
  return af;
}

/* ---------------- constructors ---------------- */

// cppcheck-suppress funcArgNamesDifferentUnnamed
CCAirfoil* cc_af_function(CCAFFunc f, void* userdata) {
  CCAirfoil* af = af_alloc(CC_AF_FUNCTION);
  if (af) {
    af->func = f;
    af->userdata = userdata;
  }
  return af;
}

CCAirfoil* cc_af_simple(double m, double alpha0, double clmax, double clmin,
                        double cd0, double cd2) {
  CCAirfoil* af = af_alloc(CC_AF_SIMPLE);
  if (af) {
    af->m = m;
    af->alpha0 = alpha0;
    af->clmax = clmax;
    af->clmin = clmin;
    af->cd0 = cd0;
    af->cd2 = cd2;
  }
  return af;
}

CCAirfoil* cc_af_alpha(const double* alpha, const double* cl, const double* cd,
                       int n, const char* info, double re, double mach) {
  CCAirfoil* af = af_alloc(CC_AF_ALPHA);
  if (!af) return NULL;
  af->nalpha = n;
  af->nre = af->nmach = 0;
  af->alpha = cc_copy(alpha, n);
  af->cl = cc_copy(cl, n);
  af->cd = cc_copy(cd, n);
  af->info = cc_strdup(info ? info : "");
  af->re_info = re;
  af->mach_info = mach;
  af->clspline = cc_akima(alpha, cl, n);
  af->cdspline = cc_akima(alpha, cd, n);
  if (!af->alpha || !af->cl || !af->cd || !af->info || !af->clspline ||
      !af->cdspline) {
    cc_af_free(af);
    return NULL;
  }
  return af;
}

CCAirfoil* cc_af_alpha_file(const char* filename, int radians) {
  char* info;
  double re, mach, *alpha, *cl, *cd;
  int n;
  if (parsefile(filename, radians, &info, &re, &mach, &alpha, &cl, &cd, &n) !=
      0)
    return NULL;
  CCAirfoil* af = cc_af_alpha(alpha, cl, cd, n, info, re, mach);
  free(info);
  free(alpha);
  free(cl);
  free(cd);
  return af;
}

/* shared constructor for the 2D/3D table types */
static CCAirfoil* af_table_alloc(CCAFKind kind, const double* alpha, int nalpha,
                                 const double* re, int nre, const double* mach,
                                 int nmach, const double* cl, const double* cd,
                                 const char* info) {
  CCAirfoil* af = af_alloc(kind);
  if (!af) return NULL;
  af->nalpha = nalpha;
  af->nre = nre;
  af->nmach = nmach;
  af->alpha = cc_copy(alpha, nalpha);
  af->re = nre > 0 ? cc_copy(re, nre) : NULL;
  af->mach = nmach > 0 ? cc_copy(mach, nmach) : NULL;
  long total = (long)nalpha * (nre > 0 ? nre : 1) * (nmach > 0 ? nmach : 1);
  af->cl = cc_copy(cl, total);
  af->cd = cc_copy(cd, total);
  af->info = cc_strdup(info ? info : "");
  if (!af->alpha || (nre > 0 && !af->re) || (nmach > 0 && !af->mach) ||
      !af->cl || !af->cd || !af->info) {
    cc_af_free(af);
    return NULL;
  }
  return af;
}

CCAirfoil* cc_af_alpha_re(const double* alpha, const double* re,
                          const double* cl, const double* cd, int nalpha,
                          int nre, const char* info, double mach) {
  CCAirfoil* af = af_table_alloc(CC_AF_ALPHA_RE, alpha, nalpha, re, nre, NULL,
                                 0, cl, cd, info);
  if (af) af->mach_info = mach;
  return af;
}

CCAirfoil* cc_af_alpha_mach(const double* alpha, const double* mach,
                            const double* cl, const double* cd, int nalpha,
                            int nmach, const char* info, double re) {
  CCAirfoil* af = af_table_alloc(CC_AF_ALPHA_MACH, alpha, nalpha, NULL, 0, mach,
                                 nmach, cl, cd, info);
  if (af) af->re_info = re;
  return af;
}

CCAirfoil* cc_af_alpha_re_mach(const double* alpha, const double* re,
                               const double* mach, const double* cl,
                               const double* cd, int nalpha, int nre, int nmach,
                               const char* info) {
  return af_table_alloc(CC_AF_ALPHA_RE_MACH, alpha, nalpha, re, nre, mach,
                        nmach, cl, cd, info);
}

/* one file per condition (Re or Mach); common alpha grid assumed */
static CCAirfoil* af_files_2d(CCAFKind kind, const char** filenames, int nfiles,
                              int radians) {
  char* info;
  double re0, mach0, *alpha, *cl0, *cd0;
  int nalpha;
  if (parsefile(filenames[0], radians, &info, &re0, &mach0, &alpha, &cl0, &cd0,
                &nalpha) != 0)
    return NULL;

  double* cond = (double*)malloc(sizeof(double) * nfiles);
  double* cl = (double*)malloc(sizeof(double) * nalpha * nfiles);
  double* cd = (double*)malloc(sizeof(double) * nalpha * nfiles);
  if (!cond || !cl || !cd) {
    free(cond);
    free(cl);
    free(cd);
    free(info);
    free(alpha);
    free(cl0);
    free(cd0);
    return NULL;
  }

  cond[0] = (kind == CC_AF_ALPHA_RE) ? re0 : mach0;
  memcpy(cl, cl0, sizeof(double) * nalpha);
  memcpy(cd, cd0, sizeof(double) * nalpha);
  free(cl0);
  free(cd0);

  int ok = 1;
  for (int j = 1; j < nfiles && ok; j++) {
    char* infoj;
    double rej, machj, *alphaj, *clj, *cdj;
    int nj;
    if (parsefile(filenames[j], radians, &infoj, &rej, &machj, &alphaj, &clj,
                  &cdj, &nj) != 0) {
      ok = 0;
      break;
    }
    if (nj != nalpha) {
      free(infoj);
      free(alphaj);
      free(clj);
      free(cdj);
      ok = 0;
      break;
    }
    cond[j] = (kind == CC_AF_ALPHA_RE) ? rej : machj;
    memcpy(cl + (size_t)nalpha * j, clj, sizeof(double) * nalpha);
    memcpy(cd + (size_t)nalpha * j, cdj, sizeof(double) * nalpha);
    free(infoj);
    free(alphaj);
    free(clj);
    free(cdj);
  }

  CCAirfoil* af = NULL;
  if (ok) {
    if (kind == CC_AF_ALPHA_RE)
      af = cc_af_alpha_re(alpha, cond, cl, cd, nalpha, nfiles, info, mach0);
    else
      af = cc_af_alpha_mach(alpha, cond, cl, cd, nalpha, nfiles, info, re0);
  }
  free(cond);
  free(cl);
  free(cd);
  free(info);
  free(alpha);
  return af;
}

CCAirfoil* cc_af_alpha_re_files(const char** filenames, int nfiles,
                                int radians) {
  return af_files_2d(CC_AF_ALPHA_RE, filenames, nfiles, radians);
}

CCAirfoil* cc_af_alpha_mach_files(const char** filenames, int nfiles,
                                  int radians) {
  return af_files_2d(CC_AF_ALPHA_MACH, filenames, nfiles, radians);
}

/* filenames[i + nre*j] corresponds to Re[i], Mach[j] (Julia: filenames[i, j])
 */
CCAirfoil* cc_af_alpha_re_mach_files(const char** filenames, int nre, int nmach,
                                     int radians) {
  char* info;
  double re0, mach0, *alpha, *cl0, *cd0;
  int nalpha;
  if (parsefile(filenames[0], radians, &info, &re0, &mach0, &alpha, &cl0, &cd0,
                &nalpha) != 0)
    return NULL;
  free(cl0);
  free(cd0);

  double* re = (double*)malloc(sizeof(double) * nre);
  double* mach = (double*)malloc(sizeof(double) * nmach);
  double* cl = (double*)malloc(sizeof(double) * nalpha * nre * nmach);
  double* cd = (double*)malloc(sizeof(double) * nalpha * nre * nmach);
  if (!re || !mach || !cl || !cd) {
    free(re);
    free(mach);
    free(cl);
    free(cd);
    free(info);
    free(alpha);
    return NULL;
  }

  int ok = 1;
  for (int j = 0; j < nmach && ok; j++) {
    for (int i = 0; i < nre && ok; i++) {
      char* infoij;
      double reij, machij, *alphaij, *clij, *cdij;
      int nij;
      if (parsefile(filenames[i + nre * j], radians, &infoij, &reij, &machij,
                    &alphaij, &clij, &cdij, &nij) != 0) {
        ok = 0;
        break;
      }
      if (nij != nalpha) {
        free(infoij);
        free(alphaij);
        free(clij);
        free(cdij);
        ok = 0;
        break;
      }
      re[i] = reij;
      mach[j] = machij;
      for (int k = 0; k < nalpha; k++) {
        cl[k + nalpha * (i + nre * j)] = clij[k];
        cd[k + nalpha * (i + nre * j)] = cdij[k];
      }
      free(infoij);
      free(alphaij);
      free(clij);
      free(cdij);
    }
  }

  CCAirfoil* af = NULL;
  if (ok)
    af = cc_af_alpha_re_mach(alpha, re, mach, cl, cd, nalpha, nre, nmach, info);
  free(re);
  free(mach);
  free(cl);
  free(cd);
  free(info);
  free(alpha);
  return af;
}

void cc_af_free(CCAirfoil* af) {
  if (!af) return;
  free(af->alpha);
  free(af->re);
  free(af->mach);
  free(af->cl);
  free(af->cd);
  free(af->info);
  cc_akima_free(af->clspline);
  cc_akima_free(af->cdspline);
  free(af);
}

/* ---------------- evaluation ---------------- */

void cc_afeval(const CCAirfoil* af, double alpha, double Re, double Mach,
               double* cl, double* cd) {
  switch (af->kind) {
    case CC_AF_FUNCTION:
      af->func(alpha, Re, Mach, af->userdata, cl, cd);
      return;
    case CC_AF_SIMPLE: {
      double c = af->m * (alpha - af->alpha0);
      if (c > af->clmax) c = af->clmax;
      if (c < af->clmin) c = af->clmin;
      *cl = c;
      *cd = af->cd0 + af->cd2 * c * c;
      return;
    }
    case CC_AF_ALPHA:
      *cl = cc_akima_eval(af->clspline, alpha);
      *cd = cc_akima_eval(af->cdspline, alpha);
      return;
    case CC_AF_ALPHA_RE:
      *cl = cc_interp2d_akima(af->alpha, af->re, af->cl, af->nalpha, af->nre,
                              alpha, Re);
      *cd = cc_interp2d_akima(af->alpha, af->re, af->cd, af->nalpha, af->nre,
                              alpha, Re);
      return;
    case CC_AF_ALPHA_MACH:
      *cl = cc_interp2d_akima(af->alpha, af->mach, af->cl, af->nalpha,
                              af->nmach, alpha, Mach);
      *cd = cc_interp2d_akima(af->alpha, af->mach, af->cd, af->nalpha,
                              af->nmach, alpha, Mach);
      return;
    case CC_AF_ALPHA_RE_MACH:
      *cl = cc_interp3d_akima(af->alpha, af->re, af->mach, af->cl, af->nalpha,
                              af->nre, af->nmach, alpha, Re, Mach);
      *cd = cc_interp3d_akima(af->alpha, af->re, af->mach, af->cd, af->nalpha,
                              af->nre, af->nmach, alpha, Re, Mach);
      return;
  }
  *cl = NAN;
  *cd = NAN;
}

/* ---------------- write ---------------- */

int cc_write_af(const char* filename, const CCAirfoil* af, int radians) {
  if (af->kind != CC_AF_ALPHA) return -1;
  return writefile(filename, af->info, af->re_info, af->mach_info, af->alpha,
                   af->cl, af->cd, af->nalpha, radians);
}

int cc_write_af_multi(const char** filenames, const CCAirfoil* af,
                      int radians) {
  if (af->kind == CC_AF_ALPHA) return cc_write_af(filenames[0], af, radians);

  if (af->kind == CC_AF_ALPHA_RE || af->kind == CC_AF_ALPHA_MACH) {
    int ncond = af->kind == CC_AF_ALPHA_RE ? af->nre : af->nmach;
    for (int j = 0; j < ncond; j++) {
      double re = af->kind == CC_AF_ALPHA_RE ? af->re[j] : af->re_info;
      double mach = af->kind == CC_AF_ALPHA_RE ? af->mach_info : af->mach[j];
      if (writefile(filenames[j], af->info, re, mach, af->alpha,
                    af->cl + (size_t)af->nalpha * j,
                    af->cd + (size_t)af->nalpha * j, af->nalpha, radians) != 0)
        return -1;
    }
    return 0;
  }

  if (af->kind == CC_AF_ALPHA_RE_MACH) {
    for (int j = 0; j < af->nmach; j++) {
      for (int i = 0; i < af->nre; i++) {
        const double* clij = af->cl + (size_t)af->nalpha * (i + af->nre * j);
        const double* cdij = af->cd + (size_t)af->nalpha * (i + af->nre * j);
        if (writefile(filenames[i + af->nre * j], af->info, af->re[i],
                      af->mach[j], af->alpha, clij, cdij, af->nalpha,
                      radians) != 0)
          return -1;
      }
    }
    return 0;
  }
  return -1;
}

/* ---------------- corrections ---------------- */

void cc_mach_correction_prandtl_glauert(double cl, double cd, double Mach,
                                        double* cl_out, double* cd_out) {
  double beta = sqrt(1.0 - Mach * Mach);
  *cl_out = cl / beta;
  *cd_out = cd;
}

void cc_re_correction_skin_friction(double re0, double p, double cl, double cd,
                                    double Re, double* cl_out, double* cd_out) {
  *cl_out = cl;
  *cd_out = cd * pow(re0 / Re, p);
}

void cc_rotation_correction_duselig_eggers(double a, double b, double d,
                                           double m, double alpha0, double cl,
                                           double cd, double cr, double rR,
                                           double tsr, double alpha, double phi,
                                           double alpha_max_corr,
                                           double* cl_out, double* cd_out) {
  /* Du-Selig correction for lift */
  double Lambda = tsr / sqrt(1.0 + tsr * tsr);
  double expon = d / (Lambda * rR);
  double fcl =
      1.0 / m *
      (1.6 * cr / 0.1267 * (a - pow(cr, expon)) / (b + pow(cr, expon)) - 1.0);

  double cl_linear = m * (alpha - alpha0);

  double amax = atan(1.0 / 0.12) - 5.0 * CC_PI / 180.0;
  double adj;
  double aa = fabs(alpha);
  if (aa >= amax)
    adj = 0.0;
  else if (aa > alpha_max_corr)
    adj = ((amax - aa) / (amax - alpha_max_corr)) *
          ((amax - aa) / (amax - alpha_max_corr));
  else
    adj = 1.0;

  double deltacl = adj * fcl * (cl_linear - cl);
  cl += deltacl;

  /* Eggers correction for drag */
  double deltacd =
      deltacl * (sin(phi) - 0.12 * cos(phi)) / (cos(phi) + 0.12 * sin(phi));
  cd += deltacd;

  *cl_out = cl;
  *cd_out = cd;
}

double cc_tip_correction_prandtl(double r, double Rtip, double phi, int B) {
  double asphi = fabs(sin(phi));
  double factortip = B / 2.0 * (Rtip / r - 1.0) / asphi;
  return 2.0 / CC_PI * acos(exp(-factortip));
}

double cc_tip_correction_prandtl_hub(double r, double Rhub, double Rtip,
                                     double phi, int B) {
  double asphi = fabs(sin(phi));
  double factortip = B / 2.0 * (Rtip / r - 1.0) / asphi;
  double Ftip = 2.0 / CC_PI * acos(exp(-factortip));
  double factorhub = B / 2.0 * (r / Rhub - 1.0) / asphi;
  double Fhub = 2.0 / CC_PI * acos(exp(-factorhub));
  return Ftip * Fhub;
}

/* ---------------- linear lift curve ---------------- */

void cc_linearliftcoeff(const double* alpha, const double* cl, int n,
                        double alphamin, double alphamax, double* m,
                        double* alpha0) {
  /* Julia: idxmin = first index with alpha > alphamin;
   *        idxmax = first index with alpha > alphamax;
   *        subrange alpha[idxmin..idxmax] inclusive (1-based) */
  int idxmin = n, idxmax = n;
  for (int i = 0; i < n; i++) {
    if (alpha[i] > alphamin) {
      idxmin = i;
      break;
    }
  }
  for (int i = 0; i < n; i++) {
    if (alpha[i] > alphamax) {
      idxmax = i;
      break;
    }
  }
  if (idxmax < idxmin) idxmax = idxmin;
  int cnt = idxmax - idxmin + 1;

  /* least squares cl = m*alpha + b */
  double sx = 0, sy = 0, sxx = 0, sxy = 0;
  for (int i = idxmin; i <= idxmax; i++) {
    sx += alpha[i];
    sy += cl[i];
    sxx += alpha[i] * alpha[i];
    sxy += alpha[i] * cl[i];
  }
  double denom = cnt * sxx - sx * sx;
  double slope = (cnt * sxy - sx * sy) / denom;
  double intercept = (sy - slope * sx) / cnt;
  *m = slope;
  *alpha0 = -intercept / slope;
}

/* ---------------- Viterna extrapolation ---------------- */

int cc_viterna(const double* alpha, const double* cl, const double* cd, int n,
               double cr75, int nalpha_out, double** alpha_out, double** cl_out,
               double** cd_out, int* n_out) {
  /* estimate cdmax */
  double AR = 1.0 / cr75;
  double cdmaxAR = 1.11 + 0.018 * AR;
  double cdmax = cd[0];
  for (int i = 1; i < n; i++)
    if (cd[i] > cdmax) cdmax = cd[i];
  if (cdmaxAR > cdmax) cdmax = cdmaxAR;

  /* positive stall: first index of max cl */
  int i_ps = 0;
  for (int i = 1; i < n; i++)
    if (cl[i] > cl[i_ps]) i_ps = i;
  double cl_ps = cl[i_ps], cd_ps = cd[i_ps], a_ps = alpha[i_ps];

  /* negative stall: argmin over cl where alpha < a_ps */
  int i_ns = -1;
  for (int i = 0; i < n; i++) {
    if (alpha[i] < a_ps) {
      if (i_ns < 0 || cl[i] < cl[i_ns]) i_ns = i;
    }
  }
  if (i_ns < 0) return -1; /* no data below stall: invalid input */
  double cl_ns = cl[i_ns], cd_ns = cd[i_ns], a_ns = alpha[i_ns];

  int na = nalpha_out;
  double* A1pos = (double*)malloc(sizeof(double) * na);
  double* B2pos = (double*)malloc(sizeof(double) * na);
  double* A2neg = (double*)malloc(sizeof(double) * na);
  double* B2neg = (double*)malloc(sizeof(double) * na);
  double* apos = (double*)malloc(sizeof(double) * na);
  double* aneg = (double*)malloc(sizeof(double) * na);
  double* adjpos = (double*)malloc(sizeof(double) * na);
  double* adjneg = (double*)malloc(sizeof(double) * na);
  if (!A1pos || !B2pos || !A2neg || !B2neg || !apos || !aneg || !adjpos ||
      !adjneg) {
    free(A1pos);
    free(B2pos);
    free(A2neg);
    free(B2neg);
    free(apos);
    free(aneg);
    free(adjpos);
    free(adjneg);
    return -1;
  }

  double B1pos = cdmax;
  double sa = sin(a_ps), ca = cos(a_ps);
  double A2pos = (cl_ps - cdmax * sa * ca) * sa / (ca * ca);
  double B2posv = (cd_ps - cdmax * sa * sa) / ca;
  for (int i = 0; i < na; i++) {
    A1pos[i] = B1pos / 2.0;
    B2pos[i] = B2posv;
  }

  double B1neg = cdmax;
  double A1neg = B1neg / 2.0;
  sa = sin(a_ns);
  ca = cos(a_ns);
  double A2negv = (cl_ns - cdmax * sa * ca) * sa / (ca * ca);
  double B2negv = (cd_ns - cdmax * sa * sa) / ca;
  for (int i = 0; i < na; i++) {
    A2neg[i] = A2negv;
    B2neg[i] = B2negv;
  }

  /* apos = range(alpha[end], pi, length=na+1)[2:]; aneg = range(-pi,
   * alpha[1], length=na+1)[1:end-1] */
  for (int i = 0; i < na; i++) {
    apos[i] = alpha[n - 1] + (CC_PI - alpha[n - 1]) * (i + 1) / na;
    aneg[i] = -CC_PI + (alpha[0] + CC_PI) * i / na;
    adjpos[i] = 1.0;
    adjneg[i] = 1.0;
  }

  for (int i = 0; i < na; i++) {
    if (apos[i] >= CC_PI / 2.0) {
      adjpos[i] = -0.7;
      A1pos[i] *= -1.0;
      B2pos[i] *= -1.0;
    }
  }
  for (int i = 0; i < na; i++) {
    if (aneg[i] <= -CC_PI / 2.0) {
      adjneg[i] = 0.7;
      A2neg[i] *= -1.0;
      B2neg[i] *= -1.0;
    }
  }

  double* clpos = (double*)malloc(sizeof(double) * na);
  double* cdpos = (double*)malloc(sizeof(double) * na);
  double* clneg = (double*)malloc(sizeof(double) * na);
  double* cdneg = (double*)malloc(sizeof(double) * na);
  if (!clpos || !cdpos || !clneg || !cdneg) {
    free(clpos);
    free(cdpos);
    free(clneg);
    free(cdneg);
    free(A1pos);
    free(B2pos);
    free(A2neg);
    free(B2neg);
    free(apos);
    free(aneg);
    free(adjpos);
    free(adjneg);
    return -1;
  }

  for (int i = 0; i < na; i++) {
    double s2 = sin(2.0 * apos[i]), c = cos(apos[i]), s = sin(apos[i]);
    clpos[i] = adjpos[i] * (A1pos[i] * s2 + A2pos * c * c / s);
    cdpos[i] = B1pos * s * s + B2pos[i] * c;
  }
  for (int i = 0; i < na; i++) {
    double s2 = sin(2.0 * aneg[i]), c = cos(aneg[i]), s = sin(aneg[i]);
    clneg[i] = adjneg[i] * (A1neg * s2 + A2neg[i] * c * c / s);
    cdneg[i] = B1neg * s * s + B2neg[i] * c;
  }

  /* linear variation at the ends */
  for (int i = 0; i < na; i++) {
    if (apos[i] >= CC_PI - a_ps)
      clpos[i] = (apos[i] - CC_PI) / a_ps * cl_ps * 0.7;
  }
  for (int i = 0; i < na; i++) {
    if (aneg[i] <= -CC_PI - a_ns)
      clneg[i] = (aneg[i] + CC_PI) / a_ns * cl_ns * 0.7;
  }

  /* concatenate: [aneg; alpha; apos] */
  int ntotal = 2 * na + n;
  double* af_ = (double*)malloc(sizeof(double) * ntotal);
  double* clf = (double*)malloc(sizeof(double) * ntotal);
  double* cdf = (double*)malloc(sizeof(double) * ntotal);
  if (!af_ || !clf || !cdf) {
    free(af_);
    free(clf);
    free(cdf);
    free(clpos);
    free(cdpos);
    free(clneg);
    free(cdneg);
    free(A1pos);
    free(B2pos);
    free(A2neg);
    free(B2neg);
    free(apos);
    free(aneg);
    free(adjpos);
    free(adjneg);
    return -1;
  }

  for (int i = 0; i < na; i++) {
    af_[i] = aneg[i];
    clf[i] = clneg[i];
    cdf[i] = cdneg[i];
  }
  for (int i = 0; i < n; i++) {
    af_[na + i] = alpha[i];
    clf[na + i] = cl[i];
    cdf[na + i] = cd[i];
  }
  for (int i = 0; i < na; i++) {
    af_[na + n + i] = apos[i];
    clf[na + n + i] = clpos[i];
    cdf[na + n + i] = cdpos[i];
  }
  for (int i = 0; i < ntotal; i++)
    if (cdf[i] < 0.0001) cdf[i] = 0.0001;

  *alpha_out = af_;
  *cl_out = clf;
  *cd_out = cdf;
  *n_out = ntotal;

  free(clpos);
  free(cdpos);
  free(clneg);
  free(cdneg);
  free(A1pos);
  free(B2pos);
  free(A2neg);
  free(B2neg);
  free(apos);
  free(aneg);
  free(adjpos);
  free(adjneg);
  return 0;
}
