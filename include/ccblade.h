/*
 * ccblade.h — C port of CCBlade.jl (Andrew Ning's BEM solver for
 * propellers/fans and wind turbines).
 *
 * Only the numerical core is ported. Julia-specific features (automatic
 * differentiation via ImplicitAD, StructArrays broadcasting) are not
 * applicable in C.
 *
 * Interpolation uses a faithful port of FLOWMath's Akima spline (including
 * its endpoint slope estimation and extrapolation-clamping behavior) and
 * the root solve is a port of FLOWMath's Brent (scipy brentq), so results
 * match the Julia implementation to within solver tolerances.
 * No external dependencies beyond libm.
 */
#ifndef CCBLADE_H
#define CCBLADE_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

#define CC_PI 3.14159265358979323846

/* ---------------- interpolation ---------------- */

typedef struct {
  int n;
  double* xdata;             /* owned copy, length n */
  double *p0, *p1, *p2, *p3; /* owned, length n-1 */
} CCAkima;

/* Build an Akima spline (FLOWMath-compatible). Requires n >= 2.
 * Returns NULL on allocation failure or n < 2. */
CCAkima* cc_akima(const double* x, const double* y, int n);
/* Evaluate; extrapolates using the end segments (like FLOWMath). */
double cc_akima_eval(const CCAkima* s, double x);
double cc_akima_deriv(const CCAkima* s, double x);
void cc_akima_free(CCAkima* s);

/* Recursive 1D-Akima interpolation over 2D/3D grids, matching
 * FLOWMath.interp2d/interp3d evaluation order.
 * f2 is indexed f2[i + nx*j] (i over x, j over y).
 * f3 is indexed f3[i + nx*(j + ny*k)]. */
double cc_interp2d_akima(const double* x, const double* y, const double* f2,
                         int nx, int ny, double xq, double yq);
double cc_interp3d_akima(const double* x, const double* y, const double* z,
                         const double* f3, int nx, int ny, int nz, double xq,
                         double yq, double zq);

double cc_trapz(const double* x, const double* y, int n);

/* ---------------- airfoil data ---------------- */

typedef enum {
  CC_AF_FUNCTION = 0, /* user callback */
  CC_AF_SIMPLE,       /* parameterized cl/cd curve */
  CC_AF_ALPHA,        /* table vs alpha (cached Akima) */
  CC_AF_ALPHA_RE,     /* table vs (alpha, Re) */
  CC_AF_ALPHA_MACH,   /* table vs (alpha, Mach) */
  CC_AF_ALPHA_RE_MACH /* table vs (alpha, Re, Mach) */
} CCAFKind;

typedef void (*CCAFFunc)(double alpha, double Re, double Mach, void* userdata,
                         double* cl, double* cd);

typedef struct {
  CCAFKind kind;

  /* CC_AF_FUNCTION */
  CCAFFunc func;
  void* userdata;

  /* CC_AF_SIMPLE */
  double m, alpha0, clmax, clmin, cd0, cd2;

  /* table data (CC_AF_ALPHA*); cl/cd layout:
   *   ALPHA:         cl[i]
   *   ALPHA_RE:      cl[i + nalpha*j]          j over Re
   *   ALPHA_MACH:    cl[i + nalpha*j]          j over Mach
   *   ALPHA_RE_MACH: cl[i + nalpha*(j + nre*k)] j over Re, k over Mach */
  int nalpha, nre, nmach;
  double *alpha, *re, *mach; /* owned */
  double *cl, *cd;           /* owned */
  char* info;                /* owned */
  double re_info, mach_info; /* informational only */

  /* cached splines (CC_AF_ALPHA only) */
  CCAkima *clspline, *cdspline;
} CCAirfoil;

CCAirfoil* cc_af_function(CCAFFunc f, void* userdata);
CCAirfoil* cc_af_simple(double m, double alpha0, double clmax, double clmin,
                        double cd0, double cd2);
CCAirfoil* cc_af_alpha(const double* alpha, const double* cl, const double* cd,
                       int n, const char* info, double re, double mach);
CCAirfoil* cc_af_alpha_file(const char* filename, int radians);
CCAirfoil* cc_af_alpha_re(const double* alpha, const double* re,
                          const double* cl, const double* cd, int nalpha,
                          int nre, const char* info, double mach);
CCAirfoil* cc_af_alpha_re_files(const char** filenames, int nfiles,
                                int radians);
CCAirfoil* cc_af_alpha_mach(const double* alpha, const double* mach,
                            const double* cl, const double* cd, int nalpha,
                            int nmach, const char* info, double re);
CCAirfoil* cc_af_alpha_mach_files(const char** filenames, int nfiles,
                                  int radians);
CCAirfoil* cc_af_alpha_re_mach(const double* alpha, const double* re,
                               const double* mach, const double* cl,
                               const double* cd, int nalpha, int nre, int nmach,
                               const char* info);
CCAirfoil* cc_af_alpha_re_mach_files(const char** filenames, int nre, int nmach,
                                     int radians);

void cc_af_free(CCAirfoil* af);

void cc_afeval(const CCAirfoil* af, double alpha, double Re, double Mach,
               double* cl, double* cd);

/* Write airfoil file(s). For CC_AF_ALPHA pass one filename; for the
 * higher-dimensional types pass nre/nmach (or nre*nmach) filenames in the
 * same order as the constructors. Returns 0 on success. */
int cc_write_af(const char* filename, const CCAirfoil* af, int radians);
int cc_write_af_multi(const char** filenames, const CCAirfoil* af, int radians);

/* Viterna post-stall extrapolation to +/-180 deg (port of CCBlade's viterna,
 * which follows the Viterna paper with NREL AirfoilPrep modifications).
 * Input arrays length n; output arrays (malloc'd, caller frees) have length
 * 2*nalpha_out + n. cr75 = chord/Rtip at 75% radius. Returns 0 on success. */
int cc_viterna(const double* alpha, const double* cl, const double* cd, int n,
               double cr75, int nalpha_out, double** alpha_out, double** cl_out,
               double** cd_out, int* n_out);

/* Linear regression on the linear portion of a lift curve (points with
 * alphamin < alpha[i] <= alphamax_next like CCBlade.linearliftcoeff:
 * subrange is alpha[idxmin..idxmax] with idxmin the first index with
 * alpha > alphamin and idxmax the first index with alpha > alphamax).
 * Returns slope m and zero-lift angle alpha0. */
void cc_linearliftcoeff(const double* alpha, const double* cl, int n,
                        double alphamin, double alphamax, double* m,
                        double* alpha0);

/* ---------------- corrections ---------------- */

typedef enum { CC_MACH_NONE = 0, CC_MACH_PRANDTL_GLAUERT } CCMachCorr;
typedef enum { CC_RE_NONE = 0, CC_RE_SKIN_FRICTION } CCReCorr;
typedef enum { CC_ROT_NONE = 0, CC_ROT_DUSELIG_EGGERS } CCRotCorr;
typedef enum { CC_TIP_NONE = 0, CC_TIP_PRANDTL, CC_TIP_PRANDTL_HUB } CCTipCorr;

void cc_mach_correction_prandtl_glauert(double cl, double cd, double Mach,
                                        double* cl_out, double* cd_out);
void cc_re_correction_skin_friction(double re0, double p, double cl, double cd,
                                    double Re, double* cl_out, double* cd_out);
/* Du-Selig (lift) + Eggers (drag) rotational stall-delay correction. */
void cc_rotation_correction_duselig_eggers(double a, double b, double d,
                                           double m, double alpha0, double cl,
                                           double cd, double cr, double rR,
                                           double tsr, double alpha, double phi,
                                           double alpha_max_corr,
                                           double* cl_out, double* cd_out);
double cc_tip_correction_prandtl(double r, double Rtip, double phi, int B);
double cc_tip_correction_prandtl_hub(double r, double Rhub, double Rtip,
                                     double phi, int B);

/* ---------------- rotor model ---------------- */

typedef struct {
  double Rhub, Rtip;
  int B;
  double precone;
  int turbine; /* nonzero -> wind turbine sign conventions */

  CCMachCorr mach;
  CCReCorr re;
  CCRotCorr rotation;
  CCTipCorr tip;

  double sf_re0, sf_p;                           /* CC_RE_SKIN_FRICTION */
  double dse_a, dse_b, dse_d, dse_m, dse_alpha0; /* CC_ROT_DUSELIG_EGGERS */
} CCRotor;

/* Defaults: precone=0, turbine=0, no Mach/Re/rotation corrections,
 * PrandtlTipHub tip loss (same as Julia Rotor(Rhub, Rtip, B)). */
CCRotor cc_rotor(double Rhub, double Rtip, int B);

typedef struct {
  double r, chord, theta;
  const CCAirfoil* af; /* not owned */
} CCSection;

typedef struct {
  double Vx, Vy; /* axial and tangential velocity along blade */
  double rho;
  double pitch;
  double mu, asound;
} CCOperatingPoint;

typedef struct {
  double Np, Tp;     /* normal/tangential load per unit length */
  double a, ap;      /* axial/tangential induction factors */
  double u, v;       /* axial/tangential induced velocity (with G correction) */
  double phi, alpha; /* inflow angle, angle of attack */
  double W;          /* relative velocity */
  double cl, cd;     /* lift/drag coefficients (as corrected) */
  double cn, ct;     /* normal/tangential force coefficients */
  double F, G;       /* hub/tip loss factor, effective factor for velocities */
} CCOutputs;

/* Convenience constructors matching the Julia versions. */
CCOperatingPoint cc_simple_op_full(double Vinf, double Omega, double r,
                                   double rho, double pitch, double mu,
                                   double asound, double precone);
CCOperatingPoint cc_simple_op(double Vinf, double Omega, double r, double rho);
CCOperatingPoint cc_windturbine_op(double Vhub, double Omega, double pitch,
                                   double r, double precone, double yaw,
                                   double tilt, double azimuth, double hubHt,
                                   double shearExp, double rho, double mu,
                                   double asound);

/* Solve the BEM residual for one section. Returns outputs; on failure
 * (no bracket found) returns all-zero outputs and, if success is non-NULL,
 * sets *success=0 (otherwise *success=1). */
CCOutputs cc_solve_ex(const CCRotor* rotor, const CCSection* section,
                      const CCOperatingPoint* op, int npts,
                      int forcebackwardsearch, int epsilon_everywhere,
                      int* success);
CCOutputs cc_solve(const CCRotor* rotor, const CCSection* section,
                   const CCOperatingPoint* op);

/* Integrate thrust/torque along the blade (trapezoidal, zero loads at
 * hub/tip) for n sections. */
void cc_thrusttorque(const CCRotor* rotor, const CCSection* sections,
                     const CCOutputs* outputs, int n, double* T, double* Q);
/* Azimuth-averaged version: outputs[i*naz + j] is section i, azimuth j. */
void cc_thrusttorque_az(const CCRotor* rotor, const CCSection* sections,
                        const CCOutputs* outputs, int nr, int naz, double* T,
                        double* Q);

/* Nondimensional coefficients. rotortype is one of "windturbine",
 * "propeller", "helicopter".
 *   windturbine: out = {CP, CT, CQ}
 *   propeller:   out = {efficiency, CT, CQ}
 *   helicopter:  out = {FM, CT, CP}  */
void cc_nondim(double T, double Q, double Vhub, double Omega, double rho,
               const CCRotor* rotor, const char* rotortype, double out[3]);

#ifdef __cplusplus
}
#endif

#endif /* CCBLADE_H */
