/*
 * ccblade.c — port of CCBlade.jl's BEM core: residual, bracketing, root
 * solve (Brent method, ported from FLOWMath's scipy-brentq style
 * implementation), inflow convenience constructors, thrust/torque
 * integration, and nondimensionalization.
 */
#include "ccblade.h"

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static double cc_sign(double x) { return (x > 0) - (x < 0); }

static void copy2(double dst[2], const double src[2]) {
  dst[0] = src[0];
  dst[1] = src[1];
}

CCRotor cc_rotor(double Rhub, double Rtip, int B) {
  CCRotor r;
  r.Rhub = Rhub;
  r.Rtip = Rtip;
  r.B = B;
  r.precone = 0.0;
  r.turbine = 0;
  r.mach = CC_MACH_NONE;
  r.re = CC_RE_NONE;
  r.rotation = CC_ROT_NONE;
  r.tip = CC_TIP_PRANDTL_HUB;
  r.sf_re0 = 0.0;
  r.sf_p = 0.0;
  r.dse_a = 1.0;
  r.dse_b = 1.0;
  r.dse_d = 1.0;
  r.dse_m = 2.0 * CC_PI;
  r.dse_alpha0 = 0.0;
  return r;
}

CCOperatingPoint cc_simple_op_full(double Vinf, double Omega, double r,
                                   double rho, double pitch, double mu,
                                   double asound, double precone) {
  CCOperatingPoint op;
  op.Vx = Vinf * cos(precone);
  op.Vy = Omega * r * cos(precone);
  op.rho = rho;
  op.pitch = pitch;
  op.mu = mu;
  op.asound = asound;
  return op;
}

CCOperatingPoint cc_simple_op(double Vinf, double Omega, double r, double rho) {
  return cc_simple_op_full(Vinf, Omega, r, rho, 0.0, 1.0, 1.0, 0.0);
}

CCOperatingPoint cc_windturbine_op(double Vhub, double Omega, double pitch,
                                   double r, double precone, double yaw,
                                   double tilt, double azimuth, double hubHt,
                                   double shearExp, double rho, double mu,
                                   double asound) {
  double sy = sin(yaw), cy = cos(yaw);
  double st = sin(tilt), ct = cos(tilt);
  double sa = sin(azimuth), ca = cos(azimuth);
  double sc = sin(precone), cc = cos(precone);

  double x_az = -r * sc;
  double z_az = r * cc;
  double y_az = 0.0;

  double heightFromHub = (y_az * sa + z_az * ca) * ct - x_az * st;

  double V = Vhub * pow(1.0 + heightFromHub / hubHt, shearExp);

  double Vwind_x = V * ((cy * st * ca + sy * sa) * sc + cy * ct * cc);
  double Vwind_y = V * (cy * st * sa - sy * ca);

  double Vrot_x = -Omega * y_az * sc;
  double Vrot_y = Omega * z_az;

  CCOperatingPoint op;
  op.Vx = Vwind_x + Vrot_x;
  op.Vy = Vwind_y + Vrot_y;
  op.rho = rho;
  op.pitch = pitch;
  op.mu = mu;
  op.asound = asound;
  return op;
}

/* ---------------- residual ---------------- */

static CCOutputs outputs_zero(void) {
  CCOutputs o;
  o.Np = o.Tp = o.a = o.ap = o.u = o.v = o.phi = o.alpha = o.W = o.cl = o.cd =
      o.cn = o.ct = o.F = o.G = 0.0;
  return o;
}

static double residual_and_outputs(double phi, const CCRotor* rotor,
                                   const CCSection* section,
                                   const CCOperatingPoint* op,
                                   CCOutputs* outp) {
  double r = section->r, chord = section->chord, theta = section->theta;
  double Rhub = rotor->Rhub, Rtip = rotor->Rtip;
  double Vx = op->Vx, Vy = op->Vy, rho = op->rho, pitch = op->pitch;
  double mu = op->mu, asound = op->asound;
  const CCAirfoil* af = section->af;
  int B = rotor->B, turbine = rotor->turbine;

  CCOutputs out = outputs_zero();

  double sigma_p = B * chord / (2.0 * CC_PI * r);
  double sphi = sin(phi);
  double cphi = cos(phi);

  double alpha = (theta + pitch) - phi;

  /* Reynolds/Mach ignoring induction (as in Julia) */
  double W0 = sqrt(Vx * Vx + Vy * Vy);
  double Re = rho * W0 * chord / mu;
  double Mach = W0 / asound;

  double cl, cd;
  if (turbine) {
    cc_afeval(af, -alpha, Re, Mach, &cl, &cd);
    cl *= -1.0;
  } else {
    cc_afeval(af, alpha, Re, Mach, &cl, &cd);
  }

  if (rotor->re == CC_RE_SKIN_FRICTION)
    cc_re_correction_skin_friction(rotor->sf_re0, rotor->sf_p, cl, cd, Re, &cl,
                                   &cd);
  if (rotor->mach == CC_MACH_PRANDTL_GLAUERT)
    cc_mach_correction_prandtl_glauert(cl, cd, Mach, &cl, &cd);
  if (rotor->rotation == CC_ROT_DUSELIG_EGGERS)
    cc_rotation_correction_duselig_eggers(
        rotor->dse_a, rotor->dse_b, rotor->dse_d, rotor->dse_m,
        rotor->dse_alpha0, cl, cd, chord / r, r / Rtip, Vy / Vx * Rtip / r,
        alpha, phi, 30.0 * CC_PI / 180.0, &cl, &cd);

  double cn = cl * cphi - cd * sphi;
  double ct = cl * sphi + cd * cphi;

  double F = 1.0;
  if (rotor->tip == CC_TIP_PRANDTL)
    F = cc_tip_correction_prandtl(r, Rtip, phi, B);
  else if (rotor->tip == CC_TIP_PRANDTL_HUB)
    F = cc_tip_correction_prandtl_hub(r, Rhub, Rtip, phi, B);

  double k = cn * sigma_p / (4.0 * F * sphi * sphi);
  double kp = ct * sigma_p / (4.0 * F * sphi * cphi);

  double a, ap, u, v, R;

  if (fabs(Vx) <= 1e-6) {
    u = cc_sign(phi) * kp * cn / ct * Vy;
    v = 0.0;
    a = 0.0;
    ap = 0.0;
    R = cc_sign(phi) - k;
  } else if (fabs(Vy) <= 1e-6) {
    u = 0.0;
    v = k * ct / cn * fabs(Vx);
    a = 0.0;
    ap = 0.0;
    R = cc_sign(Vx) + kp;
  } else {
    if (phi < 0.0) k *= -1.0;

    if (fabs(k - 1.0) <= 1e-6) {
      if (outp) *outp = outputs_zero();
      return 1.0;
    }

    if (k >= -2.0 / 3.0) {
      a = k / (1.0 - k);
    } else {
      /* empirical region (Buhl(F=1) * F, applied via loads) */
      double g1 = 2.0 * k + 1.0 / 9.0;
      double g2 = -2.0 * k - 1.0 / 3.0;
      double g3 = -2.0 * k - 7.0 / 9.0;
      a = (g1 + sqrt(g2)) / g3;
    }

    u = a * Vx;

    if (Vx < 0.0) kp *= -1.0;

    if (fabs(kp + 1.0) <= 1e-6) {
      if (outp) *outp = outputs_zero();
      return 1.0;
    }

    ap = kp / (1.0 + kp);
    v = ap * Vy;

    R = sphi / (1.0 + a) - (Vx / Vy) * cphi / (1.0 - ap);
  }

  double W = sqrt((Vx + u) * (Vx + u) + (Vy - v) * (Vy - v));
  double Np = cn * 0.5 * rho * W * W * chord;
  double Tp = ct * 0.5 * rho * W * W * chord;

  double G;
  if (fabs(Vx) <= 1e-6)
    G = sqrt(F);
  else if (fabs(Vy) <= 1e-6)
    G = F;
  else
    G = (-1.0 + sqrt(1.0 + 4.0 * a * (1.0 + a) * F)) / (2.0 * a);
  u *= G;
  v *= G;

  if (turbine) {
    out.Np = -Np;
    out.Tp = -Tp;
    out.a = -a;
    out.ap = -ap;
    out.u = -u;
    out.v = -v;
    out.phi = phi;
    out.alpha = -alpha;
    out.W = W;
    out.cl = -cl;
    out.cd = cd;
    out.cn = -cn;
    out.ct = -ct;
    out.F = F;
    out.G = G;
  } else {
    out.Np = Np;
    out.Tp = Tp;
    out.a = a;
    out.ap = ap;
    out.u = u;
    out.v = v;
    out.phi = phi;
    out.alpha = alpha;
    out.W = W;
    out.cl = cl;
    out.cd = cd;
    out.cn = cn;
    out.ct = ct;
    out.F = F;
    out.G = G;
  }

  if (outp) *outp = out;
  return R;
}

/* ---------------- bracketing + root solve ---------------- */

typedef struct {
  const CCRotor* rotor;
  const CCSection* section;
  const CCOperatingPoint* op;
} CCResidualCtx;

// cppcheck-suppress constParameterCallback
static double residual_fn(double phi, void* params) {
  const CCResidualCtx* ctx = (const CCResidualCtx*)params;
  return residual_and_outputs(phi, ctx->rotor, ctx->section, ctx->op, NULL);
}

/* Port of FLOWMath.brent (scipy brentq) with its default tolerances
 * atol=2e-12, rtol=4*eps, maxiter=100. Requires f(a)*f(b) <= 0.
 * Returns the root; on a sign error at the bracket returns 0.0 (as FLOWMath
 * does), and on non-convergence returns the last iterate. */
static double cc_brent(double (*f)(double, void*), void* params, double a,
                       double b) {
  const double atol = 2e-12;
  const double rtol = 4.0 * DBL_EPSILON;
  const int maxiter = 100;

  double xpre = a, xcur = b;
  double fpre = f(xpre, params);
  double fcur = f(xcur, params);
  double xblk = 0.0, fblk = 0.0, spre = 0.0, scur = 0.0;

  if (fpre * fcur > 0.0) return 0.0;
  if (fpre == 0.0) return xpre;
  if (fcur == 0.0) return xcur;

  for (int i = 0; i < maxiter; i++) {
    if (fpre * fcur < 0.0) {
      xblk = xpre;
      fblk = fpre;
      spre = scur = xcur - xpre;
    }
    if (fabs(fblk) < fabs(fcur)) {
      /* sequential assignments, matching scipy/FLOWMath */
      xpre = xcur;
      xcur = xblk;
      xblk = xpre;

      fpre = fcur;
      fcur = fblk;
      fblk = fpre;
    }

    double delta = (atol + rtol * fabs(xcur)) / 2.0;
    double sbis = (xblk - xcur) / 2.0;
    if (fcur == 0.0 || fabs(sbis) < delta) return xcur;

    if (fabs(spre) > delta && fabs(fcur) < fabs(fpre)) {
      double stry;
      if (xpre == xblk) {
        /* interpolate */
        stry = -fcur * (xcur - xpre) / (fcur - fpre);
      } else {
        /* extrapolate */
        double dpre = (fpre - fcur) / (xpre - xcur);
        double dblk = (fblk - fcur) / (xblk - xcur);
        stry =
            -fcur * (fblk * dblk - fpre * dpre) / (dblk * dpre * (fblk - fpre));
      }
      double bound = fabs(spre);
      double alt = 3.0 * fabs(sbis) - delta;
      if (alt < bound) bound = alt;
      if (2.0 * fabs(stry) < bound) {
        spre = scur;
        scur = stry;
      } else {
        spre = sbis;
        scur = sbis;
      }
    } else {
      spre = sbis;
      scur = sbis;
    }

    xpre = xcur;
    fpre = fcur;
    if (fabs(scur) > delta)
      xcur += scur;
    else
      xcur += (sbis > 0.0 ? delta : -delta);

    fcur = f(xcur, params);
  }
  return xcur;
}

/* port of firstbracket: subdivide (xmin, xmax) into n intervals; find the
 * first sign change (searching from xmin, or from xmax if backwardsearch).
 * Returns 1 and fills xl/xu (ordered low, high) if found. */
static int firstbracket(double (*f)(double, void*), void* params, double xmin,
                        double xmax, int n, int backwardsearch, double* xl,
                        double* xu) {
  double step = (xmax - xmin) / n;
  if (!backwardsearch) {
    double fprev = f(xmin, params);
    for (int i = 1; i <= n; i++) {
      double x = xmin + step * i;
      double fnext = f(x, params);
      if (fprev * fnext < 0.0) {
        *xl = xmin + step * (i - 1);
        *xu = x;
        return 1;
      }
      fprev = fnext;
    }
  } else {
    double fprev = f(xmax, params);
    for (int i = 1; i <= n; i++) {
      double x = xmax - step * i;
      double fnext = f(x, params);
      if (fprev * fnext < 0.0) {
        *xl = x;
        *xu = xmax - step * (i - 1);
        return 1;
      }
      fprev = fnext;
    }
  }
  return 0;
}

CCOutputs cc_solve_ex(const CCRotor* rotor, const CCSection* section,
                      const CCOperatingPoint* op, int npts,
                      int forcebackwardsearch, int epsilon_everywhere,
                      int* success) {
  if (success) *success = 1;

  /* no loads at hub/tip */
  if (fabs(section->r - rotor->Rhub) <= 1e-6 ||
      fabs(section->r - rotor->Rtip) <= 1e-6)
    return outputs_zero();

  double Vx = op->Vx, Vy = op->Vy;
  double theta = section->theta + op->pitch;

  int Vx_is_zero = fabs(Vx) <= 1e-6;
  int Vy_is_zero = fabs(Vy) <= 1e-6;

  const double epsilon = 1e-6;
  double q1[2], q2[2], q3[2], q4[2];
  if (epsilon_everywhere) {
    q1[0] = epsilon;
    q1[1] = CC_PI / 2 - epsilon;
    q2[0] = -CC_PI / 2 + epsilon;
    q2[1] = -epsilon;
    q3[0] = CC_PI / 2 + epsilon;
    q3[1] = CC_PI - epsilon;
    q4[0] = -CC_PI + epsilon;
    q4[1] = -CC_PI / 2 - epsilon;
  } else {
    q1[0] = epsilon;
    q1[1] = CC_PI / 2;
    q2[0] = -CC_PI / 2;
    q2[1] = -epsilon;
    q3[0] = CC_PI / 2;
    q3[1] = CC_PI - epsilon;
    q4[0] = -CC_PI + epsilon;
    q4[1] = -CC_PI / 2;
  }

  double order[4][2];
  int norder = 0;
  int startfrom90 = 0;

  if (Vx_is_zero && Vy_is_zero) {
    return outputs_zero();
  } else if (Vx_is_zero) {
    startfrom90 = 0;
    if (Vy > 0 && theta > 0) {
      copy2(order[norder++], q1);
      copy2(order[norder++], q2);
    } else if (Vy > 0 && theta < 0) {
      copy2(order[norder++], q2);
      copy2(order[norder++], q1);
    } else if (Vy < 0 && theta > 0) {
      copy2(order[norder++], q3);
      copy2(order[norder++], q4);
    } else {
      copy2(order[norder++], q4);
      copy2(order[norder++], q3);
    }
  } else if (Vy_is_zero) {
    startfrom90 = 1;
    if (Vx > 0 && fabs(theta) < CC_PI / 2) {
      copy2(order[norder++], q1);
      copy2(order[norder++], q3);
    } else if (Vx < 0 && fabs(theta) < CC_PI / 2) {
      copy2(order[norder++], q2);
      copy2(order[norder++], q4);
    } else if (Vx > 0 && fabs(theta) > CC_PI / 2) {
      copy2(order[norder++], q3);
      copy2(order[norder++], q1);
    } else {
      copy2(order[norder++], q4);
      copy2(order[norder++], q2);
    }
  } else {
    startfrom90 = 0;
    if (Vx > 0 && Vy > 0) {
      copy2(order[norder++], q1);
      copy2(order[norder++], q2);
      copy2(order[norder++], q3);
      copy2(order[norder++], q4);
    } else if (Vx < 0 && Vy > 0) {
      copy2(order[norder++], q2);
      copy2(order[norder++], q1);
      copy2(order[norder++], q4);
      copy2(order[norder++], q3);
    } else if (Vx > 0 && Vy < 0) {
      copy2(order[norder++], q3);
      copy2(order[norder++], q4);
      copy2(order[norder++], q1);
      copy2(order[norder++], q2);
    } else {
      copy2(order[norder++], q4);
      copy2(order[norder++], q3);
      copy2(order[norder++], q2);
      copy2(order[norder++], q1);
    }
  }

  CCResidualCtx ctx = {rotor, section, op};

  for (int j = 0; j < norder; j++) {
    double phimin = order[j][0], phimax = order[j][1];

    int backwardsearch;
    if (forcebackwardsearch) {
      backwardsearch = 1;
    } else {
      backwardsearch = 0;
      if (!startfrom90) {
        if (phimin == -CC_PI / 2 || phimax == -CC_PI / 2) backwardsearch = 1;
      } else {
        if (phimax == CC_PI / 2) backwardsearch = 1;
      }
    }

    double phiL, phiU;
    int found = firstbracket(residual_fn, &ctx, phimin, phimax, npts,
                             backwardsearch, &phiL, &phiU);
    if (!found) continue;

    double phistar = cc_brent(residual_fn, &ctx, phiL, phiU);
    CCOutputs out;
    residual_and_outputs(phistar, rotor, section, op, &out);
    return out;
  }

  fprintf(stderr,
          "ccblade: no bracket found for section (r=%g); "
          "zero loading assumed\n",
          section->r);
  if (success) *success = 0;
  return outputs_zero();
}

CCOutputs cc_solve(const CCRotor* rotor, const CCSection* section,
                   const CCOperatingPoint* op) {
  return cc_solve_ex(rotor, section, op, 10, 0, 0, NULL);
}

/* ---------------- integration ---------------- */

void cc_thrusttorque(const CCRotor* rotor, const CCSection* sections,
                     const CCOutputs* outputs, int n, double* T, double* Q) {
  double* rfull = (double*)malloc(sizeof(double) * (n + 2));
  double* Npfull = (double*)malloc(sizeof(double) * (n + 2));
  double* torque = (double*)malloc(sizeof(double) * (n + 2));
  if (!rfull || !Npfull || !torque) {
    free(rfull);
    free(Npfull);
    free(torque);
    *T = *Q = NAN;
    return;
  }

  double cpc = cos(rotor->precone);
  rfull[0] = rotor->Rhub;
  Npfull[0] = 0.0;
  torque[0] = 0.0;
  for (int i = 0; i < n; i++) {
    rfull[i + 1] = sections[i].r;
    Npfull[i + 1] = outputs[i].Np * cpc;
    torque[i + 1] = outputs[i].Tp * rfull[i + 1] * cpc;
  }
  rfull[n + 1] = rotor->Rtip;
  Npfull[n + 1] = 0.0;
  torque[n + 1] = 0.0;

  *T = rotor->B * cc_trapz(rfull, Npfull, n + 2);
  *Q = rotor->B * cc_trapz(rfull, torque, n + 2);

  free(rfull);
  free(Npfull);
  free(torque);
}

void cc_thrusttorque_az(const CCRotor* rotor, const CCSection* sections,
                        const CCOutputs* outputs, int nr, int naz, double* T,
                        double* Q) {
  /* outputs[i*naz + j] is section i, azimuth j: gather each azimuth column
   * into a contiguous array before integrating. */
  CCOutputs* col = (CCOutputs*)malloc(sizeof(CCOutputs) * nr);
  if (!col) {
    *T = *Q = NAN;
    return;
  }
  *T = 0.0;
  *Q = 0.0;
  for (int j = 0; j < naz; j++) {
    for (int i = 0; i < nr; i++) col[i] = outputs[i * naz + j];
    double Tj, Qj;
    cc_thrusttorque(rotor, sections, col, nr, &Tj, &Qj);
    *T += Tj / naz;
    *Q += Qj / naz;
  }
  free(col);
}

/* ---------------- nondim ---------------- */

void cc_nondim(double T, double Q, double Vhub, double Omega, double rho,
               const CCRotor* rotor, const char* rotortype, double out[3]) {
  double P = Q * Omega;
  double Rp = rotor->Rtip * cos(rotor->precone);

  if (strcmp(rotortype, "windturbine") == 0) {
    double q = 0.5 * rho * Vhub * Vhub;
    double A = CC_PI * Rp * Rp;
    out[0] = P / (q * A * Vhub);
    out[1] = T / (q * A);
    out[2] = Q / (q * Rp * A);
  } else if (strcmp(rotortype, "propeller") == 0) {
    double n = Omega / (2.0 * CC_PI);
    double Dp = 2.0 * Rp;
    if (T < 0.0)
      out[0] = 0.0;
    else
      out[0] = T * Vhub / P;
    out[1] = T / (rho * n * n * Dp * Dp * Dp * Dp);
    out[2] = Q / (rho * n * n * Dp * Dp * Dp * Dp * Dp);
  } else { /* helicopter */
    double A = CC_PI * Rp * Rp;
    double vtip2 = (Omega * Rp) * (Omega * Rp);
    double CT = T / (rho * A * vtip2);
    double CP = P / (rho * A * vtip2 * (Omega * Rp));
    if (CT < 0.0)
      out[0] = 0.0;
    else
      out[0] = pow(CT, 1.5) / (sqrt(2.0) * CP);
    out[1] = CT;
    out[2] = CP;
  }
}
