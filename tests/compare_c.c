/* C side of the sweep comparison: same geometry/conditions as
 * compare_julia.jl, writes sweep_c.csv. */
#include <stdio.h>
#include <stdlib.h>

#include "ccblade.h"

int main(void) {
  double Rhub = 0.0254 * 0.5;
  double Rtip = 0.0254 * 3.0;
  CCRotor rotor = cc_rotor(Rhub, Rtip, 2);

  static const double r_in[31] = {
      0.7526, 0.7928, 0.8329, 0.8731, 0.9132, 0.9586, 1.0332, 1.1128,
      1.1925, 1.2722, 1.3519, 1.4316, 1.5114, 1.5911, 1.6708, 1.7505,
      1.8302, 1.9099, 1.9896, 2.0693, 2.1490, 2.2287, 2.3084, 2.3881,
      2.4678, 2.5475, 2.6273, 2.7070, 2.7867, 2.8661, 2.9410};
  static const double chord_in[31] = {
      0.6270, 0.6255, 0.6231, 0.6199, 0.6165, 0.6125, 0.6054, 0.5973,
      0.5887, 0.5794, 0.5695, 0.5590, 0.5479, 0.5362, 0.5240, 0.5111,
      0.4977, 0.4836, 0.4689, 0.4537, 0.4379, 0.4214, 0.4044, 0.3867,
      0.3685, 0.3497, 0.3303, 0.3103, 0.2897, 0.2618, 0.1920};
  static const double theta_deg[31] = {
      40.2273, 38.7657, 37.3913, 36.0981, 34.8803, 33.5899, 31.6400, 29.7730,
      28.0952, 26.5833, 25.2155, 23.9736, 22.8421, 21.8075, 20.8586, 19.9855,
      19.1800, 18.4347, 17.7434, 17.1005, 16.5013, 15.9417, 15.4179, 14.9266,
      14.4650, 14.0306, 13.6210, 13.2343, 12.8685, 12.5233, 12.2138};

  CCAirfoil* af = cc_af_alpha_file("airfoils/NACA64_A17.dat", 0);
  if (!af) {
    fprintf(stderr, "cannot load NACA64_A17.dat\n");
    return 1;
  }

  CCSection sections[31];
  for (int i = 0; i < 31; i++) {
    sections[i].r = r_in[i] * 0.0254;
    sections[i].chord = chord_in[i] * 0.0254;
    sections[i].theta = theta_deg[i] * CC_PI / 180.0;
    sections[i].af = af;
  }

  static const double vvec[7] = {0.0, 2.5, 7.3, 13.7, 25.4, 40.9, 55.0};
  static const double rpmvec[6] = {0.0, 500.0, 1234.0, 4000.0, 7777.0, 12000.0};
  double rho = 1.225;

  FILE* f = fopen("sweep_c.csv", "w");
  if (!f) {
    fprintf(stderr, "cannot open sweep_c.csv\n");
    cc_af_free(af);
    return 1;
  }
  fprintf(f, "Vinf,RPM,T,Q\n");
  for (int iv = 0; iv < 7; iv++) {
    for (int ir = 0; ir < 6; ir++) {
      double Omega = rpmvec[ir] * CC_PI / 30.0;
      CCOutputs outs[31];
      for (int j = 0; j < 31; j++) {
        CCOperatingPoint op = cc_simple_op(vvec[iv], Omega, sections[j].r, rho);
        outs[j] = cc_solve(&rotor, &sections[j], &op);
      }
      double T, Q;
      cc_thrusttorque(&rotor, sections, outs, 31, &T, &Q);
      fprintf(f, "%.17g,%.17g,%.17g,%.17g\n", vvec[iv], rpmvec[ir], T, Q);
    }
  }
  fclose(f);
  cc_af_free(af);
  printf("sweep_c.csv written\n");
  return 0;
}
