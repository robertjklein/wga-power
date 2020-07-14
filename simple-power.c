/*
 * simple-power.c 
 * 
 * Given a number of cases, a number of controls, a GRR, a model (0=mult,
 * 1=add, 2=dom, 3=rec), a heterozygote_f, a homozygote_f (minor), a
 * qchisq value, and r^2 with marker, reports the power
 *
 * Robert J. Klein
 * September 29, 2006
 */

#include "powerfuncs.h"

int main (int argc, char **argv) {
  float n_cases, n_controls, grr, hetero_f, homo_f;
  int model;
  float q_chisq, r2;
  float prevalence;

  if (argc != 10) {
    fprintf (stderr, "Wrong usage\n");
    exit(1);
  }
  n_cases = atof(argv[1]);
  n_controls = atof(argv[2]);
  grr = atof(argv[3]);
  model = atoi(argv[4]);
  prevalence = atof(argv[5]);
  hetero_f = atof(argv[6]);
  homo_f = atof(argv[7]);
  q_chisq = atof(argv[8]);
  r2 = atof(argv[9]);

  printf ("Power is %f\n", individual_power_calc(q_chisq, 2, r2, lambda_calc(2, n_cases, n_controls, grr, hetero_f, homo_f, model, prevalence)));
}
