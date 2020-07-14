/*
 *
 * Given a value for lambda and num_markers on command line, reports power
 *
 * Robert J. Klein
 * March 13, 2006
 */

#include <stdio.h>
#include <stdlib.h>

#define MATHLIB_STANDALONE
#include <Rmath.h>

#include "powerfuncs.h"

int main (int argc, char **argv) {
  float lambda, alpha;
  double q_chisq;

  if (argc != 3) {
    fprintf (stderr, "USAGE: power-from-lambda-alpha <lambda> <alpha>\n");
    exit(2);
  }

  lambda = atof(argv[1]);
  alpha = 0.05/atoi(argv[2]);

  q_chisq = qchisq(1.0-alpha,1.,TRUE,0);
  printf ("Power is %f\n", individual_power_calc(q_chisq, 1, 1.0, lambda));
}


