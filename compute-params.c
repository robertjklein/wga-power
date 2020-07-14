/*
 * compute-params.c
 * 
 * Computes a-g parameters given input on command line
 * Wrapped by generate-command-lines.pl
 *
 * Robert J. Klein
 * April 10. 2007
 */

#include <gsl/gsl_poly.h>

#include <stdio.h>
#include <stdlib.h>

int main (int argc, char **argv) {
  double p1, p2, g1, g2, K;
  double Qcub, Acub, Bcub, Ccub;
  int res, i, good;
  double x[3];
  double a,b,c,d,e,f;

  if (argc != 6) {
    printf ("0\n0\n0\n0\n0\n0\n");
    exit(5);
  }

  p1 = atof(argv[1]);
  p2 = atof(argv[2]);
  g1 = atof(argv[3]);
  g2 = atof(argv[4]);
  K = atof(argv[5]);

  if (K > 0) {          /* Selected controls */
    Qcub = -1. + g1 + g2 - g1*g2;
    Acub = (g1 - 1.) * K + p1 - 2*g1*p1 -g1*p2+g2*(-1.+g1+K-g1*K-p1+2.*g1*p1+p2);
    Bcub = p1*(g2*(1.-K+g1*(-2.+2.*K-p1)-p2)+g1*(-1.*K+p1+p2));
    Ccub = g1*g2*p1*p1*(1.-K);
    res = gsl_poly_solve_cubic (Acub/Qcub, Bcub/Qcub, Ccub/Qcub,
				x, x+1, x+2);
    good = 0;
    for (i=0; i<res; i++) {
      b = x[i];
      
      d = x[i] * p2 / (x[i]-x[i]*g1+g1*p1);
      f = (-1.*x[i]+x[i]*p1+x[i]*p2)/(-1.*x[i]+x[i]*g2-g2*p1);
      a = p1-x[i];
      c = p2-(x[i]*p2)/(x[i]-x[i]*g1+g1*p1);
      e = 1.-p1-p2-(-1.*x[i]+x[i]*p1+x[i]*p2)/(-1.*x[i]+x[i]*g2-g2*p1);
      if (a >= 0 && a <= 1 && b >= 0 && b <= 1 && c >= 0 && c <= 1 && d >= 0 && d <= 1 && e >= 0 && e <= 1 && f >= 0 && f <= 1) {
	good = 1;
	break;
      }
    }
    if (good == 0) {
      printf ("0\n0\n0\n0\n0\n0\n");
      exit(20);
    }
    a /= K;
    c /= K;
    e /= K;
    b /= (1.-K);
    d /= (1.-K);
    f /= (1.-K);
  } else {
    b = p1;
    d = p2;
    f = 1.-p1-p2;
    a = p1/(-1.*g2-p1+g2*p1-g1*p2+g2*p2);
    c = g1*p2/(g2+p1-g2*p1+g1*p2-g2*p2);
    e = (g2-g2*p1-g2*p2)/(-1.*g2-p1+g2*p1-g1*p2+g2*p2);
  }
  printf ("%f\n%f\n%f\n%f\n%f\n%f\n", a,b,c,d,e,f);
}
