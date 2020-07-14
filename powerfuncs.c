/*
 * powerfuncs.c
 *
 * Functions for computing the power, regardless of whether we're finding best
 * SNPs or computing power for a give marker set.
 *
 * Robert J. Klein
 * May 12, 2005
 */

#define MATHLIB_STANDALONE
#include <Rmath.h>

#include <gsl/gsl_poly.h>

#include "structs.h"
#include "powerfuncs.h"

inline double individual_power_calc (double q_chisq, int df, 
				     double r2, double lambda) {
    //printf ("Computing power for %f, %d, %f * %f = %f\n", q_chisq, df, r2, lambda, r2*lambda);
  return(1. - (pnchisq(q_chisq,(double)df,r2*lambda,TRUE,0)));
}

/*
 * Function: calculate_r2
 * Date:     March 17, 2005 [RJK]
 * Purpose:  Given pointers to two SNPs, calculates the pairwise r2 between
 *           them.
 */
double calculate_r2 (snp_info_t *snp1, snp_info_t *snp2, int num_indivs) {
  unsigned char *gt1, *gt2;
  int counts[9];
  int N;
  double p,q, pq, oldpq;
  int i;
  int counter = 0;
  double r2;
  
  /* If the two SNPs are the same, force r2 to 1.0 */
  if (snp1 == snp2) return(1.0);

  gt1 = snp1->gt;
  gt2 = snp2->gt;

  for (i=0; i<9; i++) counts[i] = 0;

  for (i=0; i<num_indivs; i++) {
    if (gt1[i] < 1 || gt2[i] < 1 || gt1[i] > 3 || gt2[i] > 3)
      continue;
    counts[(gt1[i]-1)+3*(gt2[i]-1)]++;
  }
  
  N=0;
  for (i=0; i<9; i++) N+=counts[i];

  p = (2.0*(counts[0]+counts[1]+counts[2])+counts[3]+counts[4]+counts[5])/(2.0*N);
  q = (2.0*(counts[0]+counts[3]+counts[6])+counts[1]+counts[4]+counts[7])/(2.0*N);
  pq = p*q;

  do {
    oldpq = pq;
    counter++;
    if (1.0 - p - q + pq == 0.0 || pq == 0.0) {
      pq = p*q;
      break;
    }
    pq = (2*counts[0]+counts[1]+counts[3]+((counts[4]*pq*(1.0-p-q+pq))/(pq*(1.0-p-q+pq)+(p-pq)*(q-pq))))/(2.0*N);
  } while ((pq/oldpq < 0.9999 || pq/oldpq > 1.0001) && counter <= 100);
    
  if (p == 0.0 || q == 0.0 || p == 1.0 || q == 1.0) 
    r2 = 0.0;
  else
    r2 = (pq-(p*q))*(pq-(p*q))/(p*q*(1-p)*(1-q));
      
  return(r2);
}

/* 
 * Computes frequencies -- new formula from 4/9/07
 */
void compute_freqs (double *a, double *b, double *c, double *d, double *e, double *f, double p1, double p2, double g1, double g2, double K) {
  double Qcub, Acub, Bcub, Ccub;
  int res, i, good;
  double x[3];
  if (K > 0) {          /* Selected controls */
    Qcub = -1. + g1 + g2 - g1*g2;
    Acub = (g1 - 1.) * K + p1 - 2*g1*p1 -g1*p2+g2*(-1.+g1+K-g1*K-p1+2.*g1*p1+p2);
    Bcub = p1*(g2*(1.-K+g1*(-2.+2.*K-p1)-p2)+g1*(-1.*K+p1+p2));
    Ccub = g1*g2*p1*p1*(1.-K);
    res = gsl_poly_solve_cubic (Acub/Qcub, Bcub/Qcub, Ccub/Qcub,
				x, x+1, x+2);
    good = 0;
    for (i=0; i<res; i++) {
      *b = x[i];
      
      *d = x[i] * p2 / (x[i]-x[i]*g1+g1*p1);
      *f = (-1.*x[i]+x[i]*p1+x[i]*p2)/(-1.*x[i]+x[i]*g2-g2*p1);
      *a = p1-x[i];
      *c = p2-(x[i]*p2)/(x[i]-x[i]*g1+g1*p1);
      *e = 1.-p1-p2-(-1.*x[i]+x[i]*p1+x[i]*p2)/(-1.*x[i]+x[i]*g2-g2*p1);
      if (*a >= 0 && *a <= 1 && *b >= 0 && *b <= 1 && *c >= 0 && *c <= 1 && *d >= 0 && *d <= 1 && *e >= 0 && *e <= 1 && *f >= 0 && *f <= 1) {
	good = 1;
	break;
      }
    }
    if (good == 0) {
      Die ("Could not compute freqs for p1 = %f, p2 = %f, g1 = %f, g2 = %f, K = %f\n", p1, p2, g1, g2, K);
    }
    *a /= K;
    *c /= K;
    *e /= K;
    *b /= (1.-K);
    *d /= (1.-K);
    *f /= (1.-K);
  } else {
    *b = p1;
    *d = p2;
    *f = 1.-p1-p2;
    *a = p1/(-1.*g2-p1+g2*p1-g1*p2+g2*p2);
    *c = g1*p2/(g2+p1-g2*p1+g1*p2-g2*p2);
    *e = (g2-g2*p1-g2*p2)/(-1.*g2-p1+g2*p1-g1*p2+g2*p2);
  }
}

/*
 * Function: lambda_calc
 *
 * Given df and all parameters needed, calculates lambda for the given
 * sample size and a multiplicative model
 *
 */
double lambda_calc (int df, float n_cases, float n_controls, float grr, 
		    float heterozygote_f, float homozygote_f, int grr_type, float prevalence) {
  double a,b,c,d,e,f;
  double pa,pu;

  switch (grr_type) {
  case 0: 
    compute_freqs (&a, &b, &c, &d, &e, &f, 1.-heterozygote_f-homozygote_f, heterozygote_f, grr, grr*grr, prevalence);
    break;
  case 1:
    compute_freqs (&a, &b, &c, &d, &e, &f,  1.-heterozygote_f-homozygote_f, heterozygote_f, grr, grr+grr, prevalence);
    break;
  case 2:
    compute_freqs (&a, &b, &c, &d, &e, &f,  1.-heterozygote_f-homozygote_f, heterozygote_f, grr, grr, prevalence);
    break;
  case 3:
    compute_freqs (&a, &b, &c, &d, &e, &f,  1.-heterozygote_f-homozygote_f, heterozygote_f, 1, grr, prevalence);
    break;
  default:
    Die("Bad GRR type %d\n", grr_type);
    break;
  }
  if (df == 1) {
    pa = a+0.5*c;
    pu = b+0.5*d;
    return (2.*n_cases*n_controls*(pa-pu)*(pa-pu)*(n_cases+n_controls)/((n_cases*pa + n_controls*pu)*(n_cases+n_controls-n_cases*pa-n_controls*pu)));
  } else {
    return(n_cases*n_controls*(((a-b)*(a-b)/(n_cases*a+n_controls*b))+((c-d)*(c-d)/(n_cases*c+n_controls*d))+((e-f)*(e-f)/(n_cases*e+n_controls*f))));
  }
}

/*
 * Function: power
 *
 * Given an array for lambda, degrees of freedom, and snp list with best r^2 and
 * marker status, calculates the power.  Assumes total number of SNPs as
 * specified in structs.h.  If 0, calculates power for all SNPs equally.  If
 * not, assumes non-markers are representative of total number - # markers,
 * and then does markers seperately.
 */
double power (float n_cases, float n_controls, float grr, 
	      int df, double q_chisq,
	      snp_info_t **snp_info, int number_snps,
	      int total_snps, float min_r2, int grr_type, float prevalence) {
  int i;
  double power = 0.0;   /* Use first as power_0 */
  double marker_power = 0.0, non_marker_power = 0.0;
  int num_non_markers = 0, num_markers = 0;

  power = 1. - pnchisq(q_chisq, (double)(df), 0., TRUE, 0);

  for (i=0; i<number_snps; i++) {
    if (snp_info[i]->best_marker_r2 >= min_r2) {
      snp_info[i]->power = individual_power_calc(q_chisq, df, snp_info[i]->best_marker_r2, lambda_calc (df, n_cases, n_controls, grr, snp_info[i]->heterozygote_f, snp_info[i]->minor_homozygote_f, grr_type, prevalence));
    } else {
      snp_info[i]->power = power;
    }
//    printf ("POWER CALC\t%d\t%d\t%f\n", i, snp_info[i]->status, snp_info[i]->power); 
    if (snp_info[i]->status != 1) {
      num_non_markers++;
    } else {
      num_markers++;
    }
  }

  if (total_snps == 0)
    total_snps = num_markers + num_non_markers;

  power = 0.0;         /* Now use to calculate total power */

  for (i = 0; i<number_snps; i++) {
    if (snp_info[i]->status != 1) {
      non_marker_power += (snp_info[i]->power)/((float)(num_non_markers));
    } else {
      marker_power += (snp_info[i]->power)/((float)(num_markers));
    }
  }
  power = 
    (float)(num_markers)/(total_snps)*marker_power +
    (float)(total_snps-num_markers)/(total_snps)*non_marker_power;
  return(power);
}


