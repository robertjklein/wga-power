/*
 * structs.h
 *
 * Main structures for wga-power
 * 
 * Robert J. Klein
 * May 12, 2005
 */

#ifndef STRUCTSH_INCLUDED
#define STRUCTSH_INCLUDED

#include "squid.h"

#define MAX_INDIVS 90
#define MAX_NUMBER_SNPS 5000000 
#define MAX_RS_NUMBER 36223210
#define MAXDISTANCE 300000

#define BEST_N_DF 1
#define BEST_N_GRR 1.5
#define BEST_N_CASE 1000.
#define BEST_N_CONTROL 1000.

#define ALPHA 0.05

typedef struct _snp_info_t {
  unsigned int chr;
  unsigned int  pos;
  unsigned int status;
  unsigned int first_possible_marker, last_possible_marker;
  unsigned char gt[MAX_INDIVS];
  unsigned int  rs;
  int best_marker;
  float best_marker_r2;
  float *r2_values;
  double power;
  float heterozygote_f;
  float minor_homozygote_f;
} snp_info_t;

#endif
