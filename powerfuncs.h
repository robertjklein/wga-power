/*
 * powerfuncs.h
 * 
 * Headers for powerfuncs.c
 */

#ifndef POWERFUNCSH_INCLUDED
#define POWERFUNCSH_INCLUDED

#include "structs.h"

double lambda_calc (int df, float n_cases, float n_controls, float grr,
                    float heterozygote_f, float homozygote_f, int grr_type, float prevalence);

inline double individual_power_calc (double q_chisq, int df, 
				     double r2, double lambda);

double calculate_r2 (snp_info_t *snp1, snp_info_t *snp2, int num_indivs);

double power (float n_cases, float n_controls, float grr,
	      int df, double q_chisq,
	      snp_info_t **snp_info, int number_snps,
	      int total_snps, float min_r2, int grr_type, float prevalence);


#endif
