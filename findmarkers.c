/**
 * findmarkers.c
 *
 * Routines involved in finding best SNP markers
 *
 * Robert J. Klein
 * May 12, 2005
 */

#include <stdio.h>
#include <stdlib.h>

#define MATHLIB_STANDALONE
#include <Rmath.h>

#include "structs.h"
#include "powerfuncs.h"
#include "findmarkers.h"

#ifdef USE_MPI
#include "mpi.h"
#endif



/*
 * For the given SNP, calculates the change in power and returns it
 */
double try_snp (snp_info_t **snp_info, int snp, int number_indivs,
		int num_markers, int total_snps, int total_gtyped_snps,
		double q_chisq, double power_0, float min_r2, int grr_type) {
  int start, stop;
  int j;
  snp_info_t *snp1, *snp2;
  double r2;
  double power_diff = 0.;
  double local_power_diff;

  snp1 = snp_info[snp];
  start = snp1->first_possible_marker;
  stop = snp1->last_possible_marker;

  for (j=start; j<=stop; j++) {
    snp2 = snp_info[j];
    r2 = snp1->r2_values[j-start];
    //r2 = calculate_r2(snp1,snp2,number_indivs);
    if (r2 > snp2->best_marker_r2) {
      if (r2 < min_r2) {
	local_power_diff = (power_0 - snp2->power);
      } else {
	local_power_diff = (individual_power_calc(q_chisq, (double)(BEST_N_DF), r2, lambda_calc(BEST_N_DF, BEST_N_CASE, BEST_N_CONTROL, BEST_N_GRR, snp2->heterozygote_f, snp2->minor_homozygote_f, grr_type, -1.)) - snp2->power);
      }
      if (snp2->status == 1) {
	power_diff += local_power_diff/total_snps;
      } else {
	power_diff += local_power_diff*(total_snps-num_markers)/(total_snps*total_gtyped_snps);
      }
    }
  }
  return(power_diff);
}

/* 
 * Sets the best SNP as a marker and recalculates the power
 */
void set_snp (snp_info_t **snp_info, int snp, int number_indivs, double q_chisq,
	 double power_0, float min_r2, int grr_type) {
  int start, stop;
  int j;
  snp_info_t *snp1, *snp2;
  double r2;

  snp1 = snp_info[snp];
  start = snp1->first_possible_marker;
  stop = snp1->last_possible_marker;

  snp1 -> status = 1;

  for (j=start; j<=stop; j++) {
    snp2 = snp_info[j];
    r2 = calculate_r2(snp1, snp2, number_indivs);
    if (r2 > snp2->best_marker_r2) {
      snp2->best_marker_r2 = r2;
      snp2->best_marker = snp;
      if (r2 < min_r2) {
	snp2->power = power_0;
      } else {
	snp2->power = individual_power_calc(q_chisq, (double)(BEST_N_DF), r2, lambda_calc(BEST_N_DF, BEST_N_CASE, BEST_N_CONTROL, BEST_N_GRR, snp2->heterozygote_f, snp2->minor_homozygote_f, grr_type, -1.));
      }
    }
  }
}


/*
 * Repeatedly adds markers to the set of markers to be used in a greedy
 * algorithm.  For each candidate marker to add, recalculate power based on
 * the new marker.  Report the best marker and current power at specified
 * lambda (hard-coded) each time.
 */
void find_best_n_snps (snp_info_t **snp_info, int number_snps, 
		       int number_indivs, int marked_snps, int num_markers,
		       int total_snps, float min_r2, int grr_type
#ifdef USE_MPI
			 ,int mpi_my_rank, int mpi_num_procs, int mpi_master_rank
#endif
			 ) {
  int i,m, j;
  int cur_best_snp;
  double cur_best_power_change, cur_power_change;
  double final_power;
  double q_chisq, power_0;
  unsigned long long bytes_to_alloc = 0;
  float *r2_data;
#ifdef USE_MPI
  double *big_powerbuf, *small_powerbuf;
  int bigbufsize, smallbufsize;
  int k;
  struct {
    double power_change;
    int snp;
  } cur_best, overall_best;
#endif
  
  int root, stepsize;            /* 0, 1 for serial, mpi_my_rank, mpi_num_procs for MPI */

  /* 
   *First, find all r2 values 
   */
  
  /* For serial, find all values.  For MPI, only find values for some
     SNPs as i in each proc, overall all r2 is found */
#ifdef USE_MPI
  root = mpi_my_rank;
  stepsize = mpi_num_procs;
#else
  root = 0;
  stepsize = 1;
#endif

  /* Calculate number of bytes to allocate for r2 values and allocate them */
  for (i=root; i<number_snps; i+=stepsize) {
    bytes_to_alloc += sizeof(float)*(snp_info[i]->last_possible_marker-snp_info[i]->first_possible_marker+1);
  }
  r2_data = MallocOrDie(bytes_to_alloc);


  /* Now, calculate the r2 values */
  snp_info[root]->r2_values = r2_data;
  for (i=root+stepsize; i<number_snps; i+=stepsize) {
    snp_info[i]->r2_values = snp_info[i-stepsize]->r2_values + (snp_info[i-stepsize]->last_possible_marker-snp_info[i-stepsize]->first_possible_marker+1);
    for (j=snp_info[i]->first_possible_marker; j<=snp_info[i]->last_possible_marker; j++) {
      snp_info[i]->r2_values[j-snp_info[i]->first_possible_marker] = 
	(float)calculate_r2(snp_info[i], snp_info[j], number_indivs);
    }
  }

  /* If using MPI, allocate buffer space for the power calculations which
     will then be split among CPUs and Allgathered */
#ifdef USE_MPI
  smallbufsize = (number_snps/mpi_num_procs)+1;
  bigbufsize = smallbufsize * mpi_num_procs;
  big_powerbuf = MallocOrDie(sizeof(double)*bigbufsize);
  small_powerbuf = MallocOrDie(sizeof(double)*smallbufsize);
#endif
  


  for (m=marked_snps; m<num_markers; m++) {
    cur_best_power_change = 0;
    q_chisq = qchisq(1.0-(ALPHA/((double)(m+1))),(double)(BEST_N_DF),TRUE,0);
    power_0 = 1. - pnchisq(q_chisq,(double)(BEST_N_DF),0.,TRUE,0);
    
    /* First, recalculate power for all SNPs based on the current number of
       markers.  Split among procs for MPI */
#ifdef USE_MPI
    for (i=root, j=0; i<number_snps; i+=stepsize, j++) {
      if (snp_info[i]->best_marker_r2 < min_r2) {
	small_powerbuf[j] = power_0;
      } else {
	small_powerbuf[j] = individual_power_calc(q_chisq, (double)(BEST_N_DF), snp_info[i]->best_marker_r2, lambda_calc(BEST_N_DF, BEST_N_CASE, BEST_N_CONTROL, BEST_N_GRR, snp_info[i]->heterozygote_f, snp_info[i]->minor_homozygote_f, grr_type, -1.));
      }
    }
    MPI_Allgather (small_powerbuf, smallbufsize, MPI_DOUBLE, big_powerbuf, smallbufsize, MPI_DOUBLE, MPI_COMM_WORLD);
    /* Now, this loop is not the most efficient but it will get the job done */
    for (i=0; i<mpi_num_procs; i++) {
      for (j=i,k=i*smallbufsize; j<number_snps; j+=stepsize,k++) {
	snp_info[j]->power = big_powerbuf[k];
      }
    }
#else
    for (i=0; i<number_snps; i++) {
      if (snp_info[i]->best_marker_r2 < min_r2) {
	snp_info[i]->power = power_0;
      } else {
	snp_info[i]->power = individual_power_calc(q_chisq, (double)(BEST_N_DF), snp_info[i]->best_marker_r2, lambda_calc(BEST_N_DF, BEST_N_CASE, BEST_N_CONTROL, BEST_N_GRR, snp_info[i]->heterozygote_f, snp_info[i]->minor_homozygote_f, grr_type, -1.));
      }
    }
#endif
    

    /* This loop is only executed for some SNPs on each proc under MPI; again
       this is taken care of through the use of the root and stepsize variables */
    for (i=root; i<number_snps; i+=stepsize) {
      if (snp_info[i]->status != 0) {
	continue;
      }
      cur_power_change = try_snp (snp_info, i, number_indivs, m, total_snps, 
				  number_snps, q_chisq, power_0, min_r2, grr_type);
      if (cur_power_change > cur_best_power_change) {
	cur_best_power_change = cur_power_change;
	cur_best_snp = i;
      }
    }

    /* Here where it gets calculated.  Now we need to set cur_best_power_change
       and cur_best_snp to be the same based on whichever proc has the highest
       cur_best_power_change */
#ifdef USE_MPI
    cur_best.power_change = cur_best_power_change;
    cur_best.snp = cur_best_snp;
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Allreduce (&cur_best, &overall_best, 1, MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD);
    cur_best_power_change = overall_best.power_change;
    cur_best_snp = overall_best.snp;
#endif

    if (cur_best_power_change > 0.) {
      /* Do the setting algorithm on all CPUs so we have the data for later
	 This is easier than splitting up and broadcasting results to all, 
	 though I may change that later */
      set_snp (snp_info, cur_best_snp, number_indivs, q_chisq, power_0, min_r2, grr_type);
    
      /* This next command could be split over procs, but there's no need */
      final_power = power (BEST_N_CASE, BEST_N_CONTROL, BEST_N_GRR, 
			   BEST_N_DF, q_chisq, snp_info, number_snps, total_snps, min_r2, grr_type, -1.);
      /* And here, only put output on master rank proc */
#ifdef USE_MPI
      if (mpi_my_rank == mpi_master_rank) {
#endif
      printf ("%d: Added SNP rs%d for total power of %f\n", m, snp_info[cur_best_snp]->rs, final_power);
      fflush(stdout);
#ifdef USE_MPI
      }
#endif
    } else {
#ifdef USE_MPI
      if (mpi_my_rank == mpi_master_rank)
#endif
      printf ("%d: Power did not increase.  Stop here\n", m);
      break;
    }
  }

  free(r2_data);
#ifdef USE_MPI
  free(big_powerbuf);
  free(small_powerbuf);
#endif

}
      
