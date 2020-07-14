/*
 * calcpower.c
 *
 * Routines for calculating power from a set of markers
 *
 * Robert J. Klein
 * May 12, 2005
 */

#include <stdio.h>
#include <stdlib.h>

#include "structs.h"
#include "calcpower.h"
#include "powerfuncs.h"

#ifdef USE_MPI
#include "mpi.h"
#endif



/*
 * Given a file of SNP names, computes the best r2 for each SNP in the
 * initial genotyped set against the subset of SNPs in the file.
 */
void best_r2 (snp_info_t **snp_info, int number_snps, int number_indivs
#ifdef USE_MPI
	      , int mpi_my_rank, int mpi_num_procs
#endif
	      ) {
  double r2;
  float best_r2_val;
  int best_r2_snp = 0;
  int i,j;
  unsigned int last_index;
  snp_info_t *snp1, *snp2;
#ifdef USE_MPI
  struct {
    float best_r2_val_s;
    int marker;
  } float_int_data;



  for (i=mpi_my_rank; i<number_snps; i+=mpi_num_procs) {
#else
  for (i=0; i<number_snps; i++) {
#endif
    snp1 = snp_info[i];
    best_r2_val = -1.;
    if (snp1->chr != 0) {
      last_index = snp1->last_possible_marker;
      for (j=snp1->first_possible_marker; j<=last_index; j++) {
	snp2 = snp_info[j];
	if (snp2->status != 1) continue;
	r2 = calculate_r2(snp1, snp2, number_indivs);
	if ((float)r2 > best_r2_val) {
	  best_r2_val = (float)r2;
	  best_r2_snp = j;
	}
      }
    }
    snp_info[i]->best_marker_r2 = best_r2_val;
    snp_info[i]->best_marker = best_r2_snp;
  }

#ifdef USE_MPI
  /* Now, broadcast results using MPI for each SNP */
  for (i=0; i<number_snps; i++) {
    if (i % mpi_num_procs == mpi_my_rank) {
      float_int_data.best_r2_val_s = snp_info[i]->best_marker_r2;
      float_int_data.marker = snp_info[i]->best_marker;
    }
    MPI_Bcast(&float_int_data, 1, MPI_FLOAT_INT, i % mpi_num_procs, MPI_COMM_WORLD);
    if (i % mpi_num_procs != mpi_my_rank) {
      snp_info[i]->best_marker_r2 = float_int_data.best_r2_val_s;
      snp_info[i]->best_marker = float_int_data.marker;
    }
  }
#endif
}


