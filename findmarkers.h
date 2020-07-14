/*
 * findmarkers.h
 *
 * Headers for routines involved in finding the best SNP markers to use
 *
 * Robert J. Klein
 * May 12, 2005
 */

#ifndef FINDMARKERSH_INCLUDED
#define FINDMARKERSH_INCLUDED


#include "structs.h"

void find_best_n_snps (snp_info_t **snp_info, int number_snps, 
		       int number_indivs, int marked_snps, int num_markers,
		       int total_snps, float min_r2, int grr_type
#ifdef USE_MPI
			 ,int mpi_my_rank, int mpi_num_procs, int mpi_master_rank
#endif
			 );

#endif

