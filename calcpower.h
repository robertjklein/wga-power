/*
 * calcpower.h
 *
 * Routines for calculating power from a set of markers
 *
 * Robert J. Klein
 * May 12, 2005
 */

#ifndef CALCPOWERH_INCLUDED
#define CALCPOWERH_INCLUDED

#include "structs.h"

void best_r2 (snp_info_t **snp_info, int number_snps, int number_indivs
#ifdef USE_MPI
	      , int mpi_my_rank, int mpi_num_procs
#endif
	      );

#endif
