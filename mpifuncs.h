/*
 * mpifuncs.h
 *
 * Basic functions for using MPI in wga-power.  Based on mpifuncs from
 * rsearch
 *
 * Robert J. Klein
 * May 12, 2005
 */

#ifndef _MPIFUNCS_H
#define _MPIFUNCS_H

#ifdef USE_MPI

#include "mpi.h"
#include "structs.h"

/* Get rank of master process (lowest ranked one that can do I/O */
/* Also checks the version string */
int get_master_rank (MPI_Comm comm, int mpi_my_rank);

void broadcast_data (int *marked_snps, int *num_markers, int *number_snps, int *number_indivs, int *total_snps, float *min_r2, snp_info_t ***snp_info, int mpi_my_rank, int mpi_master_rank);

#endif

#endif
