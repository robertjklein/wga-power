/*
 * mpifuncs.c
 *
 * Basic functions for using MPI in rsearch.
 *
 * Robert J. Klein
 * May 28, 2002
 */

#ifdef USE_MPI

#include <string.h>
#include <unistd.h>

#include "mpi.h"
#include "mpifuncs.h"
#include "structs.h"

#define VERSION_STRING "WGA-POWER 0.1"
#define VERSION_STRING_SIZE 100

/*
 * Function: get_master_rank
 * Date:     RJK, Tue May 28, 2002 [St. Louis]
 * Purpose:  Given a communicator, returns the lowest ranked one that can
 *           do I/O.
 *           Also checks the version string -- makes sure that it's the same
 *           on all procs.
 */
int get_master_rank (MPI_Comm comm, int mpi_my_rank) {
  int *io_proc_rank_p;
  int i;
  char versionbuf[VERSION_STRING_SIZE];
  
  MPI_Attr_get (comm, MPI_IO, &io_proc_rank_p, &i);
  
  if (i == 0)                 /* Not MPI compliant */
    return(MPI_PROC_NULL);
  if (*io_proc_rank_p == MPI_PROC_NULL)
    return (MPI_PROC_NULL);
  if (*io_proc_rank_p == MPI_ANY_SOURCE)
    return (0);

  /* Take min of procs allowed to do I/O.  */
  MPI_Allreduce (io_proc_rank_p, &i, 1, MPI_INT, MPI_MIN, comm);

  /* i is now master rank */
  /* Broadcast the version from master rank */
  if (i == mpi_my_rank)
    strncpy (versionbuf, VERSION_STRING, VERSION_STRING_SIZE-1);
  MPI_Bcast (versionbuf, VERSION_STRING_SIZE, MPI_CHAR, i, comm);
  
  if (strncmp (versionbuf, VERSION_STRING, VERSION_STRING_SIZE))
    Die ("Version strings %s and %s don't match\n", versionbuf, VERSION_STRING);
  return (i);
}

/*
 * Function: broadcast_data
 * Date:     RJK, Tue May 28, 2002 [St. Louis]
 *           RJK, Thu May 12, 2005 [NYC]
 * Purpose:  Broadcasts the first set of parameters needed in all processes
 *
 */
void broadcast_data (int *marked_snps, int *num_markers, int *number_snps, int *number_indivs, int *total_snps, float *min_r2, snp_info_t ***snp_info, int mpi_my_rank, int mpi_master_rank) {

  unsigned int *buf;
  float *buf2;
  int i;
  snp_info_t *cur_snp;

  /* Broadcast the two ints to everyone */
  MPI_Bcast (marked_snps, 1, MPI_INT, mpi_master_rank, MPI_COMM_WORLD);
  MPI_Bcast (num_markers, 1, MPI_INT, mpi_master_rank, MPI_COMM_WORLD);
  MPI_Bcast (number_snps, 1, MPI_INT, mpi_master_rank, MPI_COMM_WORLD);
  MPI_Bcast (number_indivs, 1, MPI_INT, mpi_master_rank, MPI_COMM_WORLD);
  MPI_Bcast (total_snps, 1, MPI_INT, mpi_master_rank, MPI_COMM_WORLD);
  MPI_Bcast (min_r2, 1, MPI_FLOAT, mpi_master_rank, MPI_COMM_WORLD);

  /*
   * Now we need to broadcast the SNP information 
   */
  
  /* First, set up the SNP array for all but master rank */
  if (mpi_my_rank != mpi_master_rank) {
    *snp_info = MallocOrDie(sizeof(snp_info_t *)*(*number_snps));
    for (i=0; i<(*number_snps); i++) {
      (*snp_info)[i] = MallocOrDie(sizeof(snp_info_t));
      cur_snp = (*snp_info)[i];

      /* Defaults */
      cur_snp -> best_marker = -1;
      cur_snp -> best_marker_r2 = -1.;
      cur_snp -> r2_values = NULL;
      cur_snp -> power = 0.;
      cur_snp -> status = 0;
    }
  }

  //sleep(120);
  /* Now, make the buffer for unsigned ints of rs, chr, pos, first, last, status */
  buf = MallocOrDie(sizeof(unsigned int)*(*number_snps)*6);
  if (mpi_master_rank == mpi_my_rank) {
    for (i=0; i<(*number_snps)*6; i+=6) {
      buf[i] = (*snp_info)[i/6]->rs;
      buf[i+1] = (*snp_info)[i/6]->chr;
      buf[i+2] = (*snp_info)[i/6]->pos;
      buf[i+3] = (*snp_info)[i/6]->first_possible_marker;
      buf[i+4] = (*snp_info)[i/6]->last_possible_marker;
      buf[i+5] = (*snp_info)[i/6]->status;
    }
  }
  //printf ("Here in process %d\n", mpi_my_rank);
  MPI_Bcast (buf, (*number_snps)*6, MPI_UNSIGNED, mpi_master_rank, MPI_COMM_WORLD);
  if (mpi_master_rank != mpi_my_rank) {
    for (i=0; i<(*number_snps)*6; i+=6) {
      (*snp_info)[i/6]->rs = buf[i];
      (*snp_info)[i/6]->chr = buf[i+1];
      (*snp_info)[i/6]->pos = buf[i+2];
      (*snp_info)[i/6]->first_possible_marker = buf[i+3];
      (*snp_info)[i/6]->last_possible_marker = buf[i+4];
      (*snp_info)[i/6]->status = buf[i+5];
    }
  }
  free (buf);

  /* Now, do a broadcast for each gt in each snp_info */
  for (i=0; i<(*number_snps); i++) {
    MPI_Bcast ((*snp_info)[i]->gt, MAX_INDIVS, MPI_UNSIGNED_CHAR, mpi_master_rank, MPI_COMM_WORLD);
  }

  /* Now do minor_homozygote_f and heterozygote_f */
  buf2 = MallocOrDie(sizeof(float)*(*number_snps)*2);
  if (mpi_master_rank == mpi_my_rank) {
    for (i=0; i<(*number_snps)*2; i+=2) {
      buf2[i] = (*snp_info)[i/2]->minor_homozygote_f;
      buf2[i+1] = (*snp_info)[i/2]->heterozygote_f;
    }
  }
  //printf ("Here in process %d\n", mpi_my_rank);
  MPI_Bcast (buf2, (*number_snps)*2, MPI_FLOAT, mpi_master_rank, MPI_COMM_WORLD);
  if (mpi_master_rank != mpi_my_rank) {
    for (i=0; i<(*number_snps)*2; i+=2) {
      (*snp_info)[i/2]->minor_homozygote_f = buf2[i];
      (*snp_info)[i/2]->heterozygote_f = buf2[i+1];
    }
  }
  free (buf2);

  /* And we're done */
}

#endif
