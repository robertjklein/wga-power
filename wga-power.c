/* wga-power.c
 *
 * Code for dealing with whole-genome association power calculations, using
 * LD inferred from genotype data for a given population.
 *
 * Robert J. Klein
 * Started March 17, 2005
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <ctype.h>

#define MATHLIB_STANDALONE
#include <Rmath.h>

#include "squid.h"		/* general sequence analysis library    */
#include "structs.h"
#include "powerfuncs.h"
#include "findmarkers.h"
#include "calcpower.h"

#ifdef USE_MPI
#include "mpi.h"
#include "mpifuncs.h"

static int in_mpi;
#endif


static char banner[] = "wga-power -- calculates power for whole genome association\n";

static char usage[]  = "\
Usage: wga-power [-options] <LD genotype data>\n\
  Available options are:\n\
   -h          : help; print brief help on version and usage\n\
   -a <file>   : Read in list of SNPs genotyped from file\n\
                 Start with this list of SNPs in data set\n\
                 Compute best r2 with this list for each SNP in LD data\n\
   -n <num>    : Compute best set of <num> markers for genotype data\n\
                 Only happens if -n number is greater than -a number\n\
   -k <float>  : Fraction of individuals who are cases (default 0.5)\n\
   -N <string> : Start,stop,step for number of individuals\n\
                 (default 100,2000,100)\n\
";

static char experts[] = "\
   --totalsnps       : total number of SNPs assumed\n\
   --printeachmarker : prints the marker for each SNP\n\
   --grr <string>    : grr1,grr2,.... to try (default 1.5,2,4)\n\
   --qchisq <num>,<num> : chi^2 cutoff for 1 df and 2df\n\
   --minr2 <num>     : Min. r^2 value to allow (default 0.1)\n\
   --grrtype <num>   : Type of GRR to calculate (default 0=mult, 1=add, 2=dom, 3=rec\n\
   --prevalence <num> : Prevalence of the disease (implies selected controls*\n\
";

static struct opt_s OPTIONS[] = {
  { "-h", TRUE, sqdARG_NONE }, 
  { "-a", TRUE, sqdARG_STRING },
  { "-n", TRUE, sqdARG_INT },
  { "-k", TRUE, sqdARG_FLOAT },
  { "-N", TRUE, sqdARG_STRING},
  { "--grr", FALSE, sqdARG_STRING},
  { "--grrtype", FALSE, sqdARG_INT},
  { "--totalsnps", FALSE, sqdARG_INT },
  { "--printeachmarker", FALSE, sqdARG_NONE},
  { "--qchisq", FALSE, sqdARG_STRING},
  { "--minr2", FALSE, sqdARG_FLOAT},
  { "--prevalence", FALSE, sqdARG_FLOAT}
};
#define NOPTIONS (sizeof(OPTIONS) / sizeof(struct opt_s))




int snp_sort_func (const void *a, const void *b) {
  snp_info_t *x;
  snp_info_t *y;
  x = *(snp_info_t **)a;
  y = *(snp_info_t **)b;

  if (x->chr < y->chr) {
    return (-1);
  } else if (x->chr > y->chr) {
    return(1);
  } else if (x->pos < y->pos) {
    return(-1);
  } else if (x->pos > y->pos) {
    return(1);
  } else {
    return(0);
  }
}



/*
 * Function: read_genotypes
 * Date:     March 17, 2005 [RJK]
 * Purpose:  Given a file of SNP genotypes, reads into the data structure
 */
int read_genotypes (char *gt_file, snp_info_t **snp_info, int *number_indivs_p) {
  FILE *f;
  char buf[65536];
  char *cp;
  int i = 0;
  int j;
  int number_indivs = 0;
  snp_info_t *cur_snp;
  int num_snps;
  int count1, count2, count3;

  f = fopen(gt_file,"r");
  if (f==NULL) Die("Could not open %s\n", gt_file);

  while (fgets(buf,65535,f) != 0) {
    if (i > MAX_NUMBER_SNPS) {
      Die("There are more than %d SNPs\n", MAX_NUMBER_SNPS);
    }

    cur_snp = MallocOrDie(sizeof(snp_info_t));

    /* First column -- SNP number */
    cur_snp->rs = atoi(buf);
    if (cur_snp->rs >= MAX_RS_NUMBER) {
      Die("SNP has rs of %d for line %s\n", cur_snp->rs, buf);
    }
    for (cp=buf; *cp != '\0' && !isspace(*cp); cp++);
    while (*cp != '\0' && isspace(*cp)) cp++;

    if (isdigit(*cp)) {
      cur_snp->chr = atoi(cp);
    } else if (*cp == 'X') {
      cur_snp->chr = 23;
    } else {
      cur_snp->chr = 0;
    }
    while(*cp != '\0' && !isspace(*cp)) cp++;
    while (*cp != '\0' && isspace(*cp)) cp++;

    cur_snp->pos = atoi(cp);
    while (*cp != '\0' && !isspace(*cp)) cp++;
    while (*cp != '\0' && isspace(*cp)) cp++;

    number_indivs = 0;

    /* Initialize with pseudo-counts; 1 for AA and BB, 2 for AB */
    count1 = 1;
    count2 = 2;
    count3 = 1;
    do {
      if (number_indivs >= MAX_INDIVS) {
	Die("Not enough individuals allocated in .h\n");
      }
      if (cp[0] == 'A' && cp[1] == 'A') {
	cur_snp->gt[number_indivs] = 1;
	count1 ++;
      } else if (cp[0] == 'A' && cp[1] == 'B') {
	cur_snp->gt[number_indivs] = 2;
	count2++;
      } else if (cp[0] == 'B' && cp[1] == 'B') {
	cur_snp->gt[number_indivs] = 3;
	count3++;
      } else {
	cur_snp->gt[number_indivs] = 0;
      }
      while (*cp != '\0' && !isspace(*cp)) cp++;
      while (isspace(*cp) && *cp != '\0') cp++;
      number_indivs++;
    } while (*cp != '\0');

    cur_snp->best_marker = -1;
    cur_snp->best_marker_r2 = -1.0;
    cur_snp->r2_values = NULL;
    cur_snp->power = 0.;
    cur_snp->status = 0;

    cur_snp->heterozygote_f = 1.*count2/(count1+count2+count3);
    if (count1 < count3) {
      cur_snp->minor_homozygote_f = 1.*count1/(count1+count2+count3); 
    } else {
      cur_snp->minor_homozygote_f = 1.*count3/(count1+count2+count3);
    }

    if (cur_snp->pos == 0) {
      cur_snp->chr = 0;      /* If undefined position, can't get LD info */
    }

    snp_info[i] = cur_snp;
    i++;
  }
  *number_indivs_p = number_indivs;

  num_snps = i;
  for (i=0; i<num_snps; i++) {
    if (snp_info[i] == NULL) {
      Die ("WHOA! For index %d, NULL\n", i);
    }
  }

  qsort (snp_info, num_snps, sizeof(snp_info_t *), &snp_sort_func);


  /* Now, set first and last possible marker */
  for (i=0; i<num_snps; i++) {
    if (snp_info[i]->chr == 0) {
      snp_info[i]->first_possible_marker = i;
      snp_info[i]->last_possible_marker = i;
    } else {
      j = i;
      while (j >= 0) {
	if (snp_info[i] -> chr != snp_info[j]->chr) break;
	if (snp_info[i]->pos - snp_info[j]->pos > MAXDISTANCE) break;
	j--;
      }
      j++;
      snp_info[i]->first_possible_marker = j;

      j = i;
      while (j<num_snps) {
	if (snp_info[i]->chr != snp_info[j]->chr) break;
	if (snp_info[j]->pos - snp_info[i]->pos > MAXDISTANCE) break;
	j++;
      }
      j--;
      snp_info[i]->last_possible_marker = j;
      //      printf ("Marker rs%d (%d:%d), range from %d:%d to %d:%d\n", snp_info[i]->rs, snp_info[i]->chr, snp_info[i]->pos, snp_info[snp_info[i]->first_possible_marker]->chr, snp_info[snp_info[i]->first_possible_marker]->pos, snp_info[snp_info[i]->last_possible_marker]->chr, snp_info[snp_info[i]->last_possible_marker]->pos);
    }
  }

  return(num_snps);
}



/*
 * Function: mark_exact_dups
 * Date:     RJK, Wed May 18, 2005 [NYC]
 * Purpose:  Given a list of SNPs, find SNPs within MAXDISTANCE of each other
 *           that are in perfect LD and marks the lowest numbered one as status
 *           0, rest as status 2
 *
 *           Keys each SNP by its 
 */
void mark_exact_dups (snp_info_t **snp_info, int number_snps, 
		      int number_indivs) {
  int i,j,k;
  int mismatch = 0;

  for (i=0; i<number_snps; i++) {
    if (snp_info[i]->status == 2) continue;
    for (j=snp_info[i]->first_possible_marker; j<=snp_info[i]->last_possible_marker; j++) {
      if (i == j) continue;
      if (snp_info[j]->status == 2) continue;
      mismatch = 0;
      for (k=0; k<number_indivs; k++) {
	if (snp_info[i]->gt[k] != snp_info[j]->gt[k]) {
	  mismatch++;
	  break;
	}
      }
      if (mismatch == 0) {
	if (i < j) {
	  snp_info[j]->status = 2;
	} else {
	  Die("Shouldn't be here rs%d < rs%d\n", snp_info[i]->rs, snp_info[j]->rs);
	}
      }
    }
  }
}


/*
 * Given a file of SNP names, finds them in snp_info list and puts their
 * index into array that is returned.
 * as 1.
 */
int mark_snps_to_use (char *snps_file, snp_info_t **snp_info, int num_snps) {
  int *snp_indices;
  int i;
  char buf[256];
  FILE *f;
  int retval = 0;

  if (snps_file != NULL) {
 
    snp_indices = MallocOrDie(sizeof(int)*MAX_RS_NUMBER);
    for (i=0; i<MAX_RS_NUMBER; i++) {
      snp_indices[i] = -1;
    }
    
    for (i=0; i<num_snps; i++) {
      if (snp_info[i] == NULL) break;
      snp_indices[snp_info[i]->rs] = i;
    }
    
    f = fopen(snps_file, "r");
    while (fgets(buf, 255, f)) {
      if (buf[0] == 'r' && buf[1] == 's') {
	i = snp_indices[atoi(buf+2)];
      } else {
        i = snp_indices[atoi(buf)];
      }
      if (i != -1) {
	snp_info[i]->status = 1;
	retval++;
      }
    }
    fclose(f);

    free(snp_indices);
  }
  return(retval);
}


#ifdef USE_MPI
/*
 * Function: exit_from_mpi
 * Date:     RJK, Thu Jun 6, 2002 [St. Louis]
 * Purpose:  Calls MPI_Abort on exit if in_mpi flag is 1, otherwise
 *           returns
 */
void exit_from_mpi () {
  if (in_mpi)
    MPI_Abort (MPI_COMM_WORLD, -1);
}

#endif


int
main(int argc, char **argv)
{
  char            *gtfile;       /* File of genotypes for calculating LD */
  snp_info_t      **snp_info;
  int             i;

  char            *snps_file = NULL;    /* List of SNPs in genotyping set */
  int             num_markers = 0;      /* Number of markers to find, or if
					   using a set list, number there */


  int             number_indivs;   /* How many individuals in set */
  int             number_snps;     /* Total number of SNPs */
  int             marked_snps;     /* Number of marker SNPs */

  char *optname;                /* name of option found by Getopt()        */
  char *optarg;                 /* argument found by Getopt()              */
  int   optind;                 /* index in argv[]                         */


  int total_snps = 0;            /* Total number of SNPs assumed in set */
  int print_each_marker = 0;     /* If set to 1, print each marker */

  float frac_case = 0.5;
  int N_start = 100;
  int N_stop = 2000;
  int N_step = 100;
  int N;
  float GRRs[100];
  int num_GRRs = 3;
  GRRs[0] = 1.5;
  GRRs[1] = 2;
  GRRs[2] = 4;
  int GRR_index;
  float minr2 = 0.1;
  int grr_type = 0;    /* default is mult. model */
  float prevalence = -1.; /* default is unselected controls */

  char buf[65535];
  char *bufp;

  double q_chisq1 = -1., q_chisq2 = -1.;

#ifdef USE_MPI
  int mpi_my_rank;              /* My rank in MPI */
  int mpi_num_procs;            /* Total number of processes */
  int mpi_master_rank;          /* Rank of master process */

  /* Initailize MPI, get values for rank and num procs */
  MPI_Init (&argc, &argv);

  atexit (exit_from_mpi);
  in_mpi = 1;                /* Flag for exit_from_mpi() */

  MPI_Comm_rank (MPI_COMM_WORLD, &mpi_my_rank);
  MPI_Comm_size (MPI_COMM_WORLD, &mpi_num_procs);

  /*
   * Determine master process.  This is the lowest ranking one that can do I/O
   */
  mpi_master_rank = get_master_rank (MPI_COMM_WORLD, mpi_my_rank);

  /* If I'm the master, do the following set up code -- parse arguments, read
     in matrix and query, build model */
  if (mpi_my_rank == mpi_master_rank) {
#endif


  /**********************************************
   * Print header here 
   *********************************************/
#define VERSION "0.1"
#define COPYRIGHT "Development"
#define LICENSE "Do not redistribute"

  printf ("WGA-power version %s\n", VERSION);
#ifdef USE_MPI
  printf ("Running in parallel with %d processes\n", mpi_num_procs);
#endif
  printf ("%s\n%s\n\n", COPYRIGHT, LICENSE);


  /*********************************************** 
   * Parse command line
   ***********************************************/
  
  while (Getopt(argc, argv, OPTIONS, NOPTIONS, usage,
                &optind, &optname, &optarg))  {
    if (strcmp(optname, "-a") == 0) {
      snps_file = MallocOrDie(strlen(optarg)+1);
      strncpy (snps_file, optarg, strlen(optarg));
      snps_file[strlen(optarg)] = '\0';
    } else if (strcmp(optname, "-n") == 0) {
      num_markers = atoi(optarg);
    } else if (strcmp(optname, "-k") == 0) {
      frac_case = atof(optarg);
    } else if (strcmp(optname, "-N") == 0) {
      strncpy (buf, optarg, strlen(optarg));
      buf[strlen(optarg)]='\0';
      N_start = atoi(buf);
      bufp = buf;
      while (*bufp != ',') bufp++;
      bufp++;
      N_stop = atoi(bufp);
      while (*bufp!= ',') bufp++;
      bufp++;
      N_step = atoi(bufp);
    } else if (strcmp(optname, "--totalsnps") == 0) {
      total_snps = atoi(optarg);
    } else if (strcmp(optname, "--minr2") == 0) {
      minr2 = atof(optarg);
    } else if (strcmp(optname, "--grrtype") == 0) {
      grr_type = atoi(optarg);
    } else if (strcmp(optname, "--printeachmarker") == 0) {
      print_each_marker = 1;
    } else if (strcmp(optname, "--grr") == 0) {
      strncpy (buf, optarg, strlen(optarg));
      buf[strlen(optarg)] = '\0';
      GRR_index = 0;
      bufp = buf;
      do {
	GRRs[GRR_index] = atof(bufp);
	num_GRRs = GRR_index+1;
	while (*bufp != '\0' && *bufp != ',') {
	  bufp++;
	}
	if (*bufp == ',') {
	  bufp++;
	  GRR_index++;
	}
      } while (*bufp != '\0');
    } else if (strcmp(optname, "--qchisq") == 0) {
      strncpy (buf, optarg, strlen(optarg));
      buf[strlen(optarg)] = '\0';
      for (bufp = buf; *bufp != '\0' && *bufp != ','; bufp++);
      if (*bufp == ',') {
	*bufp = '\0';
	q_chisq1 = atof(buf);
	q_chisq2 = atof(bufp+1);
      }
    } else if (strcmp(optname, "--prevalence") == 0) {
      prevalence = atof(optarg);
    } else if (strcmp(optname, "-h") == 0) {
      puts(usage);
      puts(experts);
      exit(EXIT_SUCCESS);
    }
  }

  if (argc - optind != 1) 
    Die("Incorrect number of arguments.\n%s\n", usage);
  gtfile = argv[optind++];

  /* Initialize SNP info */
  snp_info = MallocOrDie(sizeof(snp_info_t *)*(MAX_NUMBER_SNPS+1));
  for (i=0; i<=MAX_NUMBER_SNPS; i++)
    snp_info[i] = NULL;

  /* Now, read in the genotypes, and then sort them by chromosome number
   and position.  Then, mark those SNPs in snps_to_use file and 
   exact dups file */
  number_snps = read_genotypes (gtfile, snp_info, &number_indivs);
  mark_exact_dups (snp_info, number_snps, number_indivs);
  marked_snps = mark_snps_to_use (snps_file, snp_info, number_snps);

#ifdef USE_MPI
  }   /* End of first block that is only done by master process */
  /* Barrier for debugging */
  MPI_Barrier(MPI_COMM_WORLD);

  /* Here, we need to broadcast the data we've read in: marked_snps,
     num_markers, number_snps, snp_info */
  broadcast_data (&marked_snps, &num_markers, &number_snps, &number_indivs, &total_snps, &minr2, &snp_info, mpi_my_rank, mpi_master_rank);

#endif

  /* Now, find best r2 for marked SNPs */
  /* best_r2 is split over several procs if using MPI */
  best_r2 (snp_info, number_snps, number_indivs
#ifdef USE_MPI
	   , mpi_my_rank, mpi_num_procs
#endif
	   );
  

  /* If we have a number of markers to find, find them now */
  if (num_markers > marked_snps)
    find_best_n_snps (snp_info, number_snps, number_indivs, marked_snps, num_markers, total_snps, minr2, grr_type
#ifdef USE_MPI
		    ,mpi_my_rank, mpi_num_procs, mpi_master_rank
#endif
		    );
  else
    num_markers = marked_snps;
  

  /* Now, find best r2 again with the new data */
  best_r2 (snp_info, number_snps, number_indivs
#ifdef USE_MPI
	   , mpi_my_rank, mpi_num_procs
#endif
	   );

  
  /* Print out the best SNP information */
#ifdef USE_MPI
  if (mpi_my_rank ==  mpi_master_rank) {
#endif

  if (print_each_marker == 1)
    for (i=0; i<number_snps; i++) {
      if (snp_info[i]->best_marker_r2 > -1.0) {
	printf ("%d: The best marker for rs%d is rs%d with r2 = %f\n", i, snp_info[i]->rs, snp_info[snp_info[i]->best_marker]->rs, snp_info[i]->best_marker_r2);
      } else {
	printf ("%d: There is no marker for rs%d\n", i, snp_info[i]->rs);
      }
    }

  /* Estimate power for different N and GRR */
  printf ("Fraction cases = %f\n", frac_case);
  printf ("Using %d total snps\n", total_snps);
  printf ("N\tGRR\tPower (1 df)\tPower (2df)\n");

  /* Calculate chisq value for the appropriate percentile based on the
     value of ALPHA, if set to -1.  Otherwise, use what was on the command
     line. */
  if (q_chisq1 < 0) {
    q_chisq1 = qchisq(1.0-(ALPHA/((double)(num_markers))),1.,TRUE,0);
  }
  if (q_chisq2 < 0) {
    q_chisq2 = qchisq(1.0-(ALPHA/((double)(num_markers))),2.,TRUE,0);
  }

  for (N = N_start; N <= N_stop; N +=N_step) {
    for (GRR_index = 0; GRR_index < num_GRRs; GRR_index++) {
      printf ("%d\t%f\t%f\t%f\n", N, GRRs[GRR_index], power(frac_case*N, (1.-frac_case)*N, GRRs[GRR_index], 1, q_chisq1, snp_info, number_snps, total_snps, minr2, grr_type, prevalence), power(frac_case*N, (1.-frac_case)*N, GRRs[GRR_index], 2, q_chisq2, snp_info, number_snps, total_snps, minr2, grr_type, prevalence));
    fflush(stdout);
    }
  }
#ifdef USE_MPI
  }

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
  in_mpi = 0;
#endif

  return EXIT_SUCCESS;
}




