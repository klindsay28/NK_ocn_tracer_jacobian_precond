#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>
#include <string.h>

#include "superlu_ddefs.h"
#include "globals.h"

#include "file_io.h"
#include "grid.h"
#include "matrix.h"
#include "memory.h"

/******************************************************************************/

int dbg_lvl;
int iam;

int_t nprow;
int_t npcol;
int_t nprocs;

char *vars = NULL;

char *matrix_fname = NULL;
char *inout_fname = NULL;

gridinfo_t grid;

int flat_len_loc;
int nnz_loc;
int fst_row;

double *nzval_loc;
int_t *colind_loc;
int_t *rowptr_loc;

/******************************************************************************/

int
parse_cmd_line (int argc, char **argv)
{
   char *subname = "parse_cmd_line";
   char *usage_msg =
      "usage: jacobian_precond [-D dbg_lvl] [-n nprow[,npcol]] [-v vars] matrix_fname inout_fname";
   extern char *optarg;
   extern int optind;
   char *optstring = "D:n:v:h";
   int opt;
   char *cp;

   while ((opt = getopt (argc, argv, optstring)) != -1) {
      switch (opt) {
      case 'D':
         dbg_lvl = atoi (optarg);
         break;
      case 'n':
         cp = strtok (optarg, ",");
         nprow = atoi (cp);
         if ((cp = strtok (NULL, ",")) == NULL)
            npcol = nprow;
         else
            npcol = atoi (cp);
         break;
      case 'v':
         /* copy optarg string, so that later use of strtok doesn't modify the argv array */
         /* this is possibly not necessary */
         if ((vars = malloc (strlen (optarg) + 1)) == NULL) {
            fprintf (stderr, "malloc failed in parse_cmd_line for vars\n");
            return 1;
         }
         strcpy (vars, optarg);
         break;
      case 'h':
      case '?':
         fprintf (stderr, "%s\n", usage_msg);
         return 1;
      default:
         fprintf (stderr, "internal error: unhandled option '-%c'\n", opt);
         return 1;
      }
   }
   if (optind != argc - 2) {
      fprintf (stderr, "unexpected number of arguments\n%s\n", usage_msg);
      return 1;
   }
   matrix_fname = argv[optind++];
   inout_fname = argv[optind++];
   return 0;
}

/******************************************************************************/
/* process 0:
 *    retrieve global row-oriented matrix entries from disk
 *    sent row-oriented segments
 *    free global row-oriented entries
 * process >0:
 *    receive row-oriented segments
 */

int
get_sparse_matrix_dist (void)
{
   MPI_Status status;
   int flat_ind;
   int nnz_ind;
   int dest;

   if (iam == 0) {
      if (get_sparse_matrix (matrix_fname))
         return 1;
      if (dbg_lvl)
         printf ("row-oriented matrix read in\n");
   }

   MPI_Bcast (&flat_len, 1, MPI_INT, 0, grid.comm);

   /* compute flat_len_loc & fst_row */
   /* adjust flat_len_loc for last task, to get remaining rows */
   flat_len_loc = flat_len / nprocs;
   fst_row = iam * flat_len_loc;
   if (iam == nprocs - 1)
      flat_len_loc = flat_len - fst_row;

   if ((rowptr_loc = intMalloc_dist (flat_len_loc + 1)) == NULL)
      ABORT ("Malloc fails for rowptr_loc.");

   if (iam == 0) {
      /* copy master task rowptr_loc entries from rowptr */
      for (flat_ind = 0; flat_ind <= flat_len_loc; flat_ind++)
         rowptr_loc[flat_ind] = rowptr[flat_ind];

      /* send global rowptr entries to non-master tasks */
      for (dest = 1; dest < nprocs; dest++) {
         int count;
         int fst_row_dest;

         count = flat_len / nprocs;
         fst_row_dest = dest * count;
         if (dest == nprocs - 1)
            count = flat_len - fst_row_dest;
         MPI_Send (&rowptr[fst_row_dest], count + 1, mpi_int_t, dest, 0, grid.comm);
      }

      if (dbg_lvl)
         printf ("rowptr segments sent\n");
   } else {
      /* receive global rowptr_loc entries from master task */
      MPI_Recv (rowptr_loc, flat_len_loc + 1, mpi_int_t, 0, 0, grid.comm, &status);

      /* localize global rowptr_loc entries */
      for (flat_ind = 1; flat_ind <= flat_len_loc; flat_ind++)
         rowptr_loc[flat_ind] -= rowptr_loc[0];
      rowptr_loc[0] = 0;
   }

   nnz_loc = rowptr_loc[flat_len_loc];

   if (dbg_lvl > 1)
      printf ("fst_row, flat_len_loc, nnz_loc = %d, %d, %d\n", fst_row, flat_len_loc, nnz_loc);

   if ((colind_loc = intMalloc_dist (nnz_loc)) == NULL)
      ABORT ("Malloc fails for colind_loc.");
   if ((nzval_loc = doubleMalloc_dist (nnz_loc)) == NULL)
      ABORT ("Malloc fails for nzval_loc.");

   if (iam == 0) {
      int count_cumsum;

      /* copy master task colind_loc and nzval_loc entries entries from colind and nzval */
      for (nnz_ind = 0; nnz_ind < nnz_loc; nnz_ind++) {
         colind_loc[nnz_ind] = colind[nnz_ind];
         nzval_loc[nnz_ind] = nzval_row_wise[nnz_ind];
      }

      count_cumsum = nnz_loc;

      /* send global colind and nzval entries to non-master tasks */
      for (dest = 1; dest < nprocs; dest++) {
         int count;
         int fst_row_dest;

         count = flat_len / nprocs;
         fst_row_dest = dest * count;
         if (dest == nprocs - 1)
            count = flat_len - fst_row_dest;

         count = rowptr[fst_row_dest + count] - rowptr[fst_row_dest];

         MPI_Send (&colind[count_cumsum], count, mpi_int_t, dest, 1, grid.comm);
         MPI_Send (&nzval_row_wise[count_cumsum], count, MPI_DOUBLE, dest, 2, grid.comm);

         count_cumsum += count;
      }

      if (dbg_lvl)
         printf ("colind and nzval segments sent\n");

      free_sparse_matrix ();
   } else {
      /* receive global colind_loc and nzval_loc entries from master task */
      MPI_Recv (colind_loc, nnz_loc, mpi_int_t, 0, 1, grid.comm, &status);
      MPI_Recv (nzval_loc, nnz_loc, MPI_DOUBLE, 0, 2, grid.comm, &status);
   }

   if (dbg_lvl > 1) {
      for (flat_ind = 0; flat_ind < flat_len_loc; flat_ind++) {
         for (nnz_ind = rowptr_loc[flat_ind]; nnz_ind < rowptr_loc[flat_ind + 1]; nnz_ind++) {
            printf (" A : %d %lld %e\n", fst_row + flat_ind, (long long) colind_loc[nnz_ind],
                    nzval_loc[nnz_ind]);
         }
      }
   }

   fflush (stdout);
   MPI_Barrier (grid.comm);

   return 0;
}

/******************************************************************************/

int
get_B_dist (char *var, double *B)
{
   char *subname = "get_B_dist";

   if (iam == 0) {
      double ***field_3d;
      double *B_global;
      int flat_ind;
      int dest;

      /* allocate needed variables */
      if ((field_3d = malloc_3d_double (km, jmt, imt)) == NULL) {
         fprintf (stderr, "malloc failed in %s for field_3d\n", subname);
         return 1;
      }
      if ((B_global = malloc ((size_t) flat_len * sizeof (double))) == NULL) {
         fprintf (stderr, "malloc failed in %s for B_global\n", subname);
         return 1;
      }

      /* read field in as 3D variable */
      if (get_var_3d_double (inout_fname, var, field_3d))
         return 1;

      /* flatten 3D variable into B_global */
      for (flat_ind = 0; flat_ind < flat_len; flat_ind++) {
         int i = flat_ind_to_int3[flat_ind].i;
         int j = flat_ind_to_int3[flat_ind].j;
         int k = flat_ind_to_int3[flat_ind].k;
         B_global[flat_ind] = field_3d[k][j][i];
      }

      /* copy master task B entries from B_global */
      for (flat_ind = 0; flat_ind <= flat_len_loc; flat_ind++)
         B[flat_ind] = B_global[flat_ind];

      /* send B_global entry segments to non-master tasks */
      for (dest = 1; dest < nprocs; dest++) {
         int count;
         int fst_row_dest;

         count = flat_len / nprocs;
         fst_row_dest = dest * count;
         if (dest == nprocs - 1)
            count = flat_len - fst_row_dest;
         MPI_Send (&B_global[fst_row_dest], count, MPI_DOUBLE, dest, 3, grid.comm);
      }

      /* free allocated memory */
      free (B_global);
      free_3d_double (field_3d);
   } else {
      MPI_Status status;
      MPI_Recv (B, flat_len_loc, MPI_DOUBLE, 0, 3, grid.comm, &status);
   }

   fflush (stdout);
   MPI_Barrier (grid.comm);

   return 0;
}

/******************************************************************************/

int
put_B_dist (char *var, double *B)
{
   char *subname = "put_B_dist";

   if (iam == 0) {
      MPI_Status status;
      double ***field_3d;
      double *B_global;
      int flat_ind;
      int source;

      /* allocate needed variables */
      if ((field_3d = malloc_3d_double (km, jmt, imt)) == NULL) {
         fprintf (stderr, "malloc failed in %s for field_3d\n", subname);
         return 1;
      }
      if ((B_global = malloc ((size_t) flat_len * sizeof (double))) == NULL) {
         fprintf (stderr, "malloc failed in %s for B_global\n", subname);
         return 1;
      }

      /* copy master task B entries from B_global */
      for (flat_ind = 0; flat_ind <= flat_len_loc; flat_ind++)
         B_global[flat_ind] = B[flat_ind];

      /* receive B_global entry segments from non-master tasks */
      for (source = 1; source < nprocs; source++) {
         int count;
         int fst_row_source;

         count = flat_len / nprocs;
         fst_row_source = source * count;
         if (source == nprocs - 1)
            count = flat_len - fst_row_source;
         MPI_Recv (&B_global[fst_row_source], count, MPI_DOUBLE, source, 4, grid.comm, &status);
      }

      /* read field in as 3D variable so that non-processed field values are preserved */
      if (get_var_3d_double (inout_fname, var, field_3d))
         return 1;

      /* convert flat B_global into 3D variable */
      for (flat_ind = 0; flat_ind < flat_len; flat_ind++) {
         int i = flat_ind_to_int3[flat_ind].i;
         int j = flat_ind_to_int3[flat_ind].j;
         int k = flat_ind_to_int3[flat_ind].k;
         field_3d[k][j][i] = B_global[flat_ind];
      }

      /* write out 3D variable */
      if (put_var_3d_double (inout_fname, var, field_3d))
         return 1;

      /* free allocated memory */
      free (B_global);
      free_3d_double (field_3d);
   } else {
      MPI_Send (B, flat_len_loc, MPI_DOUBLE, 0, 4, grid.comm);
   }

   fflush (stdout);
   MPI_Barrier (grid.comm);

   return 0;
}

/******************************************************************************/

int
main (int argc, char *argv[])
{
   SuperMatrix A;
   superlu_options_t options;
   ScalePermstruct_t ScalePermstruct;
   LUstruct_t LUstruct;
   SOLVEstruct_t SOLVEstruct;
   SuperLUStat_t stat;
   int info;

   int nrhs;
   int nrhs_tmp;

   double *B;
   double *berr;

   char *varsep = ",";
   char *var = NULL;

   /* initialize MPI */
   MPI_Init (&argc, &argv);

   /* set defaults */
   dbg_lvl = 0;
   nprow = npcol = 4;

   if (parse_cmd_line (argc, argv))
      exit (EXIT_FAILURE);

   nprocs = nprow * npcol;

   /* initialize SuperLU process grid */
   superlu_gridinit (MPI_COMM_WORLD, nprow, npcol, &grid);
   iam = grid.iam;

   if (dbg_lvl && (iam == 0)) {
      printf ("dbg_lvl            = %d\n", dbg_lvl);
      printf ("nprocs             = %lld\n", (long long) nprocs);
      printf ("nprow              = %lld\n", (long long) nprow);
      printf ("npcol              = %lld\n", (long long) npcol);
      printf ("vars               = %s\n", vars);
      printf ("matrix_fname       = %s\n", matrix_fname);
      printf ("inout_fname        = %s\n\n", inout_fname);
   }

   if (get_sparse_matrix_dist ())
      exit (EXIT_FAILURE);

   /* create distributed compressed row matrix for A */
   dCreate_CompRowLoc_Matrix_dist (&A, flat_len, flat_len, nnz_loc, flat_len_loc, fst_row,
                                   nzval_loc, colind_loc, rowptr_loc, SLU_NR_loc, SLU_D,
                                   SLU_GE);
   if (dbg_lvl && (iam == 0))
      printf ("distributed row-oriented matrix created\n");

   /* use default SuperLU options */
   set_default_options_dist (&options);

   /* set any non-standard SuperLU options */
#if 0
   options.ColPerm = PARMETIS;
#endif

   if (dbg_lvl && (iam == 0))
      print_options_dist (&options);

   /* initialize ScalePermstruct, LUstruct */
   ScalePermstructInit (flat_len, flat_len, &ScalePermstruct);
   LUstructInit (flat_len, flat_len, &LUstruct);

   /* allocate space for RHS and solution */
   nrhs = 1;
   if ((B = doubleMalloc_dist (flat_len_loc * nrhs)) == NULL)
      ABORT ("Malloc fails for B[].");
   if ((berr = doubleMalloc_dist (nrhs)) == NULL)
      ABORT ("Malloc fails for berr[].");

   fflush (stdout);
   MPI_Barrier (grid.comm);

   /* factor matrix with no RHS */
   nrhs_tmp = 0;
   PStatInit (&stat);
   printf ("calling pdgssvx\n");
   pdgssvx (&options, &A, &ScalePermstruct, B, flat_len_loc, nrhs_tmp, &grid,
            &LUstruct, &SOLVEstruct, berr, &stat, &info);
   if (dbg_lvl) {
      if (iam == 0)
         printf ("dgssvx info = %d\n", info);
      if (options.PrintStat)
         PStatPrint (&options, &stat, &grid);
   }
   PStatFree (&stat);

#if 0
   if (iam == 0) {
      int flat_ind;
      for (flat_ind = 0; flat_ind < flat_len; flat_ind++)
         printf (" perm_c[%d] = %lld\n", flat_ind,
                 (long long) (ScalePermstruct.perm_c[flat_ind]));
   }
#endif

   fflush (stdout);
   MPI_Barrier (grid.comm);

   /* solve, each tracer is a different RHS */
   options.Fact = FACTORED;
   if (iam == 0) {
      if (get_ind_maps (matrix_fname))
         exit (EXIT_FAILURE);
      if (get_grid_dims (matrix_fname))
         exit (EXIT_FAILURE);
   }
   for (var = strtok (vars, varsep); var; var = strtok (NULL, varsep)) {
      if (dbg_lvl && (iam == 0))
         printf ("processing variable %s\n", var);

      if (get_B_dist (var, B))
         exit (EXIT_FAILURE);

      PStatInit (&stat);
      printf ("calling pdgssvx\n");
      pdgssvx (&options, &A, &ScalePermstruct, B, flat_len_loc, nrhs, &grid,
               &LUstruct, &SOLVEstruct, berr, &stat, &info);
      if (dbg_lvl) {
         if (iam == 0)
            printf ("dgssvx info = %d\n", info);
         if (options.PrintStat)
            PStatPrint (&options, &stat, &grid);
      }
      PStatFree (&stat);

      fflush (stdout);
      MPI_Barrier (grid.comm);

      if (put_B_dist (var, B))
         exit (EXIT_FAILURE);
   }

   /* deallocate storage */
   Destroy_CompRowLoc_Matrix_dist (&A);
   ScalePermstructFree (&ScalePermstruct);
   Destroy_LU (flat_len_loc, &grid, &LUstruct);
   LUstructFree (&LUstruct);
   if (options.SolveInitialized)
      dSolveFinalize (&options, &SOLVEstruct);
   if (iam == 0)
      free_ind_maps ();
   free (vars);
   SUPERLU_FREE (B);
   SUPERLU_FREE (berr);

   /* release SuperLU process grid */
   superlu_gridexit (&grid);

   MPI_Finalize ();

   exit (EXIT_SUCCESS);
}
