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

char *vars = NULL;

char *matrix_fname = NULL;
char *inout_fname = NULL;

gridinfo_t grid;

double *nzval_col_wise;
int_t *rowind;
int_t *colptr;

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
         vars = optarg;
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
 *    retrieve row-oriented matrix entries from disk
 *    generate column-oriented matrix entries
 *    free row-oriented entries
 *    broadcast column-oriented entries
 * process >0:
 *    receive column-oriented entries
 */

int
get_sparse_matrix_global (void)
{
   if (iam == 0) {
      if (get_sparse_matrix (matrix_fname))
         return 1;
      if (dbg_lvl)
         printf ("row-oriented matrix read in\n");
      dCompRow_to_CompCol_dist (flat_len, flat_len, nnz, nzval_row_wise, colind, rowptr,
                                &nzval_col_wise, &rowind, &colptr);
      if (dbg_lvl)
         printf ("row-oriented matrix entries converted to column-oriented matrix entries\n");
      free_sparse_matrix ();
   }

   MPI_Bcast (&flat_len, 1, MPI_INT, 0, grid.comm);
   MPI_Bcast (&nnz, 1, MPI_INT, 0, grid.comm);
   if (iam > 0)
      dallocateA_dist (flat_len, nnz, &nzval_col_wise, &rowind, &colptr);
   MPI_Bcast (nzval_col_wise, nnz, MPI_DOUBLE, 0, grid.comm);
   MPI_Bcast (rowind, nnz, mpi_int_t, 0, grid.comm);
   MPI_Bcast (colptr, flat_len + 1, mpi_int_t, 0, grid.comm);
   if (dbg_lvl && (iam == 0))
      printf ("column-oriented matrix entries broadcast\n");

   return 0;
}

/******************************************************************************/

int
get_B_global (char *var, double *B)
{
   char *subname = "get_B_global";

   if (iam == 0) {
      double ***field_3d;
      int flat_ind;

      /* allocate needed variables */
      if ((field_3d = malloc_3d_double (km, jmt, imt)) == NULL) {
         fprintf (stderr, "malloc failed in %s for field_3d\n", subname);
         return 1;
      }

      if (get_var_3d_double (inout_fname, var, field_3d))
         return 1;

      for (flat_ind = 0; flat_ind < flat_len; flat_ind++) {
         int i = flat_ind_to_int3[flat_ind].i;
         int j = flat_ind_to_int3[flat_ind].j;
         int k = flat_ind_to_int3[flat_ind].k;
         B[flat_ind] = field_3d[k][j][i];
      }

      MPI_Bcast (B, flat_len, MPI_DOUBLE, 0, grid.comm);

      /* free allocated memory */
      free_3d_double (field_3d);
   } else {
      MPI_Bcast (B, flat_len, MPI_DOUBLE, 0, grid.comm);
   }

   return 0;
}

/******************************************************************************/

int
put_B_global (char *var, double *B)
{
   char *subname = "put_B_global";

   if (iam == 0) {
      double ***field_3d;
      int flat_ind;

      /* allocate needed variables */
      if ((field_3d = malloc_3d_double (km, jmt, imt)) == NULL) {
         fprintf (stderr, "malloc failed in %s for field_3d\n", subname);
         return 1;
      }

      /* read field in as 3D variable so that non-processed field values are preserved */
      if (get_var_3d_double (inout_fname, var, field_3d))
         return 1;

      for (flat_ind = 0; flat_ind < flat_len; flat_ind++) {
         int i = flat_ind_to_int3[flat_ind].i;
         int j = flat_ind_to_int3[flat_ind].j;
         int k = flat_ind_to_int3[flat_ind].k;
         field_3d[k][j][i] = B[flat_ind];
      }

      if (put_var_3d_double (inout_fname, var, field_3d))
         return 1;

      /* free allocated memory */
      free_3d_double (field_3d);
   }

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
   SuperLUStat_t stat;
   int info;

   int nrhs;
   int nrhs_tmp;

   double *B;
   double *berr;

   char *varsep = ",";
   char *vars_writable = NULL;
   char *var = NULL;

   /* initialize MPI */
   MPI_Init (&argc, &argv);

   /* set defaults */
   dbg_lvl = 0;
   nprow = npcol = 4;

   if (parse_cmd_line (argc, argv))
      exit (EXIT_FAILURE);

   /* initialize SuperLU process grid */
   superlu_gridinit (MPI_COMM_WORLD, nprow, npcol, &grid);
   iam = grid.iam;

   if (dbg_lvl && (iam == 0)) {
      printf ("dbg_lvl            = %d\n", dbg_lvl);
      printf ("nprow              = %lld\n", (long long) nprow);
      printf ("npcol              = %lld\n", (long long) npcol);
      printf ("vars               = %s\n", vars);
      printf ("matrix_fname       = %s\n", matrix_fname);
      printf ("inout_fname        = %s\n\n", inout_fname);
   }

   if (get_sparse_matrix_global ())
      exit (EXIT_FAILURE);

   /* create compressed column matrix for A */
   dCreate_CompCol_Matrix_dist (&A, flat_len, flat_len, nnz, nzval_col_wise, rowind, colptr,
                                SLU_NC, SLU_D, SLU_GE);
   if (dbg_lvl && (iam == 0))
      printf ("column-oriented matrix created\n");

   /* use default SuperLU options */
   set_default_options_dist (&options);
   if (dbg_lvl && (iam == 0))
      print_options_dist (&options);

   /* initialize ScalePermstruct, LUstruct */
   ScalePermstructInit (flat_len, flat_len, &ScalePermstruct);
   LUstructInit (flat_len, flat_len, &LUstruct);

   /* allocate space for RHS and solution */
   nrhs = 1;
   if ((B = doubleMalloc_dist (flat_len * nrhs)) == NULL)
      ABORT ("Malloc fails for B[].");
   if ((berr = doubleMalloc_dist (nrhs)) == NULL)
      ABORT ("Malloc fails for berr[].");

   /* factor matrix with no RHS */
   nrhs_tmp = 0;
   PStatInit (&stat);
   printf ("calling pdgssvx_ABglobal\n");
   pdgssvx_ABglobal (&options, &A, &ScalePermstruct, B, flat_len, nrhs_tmp, &grid,
                     &LUstruct, berr, &stat, &info);
   if (dbg_lvl) {
      if (iam == 0)
         printf ("dgssvx info = %d\n", info);
      if (options.PrintStat)
         PStatPrint (&options, &stat, &grid);
   }
   PStatFree (&stat);

   /* solve, each tracer is a different RHS */
   options.Fact = FACTORED;
   if (iam == 0) {
      if (get_ind_maps (matrix_fname))
         exit (EXIT_FAILURE);
      if (get_grid_dims (matrix_fname))
         exit (EXIT_FAILURE);
   }
   if ((vars_writable = malloc (strlen (vars) + 1)) == NULL) {
      fprintf (stderr, "malloc failed in main for vars_writable\n");
      exit (EXIT_FAILURE);
   }
   strcpy (vars_writable, vars);
   for (var = strtok (vars_writable, varsep); var; var = strtok (NULL, varsep)) {
      if (dbg_lvl && (iam == 0))
         printf ("processing variable %s\n", var);

      if (get_B_global (var, B))
         exit (EXIT_FAILURE);

      PStatInit (&stat);
      printf ("calling pdgssvx_ABglobal\n");
      pdgssvx_ABglobal (&options, &A, &ScalePermstruct, B, flat_len, nrhs, &grid,
                        &LUstruct, berr, &stat, &info);
      if (dbg_lvl) {
         if (iam == 0)
            printf ("dgssvx info = %d\n", info);
         if (options.PrintStat)
            PStatPrint (&options, &stat, &grid);
      }
      PStatFree (&stat);

      if (put_B_global (var, B))
         exit (EXIT_FAILURE);
   }

   /* deallocate storage */
   Destroy_CompCol_Matrix_dist (&A);
   Destroy_LU (flat_len, &grid, &LUstruct);
   ScalePermstructFree (&ScalePermstruct);
   LUstructFree (&LUstruct);
   if (iam == 0)
      free_ind_maps ();
   free (vars_writable);
   SUPERLU_FREE (B);
   SUPERLU_FREE (berr);

   /* release SuperLU process grid */
   superlu_gridexit (&grid);

   MPI_Finalize ();

   exit (EXIT_SUCCESS);
}
