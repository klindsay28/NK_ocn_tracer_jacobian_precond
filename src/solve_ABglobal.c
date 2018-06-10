#include <stdlib.h>
#include <unistd.h>
#include <getopt.h>
#include <stdio.h>
#include <string.h>

#include "superlu_ddefs.h"
#include "globals.h"

#include "file_io.h"
#include "grid.h"
#include "matrix.h"
#include "memory.h"
#include "misc.h"

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
   char *usage_msg = "usage: jacobian_precond [-D dbg_lvl] [-n nprow[,npcol]] [-v vars] matrix_fname inout_fname";
   extern char *optarg;
   extern int optind;
   char *optstring = "D:n:v:h";
   int opt;
   char *cp;
   long lval;

   while ((opt = getopt (argc, argv, optstring)) != -1) {
      switch (opt) {
      case '?':
      case 'h':
         fprintf (stderr, "(%d) %s\n", iam, usage_msg);
         return 1;
      case 'D':
         if (parse_to_int (optarg, &dbg_lvl)) {
            fprintf (stderr, "(%d) error parsing argument '%s' for option '%c'\n", iam, optarg, opt);
            return 1;
         }
         break;
      case 'n':
         cp = strtok (optarg, ",");
         if (parse_to_long (cp, &lval)) {
            fprintf (stderr, "(%d) error parsing argument '%s' for option '%c'\n", iam, cp, opt);
            return 1;
         }
         nprow = lval;
         if ((cp = strtok (NULL, ",")) == NULL) {
            npcol = nprow;
         } else {
            if (parse_to_long (cp, &lval)) {
               fprintf (stderr, "(%d) error parsing argument '%s' for option '%c'\n", iam, cp, opt);
               return 1;
            }
            npcol = lval;
         }
         break;
      case 'v':
         /* copy optarg string, so that later use of strtok doesn't modify the argv array */
         /* this is possibly not necessary */
         if ((vars = malloc (strlen (optarg) + 1)) == NULL) {
            fprintf (stderr, "(%d) malloc failed in parse_cmd_line for vars\n", iam);
            return 1;
         }
         strcpy (vars, optarg);
         break;
      default:
         fprintf (stderr, "(%d) internal error: unhandled option '-%c'\n", iam, opt);
         return 1;
      }
   }
   if (optind != argc - 2) {
      fprintf (stderr, "(%d) unexpected number of arguments\n%s\n", iam, usage_msg);
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
   char *subname = "get_sparse_matrix_global";

   if (dbg_lvl > 1) {
      printf ("(%d) entering %s\n", iam, subname);
      fflush (stdout);
   }

   if (iam == 0) {
      if (get_sparse_matrix (matrix_fname))
         return 1;
      if (dbg_lvl)
         printf ("(%d) row-oriented matrix read in\n", iam);
      dCompRow_to_CompCol_dist (flat_len, flat_len, nnz, nzval_row_wise, colind, rowptr, &nzval_col_wise, &rowind, &colptr);
      if (dbg_lvl)
         printf ("(%d) row-oriented matrix entries converted to column-oriented matrix entries\n", iam);
      free_sparse_matrix ();
   }

   MPI_Bcast (&coupled_tracer_cnt, 1, MPI_INT, 0, grid.comm);
   MPI_Bcast (&flat_len, 1, MPI_INT, 0, grid.comm);
   MPI_Bcast (&nnz, 1, MPI_INT, 0, grid.comm);
   if (iam > 0)
      dallocateA_dist (flat_len, nnz, &nzval_col_wise, &rowind, &colptr);
   MPI_Bcast (nzval_col_wise, nnz, MPI_DOUBLE, 0, grid.comm);
   MPI_Bcast (rowind, nnz, mpi_int_t, 0, grid.comm);
   MPI_Bcast (colptr, flat_len + 1, mpi_int_t, 0, grid.comm);
   if (dbg_lvl && (iam == 0))
      printf ("(%d) column-oriented matrix entries broadcast\n", iam);

   if (dbg_lvl > 1) {
      printf ("(%d) exiting %s\n", iam, subname);
      fflush (stdout);
   }

   return 0;
}

/******************************************************************************/

int
get_B_global (char **vars_per_solve, double *B)
{
   char *subname = "get_B_global";

   if (dbg_lvl > 1) {
      printf ("(%d) entering %s\n", iam, subname);
      fflush (stdout);
   }

   if (iam == 0) {
      double ***field_3d;
      int tracer_ind;
      int flat_ind_offset;
      int tracer_state_ind;
      int flat_ind;

      /* allocate needed variables */
      if ((field_3d = malloc_3d_double (km, jmt, imt)) == NULL) {
         fprintf (stderr, "(%d) malloc failed in %s for field_3d\n", iam, subname);
         return 1;
      }

      for (tracer_ind = 0; tracer_ind < coupled_tracer_cnt; tracer_ind++) {
         /* read single tracer in as 3D variable */
         if (dbg_lvl)
            printf ("(%d) reading %s from %s\n", iam, vars_per_solve[tracer_ind], inout_fname);
         if (get_var_3d_double (inout_fname, vars_per_solve[tracer_ind], field_3d))
            return 1;

         /* flatten 3D variable into B */
         flat_ind_offset = tracer_ind * tracer_state_len;
         for (tracer_state_ind = 0; tracer_state_ind < tracer_state_len; tracer_state_ind++) {
            int i = tracer_state_ind_to_int3[tracer_state_ind].i;
            int j = tracer_state_ind_to_int3[tracer_state_ind].j;
            int k = tracer_state_ind_to_int3[tracer_state_ind].k;
            flat_ind = flat_ind_offset + tracer_state_ind;
            B[flat_ind] = field_3d[k][j][i];
         }
      }

      MPI_Bcast (B, flat_len, MPI_DOUBLE, 0, grid.comm);

      /* free allocated memory */
      free_3d_double (field_3d);
   } else {
      MPI_Bcast (B, flat_len, MPI_DOUBLE, 0, grid.comm);
   }

   if (dbg_lvl > 1) {
      printf ("(%d) exiting %s\n", iam, subname);
      fflush (stdout);
   }

   return 0;
}

/******************************************************************************/

int
put_B_global (char **vars_per_solve, double *B)
{
   char *subname = "put_B_global";

   if (dbg_lvl > 1) {
      printf ("(%d) entering %s\n", iam, subname);
      fflush (stdout);
   }

   if (iam == 0) {
      double ***field_3d;
      int tracer_ind;
      int flat_ind_offset;
      int tracer_state_ind;
      int flat_ind;

      /* allocate needed variables */
      if ((field_3d = malloc_3d_double (km, jmt, imt)) == NULL) {
         fprintf (stderr, "(%d) malloc failed in %s for field_3d\n", iam, subname);
         return 1;
      }

      for (tracer_ind = 0; tracer_ind < coupled_tracer_cnt; tracer_ind++) {
         /* read field in as 3D variable so that non-processed field values are preserved */
         if (get_var_3d_double (inout_fname, vars_per_solve[tracer_ind], field_3d))
            return 1;

         /* convert flat B in to 3D variable B */
         flat_ind_offset = tracer_ind * tracer_state_len;
         for (tracer_state_ind = 0; tracer_state_ind < tracer_state_len; tracer_state_ind++) {
            int i = tracer_state_ind_to_int3[tracer_state_ind].i;
            int j = tracer_state_ind_to_int3[tracer_state_ind].j;
            int k = tracer_state_ind_to_int3[tracer_state_ind].k;
            flat_ind = flat_ind_offset + tracer_state_ind;
            field_3d[k][j][i] = B[flat_ind];
         }

         /* write out 3D variable */
         if (dbg_lvl)
            printf ("(%d) writing %s to %s\n", iam, vars_per_solve[tracer_ind], inout_fname);
         if (put_var_3d_double (inout_fname, vars_per_solve[tracer_ind], field_3d))
            return 1;
      }

      /* free allocated memory */
      free_3d_double (field_3d);
   }

   if (dbg_lvl > 1) {
      printf ("(%d) exiting %s\n", iam, subname);
      fflush (stdout);
   }

   return 0;
}

/******************************************************************************/

int
main (int argc, char *argv[])
{
   SuperMatrix A;
   superlu_dist_options_t options;
   ScalePermstruct_t ScalePermstruct;
   LUstruct_t LUstruct;
   SuperLUStat_t stat;
   int info;

   int nrhs;
   int nrhs_tmp;

   double *B;
   double *berr;

   char **vars_per_solve = NULL;
   char *varsep = ",";
   char *var = NULL;

   /* initialize MPI */
   MPI_Init (&argc, &argv);

   /* set defaults */
   dbg_lvl = 0;
   nprow = npcol = 4;

   MPI_Comm_rank (MPI_COMM_WORLD, &iam);

   if (parse_cmd_line (argc, argv))
      exit (EXIT_FAILURE);

   /* only run remainder of program for nprow * npcol tasks */
   if (iam < nprow * npcol) {

      /* initialize SuperLU process grid */
      superlu_gridinit (MPI_COMM_WORLD, nprow, npcol, &grid);

      if (dbg_lvl && (iam == 0)) {
         printf ("(%d) dbg_lvl            = %d\n", iam, dbg_lvl);
         printf ("(%d) nprow              = %lld\n", iam, (long long) nprow);
         printf ("(%d) npcol              = %lld\n", iam, (long long) npcol);
         printf ("(%d) vars               = %s\n", iam, vars);
         printf ("(%d) matrix_fname       = %s\n", iam, matrix_fname);
         printf ("(%d) inout_fname        = %s\n\n", iam, inout_fname);
      }

      if (get_sparse_matrix_global ())
         exit (EXIT_FAILURE);

      if ((vars_per_solve = malloc ((size_t) coupled_tracer_cnt * sizeof (char *))) == NULL) {
         fprintf (stderr, "(%d) malloc failed in %s for vars_per_solve\n", iam, argv[0]);
         exit (EXIT_FAILURE);
      }

      /* create compressed column matrix for A */
      dCreate_CompCol_Matrix_dist (&A, flat_len, flat_len, nnz, nzval_col_wise, rowind, colptr, SLU_NC, SLU_D, SLU_GE);
      if (dbg_lvl && (iam == 0))
         printf ("(%d) column-oriented matrix created\n", iam);

      /* use default SuperLU options */
      set_default_options_dist (&options);
      options.ParSymbFact = YES;
      options.ColPerm = PARMETIS;
      if (dbg_lvl && (iam == 0))
         print_options_dist (&options);

      /* initialize ScalePermstruct, LUstruct */
      ScalePermstructInit (flat_len, flat_len, &ScalePermstruct);
      LUstructInit (flat_len, &LUstruct);

      /* allocate space for RHS and solution */
      nrhs = 1;
      if ((B = doubleMalloc_dist (flat_len * nrhs)) == NULL)
         ABORT ("Malloc fails for B[].");
      if ((berr = doubleMalloc_dist (nrhs)) == NULL)
         ABORT ("Malloc fails for berr[].");

      /* factor matrix with no RHS */
      nrhs_tmp = 0;
      PStatInit (&stat);
      printf ("(%d) calling pdgssvx_ABglobal\n", iam);
      pdgssvx_ABglobal (&options, &A, &ScalePermstruct, B, flat_len, nrhs_tmp, &grid, &LUstruct, berr, &stat, &info);
      if (dbg_lvl) {
         if (iam == 0)
            printf ("(%d) dgssvx info = %d\n", iam, info);
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
      for (var = strtok (vars, varsep); var; var = strtok (NULL, varsep)) {
         int tracer_ind;

         for (tracer_ind = 0; tracer_ind < coupled_tracer_cnt; tracer_ind++) {
            if (tracer_ind > 0) {
               var = strtok (NULL, varsep);
               if (var == NULL) {
                  fprintf (stderr, "(%d) error extracting tracer_ind=%d, ran out of var names\n", iam, tracer_ind);
                  exit (EXIT_FAILURE);
               }
            }
            if (dbg_lvl && (iam == 0))
               printf ("(%d) processing variable %s\n", iam, var);
            if ((vars_per_solve[tracer_ind] = malloc (1 + strlen (var))) == NULL) {
               fprintf (stderr, "(%d) malloc failed in %s for vars_per_solve[%d]\n", iam, argv[0], tracer_ind);
               exit (EXIT_FAILURE);
            }
            strcpy (vars_per_solve[tracer_ind], var);
         }

         if (get_B_global (vars_per_solve, B))
            exit (EXIT_FAILURE);

         PStatInit (&stat);
         printf ("(%d) calling pdgssvx_ABglobal\n", iam);
         pdgssvx_ABglobal (&options, &A, &ScalePermstruct, B, flat_len, nrhs, &grid, &LUstruct, berr, &stat, &info);
         if (dbg_lvl) {
            if (iam == 0)
               printf ("(%d) dgssvx info = %d\n", iam, info);
            if (options.PrintStat)
               PStatPrint (&options, &stat, &grid);
         }
         PStatFree (&stat);

         if (put_B_global (vars_per_solve, B))
            exit (EXIT_FAILURE);

         for (tracer_ind = 0; tracer_ind < coupled_tracer_cnt; tracer_ind++)
            free (vars_per_solve[tracer_ind]);
      }

      /* deallocate storage */
      Destroy_CompCol_Matrix_dist (&A);
      Destroy_LU (flat_len, &grid, &LUstruct);
      ScalePermstructFree (&ScalePermstruct);
      LUstructFree (&LUstruct);
      if (iam == 0)
         free_ind_maps ();
      free (vars_per_solve);
      free (vars);
      SUPERLU_FREE (B);
      SUPERLU_FREE (berr);

      /* release SuperLU process grid */
      superlu_gridexit (&grid);

   }

   MPI_Finalize ();

   exit (EXIT_SUCCESS);
}
