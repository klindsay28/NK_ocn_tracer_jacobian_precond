/* functions to support handling sparse approximation to NK Jacobian
 *
 * primary functions for matrix handling
 *    gen_sparse_matrix (double day_cnt) : generate the sparse matrix
 *    put_sparse_matrix (char *fname) : write sparse matrix to a file
 *    get_sparse_matrix (char *fname) : read sparse matrix from a file
 *    free_sparse_matrix () : free memory associated with matrix
 *
 * The sparse matrix is stored in a compressed row storage format.
 * Related variables:
 *    coupled_tracer_cnt : number of tracers
 *    flat_len : size of matrix, i.e. number of tracers x tracer_state_len
 *    nnz : number of non-zero values in matrix
 *    nzval_row_wise : matrix values, stored row-by-row
 *    colind : column index for each matrix value in nzval_row_wise
 *    rowptr : index into nzval_row_wise and colind of starting location of entries for jth row
 *
 *    The following arrays are indices into nzval_row_wise and colind of specific categories of entries.
 *
 *    coef_ind_self : index of term for tracer cell value's dependence on itself
 *    coef_ind_adv_non_nbr : index of terms for non-neighbor advective dependencies
 *    coef_ind_hmix_non_nbr : index of terms for non-neighbor hmix dependencies
 *    coef_ind_vmix_non_nbr : index of terms for non-neighbor vmix dependencies
 *    coef_ind_sink_non_nbr : index of terms for non-neighbor sink dependencies
 *    coef_ind_sink_other_tracers : index of terms for dependencies on other tracers
 *
 * related functions
 *    comp_flat_len () : compute flat_len
 *    comp_nnz () : compute nnz
 *       this uses the functions adv_non_nbr_cnt, hmix_non_nbr_cnt, vmix_non_nbr_cnt, sink_non_nbr_cnt
 *    allocate_matrix_arrays() : allocate matrix variables
 *    init_matrix () : generate matrix index arrays, initialize nzval_row_wise to zero
 *
 *    add_adv () : add advection related terms to matrix
 *    add_hmix () : add lateral mixing related terms to matrix
 *    add_vmix () : add vertical mixing related terms to matrix
 *    add_sink_pure_diag () : add pure diagonal, i.e. non-generic, source-sink terms to matrix
 *    add_sink_generic_tracer () : add single generic tracer source-sink terms to matrix
 *    add_pv () : add piston velocity terms to matrix
 *    add_d_SF_d_TRACER () : add generic surface flux terms to matrix
 *
 *    sum_dup_vals () : combine matrix entries in a row that have the same column index
 *    strip_matrix_zeros () : remove zeros from matrix
 *    check_matrix_diag () : verify that matrix diagonal entries are non-zero
 *    sort_cols_all_rows () : sort entries in a row into increasing column order
 *
 * Access to GCM output fields is natural with index triplets, while access to matrix
 * is convenient with flat tracer indexing.
 *
 * variables related to mapping between GCM index triplets to flat tracer index space
 *    tracer_state_len : length of state vector for one tracer, i.e. number of active grid points in GCM
 *    int3_to_tracer_state_ind : mapping from GCM index triplets to flat tracer index space
 *    tracer_state_ind_to_int3 : mapping from flat tracer index space to GCM index triplets
 *
 * related functions
 *    comp_tracer_state_len () : compute tracer_state_len
 *    gen_ind_maps () : generate index mapping variables, int3_to_tracer_state_ind and tracer_state_ind_to_int3
 *    put_ind_maps (char *fname) : write index mapping variables to a file
 *    get_ind_maps (char *fname) : read index mapping variables from a file
 *    free_ind_maps () : free memory associated with mapping variables
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "globals.h"

#include "file_io.h"
#include "grid.h"
#include "memory.h"
#include "netcdf.h"
#include "superlu_ddefs.h"

#include "matrix.h"

int tracer_state_len;
int ***int3_to_tracer_state_ind;
int3 *tracer_state_ind_to_int3;

int coupled_tracer_cnt;

int flat_len;
int nnz;
double *nzval_row_wise;
int_t *colind;
int_t *rowptr;
int_t **coef_ind_self;
int_t **coef_ind_adv_non_nbr;
int_t **coef_ind_hmix_non_nbr;
int_t **coef_ind_vmix_non_nbr;
int_t **coef_ind_sink_non_nbr;
int_t **coef_ind_sink_other_tracers;

double delta_t;
double year_cnt;

adv_opt_t adv_opt;
int l_adv_enforce_divfree;

hmix_opt_t hmix_opt;

vmix_opt_t vmix_opt;

char *tracer_fname = NULL;

per_tracer_opt_t *per_tracer_opt = NULL;

coupled_tracer_opt_t coupled_tracer_opt;

char *OCMIP_BGC_PO4_DOP_names[] = { "OCMIP_BGC_PO4", "OCMIP_BGC_DOP" };
char *DIC_SHADOW_ALK_SHADOW_names[] = { "DIC_SHADOW", "ALK_SHADOW" };

/******************************************************************************/

void
set_3d_double (double val, double ***FIELD)
{
   int i;
   int j;
   int k;

   for (k = 0; k < km; k++)
      for (j = 0; j < jmt; j++)
         for (i = 0; i < imt; i++)
            FIELD[k][j][i] = val;
}

/******************************************************************************/

void
set_fv_2d_double (double fv, double val, double **FIELD)
{
   int i;
   int j;

   for (j = 0; j < jmt; j++)
      for (i = 0; i < imt; i++)
         if (FIELD[j][i] == fv)
            FIELD[j][i] = val;
}

/******************************************************************************/

void
set_fv_3d_double (double fv, double val, double ***FIELD)
{
   int i;
   int j;
   int k;

   for (k = 0; k < km; k++)
      for (j = 0; j < jmt; j++)
         for (i = 0; i < imt; i++)
            if (FIELD[k][j][i] == fv)
               FIELD[k][j][i] = val;
}

/******************************************************************************/

int
comp_tracer_state_len (void)
{
   char *subname = "comp_tracer_state_len";
   int south_flag;
   int north_flag;
   int i;
   int j;

   if (dbg_lvl > 1) {
      printf ("(%d) entering %s\n", iam, subname);
      fflush (stdout);
   }

   /* verify that KMT is 0 on southern- and northern-most rows */
   south_flag = north_flag = 0;
   for (i = 0; i < imt; i++) {
      if (KMT[0][i])
         south_flag = 1;
      if (KMT[jmt - 1][i])
         north_flag = 1;
   }
   if (south_flag)
      fprintf (stderr, "(%d) non-land found on southern-most row in %s\n", iam, subname);
   if (north_flag)
      fprintf (stderr, "(%d) non-land found on northern-most row in %s\n", iam, subname);
   if (south_flag || north_flag)
      return 1;

   tracer_state_len = 0;
   for (j = 0; j < jmt; j++)
      for (i = 0; i < imt; i++)
         tracer_state_len += KMT[j][i];
   if (dbg_lvl)
      printf ("(%d) tracer_state_len = %d\n\n", iam, tracer_state_len);

   if (dbg_lvl > 1) {
      printf ("(%d) exiting %s\n", iam, subname);
      fflush (stdout);
   }

   return 0;
}

/******************************************************************************/

/* generate mappings between GCM index triplets and flat tracer index space */

int
gen_ind_maps (void)
{
   char *subname = "gen_ind_maps";
   int i;
   int j;
   int k;
   int tracer_state_ind;

   if (dbg_lvl > 1) {
      printf ("(%d) entering %s\n", iam, subname);
      fflush (stdout);
   }

   if (comp_tracer_state_len ()) {
      fprintf (stderr, "(%d) comp_tracer_state_len call failed in %s\n", iam, subname);
      return 1;
   }
   if ((int3_to_tracer_state_ind = malloc_3d_int (km, jmt, imt)) == NULL) {
      fprintf (stderr, "(%d) malloc failed in %s for int3_to_tracer_state_ind\n", iam, subname);
      return 1;
   }
   if ((tracer_state_ind_to_int3 = malloc ((size_t) tracer_state_len * sizeof (int3))) == NULL) {
      fprintf (stderr, "(%d) malloc failed in %s for tracer_state_ind_to_int3\n", iam, subname);
      return 1;
   }
   tracer_state_ind = 0;
   if (dbg_lvl > 2)
      printf ("(%d) mappings between flat and 3d indices\n", iam);
   for (j = 0; j < jmt; j++)
      for (i = 0; i < imt; i++)
         for (k = 0; k < km; k++)
            if (k < KMT[j][i]) {
               if (dbg_lvl > 2)
                  printf ("(%d) i = %3d, j = %3d, k = %2d, tracer_state_ind = %d\n", iam, i, j, k, tracer_state_ind);
               int3_to_tracer_state_ind[k][j][i] = tracer_state_ind;
               tracer_state_ind_to_int3[tracer_state_ind].i = i;
               tracer_state_ind_to_int3[tracer_state_ind].j = j;
               tracer_state_ind_to_int3[tracer_state_ind].k = k;
               tracer_state_ind++;
            } else
               int3_to_tracer_state_ind[k][j][i] = -1;

   if (dbg_lvl > 1) {
      printf ("(%d) exiting %s\n", iam, subname);
      fflush (stdout);
   }

   return 0;
}

/******************************************************************************/

int
put_ind_maps (char *fname)
{
   char *subname = "put_ind_maps";
   int status;
   int ncid;
   int dimids[3];
   int tracer_state_len_dimid;
   int nlon_dimid;
   int nlat_dimid;
   int z_t_dimid;
   int varid;
   char *string;
   int attval_int;

   if (dbg_lvl > 1) {
      printf ("(%d) entering %s\n", iam, subname);
      fflush (stdout);
   }

   if ((status = nc_open (fname, NC_WRITE, &ncid)))
      return handle_nc_error (subname, "nc_open", fname, status);

   /* define new dimensions and get dimids for existing ones */

   if ((status = nc_redef (ncid)))
      return handle_nc_error (subname, "nc_redef", fname, status);

   if ((status = nc_def_dim (ncid, "tracer_state_len", tracer_state_len, &tracer_state_len_dimid)))
      return handle_nc_error (subname, "nc_def_dim", "tracer_state_len", status);

   if ((status = nc_inq_dimid (ncid, "nlon", &nlon_dimid)))
      return handle_nc_error (subname, "nc_inq_dimid", "nlon", status);

   if ((status = nc_inq_dimid (ncid, "nlat", &nlat_dimid)))
      return handle_nc_error (subname, "nc_inq_dimid", "nlat", status);

   if ((status = nc_inq_dimid (ncid, "z_t", &z_t_dimid)))
      return handle_nc_error (subname, "nc_inq_dimid", "z_t", status);

   /* define new variables */

   dimids[0] = z_t_dimid;
   dimids[1] = nlat_dimid;
   dimids[2] = nlon_dimid;

   if ((status = nc_def_var (ncid, "int3_to_tracer_state_ind", NC_INT, 3, dimids, &varid)))
      return handle_nc_error (subname, "nc_def_var", "int3_to_tracer_state_ind", status);
   string = "TLONG TLAT";
   if ((status = nc_put_att_text (ncid, varid, "coordinates", strlen (string), string)))
      return handle_nc_error (subname, "nc_put_att_text", "int3_to_tracer_state_ind", status);
   attval_int = -1;
   if ((status = nc_put_att_int (ncid, varid, "_FillValue", NC_INT, 1, &attval_int)))
      return handle_nc_error (subname, "nc_put_att_int", "int3_to_tracer_state_ind", status);
   if ((status = nc_put_att_int (ncid, varid, "missing_value", NC_INT, 1, &attval_int)))
      return handle_nc_error (subname, "nc_put_att_int", "int3_to_tracer_state_ind", status);

   dimids[0] = tracer_state_len_dimid;

   if ((status = nc_def_var (ncid, "tracer_state_ind_to_i", NC_INT, 1, dimids, &varid)))
      return handle_nc_error (subname, "nc_def_var", "int3_to_tracer_state_i", status);

   if ((status = nc_def_var (ncid, "tracer_state_ind_to_j", NC_INT, 1, dimids, &varid)))
      return handle_nc_error (subname, "nc_def_var", "int3_to_tracer_state_j", status);

   if ((status = nc_def_var (ncid, "tracer_state_ind_to_k", NC_INT, 1, dimids, &varid)))
      return handle_nc_error (subname, "nc_def_var", "int3_to_tracer_state_k", status);

   if ((status = nc_close (ncid)))
      return handle_nc_error (subname, "nc_close", fname, status);

   /* write out new variables */

   if (put_var_3d_int (fname, "int3_to_tracer_state_ind", int3_to_tracer_state_ind))
      return 1;

   /* write out tracer_state_ind_to_[ijk] variables, first copying to a temporary contiguous array tracer_state_ind_to_ijk */
   {
      int *tracer_state_ind_to_ijk;
      int tracer_state_ind;

      if ((tracer_state_ind_to_ijk = malloc ((size_t) tracer_state_len * sizeof (int))) == NULL) {
         fprintf (stderr, "(%d) malloc failed in %s for tracer_state_ind_to_ijk\n", iam, subname);
         return 1;
      }
      for (tracer_state_ind = 0; tracer_state_ind < tracer_state_len; tracer_state_ind++)
         tracer_state_ind_to_ijk[tracer_state_ind] = tracer_state_ind_to_int3[tracer_state_ind].i;
      if (put_var_1d_int (fname, "tracer_state_ind_to_i", tracer_state_ind_to_ijk))
         return 1;
      for (tracer_state_ind = 0; tracer_state_ind < tracer_state_len; tracer_state_ind++)
         tracer_state_ind_to_ijk[tracer_state_ind] = tracer_state_ind_to_int3[tracer_state_ind].j;
      if (put_var_1d_int (fname, "tracer_state_ind_to_j", tracer_state_ind_to_ijk))
         return 1;
      for (tracer_state_ind = 0; tracer_state_ind < tracer_state_len; tracer_state_ind++)
         tracer_state_ind_to_ijk[tracer_state_ind] = tracer_state_ind_to_int3[tracer_state_ind].k;
      if (put_var_1d_int (fname, "tracer_state_ind_to_k", tracer_state_ind_to_ijk))
         return 1;
      free (tracer_state_ind_to_ijk);
   }

   if (dbg_lvl > 1) {
      printf ("(%d) exiting %s\n", iam, subname);
      fflush (stdout);
   }

   return 0;
}

/******************************************************************************/

int
get_ind_maps (char *fname)
{
   char *subname = "get_ind_maps";
   int status;
   int ncid;
   int dimid;
   size_t dimlen;

   if (dbg_lvl > 1) {
      printf ("(%d) entering %s\n", iam, subname);
      fflush (stdout);
   }

   if (get_grid_dims (fname))
      return (1);

   if ((status = nc_open (fname, NC_NOWRITE, &ncid)))
      return handle_nc_error (subname, "nc_open", fname, status);

   if ((status = nc_inq_dimid (ncid, "tracer_state_len", &dimid)))
      return handle_nc_error (subname, "nc_inq_dimid", "tracer_state_len", status);

   if ((status = nc_inq_dimlen (ncid, dimid, &dimlen)))
      return handle_nc_error (subname, "nc_inq_dimlen", "tracer_state_len", status);
   tracer_state_len = dimlen;

   if ((status = nc_close (ncid)))
      return handle_nc_error (subname, "nc_close", fname, status);

   if (dbg_lvl && (iam == 0)) {
      printf ("(%d) %s: tracer_state_len = %d\n", iam, subname, tracer_state_len);
   }

   /* allocate space for ind_maps */
   if ((int3_to_tracer_state_ind = malloc_3d_int (km, jmt, imt)) == NULL) {
      fprintf (stderr, "(%d) malloc failed in %s for int3_to_tracer_state_ind\n", iam, subname);
      return 1;
   }
   if ((tracer_state_ind_to_int3 = malloc ((size_t) tracer_state_len * sizeof (int3))) == NULL) {
      fprintf (stderr, "(%d) malloc failed in %s for tracer_state_ind_to_int3\n", iam, subname);
      return 1;
   }

   /* read variables */
   /* read tracer_state_ind_to_[ijk] variables to a temporary contiguous array tracer_state_ind_to_ijk */

   if (get_var_3d_int (fname, "int3_to_tracer_state_ind", int3_to_tracer_state_ind))
      return 1;

   {
      int *tracer_state_ind_to_ijk;
      int tracer_state_ind;

      if ((tracer_state_ind_to_ijk = malloc ((size_t) tracer_state_len * sizeof (int))) == NULL) {
         fprintf (stderr, "(%d) malloc failed in %s for tracer_state_ind_to_ijk\n", iam, subname);
         return 1;
      }
      if (get_var_1d_int (fname, "tracer_state_ind_to_i", tracer_state_ind_to_ijk))
         return 1;
      for (tracer_state_ind = 0; tracer_state_ind < tracer_state_len; tracer_state_ind++)
         tracer_state_ind_to_int3[tracer_state_ind].i = tracer_state_ind_to_ijk[tracer_state_ind];

      if (get_var_1d_int (fname, "tracer_state_ind_to_j", tracer_state_ind_to_ijk))
         return 1;
      for (tracer_state_ind = 0; tracer_state_ind < tracer_state_len; tracer_state_ind++)
         tracer_state_ind_to_int3[tracer_state_ind].j = tracer_state_ind_to_ijk[tracer_state_ind];

      if (get_var_1d_int (fname, "tracer_state_ind_to_k", tracer_state_ind_to_ijk))
         return 1;
      for (tracer_state_ind = 0; tracer_state_ind < tracer_state_len; tracer_state_ind++)
         tracer_state_ind_to_int3[tracer_state_ind].k = tracer_state_ind_to_ijk[tracer_state_ind];

      free (tracer_state_ind_to_ijk);
   }

   if (dbg_lvl > 1) {
      printf ("(%d) exiting %s\n", iam, subname);
      fflush (stdout);
   }

   return 0;
}

/******************************************************************************/

void
free_ind_maps (void)
{
   free_3d_int (int3_to_tracer_state_ind);
   free (tracer_state_ind_to_int3);
}

/******************************************************************************/

void
comp_flat_len (void)
{
   flat_len = coupled_tracer_cnt * tracer_state_len;
   if (dbg_lvl)
      printf ("(%d) flat_len = %d\n\n", iam, flat_len);
}

/******************************************************************************/

int
adv_non_nbr_cnt (int k, int j, int i)
{
   int cnt;
   int ip1;
   int im1;
   int ip2;
   int im2;

   cnt = 0;
   ip1 = (i < imt - 1) ? i + 1 : 0;
   im1 = (i > 0) ? i - 1 : imt - 1;
   ip2 = (ip1 < imt - 1) ? ip1 + 1 : 0;
   im2 = (im1 > 0) ? im1 - 1 : imt - 1;

   if (adv_opt == adv_upwind3) {
      /* cell 2 level shallower */
      if (k - 2 >= 0)
         cnt++;
      /* cell 2 level deeper */
      if (k + 2 < KMT[j][i])
         cnt++;
      /* cell 2 unit east */
      if (k < KMT[j][ip2])
         cnt++;
      /* cell 2 unit west */
      if (k < KMT[j][im2])
         cnt++;
      /* cell 2 unit north */
      if ((j + 2 < jmt) && (k < KMT[j + 2][i]))
         cnt++;
      /* cell 2 unit south */
      if ((j - 2 >= 0) && (k < KMT[j - 2][i]))
         cnt++;
   }

   return cnt;
}

/******************************************************************************/

int
hmix_non_nbr_cnt (int k, int j, int i)
{
   int cnt;
   int ip1;
   int im1;

   cnt = 0;
   ip1 = (i < imt - 1) ? i + 1 : 0;
   im1 = (i > 0) ? i - 1 : imt - 1;

   if (hmix_opt == hmix_isop_file) {
      /* shallower & east */
      if ((k - 1 >= 0) && (k - 1 < KMT[j][ip1]))
         cnt++;
      /* deeper & east */
      if (k + 1 < KMT[j][ip1])
         cnt++;
      /* shallower & west */
      if ((k - 1 >= 0) && (k - 1 < KMT[j][im1]))
         cnt++;
      /* deeper & west */
      if (k + 1 < KMT[j][im1])
         cnt++;
      /* shallower & north */
      if ((k - 1 >= 0) && (k - 1 < KMT[j + 1][i]))
         cnt++;
      /* deeper & north */
      if (k + 1 < KMT[j + 1][i])
         cnt++;
      /* shallower & south */
      if ((k - 1 >= 0) && (k - 1 < KMT[j - 1][i]))
         cnt++;
      /* deeper & south */
      if (k + 1 < KMT[j - 1][i])
         cnt++;
   }

   return cnt;
}

/******************************************************************************/

int
vmix_non_nbr_cnt (int k, int j, int i)
{
   int cnt;

   if (vmix_opt == vmix_matrix_file)
      cnt = KMT[j][i];
   else
      cnt = 0;

   return cnt;
}

/******************************************************************************/

int
sink_non_nbr_cnt (int tracer_ind, int k, int j, int i)
{
   int cnt;

   cnt = 0;
   if (per_tracer_opt[tracer_ind].sink_opt == sink_generic_tracer) {
      int kmax =
         (per_tracer_opt[tracer_ind].sink_generic_tracer_depends_layer_cnt ==
          -1) ? km - 1 : per_tracer_opt[tracer_ind].sink_generic_tracer_depends_layer_cnt - 1;
      cnt = (k <= kmax) ? k + 1 : kmax + 1;
   }

   return cnt;
}

/******************************************************************************/

void
comp_nnz (void)
{
   char *subname = "comp_nnz";
   int tracer_ind;
   int tracer_state_ind;
   int i;
   int ip1;
   int im1;
   int j;
   int k;

   if (dbg_lvl > 1) {
      printf ("(%d) entering %s\n", iam, subname);
      fflush (stdout);
   }

   nnz = 0;
   for (tracer_ind = 0; tracer_ind < coupled_tracer_cnt; tracer_ind++) {
      for (tracer_state_ind = 0; tracer_state_ind < tracer_state_len; tracer_state_ind++) {
         i = tracer_state_ind_to_int3[tracer_state_ind].i;
         j = tracer_state_ind_to_int3[tracer_state_ind].j;
         k = tracer_state_ind_to_int3[tracer_state_ind].k;
         ip1 = (i < imt - 1) ? i + 1 : 0;
         im1 = (i > 0) ? i - 1 : imt - 1;

         /* cell itself */
         (nnz)++;
         /* cell 1 level shallower */
         if (k - 1 >= 0)
            (nnz)++;
         /* cell 1 level deeper */
         if (k + 1 < KMT[j][i])
            (nnz)++;
         /* cell 1 unit east */
         if (k < KMT[j][ip1])
            (nnz)++;
         /* cell 1 unit west */
         if (k < KMT[j][im1])
            (nnz)++;
         /* cell 1 unit north */
         if (k < KMT[j + 1][i])
            (nnz)++;
         /* cell 1 unit south */
         if (k < KMT[j - 1][i])
            (nnz)++;

         nnz += adv_non_nbr_cnt (k, j, i);

         nnz += hmix_non_nbr_cnt (k, j, i);

         nnz += vmix_non_nbr_cnt (k, j, i);

         nnz += sink_non_nbr_cnt (tracer_ind, k, j, i);

         nnz += coupled_tracer_cnt - 1;
      }
   }

   if (dbg_lvl)
      printf ("(%d) nnz       = %d\n\n", iam, nnz);

   if (dbg_lvl > 1) {
      printf ("(%d) exiting %s\n", iam, subname);
      fflush (stdout);
   }
}

/******************************************************************************/

int
allocate_matrix_arrays (void)
{
   char *subname = "allocate_matrix_arrays";
   int tracer_ind;

   if (dbg_lvl > 1) {
      printf ("(%d) entering %s\n", iam, subname);
      fflush (stdout);
   }

   if ((nzval_row_wise = malloc ((size_t) nnz * sizeof (double))) == NULL) {
      fprintf (stderr, "(%d) malloc failed in %s for nzval_row_wise\n", iam, subname);
      return 1;
   }
   if ((colind = malloc ((size_t) nnz * sizeof (int_t))) == NULL) {
      fprintf (stderr, "(%d) malloc failed in %s for colind\n", iam, subname);
      return 1;
   }
   if ((rowptr = malloc ((size_t) (flat_len + 1) * sizeof (int_t))) == NULL) {
      fprintf (stderr, "(%d) malloc failed in %s for rowptr\n", iam, subname);
      return 1;
   }

   if ((coef_ind_self = malloc ((size_t) coupled_tracer_cnt * sizeof (int_t *))) == NULL) {
      fprintf (stderr, "(%d) malloc failed in %s for coef_ind_self\n", iam, subname);
      return 1;
   }
   if ((coef_ind_adv_non_nbr = malloc ((size_t) coupled_tracer_cnt * sizeof (int_t *))) == NULL) {
      fprintf (stderr, "(%d) malloc failed in %s for coef_ind_adv_non_nbr\n", iam, subname);
      return 1;
   }
   if ((coef_ind_hmix_non_nbr = malloc ((size_t) coupled_tracer_cnt * sizeof (int_t *))) == NULL) {
      fprintf (stderr, "(%d) malloc failed in %s for coef_ind_hmix_non_nbr\n", iam, subname);
      return 1;
   }
   if ((coef_ind_vmix_non_nbr = malloc ((size_t) coupled_tracer_cnt * sizeof (int_t *))) == NULL) {
      fprintf (stderr, "(%d) malloc failed in %s for coef_ind_vmix_non_nbr\n", iam, subname);
      return 1;
   }
   if ((coef_ind_sink_non_nbr = malloc ((size_t) coupled_tracer_cnt * sizeof (int_t *))) == NULL) {
      fprintf (stderr, "(%d) malloc failed in %s for coef_ind_sink_non_nbr\n", iam, subname);
      return 1;
   }
   if ((coef_ind_sink_other_tracers = malloc ((size_t) coupled_tracer_cnt * sizeof (int_t *))) == NULL) {
      fprintf (stderr, "(%d) malloc failed in %s for coef_ind_sink_other_tracers\n", iam, subname);
      return 1;
   }

   for (tracer_ind = 0; tracer_ind < coupled_tracer_cnt; tracer_ind++) {
      if ((coef_ind_self[tracer_ind] = malloc ((size_t) tracer_state_len * sizeof (int_t))) == NULL) {
         fprintf (stderr, "(%d) malloc failed in %s for coef_ind_self[%d]\n", iam, subname, tracer_ind);
         return 1;
      }
      if ((coef_ind_adv_non_nbr[tracer_ind] = malloc ((size_t) tracer_state_len * sizeof (int_t))) == NULL) {
         fprintf (stderr, "(%d) malloc failed in %s for coef_ind_adv_non_nbr[%d]\n", iam, subname, tracer_ind);
         return 1;
      }
      if ((coef_ind_hmix_non_nbr[tracer_ind] = malloc ((size_t) tracer_state_len * sizeof (int_t))) == NULL) {
         fprintf (stderr, "(%d) malloc failed in %s for coef_ind_hmix_non_nbr[%d]\n", iam, subname, tracer_ind);
         return 1;
      }
      if ((coef_ind_vmix_non_nbr[tracer_ind] = malloc ((size_t) tracer_state_len * sizeof (int_t))) == NULL) {
         fprintf (stderr, "(%d) malloc failed in %s for coef_ind_vmix_non_nbr[%d]\n", iam, subname, tracer_ind);
         return 1;
      }
      if ((coef_ind_sink_non_nbr[tracer_ind] = malloc ((size_t) tracer_state_len * sizeof (int_t))) == NULL) {
         fprintf (stderr, "(%d) malloc failed in %s for coef_ind_sink_non_nbr[%d]\n", iam, subname, tracer_ind);
         return 1;
      }
      if ((coef_ind_sink_other_tracers[tracer_ind] = malloc ((size_t) tracer_state_len * sizeof (int_t))) == NULL) {
         fprintf (stderr, "(%d) malloc failed in %s for coef_ind_sink_other_tracers[%d]\n", iam, subname, tracer_ind);
         return 1;
      }
   }

   if (dbg_lvl > 1) {
      printf ("(%d) exiting %s\n", iam, subname);
      fflush (stdout);
   }

   return 0;
}

/******************************************************************************/

/* initialize matrix values to zero and set up sparsity pattern arrays */

int
init_matrix (void)
{
   char *subname = "init_matrix";
   int coef_ind;
   int tracer_ind;
   int tracer_ind_2;
   int flat_ind_offset;
   int tracer_state_ind;
   int flat_ind;
   int i;
   int ip1;
   int im1;
   int ip2;
   int im2;
   int j;
   int k;

   if (dbg_lvl > 1) {
      printf ("(%d) entering %s\n", iam, subname);
      fflush (stdout);
   }

   coef_ind = 0;
   for (tracer_ind = 0; tracer_ind < coupled_tracer_cnt; tracer_ind++) {
      flat_ind_offset = tracer_ind * tracer_state_len;
      if (dbg_lvl > 1) {
         printf ("(%d) tracer_ind = %d, flat_ind_offset = %d\n", iam, tracer_ind, flat_ind_offset);
         fflush (stdout);
      }
      for (tracer_state_ind = 0; tracer_state_ind < tracer_state_len; tracer_state_ind++) {
         flat_ind = flat_ind_offset + tracer_state_ind;
         rowptr[flat_ind] = coef_ind;

         if (dbg_lvl > 2) {
            printf ("(%d) rowptr[%d]=%d\n", iam, flat_ind, coef_ind);
            fflush (stdout);
         }

         i = tracer_state_ind_to_int3[tracer_state_ind].i;
         j = tracer_state_ind_to_int3[tracer_state_ind].j;
         k = tracer_state_ind_to_int3[tracer_state_ind].k;
         ip1 = (i < imt - 1) ? i + 1 : 0;
         im1 = (i > 0) ? i - 1 : imt - 1;
         ip2 = (ip1 < imt - 1) ? ip1 + 1 : 0;
         im2 = (im1 > 0) ? im1 - 1 : imt - 1;

         /* cell itself */
         coef_ind_self[tracer_ind][tracer_state_ind] = coef_ind;
         nzval_row_wise[coef_ind] = 0.0;
         colind[coef_ind] = flat_ind_offset + int3_to_tracer_state_ind[k][j][i];
         coef_ind++;
         /* cell 1 level shallower */
         if (k - 1 >= 0) {
            nzval_row_wise[coef_ind] = 0.0;
            colind[coef_ind] = flat_ind_offset + int3_to_tracer_state_ind[k - 1][j][i];
            coef_ind++;
         }
         /* cell 1 level deeper */
         if (k + 1 < KMT[j][i]) {
            nzval_row_wise[coef_ind] = 0.0;
            colind[coef_ind] = flat_ind_offset + int3_to_tracer_state_ind[k + 1][j][i];
            coef_ind++;
         }
         /* cell 1 unit east */
         if (k < KMT[j][ip1]) {
            nzval_row_wise[coef_ind] = 0.0;
            colind[coef_ind] = flat_ind_offset + int3_to_tracer_state_ind[k][j][ip1];
            coef_ind++;
         }
         /* cell 1 unit west */
         if (k < KMT[j][im1]) {
            nzval_row_wise[coef_ind] = 0.0;
            colind[coef_ind] = flat_ind_offset + int3_to_tracer_state_ind[k][j][im1];
            coef_ind++;
         }
         /* cell 1 unit north */
         if (k < KMT[j + 1][i]) {
            nzval_row_wise[coef_ind] = 0.0;
            colind[coef_ind] = flat_ind_offset + int3_to_tracer_state_ind[k][j + 1][i];
            coef_ind++;
         }
         /* cell 1 unit south */
         if (k < KMT[j - 1][i]) {
            nzval_row_wise[coef_ind] = 0.0;
            colind[coef_ind] = flat_ind_offset + int3_to_tracer_state_ind[k][j - 1][i];
            coef_ind++;
         }
         coef_ind_adv_non_nbr[tracer_ind][tracer_state_ind] = coef_ind;
         if (adv_opt == adv_upwind3) {
            /* cell 2 level shallower */
            if (k - 2 >= 0) {
               nzval_row_wise[coef_ind] = 0.0;
               colind[coef_ind] = flat_ind_offset + int3_to_tracer_state_ind[k - 2][j][i];
               coef_ind++;
            }
            /* cell 2 level deeper */
            if (k + 2 < KMT[j][i]) {
               nzval_row_wise[coef_ind] = 0.0;
               colind[coef_ind] = flat_ind_offset + int3_to_tracer_state_ind[k + 2][j][i];
               coef_ind++;
            }
            /* cell 2 unit east */
            if (k < KMT[j][ip2]) {
               nzval_row_wise[coef_ind] = 0.0;
               colind[coef_ind] = flat_ind_offset + int3_to_tracer_state_ind[k][j][ip2];
               coef_ind++;
            }
            /* cell 2 unit west */
            if (k < KMT[j][im2]) {
               nzval_row_wise[coef_ind] = 0.0;
               colind[coef_ind] = flat_ind_offset + int3_to_tracer_state_ind[k][j][im2];
               coef_ind++;
            }
            /* cell 2 unit north */
            if ((j + 2 < jmt) && (k < KMT[j + 2][i])) {
               nzval_row_wise[coef_ind] = 0.0;
               colind[coef_ind] = flat_ind_offset + int3_to_tracer_state_ind[k][j + 2][i];
               coef_ind++;
            }
            /* cell 2 unit south */
            if ((j - 2 >= 0) && (k < KMT[j - 2][i])) {
               nzval_row_wise[coef_ind] = 0.0;
               colind[coef_ind] = flat_ind_offset + int3_to_tracer_state_ind[k][j - 2][i];
               coef_ind++;
            }
         }
         coef_ind_hmix_non_nbr[tracer_ind][tracer_state_ind] = coef_ind;
         if (hmix_opt == hmix_isop_file) {
            /* shallower & east */
            if ((k - 1 >= 0) && (k - 1 < KMT[j][ip1])) {
               nzval_row_wise[coef_ind] = 0.0;
               colind[coef_ind] = flat_ind_offset + int3_to_tracer_state_ind[k - 1][j][ip1];
               coef_ind++;
            }
            /* deeper & east */
            if (k + 1 < KMT[j][ip1]) {
               nzval_row_wise[coef_ind] = 0.0;
               colind[coef_ind] = flat_ind_offset + int3_to_tracer_state_ind[k + 1][j][ip1];
               coef_ind++;
            }
            /* shallower & west */
            if ((k - 1 >= 0) && (k - 1 < KMT[j][im1])) {
               nzval_row_wise[coef_ind] = 0.0;
               colind[coef_ind] = flat_ind_offset + int3_to_tracer_state_ind[k - 1][j][im1];
               coef_ind++;
            }
            /* deeper & west */
            if (k + 1 < KMT[j][im1]) {
               nzval_row_wise[coef_ind] = 0.0;
               colind[coef_ind] = flat_ind_offset + int3_to_tracer_state_ind[k + 1][j][im1];
               coef_ind++;
            }
            /* shallower & north */
            if ((k - 1 >= 0) && (k - 1 < KMT[j + 1][i])) {
               nzval_row_wise[coef_ind] = 0.0;
               colind[coef_ind] = flat_ind_offset + int3_to_tracer_state_ind[k - 1][j + 1][i];
               coef_ind++;
            }
            /* deeper & north */
            if (k + 1 < KMT[j + 1][i]) {
               nzval_row_wise[coef_ind] = 0.0;
               colind[coef_ind] = flat_ind_offset + int3_to_tracer_state_ind[k + 1][j + 1][i];
               coef_ind++;
            }
            /* shallower & south */
            if ((k - 1 >= 0) && (k - 1 < KMT[j - 1][i])) {
               nzval_row_wise[coef_ind] = 0.0;
               colind[coef_ind] = flat_ind_offset + int3_to_tracer_state_ind[k - 1][j - 1][i];
               coef_ind++;
            }
            /* deeper & south */
            if (k + 1 < KMT[j - 1][i]) {
               nzval_row_wise[coef_ind] = 0.0;
               colind[coef_ind] = flat_ind_offset + int3_to_tracer_state_ind[k + 1][j - 1][i];
               coef_ind++;
            }
         }
         coef_ind_vmix_non_nbr[tracer_ind][tracer_state_ind] = coef_ind;
         if (vmix_opt == vmix_matrix_file) {
            int k2;

            for (k2 = 0; k2 < KMT[j][i]; k2++) {
               nzval_row_wise[coef_ind] = 0.0;
               colind[coef_ind] = flat_ind_offset + int3_to_tracer_state_ind[k2][j][i];
               coef_ind++;
            }
         }
         coef_ind_sink_non_nbr[tracer_ind][tracer_state_ind] = coef_ind;
         if (per_tracer_opt[tracer_ind].sink_opt == sink_generic_tracer) {
            int k2;
            int kmax =
               (per_tracer_opt[tracer_ind].sink_generic_tracer_depends_layer_cnt ==
                -1) ? km - 1 : per_tracer_opt[tracer_ind].sink_generic_tracer_depends_layer_cnt - 1;

            for (k2 = (k <= kmax) ? k : kmax; k2 >= 0; k2--) {
               nzval_row_wise[coef_ind] = 0.0;
               colind[coef_ind] = flat_ind_offset + int3_to_tracer_state_ind[k2][j][i];
               coef_ind++;
            }
         }
         coef_ind_sink_other_tracers[tracer_ind][tracer_state_ind] = coef_ind;
         for (tracer_ind_2 = 0; tracer_ind_2 < coupled_tracer_cnt; tracer_ind_2++) {
            if (tracer_ind_2 == tracer_ind)
               continue;
            nzval_row_wise[coef_ind] = 0.0;
            colind[coef_ind] = tracer_ind_2 * tracer_state_len + int3_to_tracer_state_ind[k][j][i];
            coef_ind++;
         }
      }
   }
   if (coef_ind != nnz) {
      fprintf (stderr, "(%d) internal error in %s, coef_ind != nnz after setting sparsity pattern\n", iam, subname);
      fprintf (stderr, "(%d) coef_ind = %d\nnnz      = %d\n", iam, coef_ind, nnz);
      return 1;
   }
   rowptr[flat_len] = coef_ind;
   if (dbg_lvl > 2) {
      printf ("(%d) rowptr[%d]=%d\n", iam, flat_len, coef_ind);
      fflush (stdout);
   }

   if (dbg_lvl > 1) {
      printf ("(%d) exiting %s\n", iam, subname);
      fflush (stdout);
   }

   return 0;
}

/******************************************************************************/

int
load_UTE (double ***UTE)
{
   char *subname = "load_UTE";
   double ***WORK;
   double **DY;
   double fv;
   int i;
   int ip1;
   int j;
   int k;

   if (dbg_lvl > 1) {
      printf ("(%d) entering %s\n", iam, subname);
      fflush (stdout);
   }

   set_3d_double (0.0, UTE);

   if ((WORK = malloc_3d_double (km, jmt, imt)) == NULL) {
      fprintf (stderr, "(%d) malloc failed in %s for WORK\n", iam, subname);
      return 1;
   }
   if ((DY = malloc_2d_double (jmt, imt)) == NULL) {
      fprintf (stderr, "(%d) malloc failed in %s for DY\n", iam, subname);
      return 1;
   }
   if (dbg_lvl)
      printf ("(%d) %s: reading UVEL,DYU from %s\n", iam, subname, circ_fname);
   if (get_var_3d_double (circ_fname, "UVEL", WORK))
      return 1;
   if (get_att_double (circ_fname, "UVEL", "_FillValue", &fv))
      return 1;
   set_fv_3d_double (fv, 0.0, WORK);
   if (get_var_2d_double (circ_fname, "DYU", DY))
      return 1;
   if (get_att_double (circ_fname, "DYU", "_FillValue", &fv))
      return 1;
   set_fv_2d_double (fv, 0.0, DY);
   for (k = 0; k < km; k++)
      for (j = 1; j < jmt - 1; j++)
         for (i = 0; i < imt; i++) {
            if (k < KMU[j][i])
               UTE[k][j][i] += 0.5 * WORK[k][j][i] * DY[j][i];
            if (k < KMU[j - 1][i])
               UTE[k][j][i] += 0.5 * WORK[k][j - 1][i] * DY[j - 1][i];
         }

   if (hmix_opt == hmix_hor_file) {
      if (dbg_lvl)
         printf ("(%d) %s: reading UISOP,HTE from %s\n", iam, subname, circ_fname);
      if (get_var_3d_double (circ_fname, "UISOP", WORK))
         return 1;
      if (get_var_2d_double (circ_fname, "HTE", DY))
         return 1;
      if (get_att_double (circ_fname, "HTE", "_FillValue", &fv))
         return 1;
      set_fv_2d_double (fv, 0.0, DY);
      for (k = 0; k < km; k++)
         for (j = 1; j < jmt - 1; j++)
            for (i = 0; i < imt; i++) {
               ip1 = (i < imt - 1) ? i + 1 : 0;
               if ((k < KMT[j][i]) && (k < KMT[j][ip1]))
                  UTE[k][j][i] += WORK[k][j][i] * DY[j][i];
            }
   }
   free_2d_double (DY);
   free_3d_double (WORK);

   if (dbg_lvl > 1) {
      printf ("(%d) exiting %s\n", iam, subname);
      fflush (stdout);
   }

   return 0;
}

/******************************************************************************/

int
load_VTN (double ***VTN)
{
   char *subname = "load_VTN";
   double ***WORK;
   double **DX;
   double fv;
   int i;
   int im1;
   int j;
   int k;

   if (dbg_lvl > 1) {
      printf ("(%d) entering %s\n", iam, subname);
      fflush (stdout);
   }

   set_3d_double (0.0, VTN);

   if ((WORK = malloc_3d_double (km, jmt, imt)) == NULL) {
      fprintf (stderr, "(%d) malloc failed in %s for WORK\n", iam, subname);
      return 1;
   }
   if ((DX = malloc_2d_double (jmt, imt)) == NULL) {
      fprintf (stderr, "(%d) malloc failed in %s for DX\n", iam, subname);
      return 1;
   }
   if (dbg_lvl)
      printf ("(%d) %s: reading VVEL,DXU from %s\n", iam, subname, circ_fname);
   if (get_var_3d_double (circ_fname, "VVEL", WORK))
      return 1;
   if (get_att_double (circ_fname, "VVEL", "_FillValue", &fv))
      return 1;
   set_fv_3d_double (fv, 0.0, WORK);
   if (get_var_2d_double (circ_fname, "DXU", DX))
      return 1;
   if (get_att_double (circ_fname, "DXU", "_FillValue", &fv))
      return 1;
   set_fv_2d_double (fv, 0.0, DX);
   for (k = 0; k < km; k++)
      for (j = 1; j < jmt - 1; j++)
         for (i = 0; i < imt; i++) {
            im1 = (i > 0) ? i - 1 : imt - 1;
            if (k < KMU[j][i])
               VTN[k][j][i] += 0.5 * WORK[k][j][i] * DX[j][i];
            if (k < KMU[j][im1])
               VTN[k][j][i] += 0.5 * WORK[k][j][im1] * DX[j][im1];
         }

   if (hmix_opt == hmix_hor_file) {
      if (dbg_lvl)
         printf ("(%d) %s: reading VISOP,HTN from %s\n", iam, subname, circ_fname);
      if (get_var_3d_double (circ_fname, "VISOP", WORK))
         return 1;
      if (get_att_double (circ_fname, "VISOP", "_FillValue", &fv))
         return 1;
      set_fv_3d_double (fv, 0.0, WORK);
      if (get_var_2d_double (circ_fname, "HTN", DX))
         return 1;
      if (get_att_double (circ_fname, "HTN", "_FillValue", &fv))
         return 1;
      set_fv_2d_double (fv, 0.0, DX);
      for (k = 0; k < km; k++)
         for (j = 1; j < jmt - 1; j++)
            for (i = 0; i < imt; i++)
               if ((k < KMT[j][i]) && (k < KMT[j + 1][i]))
                  VTN[k][j][i] += WORK[k][j][i] * DX[j][i];
   }
   free_2d_double (DX);
   free_3d_double (WORK);

   if (dbg_lvl > 1) {
      printf ("(%d) exiting %s\n", iam, subname);
      fflush (stdout);
   }

   return 0;
}

/******************************************************************************/

int
load_WVEL (double ***WVEL)
{
   char *subname = "load_WVEL";
   double ***WORK;
   double fv;
   int i;
   int j;
   int k;

   if (dbg_lvl > 1) {
      printf ("(%d) entering %s\n", iam, subname);
      fflush (stdout);
   }

   set_3d_double (0.0, WVEL);

   if ((WORK = malloc_3d_double (km, jmt, imt)) == NULL) {
      fprintf (stderr, "(%d) malloc failed in %s for WORK\n", iam, subname);
      return 1;
   }
   if (dbg_lvl)
      printf ("(%d) %s: reading WVEL from %s\n", iam, subname, circ_fname);
   if (get_var_3d_double (circ_fname, "WVEL", WORK))
      return 1;
   if (get_att_double (circ_fname, "WVEL", "_FillValue", &fv))
      return 1;
   set_fv_3d_double (fv, 0.0, WORK);
   for (k = 0; k < km; k++)
      for (j = 1; j < jmt - 1; j++)
         for (i = 0; i < imt; i++)
            if (k < KMT[j][i])
               WVEL[k][j][i] += WORK[k][j][i];

   if (hmix_opt == hmix_hor_file) {
      if (dbg_lvl)
         printf ("(%d) %s: reading WISOP from %s\n", iam, subname, circ_fname);
      if (get_var_3d_double (circ_fname, "WISOP", WORK))
         return 1;
      if (get_att_double (circ_fname, "WISOP", "_FillValue", &fv))
         return 1;
      set_fv_3d_double (fv, 0.0, WORK);
      for (k = 0; k < km; k++)
         for (j = 1; j < jmt - 1; j++)
            for (i = 0; i < imt; i++)
               if (k < KMT[j][i])
                  WVEL[k][j][i] += WORK[k][j][i];
   }
   free_3d_double (WORK);

   /* explicitly set surface velocity to zero */
   for (j = 1; j < jmt - 1; j++)
      for (i = 0; i < imt; i++)
         WVEL[0][j][i] = 0.0;

   if (dbg_lvl > 1) {
      printf ("(%d) exiting %s\n", iam, subname);
      fflush (stdout);
   }

   return 0;
}

/******************************************************************************/

void
add_UTE_coeffs (double ***UTE)
{
   char *subname = "add_UTE_coeffs";
   int tracer_ind;
   int tracer_state_ind;
   double east_self_interp_w;
   double west_self_interp_w;

   if (dbg_lvl > 1) {
      printf ("(%d) entering %s\n", iam, subname);
      fflush (stdout);
   }

   for (tracer_ind = 0; tracer_ind < coupled_tracer_cnt; tracer_ind++) {
      for (tracer_state_ind = 0; tracer_state_ind < tracer_state_len; tracer_state_ind++) {
         int i;
         int ip1;
         int im1;
         int j;
         int k;
         int coef_ind;

         i = tracer_state_ind_to_int3[tracer_state_ind].i;
         j = tracer_state_ind_to_int3[tracer_state_ind].j;
         k = tracer_state_ind_to_int3[tracer_state_ind].k;
         ip1 = (i < imt - 1) ? i + 1 : 0;
         im1 = (i > 0) ? i - 1 : imt - 1;

         switch (adv_opt) {
         case adv_donor:
            east_self_interp_w = (UTE[k][j][i] > 0.0) ? 1.0 : 0.0;
            west_self_interp_w = (UTE[k][j][im1] < 0.0) ? 1.0 : 0.0;
            break;
         case adv_cent:
            east_self_interp_w = 0.5;
            west_self_interp_w = 0.5;
            break;
         }

         coef_ind = coef_ind_self[tracer_ind][tracer_state_ind];

         /* cell itself */
         if (k < KMT[j][ip1])
            nzval_row_wise[coef_ind] -= east_self_interp_w * UTE[k][j][i] / TAREA[j][i] * delta_t;
         if (k < KMT[j][im1])
            nzval_row_wise[coef_ind] += west_self_interp_w * UTE[k][j][im1] / TAREA[j][i] * delta_t;
         coef_ind++;
         /* cell 1 level shallower */
         if (k - 1 >= 0)
            coef_ind++;
         /* cell 1 level deeper */
         if (k + 1 < KMT[j][i])
            coef_ind++;
         /* cell 1 unit east */
         if (k < KMT[j][ip1]) {
            nzval_row_wise[coef_ind] -= (1.0 - east_self_interp_w) * UTE[k][j][i] / TAREA[j][i] * delta_t;
            coef_ind++;
         }
         /* cell 1 unit west */
         if (k < KMT[j][im1]) {
            nzval_row_wise[coef_ind] += (1.0 - west_self_interp_w) * UTE[k][j][im1] / TAREA[j][i] * delta_t;
            coef_ind++;
         }
         /* cell 1 unit north */
         if (k < KMT[j + 1][i])
            coef_ind++;
         /* cell 1 unit south */
         if (k < KMT[j - 1][i])
            coef_ind++;
      }
   }

   if (dbg_lvl > 1) {
      printf ("(%d) exiting %s\n", iam, subname);
      fflush (stdout);
   }
}

/******************************************************************************/

void
add_VTN_coeffs (double ***VTN)
{
   char *subname = "add_VTN_coeffs";
   int tracer_ind;
   int tracer_state_ind;
   double north_self_interp_w;
   double south_self_interp_w;

   if (dbg_lvl > 1) {
      printf ("(%d) entering %s\n", iam, subname);
      fflush (stdout);
   }

   for (tracer_ind = 0; tracer_ind < coupled_tracer_cnt; tracer_ind++) {
      for (tracer_state_ind = 0; tracer_state_ind < tracer_state_len; tracer_state_ind++) {
         int i;
         int ip1;
         int im1;
         int j;
         int k;
         int coef_ind;

         i = tracer_state_ind_to_int3[tracer_state_ind].i;
         j = tracer_state_ind_to_int3[tracer_state_ind].j;
         k = tracer_state_ind_to_int3[tracer_state_ind].k;
         ip1 = (i < imt - 1) ? i + 1 : 0;
         im1 = (i > 0) ? i - 1 : imt - 1;

         switch (adv_opt) {
         case adv_donor:
            north_self_interp_w = (VTN[k][j][i] > 0.0) ? 1.0 : 0.0;
            south_self_interp_w = (VTN[k][j - 1][i] < 0.0) ? 1.0 : 0.0;
            break;
         case adv_cent:
            north_self_interp_w = 0.5;
            south_self_interp_w = 0.5;
            break;
         }

         coef_ind = coef_ind_self[tracer_ind][tracer_state_ind];

         /* cell itself */
         if (k < KMT[j + 1][i])
            nzval_row_wise[coef_ind] -= north_self_interp_w * VTN[k][j][i] / TAREA[j][i] * delta_t;
         if (k < KMT[j - 1][i])
            nzval_row_wise[coef_ind] += south_self_interp_w * VTN[k][j - 1][i] / TAREA[j][i] * delta_t;
         coef_ind++;
         /* cell 1 level shallower */
         if (k - 1 >= 0)
            coef_ind++;
         /* cell 1 level deeper */
         if (k + 1 < KMT[j][i])
            coef_ind++;
         /* cell 1 unit east */
         if (k < KMT[j][ip1])
            coef_ind++;
         /* cell 1 unit west */
         if (k < KMT[j][im1])
            coef_ind++;
         /* cell 1 unit north */
         if (k < KMT[j + 1][i]) {
            nzval_row_wise[coef_ind] -= (1.0 - north_self_interp_w) * VTN[k][j][i] / TAREA[j][i] * delta_t;
            coef_ind++;
         }
         /* cell 1 unit south */
         if (k < KMT[j - 1][i]) {
            nzval_row_wise[coef_ind] += (1.0 - south_self_interp_w) * VTN[k][j - 1][i] / TAREA[j][i] * delta_t;
            coef_ind++;
         }
      }
   }

   if (dbg_lvl > 1) {
      printf ("(%d) exiting %s\n", iam, subname);
      fflush (stdout);
   }
}

/******************************************************************************/

void
add_WVEL_coeffs (double ***WVEL)
{
   char *subname = "add_WVEL_coeffs";
   int tracer_ind;
   int tracer_state_ind;
   double top_self_interp_w;
   double bot_self_interp_w;

   if (dbg_lvl > 1) {
      printf ("(%d) entering %s\n", iam, subname);
      fflush (stdout);
   }

   for (tracer_ind = 0; tracer_ind < coupled_tracer_cnt; tracer_ind++) {
      for (tracer_state_ind = 0; tracer_state_ind < tracer_state_len; tracer_state_ind++) {
         int i;
         int ip1;
         int im1;
         int j;
         int k;
         int coef_ind;

         i = tracer_state_ind_to_int3[tracer_state_ind].i;
         j = tracer_state_ind_to_int3[tracer_state_ind].j;
         k = tracer_state_ind_to_int3[tracer_state_ind].k;
         ip1 = (i < imt - 1) ? i + 1 : 0;
         im1 = (i > 0) ? i - 1 : imt - 1;

         switch (adv_opt) {
         case adv_donor:
            top_self_interp_w = (WVEL[k][j][i] > 0.0) ? 1.0 : 0.0;
            if (k + 1 < KMT[j][i])
               bot_self_interp_w = (WVEL[k + 1][j][i] < 0.0) ? 1.0 : 0.0;
            break;
         case adv_cent:
            top_self_interp_w = 0.5;
            bot_self_interp_w = 0.5;
            break;
         }

         coef_ind = coef_ind_self[tracer_ind][tracer_state_ind];

         /* cell itself */
         if (k - 1 >= 0)
            nzval_row_wise[coef_ind] -= top_self_interp_w * WVEL[k][j][i] / dz[k] * delta_t;
         if (k + 1 < KMT[j][i])
            nzval_row_wise[coef_ind] += bot_self_interp_w * WVEL[k + 1][j][i] / dz[k] * delta_t;
         coef_ind++;
         /* cell 1 level shallower */
         if (k - 1 >= 0) {
            nzval_row_wise[coef_ind] -= (1.0 - top_self_interp_w) * WVEL[k][j][i] / dz[k] * delta_t;
            coef_ind++;
         }
         /* cell 1 level deeper */
         if (k + 1 < KMT[j][i]) {
            nzval_row_wise[coef_ind] += (1.0 - bot_self_interp_w) * WVEL[k + 1][j][i] / dz[k] * delta_t;
            coef_ind++;
         }
         /* cell 1 unit east */
         if (k < KMT[j][ip1])
            coef_ind++;
         /* cell 1 unit west */
         if (k < KMT[j][im1])
            coef_ind++;
         /* cell 1 unit north */
         if (k < KMT[j + 1][i])
            coef_ind++;
         /* cell 1 unit south */
         if (k < KMT[j - 1][i])
            coef_ind++;
      }
   }

   if (dbg_lvl > 1) {
      printf ("(%d) exiting %s\n", iam, subname);
      fflush (stdout);
   }
}

/******************************************************************************/

int
load_UTE_upwind3 (double ***UTE_POS, double ***UTE_NEG)
{
   char *subname = "load_UTE_upwind3";
   double fv;

   if (dbg_lvl > 1) {
      printf ("(%d) entering %s\n", iam, subname);
      fflush (stdout);
   }

   set_3d_double (0.0, UTE_POS);
   set_3d_double (0.0, UTE_NEG);

   if (dbg_lvl)
      printf ("(%d) %s: reading UTE_{POS,NEG} from %s\n", iam, subname, circ_fname);

   if (get_var_3d_double (circ_fname, "UTE_POS", UTE_POS))
      return 1;
   if (get_att_double (circ_fname, "UTE_POS", "_FillValue", &fv))
      return 1;
   set_fv_3d_double (fv, 0.0, UTE_POS);
   if (get_var_3d_double (circ_fname, "UTE_NEG", UTE_NEG))
      return 1;
   if (get_att_double (circ_fname, "UTE_NEG", "_FillValue", &fv))
      return 1;
   set_fv_3d_double (fv, 0.0, UTE_NEG);

   if (dbg_lvl > 1) {
      printf ("(%d) exiting %s\n", iam, subname);
      fflush (stdout);
   }

   return 0;
}

/******************************************************************************/

int
load_VTN_upwind3 (double ***VTN_POS, double ***VTN_NEG)
{
   char *subname = "load_VTN_upwind3";
   double fv;

   if (dbg_lvl > 1) {
      printf ("(%d) entering %s\n", iam, subname);
      fflush (stdout);
   }

   set_3d_double (0.0, VTN_POS);
   set_3d_double (0.0, VTN_NEG);

   if (dbg_lvl)
      printf ("(%d) %s: reading VTN_{POS,NEG} from %s\n", iam, subname, circ_fname);

   if (get_var_3d_double (circ_fname, "VTN_POS", VTN_POS))
      return 1;
   if (get_att_double (circ_fname, "VTN_POS", "_FillValue", &fv))
      return 1;
   set_fv_3d_double (fv, 0.0, VTN_POS);
   if (get_var_3d_double (circ_fname, "VTN_NEG", VTN_NEG))
      return 1;
   if (get_att_double (circ_fname, "VTN_NEG", "_FillValue", &fv))
      return 1;
   set_fv_3d_double (fv, 0.0, VTN_NEG);

   if (dbg_lvl > 1) {
      printf ("(%d) exiting %s\n", iam, subname);
      fflush (stdout);
   }

   return 0;
}

/******************************************************************************/

int
load_WVEL_upwind3 (double ***WVEL_POS, double ***WVEL_NEG)
{
   char *subname = "load_WVEL_upwind3";
   double fv;
   int i;
   int j;

   if (dbg_lvl > 1) {
      printf ("(%d) entering %s\n", iam, subname);
      fflush (stdout);
   }

   set_3d_double (0.0, WVEL_POS);
   set_3d_double (0.0, WVEL_NEG);

   if (dbg_lvl)
      printf ("(%d) %s: reading WTK_{POS,NEG} from %s\n", iam, subname, circ_fname);

   if (get_var_3d_double (circ_fname, "WTK_POS", WVEL_POS))
      return 1;
   if (get_att_double (circ_fname, "WTK_POS", "_FillValue", &fv))
      return 1;
   set_fv_3d_double (fv, 0.0, WVEL_POS);
   if (get_var_3d_double (circ_fname, "WTK_NEG", WVEL_NEG))
      return 1;
   if (get_att_double (circ_fname, "WTK_NEG", "_FillValue", &fv))
      return 1;
   set_fv_3d_double (fv, 0.0, WVEL_NEG);

   /* explicitly set surface velocity to zero */
   for (j = 1; j < jmt - 1; j++)
      for (i = 0; i < imt; i++) {
         WVEL_POS[0][j][i] = 0.0;
         WVEL_NEG[0][j][i] = 0.0;
      }

   if (dbg_lvl > 1) {
      printf ("(%d) exiting %s\n", iam, subname);
      fflush (stdout);
   }

   return 0;
}

/******************************************************************************/

void
add_UTE_coeffs_upwind3 (double ***UTE_POS, double ***UTE_NEG)
{
   char *subname = "add_UTE_coeffs_upwind3";
   int tracer_ind;
   int tracer_state_ind;

   if (dbg_lvl > 1) {
      printf ("(%d) entering %s\n", iam, subname);
      fflush (stdout);
   }

   for (tracer_ind = 0; tracer_ind < coupled_tracer_cnt; tracer_ind++) {
      for (tracer_state_ind = 0; tracer_state_ind < tracer_state_len; tracer_state_ind++) {
         int i;
         int ip1;
         int im1;
         int ip2;
         int im2;
         int j;
         int k;
         int coef_ind;

         i = tracer_state_ind_to_int3[tracer_state_ind].i;
         j = tracer_state_ind_to_int3[tracer_state_ind].j;
         k = tracer_state_ind_to_int3[tracer_state_ind].k;
         ip1 = (i < imt - 1) ? i + 1 : 0;
         im1 = (i > 0) ? i - 1 : imt - 1;
         ip2 = (ip1 < imt - 1) ? ip1 + 1 : 0;
         im2 = (im1 > 0) ? im1 - 1 : imt - 1;

         coef_ind = coef_ind_self[tracer_ind][tracer_state_ind];

         /* cell itself */
         /* advection through east face */
         if (k < KMT[j][im1])
            nzval_row_wise[coef_ind] -= 0.75 * UTE_POS[k][j][i] / TAREA[j][i] * delta_t;
         else
            nzval_row_wise[coef_ind] -= (0.75 - 0.125) * UTE_POS[k][j][i] / TAREA[j][i] * delta_t;
         nzval_row_wise[coef_ind] -= 0.375 * UTE_NEG[k][j][i] / TAREA[j][i] * delta_t;
         /* advection through west face */
         nzval_row_wise[coef_ind] += 0.375 * UTE_POS[k][j][im1] / TAREA[j][i] * delta_t;
         if (k < KMT[j][ip1])
            nzval_row_wise[coef_ind] += 0.75 * UTE_NEG[k][j][im1] / TAREA[j][i] * delta_t;
         else
            nzval_row_wise[coef_ind] += (0.75 - 0.125) * UTE_NEG[k][j][im1] / TAREA[j][i] * delta_t;
         coef_ind++;
         /* cell 1 level shallower */
         if (k - 1 >= 0)
            coef_ind++;
         /* cell 1 level deeper */
         if (k + 1 < KMT[j][i])
            coef_ind++;
         /* cell 1 unit east */
         if (k < KMT[j][ip1]) {
            /* advection through east face */
            nzval_row_wise[coef_ind] -= 0.375 * UTE_POS[k][j][i] / TAREA[j][i] * delta_t;
            if (k < KMT[j][ip2])
               nzval_row_wise[coef_ind] -= 0.75 * UTE_NEG[k][j][i] / TAREA[j][i] * delta_t;
            else
               nzval_row_wise[coef_ind] -= (0.75 - 0.125) * UTE_NEG[k][j][i] / TAREA[j][i] * delta_t;
            /* advection through west face */
            nzval_row_wise[coef_ind] += (-0.125) * UTE_NEG[k][j][im1] / TAREA[j][i] * delta_t;
            coef_ind++;
         }
         /* cell 1 unit west */
         if (k < KMT[j][im1]) {
            /* advection through east face */
            nzval_row_wise[coef_ind] -= (-0.125) * UTE_POS[k][j][i] / TAREA[j][i] * delta_t;
            /* advection through west face */
            if (k < KMT[j][im2])
               nzval_row_wise[coef_ind] += 0.75 * UTE_POS[k][j][im1] / TAREA[j][i] * delta_t;
            else
               nzval_row_wise[coef_ind] += (0.75 - 0.125) * UTE_POS[k][j][im1] / TAREA[j][i] * delta_t;
            nzval_row_wise[coef_ind] += 0.375 * UTE_NEG[k][j][im1] / TAREA[j][i] * delta_t;
            coef_ind++;
         }
         /* cell 1 unit north */
         if (k < KMT[j + 1][i])
            coef_ind++;
         /* cell 1 unit south */
         if (k < KMT[j - 1][i])
            coef_ind++;

         coef_ind = coef_ind_adv_non_nbr[tracer_ind][tracer_state_ind];

         /* cell 2 level shallower */
         if (k - 2 >= 0)
            coef_ind++;
         /* cell 2 level deeper */
         if (k + 2 < KMT[j][i])
            coef_ind++;
         /* cell 2 unit east */
         if (k < KMT[j][ip2]) {
            /* advection through east face */
            nzval_row_wise[coef_ind] -= (-0.125) * UTE_NEG[k][j][i] / TAREA[j][i] * delta_t;
            coef_ind++;
         }
         /* cell 2 unit west */
         if (k < KMT[j][im2]) {
            /* advection through west face */
            nzval_row_wise[coef_ind] += (-0.125) * UTE_POS[k][j][im1] / TAREA[j][i] * delta_t;
            coef_ind++;
         }
         /* cell 2 unit north */
         if ((j + 2 < jmt) && (k < KMT[j + 2][i]))
            coef_ind++;
         /* cell 2 unit south */
         if ((j - 2 >= 0) && (k < KMT[j - 2][i]))
            coef_ind++;
      }
   }

   if (dbg_lvl > 1) {
      printf ("(%d) exiting %s\n", iam, subname);
      fflush (stdout);
   }
}

/******************************************************************************/

void
add_VTN_coeffs_upwind3 (double ***VTN_POS, double ***VTN_NEG)
{
   char *subname = "add_VTN_coeffs_upwind3";
   int tracer_ind;
   int tracer_state_ind;

   if (dbg_lvl > 1) {
      printf ("(%d) entering %s\n", iam, subname);
      fflush (stdout);
   }

   for (tracer_ind = 0; tracer_ind < coupled_tracer_cnt; tracer_ind++) {
      for (tracer_state_ind = 0; tracer_state_ind < tracer_state_len; tracer_state_ind++) {
         int i;
         int ip1;
         int im1;
         int ip2;
         int im2;
         int j;
         int k;
         int coef_ind;

         i = tracer_state_ind_to_int3[tracer_state_ind].i;
         j = tracer_state_ind_to_int3[tracer_state_ind].j;
         k = tracer_state_ind_to_int3[tracer_state_ind].k;
         ip1 = (i < imt - 1) ? i + 1 : 0;
         im1 = (i > 0) ? i - 1 : imt - 1;
         ip2 = (ip1 < imt - 1) ? ip1 + 1 : 0;
         im2 = (im1 > 0) ? im1 - 1 : imt - 1;

         coef_ind = coef_ind_self[tracer_ind][tracer_state_ind];

         /* cell itself */
         /* advection through north face */
         if (k < KMT[j - 1][i])
            nzval_row_wise[coef_ind] -= 0.75 * VTN_POS[k][j][i] / TAREA[j][i] * delta_t;
         else
            nzval_row_wise[coef_ind] -= (0.75 - 0.125) * VTN_POS[k][j][i] / TAREA[j][i] * delta_t;
         nzval_row_wise[coef_ind] -= 0.375 * VTN_NEG[k][j][i] / TAREA[j][i] * delta_t;
         /* advection through south face */
         nzval_row_wise[coef_ind] += 0.375 * VTN_POS[k][j - 1][i] / TAREA[j][i] * delta_t;
         if (k < KMT[j + 1][i])
            nzval_row_wise[coef_ind] += 0.75 * VTN_NEG[k][j - 1][i] / TAREA[j][i] * delta_t;
         else
            nzval_row_wise[coef_ind] += (0.75 - 0.125) * VTN_NEG[k][j - 1][i] / TAREA[j][i] * delta_t;
         coef_ind++;
         /* cell 1 level shallower */
         if (k - 1 >= 0)
            coef_ind++;
         /* cell 1 level deeper */
         if (k + 1 < KMT[j][i])
            coef_ind++;
         /* cell 1 unit east */
         if (k < KMT[j][ip1])
            coef_ind++;
         /* cell 1 unit west */
         if (k < KMT[j][im1])
            coef_ind++;
         /* cell 1 unit north */
         if (k < KMT[j + 1][i]) {
            /* advection through north face */
            nzval_row_wise[coef_ind] -= 0.375 * VTN_POS[k][j][i] / TAREA[j][i] * delta_t;
            if ((j + 2 < jmt) && (k < KMT[j + 2][i]))
               nzval_row_wise[coef_ind] -= 0.75 * VTN_NEG[k][j][i] / TAREA[j][i] * delta_t;
            else
               nzval_row_wise[coef_ind] -= (0.75 - 0.125) * VTN_NEG[k][j][i] / TAREA[j][i] * delta_t;
            /* advection through south face */
            nzval_row_wise[coef_ind] += (-0.125) * VTN_NEG[k][j - 1][i] / TAREA[j][i] * delta_t;
            coef_ind++;
         }
         /* cell 1 unit south */
         if (k < KMT[j - 1][i]) {
            /* advection through north face */
            nzval_row_wise[coef_ind] -= (-0.125) * VTN_POS[k][j][i] / TAREA[j][i] * delta_t;
            /* advection through south face */
            if ((j - 2 >= 0) && (k < KMT[j - 2][i]))
               nzval_row_wise[coef_ind] += 0.75 * VTN_POS[k][j - 1][i] / TAREA[j][i] * delta_t;
            else
               nzval_row_wise[coef_ind] += (0.75 - 0.125) * VTN_POS[k][j - 1][i] / TAREA[j][i] * delta_t;
            nzval_row_wise[coef_ind] += 0.375 * VTN_NEG[k][j - 1][i] / TAREA[j][i] * delta_t;
            coef_ind++;
         }

         coef_ind = coef_ind_adv_non_nbr[tracer_ind][tracer_state_ind];

         /* cell 2 level shallower */
         if (k - 2 >= 0)
            coef_ind++;
         /* cell 2 level deeper */
         if (k + 2 < KMT[j][i])
            coef_ind++;
         /* cell 2 unit east */
         if (k < KMT[j][ip2])
            coef_ind++;
         /* cell 2 unit west */
         if (k < KMT[j][im2])
            coef_ind++;
         /* cell 2 unit north */
         if ((j + 2 < jmt) && (k < KMT[j + 2][i])) {
            /* advection through north face */
            nzval_row_wise[coef_ind] -= (-0.125) * VTN_NEG[k][j][i] / TAREA[j][i] * delta_t;
            coef_ind++;
         }
         /* cell 2 unit south */
         if ((j - 2 >= 0) && (k < KMT[j - 2][i])) {
            /* advection through south face */
            nzval_row_wise[coef_ind] += (-0.125) * VTN_POS[k][j - 1][i] / TAREA[j][i] * delta_t;
            coef_ind++;
         }
      }
   }

   if (dbg_lvl > 1) {
      printf ("(%d) exiting %s\n", iam, subname);
      fflush (stdout);
   }
}

/******************************************************************************/

int
add_WVEL_coeffs_upwind3 (double ***WVEL_POS, double ***WVEL_NEG)
{
   char *subname = "add_WVEL_coeffs_upwind3";
   int tracer_ind;
   int tracer_state_ind;
   double *dzc_tmp;
   double *dzc;
   double *talfzp;
   double *tbetzp;
   double *tgamzp;
   double *talfzm;
   double *tbetzm;
   double *tdelzm;
   int k;

   if (dbg_lvl > 1) {
      printf ("(%d) entering %s\n", iam, subname);
      fflush (stdout);
   }

   if ((dzc_tmp = malloc ((size_t) (km + 2) * sizeof (double))) == NULL) {
      fprintf (stderr, "(%d) malloc failed in %s for dzc\n", iam, subname);
      return 1;
   }
   dzc = dzc_tmp + 1;
   if ((talfzp = malloc ((size_t) km * sizeof (double))) == NULL) {
      fprintf (stderr, "(%d) malloc failed in %s for talfzp\n", iam, subname);
      return 1;
   }
   if ((tbetzp = malloc ((size_t) km * sizeof (double))) == NULL) {
      fprintf (stderr, "(%d) malloc failed in %s for tbetzp\n", iam, subname);
      return 1;
   }
   if ((tgamzp = malloc ((size_t) km * sizeof (double))) == NULL) {
      fprintf (stderr, "(%d) malloc failed in %s for tgamzp\n", iam, subname);
      return 1;
   }
   if ((talfzm = malloc ((size_t) km * sizeof (double))) == NULL) {
      fprintf (stderr, "(%d) malloc failed in %s for talfzm\n", iam, subname);
      return 1;
   }
   if ((tbetzm = malloc ((size_t) km * sizeof (double))) == NULL) {
      fprintf (stderr, "(%d) malloc failed in %s for tbetzm\n", iam, subname);
      return 1;
   }
   if ((tdelzm = malloc ((size_t) km * sizeof (double))) == NULL) {
      fprintf (stderr, "(%d) malloc failed in %s for tdelzm\n", iam, subname);
      return 1;
   }

   /* k indices shifted wrt index in POP, but relative */
   /* indices between dz and dzc is the same           */
   /*                                                  */
   /*       index range   index range                  */
   /* VAR     in POP         here                      */
   /*  dz      1:km          0:km-1                    */
   /* dzc      0:km+1       -1:km                      */

   dzc[-1] = dz[0];
   for (k = 0; k < km; k++)
      dzc[k] = dz[k];
   dzc[km] = dzc[km - 1];

   for (k = 0; k < km - 1; k++) {
      talfzp[k] = dz[k] * (2.0 * dz[k] + dzc[k - 1]) / (dz[k] + dz[k + 1]) / (dzc[k - 1] + 2.0 * dz[k] + dz[k + 1]);
      tbetzp[k] = dz[k + 1] * (2.0 * dz[k] + dzc[k - 1]) / (dz[k] + dz[k + 1]) / (dz[k] + dzc[k - 1]);
      tgamzp[k] = -(dz[k] * dz[k + 1]) / (dz[k] + dzc[k - 1]) / (dz[k + 1] + dzc[k - 1] + 2.0 * dz[k]);
   }
   tbetzp[0] = tbetzp[0] + tgamzp[0];
   tgamzp[0] = 0.0;
   talfzp[km - 1] = 0.0;
   tbetzp[km - 1] = 0.0;
   tgamzp[km - 1] = 0.0;

   for (k = 0; k < km - 1; k++) {
      talfzm[k] = dz[k] * (2.0 * dz[k + 1] + dzc[k + 2]) / (dz[k] + dz[k + 1]) / (dz[k + 1] + dzc[k + 2]);
      tbetzm[k] = dz[k + 1] * (2.0 * dz[k + 1] + dzc[k + 2]) / (dz[k] + dz[k + 1]) / (dz[k] + dzc[k + 2] + 2.0 * dz[k + 1]);
      tdelzm[k] = -(dz[k] * dz[k + 1]) / (dz[k + 1] + dzc[k + 2]) / (dz[k] + dzc[k + 2] + 2.0 * dz[k + 1]);
   }
   talfzm[km - 1] = 0.0;
   tbetzm[km - 1] = 0.0;
   tdelzm[km - 1] = 0.0;

   for (tracer_ind = 0; tracer_ind < coupled_tracer_cnt; tracer_ind++) {
      for (tracer_state_ind = 0; tracer_state_ind < tracer_state_len; tracer_state_ind++) {
         int i;
         int ip1;
         int im1;
         int ip2;
         int im2;
         int j;
         int coef_ind;

         i = tracer_state_ind_to_int3[tracer_state_ind].i;
         j = tracer_state_ind_to_int3[tracer_state_ind].j;
         k = tracer_state_ind_to_int3[tracer_state_ind].k;
         ip1 = (i < imt - 1) ? i + 1 : 0;
         im1 = (i > 0) ? i - 1 : imt - 1;
         ip2 = (ip1 < imt - 1) ? ip1 + 1 : 0;
         im2 = (im1 > 0) ? im1 - 1 : imt - 1;

         coef_ind = coef_ind_self[tracer_ind][tracer_state_ind];

         /* cell itself */
         /* advection through top face */
         if (k - 1 >= 0) {
            if (k + 1 < KMT[j][i])
               nzval_row_wise[coef_ind] -= talfzm[k - 1] * WVEL_POS[k][j][i] / dz[k] * delta_t;
            else
               nzval_row_wise[coef_ind] -= (talfzm[k - 1] + tdelzm[k - 1]) * WVEL_POS[k][j][i] / dz[k] * delta_t;
            nzval_row_wise[coef_ind] -= talfzp[k - 1] * WVEL_NEG[k][j][i] / dz[k] * delta_t;
         }
         /* advection through bottom face */
         if (k + 1 < KMT[j][i]) {
            nzval_row_wise[coef_ind] += tbetzm[k] * WVEL_POS[k + 1][j][i] / dz[k] * delta_t;
            nzval_row_wise[coef_ind] += tbetzp[k] * WVEL_NEG[k + 1][j][i] / dz[k] * delta_t;
         }
         coef_ind++;
         /* cell 1 level shallower */
         if (k - 1 >= 0) {
            /* advection through top face */
            nzval_row_wise[coef_ind] -= tbetzm[k - 1] * WVEL_POS[k][j][i] / dz[k] * delta_t;
            nzval_row_wise[coef_ind] -= tbetzp[k - 1] * WVEL_NEG[k][j][i] / dz[k] * delta_t;
            /* advection through bottom face */
            if (k + 1 < KMT[j][i])
               nzval_row_wise[coef_ind] += tgamzp[k] * WVEL_NEG[k + 1][j][i] / dz[k] * delta_t;
            coef_ind++;
         }
         /* cell 1 level deeper */
         if (k + 1 < KMT[j][i]) {
            /* advection through top face */
            if (k - 1 >= 0)
               nzval_row_wise[coef_ind] -= tdelzm[k - 1] * WVEL_POS[k][j][i] / dz[k] * delta_t;
            /* advection through bottom face */
            if (k + 2 < KMT[j][i])
               nzval_row_wise[coef_ind] += talfzm[k] * WVEL_POS[k + 1][j][i] / dz[k] * delta_t;
            else
               nzval_row_wise[coef_ind] += (talfzm[k] + tdelzm[k]) * WVEL_POS[k + 1][j][i] / dz[k] * delta_t;
            nzval_row_wise[coef_ind] += talfzp[k] * WVEL_NEG[k + 1][j][i] / dz[k] * delta_t;
            coef_ind++;
         }
         /* cell 1 unit east */
         if (k < KMT[j][ip1])
            coef_ind++;
         /* cell 1 unit west */
         if (k < KMT[j][im1])
            coef_ind++;
         /* cell 1 unit north */
         if (k < KMT[j + 1][i])
            coef_ind++;
         /* cell 1 unit south */
         if (k < KMT[j - 1][i])
            coef_ind++;

         coef_ind = coef_ind_adv_non_nbr[tracer_ind][tracer_state_ind];

         /* cell 2 level shallower */
         if (k - 2 >= 0) {
            /* advection through top face */
            nzval_row_wise[coef_ind] -= tgamzp[k - 1] * WVEL_NEG[k][j][i] / dz[k] * delta_t;
            coef_ind++;
         }
         /* cell 2 level deeper */
         if (k + 2 < KMT[j][i]) {
            /* advection through bottom face */
            nzval_row_wise[coef_ind] += tdelzm[k] * WVEL_POS[k + 1][j][i] / dz[k] * delta_t;
            coef_ind++;
         }
         /* cell 2 unit east */
         if (k < KMT[j][ip2])
            coef_ind++;
         /* cell 2 unit west */
         if (k < KMT[j][im2])
            coef_ind++;
         /* cell 2 unit north */
         if ((j + 2 < jmt) && (k < KMT[j + 2][i]))
            coef_ind++;
         /* cell 2 unit south */
         if ((j - 2 >= 0) && (k < KMT[j - 2][i]))
            coef_ind++;
      }
   }
   free (dzc_tmp);
   free (talfzp);
   free (tbetzp);
   free (tgamzp);
   free (talfzm);
   free (tbetzm);
   free (tdelzm);

   if (dbg_lvl > 1) {
      printf ("(%d) exiting %s\n", iam, subname);
      fflush (stdout);
   }

   return 0;
}

/******************************************************************************/

int
add_adv (void)
{
   char *subname = "add_adv";

   double ***VEL_WIDTH;
   double ***VEL_WIDTH_POS;
   double ***VEL_WIDTH_NEG;

   if (dbg_lvl > 1) {
      printf ("(%d) entering %s\n", iam, subname);
      fflush (stdout);
   }

   switch (adv_opt) {
   case adv_none:
      break;
   case adv_donor:
   case adv_cent:
      if ((VEL_WIDTH = malloc_3d_double (km, jmt, imt)) == NULL) {
         fprintf (stderr, "(%d) malloc failed in %s for VEL_WIDTH\n", iam, subname);
         return 1;
      }
      if (load_UTE (VEL_WIDTH))
         return 1;
      add_UTE_coeffs (VEL_WIDTH);
      if (load_VTN (VEL_WIDTH))
         return 1;
      add_VTN_coeffs (VEL_WIDTH);
      if (load_WVEL (VEL_WIDTH))
         return 1;
      add_WVEL_coeffs (VEL_WIDTH);
      free_3d_double (VEL_WIDTH);
      break;
   case adv_upwind3:
      if ((VEL_WIDTH_POS = malloc_3d_double (km, jmt, imt)) == NULL) {
         fprintf (stderr, "(%d) malloc failed in %s for VEL_WIDTH_POS\n", iam, subname);
         return 1;
      }
      if ((VEL_WIDTH_NEG = malloc_3d_double (km, jmt, imt)) == NULL) {
         fprintf (stderr, "(%d) malloc failed in %s for VEL_WIDTH_NEG\n", iam, subname);
         return 1;
      }
      if (load_UTE_upwind3 (VEL_WIDTH_POS, VEL_WIDTH_NEG))
         return 1;
      add_UTE_coeffs_upwind3 (VEL_WIDTH_POS, VEL_WIDTH_NEG);
      if (load_VTN_upwind3 (VEL_WIDTH_POS, VEL_WIDTH_NEG))
         return 1;
      add_VTN_coeffs_upwind3 (VEL_WIDTH_POS, VEL_WIDTH_NEG);
      if (load_WVEL_upwind3 (VEL_WIDTH_POS, VEL_WIDTH_NEG))
         return 1;
      if (add_WVEL_coeffs_upwind3 (VEL_WIDTH_POS, VEL_WIDTH_NEG))
         return 1;
      free_3d_double (VEL_WIDTH_POS);
      free_3d_double (VEL_WIDTH_NEG);
      break;
   }

   if (dbg_lvl > 1) {
      printf ("(%d) adv terms added\n\n", iam);
      fflush (stdout);
   }

   if (dbg_lvl > 1) {
      printf ("(%d) exiting %s\n", iam, subname);
      fflush (stdout);
   }

   return 0;
}

/******************************************************************************/

int
adv_enforce_divfree (void)
{
   char *subname = "adv_enforce_divfree";
   int tracer_ind;
   int tracer_state_ind;

   if (dbg_lvl > 1) {
      printf ("(%d) entering %s\n", iam, subname);
      fflush (stdout);
   }

   for (tracer_ind = 0; tracer_ind < coupled_tracer_cnt; tracer_ind++) {
      for (tracer_state_ind = 0; tracer_state_ind < tracer_state_len; tracer_state_ind++) {
         double nzval_sum_non_self;
         int i;
         int ip1;
         int im1;
         int ip2;
         int im2;
         int j;
         int k;
         int coef_ind;

         i = tracer_state_ind_to_int3[tracer_state_ind].i;
         j = tracer_state_ind_to_int3[tracer_state_ind].j;
         k = tracer_state_ind_to_int3[tracer_state_ind].k;
         ip1 = (i < imt - 1) ? i + 1 : 0;
         im1 = (i > 0) ? i - 1 : imt - 1;
         ip2 = (ip1 < imt - 1) ? ip1 + 1 : 0;
         im2 = (im1 > 0) ? im1 - 1 : imt - 1;

         coef_ind = coef_ind_self[tracer_ind][tracer_state_ind];

         /* cell itself */
         nzval_sum_non_self = 0.0;
         coef_ind++;
         /* cell 1 level shallower */
         if (k - 1 >= 0) {
            nzval_sum_non_self += nzval_row_wise[coef_ind];
            coef_ind++;
         }
         /* cell 1 level deeper */
         if (k + 1 < KMT[j][i]) {
            nzval_sum_non_self += nzval_row_wise[coef_ind];
            coef_ind++;
         }
         /* cell 1 unit east */
         if (k < KMT[j][ip1]) {
            nzval_sum_non_self += nzval_row_wise[coef_ind];
            coef_ind++;
         }
         /* cell 1 unit west */
         if (k < KMT[j][im1]) {
            nzval_sum_non_self += nzval_row_wise[coef_ind];
            coef_ind++;
         }
         /* cell 1 unit north */
         if (k < KMT[j + 1][i]) {
            nzval_sum_non_self += nzval_row_wise[coef_ind];
            coef_ind++;
         }
         /* cell 1 unit south */
         if (k < KMT[j - 1][i]) {
            nzval_sum_non_self += nzval_row_wise[coef_ind];
            coef_ind++;
         }

         coef_ind = coef_ind_adv_non_nbr[tracer_ind][tracer_state_ind];
         if (adv_opt == adv_upwind3) {
            /* cell 2 level shallower */
            if (k - 2 >= 0) {
               nzval_sum_non_self += nzval_row_wise[coef_ind];
               coef_ind++;
            }
            /* cell 2 level deeper */
            if (k + 2 < KMT[j][i]) {
               nzval_sum_non_self += nzval_row_wise[coef_ind];
               coef_ind++;
            }
            /* cell 2 unit east */
            if (k < KMT[j][ip2]) {
               nzval_sum_non_self += nzval_row_wise[coef_ind];
               coef_ind++;
            }
            /* cell 2 unit west */
            if (k < KMT[j][im2]) {
               nzval_sum_non_self += nzval_row_wise[coef_ind];
               coef_ind++;
            }
            /* cell 2 unit north */
            if ((j + 2 < jmt) && (k < KMT[j + 2][i])) {
               nzval_sum_non_self += nzval_row_wise[coef_ind];
               coef_ind++;
            }
            /* cell 2 unit south */
            if ((j - 2 >= 0) && (k < KMT[j - 2][i])) {
               nzval_sum_non_self += nzval_row_wise[coef_ind];
               coef_ind++;
            }
         }

         coef_ind = coef_ind_self[tracer_ind][tracer_state_ind];
         nzval_sum_non_self = -nzval_sum_non_self;
      }
   }

   if (dbg_lvl > 1) {
      printf ("(%d) exiting %s\n", iam, subname);
      fflush (stdout);
   }

   return 0;
}

/******************************************************************************/

int
add_hmix_isop_file (void)
{
   char *subname = "add_hmix_isop_file";
   int iprime;
   int jprime;
   int kprime;
   int var_exists;
   double ***IRF;
   char IRF_name[64];

   int tracer_ind;
   int tracer_state_ind;

   if (dbg_lvl > 1) {
      printf ("(%d) entering %s\n", iam, subname);
      fflush (stdout);
   }

   if ((IRF = malloc_3d_double (km, jmt, imt)) == NULL) {
      fprintf (stderr, "(%d) malloc failed in %s for IRF\n", iam, subname);
      return 1;
   }
   for (iprime = 0; iprime < 4; iprime++)
      for (jprime = 0; jprime < 3; jprime++)
         for (kprime = 0; kprime < 3; kprime++) {
            sprintf (IRF_name, "HDIF_EXPLICIT_3D_IRF_%d_%d_%d", iprime + 1, jprime + 1, kprime + 1);
            if (var_exists_in_file (circ_fname, IRF_name, &var_exists)) {
               fprintf (stderr, "(%d) var_exists_in_file failed in %s for field_name %s in file %s\n", iam, subname,
                        IRF_name, circ_fname);
               return 1;
            }
            if (!var_exists) {
               if (dbg_lvl)
                  printf ("(%d) %s: %s not found in %s\n", iam, subname, IRF_name, circ_fname);
               sprintf (IRF_name, "HDIF_EXPLICIT_3D_IRF_NK_%d_%d_%d", iprime + 1, jprime + 1, kprime + 1);
               if (var_exists_in_file (circ_fname, IRF_name, &var_exists)) {
                  fprintf (stderr, "(%d) var_exists_in_file failed in %s for field_name %s in file %s\n", iam, subname,
                           IRF_name, circ_fname);
                  return 1;
               }
               if (!var_exists) {
                  if (dbg_lvl)
                     printf ("(%d) %s: %s not found in %s\n", iam, subname, IRF_name, circ_fname);
                  return 1;
               }
            }
            if (dbg_lvl)
               printf ("(%d) %s: reading %s from %s\n", iam, subname, IRF_name, circ_fname);
            if (get_var_3d_double (circ_fname, IRF_name, IRF))
               return 1;

            for (tracer_ind = 0; tracer_ind < coupled_tracer_cnt; tracer_ind++) {
               for (tracer_state_ind = 0; tracer_state_ind < tracer_state_len; tracer_state_ind++) {
                  int i;
                  int ip1;
                  int im1;
                  int j;
                  int k;
                  int coef_ind;

                  i = tracer_state_ind_to_int3[tracer_state_ind].i;
                  j = tracer_state_ind_to_int3[tracer_state_ind].j;
                  k = tracer_state_ind_to_int3[tracer_state_ind].k;
                  ip1 = (i < imt - 1) ? i + 1 : 0;
                  im1 = (i > 0) ? i - 1 : imt - 1;

                  coef_ind = coef_ind_self[tracer_ind][tracer_state_ind];

                  /* cell itself */
                  if ((i % 4 == iprime) && (j % 3 == jprime) && (k % 3 == kprime))
                     nzval_row_wise[coef_ind] += IRF[k][j][i] * delta_t;
                  coef_ind++;
                  /* cell 1 level shallower */
                  if (k - 1 >= 0) {
                     if ((i % 4 == iprime) && (j % 3 == jprime) && ((k - 1) % 3 == kprime))
                        nzval_row_wise[coef_ind] += IRF[k][j][i] * delta_t;
                     coef_ind++;
                  }
                  /* cell 1 level deeper */
                  if (k + 1 < KMT[j][i]) {
                     if ((i % 4 == iprime) && (j % 3 == jprime) && ((k + 1) % 3 == kprime))
                        nzval_row_wise[coef_ind] += IRF[k][j][i] * delta_t;
                     coef_ind++;
                  }
                  /* cell 1 unit east */
                  if (k < KMT[j][ip1]) {
                     if ((ip1 % 4 == iprime) && (j % 3 == jprime) && (k % 3 == kprime))
                        nzval_row_wise[coef_ind] += IRF[k][j][i] * delta_t;
                     coef_ind++;
                  }
                  /* cell 1 unit west */
                  if (k < KMT[j][im1]) {
                     if ((im1 % 4 == iprime) && (j % 3 == jprime) && (k % 3 == kprime))
                        nzval_row_wise[coef_ind] += IRF[k][j][i] * delta_t;
                     coef_ind++;
                  }
                  /* cell 1 unit north */
                  if (k < KMT[j + 1][i]) {
                     if ((i % 4 == iprime) && ((j + 1) % 3 == jprime) && (k % 3 == kprime))
                        nzval_row_wise[coef_ind] += IRF[k][j][i] * delta_t;
                     coef_ind++;
                  }
                  /* cell 1 unit south */
                  if (k < KMT[j - 1][i]) {
                     if ((i % 4 == iprime) && ((j - 1) % 3 == jprime) && (k % 3 == kprime))
                        nzval_row_wise[coef_ind] += IRF[k][j][i] * delta_t;
                     coef_ind++;
                  }

                  coef_ind = coef_ind_hmix_non_nbr[tracer_ind][tracer_state_ind];

                  /* shallower & east */
                  if ((k - 1 >= 0) && (k - 1 < KMT[j][ip1])) {
                     if ((ip1 % 4 == iprime) && (j % 3 == jprime) && ((k - 1) % 3 == kprime))
                        nzval_row_wise[coef_ind] += IRF[k][j][i] * delta_t;
                     coef_ind++;
                  }
                  /* deeper & east */
                  if (k + 1 < KMT[j][ip1]) {
                     if ((ip1 % 4 == iprime) && (j % 3 == jprime) && ((k + 1) % 3 == kprime))
                        nzval_row_wise[coef_ind] += IRF[k][j][i] * delta_t;
                     coef_ind++;
                  }
                  /* shallower & west */
                  if ((k - 1 >= 0) && (k - 1 < KMT[j][im1])) {
                     if ((im1 % 4 == iprime) && (j % 3 == jprime) && ((k - 1) % 3 == kprime))
                        nzval_row_wise[coef_ind] += IRF[k][j][i] * delta_t;
                     coef_ind++;
                  }
                  /* deeper & west */
                  if (k + 1 < KMT[j][im1]) {
                     if ((im1 % 4 == iprime) && (j % 3 == jprime) && ((k + 1) % 3 == kprime))
                        nzval_row_wise[coef_ind] += IRF[k][j][i] * delta_t;
                     coef_ind++;
                  }
                  /* shallower & north */
                  if ((k - 1 >= 0) && (k - 1 < KMT[j + 1][i])) {
                     if ((i % 4 == iprime) && ((j + 1) % 3 == jprime)
                         && ((k - 1) % 3 == kprime))
                        nzval_row_wise[coef_ind] += IRF[k][j][i] * delta_t;
                     coef_ind++;
                  }
                  /* deeper & north */
                  if (k + 1 < KMT[j + 1][i]) {
                     if ((i % 4 == iprime) && ((j + 1) % 3 == jprime)
                         && ((k + 1) % 3 == kprime))
                        nzval_row_wise[coef_ind] += IRF[k][j][i] * delta_t;
                     coef_ind++;
                  }
                  /* shallower & south */
                  if ((k - 1 >= 0) && (k - 1 < KMT[j - 1][i])) {
                     if ((i % 4 == iprime) && ((j - 1) % 3 == jprime)
                         && ((k - 1) % 3 == kprime))
                        nzval_row_wise[coef_ind] += IRF[k][j][i] * delta_t;
                     coef_ind++;
                  }
                  /* deeper & south */
                  if (k + 1 < KMT[j - 1][i]) {
                     if ((i % 4 == iprime) && ((j - 1) % 3 == jprime)
                         && ((k + 1) % 3 == kprime))
                        nzval_row_wise[coef_ind] += IRF[k][j][i] * delta_t;
                     coef_ind++;
                  }
               }
            }
         }

   free_3d_double (IRF);

   if (dbg_lvl > 1) {
      printf ("(%d) exiting %s\n", iam, subname);
      fflush (stdout);
   }

   return 0;
}

/******************************************************************************/

int
add_hmix_hor_file (void)
{
   char *subname = "add_hmix_hor_file";
   double ***WORK;
   double ***KAPPA;
   double **HUS;
   double **HTE;
   double **HUW;
   double **HTN;
   double fv;

   int tracer_ind;
   int tracer_state_ind;

   if (dbg_lvl > 1) {
      printf ("(%d) entering %s\n", iam, subname);
      fflush (stdout);
   }

   if ((KAPPA = malloc_3d_double (km, jmt, imt)) == NULL) {
      fprintf (stderr, "(%d) malloc failed in %s for KAPPA\n", iam, subname);
      return 1;
   }
   if ((WORK = malloc_3d_double (km, jmt, imt)) == NULL) {
      fprintf (stderr, "(%d) malloc failed in %s for WORK\n", iam, subname);
      return 1;
   }
   if (dbg_lvl)
      printf ("(%d) %s: reading KAPPA_ISOP,HOR_DIFF from %s\n", iam, subname, circ_fname);
   if (get_var_3d_double (circ_fname, "KAPPA_ISOP", KAPPA))
      return 1;
   if (get_att_double (circ_fname, "KAPPA_ISOP", "_FillValue", &fv))
      return 1;
   set_fv_3d_double (fv, 0.0, KAPPA);
   if (get_var_3d_double (circ_fname, "HOR_DIFF", WORK))
      return 1;
   if (get_att_double (circ_fname, "HOR_DIFF", "_FillValue", &fv))
      return 1;
   set_fv_3d_double (fv, 0.0, WORK);
   {
      int i;
      int j;
      int k;
      for (k = 0; k < km; k++)
         for (j = 1; j < jmt - 1; j++)
            for (i = 0; i < imt; i++)
               if (k < KMT[j][i])
                  KAPPA[k][j][i] += WORK[k][j][i];
   }
   free_3d_double (WORK);

   if ((HUS = malloc_2d_double (jmt, imt)) == NULL) {
      fprintf (stderr, "(%d) malloc failed in %s for HUS\n", iam, subname);
      return 1;
   }
   if ((HTE = malloc_2d_double (jmt, imt)) == NULL) {
      fprintf (stderr, "(%d) malloc failed in %s for HTE\n", iam, subname);
      return 1;
   }
   if ((HUW = malloc_2d_double (jmt, imt)) == NULL) {
      fprintf (stderr, "(%d) malloc failed in %s for HUW\n", iam, subname);
      return 1;
   }
   if ((HTN = malloc_2d_double (jmt, imt)) == NULL) {
      fprintf (stderr, "(%d) malloc failed in %s for HTN\n", iam, subname);
      return 1;
   }
   if (dbg_lvl)
      printf ("(%d) %s: reading HUS,HTE,HUW,HTN from %s\n", iam, subname, circ_fname);
   if (get_var_2d_double (circ_fname, "HUS", HUS))
      return 1;
   if (get_att_double (circ_fname, "HUS", "_FillValue", &fv))
      return 1;
   set_fv_2d_double (fv, 0.0, HUS);
   if (get_var_2d_double (circ_fname, "HTE", HTE))
      return 1;
   if (get_att_double (circ_fname, "HTE", "_FillValue", &fv))
      return 1;
   set_fv_2d_double (fv, 0.0, HTE);
   if (get_var_2d_double (circ_fname, "HUW", HUW))
      return 1;
   if (get_att_double (circ_fname, "HUW", "_FillValue", &fv))
      return 1;
   set_fv_2d_double (fv, 0.0, HUW);
   if (get_var_2d_double (circ_fname, "HTN", HTN))
      return 1;
   if (get_att_double (circ_fname, "HTN", "_FillValue", &fv))
      return 1;
   set_fv_2d_double (fv, 0.0, HTN);

   for (tracer_ind = 0; tracer_ind < coupled_tracer_cnt; tracer_ind++) {
      for (tracer_state_ind = 0; tracer_state_ind < tracer_state_len; tracer_state_ind++) {
         int i;
         int ip1;
         int im1;
         int j;
         int k;
         int coef_ind;

         double ce;
         double cw;
         double cn;
         double cs;

         i = tracer_state_ind_to_int3[tracer_state_ind].i;
         j = tracer_state_ind_to_int3[tracer_state_ind].j;
         k = tracer_state_ind_to_int3[tracer_state_ind].k;
         ip1 = (i < imt - 1) ? i + 1 : 0;
         im1 = (i > 0) ? i - 1 : imt - 1;

         /* cell 1 unit east */
         if (k < KMT[j][ip1])
            ce = 0.5 * (KAPPA[k][j][i] + KAPPA[k][j][ip1]) * HTE[j][i] / HUS[j][i] / TAREA[j][i] * delta_t;
         else
            ce = 0.0;

         /* cell 1 unit west */
         if (k < KMT[j][im1])
            cw = 0.5 * (KAPPA[k][j][im1] + KAPPA[k][j][i]) * HTE[j][im1] / HUS[j][im1] / TAREA[j][i] * delta_t;
         else
            cw = 0.0;

         /* cell 1 unit north */
         if (k < KMT[j + 1][i])
            cn = 0.5 * (KAPPA[k][j][i] + KAPPA[k][j + 1][i]) * HTN[j][i] / HUW[j][i] / TAREA[j][i] * delta_t;
         else
            cn = 0.0;

         /* cell 1 unit south */
         if (k < KMT[j - 1][i])
            cs = 0.5 * (KAPPA[k][j - 1][i] + KAPPA[k][j][i]) * HTN[j - 1][i] / HUW[j - 1]
               [i] / TAREA[j][i] * delta_t;
         else
            cs = 0.0;

         coef_ind = coef_ind_self[tracer_ind][tracer_state_ind];

         /* cell itself */
         nzval_row_wise[coef_ind] -= (ce + cw + cn + cs);
         coef_ind++;
         /* cell 1 level shallower */
         if (k - 1 >= 0)
            coef_ind++;
         /* cell 1 level deeper */
         if (k + 1 < KMT[j][i])
            coef_ind++;
         /* cell 1 unit east */
         if (k < KMT[j][ip1]) {
            nzval_row_wise[coef_ind] += ce;
            coef_ind++;
         }
         /* cell 1 unit west */
         if (k < KMT[j][im1]) {
            nzval_row_wise[coef_ind] += cw;
            coef_ind++;
         }
         /* cell 1 unit north */
         if (k < KMT[j + 1][i]) {
            nzval_row_wise[coef_ind] += cn;
            coef_ind++;
         }
         /* cell 1 unit south */
         if (k < KMT[j - 1][i]) {
            nzval_row_wise[coef_ind] += cs;
            coef_ind++;
         }
      }
   }

   free_2d_double (HTN);
   free_2d_double (HUW);
   free_2d_double (HTE);
   free_2d_double (HUS);
   free_3d_double (KAPPA);

   if (dbg_lvl > 1) {
      printf ("(%d) exiting %s\n", iam, subname);
      fflush (stdout);
   }

   return 0;
}

/******************************************************************************/

int
add_hmix_const (void)
{
   char *subname = "add_hmix_const";
   double **HUS;
   double **HTE;
   double **HUW;
   double **HTN;
   double fv;

   int tracer_ind;
   int tracer_state_ind;
   double ah;

   if (dbg_lvl > 1) {
      printf ("(%d) entering %s\n", iam, subname);
      fflush (stdout);
   }

   ah = 4.0e6;

   if ((HUS = malloc_2d_double (jmt, imt)) == NULL) {
      fprintf (stderr, "(%d) malloc failed in %s for HUS\n", iam, subname);
      return 1;
   }
   if ((HTE = malloc_2d_double (jmt, imt)) == NULL) {
      fprintf (stderr, "(%d) malloc failed in %s for HTE\n", iam, subname);
      return 1;
   }
   if ((HUW = malloc_2d_double (jmt, imt)) == NULL) {
      fprintf (stderr, "(%d) malloc failed in %s for HUW\n", iam, subname);
      return 1;
   }
   if ((HTN = malloc_2d_double (jmt, imt)) == NULL) {
      fprintf (stderr, "(%d) malloc failed in %s for HTN\n", iam, subname);
      return 1;
   }
   if (dbg_lvl)
      printf ("(%d) %s: reading HUS,HTE,HUW,HTN from %s\n", iam, subname, circ_fname);
   if (get_var_2d_double (circ_fname, "HUS", HUS))
      return 1;
   if (get_att_double (circ_fname, "HUS", "_FillValue", &fv))
      return 1;
   set_fv_2d_double (fv, 0.0, HUS);
   if (get_var_2d_double (circ_fname, "HTE", HTE))
      return 1;
   if (get_att_double (circ_fname, "HTE", "_FillValue", &fv))
      return 1;
   set_fv_2d_double (fv, 0.0, HTE);
   if (get_var_2d_double (circ_fname, "HUW", HUW))
      return 1;
   if (get_att_double (circ_fname, "HUW", "_FillValue", &fv))
      return 1;
   set_fv_2d_double (fv, 0.0, HUW);
   if (get_var_2d_double (circ_fname, "HTN", HTN))
      return 1;
   if (get_att_double (circ_fname, "HTN", "_FillValue", &fv))
      return 1;
   set_fv_2d_double (fv, 0.0, HTN);

   for (tracer_ind = 0; tracer_ind < coupled_tracer_cnt; tracer_ind++) {
      for (tracer_state_ind = 0; tracer_state_ind < tracer_state_len; tracer_state_ind++) {
         int i;
         int ip1;
         int im1;
         int j;
         int k;
         int coef_ind;

         double ce;
         double cw;
         double cn;
         double cs;

         i = tracer_state_ind_to_int3[tracer_state_ind].i;
         j = tracer_state_ind_to_int3[tracer_state_ind].j;
         k = tracer_state_ind_to_int3[tracer_state_ind].k;
         ip1 = (i < imt - 1) ? i + 1 : 0;
         im1 = (i > 0) ? i - 1 : imt - 1;

         /* cell 1 unit east */
         if (k < KMT[j][ip1])
            ce = ah * HTE[j][i] / HUS[j][i] / TAREA[j][i] * delta_t;
         else
            ce = 0.0;

         /* cell 1 unit west */
         if (k < KMT[j][im1])
            cw = ah * HTE[j][im1] / HUS[j][im1] / TAREA[j][i] * delta_t;
         else
            cw = 0.0;

         /* cell 1 unit north */
         if (k < KMT[j + 1][i])
            cn = ah * HTN[j][i] / HUW[j][i] / TAREA[j][i] * delta_t;
         else
            cn = 0.0;

         /* cell 1 unit south */
         if (k < KMT[j - 1][i])
            cs = ah * HTN[j - 1][i] / HUW[j - 1][i] / TAREA[j][i] * delta_t;
         else
            cs = 0.0;

         coef_ind = coef_ind_self[tracer_ind][tracer_state_ind];

         /* cell itself */
         nzval_row_wise[coef_ind] -= (ce + cw + cn + cs);
         coef_ind++;
         /* cell 1 level shallower */
         if (k - 1 >= 0)
            coef_ind++;
         /* cell 1 level deeper */
         if (k + 1 < KMT[j][i])
            coef_ind++;
         /* cell 1 unit east */
         if (k < KMT[j][ip1]) {
            nzval_row_wise[coef_ind] += ce;
            coef_ind++;
         }
         /* cell 1 unit west */
         if (k < KMT[j][im1]) {
            nzval_row_wise[coef_ind] += cw;
            coef_ind++;
         }
         /* cell 1 unit north */
         if (k < KMT[j + 1][i]) {
            nzval_row_wise[coef_ind] += cn;
            coef_ind++;
         }
         /* cell 1 unit south */
         if (k < KMT[j - 1][i]) {
            nzval_row_wise[coef_ind] += cs;
            coef_ind++;
         }
      }
   }

   free_2d_double (HTN);
   free_2d_double (HUW);
   free_2d_double (HTE);
   free_2d_double (HUS);

   if (dbg_lvl > 1) {
      printf ("(%d) exiting %s\n", iam, subname);
      fflush (stdout);
   }

   return 0;
}

/******************************************************************************/

int
add_hmix (void)
{
   char *subname = "add_hmix";

   if (dbg_lvl > 1) {
      printf ("(%d) entering %s\n", iam, subname);
      fflush (stdout);
   }

   switch (hmix_opt) {
   case hmix_none:
      break;
   case hmix_const:
      if (add_hmix_const ())
         return 1;
      break;
   case hmix_hor_file:
      if (adv_opt == adv_upwind3) {
         fprintf (stderr, "(%d) cannot use hmix_hor_file with adv_upwind3\n", iam);
         return 1;
      }
      if (add_hmix_hor_file ())
         return 1;
      break;
   case hmix_isop_file:
      if (add_hmix_isop_file ())
         return 1;
      break;
   }

   if (dbg_lvl > 1) {
      printf ("(%d) hmix terms added\n\n", iam);
      fflush (stdout);
   }

   if (dbg_lvl > 1) {
      printf ("(%d) exiting %s\n", iam, subname);
      fflush (stdout);
   }

   return 0;
}

/******************************************************************************/

int
add_vmix_matrix_file (void)
{
   char *subname = "add_vmix_matrix_file";
   char varname[64];
   double ***vmix_matrix_var;

   int kprime;
   int tracer_ind;
   int tracer_state_ind;

   if (dbg_lvl > 1) {
      printf ("(%d) entering %s\n", iam, subname);
      fflush (stdout);
   }

   if ((vmix_matrix_var = malloc_3d_double (km, jmt, imt)) == NULL) {
      fprintf (stderr, "(%d) malloc failed in %s for vmix_matrix_var\n", iam, subname);
      return 1;
   }

   if (dbg_lvl)
      printf ("(%d) %s: reading vmix_matrix vars from %s\n", iam, subname, circ_fname);

   for (kprime = 0; kprime < km; kprime++) {
      sprintf (varname, "vmix_matrix_%03d_CUR", kprime + 1);
      if (dbg_lvl)
         printf ("(%d) %s: reading %s from %s\n", iam, subname, varname, circ_fname);
      if (get_var_3d_double (circ_fname, varname, vmix_matrix_var))
         return 1;

      for (tracer_ind = 0; tracer_ind < coupled_tracer_cnt; tracer_ind++) {
         for (tracer_state_ind = 0; tracer_state_ind < tracer_state_len; tracer_state_ind++) {
            int i;
            int j;
            int k;
            int k2;
            int coef_ind;

            i = tracer_state_ind_to_int3[tracer_state_ind].i;
            j = tracer_state_ind_to_int3[tracer_state_ind].j;
            k = tracer_state_ind_to_int3[tracer_state_ind].k;

            coef_ind = coef_ind_vmix_non_nbr[tracer_ind][tracer_state_ind];

            for (k2 = 0; k2 < KMT[j][i]; k2++) {
               if (k2 == kprime)
                  nzval_row_wise[coef_ind] += vmix_matrix_var[k][j][i] * delta_t;
               coef_ind++;
            }
         }
      }
   }

   free_3d_double (vmix_matrix_var);

   if (dbg_lvl > 1) {
      printf ("(%d) exiting %s\n", iam, subname);
      fflush (stdout);
   }

   return 0;
}

/******************************************************************************/

int
add_vmix_file (void)
{
   char *subname = "add_vmix_file";
   char *varname;
   double ***VDC_TOTAL, ***VDC_READ;
   double fv;
   int i;
   int j;
   int k;

   int tracer_ind;
   int tracer_state_ind;

   if (dbg_lvl > 1) {
      printf ("(%d) entering %s\n", iam, subname);
      fflush (stdout);
   }

   if ((VDC_TOTAL = malloc_3d_double (km, jmt, imt)) == NULL) {
      fprintf (stderr, "(%d) malloc failed in %s for VDC_TOTAL\n", iam, subname);
      return 1;
   }
   if ((VDC_READ = malloc_3d_double (km, jmt, imt)) == NULL) {
      fprintf (stderr, "(%d) malloc failed in %s for VDC_READ\n", iam, subname);
      return 1;
   }

   varname = "VDC_S";
   if (dbg_lvl)
      printf ("(%d) %s: reading %s from %s for VDC\n", iam, subname, varname, circ_fname);
   if (get_var_3d_double (circ_fname, varname, VDC_TOTAL))
      return 1;
   if (get_att_double (circ_fname, varname, "_FillValue", &fv))
      return 1;
   set_fv_3d_double (fv, 0.0, VDC_TOTAL);

   varname = "VDC_GM";
   if (dbg_lvl)
      printf ("(%d) %s: reading %s from %s for VDC\n", iam, subname, varname, circ_fname);
   if (get_var_3d_double (circ_fname, varname, VDC_READ))
      return 1;
   if (get_att_double (circ_fname, varname, "_FillValue", &fv))
      return 1;
   set_fv_3d_double (fv, 0.0, VDC_READ);

   for (k = 0; k < km; k++)
      for (j = 1; j < jmt - 1; j++)
         for (i = 0; i < imt; i++)
            VDC_TOTAL[k][j][i] += VDC_READ[k][j][i];

   for (tracer_ind = 0; tracer_ind < coupled_tracer_cnt; tracer_ind++) {
      for (tracer_state_ind = 0; tracer_state_ind < tracer_state_len; tracer_state_ind++) {
         int i;
         int j;
         int k;
         int coef_ind;

         double ct;
         double cb;

         i = tracer_state_ind_to_int3[tracer_state_ind].i;
         j = tracer_state_ind_to_int3[tracer_state_ind].j;
         k = tracer_state_ind_to_int3[tracer_state_ind].k;

         /* cell 1 level shallower */
         if (k - 1 >= 0)
            ct = VDC_TOTAL[k - 1][j][i] / (0.5 * (dz[k - 1] + dz[k])) / dz[k] * delta_t;
         else
            ct = 0.0;

         /* cell 1 level deeper */
         if (k + 1 < KMT[j][i])
            cb = VDC_TOTAL[k][j][i] / (0.5 * (dz[k] + dz[k + 1])) / dz[k] * delta_t;
         else
            cb = 0.0;

         coef_ind = coef_ind_self[tracer_ind][tracer_state_ind];

         /* cell itself */
         nzval_row_wise[coef_ind] -= (ct + cb);
         coef_ind++;
         /* cell 1 level shallower */
         if (k - 1 >= 0) {
            nzval_row_wise[coef_ind] += ct;
            coef_ind++;
         }
         /* cell 1 level deeper */
         if (k + 1 < KMT[j][i]) {
            nzval_row_wise[coef_ind] += cb;
            coef_ind++;
         }
      }
   }

   free_3d_double (VDC_READ);
   free_3d_double (VDC_TOTAL);

   if (dbg_lvl > 1) {
      printf ("(%d) exiting %s\n", iam, subname);
      fflush (stdout);
   }

   return 0;
}

/******************************************************************************/

void
add_vmix_const (void)
{
   char *subname = "add_vmix_const";
   int tracer_ind;
   int tracer_state_ind;
   double vdc;

   if (dbg_lvl > 1) {
      printf ("(%d) entering %s\n", iam, subname);
      fflush (stdout);
   }

   vdc = 0.1;

   for (tracer_ind = 0; tracer_ind < coupled_tracer_cnt; tracer_ind++) {
      for (tracer_state_ind = 0; tracer_state_ind < tracer_state_len; tracer_state_ind++) {
         int i;
         int j;
         int k;
         int coef_ind;

         double ct;
         double cb;

         i = tracer_state_ind_to_int3[tracer_state_ind].i;
         j = tracer_state_ind_to_int3[tracer_state_ind].j;
         k = tracer_state_ind_to_int3[tracer_state_ind].k;

         /* cell 1 level shallower */
         if (k - 1 >= 0)
            ct = vdc / (0.5 * (dz[k - 1] + dz[k])) / dz[k] * delta_t;
         else
            ct = 0.0;

         /* cell 1 level deeper */
         if (k + 1 < KMT[j][i])
            cb = vdc / (0.5 * (dz[k] + dz[k + 1])) / dz[k] * delta_t;
         else
            cb = 0.0;

         coef_ind = coef_ind_self[tracer_ind][tracer_state_ind];

         /* cell itself */
         nzval_row_wise[coef_ind] -= (ct + cb);
         coef_ind++;
         /* cell 1 level shallower */
         if (k - 1 >= 0) {
            nzval_row_wise[coef_ind] += ct;
            coef_ind++;
         }
         /* cell 1 level deeper */
         if (k + 1 < KMT[j][i]) {
            nzval_row_wise[coef_ind] += cb;
            coef_ind++;
         }
      }
   }

   if (dbg_lvl > 1) {
      printf ("(%d) exiting %s\n", iam, subname);
      fflush (stdout);
   }
}

/******************************************************************************/


int
add_vmix (void)
{
   char *subname = "add_vmix";

   if (dbg_lvl > 1) {
      printf ("(%d) entering %s\n", iam, subname);
      fflush (stdout);
   }

   switch (vmix_opt) {
   case vmix_matrix_file:
      if (add_vmix_matrix_file ())
         return 1;
      break;
   case vmix_file:
      if (add_vmix_file ())
         return 1;
      break;
   case vmix_const:
      add_vmix_const ();
      break;
   case vmix_none:
      break;
   }

   if (dbg_lvl > 1) {
      printf ("(%d) vmix terms added\n\n", iam);
      fflush (stdout);
   }

   if (dbg_lvl > 1) {
      printf ("(%d) exiting %s\n", iam, subname);
      fflush (stdout);
   }

   return 0;
}

/******************************************************************************/

int
add_sink_pure_diag (void)
{
   char *subname = "add_sink_pure_diag";
   int tracer_ind;
   int tracer_state_ind;
   double ***SINK_RATE_FIELD;

   if (dbg_lvl > 1) {
      printf ("(%d) entering %s\n", iam, subname);
      fflush (stdout);
   }

   for (tracer_ind = 0; tracer_ind < coupled_tracer_cnt; tracer_ind++) {
      switch (per_tracer_opt[tracer_ind].sink_opt) {
      case sink_none:
         break;
      case sink_const:
         for (tracer_state_ind = 0; tracer_state_ind < tracer_state_len; tracer_state_ind++) {
            int coef_ind = coef_ind_self[tracer_ind][tracer_state_ind];
            nzval_row_wise[coef_ind] += -year_cnt * per_tracer_opt[tracer_ind].sink_rate;
         }
         if (dbg_lvl > 1) {
            printf ("(%d) sink const (%e) added for tracer %d\n\n", iam, per_tracer_opt[tracer_ind].sink_rate, tracer_ind);
            fflush (stdout);
         }
         break;
      case sink_const_shallow:
         for (tracer_state_ind = 0; tracer_state_ind < tracer_state_len; tracer_state_ind++) {
            int k = tracer_state_ind_to_int3[tracer_state_ind].k;
            int coef_ind = coef_ind_self[tracer_ind][tracer_state_ind];
            if (z_t[k] < per_tracer_opt[tracer_ind].sink_depth)
               nzval_row_wise[coef_ind] += -year_cnt * per_tracer_opt[tracer_ind].sink_rate;
         }
         if (dbg_lvl > 1) {
            printf ("(%d) sink const shallow (%e,%e) added for tracer %d\n\n", iam, per_tracer_opt[tracer_ind].sink_depth,
                    per_tracer_opt[tracer_ind].sink_rate, tracer_ind);
            fflush (stdout);
         }
         break;
      case sink_file:
         if ((SINK_RATE_FIELD = malloc_3d_double (km, jmt, imt)) == NULL) {
            fprintf (stderr, "(%d) malloc failed in %s for SINK_RATE_FIELD\n", iam, subname);
            return 1;
         }
         if (dbg_lvl)
            printf ("(%d) %s: reading %s from %s\n", iam, subname, per_tracer_opt[tracer_ind].sink_field_name, tracer_fname);
         if (get_var_3d_double (tracer_fname, per_tracer_opt[tracer_ind].sink_field_name, SINK_RATE_FIELD))
            return 1;
         for (tracer_state_ind = 0; tracer_state_ind < tracer_state_len; tracer_state_ind++) {
            int i = tracer_state_ind_to_int3[tracer_state_ind].i;
            int j = tracer_state_ind_to_int3[tracer_state_ind].j;
            int k = tracer_state_ind_to_int3[tracer_state_ind].k;
            int coef_ind = coef_ind_self[tracer_ind][tracer_state_ind];
            nzval_row_wise[coef_ind] += -year_cnt * SINK_RATE_FIELD[k][j][i];
         }
         free_3d_double (SINK_RATE_FIELD);
         if (dbg_lvl > 1) {
            printf ("(%d) file sink (%s,%s) added for tracer %d\n\n", iam, tracer_fname,
                    per_tracer_opt[tracer_ind].sink_field_name, tracer_ind);
            fflush (stdout);
         }
         break;
      }
   }

   if (dbg_lvl > 1) {
      printf ("(%d) exiting %s\n", iam, subname);
      fflush (stdout);
   }

   return 0;
}

/******************************************************************************/

int
add_sink_generic_tracer (void)
{
   char *subname = "add_sink_generic_tracer";
   int k2;
   int tracer_ind;
   int sink_var_exists;
   char *field_name;
   double ***SINK_RATE_FIELD_SAME_LEVEL;
   double ****SINK_RATE_FIELDS_SHALLOWER;
   int tracer_state_ind;

   if (dbg_lvl > 1) {
      printf ("(%d) entering %s\n", iam, subname);
      fflush (stdout);
   }

   if ((SINK_RATE_FIELD_SAME_LEVEL = malloc_3d_double (km, jmt, imt)) == NULL) {
      fprintf (stderr, "(%d) malloc failed in %s for SINK_RATE_FIELD_SAME_LEVEL\n", iam, subname);
      return 1;
   }

   if ((SINK_RATE_FIELDS_SHALLOWER = malloc ((size_t) km * sizeof (double ***))) == NULL) {
      fprintf (stderr, "(%d) malloc failed in %s for SINK_RATE_FIELD_SAME_LEVEL\n", iam, subname);
      return 1;
   }
   for (k2 = 0; k2 < km; k2++)
      SINK_RATE_FIELDS_SHALLOWER[k2] = NULL;

   for (tracer_ind = 0; tracer_ind < coupled_tracer_cnt; tracer_ind++) {
      if (per_tracer_opt[tracer_ind].sink_opt == sink_generic_tracer) {
         int kmax =
            (per_tracer_opt[tracer_ind].sink_generic_tracer_depends_layer_cnt ==
             -1) ? km - 1 : per_tracer_opt[tracer_ind].sink_generic_tracer_depends_layer_cnt - 1;

         /* read and add pure diagonal term, if present in input file */

         if ((field_name = malloc (13 + 2 * strlen (per_tracer_opt[tracer_ind].sink_generic_tracer_name))) == NULL) {
            fprintf (stderr, "(%d) malloc failed in %s for field_name for tracer %s\n", iam, subname,
                     per_tracer_opt[tracer_ind].sink_generic_tracer_name);
            return 1;
         }
         sprintf (field_name, "d_J_%s_d_%s", per_tracer_opt[tracer_ind].sink_generic_tracer_name,
                  per_tracer_opt[tracer_ind].sink_generic_tracer_name);
         if (var_exists_in_file (tracer_fname, field_name, &sink_var_exists)) {
            fprintf (stderr, "(%d) var_exists_in_file failed in %s for field_name %s for tracer %s\n", iam, subname,
                     field_name, per_tracer_opt[tracer_ind].sink_generic_tracer_name);
            return 1;
         }
         if (sink_var_exists) {
            if (dbg_lvl)
               printf ("(%d) %s: reading %s from %s\n", iam, subname, field_name, tracer_fname);
            if (get_var_3d_double (tracer_fname, field_name, SINK_RATE_FIELD_SAME_LEVEL))
               return 1;
            for (tracer_state_ind = 0; tracer_state_ind < tracer_state_len; tracer_state_ind++) {
               int i = tracer_state_ind_to_int3[tracer_state_ind].i;
               int j = tracer_state_ind_to_int3[tracer_state_ind].j;
               int k = tracer_state_ind_to_int3[tracer_state_ind].k;
               int coef_ind = coef_ind_self[tracer_ind][tracer_state_ind];
               nzval_row_wise[coef_ind] += delta_t * SINK_RATE_FIELD_SAME_LEVEL[k][j][i];
            }
         } else {
            if (dbg_lvl)
               printf ("(%d) %s: %s does not exist in %s\n", iam, subname, field_name, tracer_fname);
         }

         /* process levels shallower than each specific level */

         /* allocate space for and read in corresponding file variables, if present in input file */

         for (k2 = 0; k2 <= kmax; k2++) {
            sprintf (field_name, "d_J_%s_d_%s_k_%02d", per_tracer_opt[tracer_ind].sink_generic_tracer_name,
                     per_tracer_opt[tracer_ind].sink_generic_tracer_name, k2 + 1);
            if (var_exists_in_file (tracer_fname, field_name, &sink_var_exists)) {
               fprintf (stderr, "(%d) var_exists_in_file failed in %s for field_name %s for tracer %s\n", iam, subname,
                        field_name, per_tracer_opt[tracer_ind].sink_generic_tracer_name);
               return 1;
            }
            if (sink_var_exists) {
               if ((SINK_RATE_FIELDS_SHALLOWER[k2] = malloc_3d_double (km, jmt, imt)) == NULL) {
                  fprintf (stderr, "(%d) malloc failed in %s for SINK_RATE_FIELDS_SHALLOWER[%d]\n", iam, subname, k2);
                  return 1;
               }
               if (dbg_lvl)
                  printf ("(%d) %s: reading %s from %s\n", iam, subname, field_name, tracer_fname);
               if (get_var_3d_double (tracer_fname, field_name, SINK_RATE_FIELDS_SHALLOWER[k2]))
                  return 1;
            } else {
               if (dbg_lvl)
                  printf ("(%d) %s: %s does not exist in %s\n", iam, subname, field_name, tracer_fname);
               SINK_RATE_FIELDS_SHALLOWER[k2] = NULL;
            }
         }

         for (tracer_state_ind = 0; tracer_state_ind < tracer_state_len; tracer_state_ind++) {
            int i = tracer_state_ind_to_int3[tracer_state_ind].i;
            int j = tracer_state_ind_to_int3[tracer_state_ind].j;
            int k = tracer_state_ind_to_int3[tracer_state_ind].k;
            int coef_ind = coef_ind_sink_non_nbr[tracer_ind][tracer_state_ind];
            for (k2 = (k <= kmax) ? k : kmax; k2 >= 0; k2--) {
               if (SINK_RATE_FIELDS_SHALLOWER[k2] != NULL) {
                  nzval_row_wise[coef_ind] += delta_t * SINK_RATE_FIELDS_SHALLOWER[k2][k][j][i];
               }
               coef_ind++;
            }
         }

         /* deallocate allocated space */

         for (k2 = 0; k2 < kmax; k2++) {
            if (SINK_RATE_FIELDS_SHALLOWER[k2] != NULL) {
               free_3d_double (SINK_RATE_FIELDS_SHALLOWER[k2]);
               SINK_RATE_FIELDS_SHALLOWER[k2] = NULL;
            }
         }

         free (field_name);

         if (dbg_lvl > 1) {
            printf ("(%d) generic tracer sink added for tracer %d, %s\n\n", iam, tracer_ind,
                    per_tracer_opt[tracer_ind].sink_generic_tracer_name);
            fflush (stdout);
         }
      }
   }

   free (SINK_RATE_FIELDS_SHALLOWER);
   free_3d_double (SINK_RATE_FIELD_SAME_LEVEL);

   if (dbg_lvl > 1) {
      printf ("(%d) exiting %s\n", iam, subname);
      fflush (stdout);
   }

   return 0;
}

/******************************************************************************/

int
add_sink_coupled_tracers (void)
{
   char *subname = "add_sink_coupled_tracers";
   int tracer_ind;
   int tracer_ind_2;
   int sink_var_exists;
   char *field_name;
   double ****SINK_RATE_FIELDS;
   int tracer_state_ind;

   char **tracer_names = NULL;

   if (dbg_lvl > 1) {
      printf ("(%d) entering %s\n", iam, subname);
      fflush (stdout);
   }

   switch (coupled_tracer_opt) {
   case coupled_tracer_none:
      if (dbg_lvl > 1) {
         printf ("(%d) exiting %s\n", iam, subname);
         fflush (stdout);
      }
      return 0;
   case coupled_tracer_OCMIP_BGC_PO4_DOP:
      tracer_names = OCMIP_BGC_PO4_DOP_names;
      break;
   case coupled_tracer_DIC_SHADOW_ALK_SHADOW:
      tracer_names = DIC_SHADOW_ALK_SHADOW_names;
      break;
   }

   if (tracer_names != NULL) {
      if ((SINK_RATE_FIELDS = malloc ((size_t) coupled_tracer_cnt * sizeof (double ***))) == NULL) {
         fprintf (stderr, "(%d) malloc failed in %s for SINK_RATE_FIELDS\n", iam, subname);
         return 1;
      }
      for (tracer_ind = 0; tracer_ind < coupled_tracer_cnt; tracer_ind++)
         SINK_RATE_FIELDS[tracer_ind] = NULL;

      for (tracer_ind = 0; tracer_ind < coupled_tracer_cnt; tracer_ind++) {

         /* allocate space for and read in corresponding file variables, if present in input file */

         for (tracer_ind_2 = 0; tracer_ind_2 < coupled_tracer_cnt; tracer_ind_2++) {
            if (tracer_ind_2 == tracer_ind)
               continue;

            if ((field_name = malloc (8 + strlen (tracer_names[tracer_ind]) + strlen (tracer_names[tracer_ind_2]))) == NULL) {
               fprintf (stderr, "(%d) malloc failed in %s for field_name\n", iam, subname);
               return 1;
            }
            sprintf (field_name, "d_J_%s_d_%s", tracer_names[tracer_ind], tracer_names[tracer_ind_2]);
            if (var_exists_in_file (tracer_fname, field_name, &sink_var_exists)) {
               fprintf (stderr, "(%d) var_exists_in_file failed in %s for field_name %s for tracer %s\n", iam, subname,
                        field_name, per_tracer_opt[tracer_ind].sink_generic_tracer_name);
               return 1;
            }
            if (sink_var_exists) {
               if ((SINK_RATE_FIELDS[tracer_ind_2] = malloc_3d_double (km, jmt, imt)) == NULL) {
                  fprintf (stderr, "(%d) malloc failed in %s for SINK_RATE_FIELDS[%d]\n", iam, subname, tracer_ind_2);
                  return 1;
               }
               if (dbg_lvl)
                  printf ("(%d) %s: reading %s from %s\n", iam, subname, field_name, tracer_fname);
               if (get_var_3d_double (tracer_fname, field_name, SINK_RATE_FIELDS[tracer_ind_2]))
                  return 1;
            } else {
               if (dbg_lvl)
                  printf ("(%d) %s: %s does not exist in %s\n", iam, subname, field_name, tracer_fname);
               SINK_RATE_FIELDS[tracer_ind_2] = NULL;
            }

            free (field_name);
         }

         for (tracer_state_ind = 0; tracer_state_ind < tracer_state_len; tracer_state_ind++) {
            int i = tracer_state_ind_to_int3[tracer_state_ind].i;
            int j = tracer_state_ind_to_int3[tracer_state_ind].j;
            int k = tracer_state_ind_to_int3[tracer_state_ind].k;
            int coef_ind = coef_ind_sink_other_tracers[tracer_ind][tracer_state_ind];
            for (tracer_ind_2 = 0; tracer_ind_2 < coupled_tracer_cnt; tracer_ind_2++) {
               if (tracer_ind_2 == tracer_ind)
                  continue;
               if (SINK_RATE_FIELDS[tracer_ind_2] != NULL) {
                  nzval_row_wise[coef_ind] += delta_t * SINK_RATE_FIELDS[tracer_ind_2][k][j][i];
               }
               coef_ind++;
            }
         }

         /* deallocate allocated space */

         for (tracer_ind_2 = 0; tracer_ind_2 < coupled_tracer_cnt; tracer_ind_2++) {
            if (SINK_RATE_FIELDS[tracer_ind_2] != NULL) {
               free_3d_double (SINK_RATE_FIELDS[tracer_ind_2]);
               SINK_RATE_FIELDS[tracer_ind_2] = NULL;
            }
         }
      }
      free (SINK_RATE_FIELDS);
   }

   if (dbg_lvl > 1) {
      printf ("(%d) exiting %s\n", iam, subname);
      fflush (stdout);
   }

   return 0;
}

/******************************************************************************/

int
add_pv (void)
{
   char *subname = "add_pv";
   double **pv = NULL;
   int tracer_ind;
   int tracer_state_ind;

   if (dbg_lvl > 1) {
      printf ("(%d) entering %s\n", iam, subname);
      fflush (stdout);
   }

   for (tracer_ind = 0; tracer_ind < coupled_tracer_cnt; tracer_ind++) {
      if (per_tracer_opt[tracer_ind].pv_field_name != NULL) {
         if (tracer_fname == NULL) {
            fprintf (stderr, "(%d) %s:tracer_fname not specified for tracer pv %s\n", iam, subname,
                     per_tracer_opt[tracer_ind].pv_field_name);
            return 1;
         }
         if (pv == NULL) {
            if ((pv = malloc_2d_double (jmt, imt)) == NULL) {
               fprintf (stderr, "(%d) malloc failed in %s for pv\n", iam, subname);
               return 1;
            }
         }
         if (dbg_lvl)
            printf ("(%d) %s: reading %s for piston velocity from %s\n", iam, subname,
                    per_tracer_opt[tracer_ind].pv_field_name, tracer_fname);
         if (get_var_2d_double (tracer_fname, per_tracer_opt[tracer_ind].pv_field_name, pv))
            return 1;
         for (tracer_state_ind = 0; tracer_state_ind < tracer_state_len; tracer_state_ind++) {
            int i = tracer_state_ind_to_int3[tracer_state_ind].i;
            int j = tracer_state_ind_to_int3[tracer_state_ind].j;
            int k = tracer_state_ind_to_int3[tracer_state_ind].k;
            int coef_ind = coef_ind_self[tracer_ind][tracer_state_ind];
            if (k == 0)
               nzval_row_wise[coef_ind] -= pv[j][i] / dz[0] * delta_t;
         }
      }
   }

   if (dbg_lvl > 1) {
      printf ("(%d) pv terms added\n\n", iam);
      fflush (stdout);
   }

   if (pv != NULL)
      free_2d_double (pv);

   if (dbg_lvl > 1) {
      printf ("(%d) exiting %s\n", iam, subname);
      fflush (stdout);
   }

   return 0;
}

/******************************************************************************/

int
add_d_SF_d_TRACER (void)
{
   char *subname = "add_d_SF_d_TRACER";
   double **d_SF_d_TRACER = NULL;
   int tracer_ind;
   int tracer_state_ind;

   if (dbg_lvl > 1) {
      printf ("(%d) entering %s\n", iam, subname);
      fflush (stdout);
   }

   for (tracer_ind = 0; tracer_ind < coupled_tracer_cnt; tracer_ind++) {
      if (per_tracer_opt[tracer_ind].d_SF_d_TRACER_field_name != NULL) {
         if (tracer_fname == NULL) {
            fprintf (stderr, "(%d) %s:tracer_fname not specified for tracer d_SF_d_TRACER %s\n", iam, subname,
                     per_tracer_opt[tracer_ind].d_SF_d_TRACER_field_name);
            return 1;
         }
         if (d_SF_d_TRACER == NULL) {
            if ((d_SF_d_TRACER = malloc_2d_double (jmt, imt)) == NULL) {
               fprintf (stderr, "(%d) malloc failed in %s for d_SF_d_TRACER\n", iam, subname);
               return 1;
            }
         }
         if (dbg_lvl)
            printf ("(%d) %s: reading %s from %s\n", iam, subname, per_tracer_opt[tracer_ind].d_SF_d_TRACER_field_name,
                    tracer_fname);
         if (get_var_2d_double (tracer_fname, per_tracer_opt[tracer_ind].d_SF_d_TRACER_field_name, d_SF_d_TRACER))
            return 1;
         for (tracer_state_ind = 0; tracer_state_ind < tracer_state_len; tracer_state_ind++) {
            int i = tracer_state_ind_to_int3[tracer_state_ind].i;
            int j = tracer_state_ind_to_int3[tracer_state_ind].j;
            int k = tracer_state_ind_to_int3[tracer_state_ind].k;
            int coef_ind = coef_ind_self[tracer_ind][tracer_state_ind];
            if (k == 0)
               nzval_row_wise[coef_ind] += d_SF_d_TRACER[j][i] / dz[0] * delta_t;
         }
      }
   }

   if (dbg_lvl > 1) {
      printf ("(%d) d_SF_d_TRACER terms added\n\n", iam);
      fflush (stdout);
   }

   if (d_SF_d_TRACER != NULL)
      free_2d_double (d_SF_d_TRACER);

   if (dbg_lvl > 1) {
      printf ("(%d) exiting %s\n", iam, subname);
      fflush (stdout);
   }

   return 0;
}

/******************************************************************************/

int
add_sf_coupled_tracers (void)
{
   char *subname = "add_sf_coupled_tracers";
   int tracer_ind;
   int tracer_ind_2;
   int d_SF_d_TRACER_var_exists;
   char *field_name;
   double ***d_SF_d_TRACER_FIELDS;
   int tracer_state_ind;

   char **tracer_names = NULL;

   if (dbg_lvl > 1) {
      printf ("(%d) entering %s\n", iam, subname);
      fflush (stdout);
   }

   switch (coupled_tracer_opt) {
   case coupled_tracer_none:
      if (dbg_lvl > 1) {
         printf ("(%d) exiting %s\n", iam, subname);
         fflush (stdout);
      }
      return 0;
   case coupled_tracer_DIC_SHADOW_ALK_SHADOW:
      tracer_names = DIC_SHADOW_ALK_SHADOW_names;
      break;
   }

   if (tracer_names != NULL) {
      if ((d_SF_d_TRACER_FIELDS = malloc ((size_t) coupled_tracer_cnt * sizeof (double **))) == NULL) {
         fprintf (stderr, "(%d) malloc failed in %s for d_SF_d_TRACER_FIELDS\n", iam, subname);
         return 1;
      }
      for (tracer_ind = 0; tracer_ind < coupled_tracer_cnt; tracer_ind++)
         d_SF_d_TRACER_FIELDS[tracer_ind] = NULL;

      for (tracer_ind = 0; tracer_ind < coupled_tracer_cnt; tracer_ind++) {

         /* allocate space for and read in corresponding file variables, if present in input file */

         for (tracer_ind_2 = 0; tracer_ind_2 < coupled_tracer_cnt; tracer_ind_2++) {
            if (tracer_ind_2 == tracer_ind)
               continue;

            if ((field_name = malloc (9 + strlen (tracer_names[tracer_ind]) + strlen (tracer_names[tracer_ind_2]))) == NULL) {
               fprintf (stderr, "(%d) malloc failed in %s for field_name\n", iam, subname);
               return 1;
            }
            sprintf (field_name, "d_SF_%s_d_%s", tracer_names[tracer_ind], tracer_names[tracer_ind_2]);
            if (var_exists_in_file (tracer_fname, field_name, &d_SF_d_TRACER_var_exists)) {
               fprintf (stderr, "(%d) var_exists_in_file failed in %s for field_name %s for tracer %s\n", iam, subname,
                        field_name, per_tracer_opt[tracer_ind].sink_generic_tracer_name);
               return 1;
            }
            if (d_SF_d_TRACER_var_exists) {
               if ((d_SF_d_TRACER_FIELDS[tracer_ind_2] = malloc_2d_double (jmt, imt)) == NULL) {
                  fprintf (stderr, "(%d) malloc failed in %s for d_SF_d_TRACER_FIELDS[%d]\n", iam, subname, tracer_ind_2);
                  return 1;
               }
               if (dbg_lvl)
                  printf ("(%d) %s: reading %s from %s\n", iam, subname, field_name, tracer_fname);
               if (get_var_2d_double (tracer_fname, field_name, d_SF_d_TRACER_FIELDS[tracer_ind_2]))
                  return 1;
            } else {
               if (dbg_lvl)
                  printf ("(%d) %s: %s does not exist in %s\n", iam, subname, field_name, tracer_fname);
               d_SF_d_TRACER_FIELDS[tracer_ind_2] = NULL;
            }

            free (field_name);
         }

         for (tracer_state_ind = 0; tracer_state_ind < tracer_state_len; tracer_state_ind++) {
            int i = tracer_state_ind_to_int3[tracer_state_ind].i;
            int j = tracer_state_ind_to_int3[tracer_state_ind].j;
            int k = tracer_state_ind_to_int3[tracer_state_ind].k;
            int coef_ind = coef_ind_sink_other_tracers[tracer_ind][tracer_state_ind];
            for (tracer_ind_2 = 0; tracer_ind_2 < coupled_tracer_cnt; tracer_ind_2++) {
               if (tracer_ind_2 == tracer_ind)
                  continue;
               if ((d_SF_d_TRACER_FIELDS[tracer_ind_2] != NULL) && (k == 0)) {
                  nzval_row_wise[coef_ind] += delta_t * d_SF_d_TRACER_FIELDS[tracer_ind_2][j][i] / dz[0];
               }
               coef_ind++;
            }
         }

         /* deallocate allocated space */

         for (tracer_ind_2 = 0; tracer_ind_2 < coupled_tracer_cnt; tracer_ind_2++) {
            if (d_SF_d_TRACER_FIELDS[tracer_ind_2] != NULL) {
               free_2d_double (d_SF_d_TRACER_FIELDS[tracer_ind_2]);
               d_SF_d_TRACER_FIELDS[tracer_ind_2] = NULL;
            }
         }
      }
      free (d_SF_d_TRACER_FIELDS);
   }

   if (dbg_lvl > 1) {
      printf ("(%d) exiting %s\n", iam, subname);
      fflush (stdout);
   }

   return 0;
}

/******************************************************************************/

/* sum entries with multiple values */

void
sum_dup_vals (void)
{
   char *subname = "sum_dup_vals";
   int flat_ind;
   int coef_ind;
   int coef_ind_dup;
   int dup_cnt;

   if (dbg_lvl > 1) {
      printf ("(%d) entering %s\n", iam, subname);
      fflush (stdout);
   }

   dup_cnt = 0;
   for (flat_ind = 0; flat_ind < flat_len; flat_ind++)
      for (coef_ind = rowptr[flat_ind]; coef_ind < rowptr[flat_ind + 1]; coef_ind++)
         for (coef_ind_dup = coef_ind + 1; coef_ind_dup < rowptr[flat_ind + 1]; coef_ind_dup++)
            if (colind[coef_ind_dup] == colind[coef_ind]) {
               nzval_row_wise[coef_ind] += nzval_row_wise[coef_ind_dup];
               nzval_row_wise[coef_ind_dup] = 0.0;
               dup_cnt++;
            }
   if (dbg_lvl)
      printf ("(%d) subname = %s, dup_cnt = %d\n", iam, subname, dup_cnt);

   if (dbg_lvl > 1) {
      printf ("(%d) exiting %s\n", iam, subname);
      fflush (stdout);
   }
}

/******************************************************************************/

/* remove zeros from matrix */

void
strip_matrix_zeros (void)
{
   char *subname = "strip_matrix_zeros";
   int flat_ind;
   int coef_ind_all;
   int coef_ind_nonzero;

   if (dbg_lvl > 1) {
      printf ("(%d) entering %s\n", iam, subname);
      fflush (stdout);
   }

   coef_ind_all = 0;
   coef_ind_nonzero = 0;
   for (flat_ind = 0; flat_ind < flat_len; flat_ind++) {
      for (; coef_ind_all < rowptr[flat_ind + 1]; coef_ind_all++)
         if (nzval_row_wise[coef_ind_all] != 0.0) {
            nzval_row_wise[coef_ind_nonzero] = nzval_row_wise[coef_ind_all];
            colind[coef_ind_nonzero] = colind[coef_ind_all];
            coef_ind_nonzero++;
         }
      rowptr[flat_ind + 1] = coef_ind_nonzero;
   }
   if (dbg_lvl)
      printf ("(%d) subname = %s, nnz_pre = %d, nnz_new = %d\n", iam, subname, nnz, coef_ind_nonzero);
   nnz = coef_ind_nonzero;

   if (dbg_lvl > 1) {
      printf ("(%d) exiting %s\n", iam, subname);
      fflush (stdout);
   }
}

/******************************************************************************/

void
check_matrix_diag (void)
{
   char *subname = "check_matrix_diag";
   int flat_ind;
   int coef_ind;

   if (dbg_lvl > 1) {
      printf ("(%d) entering %s\n", iam, subname);
      fflush (stdout);
   }

   for (flat_ind = 0; flat_ind < flat_len; flat_ind++) {
      int diag_found = 0;
      for (coef_ind = rowptr[flat_ind]; coef_ind < rowptr[flat_ind + 1]; coef_ind++)
         if (colind[coef_ind] == flat_ind) {
            diag_found = 1;
            if (nzval_row_wise[coef_ind] == 0.0)
               printf
                  ("(%d) subname = %s, zero on diagonal, flat_ind = %d, flat_ind = %d, colind = %lld\n", iam,
                   subname, flat_ind, flat_ind, (long long) colind[coef_ind]);
         }
      if (!diag_found) {
         printf ("(%d) subname = %s, no diagonal found, flat_ind = %d, flat_ind = %d, colind = ", iam, subname, flat_ind,
                 flat_ind);
         for (coef_ind = rowptr[flat_ind]; coef_ind < rowptr[flat_ind + 1]; coef_ind++)
            printf ("(%d)  %lld", iam, (long long) colind[coef_ind]);
         putchar ('\n');
      }
   }

   if (dbg_lvl > 1) {
      printf ("(%d) exiting %s\n", iam, subname);
      fflush (stdout);
   }
}

/******************************************************************************/

void
sort_cols_one_row (int_t len, double *nzval_row_wise_loc, int_t * colind_loc)
{
   int i, j;

   /* use insertion sort */

   for (i = 1; i < len; i++) {
      int_t key = colind_loc[i];
      double nzval = nzval_row_wise_loc[i];
      for (j = i - 1; j >= 0 && colind_loc[j] > key; --j) {
         colind_loc[j + 1] = colind_loc[j];
         nzval_row_wise_loc[j + 1] = nzval_row_wise_loc[j];
         colind_loc[j] = key;
         nzval_row_wise_loc[j] = nzval;
      }
   }
}

/******************************************************************************/

void
sort_cols_all_rows (void)
{
   char *subname = "sort_cols_all_rows";
   int flat_ind;

   if (dbg_lvl > 1) {
      printf ("(%d) entering %s\n", iam, subname);
      fflush (stdout);
   }

   for (flat_ind = 0; flat_ind < flat_len; flat_ind++)
      sort_cols_one_row (rowptr[flat_ind + 1] - rowptr[flat_ind], nzval_row_wise + rowptr[flat_ind], colind + rowptr[flat_ind]);

   if (dbg_lvl > 1) {
      printf ("(%d) exiting %s\n", iam, subname);
      fflush (stdout);
   }
}

/******************************************************************************/

int
gen_sparse_matrix (double day_cnt)
{
   char *subname = "gen_sparse_matrix";
   delta_t = 60.0 * 60.0 * 24.0 * day_cnt;
   year_cnt = day_cnt / 365.0;

   if (dbg_lvl > 1) {
      printf ("(%d) entering %s\n", iam, subname);
      fflush (stdout);
   }

   comp_flat_len ();
   comp_nnz ();

   if (allocate_matrix_arrays ())
      return 1;

   if (init_matrix ())
      return 1;

   /* add_adv must be called first, so that adv_enforce_divfree can work properly */
   if (add_adv ())
      return 1;

   if (l_adv_enforce_divfree)
      adv_enforce_divfree ();

   if (add_hmix ())
      return 1;

   if (add_vmix ())
      return 1;

   if (add_sink_pure_diag ())
      return 1;

   if (add_sink_generic_tracer ())
      return 1;

   if (add_sink_coupled_tracers ())
      return 1;

   if (add_pv ())
      return 1;

   if (add_d_SF_d_TRACER ())
      return 1;

   if (add_sf_coupled_tracers ())
      return 1;

   sum_dup_vals ();

   strip_matrix_zeros ();

   check_matrix_diag ();

   sort_cols_all_rows ();

   if (dbg_lvl > 1) {
      printf ("(%d) exiting %s\n", iam, subname);
      fflush (stdout);
   }

   return 0;
}

/******************************************************************************/

int
put_sparse_matrix (char *fname)
{
   char *subname = "put_sparse_matrix";
   int status;
   int ncid;
   int dimids[1];
   int flat_len_p1_dimid;
   int nnz_dimid;
   int varid;

   if (dbg_lvl > 1) {
      printf ("(%d) entering %s\n", iam, subname);
      fflush (stdout);
   }

   if ((status = nc_open (fname, NC_WRITE, &ncid)))
      return handle_nc_error (subname, "nc_open", fname, status);

   /* define new dimensions and get dimids for existing ones */

   if ((status = nc_redef (ncid)))
      return handle_nc_error (subname, "nc_redef", fname, status);

   if ((status = nc_def_dim (ncid, "nnz", nnz, &nnz_dimid)))
      return handle_nc_error (subname, "nc_def_dim", "nnz", status);

   if ((status = nc_def_dim (ncid, "flat_len_p1", flat_len + 1, &flat_len_p1_dimid)))
      return handle_nc_error (subname, "nc_def_dim", "flat_len_p1", status);

   /* define new variables */

   if ((status = nc_def_var (ncid, "coupled_tracer_cnt", NC_INT, 0, dimids, &varid)))
      return handle_nc_error (subname, "nc_def_var", "coupled_tracer_cnt", status);

   dimids[0] = nnz_dimid;
   if ((status = nc_def_var (ncid, "nzval_row_wise", NC_DOUBLE, 1, dimids, &varid)))
      return handle_nc_error (subname, "nc_def_var", "nzval_row_wise", status);

   dimids[0] = nnz_dimid;
   if ((status = nc_def_var (ncid, "colind", NC_INT, 1, dimids, &varid)))
      return handle_nc_error (subname, "nc_def_var", "colind", status);

   dimids[0] = flat_len_p1_dimid;
   if ((status = nc_def_var (ncid, "rowptr", NC_INT, 1, dimids, &varid)))
      return handle_nc_error (subname, "nc_def_var", "rowptr", status);

   if ((status = nc_close (ncid)))
      return handle_nc_error (subname, "nc_close", fname, status);

   /* write out new variables */
   /* copy int_t variables to an int array before writing */

   if (put_var_1d_int (fname, "coupled_tracer_cnt", &coupled_tracer_cnt))
      return 1;

   if (put_var_1d_double (fname, "nzval_row_wise", nzval_row_wise))
      return 1;

   /* colind is an array of type int_t, write it to a temporary arrays of type int, for safe writing to NetCDF */
   {
      int *tmp_int_array;
      int ind;
      if ((tmp_int_array = malloc ((size_t) nnz * sizeof (int))) == NULL) {
         fprintf (stderr, "(%d) malloc failed in %s for temp nnz int array\n", iam, subname);
         return 1;
      }
      for (ind = 0; ind < nnz; ind++)
         tmp_int_array[ind] = colind[ind];
      if (put_var_1d_int (fname, "colind", tmp_int_array))
         return 1;
      free (tmp_int_array);
   }

   /* rowptr is an array of type int_t, write it to a temporary arrays of type int, for safe writing to NetCDF */
   {
      int *tmp_int_array;
      int ind;
      if ((tmp_int_array = malloc ((size_t) (flat_len + 1) * sizeof (int))) == NULL) {
         fprintf (stderr, "(%d) malloc failed in %s for (flat_len+1) nnz int array\n", iam, subname);
         return 1;
      }
      for (ind = 0; ind < flat_len + 1; ind++)
         tmp_int_array[ind] = rowptr[ind];
      if (put_var_1d_int (fname, "rowptr", tmp_int_array))
         return 1;
      free (tmp_int_array);
   }

   if (dbg_lvl > 1) {
      printf ("(%d) exiting %s\n", iam, subname);
      fflush (stdout);
   }

   return 0;
}

/******************************************************************************/

int
get_sparse_matrix (char *fname)
{
   char *subname = "get_sparse_matrix";
   int status;
   int ncid;
   int dimid;
   size_t dimlen;

   if (dbg_lvl > 1) {
      printf ("(%d) entering %s\n", iam, subname);
      fflush (stdout);
   }

   if ((status = nc_open (fname, NC_NOWRITE, &ncid)))
      return handle_nc_error (subname, "nc_open", fname, status);

   if ((status = nc_inq_dimid (ncid, "nnz", &dimid)))
      return handle_nc_error (subname, "nc_inq_dimid", "nnz", status);

   if ((status = nc_inq_dimlen (ncid, dimid, &dimlen)))
      return handle_nc_error (subname, "nc_inq_dimlen", "nnz", status);
   nnz = dimlen;

   if ((status = nc_inq_dimid (ncid, "flat_len_p1", &dimid)))
      return handle_nc_error (subname, "nc_inq_dimid", "flat_len_p1", status);

   if ((status = nc_inq_dimlen (ncid, dimid, &dimlen)))
      return handle_nc_error (subname, "nc_inq_dimlen", "flat_len_p1", status);
   flat_len = dimlen - 1;

   if ((status = nc_close (ncid)))
      return handle_nc_error (subname, "nc_close", fname, status);

   if (get_var_1d_int (fname, "coupled_tracer_cnt", &coupled_tracer_cnt))
      return 1;

   if (dbg_lvl && (iam == 0)) {
      printf ("(%d) %s: coupled_tracer_cnt = %d\n", iam, subname, coupled_tracer_cnt);
      printf ("(%d) %s: nnz = %d\n", iam, subname, nnz);
      printf ("(%d) %s: flat_len = %d\n", iam, subname, flat_len);
   }

   /* allocate arrays for sparse matrix */

   if (allocate_matrix_arrays ())
      return 1;

   /* read variables */
   /* read int_t variables to an int array and copy to int_t arrays */

   if (get_var_1d_double (fname, "nzval_row_wise", nzval_row_wise))
      return 1;

   {
      int *tmp_int_array;
      int ind;
      if ((tmp_int_array = malloc ((size_t) nnz * sizeof (int))) == NULL) {
         fprintf (stderr, "(%d) malloc failed in %s for nnz length int array\n", iam, subname);
         return 1;
      }
      if (get_var_1d_int (fname, "colind", tmp_int_array))
         return 1;
      for (ind = 0; ind < nnz; ind++)
         colind[ind] = tmp_int_array[ind];
      free (tmp_int_array);
   }

   {
      int *tmp_int_array;
      int ind;
      if ((tmp_int_array = malloc ((size_t) (flat_len + 1) * sizeof (int))) == NULL) {
         fprintf (stderr, "(%d) malloc failed in %s for (flat_len+1) length int array\n", iam, subname);
         return 1;
      }
      if (get_var_1d_int (fname, "rowptr", tmp_int_array))
         return 1;
      for (ind = 0; ind < flat_len + 1; ind++)
         rowptr[ind] = tmp_int_array[ind];
      free (tmp_int_array);
   }

   if (dbg_lvl > 1) {
      printf ("(%d) exiting %s\n", iam, subname);
      fflush (stdout);
   }

   return 0;
}

/******************************************************************************/

void
free_sparse_matrix (void)
{
   char *subname = "free_sparse_matrix";
   int tracer_ind;

   if (dbg_lvl > 1) {
      printf ("(%d) entering %s\n", iam, subname);
      fflush (stdout);
   }

   free (rowptr);
   free (colind);
   free (nzval_row_wise);

   for (tracer_ind = 0; tracer_ind < coupled_tracer_cnt; tracer_ind++) {
      free (coef_ind_self[tracer_ind]);
      free (coef_ind_adv_non_nbr[tracer_ind]);
      free (coef_ind_hmix_non_nbr[tracer_ind]);
      free (coef_ind_vmix_non_nbr[tracer_ind]);
      free (coef_ind_sink_non_nbr[tracer_ind]);
      free (coef_ind_sink_other_tracers[tracer_ind]);
   }

   free (coef_ind_self);
   free (coef_ind_adv_non_nbr);
   free (coef_ind_hmix_non_nbr);
   free (coef_ind_vmix_non_nbr);
   free (coef_ind_sink_non_nbr);
   free (coef_ind_sink_other_tracers);

   if (dbg_lvl > 1) {
      printf ("(%d) exiting %s\n", iam, subname);
      fflush (stdout);
   }
}
