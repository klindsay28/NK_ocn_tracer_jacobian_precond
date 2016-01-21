#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "globals.h"

#include "file_io.h"
#include "memory.h"
#include "netcdf.h"

#include "grid.h"

char *circ_fname = NULL;

char *reg_fname = NULL;

int imt;
int jmt;
int km;

double *z_t;
double *dz;
double **TLONG;
double **TLAT;

int **KMT;
int **KMU;

double **TAREA;

/******************************************************************************/

int
get_grid_dims (char *fname)
{
   char *subname = "get_grid_dims";
   int status;
   int ncid;
   int dimid;
   size_t dimlen;

   if (dbg_lvl > 1) {
      printf ("entering %s\n", subname);
      fflush (stdout);
   }

   if ((status = nc_open (fname, NC_NOWRITE, &ncid)))
      return handle_nc_error (subname, "nc_open", fname, status);

   if ((status = nc_inq_dimid (ncid, "nlon", &dimid)))
      return handle_nc_error (subname, "nc_inq_dimid", "nlon", status);

   if ((status = nc_inq_dimlen (ncid, dimid, &dimlen)))
      return handle_nc_error (subname, "nc_inq_dimlen", "nlon", status);
   imt = dimlen;

   if ((status = nc_inq_dimid (ncid, "nlat", &dimid)))
      return handle_nc_error (subname, "nc_inq_dimid", "nlat", status);

   if ((status = nc_inq_dimlen (ncid, dimid, &dimlen)))
      return handle_nc_error (subname, "nc_inq_dimlen", "nlat", status);
   jmt = dimlen;

   if ((status = nc_inq_dimid (ncid, "z_t", &dimid)))
      return handle_nc_error (subname, "nc_inq_dimid", "z_t", status);

   if ((status = nc_inq_dimlen (ncid, dimid, &dimlen)))
      return handle_nc_error (subname, "nc_inq_dimlen", "z_t", status);
   km = dimlen;

   if ((status = nc_close (ncid)))
      return handle_nc_error (subname, "nc_close", fname, status);

   if (dbg_lvl && (iam == 0)) {
      printf ("imt = %d\n", imt);
      printf ("jmt = %d\n", jmt);
      printf ("km  = %d\n", km);
   }

   if (dbg_lvl > 1) {
      printf ("exiting %s\n", subname);
      fflush (stdout);
   }

   return 0;
}

/******************************************************************************/

int
get_grid_info (char *circ_fname, char *reg_fname)
{
   char *subname = "get_grid_info";
   int **DYN_REGMASK;
   int j;
   int i;
   int ip1;

   if (dbg_lvl > 1) {
      printf ("entering %s\n", subname);
      fflush (stdout);
   }

   if (get_grid_dims (circ_fname))
      return 1;

   if ((z_t = malloc ((size_t) km * sizeof (double))) == NULL) {
      fprintf (stderr, "malloc failed in %s for z_t\n", subname);
      return 1;
   }
   if (get_var_1d_double (circ_fname, "z_t", z_t))
      return 1;

   if ((dz = malloc ((size_t) km * sizeof (double))) == NULL) {
      fprintf (stderr, "malloc failed in %s for dz\n", subname);
      return 1;
   }
   if (get_var_1d_double (circ_fname, "dz", dz))
      return 1;

   if ((TLONG = malloc_2d_double (jmt, imt)) == NULL) {
      fprintf (stderr, "malloc failed in %s for TLONG\n", subname);
      return 1;
   }
   if (get_var_2d_double (circ_fname, "TLONG", TLONG))
      return 1;

   if ((TLAT = malloc_2d_double (jmt, imt)) == NULL) {
      fprintf (stderr, "malloc failed in %s for TLAT\n", subname);
      return 1;
   }
   if (get_var_2d_double (circ_fname, "TLAT", TLAT))
      return 1;

   if ((KMT = malloc_2d_int (jmt, imt)) == NULL) {
      fprintf (stderr, "malloc failed in %s for KMT\n", subname);
      return 1;
   }
   if (get_var_2d_int (circ_fname, "KMT", KMT))
      return 1;

   /* set KMT to 0 in ignored regions */
   if (reg_fname != NULL) {
      if ((DYN_REGMASK = malloc_2d_int (jmt, imt)) == NULL) {
         fprintf (stderr, "malloc failed in %s for DYN_REGMASK\n", subname);
         return 1;
      }
      if (get_var_2d_int (reg_fname, "DYN_REGMASK", DYN_REGMASK))
         return 1;
      for (j = 1; j < jmt - 1; j++)
         for (i = 0; i < imt; i++)
            if (DYN_REGMASK[j][i] < 0)
               KMT[j][i] = 0;
      free_2d_int (DYN_REGMASK);
   }

   /* verify that KMT is 0 on southern- and northern-most rows */
   {
      int south_flag;
      int north_flag;

      south_flag = north_flag = 0;
      for (i = 0; i < imt; i++) {
         if (KMT[0][i])
            south_flag = 1;
         if (KMT[jmt - 1][i])
            north_flag = 1;
      }
      if (south_flag)
         fprintf (stderr, "non-land found on southern-most row in %s\n", subname);
      if (north_flag)
         fprintf (stderr, "non-land found on northern-most row in %s\n", subname);
      if (south_flag || north_flag)
         return 1;
   }

   /* construct KMU from (potentially modified) KMT */
   if ((KMU = malloc_2d_int (jmt, imt)) == NULL) {
      fprintf (stderr, "malloc failed in %s for KMU\n", subname);
      return 1;
   }
   for (j = 0; j < jmt - 1; j++)
      for (i = 0; i < imt; i++) {
         ip1 = (i < imt - 1) ? i + 1 : 0;
         KMU[j][i] = KMT[j][i];
         KMU[j][i] = (KMT[j + 1][i] < KMU[j][i]) ? KMT[j + 1][i] : KMU[j][i];
         KMU[j][i] = (KMT[j][ip1] < KMU[j][i]) ? KMT[j][ip1] : KMU[j][i];
         KMU[j][i] = (KMT[j + 1][ip1] < KMU[j][i]) ? KMT[j + 1][ip1] : KMU[j][i];
      }
   j = jmt - 1;
   for (i = 0; i < imt; i++) {
      KMU[j][i] = 0;
   }

   if ((TAREA = malloc_2d_double (jmt, imt)) == NULL) {
      fprintf (stderr, "malloc failed in %s for TAREA\n", subname);
      return 1;
   }
   if (get_var_2d_double (circ_fname, "TAREA", TAREA))
      return 1;

   if (dbg_lvl > 1) {
      printf ("exiting %s\n", subname);
      fflush (stdout);
   }

   return 0;
}

/******************************************************************************/

int
put_grid_info (char *fname)
{
   char *subname = "put_grid_info";
   int status;
   int ncid;
   int dimids[2];
   int nlon_dimid;
   int nlat_dimid;
   int z_t_dimid;
   int varid;
   char *string;

   if (dbg_lvl > 1) {
      printf ("entering %s\n", subname);
      fflush (stdout);
   }

   if ((status = nc_create (fname, NC_64BIT_OFFSET, &ncid)))
      return handle_nc_error (subname, "nc_create", fname, status);

   /* define dimensions */

   if ((status = nc_def_dim (ncid, "nlon", imt, &nlon_dimid)))
      return handle_nc_error (subname, "nc_def_dimid", "nlon", status);

   if ((status = nc_def_dim (ncid, "nlat", jmt, &nlat_dimid)))
      return handle_nc_error (subname, "nc_def_dimid", "nlat", status);

   if ((status = nc_def_dim (ncid, "z_t", km, &z_t_dimid)))
      return handle_nc_error (subname, "nc_def_dimid", "z_t", status);

   /* define variables */

   dimids[0] = z_t_dimid;

   if ((status = nc_def_var (ncid, "z_t", NC_DOUBLE, 1, dimids, &varid)))
      return handle_nc_error (subname, "nc_def_var", "z_t", status);
   string = "depth from surface to midpoint of layer";
   if ((status = nc_put_att_text (ncid, varid, "long_name", strlen (string), string)))
      return handle_nc_error (subname, "nc_put_att_text", "z_t", status);
   string = "centimeters";
   if ((status = nc_put_att_text (ncid, varid, "units", strlen (string), string)))
      return handle_nc_error (subname, "nc_put_att_text", "z_t", status);
   string = "down";
   if ((status = nc_put_att_text (ncid, varid, "positive", strlen (string), string)))
      return handle_nc_error (subname, "nc_put_att_text", "z_t", status);

   dimids[0] = nlat_dimid;
   dimids[1] = nlon_dimid;

   if ((status = nc_def_var (ncid, "TLONG", NC_DOUBLE, 2, dimids, &varid)))
      return handle_nc_error (subname, "nc_def_var", "TLONG", status);
   string = "array of t-grid longitudes";
   if ((status = nc_put_att_text (ncid, varid, "long_name", strlen (string), string)))
      return handle_nc_error (subname, "nc_put_att_text", "TLONG", status);
   string = "degrees_east";
   if ((status = nc_put_att_text (ncid, varid, "units", strlen (string), string)))
      return handle_nc_error (subname, "nc_put_att_text", "TLONG", status);

   if ((status = nc_def_var (ncid, "TLAT", NC_DOUBLE, 2, dimids, &varid)))
      return handle_nc_error (subname, "nc_def_var", "TLAT", status);
   string = "array of t-grid latitudes";
   if ((status = nc_put_att_text (ncid, varid, "long_name", strlen (string), string)))
      return handle_nc_error (subname, "nc_put_att_text", "TLAT", status);
   string = "degrees_north";
   if ((status = nc_put_att_text (ncid, varid, "units", strlen (string), string)))
      return handle_nc_error (subname, "nc_put_att_text", "TLAT", status);

   if ((status = nc_def_var (ncid, "KMT", NC_INT, 2, dimids, &varid)))
      return handle_nc_error (subname, "nc_def_var", "KMT", status);
   string = "k Index of Deepest Grid Cell on T Grid";
   if ((status = nc_put_att_text (ncid, varid, "long_name", strlen (string), string)))
      return handle_nc_error (subname, "nc_put_att_text", "KMT", status);
   string = "TLONG TLAT";
   if ((status = nc_put_att_text (ncid, varid, "coordinates", strlen (string), string)))
      return handle_nc_error (subname, "nc_put_att_text", "KMT", status);

   if ((status = nc_close (ncid)))
      return handle_nc_error (subname, "nc_close", fname, status);

   /* write out variables */

   if (put_var_1d_double (fname, "z_t", z_t))
      return 1;

   if (put_var_2d_double (fname, "TLONG", TLONG))
      return 1;
   if (put_var_2d_double (fname, "TLAT", TLAT))
      return 1;
   if (put_var_2d_int (fname, "KMT", KMT))
      return 1;

   if (dbg_lvl > 1) {
      printf ("exiting %s\n", subname);
      fflush (stdout);
   }

   return 0;
}

/******************************************************************************/

void
free_grid_info (void)
{
   free (z_t);
   free (dz);
   free_2d_double (TLONG);
   free_2d_double (TLAT);
   free_2d_int (KMT);
   free_2d_int (KMU);
   free_2d_double (TAREA);
}
