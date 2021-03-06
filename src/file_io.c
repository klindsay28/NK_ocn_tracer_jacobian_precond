#include <stdio.h>

#include "netcdf.h"

#include "file_io.h"
#include "globals.h"

/******************************************************************************/

int
handle_nc_error (char *subname, char *cdf_subname, char *msg, int status)
{
   fprintf (stderr, "(%d) ERROR returned from netCDF routine\n"
            "\tsubname     : %s\n\tcdf_subname : %s\n\tmsg         : %s\n\tnetCDF msg  : %s\n", iam,
            subname, cdf_subname, msg, nc_strerror (status));
   return status;
}

/******************************************************************************/

int
var_exists_in_file (char *fname, char *varname, int *retval)
{
   char *subname = "var_exists_in_file";
   int status;
   int ncid;
   int varid;

   if ((status = nc_open (fname, NC_NOWRITE, &ncid)))
      return handle_nc_error (subname, "nc_open", fname, status);

   status = nc_inq_varid (ncid, varname, &varid);

   if ((status != NC_NOERR) && (status != NC_ENOTVAR))
      return handle_nc_error (subname, "nc_inq_varid", varname, status);

   *retval = (status == NC_NOERR);

   if ((status = nc_close (ncid)))
      return handle_nc_error (subname, "nc_close", fname, status);

   return 0;
}

/******************************************************************************/

int
get_att_double (char *fname, char *varname, char *attname, double *val)
{
   char *subname = "get_att_double";
   int status;
   int ncid;
   int varid;

   if ((status = nc_open (fname, NC_NOWRITE, &ncid)))
      return handle_nc_error (subname, "nc_open", fname, status);

   if ((status = nc_inq_varid (ncid, varname, &varid)))
      return handle_nc_error (subname, "nc_inq_varid", varname, status);

   if ((status = nc_get_att_double (ncid, varid, attname, val)))
      return handle_nc_error (subname, "nc_get_att_double", varname, status);

   if ((status = nc_close (ncid)))
      return handle_nc_error (subname, "nc_close", fname, status);

   return 0;
}

/******************************************************************************/

int
get_var_1d_int (char *fname, char *varname, int *field)
{
   char *subname = "get_var_1d_int";
   int status;
   int ncid;
   int varid;

   if ((status = nc_open (fname, NC_NOWRITE, &ncid)))
      return handle_nc_error (subname, "nc_open", fname, status);

   if ((status = nc_inq_varid (ncid, varname, &varid)))
      return handle_nc_error (subname, "nc_inq_varid", varname, status);

   if ((status = nc_get_var_int (ncid, varid, field)))
      return handle_nc_error (subname, "nc_get_var_int", varname, status);

   if ((status = nc_close (ncid)))
      return handle_nc_error (subname, "nc_close", fname, status);

   return 0;
}

/******************************************************************************/

int
get_var_2d_int (char *fname, char *varname, int **field)
{
   char *subname = "get_var_2d_int";
   int status;
   int ncid;
   int varid;

   if ((status = nc_open (fname, NC_NOWRITE, &ncid)))
      return handle_nc_error (subname, "nc_open", fname, status);

   if ((status = nc_inq_varid (ncid, varname, &varid)))
      return handle_nc_error (subname, "nc_inq_varid", varname, status);

   if ((status = nc_get_var_int (ncid, varid, field[0])))
      return handle_nc_error (subname, "nc_get_var_int", varname, status);

   if ((status = nc_close (ncid)))
      return handle_nc_error (subname, "nc_close", fname, status);

   return 0;
}

/******************************************************************************/

int
get_var_3d_int (char *fname, char *varname, int ***field)
{
   char *subname = "get_var_3d_int";
   int status;
   int ncid;
   int varid;

   if ((status = nc_open (fname, NC_NOWRITE, &ncid)))
      return handle_nc_error (subname, "nc_open", fname, status);

   if ((status = nc_inq_varid (ncid, varname, &varid)))
      return handle_nc_error (subname, "nc_inq_varid", varname, status);

   if ((status = nc_get_var_int (ncid, varid, field[0][0])))
      return handle_nc_error (subname, "nc_get_var_int", varname, status);

   if ((status = nc_close (ncid)))
      return handle_nc_error (subname, "nc_close", fname, status);

   return 0;
}

/******************************************************************************/

int
put_var_1d_int (char *fname, char *varname, int *field)
{
   char *subname = "put_var_1d_int";
   int status;
   int ncid;
   int varid;

   if ((status = nc_open (fname, NC_WRITE, &ncid)))
      return handle_nc_error (subname, "nc_open", fname, status);

   if ((status = nc_inq_varid (ncid, varname, &varid)))
      return handle_nc_error (subname, "nc_inq_varid", varname, status);

   if ((status = nc_put_var_int (ncid, varid, field)))
      return handle_nc_error (subname, "nc_put_var_int", varname, status);

   if ((status = nc_close (ncid)))
      return handle_nc_error (subname, "nc_close", fname, status);

   return 0;
}

/******************************************************************************/

int
put_var_2d_int (char *fname, char *varname, int **field)
{
   char *subname = "put_var_2d_int";
   int status;
   int ncid;
   int varid;

   if ((status = nc_open (fname, NC_WRITE, &ncid)))
      return handle_nc_error (subname, "nc_open", fname, status);

   if ((status = nc_inq_varid (ncid, varname, &varid)))
      return handle_nc_error (subname, "nc_inq_varid", varname, status);

   if ((status = nc_put_var_int (ncid, varid, field[0])))
      return handle_nc_error (subname, "nc_put_var_int", varname, status);

   if ((status = nc_close (ncid)))
      return handle_nc_error (subname, "nc_close", fname, status);

   return 0;
}

/******************************************************************************/

int
put_var_3d_int (char *fname, char *varname, int ***field)
{
   char *subname = "put_var_3d_int";
   int status;
   int ncid;
   int varid;

   if ((status = nc_open (fname, NC_WRITE, &ncid)))
      return handle_nc_error (subname, "nc_open", fname, status);

   if ((status = nc_inq_varid (ncid, varname, &varid)))
      return handle_nc_error (subname, "nc_inq_varid", varname, status);

   if ((status = nc_put_var_int (ncid, varid, field[0][0])))
      return handle_nc_error (subname, "nc_put_var_int", varname, status);

   if ((status = nc_close (ncid)))
      return handle_nc_error (subname, "nc_close", fname, status);

   return 0;
}

/******************************************************************************/

int
get_var_1d_double (char *fname, char *varname, double *field)
{
   char *subname = "get_var_1d_double";
   int status;
   int ncid;
   int varid;

   if ((status = nc_open (fname, NC_NOWRITE, &ncid)))
      return handle_nc_error (subname, "nc_open", fname, status);

   if ((status = nc_inq_varid (ncid, varname, &varid)))
      return handle_nc_error (subname, "nc_inq_varid", varname, status);

   if ((status = nc_get_var_double (ncid, varid, field)))
      return handle_nc_error (subname, "nc_get_var_double", varname, status);

   if ((status = nc_close (ncid)))
      return handle_nc_error (subname, "nc_close", fname, status);

   return 0;
}

/******************************************************************************/

int
get_var_2d_double (char *fname, char *varname, double **field)
{
   char *subname = "get_var_2d_double";
   int status;
   int ncid;
   int varid;

   if ((status = nc_open (fname, NC_NOWRITE, &ncid)))
      return handle_nc_error (subname, "nc_open", fname, status);

   if ((status = nc_inq_varid (ncid, varname, &varid)))
      return handle_nc_error (subname, "nc_inq_varid", varname, status);

   if ((status = nc_get_var_double (ncid, varid, field[0])))
      return handle_nc_error (subname, "nc_get_var_double", varname, status);

   if ((status = nc_close (ncid)))
      return handle_nc_error (subname, "nc_close", fname, status);

   return 0;
}

/******************************************************************************/

int
get_var_3d_double (char *fname, char *varname, double ***field)
{
   char *subname = "get_var_3d_double";
   int status;
   int ncid;
   int varid;

   if ((status = nc_open (fname, NC_NOWRITE, &ncid)))
      return handle_nc_error (subname, "nc_open", fname, status);

   if ((status = nc_inq_varid (ncid, varname, &varid)))
      return handle_nc_error (subname, "nc_inq_varid", varname, status);

   if ((status = nc_get_var_double (ncid, varid, field[0][0])))
      return handle_nc_error (subname, "nc_get_var_double", varname, status);

   if ((status = nc_close (ncid)))
      return handle_nc_error (subname, "nc_close", fname, status);

   return 0;
}

/******************************************************************************/

int
put_var_1d_double (char *fname, char *varname, double *field)
{
   char *subname = "put_var_1d_double";
   int status;
   int ncid;
   int varid;

   if ((status = nc_open (fname, NC_WRITE, &ncid)))
      return handle_nc_error (subname, "nc_open", fname, status);

   if ((status = nc_inq_varid (ncid, varname, &varid)))
      return handle_nc_error (subname, "nc_inq_varid", varname, status);

   if ((status = nc_put_var_double (ncid, varid, field)))
      return handle_nc_error (subname, "nc_put_var_double", varname, status);

   if ((status = nc_close (ncid)))
      return handle_nc_error (subname, "nc_close", fname, status);

   return 0;
}

/******************************************************************************/

int
put_var_2d_double (char *fname, char *varname, double **field)
{
   char *subname = "put_var_2d_double";
   int status;
   int ncid;
   int varid;

   if ((status = nc_open (fname, NC_WRITE, &ncid)))
      return handle_nc_error (subname, "nc_open", fname, status);

   if ((status = nc_inq_varid (ncid, varname, &varid)))
      return handle_nc_error (subname, "nc_inq_varid", varname, status);

   if ((status = nc_put_var_double (ncid, varid, field[0])))
      return handle_nc_error (subname, "nc_put_var_double", varname, status);

   if ((status = nc_close (ncid)))
      return handle_nc_error (subname, "nc_close", fname, status);

   return 0;
}

/******************************************************************************/

int
put_var_3d_double (char *fname, char *varname, double ***field)
{
   char *subname = "put_var_3d_double";
   int status;
   int ncid;
   int varid;

   if ((status = nc_open (fname, NC_WRITE, &ncid)))
      return handle_nc_error (subname, "nc_open", fname, status);

   if ((status = nc_inq_varid (ncid, varname, &varid)))
      return handle_nc_error (subname, "nc_inq_varid", varname, status);

   if ((status = nc_put_var_double (ncid, varid, field[0][0])))
      return handle_nc_error (subname, "nc_put_var_double", varname, status);

   if ((status = nc_close (ncid)))
      return handle_nc_error (subname, "nc_close", fname, status);

   return 0;
}
