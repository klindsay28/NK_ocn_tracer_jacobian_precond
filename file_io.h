
/******************************************************************************/
/* function prototypes                                                        */
/******************************************************************************/

int handle_nc_error (char *subname, char *cdf_subname, char *msg, int status);

int get_var_1d_int (char *fname, char *varname, int *field);
int get_var_2d_int (char *fname, char *varname, int **field);
int get_var_3d_int (char *fname, char *varname, int ***field);

int put_var_1d_int (char *fname, char *varname, int *field);
int put_var_2d_int (char *fname, char *varname, int **field);
int put_var_3d_int (char *fname, char *varname, int ***field);

int get_var_1d_double (char *fname, char *varname, double *field);
int get_var_2d_double (char *fname, char *varname, double **field);
int get_var_3d_double (char *fname, char *varname, double ***field);

int put_var_1d_double (char *fname, char *varname, double *field);
int put_var_2d_double (char *fname, char *varname, double **field);
int put_var_3d_double (char *fname, char *varname, double ***field);
