
/******************************************************************************/
/* function prototypes                                                        */
/******************************************************************************/

int **malloc_2d_int (int jmt, int imt);
void free_2d_int (int **ptr);
int ***malloc_3d_int (int km, int jmt, int imt);
void free_3d_int (int ***ptr);

double **malloc_2d_double (int jmt, int imt);
void free_2d_double (double **ptr);
double ***malloc_3d_double (int km, int jmt, int imt);
void free_3d_double (double ***ptr);
