
/******************************************************************************/
/* function prototypes                                                        */
/******************************************************************************/

int get_grid_dims (char *fname);
int get_grid_info (char *circ_fname, char *reg_fname);
int put_grid_info (char *fname);
void free_grid_info (void);

/******************************************************************************/
/* external variable declarations                                             */
/******************************************************************************/

extern char *circ_fname;

extern char *reg_fname;

extern int imt;
extern int jmt;
extern int km;

extern double *z_t;
extern double *dz;
extern double **TLONG;
extern double **TLAT;

extern int **KMT;
extern int **KMU;

extern double **TAREA;
