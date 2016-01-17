
/******************************************************************************/
/* function prototypes                                                        */
/******************************************************************************/

int gen_ind_maps (void);
int put_ind_maps (char *fname);
int get_ind_maps (char *fname);
void free_ind_maps (void);
int gen_sparse_matrix (int day_cnt);
int put_sparse_matrix (char *fname);
int get_sparse_matrix (char *fname);
void free_sparse_matrix (void);

/******************************************************************************/
/* typedef's                                                                  */
/******************************************************************************/

typedef struct
{
   int i;
   int j;
   int k;
} int3;

typedef enum
{ adv_none, adv_donor, adv_cent, adv_upwind3 } adv_opt_t;

typedef enum
{ hmix_none, hmix_const, hmix_hor_file, hmix_isop_file } hmix_opt_t;

typedef enum
{ vmix_const, vmix_file, vmix_matrix_file } vmix_opt_t;

typedef enum
{ sink_none, sink_const, sink_const_shallow, sink_file, sink_tracer } sink_opt_t;

/******************************************************************************/
/* external variable declarations                                             */
/******************************************************************************/

extern int tracer_state_len;
extern int ***int3_to_tracer_state_ind;
extern int3 *tracer_state_ind_to_int3;

extern int coupled_tracer_cnt;

extern int flat_len;
extern int nnz;
extern double *nzval_row_wise;
extern int_t *colind;
extern int_t *rowptr;

extern adv_opt_t adv_opt;

extern hmix_opt_t hmix_opt;

extern vmix_opt_t vmix_opt;

extern sink_opt_t sink_opt;
extern double sink_rate;
extern double sink_depth;
extern char *sink_file_name;
extern char *sink_field_name;
extern char *sink_tracer_name;
extern char *pv_file_name;
extern char *pv_field_name;
extern char *d_SF_d_TRACER_file_name;
extern char *d_SF_d_TRACER_field_name;
