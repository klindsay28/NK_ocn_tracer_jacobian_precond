
/******************************************************************************/
/* function prototypes                                                        */
/******************************************************************************/

int gen_ind_maps (void);
int put_ind_maps (char *fname);
int get_ind_maps (char *fname);
void free_ind_maps (void);
int gen_sparse_matrix (double day_cnt);
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

typedef struct
{
   sink_opt_t sink_opt;
   double sink_rate;            /* loss rate, units = 1/yr */
   double sink_depth;           /* depth threshold for sink_const_shallow, units = cm (same as model's z_t) */
   char *sink_file_name;
   char *sink_field_name;
   char *sink_tracer_name;

   char *pv_file_name;
   char *pv_field_name;
   char *d_SF_d_TRACER_file_name;
   char *d_SF_d_TRACER_field_name;
} per_tracer_opt_t;

typedef enum
{ coupled_tracer_none, coupled_tracer_OCMIP_BGC_PO4_DOP } coupled_tracer_opt_t;

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

extern per_tracer_opt_t *per_tracer_opt;

extern coupled_tracer_opt_t coupled_tracer_opt;
