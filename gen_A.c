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

char *opt_fname = NULL;
double day_cnt;
char *matrix_fname = NULL;

/******************************************************************************/

int
parse_cmd_line (int argc, char **argv)
{
   char *usage_msg = "usage: gen_matrix_file [-h] [-D dbg_lvl] [-o opt_fname] matrix_fname";
   extern char *optarg;
   extern int optind;
   char *optstring = "D:o:h";
   int opt;

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
      case 'o':
         opt_fname = optarg;
         break;
      default:
         fprintf (stderr, "(%d) internal error: unhandled option '-%c'\n", iam, opt);
         return 1;
      }
   }
   if (optind != argc - 1) {
      fprintf (stderr, "(%d) unexpected number of arguments\n%s\n", iam, usage_msg);
      return 1;
   }
   matrix_fname = argv[optind++];
   return 0;
}

/******************************************************************************/

int
grow_per_tracer_opt (int prev_tracer_cnt, int new_tracer_cnt)
{
   char *subname = "grow_per_tracer_opt";
   int tracer_ind;

   if ((per_tracer_opt = realloc (per_tracer_opt, (size_t) new_tracer_cnt * sizeof (per_tracer_opt_t))) == NULL) {
      fprintf (stderr, "(%d) realloc failed in %s for grow_per_tracer_opt\n", iam, subname);
      return 1;
   }

   /* initialize new array elements */
   for (tracer_ind = prev_tracer_cnt; tracer_ind < new_tracer_cnt; tracer_ind++) {
      per_tracer_opt[tracer_ind].sink_opt = sink_none;
      per_tracer_opt[tracer_ind].sink_rate = 1.21e-4;   /* radiocarbon decay rate */
      per_tracer_opt[tracer_ind].sink_depth = 10.0e2;   /* 10m */
      per_tracer_opt[tracer_ind].sink_field_name = NULL;
      per_tracer_opt[tracer_ind].sink_generic_tracer_name = NULL;
      per_tracer_opt[tracer_ind].sink_generic_tracer_depends_layer_cnt = -1;

      per_tracer_opt[tracer_ind].pv_field_name = NULL;
      per_tracer_opt[tracer_ind].d_SF_d_TRACER_field_name = NULL;
   }

   return 0;
}

/******************************************************************************/

int
set_opt_defaults (void)
{
   day_cnt = 365.0;
   adv_opt = adv_cent;
   l_adv_enforce_divfree = 1;
   hmix_opt = hmix_isop_file;
   vmix_opt = vmix_file;
   coupled_tracer_cnt = 1;
   if (grow_per_tracer_opt (0, 1)) {
      fprintf (stderr, "(%d) error from grow_per_tracer_opt\n", iam);
      return 1;
   }
   coupled_tracer_opt = coupled_tracer_none;
   return 0;
}

/******************************************************************************/

int
read_opt_file ()
{
   char *subname = "read_opt_file";
   FILE *fp;
   char *optname;
   char *optval;
   int new_coupled_tracer_cnt;
   int tracer_ind = 0;

#define MAX_LINE_LEN 256
   char line[MAX_LINE_LEN];

   if (opt_fname == NULL)
      return 0;

   if ((fp = fopen (opt_fname, "r")) == NULL) {
      fprintf (stderr, "(%d) fopen failed in %s for %s\n", iam, subname, opt_fname);
      return 1;
   }

   while (fgets (line, MAX_LINE_LEN, fp) != NULL) {
      optname = strtok (line, " \n");
      if ((optval = strtok (NULL, " \n")) == NULL) {
         fprintf (stderr, "(%d) unspecified value for %s\n", iam, optname);
         return 1;
      }
      if (strcmp (optname, "day_cnt") == 0) {
         if (parse_to_double (optval, &day_cnt)) {
            fprintf (stderr, "(%d) error parsing argument '%s' for option '%s'\n", iam, optval, optname);
            return 1;
         }
      } else if (strcmp (optname, "reg_fname") == 0) {
         if ((reg_fname = malloc (1 + strlen (optval))) == NULL) {
            fprintf (stderr, "(%d) malloc failed in %s for reg_fname\n", iam, subname);
            return 1;
         }
         strcpy (reg_fname, optval);
      } else if (strcmp (optname, "circ_fname") == 0) {
         if ((circ_fname = malloc (1 + strlen (optval))) == NULL) {
            fprintf (stderr, "(%d) malloc failed in %s for circ_fname\n", iam, subname);
            return 1;
         }
         strcpy (circ_fname, optval);
      } else if (strcmp (optname, "adv_type") == 0) {
         if (strcmp (optval, "none") == 0)
            adv_opt = adv_none;
         else if (strcmp (optval, "donor") == 0)
            adv_opt = adv_donor;
         else if (strncmp (optval, "centered", 4) == 0)
            adv_opt = adv_cent;
         else if (strcmp (optval, "upwind3") == 0)
            adv_opt = adv_upwind3;
         else {
            fprintf (stderr, "(%d) unknown %s: %s\n", iam, optname, optval);
            return 1;
         }
      } else if (strcmp (optname, "l_adv_enforce_divfree") == 0) {
         if (strcmp (optval, "0") == 0)
            l_adv_enforce_divfree = 0;
         else if (strcmp (optval, "1") == 0)
            l_adv_enforce_divfree = 1;
         else {
            fprintf (stderr, "(%d) unknown %s: %s\n", iam, optname, optval);
            return 1;
         }
      } else if (strcmp (optname, "hmix_type") == 0) {
         if (strcmp (optval, "none") == 0)
            hmix_opt = hmix_none;
         else if (strcmp (optval, "const") == 0)
            hmix_opt = hmix_const;
         else if (strcmp (optval, "hor_file") == 0)
            hmix_opt = hmix_hor_file;
         else if (strcmp (optval, "isop_file") == 0)
            hmix_opt = hmix_isop_file;
         else {
            fprintf (stderr, "(%d) unknown %s: %s\n", iam, optname, optval);
            return 1;
         }
      } else if (strcmp (optname, "vmix_type") == 0) {
         if (strcmp (optval, "none") == 0)
            vmix_opt = vmix_none;
         else if (strcmp (optval, "const") == 0)
            vmix_opt = vmix_const;
         else if (strcmp (optval, "file") == 0)
            vmix_opt = vmix_file;
         else if (strcmp (optval, "matrix_file") == 0)
            vmix_opt = vmix_matrix_file;
         else {
            fprintf (stderr, "(%d) unknown %s: %s\n", iam, optname, optval);
            return 1;
         }
      } else if (strcmp (optname, "tracer_fname") == 0) {
         if ((tracer_fname = malloc (1 + strlen (optval))) == NULL) {
            fprintf (stderr, "(%d) malloc failed in %s for tracer_fname\n", iam, subname);
            return 1;
         }
         strcpy (tracer_fname, optval);
      } else if (strcmp (optname, "coupled_tracer_cnt") == 0) {
         if (parse_to_int (optval, &new_coupled_tracer_cnt)) {
            fprintf (stderr, "(%d) error parsing argument '%s' for option '%s'\n", iam, optval, optname);
            return 1;
         }
         if (grow_per_tracer_opt (coupled_tracer_cnt, new_coupled_tracer_cnt)) {
            fprintf (stderr, "(%d) error from grow_per_tracer_opt\n", iam);
            return 1;
         }
         coupled_tracer_cnt = new_coupled_tracer_cnt;
         if ((coupled_tracer_cnt < 1) || (coupled_tracer_cnt > 2)) {
            fprintf (stderr, "(%d) coupled_tracer_cnt = %d not supported\n", iam, coupled_tracer_cnt);
            return 1;
         }
      } else if (strcmp (optname, "tracer_ind") == 0) {
         if (parse_to_int (optval, &tracer_ind)) {
            fprintf (stderr, "(%d) error parsing argument '%s' for option '%s'\n", iam, optval, optname);
            return 1;
         }
         if ((tracer_ind < 0) || (tracer_ind >= coupled_tracer_cnt)) {
            fprintf (stderr, "(%d) tracer_ind = %d out of bounds for coupled_tracer_cnt = %d\n", iam, tracer_ind,
                     coupled_tracer_cnt);
            return 1;
         }
      } else if (strcmp (optname, "sink_type") == 0) {
         if (strcmp (optval, "none") == 0)
            per_tracer_opt[tracer_ind].sink_opt = sink_none;
         else if (strcmp (optval, "const") == 0)
            per_tracer_opt[tracer_ind].sink_opt = sink_const;
         else if (strcmp (optval, "const_shallow") == 0)
            per_tracer_opt[tracer_ind].sink_opt = sink_const_shallow;
         else if (strcmp (optval, "file") == 0)
            per_tracer_opt[tracer_ind].sink_opt = sink_file;
         else if (strcmp (optval, "generic_tracer") == 0)
            per_tracer_opt[tracer_ind].sink_opt = sink_generic_tracer;
         else {
            fprintf (stderr, "(%d) unknown %s: %s\n", iam, optname, optval);
            return 1;
         }
         if ((per_tracer_opt[tracer_ind].sink_opt == sink_const)
             || (per_tracer_opt[tracer_ind].sink_opt == sink_const_shallow)) {
            if ((optval = strtok (NULL, " \n")) == NULL) {
               fprintf (stderr, "(%d) unspecified sink_rate\n", iam);
               return 1;
            }
            if (parse_to_double (optval, &per_tracer_opt[tracer_ind].sink_rate)) {
               fprintf (stderr, "(%d) error parsing argument '%s' for option '%s'\n", iam, optval, optname);
               return 1;
            }
            if (per_tracer_opt[tracer_ind].sink_opt == sink_const_shallow) {
               if ((optval = strtok (NULL, " \n")) == NULL) {
                  fprintf (stderr, "(%d) unspecified sink_depth\n", iam);
                  return 1;
               }
               if (parse_to_double (optval, &per_tracer_opt[tracer_ind].sink_depth)) {
                  fprintf (stderr, "(%d) error parsing argument '%s' for option '%s'\n", iam, optval, optname);
                  return 1;
               }
            }
         }
         if (per_tracer_opt[tracer_ind].sink_opt == sink_file) {
            if ((optval = strtok (NULL, " \n")) == NULL) {
               fprintf (stderr, "(%d) unspecified sink_field_name\n", iam);
               return 1;
            }
            if ((per_tracer_opt[tracer_ind].sink_field_name = malloc (1 + strlen (optval))) == NULL) {
               fprintf (stderr, "(%d) malloc failed in %s for sink_field_name\n", iam, subname);
               return 1;
            }
            strcpy (per_tracer_opt[tracer_ind].sink_field_name, optval);
         }
         if (per_tracer_opt[tracer_ind].sink_opt == sink_generic_tracer) {
            if ((optval = strtok (NULL, " \n")) == NULL) {
               fprintf (stderr, "(%d) unspecified sink_generic_tracer_name\n", iam);
               return 1;
            }
            if ((per_tracer_opt[tracer_ind].sink_generic_tracer_name = malloc (1 + strlen (optval))) == NULL) {
               fprintf (stderr, "(%d) malloc failed in %s for sink_generic_tracer_name\n", iam, subname);
               return 1;
            }
            strcpy (per_tracer_opt[tracer_ind].sink_generic_tracer_name, optval);
            if ((optval = strtok (NULL, " \n")) != NULL) {
               if (parse_to_int (optval, &per_tracer_opt[tracer_ind].sink_generic_tracer_depends_layer_cnt)) {
                  fprintf (stderr, "(%d) error parsing sink_generic_tracer_depends_layer_cnt\n", iam);
                  return 1;
               }
            }
         }
      } else if (strcmp (optname, "pv") == 0) {
         if ((per_tracer_opt[tracer_ind].pv_field_name = malloc (1 + strlen (optval))) == NULL) {
            fprintf (stderr, "(%d) malloc failed in %s for pv_field_name\n", iam, subname);
            return 1;
         }
         strcpy (per_tracer_opt[tracer_ind].pv_field_name, optval);
      } else if (strcmp (optname, "sf") == 0) {
         if ((per_tracer_opt[tracer_ind].d_SF_d_TRACER_field_name = malloc (1 + strlen (optval))) == NULL) {
            fprintf (stderr, "(%d) malloc failed in %s for d_SF_d_TRACER_field_name\n", iam, subname);
            return 1;
         }
         strcpy (per_tracer_opt[tracer_ind].d_SF_d_TRACER_field_name, optval);
      } else if (strcmp (optname, "coupled_tracer_type") == 0) {
         if (strcmp (optval, "none") == 0)
            coupled_tracer_opt = coupled_tracer_none;
         else if (strcmp (optval, "OCMIP_BGC_PO4_DOP") == 0)
            coupled_tracer_opt = coupled_tracer_OCMIP_BGC_PO4_DOP;
         else if (strcmp (optval, "DIC_SHADOW_ALK_SHADOW") == 0)
            coupled_tracer_opt = coupled_tracer_DIC_SHADOW_ALK_SHADOW;
         else {
            fprintf (stderr, "(%d) unknown %s: %s\n", iam, optname, optval);
            return 1;
         }
      } else {
         fprintf (stderr, "(%d) unknown option name: %s\n", iam, optname);
         return 1;
      }
   }

   fclose (fp);

   if (coupled_tracer_cnt == 2) {
      if ((coupled_tracer_opt != coupled_tracer_OCMIP_BGC_PO4_DOP)
          && (coupled_tracer_opt != coupled_tracer_DIC_SHADOW_ALK_SHADOW)) {
         fprintf (stderr,
                  "(%d) coupled_tracer_cnt = 2 only supported for "
                  "coupled_tracer_type = OCMIP_BGC_PO4_DOP, DIC_SHADOW_ALK_SHADOW\n", iam);
         return 1;
      }
   }

   return 0;
}

/******************************************************************************/

void
write_opts (void)
{
   int tracer_ind;

   if (dbg_lvl) {
      printf ("(%d) dbg_lvl                    = %d\n", iam, dbg_lvl);
      printf ("(%d) day_cnt                    = %e\n", iam, day_cnt);
      printf ("(%d) reg_fname                  = %s\n", iam, reg_fname ? reg_fname : "none");
      printf ("(%d) circ_fname                 = %s\n", iam, circ_fname);
      switch (adv_opt) {
      case adv_none:
         printf ("(%d) adv_opt                    = %s\n", iam, "none");
         break;
      case adv_donor:
         printf ("(%d) adv_opt                    = %s\n", iam, "donor");
         break;
      case adv_cent:
         printf ("(%d) adv_opt                    = %s\n", iam, "centered");
         break;
      case adv_upwind3:
         printf ("(%d) adv_opt                    = %s\n", iam, "upwind3");
         break;
      }
      printf ("(%d) l_adv_enforce_divfree      = %d\n", iam, l_adv_enforce_divfree);
      switch (hmix_opt) {
      case hmix_none:
         printf ("(%d) hmix_opt                   = %s\n", iam, "none");
         break;
      case hmix_const:
         printf ("(%d) hmix_opt                   = %s\n", iam, "const");
         break;
      case hmix_hor_file:
         printf ("(%d) hmix_opt                   = %s\n", iam, "hor_file");
         break;
      case hmix_isop_file:
         printf ("(%d) hmix_opt                   = %s\n", iam, "isop_file");
         break;
      }
      switch (vmix_opt) {
      case vmix_none:
         printf ("(%d) vmix_opt                   = %s\n", iam, "none");
         break;
      case vmix_const:
         printf ("(%d) vmix_opt                   = %s\n", iam, "const");
         break;
      case vmix_file:
         printf ("(%d) vmix_opt                   = %s\n", iam, "file");
         break;
      case vmix_matrix_file:
         printf ("(%d) vmix_opt                   = %s\n", iam, "matrix_file");
         break;
      }
      printf ("(%d) tracer_fname               = %s\n", iam, tracer_fname ? tracer_fname : "none");
      printf ("(%d) coupled_tracer_cnt         = %d\n", iam, coupled_tracer_cnt);
      for (tracer_ind = 0; tracer_ind < coupled_tracer_cnt; tracer_ind++) {
         printf ("(%d) options for tracer %d\n", iam, tracer_ind);
         switch (per_tracer_opt[tracer_ind].sink_opt) {
         case sink_none:
            printf ("(%d)    sink_opt                = %s\n", iam, "none");
            break;
         case sink_const:
            printf ("(%d)    sink_opt                = %s\n", iam, "const");
            printf ("(%d)    sink_rate               = %e\n", iam, per_tracer_opt[tracer_ind].sink_rate);
            break;
         case sink_const_shallow:
            printf ("(%d)    sink_opt                = %s\n", iam, "const_shallow");
            printf ("(%d)    sink_rate               = %e\n", iam, per_tracer_opt[tracer_ind].sink_rate);
            printf ("(%d)    sink_depth              = %e\n", iam, per_tracer_opt[tracer_ind].sink_depth);
            break;
         case sink_file:
            printf ("(%d)    sink_opt                = %s\n", iam, "file");
            printf ("(%d)    sink_field_name         = %s\n", iam, per_tracer_opt[tracer_ind].sink_field_name);
            break;
         case sink_generic_tracer:
            printf ("(%d)    sink_opt                = %s\n", iam, "generic_tracer");
            printf ("(%d)    sink_generic_tracer_name= %s\n", iam, per_tracer_opt[tracer_ind].sink_generic_tracer_name);
            printf ("(%d)    depends_layer_cnt       = %d\n", iam,
                    per_tracer_opt[tracer_ind].sink_generic_tracer_depends_layer_cnt);
            break;
         }
         printf ("(%d)    pv_field_name           = %s\n", iam,
                 per_tracer_opt[tracer_ind].pv_field_name ? per_tracer_opt[tracer_ind].pv_field_name : "none");
         printf ("(%d)    d_SF_d_TRACER_field_name= %s\n", iam,
                 per_tracer_opt[tracer_ind].d_SF_d_TRACER_field_name ?
                 per_tracer_opt[tracer_ind].d_SF_d_TRACER_field_name : "none");
      }
      switch (coupled_tracer_opt) {
      case coupled_tracer_none:
         printf ("(%d) coupled_tracer_opt         = %s\n", iam, "none");
         break;
      case coupled_tracer_OCMIP_BGC_PO4_DOP:
         printf ("(%d) coupled_tracer_opt         = %s\n", iam, "OCMIP_BGC_PO4_DOP");
         break;
      case coupled_tracer_DIC_SHADOW_ALK_SHADOW:
         printf ("(%d) coupled_tracer_opt         = %s\n", iam, "DIC_SHADOW_ALK_SHADOW");
         break;
      }
      printf ("(%d) matrix_fname               = %s\n\n", iam, matrix_fname);
   }
}

/******************************************************************************/

int
main (int argc, char *argv[])
{
   iam = 0;
   dbg_lvl = 0;

   if (parse_cmd_line (argc, argv))
      exit (EXIT_FAILURE);

   if (set_opt_defaults ())
      exit (EXIT_FAILURE);

   if (read_opt_file ())
      exit (EXIT_FAILURE);

   write_opts ();

   if (get_grid_info (circ_fname, reg_fname))
      exit (EXIT_FAILURE);

   if (put_grid_info (matrix_fname))
      exit (EXIT_FAILURE);

   if (gen_ind_maps ())
      exit (EXIT_FAILURE);

   if (put_ind_maps (matrix_fname))
      exit (EXIT_FAILURE);

   if (gen_sparse_matrix (day_cnt))
      exit (EXIT_FAILURE);

   if (put_sparse_matrix (matrix_fname))
      exit (EXIT_FAILURE);

   free_sparse_matrix ();
   free_ind_maps ();
   free_grid_info ();

   free (per_tracer_opt);

   exit (EXIT_SUCCESS);
}
