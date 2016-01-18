#include <stdlib.h>
#include <unistd.h>
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
   char *usage_msg =
      "usage: gen_matrix_file [-h] [-D dbg_lvl] [-o opt_fname] circ_fname matrix_fname";
   extern char *optarg;
   extern int optind;
   char *optstring = "D:o:h";
   int opt;

   while ((opt = getopt (argc, argv, optstring)) != -1) {
      switch (opt) {
      case '?':
      case 'h':
         fprintf (stderr, "%s\n", usage_msg);
         return 1;
      case 'D':
         if (parse_to_int (optarg, &dbg_lvl)) {
            fprintf (stderr, "error parsing argument '%s' for option '%c'\n", optarg, opt);
            return 1;
         }
         break;
      case 'o':
         opt_fname = optarg;
         break;
      default:
         fprintf (stderr, "internal error: unhandled option '-%c'\n", opt);
         return 1;
      }
   }
   if (optind != argc - 2) {
      fprintf (stderr, "unexpected number of arguments\n%s\n", usage_msg);
      return 1;
   }
   circ_fname = argv[optind++];
   matrix_fname = argv[optind++];
   return 0;
}

/******************************************************************************/

void
set_opt_defaults (void)
{
   day_cnt = 365.0;
   adv_opt = adv_cent;
   hmix_opt = hmix_isop_file;
   vmix_opt = vmix_file;
   coupled_tracer_cnt = 1;
   sink_opt = sink_none;
   sink_rate = 1.0e-5;
   sink_depth = 10.0e2;
}

/******************************************************************************/

int
read_opt_file ()
{
   char *subname = "read_opt_file";
   FILE *fp;
   char *optname;
   char *optval;

#define MAX_LINE_LEN 256
   char line[MAX_LINE_LEN];

   if (opt_fname == NULL)
      return 0;

   if ((fp = fopen (opt_fname, "r")) == NULL) {
      fprintf (stderr, "fopen failed in %s for %s\n", subname, opt_fname);
      return 1;
   }

   while (fgets (line, MAX_LINE_LEN, fp) != NULL) {
      optname = strtok (line, " \n");
      if ((optval = strtok (NULL, " \n")) == NULL) {
         fprintf (stderr, "unspecified value for %s\n", optname);
         return 1;
      }
      if (strcmp (optname, "day_cnt") == 0) {
         if (parse_to_double (optval, &day_cnt)) {
            fprintf (stderr, "error parsing argument '%s' for option '%s'\n", optval, optname);
            return 1;
         }
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
            fprintf (stderr, "unknown adv_type: %s\n", optval);
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
            fprintf (stderr, "unknown hmix_type: %s\n", optval);
            return 1;
         }
      } else if (strcmp (optname, "vmix_type") == 0) {
         if (strcmp (optval, "const") == 0)
            vmix_opt = vmix_const;
         else if (strcmp (optval, "file") == 0)
            vmix_opt = vmix_file;
         else if (strcmp (optval, "matrix_file") == 0)
            vmix_opt = vmix_matrix_file;
         else {
            fprintf (stderr, "unknown vmix_type: %s\n", optval);
            return 1;
         }
      } else if (strcmp (optname, "coupled_tracer_cnt") == 0) {
         if (parse_to_int (optarg, &coupled_tracer_cnt)) {
            fprintf (stderr, "error parsing argument '%s' for option '%s'\n", optval, optname);
            return 1;
         }
         if (coupled_tracer_cnt != 1) {
            fprintf (stderr, "coupled_tracer_cnt = %d not supported\n", coupled_tracer_cnt);
            return 1;
         }
      } else if (strcmp (optname, "sink_type") == 0) {
         if (strcmp (optval, "none") == 0)
            sink_opt = sink_none;
         else if (strcmp (optval, "const") == 0)
            sink_opt = sink_const;
         else if (strcmp (optval, "const_shallow") == 0)
            sink_opt = sink_const_shallow;
         else if (strcmp (optval, "file") == 0)
            sink_opt = sink_file;
         else if (strcmp (optval, "tracer") == 0)
            sink_opt = sink_tracer;
         else {
            fprintf (stderr, "unknown sink_type: %s\n", optval);
            return 1;
         }
         if ((sink_opt == sink_const) || (sink_opt == sink_const_shallow)) {
            if ((optval = strtok (NULL, " \n")) == NULL) {
               fprintf (stderr, "unspecified sink_rate\n");
               return 1;
            }
            if (parse_to_double (optval, &sink_rate)) {
               fprintf (stderr, "error parsing argument '%s' for option '%s'\n", optval,
                        optname);
               return 1;
            }
            if (sink_opt == sink_const_shallow) {
               if ((optval = strtok (NULL, " \n")) == NULL) {
                  fprintf (stderr, "unspecified sink_depth\n");
                  return 1;
               }
               if (parse_to_double (optval, &sink_depth)) {
                  fprintf (stderr, "error parsing argument '%s' for option '%s'\n", optval,
                           optname);
                  return 1;
               }
            }
         }
         if (sink_opt == sink_file) {
            if ((optval = strtok (NULL, " \n")) == NULL) {
               fprintf (stderr, "unspecified sink_file_name\n");
               return 1;
            }
            if ((sink_file_name = malloc (1 + strlen (optval))) == NULL) {
               fprintf (stderr, "malloc failed in %s for sink_file_name\n", subname);
               return 1;
            }
            strcpy (sink_file_name, optval);
            if ((optval = strtok (NULL, " \n")) == NULL) {
               fprintf (stderr, "unspecified sink_field_name\n");
               return 1;
            }
            if ((sink_field_name = malloc (1 + strlen (optval))) == NULL) {
               fprintf (stderr, "malloc failed in %s for sink_field_name\n", subname);
               return 1;
            }
            strcpy (sink_field_name, optval);
         }
         if (sink_opt == sink_tracer) {
            if ((optval = strtok (NULL, " \n")) == NULL) {
               fprintf (stderr, "unspecified sink_file_name\n");
               return 1;
            }
            if ((sink_file_name = malloc (1 + strlen (optval))) == NULL) {
               fprintf (stderr, "malloc failed in %s for sink_file_name\n", subname);
               return 1;
            }
            strcpy (sink_file_name, optval);
            if ((optval = strtok (NULL, " \n")) == NULL) {
               fprintf (stderr, "unspecified sink_tracer_name\n");
               return 1;
            }
            if ((sink_tracer_name = malloc (1 + strlen (optval))) == NULL) {
               fprintf (stderr, "malloc failed in %s for sink_tracer_name\n", subname);
               return 1;
            }
            strcpy (sink_tracer_name, optval);
            if (strcmp (sink_tracer_name, "Fe")
                && strcmp (sink_tracer_name, "OCMIP_BGC_PO4_DOP")) {
               fprintf (stderr, "unknown tracer name %s for sink_type == sink_tracer in %s\n",
                        sink_tracer_name, subname);
               return 1;
            }
         }
      } else if (strcmp (optname, "reg_fname") == 0) {
         if ((reg_fname = malloc (1 + strlen (optval))) == NULL) {
            fprintf (stderr, "malloc failed in %s for reg_fname\n", subname);
            return 1;
         }
         strcpy (reg_fname, optval);
      } else if (strcmp (optname, "pv") == 0) {
         if ((pv_file_name = malloc (1 + strlen (optval))) == NULL) {
            fprintf (stderr, "malloc failed in %s for pv_file_name\n", subname);
            return 1;
         }
         strcpy (pv_file_name, optval);
         if ((optval = strtok (NULL, " \n")) == NULL) {
            fprintf (stderr, "unspecified pv_field_name\n");
            return 1;
         }
         if ((pv_field_name = malloc (1 + strlen (optval))) == NULL) {
            fprintf (stderr, "malloc failed in %s for pv_field_name\n", subname);
            return 1;
         }
         strcpy (pv_field_name, optval);
      } else if (strcmp (optname, "sf") == 0) {
         if ((d_SF_d_TRACER_file_name = malloc (1 + strlen (optval))) == NULL) {
            fprintf (stderr, "malloc failed in %s for d_SF_d_TRACER_file_name\n", subname);
            return 1;
         }
         strcpy (d_SF_d_TRACER_file_name, optval);
         if ((optval = strtok (NULL, " \n")) == NULL) {
            fprintf (stderr, "unspecified d_SF_d_TRACER_field_name\n");
            return 1;
         }
         if ((d_SF_d_TRACER_field_name = malloc (1 + strlen (optval))) == NULL) {
            fprintf (stderr, "malloc failed in %s for d_SF_d_TRACER_field_name\n", subname);
            return 1;
         }
         strcpy (d_SF_d_TRACER_field_name, optval);
      } else {
         fprintf (stderr, "unknown option name: %s\n", optname);
         return 1;
      }
   }

   fclose (fp);

   return 0;
}

/******************************************************************************/

void
write_opts (void)
{

   if (dbg_lvl) {
      printf ("dbg_lvl                    = %d\n", dbg_lvl);
      printf ("day_cnt                    = %e\n", day_cnt);
      switch (adv_opt) {
      case adv_none:
         printf ("adv_opt                    = %s\n", "none");
         break;
      case adv_donor:
         printf ("adv_opt                    = %s\n", "donor");
         break;
      case adv_cent:
         printf ("adv_opt                    = %s\n", "centered");
         break;
      case adv_upwind3:
         printf ("adv_opt                    = %s\n", "upwind3");
         break;
      }
      switch (hmix_opt) {
      case hmix_none:
         printf ("hmix_opt                   = %s\n", "none");
         break;
      case hmix_const:
         printf ("hmix_opt                   = %s\n", "const");
         break;
      case hmix_hor_file:
         printf ("hmix_opt                   = %s\n", "hor_file");
         break;
      case hmix_isop_file:
         printf ("hmix_opt                   = %s\n", "isop_file");
         break;
      }
      switch (vmix_opt) {
      case vmix_const:
         printf ("vmix_opt                   = %s\n", "const");
         break;
      case vmix_file:
         printf ("vmix_opt                   = %s\n", "file");
         break;
      case vmix_matrix_file:
         printf ("vmix_opt                   = %s\n", "matrix_file");
         break;
      }
      switch (sink_opt) {
      case sink_none:
         printf ("sink_opt                   = %s\n", "none");
         break;
      case sink_const:
         printf ("sink_opt                   = %s\n", "const");
         printf ("sink_rate                  = %e\n", sink_rate);
         break;
      case sink_const_shallow:
         printf ("sink_opt                   = %s\n", "const_shallow");
         printf ("sink_rate                  = %e\n", sink_rate);
         printf ("sink_depth                 = %e\n", sink_depth);
         break;
      case sink_file:
         printf ("sink_opt                   = %s\n", "file");
         printf ("sink_file_name             = %s\n", sink_file_name);
         printf ("sink_field_name            = %s\n", sink_field_name);
         break;
      case sink_tracer:
         printf ("sink_opt                   = %s\n", "tracer");
         printf ("sink_file_name             = %s\n", sink_file_name);
         printf ("sink_tracer_name           = %s\n", sink_tracer_name);
         break;
      }
      printf ("reg_fname                  = %s\n", reg_fname ? reg_fname : "none");
      printf ("pv_file_name               = %s\n", pv_file_name ? pv_file_name : "none");
      printf ("pv_field_name              = %s\n", pv_field_name ? pv_field_name : "none");
      printf ("d_SF_d_TRACER_file_name    = %s\n",
              d_SF_d_TRACER_file_name ? d_SF_d_TRACER_file_name : "none");
      printf ("d_SF_d_TRACER_field_name   = %s\n",
              d_SF_d_TRACER_field_name ? d_SF_d_TRACER_field_name : "none");
      printf ("circ_fname                 = %s\n", circ_fname);
      printf ("matrix_fname               = %s\n\n", matrix_fname);
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

   set_opt_defaults ();

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

   exit (EXIT_SUCCESS);
}
