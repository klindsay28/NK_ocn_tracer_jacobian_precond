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

/******************************************************************************/

int dbg_lvl;
int iam;

/******************************************************************************/

int
parse_cmd_line (int argc, char **argv, double *day_cnt, char **matrix_fname)
{
   char *subname = "parse_cmd_line";
   char *usage_msg =
      "usage: gen_matrix_file [-D dbg_lvl] [-c day_cnt] [-o adv_type,[none|donor|centered|upwind3]] [-o hmix_type,[none|const|hor_file|isop_file]] [-o vmix_type,[const|file|matrix_file]] [-o sink_type,[none|const|const_shallow|file|tracer][,rate[,depth]|,file_name,[field_name|tracer_name]]] [-o reg_fname,reg_fname] [-p pv_file_name,pv_field_name] [-s d_SF_d_TRACER_file_name,d_SF_d_TRACER_field_name] circ_fname matrix_fname";
   extern char *optarg;
   extern int optind;
   char *optstring = "D:c:o:p:s:h";
   int opt;
   char *cp;

   while ((opt = getopt (argc, argv, optstring)) != -1) {
      switch (opt) {
      case 'D':
         dbg_lvl = atoi (optarg);
         break;
      case 'c':
         *day_cnt = atof (optarg);
         break;
      case 'o':
         cp = strtok (optarg, ",");
         if (strcmp (cp, "adv_type") == 0) {
            if ((cp = strtok (NULL, ",")) == NULL) {
               fprintf (stderr, "unspecified advection_type\n");
               return 1;
            }
            if (strcmp (cp, "none") == 0)
               adv_opt = adv_none;
            else if (strcmp (cp, "donor") == 0)
               adv_opt = adv_donor;
            else if (strncmp (cp, "centered", 4) == 0)
               adv_opt = adv_cent;
            else if (strcmp (cp, "upwind3") == 0)
               adv_opt = adv_upwind3;
            else {
               fprintf (stderr, "unknown advection_type: %s\n", cp);
               return 1;
            }
         } else if (strcmp (cp, "hmix_type") == 0) {
            if ((cp = strtok (NULL, ",")) == NULL) {
               fprintf (stderr, "unspecified hmix_type\n");
               return 1;
            }
            if (strcmp (cp, "none") == 0)
               hmix_opt = hmix_none;
            else if (strcmp (cp, "const") == 0)
               hmix_opt = hmix_const;
            else if (strcmp (cp, "hor_file") == 0)
               hmix_opt = hmix_hor_file;
            else if (strcmp (cp, "isop_file") == 0)
               hmix_opt = hmix_isop_file;
            else {
               fprintf (stderr, "unknown hmix_type: %s\n", cp);
               return 1;
            }
         } else if (strcmp (cp, "vmix_type") == 0) {
            if ((cp = strtok (NULL, ",")) == NULL) {
               fprintf (stderr, "unspecified vmix_type\n");
               return 1;
            }
            if (strcmp (cp, "const") == 0)
               vmix_opt = vmix_const;
            else if (strcmp (cp, "file") == 0)
               vmix_opt = vmix_file;
            else if (strcmp (cp, "matrix_file") == 0)
               vmix_opt = vmix_matrix_file;
            else {
               fprintf (stderr, "unknown vmix_type: %s\n", cp);
               return 1;
            }
         } else if (strcmp (cp, "sink_type") == 0) {
            if ((cp = strtok (NULL, ",")) == NULL) {
               fprintf (stderr, "unspecified sink_type\n");
               return 1;
            }
            if (strcmp (cp, "none") == 0)
               sink_opt = sink_none;
            else if (strcmp (cp, "const") == 0)
               sink_opt = sink_const;
            else if (strcmp (cp, "const_shallow") == 0)
               sink_opt = sink_const_shallow;
            else if (strcmp (cp, "file") == 0)
               sink_opt = sink_file;
            else if (strcmp (cp, "tracer") == 0)
               sink_opt = sink_tracer;
            else {
               fprintf (stderr, "unknown sink_type: %s\n", cp);
               return 1;
            }
            if ((sink_opt == sink_const) || (sink_opt == sink_const_shallow)) {
               if ((cp = strtok (NULL, ",")) == NULL) {
                  fprintf (stderr, "unspecified sink_rate\n");
                  return 1;
               }
               sink_rate = atof (cp);
               if (sink_opt == sink_const_shallow) {
                  if ((cp = strtok (NULL, ",")) == NULL) {
                     fprintf (stderr, "unspecified sink_depth\n");
                     return 1;
                  }
                  sink_depth = atof (cp);
               }
            }
            if (sink_opt == sink_file) {
               if ((cp = strtok (NULL, ",")) == NULL) {
                  fprintf (stderr, "unspecified sink_file_name\n");
                  return 1;
               }
               if ((sink_file_name = malloc (1 + strlen (cp))) == NULL) {
                  fprintf (stderr, "malloc failed in %s for sink_file_name\n", subname);
                  return 1;
               }
               strcpy (sink_file_name, cp);
               if ((cp = strtok (NULL, ",")) == NULL) {
                  fprintf (stderr, "unspecified sink_field_name\n");
                  return 1;
               }
               if ((sink_field_name = malloc (1 + strlen (cp))) == NULL) {
                  fprintf (stderr, "malloc failed in %s for sink_field_name\n", subname);
                  return 1;
               }
               strcpy (sink_field_name, cp);
            }
            if (sink_opt == sink_tracer) {
               if ((cp = strtok (NULL, ",")) == NULL) {
                  fprintf (stderr, "unspecified sink_file_name\n");
                  return 1;
               }
               if ((sink_file_name = malloc (1 + strlen (cp))) == NULL) {
                  fprintf (stderr, "malloc failed in %s for sink_file_name\n", subname);
                  return 1;
               }
               strcpy (sink_file_name, cp);
               if ((cp = strtok (NULL, ",")) == NULL) {
                  fprintf (stderr, "unspecified sink_tracer_name\n");
                  return 1;
               }
               if ((sink_tracer_name = malloc (1 + strlen (cp))) == NULL) {
                  fprintf (stderr, "malloc failed in %s for sink_tracer_name\n", subname);
                  return 1;
               }
               strcpy (sink_tracer_name, cp);
               if (strcmp (sink_tracer_name, "Fe")) {
                  fprintf (stderr,
                           "unknown tracer name %s for sink_type == sink_tracer in %s\n",
                           sink_tracer_name, subname);
                  return 1;
               }
            }
         } else if (strcmp (cp, "reg_fname") == 0) {
            if ((cp = strtok (NULL, ",")) == NULL) {
               fprintf (stderr, "unspecified reg_fname\n");
               return 1;
            }
            if ((reg_fname = malloc (1 + strlen (cp))) == NULL) {
               fprintf (stderr, "malloc failed in %s for reg_fname\n", subname);
               return 1;
            }
            strcpy (reg_fname, cp);
         } else {
            fprintf (stderr, "unknown option string: %s\n", cp);
            return 1;
         }
         break;
      case 'p':
         cp = strtok (optarg, ",");
         if ((pv_file_name = malloc (1 + strlen (cp))) == NULL) {
            fprintf (stderr, "malloc failed in %s for pv_file_name\n", subname);
            return 1;
         }
         strcpy (pv_file_name, cp);
         if ((cp = strtok (NULL, ",")) == NULL) {
            fprintf (stderr, "unspecified pv_field_name\n");
            return 1;
         }
         if ((pv_field_name = malloc (1 + strlen (cp))) == NULL) {
            fprintf (stderr, "malloc failed in %s for pv_field_name\n", subname);
            return 1;
         }
         strcpy (pv_field_name, cp);
         break;
      case 's':
         cp = strtok (optarg, ",");
         if ((d_SF_d_TRACER_file_name = malloc (1 + strlen (cp))) == NULL) {
            fprintf (stderr, "malloc failed in %s for d_SF_d_TRACER_file_name\n", subname);
            return 1;
         }
         strcpy (d_SF_d_TRACER_file_name, cp);
         if ((cp = strtok (NULL, ",")) == NULL) {
            fprintf (stderr, "unspecified d_SF_d_TRACER_field_name\n");
            return 1;
         }
         if ((d_SF_d_TRACER_field_name = malloc (1 + strlen (cp))) == NULL) {
            fprintf (stderr, "malloc failed in %s for d_SF_d_TRACER_field_name\n", subname);
            return 1;
         }
         strcpy (d_SF_d_TRACER_field_name, cp);
         break;
      case 'h':
      case '?':
         fprintf (stderr, "%s\n", usage_msg);
         return 1;
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
   *matrix_fname = argv[optind++];
   return 0;
}

/******************************************************************************/

int
main (int argc, char *argv[])
{
   double day_cnt;
   char *matrix_fname = NULL;

   int i;
   int j;
   int k;

   iam = 0;

   /* set defaults */
   dbg_lvl = 0;
   day_cnt = 365.0;
   adv_opt = adv_cent;
   hmix_opt = hmix_isop_file;
   vmix_opt = vmix_file;
   sink_opt = sink_none;
   sink_rate = 1.0e-5;
   sink_depth = 10.0e2;

   if (parse_cmd_line (argc, argv, &day_cnt, &matrix_fname))
      exit (EXIT_FAILURE);

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

   if (get_grid_info (circ_fname, reg_fname))
      exit (EXIT_FAILURE);

   if (put_grid_info (matrix_fname))
      exit (EXIT_FAILURE);

   if (comp_sparse_matrix_size ())
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
