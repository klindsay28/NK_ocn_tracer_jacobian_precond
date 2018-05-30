#include "misc.h"
#include "globals.h"

#include <errno.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>

/******************************************************************************/

int
parse_to_long (char *str, long *val)
{
   char *subname = "parse_to_long";
   char *cp;

   if ((str == NULL) || (*str == '\0')) {
      fprintf (stderr, "(%d) %s:nothing to parse\n", iam, subname);
      return 1;
   }

   errno = 0;

   *val = strtol (str, &cp, 0);

   if (errno == ERANGE) {
      fprintf (stderr, "(%d) %s:ERANGE error parsing '%s'\n", iam, subname, str);
      return 1;
   }

   if (*cp != '\0') {
      fprintf (stderr, "(%d) %s:unexpected character '%c' parsing '%s'\n", iam, subname, *cp, str);
      return 1;
   }

   return 0;
}

/******************************************************************************/

int
parse_to_int (char *str, int *val)
{
   char *subname = "parse_to_int";
   long lval;

   if ((str == NULL) || (*str == '\0')) {
      fprintf (stderr, "(%d) %s:nothing to parse\n", iam, subname);
      return 1;
   }

   if (parse_to_long (str, &lval)) {
      fprintf (stderr, "(%d) %s:error from parse_to_long\n", iam, subname);
      return 1;
   }

   if ((lval > INT_MAX) || (lval < INT_MIN)) {
      fprintf (stderr, "(%d) %s:value %ld out of int range\n", iam, subname, lval);
      return 1;
   }

   *val = (int) lval;

   return 0;
}

/******************************************************************************/

int
parse_to_double (char *str, double *val)
{
   char *subname = "parse_to_double";
   char *cp;

   if ((str == NULL) || (*str == '\0')) {
      fprintf (stderr, "(%d) %s:nothing to parse\n", iam, subname);
      return 1;
   }

   errno = 0;

   *val = strtod (str, &cp);

   if (errno == ERANGE) {
      fprintf (stderr, "(%d) %s:ERANGE error parsing '%s'\n", iam, subname, str);
      return 1;
   }

   if (*cp != '\0') {
      fprintf (stderr, "(%d) %s:unexpected character '%c' parsing '%s'\n", iam, subname, *cp, str);
      return 1;
   }

   return 0;
}
