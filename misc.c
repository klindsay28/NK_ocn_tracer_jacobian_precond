#include "misc.h"

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
      fprintf (stderr, "%s:nothing to parse\n", subname);
      return 1;
   }

   errno = 0;

   *val = strtol (str, &cp, 0);

   if (errno == ERANGE) {
      fprintf (stderr, "%s:ERANGE error parsing '%s'\n", subname, str);
      return 1;
   }

   if (*cp != '\0') {
      fprintf (stderr, "%s:unexpected character '%c' parsing '%s'\n", subname, *cp, str);
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
      fprintf (stderr, "%s:nothing to parse\n", subname);
      return 1;
   }

   if (parse_to_long (str, &lval)) {
      fprintf (stderr, "%s:error from parse_to_long\n", subname);
      return 1;
   }

   if ((lval > INT_MAX) || (lval < INT_MIN)) {
      fprintf (stderr, "%s:value %ld out of int range\n", subname, lval);
      return 1;
   }

   *val = (int) lval;

   return 0;
}
