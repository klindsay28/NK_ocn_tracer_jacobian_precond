#include "memory.h"
#include "globals.h"

#include <stdio.h>
#include <stdlib.h>

/******************************************************************************/

/* allocate a 2D integer array */
/* the underlying memory for the 2D array is contiguous */
int **
malloc_2d_int (int jmt, int imt)
{
   char *subname = "malloc_2d_int";
   int **retval;
   int j;

   if ((retval = malloc ((size_t) jmt * sizeof (int *))) == NULL) {
      fprintf (stderr, "(%d) malloc failure for retval in %s\n", iam, subname);
      return NULL;
   }
   if ((retval[0] = malloc ((size_t) jmt * (size_t) imt * sizeof (int))) == NULL) {
      fprintf (stderr, "(%d) malloc failure for retval[0] in %s\n", iam, subname);
      return NULL;
   }
   for (j = 0; j < jmt; j++)
      retval[j] = retval[0] + j * imt;
   return retval;
}

/******************************************************************************/

void
free_2d_int (int **ptr)
{
   free (ptr[0]);
   free (ptr);
}

/******************************************************************************/

/* allocate a 3D integer array */
/* the underlying memory for the 3D array is contiguous */
int ***
malloc_3d_int (int km, int jmt, int imt)
{
   char *subname = "malloc_3d_int";
   int ***retval;
   int j;
   int k;

   if ((retval = malloc ((size_t) km * sizeof (int **))) == NULL) {
      fprintf (stderr, "(%d) malloc failure for retval in %s\n", iam, subname);
      return NULL;
   }
   if ((retval[0] = malloc ((size_t) km * (size_t) jmt * sizeof (int *))) == NULL) {
      fprintf (stderr, "(%d) malloc failure for retval[0] in %s\n", iam, subname);
      return NULL;
   }
   if ((retval[0][0] = malloc ((size_t) km * (size_t) jmt * (size_t) imt * sizeof (int))) == NULL) {
      fprintf (stderr, "(%d) malloc failure for retval[0][0] in %s\n", iam, subname);
      return NULL;
   }
   for (k = 0; k < km; k++) {
      retval[k] = retval[0] + k * jmt;
      for (j = 0; j < jmt; j++)
         retval[k][j] = retval[0][0] + k * jmt * imt + j * imt;
   }
   return retval;
}

/******************************************************************************/

void
free_3d_int (int ***ptr)
{
   free (ptr[0][0]);
   free (ptr[0]);
   free (ptr);
}

/******************************************************************************/

/* allocate a 2D double array */
/* the underlying memory for the 2D array is contiguous */
double **
malloc_2d_double (int jmt, int imt)
{
   char *subname = "malloc_2d_double";
   double **retval;
   int j;

   if ((retval = malloc ((size_t) jmt * sizeof (double *))) == NULL) {
      fprintf (stderr, "(%d) malloc failure for retval in %s\n", iam, subname);
      return NULL;
   }
   if ((retval[0] = malloc ((size_t) jmt * (size_t) imt * sizeof (double))) == NULL) {
      fprintf (stderr, "(%d) malloc failure for retval[0] in %s\n", iam, subname);
      return NULL;
   }
   for (j = 0; j < jmt; j++)
      retval[j] = retval[0] + j * imt;
   return retval;
}

/******************************************************************************/

void
free_2d_double (double **ptr)
{
   free (ptr[0]);
   free (ptr);
}

/******************************************************************************/

/* allocate a 3D double array */
/* the underlying memory for the 3D array is contiguous */
double ***
malloc_3d_double (int km, int jmt, int imt)
{
   char *subname = "malloc_3d_double";
   double ***retval;
   int j;
   int k;

   if ((retval = malloc ((size_t) km * sizeof (double **))) == NULL) {
      fprintf (stderr, "(%d) malloc failure for retval in %s\n", iam, subname);
      return NULL;
   }
   if ((retval[0] = malloc ((size_t) km * (size_t) jmt * sizeof (double *))) == NULL) {
      fprintf (stderr, "(%d) malloc failure for retval[0] in %s\n", iam, subname);
      return NULL;
   }
   if ((retval[0][0] = malloc ((size_t) km * (size_t) jmt * (size_t) imt * sizeof (double))) == NULL) {
      fprintf (stderr, "(%d) malloc failure for retval[0][0] in %s\n", iam, subname);
      return NULL;
   }
   for (k = 0; k < km; k++) {
      retval[k] = retval[0] + k * jmt;
      for (j = 0; j < jmt; j++)
         retval[k][j] = retval[0][0] + k * jmt * imt + j * imt;
   }
   return retval;
}

/******************************************************************************/

void
free_3d_double (double ***ptr)
{
   free (ptr[0][0]);
   free (ptr[0]);
   free (ptr);
}
