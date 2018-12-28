/**************************************************************
* su3matrices.h -- The header file for the su3matrices module
*	which contains functions for working with SU(3) matrices.
*
* Author: Leon Hostetler, Jul 16, 2018
*
**************************************************************/
#include "complex.h"
#include <stdio.h>

/* Structures */

/*************************************************************
One benefit of using an su3 struct instead of just a complex
matrix, is that structs can be passed into and returned from
functions.
*************************************************************/
struct su3_struct {
    complex a[3][3];
};
typedef struct su3_struct su3;

/* Function Prototypes */

/***********************************************************
su3_pri -- A function that prints an SU(3) array to the screen.
	May be useful for debugging or verbose output.

parameters
  array -- the array to print

returns
  void
***********************************************************/
void su3_pri(su3 su3matrix); 

/***********************************************************
su3_add -- A function that adds two SU(3) matrices.

parameters
  arr1, arr2 -- the arrays to add
  arrsum -- the array to put the sum into

returns
  void
***********************************************************/
su3 su3_add(su3 arr1, su3 arr2);

/***********************************************************
su3_sub -- A function that subtracts two SU(3) matrices.

parameters
  arr1, arr2 -- the matrices to subtract (2nd is subtracted from first)
  arrsub -- the array to put the difference into

returns
  an SU(3) matrix
***********************************************************/
su3 su3_sub(su3 arr1, su3 arr2);

/***********************************************************
su3_mult -- A function that multiplies two SU(3) matrices.

parameters
  arr1, arr2 -- the matrices to multiply
  arrprod -- the matrix to put the product into

returns
  void
***********************************************************/
/*void su3_mul(su3 arr1, su3 arr2, su3 arrprod);*/
su3 su3_mul(su3 arr1, su3 arr2);

/***********************************************************
su3_hco -- A function that returns the Hermitian conjugate of
	an  SU(3) matrix

parameters
  arr -- the SU(3) matrix to take the Hermitian conjugate of

returns
  the Hermitian conjugate matrix
***********************************************************/
su3 su3_hco(su3 arr);

/***********************************************************
su3_rtr -- A function that returns the real part of the trace
	of an SU(3) matrix

parameters
  arr -- the SU(3) matrix to take the real trace of

returns
  the real part of the trace
***********************************************************/
double su3_rtr(su3 arr);
