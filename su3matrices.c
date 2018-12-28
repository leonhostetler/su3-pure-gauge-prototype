/**************************************************************
* su3matrices.c -- This is a module containing functions for
*	working with complex 3x3 SU(3) matrices.
*
* Routines
*
*     su3_pri -- Prints the SU(3) matrices
*     su3_add -- Adds two SU(3) matrices
*     su3_sub -- Subtracts two SU(3) matrices
*     su3_mul -- Multiplies two SU(3) matrices
*     su3_hco -- Takes the Hermitian conjugate of an SU(3) matrix
*     su3_rtr -- Takes the (real) trace of an SU(3) matrix
*
* Author: Leon Hostetler, Jul 16, 2018
*
**************************************************************/
#include "su3matrices.h"
#include "complex.h"

/***********************************************************
su3_pri -- A function that prints an SU(3) array to the screen.
	May be useful for debugging or verbose output.

parameters
  su3matrix -- the array to print

returns
  void
***********************************************************/
void su3_pri(su3 su3matrix)
{
    int i, j;
    for(i=0; i<3; ++i) {
        for(j=0; j<3; ++j) {
            cpri(su3matrix.a[i][j]);
        }
        printf("\n");
    }
}


/***********************************************************
su3_add -- A function that adds two SU(3) matrices.

parameters
  arr1, arr2 -- the arrays to add
  arrsum -- the array to put the sum into

returns
  void
***********************************************************/
su3 su3_add(su3 arr1, su3 arr2)
{
    
    su3 arrsum;

    /* Add the elements */
    int i, j;
    for(i=0; i<3; ++i) {
        for(j=0; j<3; ++j) {
            arrsum.a[i][j] = cadd(arr1.a[i][j],arr2.a[i][j]);
        }
    }
    return (arrsum);
}


/***********************************************************
su3_sub -- A function that subtracts two SU(3) matrices.

parameters
  arr1, arr2 -- the matrices to subtract (2nd is subtracted from first)
  arrsub -- the array to put the difference into

returns
  an SU(3) matrix
***********************************************************/
su3 su3_sub(su3 arr1, su3 arr2)
{
    
    su3 arrsub;
    complex neg1 = {-1,0};

    /* Subtract the elements */
    int i, j;
    for(i=0; i<3; ++i) {
        for(j=0; j<3; ++j) {
            arrsub.a[i][j] = cadd(arr1.a[i][j],cmul(neg1,arr2.a[i][j]));
        }
    }
    return (arrsub);
}


/***********************************************************
su3_mult -- A function that multiplies two SU(3) matrices.

parameters
  arr1, arr2 -- the matrices to multiply
  arrprod -- the matrix to put the product into

returns
  void
***********************************************************/
su3 su3_mul(su3 arr1, su3 arr2)
{
    su3 arrprod;
    complex sum;
    int i,j,k;

    /* Multiply the elements */
    for(i=0; i<3; ++i) {
        for(j=0; j<3; ++j) {
            sum.real = 0.0;
            sum.imag = 0.0;
            for(k=0; k<3; ++k) {
                sum = cadd(sum,cmul(arr1.a[i][k],arr2.a[k][j]));
            }
            arrprod.a[i][j] = sum;
        }
    }
    return (arrprod);
}


/***********************************************************
su3_hco -- A function that returns the Hermitian conjugate of
	an  SU(3) matrix

parameters
  arr -- the SU(3) matrix to take the Hermitian conjugate of

returns
  the Hermitian conjugate matrix
***********************************************************/
su3 su3_hco(su3 arr)
{
    su3 arrconj;
    int i,j;

    /* Transpose and conjugate */
    for(i=0; i<3; ++i) {
        for(j=0; j<3; ++j) {
            arrconj.a[j][i] = ccon(arr.a[i][j]);
        }
    }
    return (arrconj);
}


/***********************************************************
su3_rtr -- A function that returns the real part of the trace
	of an SU(3) matrix

parameters
  arr -- the SU(3) matrix to take the real trace of

returns
  the real part of the trace
***********************************************************/
double su3_rtr(su3 arr)
{
    double real_trace;

    real_trace = arr.a[0][0].real + arr.a[1][1].real + arr.a[2][2].real;

    return (real_trace);
}
