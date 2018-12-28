/**************************************************************
* complex.c -- This is a module containing useful functions
* 	functions for working with complex numbers
*
* Routines
*	cadd -- Add two complex numbers
*	cmul -- Multiply two complex numbers
*       cpri -- Print a complex number
*       ccon -- Complex conjugate of a complex number
* 
* Author: Leon Hostetler, Jul 16, 2018
*
**************************************************************/
#include "complex.h"
#include <stdio.h>

/***********************************************************
cadd -- A function that adds two complex numbers

parameters
  num1, num2 -- The two complex numbers to add

returns
  the sum as a complex struct
***********************************************************/
complex cadd(complex num1, complex num2)
{
    
    complex sum;

    sum.real = num1.real + num2.real;
    sum.imag = num1.imag + num2.imag;

    return (sum);
}

/***********************************************************
cmul -- A function that multiplies two complex numbers

parameters
  num1, num2 -- The two complex numbers to multiply

returns
  the product as a complex struct
***********************************************************/
complex cmul(complex num1, complex num2)
{
    
    complex product;

    product.real = num1.real*num2.real - num1.imag*num2.imag;
    product.imag = num1.real*num2.imag + num2.real*num1.imag;

    return (product);
}

/***********************************************************
cpri -- A function that prints a complex number to the screen
  May be useful for debugging or verbose output

parameters
  num -- The complex number to print

returns
  void
***********************************************************/
void cpri(complex num)
{
    printf(" %lf%+lfi ", num.real, num.imag);
}


/***********************************************************
ccon -- A function that takes a given complex number and
	returns its complex conjugate.

parameters
  num -- The complex number

returns
  the complex conjugate
***********************************************************/
complex ccon(complex num)
{
    complex conj;

    conj.real = num.real;
    conj.imag = -num.imag;

    return (conj);
}
