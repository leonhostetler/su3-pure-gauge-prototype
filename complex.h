/**************************************************************
* complex.h -- This is the header file for the complex module
*	which deals with complex numbers.
* 
* The constants are wrapped in a preprocessor directive that
* ensures that the contents of this header file are included
* by the compiler only once, even if perhaps nested modules
* in a modular program instruct the compiler to include
* this header file.
*
* Author: Leon Hostetler, Jul 16, 2018
*
**************************************************************/
#ifndef _COMPLEX_H_INCLUDED_

/* Structures */
struct complex_struct {
    double real;	/* The real part */
    double imag;	/* The imaginary part */
};
typedef struct complex_struct complex;

/* Function Prototypes */
/***********************************************************
cadd -- A function that adds two complex numbers

parameters
  num1, num2 -- The two complex numbers to add

returns
  the sum as a complex struct
***********************************************************/
complex cadd(complex num1, complex num2);

/***********************************************************
cmul -- A function that multiplies two complex numbers

parameters
  num1, num2 -- The two complex numbers to multiply

returns
  the product as a complex struct
***********************************************************/
complex cmul(complex num1, complex num2);

/***********************************************************
cpri -- A function that prints a complex number to the screen
  May be useful for debugging or verbose output

parameters
  num -- The complex number to print

returns
  void
***********************************************************/
void cpri(complex num);

/***********************************************************
ccon -- A function that takes a given complex number and
	returns its complex conjugate.

parameters
  num -- The complex number

returns
  the complex conjugate
***********************************************************/
complex ccon(complex num);

#define _COMPLEX_H_INCLUDED_
#endif /* _COMPLEX_H_INCLUDED_ */

