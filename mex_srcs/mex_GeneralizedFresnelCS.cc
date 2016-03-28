/****************************************************************************\
  Copyright (c) Enrico Bertolazzi 2014
  All Rights Reserved.

  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License
  as published by the Free Software Foundation;

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
\****************************************************************************/

#include "Clothoid.hh"
#include "mex.h"

#define MEX_ERROR_MESSAGE \
"%===================================================================%\n" \
"%  Compute Fresnel sine and cosine integrals momenta                %\n" \
"%                                                                   %\n" \
"%  USAGE: [X,Y] = GeneralizedFresnelCS( nk, a, b, c ) ;             %\n" \
"%                                                                   %\n" \
"%  Integrals are defined as:                                        %\n" \
"%                                                                   %\n" \
"%  X_k(a,b,c) = int_0^1 t^k * cos( (a/2)*t^2 + b*t + c ) dt         %\n" \
"%  Y_k(a,b,c) = int_0^1 t^k * sin( (a/2)*t^2 + b*t + c ) dt         %\n" \
"%                                                                   %\n" \
"%  On input:                                                        %\n" \
"%                                                                   %\n" \
"%    nk      = number of momentae to be computed (1 <= nk <= 3)     %\n" \
"%    a, b, c = the parameters of the integrals                      %\n" \
"%                                                                   %\n" \
"%  On output:                                                       %\n" \
"%                                                                   %\n" \
"%    X = vector with Fresnel cosine momenta [X_0,X_1,...,X_{nk-1}]  %\n" \
"%    Y = vector with Fresnel sine momenta   [Y_0,Y_1,...,Y_{nk-1}]  %\n" \
"%                                                                   %\n" \
"%===================================================================%\n" \
"%                                                                   %\n" \
"%  Autor: Enrico Bertolazzi                                         %\n" \
"%         Department of Industrial Engineering                      %\n" \
"%         University of Trento                                      %\n" \
"%         enrico.bertolazzi@unitn.it                                %\n" \
"%                                                                   %\n" \
"%===================================================================%\n"


#define arg_nk prhs[0]
#define arg_a  prhs[1]
#define arg_b  prhs[2]
#define arg_c  prhs[3]
#define arg_C  plhs[0]
#define arg_S  plhs[1]

extern "C"
void
mexFunction( int nlhs, mxArray       *plhs[],
             int nrhs, mxArray const *prhs[] ) {

  // Check for proper number of arguments, etc
  if ( nrhs != 4 ) {
	  mexErrMsgTxt(MEX_ERROR_MESSAGE) ;
    return ;
  } else if ( mxGetClassID(arg_a) != mxDOUBLE_CLASS ||
              mxGetClassID(arg_b) != mxDOUBLE_CLASS ||
              mxGetClassID(arg_c) != mxDOUBLE_CLASS ) {
	  mexErrMsgTxt("Input arguments should be double");
    return ;
  } else if ( mxIsComplex(arg_nk) ||
              mxIsComplex(arg_a)  ||
              mxIsComplex(arg_b)  ||
              mxIsComplex(arg_c) ) {
	  mexErrMsgTxt("Input arguments should be real (not complex)");
    return ;
  }

  for ( int kk = 0 ; kk < 4 ; ++kk )
    if ( mxGetM(prhs[kk]) != 1 || mxGetN(prhs[kk]) != 1 )
	    mexErrMsgTxt("Input argument must be a scalar");

  // Output array
  int nk = int(mxGetScalar(arg_nk)) ;
  if ( nk < 1 || nk > 3 )
    mexErrMsgTxt("First argument must be an integer in the range [1..3]");
  
  Clothoid::valueType a = mxGetScalar(arg_a) ;
  Clothoid::valueType b = mxGetScalar(arg_b) ;
  Clothoid::valueType c = mxGetScalar(arg_c) ;

  if ( nk == 1 ) {
    Clothoid::valueType C, S ;
    Clothoid::GeneralizedFresnelCS( a, b, c, C, S ) ;
    arg_C = mxCreateDoubleScalar(C) ;
    arg_S = mxCreateDoubleScalar(S) ;
  } else {
    arg_C = mxCreateDoubleMatrix(1,nk,mxREAL) ;
    arg_S = mxCreateDoubleMatrix(1,nk,mxREAL) ;
    double * C = mxGetPr(arg_C) ;
	  double * S = mxGetPr(arg_S) ;
    Clothoid::GeneralizedFresnelCS( nk, a, b, c, C, S ) ;
  }
}
