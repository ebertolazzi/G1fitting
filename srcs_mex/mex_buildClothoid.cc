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

#include <sstream>
#include <stdexcept>

#define MEX_ERROR_MESSAGE \
"%=============================================================================%\n" \
"%  buildClothoid:  Compute parameters of the G1 Hermite clothoid fitting      %\n" \
"%                                                                             %\n" \
"%  USAGE: [k,dk,L,iter] = buildClothoid( x0, y0, theta0, x1, y1, theta1 ) ;   %\n" \
"%         [k,dk,L,iter,k_1,dk_1,L_1,k_2,dk_2,L_2] = ...                       %\n" \
"%                         buildClothoid( x0, y0, theta0, x1, y1, theta1 ) ;   %\n" \
"%                                                                             %\n" \
"%  On input:                                                                  %\n" \
"%                                                                             %\n" \
"%       x0, y0  = coodinate of initial point                                  %\n" \
"%       theta0  = orientation (angle) of the clothoid at initial point        %\n" \
"%       x1, y1  = coodinate of final point                                    %\n" \
"%       theta1  = orientation (angle) of the clothoid at final point          %\n" \
"%                                                                             %\n" \
"%  On output:                                                                 %\n" \
"%                                                                             %\n" \
"%       L  = the lenght of the clothoid curve from initial to final point     %\n" \
"%       k  = curvature at initial point                                       %\n" \
"%       dk = derivative of curvature respect to arclength,                    %\n" \
"%            notice that curvature at final point is k+dk*L                   %\n" \
"%       iter = Newton Iterations used to solve the interpolation problem      %\n" \
"%                                                                             %\n" \
"%  optional output                                                            %\n" \
"%                                                                             %\n" \
"%       k_1  = partial derivative of the solution respect to theta0           %\n" \
"%       dk_1 = partial derivative of the solution respect to theta0           %\n" \
"%       L_1  = partial derivative of the solution respect to theta0           %\n" \
"%       k_2  = partial derivative of the solution respect to theta1           %\n" \
"%       dk_2 = partial derivative of the solution respect to theta1           %\n" \
"%       L_2  = partial derivative of the solution respect to theta1           %\n" \
"%                                                                             %\n" \
"%=============================================================================%\n" \
"%                                                                             %\n" \
"%  Autors: Enrico Bertolazzi and Marco Frego                                  %\n" \
"%          Department of Industrial Engineering                               %\n" \
"%          University of Trento                                               %\n" \
"%          enrico.bertolazzi@unitn.it                                         %\n" \
"%          m.fregox@gmail.com                                                 %\n" \
"%                                                                             %\n" \
"%=============================================================================%\n"

#define ASSERT(COND,MSG)                      \
  if ( !(COND) ) {                            \
    std::ostringstream ost ;                  \
    ost << "buildClothoid: " << MSG << '\n' ; \
    mexErrMsgTxt(ost.str().c_str()) ;         \
  }

#define arg_x0     prhs[0]
#define arg_y0     prhs[1]
#define arg_theta0 prhs[2]
#define arg_x1     prhs[3]
#define arg_y1     prhs[4]
#define arg_theta1 prhs[5]

#define arg_k      plhs[0]
#define arg_dk     plhs[1]
#define arg_L      plhs[2]
#define arg_iter   plhs[3]
#define arg_k_1    plhs[4]
#define arg_dk_1   plhs[5]
#define arg_L_1    plhs[6]
#define arg_k_2    plhs[7]
#define arg_dk_2   plhs[8]
#define arg_L_2    plhs[9]

extern "C"
void
mexFunction( int nlhs, mxArray       *plhs[],
             int nrhs, mxArray const *prhs[] ) {

  // Check for proper number of arguments, etc
  if ( nrhs < 6 ) {
	  mexErrMsgTxt(MEX_ERROR_MESSAGE) ;
    return ;
  } else if ( mxGetClassID(arg_x0)     != mxDOUBLE_CLASS ||
              mxGetClassID(arg_y0)     != mxDOUBLE_CLASS ||
              mxGetClassID(arg_theta0) != mxDOUBLE_CLASS ||
              mxGetClassID(arg_x1)     != mxDOUBLE_CLASS ||
              mxGetClassID(arg_y1)     != mxDOUBLE_CLASS ||
              mxGetClassID(arg_theta1) != mxDOUBLE_CLASS ) {
	  mexErrMsgTxt("Input arguments should be double");
  } else if ( mxIsComplex(arg_x0)     ||
              mxIsComplex(arg_y0)     ||
              mxIsComplex(arg_theta0) ||
              mxIsComplex(arg_x1)     ||
              mxIsComplex(arg_y1)     ||
              mxIsComplex(arg_theta1) ) {
	  mexErrMsgTxt("Input arguments should be real (not complex)");
  }

  for ( int kk = 0 ; kk < 6 ; ++kk )
    if ( mxGetM(prhs[kk]) != 1 || mxGetN(prhs[kk]) != 1 )
	    mexErrMsgTxt("Input arguments must be scalars");

  ASSERT( nlhs == 3 || nlhs == 4 || nlhs == 10,
          "wrong number of outout arguments\n"
          "expected 3 or 4 or 10, found " << nlhs ) ;

  int iter ;
  if ( nlhs == 10 ) {
    Clothoid::valueType k, dk, L, k_1, dk_1, L_1, k_2, dk_2, L_2 ;
    iter = Clothoid::buildClothoid( mxGetScalar(arg_x0),
                                    mxGetScalar(arg_y0),
                                    mxGetScalar(arg_theta0),
                                    mxGetScalar(arg_x1),
                                    mxGetScalar(arg_y1),
                                    mxGetScalar(arg_theta1),
                                    k, dk, L, k_1, dk_1, L_1, k_2, dk_2, L_2 ) ;

    arg_k    = mxCreateDoubleScalar(k) ;
    arg_dk   = mxCreateDoubleScalar(dk) ;
    arg_L    = mxCreateDoubleScalar(L) ;

    arg_k_1  = mxCreateDoubleScalar(k_1) ;
    arg_dk_1 = mxCreateDoubleScalar(dk_1) ;
    arg_L_1  = mxCreateDoubleScalar(L_1) ;

    arg_k_2  = mxCreateDoubleScalar(k_2) ;
    arg_dk_2 = mxCreateDoubleScalar(dk_2) ;
    arg_L_2  = mxCreateDoubleScalar(L_2) ;
  
  } else {
    Clothoid::valueType k, dk, L ;
    iter = Clothoid::buildClothoid( mxGetScalar(arg_x0),
                                    mxGetScalar(arg_y0),
                                    mxGetScalar(arg_theta0),
                                    mxGetScalar(arg_x1),
                                    mxGetScalar(arg_y1),
                                    mxGetScalar(arg_theta1),
                                    k, dk, L ) ;
    arg_k  = mxCreateDoubleScalar(k) ;
    arg_dk = mxCreateDoubleScalar(dk) ;
    arg_L  = mxCreateDoubleScalar(L) ;
  }
  if ( nlhs >= 4 ) {
    arg_iter = mxCreateNumericMatrix(1,1,mxINT32_CLASS,mxREAL) ;
    *((int*)mxGetData(arg_iter)) = iter ;
  }

}
