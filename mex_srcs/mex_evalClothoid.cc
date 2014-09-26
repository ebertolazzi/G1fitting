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

#define arg_x0     prhs[0]
#define arg_y0     prhs[1]
#define arg_theta0 prhs[2]
#define arg_k      prhs[3]
#define arg_dk     prhs[4]
#define arg_s      prhs[5]

#define arg_X      plhs[0]
#define arg_Y      plhs[1]
#define arg_TH     plhs[2]
#define arg_CURV   plhs[3]

#define ASSERT(COND,MSG)                         \
  if ( !(COND) ) {                               \
    std::ostringstream ost ;                     \
    ost << "pointsOnClothoid: " << MSG << '\n' ; \
    mexErrMsgTxt(ost.str().c_str()) ;            \
  }

extern "C"
void
mexFunction( int nlhs, mxArray       *plhs[],
             int nrhs, mxArray const *prhs[] ) {

  // Check for proper number of arguments, etc
  if ( nrhs != 6 ) {
	  mexErrMsgTxt(
"%======================================================================%\n"
"% evalClothoid:  Compute clothoid parameters along a Clothoid curve    %\n"
"%                                                                      %\n"
"% USAGE: [X,Y,TH,K] = evalClothoid( x0, y0, theta0, k, dk, s ) ;       %\n"
"%                                                                      %\n"
"% On input:                                                            %\n"
"%                                                                      %\n"
"%  x0, y0 = coodinate of initial point                                 %\n"
"%  theta0 = orientation (angle) of the clothoid at initial point       %\n"
"%  k      = curvature at initial point                                 %\n"
"%  dk     = derivative of curvature respect to arclength               %\n"
"%  s      = vector of curvilinear coordinate where to compute clothoid %\n"
"%                                                                      %\n"
"% On output:                                                           %\n"
"%                                                                      %\n"
"%  X     = X coordinate of the points of the clothoid at s coordinate  %\n"
"%  Y     = Y coordinate of the points of the clothoid at s coordinate  %\n"
"%  TH    = angle of the clothoid at s coordinate                       %\n"
"%  K     = curvature of the clothoid at s coordinate                   %\n"
"%======================================================================%\n"
"%                                                                      %\n"
"%  Autor: Enrico Bertolazzi                                            %\n"
"%         Department of Industrial Engineering                         %\n"
"%         University of Trento                                         %\n"
"%         enrico.bertolazzi@unitn.it                                   %\n"
"%                                                                      %\n"
"%======================================================================%\n" ) ;
  }

  for ( int kk = 0 ; kk < 5 ; ++kk )
    if ( mxGetM(prhs[kk]) != 1 || mxGetN(prhs[kk]) != 1 )
	    mexErrMsgTxt("First 5 input arguments must be scalars");

  for ( int kk = 0 ; kk < 6 ; ++kk )
    if ( mxGetClassID(prhs[kk]) != mxDOUBLE_CLASS || mxIsComplex(prhs[kk]) )
	    mexErrMsgTxt("Input argument should be real double (not complex)");

  Clothoid::valueType x0     = mxGetScalar(arg_x0) ;
  Clothoid::valueType y0     = mxGetScalar(arg_y0) ;
  Clothoid::valueType theta0 = mxGetScalar(arg_theta0) ;
  Clothoid::valueType k      = mxGetScalar(arg_k) ;
  Clothoid::valueType dk     = mxGetScalar(arg_dk) ;
  
  // Output array
  mwSize         nDimNum = mxGetNumberOfDimensions(arg_s);
  mwSize const * pDims   = mxGetDimensions(arg_s) ;

  double *pX, *pY, *pTH, *pCURV ;

  if ( nlhs > 0 ) {
    arg_X = mxCreateNumericArray(nDimNum, pDims, mxDOUBLE_CLASS, mxREAL);
    pX    = mxGetPr(arg_X) ;
  }
  if ( nlhs > 1 ) {
    arg_Y = mxCreateNumericArray(nDimNum, pDims, mxDOUBLE_CLASS, mxREAL);
    pY    = mxGetPr(arg_Y) ;
  }
  if ( nlhs > 2 ) {
    arg_TH = mxCreateNumericArray(nDimNum, pDims, mxDOUBLE_CLASS, mxREAL);
    pTH    = mxGetPr(arg_TH) ;
  }
  if ( nlhs > 3 ) {
    arg_CURV = mxCreateNumericArray(nDimNum, pDims, mxDOUBLE_CLASS, mxREAL);
    pCURV    = mxGetPr(arg_CURV) ;
  }

  int nElemNum = mxGetNumberOfElements(arg_s);
  double * pS = mxGetPr(arg_s) ;
  for ( int i = 0 ; i < nElemNum ; ++i ) {
    Clothoid::valueType C, S, t = *pS++ ;
    Clothoid::GeneralizedFresnelCS( dk*t*t, k*t, theta0, C, S ) ;
    if ( nlhs > 0 ) *pX++    = x0 + t*C ;
    if ( nlhs > 1 ) *pY++    = y0 + t*S ;
    if ( nlhs > 3 ) *pTH++   = theta0 + t*(k+t*(dk/2)) ;
    if ( nlhs > 4 ) *pCURV++ = k+t*dk ;
  }

}
