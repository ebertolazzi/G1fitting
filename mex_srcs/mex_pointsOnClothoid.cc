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
#define arg_L      prhs[5]
#define arg_npts   prhs[6]

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
  if ( nrhs != 7 ) {
	  mexErrMsgTxt(
"%======================================================================%\n"
"% pointsOnClothoid:  Compute points on a clothoid curve.               %\n"
"%                    Used for plotting purpose.                        %\n"
"%                                                                      %\n"
"% USAGE: XY    = pointsOnClothoid( x0, y0, theta0, k, dk, L, npts ) ;  %\n"
"% USAGE: [X,Y] = pointsOnClothoid( x0, y0, theta0, k, dk, L, npts ) ;  %\n"
"%                                                                      %\n"
"% On input:                                                            %\n"
"%                                                                      %\n"
"%  x0, y0  = coodinate of initial point                                %\n"
"%  theta0  = orientation (angle) of the clothoid at initial point      %\n"
"%  k       = curvature at initial point                                %\n"
"%  dk      = derivative of curvature respect to arclength              %\n"
"%  L       = the lenght of the clothoid curve                          %\n"
"%  npts    = number of points along the clothoid                       %\n"
"%                                                                      %\n"
"% On output: (1 argument)                                              %\n"
"%                                                                      %\n"
"%  XY = matrix 2 x NPTS whose column are the points of the clothoid    %\n"
"%                                                                      %\n"
"% On output: (2 argument)                                              %\n"
"%                                                                      %\n"
"%   X  = matrix 1 x NPTS X coordinate of points of the clothoid        %\n"
"%   Y  = matrix 1 x NPTS Y coordinate of points of the clothoid        %\n"
"%                                                                      %\n"
"%======================================================================%\n"
"%                                                                      %\n"
"%  Autor: Enrico Bertolazzi                                            %\n"
"%         Department of Industrial Engineering                         %\n"
"%         University of Trento                                         %\n"
"%         enrico.bertolazzi@unitn.it                                   %\n"
"%                                                                      %\n"
"%======================================================================%\n" ) ;
  }

  for ( int kk = 0 ; kk < 7 ; ++kk )
    if ( mxGetM(prhs[kk]) != 1 || mxGetN(prhs[kk]) != 1 )
	    mexErrMsgTxt("Input arguments must be scalars");

  for ( int kk = 0 ; kk < 6 ; ++kk )
    if ( mxGetClassID(prhs[kk]) != mxDOUBLE_CLASS || mxIsComplex(prhs[kk]) )
	    mexErrMsgTxt("First 6 argument must be real double scalars (not complex)");

  Clothoid::valueType x0     = mxGetScalar(arg_x0) ;
  Clothoid::valueType y0     = mxGetScalar(arg_y0) ;
  Clothoid::valueType theta0 = mxGetScalar(arg_theta0) ;
  Clothoid::valueType k      = mxGetScalar(arg_k) ;
  Clothoid::valueType dk     = mxGetScalar(arg_dk) ;
  Clothoid::valueType L      = mxGetScalar(arg_L) ;
  int                 npts   = int(mxGetScalar(arg_npts)) ;

  ASSERT( npts > 1, "7th arguments (npts) must be > 1, found " << npts ) ;
  ASSERT( L    > 0, "6th arguments (length) must be > 0, found " << L ) ;

  Clothoid::valueType C, S, dt = L/npts ;

  if ( nlhs < 2 ) {
	  plhs[0] = mxCreateDoubleMatrix(2, npts, mxREAL);
	  double * pXY = mxGetPr(plhs[0]);
    for ( Clothoid::valueType t = 0 ; t <= L ; t += dt ) {
      Clothoid::GeneralizedFresnelCS( dk*t*t, k*t, theta0, C, S ) ;
      *pXY++ = x0 + t*C ;
      *pXY++ = y0 + t*S ;
    }
  } else {
	  plhs[0] = mxCreateDoubleMatrix(1, npts, mxREAL);
	  plhs[1] = mxCreateDoubleMatrix(1, npts, mxREAL);
	  double * pX = mxGetPr(plhs[0]);
	  double * pY = mxGetPr(plhs[1]);
    for ( Clothoid::valueType t = 0 ; t <= L ; t += dt ) {
      Clothoid::GeneralizedFresnelCS( dk*t*t, k*t, theta0, C, S ) ;
      *pX++ = x0 + t*C ;
      *pY++ = y0 + t*S ;
    }
  }

}
