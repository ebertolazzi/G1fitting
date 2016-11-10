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

#define ASSERT(COND,MSG)                         \
  if ( !(COND) ) {                               \
    std::ostringstream ost ;                     \
    ost << "pointsOnClothoid: " << MSG << '\n' ; \
    mexErrMsgTxt(ost.str().c_str()) ;            \
  }

#define arg_x0     prhs[0]
#define arg_y0     prhs[1]
#define arg_theta0 prhs[2]
#define arg_k      prhs[3]
#define arg_dk     prhs[4]
#define arg_L      prhs[5]
#define arg_npts   prhs[6]

#define MEX_ERROR_MESSAGE \
"%======================================================================%\n" \
"% pointsOnClothoid:  Compute points on a clothoid curve.               %\n" \
"%                    Used for plotting purpose.                        %\n" \
"%                                                                      %\n" \
"% USAGE: XY    = pointsOnClothoid( x0, y0, theta0, k, dk, L, npts ) ;  %\n" \
"% USAGE: [X,Y] = pointsOnClothoid( x0, y0, theta0, k, dk, L, npts ) ;  %\n" \
"% USAGE: XY    = pointsOnClothoid( x0, y0, theta0, k, dk, L ) ;        %\n" \
"% USAGE: [X,Y] = pointsOnClothoid( x0, y0, theta0, k, dk, L ) ;        %\n" \
"%                                                                      %\n" \
"% On input:                                                            %\n" \
"%                                                                      %\n" \
"%  x0, y0  = coodinate of initial point                                %\n" \
"%  theta0  = orientation (angle) of the clothoid at initial point      %\n" \
"%  k       = curvature at initial point                                %\n" \
"%  dk      = derivative of curvature respect to arclength              %\n" \
"%  L       = the lenght of the clothoid curve or a vector of length    %\n" \
"%            where to compute the clothoid values                      %\n" \
"%  npts    = number of points along the clothoid                       %\n" \
"%                                                                      %\n" \
"% On output: (1 argument)                                              %\n" \
"%                                                                      %\n" \
"%  XY = matrix 2 x NPTS whose column are the points of the clothoid    %\n" \
"%                                                                      %\n" \
"% On output: (2 argument)                                              %\n" \
"%                                                                      %\n" \
"%   X  = matrix 1 x NPTS X coordinate of points of the clothoid        %\n" \
"%   Y  = matrix 1 x NPTS Y coordinate of points of the clothoid        %\n" \
"%                                                                      %\n" \
"%======================================================================%\n" \
"%                                                                      %\n" \
"%  Autor: Enrico Bertolazzi                                            %\n" \
"%         Department of Industrial Engineering                         %\n" \
"%         University of Trento                                         %\n" \
"%         enrico.bertolazzi@unitn.it                                   %\n" \
"%                                                                      %\n" \
"%======================================================================%\n"

extern "C"
void
mexFunction( int nlhs, mxArray       *plhs[],
             int nrhs, mxArray const *prhs[] ) {

  // Check for proper number of arguments, etc
  if ( nrhs < 6 || nrhs > 7 ) {
    mexErrMsgTxt(MEX_ERROR_MESSAGE) ;
    return ;
  }

  for ( int kk = 0 ; kk < 5 ; ++kk )
    if ( mxGetM(prhs[kk]) != 1 || mxGetN(prhs[kk]) != 1 ) {
      mexErrMsgTxt("First 5 input arguments must be scalars");
      return ;
    }

  for ( int kk = 0 ; kk < 6 ; ++kk )
    if ( mxGetClassID(prhs[kk]) != mxDOUBLE_CLASS || mxIsComplex(prhs[kk]) ) {
      mexErrMsgTxt("First 6 argument must be real double scalars (not complex)");
      return ;
    }

  Clothoid::valueType x0     = mxGetScalar(arg_x0) ;
  Clothoid::valueType y0     = mxGetScalar(arg_y0) ;
  Clothoid::valueType theta0 = mxGetScalar(arg_theta0) ;
  Clothoid::valueType k      = mxGetScalar(arg_k) ;
  Clothoid::valueType dk     = mxGetScalar(arg_dk) ;

  Clothoid::valueType C, S ;

  if ( nrhs == 7 ) {

    Clothoid::valueType L    = mxGetScalar(arg_L) ;
    int                 npts = int(mxGetScalar(arg_npts)) ;

    //ASSERT( L    > 0, "6th arguments (L) must be > 0, found " << L ) ;
    ASSERT( npts > 1, "7th arguments (npts) must be > 1, found " << npts ) ;

    Clothoid::valueType dt = L/(npts-1) ;

    if ( nlhs == 1 ) {
  	  plhs[0] = mxCreateDoubleMatrix(2, npts, mxREAL);
  	  double * pXY = mxGetPr(plhs[0]);
      *pXY++ = x0 ;
      *pXY++ = y0 ;
      for ( int i = 1 ; i < npts-1 ; ++i ) {
        Clothoid::valueType t = i*dt ;
        Clothoid::GeneralizedFresnelCS( dk*(t*t), k*t, theta0, C, S ) ;
        *pXY++ = x0 + t*C ;
        *pXY++ = y0 + t*S ;
      }
      Clothoid::GeneralizedFresnelCS( dk*(L*L), k*L, theta0, C, S ) ;
      *pXY++ = x0 + L*C ;
      *pXY++ = y0 + L*S ;
    } else if ( nlhs == 2 ) {
      plhs[0] = mxCreateDoubleMatrix(1, npts, mxREAL);
      plhs[1] = mxCreateDoubleMatrix(1, npts, mxREAL);
      double * pX = mxGetPr(plhs[0]);
      double * pY = mxGetPr(plhs[1]);
      *pX++ = x0 ;
      *pY++ = y0 ;
      for ( int i = 1 ; i < npts-1 ; ++i ) {
        Clothoid::valueType t = i*dt ;
        Clothoid::GeneralizedFresnelCS( dk*(t*t), k*t, theta0, C, S ) ;
        *pX++ = x0 + t*C ;
        *pY++ = y0 + t*S ;
      }
      Clothoid::GeneralizedFresnelCS( dk*(L*L), k*L, theta0, C, S ) ;
      *pX++ = x0 + L*C ;
      *pY++ = y0 + L*S ;
    } else {
      mexErrMsgTxt("Output argument must be 1 or 2");
    }

  } else {

    double * L = mxGetPr(arg_L) ;
    int   npts = mxGetN(arg_L)*mxGetM(arg_L) ;

    if ( nlhs == 1 ) {
  	  plhs[0] = mxCreateDoubleMatrix(2, npts, mxREAL);
  	  double * pXY = mxGetPr(plhs[0]);
      for ( int i = 0 ; i < npts ; ++i ) {
        double t = L[i] ;
        Clothoid::GeneralizedFresnelCS( dk*(t*t), k*t, theta0, C, S ) ;
        *pXY++ = x0 + t*C ;
        *pXY++ = y0 + t*S ;
      }
    } else if ( nlhs == 2 ) {
      plhs[0] = mxCreateDoubleMatrix(1, npts, mxREAL);
      plhs[1] = mxCreateDoubleMatrix(1, npts, mxREAL);
      double * pX = mxGetPr(plhs[0]);
      double * pY = mxGetPr(plhs[1]);
      for ( int i = 0 ; i < npts ; ++i ) {
        double t = L[i] ;
        Clothoid::GeneralizedFresnelCS( dk*(t*t), k*t, theta0, C, S ) ;
        *pX++ = x0 + t*C ;
        *pY++ = y0 + t*S ;
      }
    } else {
      mexErrMsgTxt("Output argument must be 1 or 2");
    }
  }
}
