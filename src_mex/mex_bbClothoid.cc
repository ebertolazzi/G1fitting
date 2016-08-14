/****************************************************************************\
  Copyright (c) Enrico Bertolazzi 2016
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

#include <vector>
#include <sstream>
#include <stdexcept>

#define ASSERT(COND,MSG)                   \
  if ( !(COND) ) {                         \
    std::ostringstream ost ;               \
    ost << "bbClothoid: " << MSG << '\n' ; \
    mexErrMsgTxt(ost.str().c_str()) ;      \
  }

#define arg_x0     prhs[0]
#define arg_y0     prhs[1]
#define arg_theta0 prhs[2]
#define arg_k      prhs[3]
#define arg_dk     prhs[4]
#define arg_L      prhs[5]
#define arg_angle  prhs[6]
#define arg_size   prhs[7]
#define arg_offs   prhs[8]

#define MEX_ERROR_MESSAGE \
"%======================================================================%\n" \
"%  bbClothoid:  Compute a series of bounding triangles for             %\n" \
"%               a clothoid curve.                                      %\n" \
"%                                                                      %\n" \
"%  USAGE:                                                              %\n" \
"%    TT = bbClothoid( x0, y0, theta0, k, dk, L, angle, size ) ;        %\n" \
"%    TT = bbClothoid( x0, y0, theta0, k, dk, L, angle, size, offs ) ;  %\n" \
"%                                                                      %\n" \
"%  On input:                                                           %\n" \
"%                                                                      %\n" \
"%  x0, y0 = coodinate of initial point                                 %\n" \
"%  theta0 = orientation (angle) of the clothoid at initial point       %\n" \
"%  k      = curvature at initial point                                 %\n" \
"%  dk     = derivative of curvature respect to arclength               %\n" \
"%  L      = the lenght of the clothoid curve or a vector of length     %\n" \
"%           where to compute the clothoid values                       %\n" \
"%  angle  = maximum variation of angle in the bounding box triangle    %\n" \
"%  size   = maximum height of the bounding box triangle                %\n" \
"%                                                                      %\n" \
"%  On output:                                                          %\n" \
"%                                                                      %\n" \
"%  TT = matrix 6 x N whose column are the coordinates of bounding      %\n" \
"%       triangles.                                                     %\n" \
"%       TT(:,i) = [ x0, y0, x1, y1, x2, y2 ].'                         %\n" \
"%                                                                      %\n" \
"%======================================================================%\n" \
"%                                                                      %\n" \
"%  Autor: Enrico Bertolazzi                                            %\n" \
"%         Department of Industrial Engineering                         %\n" \
"%         University of Trento                                         %\n" \
"%         enrico.bertolazzi@unitn.it                                   %\n" \
"%                                                                      %\n" \
"%======================================================================%\n"

#include <iostream> 

extern "C"
void
mexFunction( int nlhs, mxArray       *plhs[],
             int nrhs, mxArray const *prhs[] ) {

  // Check for proper number of arguments, etc
  if ( ! (nrhs == 8 || nrhs == 9) ) {
    mexErrMsgTxt(MEX_ERROR_MESSAGE) ;
    return ;
  }

  for ( int kk = 0 ; kk < nrhs ; ++kk )
    if ( mxGetM(prhs[kk]) != 1 || mxGetN(prhs[kk]) != 1 ) {
      mexErrMsgTxt("Input arguments must be scalars");
      return ;
    }

  Clothoid::valueType x0        = mxGetScalar(arg_x0) ;
  Clothoid::valueType y0        = mxGetScalar(arg_y0) ;
  Clothoid::valueType theta0    = mxGetScalar(arg_theta0) ;
  Clothoid::valueType k         = mxGetScalar(arg_k) ;
  Clothoid::valueType dk        = mxGetScalar(arg_dk) ;
  Clothoid::valueType L         = mxGetScalar(arg_L) ;
  Clothoid::valueType max_angle = mxGetScalar(arg_angle) ;
  Clothoid::valueType max_size  = mxGetScalar(arg_size) ;

  Clothoid::valueType offs = 0 ;
  if ( nrhs == 9 ) offs = mxGetScalar(arg_offs) ;

  // costruisco bb
  Clothoid::ClothoidCurve clot( x0, y0, theta0, k, dk, L ) ;
  std::vector<Clothoid::ClothoidCurve> c ;
  std::vector<Clothoid::Triangle2D>    t ;

  clot.bbSplit( max_angle, max_size, offs, c, t ) ;

  if ( nlhs == 1 ) {
    plhs[0] = mxCreateDoubleMatrix(6, t.size(), mxREAL);
    double * pT = mxGetPr(plhs[0]);
    for ( int i = 0 ; i < t.size() ; ++i ) {
      *pT++ = t[i].x1() ; *pT++ = t[i].y1() ;
      *pT++ = t[i].x2() ; *pT++ = t[i].y2() ;
      *pT++ = t[i].x3() ; *pT++ = t[i].y3() ;
    }
  } else {
    mexErrMsgTxt("There must be an output argument\n");
  }
}
