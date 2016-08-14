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
"%======================================================================%\n" \
"% pointsOnClothoid:  Compute points on a clothoid curve.               %\n" \
"%                    Used for plotting purpose.                        %\n" \
"%                                                                      %\n" \
"% USAGE: XY    = pointsOnClothoid( x0, y0, theta0, k, dk, L, npts ) ;  %\n" \
"% USAGE: [X,Y] = pointsOnClothoid( x0, y0, theta0, k, dk, L, npts ) ;  %\n" \
"% USAGE: XY    = pointsOnClothoid( x0, y0, theta0, k, dk, L ) ;        %\n" \
"% USAGE: [X,Y] = pointsOnClothoid( x0, y0, theta0, k, dk, L ) ;        %\n" \
"% USAGE: XY    = pointsOnClothoid( clot, npts ) ;                      %\n" \
"% USAGE: [X,Y] = pointsOnClothoid( clot, npts ) ;                      %\n" \
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
"% In alternative                                                       %\n" \
"%  clot    = structure with the field x0, y0, kappa, dkappa, L         %\n" \
"%                                                                      %\n" \
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

#define arg_x0     prhs[0]
#define arg_y0     prhs[1]
#define arg_theta0 prhs[2]
#define arg_k      prhs[3]
#define arg_dk     prhs[4]
#define arg_L      prhs[5]
#define arg_npts   prhs[6]
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
  Clothoid::valueType x0, y0, theta0, k, dk, L ;
  Clothoid::indexType npts ;

  if ( nrhs == 2 ) {
    if ( mxGetClassID(prhs[0]) == mxSTRUCT_CLASS ) {
      mxArray * mxX0     = mxGetField( prhs[0], 0, "x0"     );
      mxArray * mxY0     = mxGetField( prhs[0], 0, "y0"     );
      mxArray * mxTheta0 = mxGetField( prhs[0], 0, "theta0" );
      mxArray * mxK      = mxGetField( prhs[0], 0, "kappa"  );
      mxArray * mxDK     = mxGetField( prhs[0], 0, "dkappa" );
      mxArray * mxL      = mxGetField( prhs[0], 0, "L"      );
      if ( mxX0 == nullptr ) {
    	  mexErrMsgTxt("missing field `x0` in first arguments\n") ;
        mexErrMsgTxt(MEX_ERROR_MESSAGE) ;
        return ;
      } else {
        x0 = mxGetScalar(mxX0) ;        
      }
      if ( mxY0 == nullptr ) {
    	  mexErrMsgTxt("missing field `y0` in first arguments\n") ;
        mexErrMsgTxt(MEX_ERROR_MESSAGE) ;
        return ;
      } else {
        y0 = mxGetScalar(mxY0) ;        
      }
      if ( mxTheta0 == nullptr ) {
    	  mexErrMsgTxt("missing field `theta0` in first arguments\n") ;
        mexErrMsgTxt(MEX_ERROR_MESSAGE) ;
        return ;        
      } else {
        theta0 = mxGetScalar(mxTheta0) ;        
      }
      if ( mxK == nullptr ) {
    	  mexErrMsgTxt("missing field `kappa` in first arguments\n") ;
        mexErrMsgTxt(MEX_ERROR_MESSAGE) ;
        return ;
      } else {
        k = mxGetScalar(mxK) ;                
      }
      if ( mxDK == nullptr ) {
    	  mexErrMsgTxt("missing field `dkappa` in first arguments\n") ;
        mexErrMsgTxt(MEX_ERROR_MESSAGE) ;
        return ;
      } else {
        dk = mxGetScalar(mxDK) ;                        
      }
      if ( mxL == nullptr ) {
    	  mexErrMsgTxt("missing field `L` in first arguments\n") ;
        mexErrMsgTxt(MEX_ERROR_MESSAGE) ;
        return ;
      } else {
        L = mxGetScalar(mxL) ;                                
      }
      npts = mxGetScalar(prhs[1]) ;  
    } else {
  	  mexErrMsgTxt("first argument expected to be a STRUCT\n") ;
      mexErrMsgTxt(MEX_ERROR_MESSAGE) ;
      return ;
    }
  } else if ( nrhs == 6 ) {
    for ( int kk = 0 ; kk < 5 ; ++kk )
      if ( mxGetM(prhs[kk]) != 1 || mxGetN(prhs[kk]) != 1 )
  	    mexErrMsgTxt("First 7 input arguments must be scalars");
    for ( int kk = 0 ; kk < 6 ; ++kk )
      if ( mxGetClassID(prhs[kk]) != mxDOUBLE_CLASS || mxIsComplex(prhs[kk]) )
  	    mexErrMsgTxt("Input argument should be real double (not complex)");
    x0     = mxGetScalar(arg_x0) ;
    y0     = mxGetScalar(arg_y0) ;
    theta0 = mxGetScalar(arg_theta0) ;
    k      = mxGetScalar(arg_k) ;
    dk     = mxGetScalar(arg_dk) ;
  } else if ( nrhs == 7 ) {
    for ( int kk = 0 ; kk < 7 ; ++kk )
      if ( mxGetM(prhs[kk]) != 1 || mxGetN(prhs[kk]) != 1 )
  	    mexErrMsgTxt("First 7 input arguments must be scalars");
    for ( int kk = 0 ; kk < 6 ; ++kk )
      if ( mxGetClassID(prhs[kk]) != mxDOUBLE_CLASS || mxIsComplex(prhs[kk]) )
  	    mexErrMsgTxt("Input argument should be real double (not complex)");
    x0     = mxGetScalar(arg_x0) ;
    y0     = mxGetScalar(arg_y0) ;
    theta0 = mxGetScalar(arg_theta0) ;
    k      = mxGetScalar(arg_k) ;
    dk     = mxGetScalar(arg_dk) ;
    L      = mxGetScalar(arg_L) ;
    npts   = mxGetScalar(arg_npts) ;
  } else {
	  mexErrMsgTxt("expexted 2,6 or 7 input arguments\n") ;
    mexErrMsgTxt(MEX_ERROR_MESSAGE) ;
    return ;
  }
  
  // Output array

  double *pX, *pY, *pTH, *pCURV ;
  if ( nrhs == 6 ) {
    mwSize         nDimNum = mxGetNumberOfDimensions(arg_s);
    mwSize const * pDims   = mxGetDimensions(arg_s) ;

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
  } else {
    if ( nlhs > 0 ) {
      arg_X = mxCreateNumericMatrix(npts,1, mxDOUBLE_CLASS, mxREAL);
      pX    = mxGetPr(arg_X) ;
    }
    if ( nlhs > 1 ) {
      arg_Y = mxCreateNumericMatrix(npts,1, mxDOUBLE_CLASS, mxREAL);
      pY    = mxGetPr(arg_Y) ;
    }
    if ( nlhs > 2 ) {
      arg_TH = mxCreateNumericMatrix(npts,1, mxDOUBLE_CLASS, mxREAL);
      pTH    = mxGetPr(arg_TH) ;
    }
    if ( nlhs > 3 ) {
      arg_CURV = mxCreateNumericMatrix(npts,1, mxDOUBLE_CLASS, mxREAL);
      pCURV    = mxGetPr(arg_CURV) ;
    }
    for ( int i = 0 ; i < npts ; ++i ) {
      Clothoid::valueType t = (i*L)/(npts-1) ;
      Clothoid::valueType C, S ;
      Clothoid::GeneralizedFresnelCS( dk*t*t, k*t, theta0, C, S ) ;
      if ( nlhs > 0 ) *pX++    = x0 + t*C ;
      if ( nlhs > 1 ) *pY++    = y0 + t*S ;
      if ( nlhs > 3 ) *pTH++   = theta0 + t*(k+t*(dk/2)) ;
      if ( nlhs > 4 ) *pCURV++ = k+t*dk ;
    }
  }

}
