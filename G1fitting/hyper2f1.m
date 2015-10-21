%=============================================================================%
%  hyper2f1:  Compute Gauss Hypergeometric function                           %
%                                                                             %
%  USAGE: res = hyper2f1( a, b, c, z ) ;                                      %
%                                                                             %
%=============================================================================%
%                                                                             %
%  Autors: Enrico Bertolazzi                                                  %
%          Department of Industrial Engineering                               %
%          University of Trento                                               %
%          enrico.bertolazzi@unitn.it                                         %
%                                                                             %
%=============================================================================%
function res = hyper2f1(a,b,c,z)
  %-------------------------------------------------------------------------%
  % based on the work of John Pearson, University of Oxford                 %
  % in the MSc dissertation 'Computation of Hypergeometric Functions'       %
  %-------------------------------------------------------------------------%
  % Applies transformation formulae detailed in Section 4.6.                %
  %-------------------------------------------------------------------------%
  if z < -1
    mz  = 1-z ;
    z1  = 1/mz ;
    bma = b-a ;
    cma = c-a ;
    cmb = c-b ;
    cf1 = gamma(bma)/gamma(b)/gamma(cma)*eval2f1(a,cmb,1-bma,z1) ;
    cf2 = gamma(-bma)/gamma(a)/gamma(cmb)*eval2f1(b,cma,1+bma,z1) ;
    res = gamma(c) * ( mz^(-a)*cf1+mz^(-b)*cf2 ) ;
  elseif z < 0
    res = (1-z)^(-a)*eval2f1(a,c-b,c,z/(z-1)) ;
  elseif z <= 0.5
    res = eval2f1(a,b,c,z) ;
  elseif z < 1
    z1    = 1-z ;
    cma   = c-a ;
    cmb   = c-b ;
    cmamb = c-a-b ;
    cf1   = gamma(cmamb)/gamma(cma)/gamma(cmb)*eval2f1(a,b,1-cmamb,z1) ;
    cf2   = gamma(-cmamb)/gamma(a)/gamma(b)*eval2f1(cma,cmb,cmamb+1,z1) ;
    res   = gamma(c)*( cf1 +z1^cmamb*cf2 ) ;
  elseif z < 2
    z1    = 1-1/z ;
    cma   = c-a ;
    cmamb = cma-b ;
    cf1   = gamma(cmamb)/gamma(cma)/gamma(c-b)*eval2f1(a,1-cma,1-cmamb,z1) ;
    cf2   = gamma(-cmamb)/gamma(a)/gamma(b)*eval2f1(cma,1-a,1+cmamb,z1) ;
    res   = gamma(c) * ( z^(-a)*cf1 +z^(-cma)*(1-z)^cmamb*cf2 ) ;
  else
    z1  = 1/z ;
    cma = c-a ;
    cmb = c-b ;
    bma = b-a ;
    cf1 = gamma(bma)/gamma(b)/gamma(cma)*eval2f1(a,1-cma,-bma,z1) ;
    cf2 = gamma(-bma)/gamma(a)/gamma(cmb)*eval2f1(1-cmb,b,bma+1,z1) ;
    res = gamma(c)*( (-z)^(-a)*cf1 +(-z)^(-b)*cf2 ) ;
  end
end
%
%
%
function res = eval2f1(a,b,c,z)
  tol = 1e-20 ;
  % Initialise vector of individual terms and sum of terms computed thus far
  a1 = 1 ;
  b1 = 1 ;
  aj = a ;
  bj = b ;
  cj = c ;
  for j=1:500
    % Update value of a1, current term, and b1, sum of all terms in terms of previous values
    a1new = (aj*bj/(cj*j))*z*a1 ;
    b1    = b1 + a1new ;
    % Terminate summation if stopping criterion is satisfied
    tol1 = abs(b1)*tol ;
    if abs(a1) < tol1 && abs(a1new) < tol1
      break ;
    end;
    aj = aj+1 ;
    bj = bj+1 ;
    cj = cj+1 ;
    a1 = a1new ;
  end
  if j >= 500
    error( ['eval2f1(' a ',' b ',' c ') failed to converge'] );
  end
  res = b1 ;
end
