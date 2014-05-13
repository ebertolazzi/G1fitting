%=============================================================================%
%  buildClothoid:  Compute parameters of the G1 Hermite clothoid fitting      %
%                                                                             %
%  USAGE: [k,dk,L] = buildClothoid( x0, y0, theta0, x1, y1, theta1, tol ) ;   %
%                                                                             %
%  On input:                                                                  %
%                                                                             %
%       x0, y0  = coodinate of initial point                                  %
%       theta0  = orientation (angle) of the clothoid at initial point        %
%       x1, y1  = coodinate of final point                                    %
%       theta1  = orientation (angle) of the clothoid at final point          %
%       tol     = tolerance used to stop Newton Iterations                    %
%                                                                             %
%  On output:                                                                 %
%                                                                             %
%       L  = the lenght of the clothoid curve from initial to final point     %
%       k  = curvature at initial point                                       %
%       dk = derivative of curvature respect to arclength,                    %
%            notice that curvature at final point is k+dk*L                   %
%       iter = Newton Iterations used to solve the interpolation problem      %
%                                                                             %
%=============================================================================%
%                                                                             %
%  Autors: Enrico Bertolazzi and Marco Frego                                  %
%          Department of Industrial Engineering                               %
%          University of Trento                                               %
%          enrico.bertolazzi@unitn.it                                         %
%          m.fregox@gmail.com                                                 %
%                                                                             %
%=============================================================================%
function [ k, dk, L, iter ] = buildClothoidMW( x0, y0, theta0, x1, y1, theta1, tol )
  global MW_fails MW_calls MW_calls_deg ;
  [P0,T0,N0,a,t1,t2,iter,failFlag,reflectFlag,reverseFlag,kind] = completeShape( [x0;y0], [cos(theta0),sin(theta0)], [x1;y1], [cos(theta1),sin(theta1)], tol, 100 ) ;
  MW_calls = MW_calls+1 ;
  if failFlag
    MW_fails = MW_fails+1 ;
    %warning('non converge') ;
  end
  
  if kind == 0
    dk = pi/a^2 ;
    L = a*(t2-t1) ;
    if reverseFlag
      k = -pi*t2/a ;
    else
      k = pi*t1/a ;
    end
    if reflectFlag
      k  = -k ;
      dk = -dk ;
    end
  elseif kind == 1 % circle
    MW_calls_deg = MW_calls_deg+1 ;
    dk = 0 ;
    d  = sqrt((x1-x0)^2+(y1-y0)^2) ;
    k  = 2*sin(a)/d ;
    L  = a*d/sin(a) ;
    if reverseFlag
      k = -k ;
    end
    if reflectFlag
      k = -k ;
    end
  else % straight line
    MW_calls_deg = MW_calls_deg+1 ;
    dk = 0 ;
    k  = 0 ;
    L  = sqrt((x1-x0)^2+(y1-y0)^2) ;
  end
end