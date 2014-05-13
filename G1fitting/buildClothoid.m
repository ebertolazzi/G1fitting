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
function [ k, dk, L, iter ] = buildClothoid( x0, y0, theta0, x1, y1, theta1, tol )

  dx  = x1 - x0 ;
  dy  = y1 - y0 ;
  r   = sqrt( dx^2 + dy^2 ) ;
  phi = atan2( dy, dx ) ;

  phi0 = normalizeAngle(theta0 - phi) ;
  phi1 = normalizeAngle(theta1 - phi) ;

  % % check if we must reverse the problem
  % reverse = false ;
  % if abs(phi1) < abs(phi0)
  %   reverse = true ;
  %   tmp  = phi1 ;
  %   phi1 = phi0 ;
  %   phi0 = tmp  ;
  % end
  % 
  % % check if we must mirror the problem
  % mirror = false ;
  % if phi1 < 0
  %   mirror = true ;
  %   phi1 = -phi1 ;
  %   phi0 = -phi0 ;
  % end

  delta = phi1 - phi0 ;

  % initial point
  Aguess = guessA( phi0, phi1 ) ;

  % Newton iteration
  [A,iter] = findA( Aguess, delta, phi0, tol ) ;

  % final operation
  [h,g] = GeneralizedFresnelCS( 1, 2*A, delta-A, phi0 ) ;
  L = r/h ;
  
  if L > 0
    k  = (delta - A)/L ;
    dk = 2*A/L^2 ;
  else
    error('negative length') ;
  end
  
  % % check if we must mirror the problem
  % if mirror
  %   k  = -k  ;
  %   dk = -dk ;
  % end
  % 
  % % reverse the problem
  % if reverse
  %   k  = k+L*dk ;
  %   dk = -dk ;
  % end

end
%=============================================================================%
%  normalizeAngle:  normalize angle in the range [-pi,pi]                     %
%=============================================================================%
function phi = normalizeAngle( phi_in )
  phi = phi_in ;
  while ( phi > pi )
    phi = phi - 2*pi ;
  end
  while ( phi < -pi )
    phi = phi + 2*pi ;
  end
end
%=============================================================================%
%  findA:  Find a zero of function g(A) defined as                            %
%  g(A) = \int_0^1 \sin( A*t^2+(delta-A)*t+phi0 ) dt                          %
%                                                                             %
%  USAGE:  A = findA( Aguess, delta, phi0, tol );                             %
%                                                                             %
%  Given an initial guess Aguess find the closest zero of equation g(A)       %
%                                                                             %
%  On input:                                                                  %
%       Aguess      = initial guess.                                          %
%       delta, phi0 = Angles used in the clothoid fitting problem.            %
%       tol         = Tolerance for stopping criterium of Newton iteration.   %
%                                                                             %
%  On output:                                                                 %
%       A           = the zero of function g(A) closest to Aguess.            %
%       iter        = iteration performed                                     %
%                                                                             %
%=============================================================================%
function [A,iter] = findA( Aguess, delta, phi0, tol )
  A = Aguess ;
  for iter=1:100
    [intC,intS] = GeneralizedFresnelCS( 3, 2*A, delta-A, phi0 ) ;
    f  = intS(1) ;
    df = intC(3)-intC(2) ;
    A  = A - f/df ;
    if abs(f) < tol
      break ;
    end
  end
  if abs(f) > tol*10
    fprintf( 1, 'Newton iteration fails, f = %g\n', f ) ;
    fprintf( 1, 'Aguess = %g, A = %g, delta = %g , phi0 = %g\n', Aguess, A, delta, phi0 ) ;
  end
end
%=============================================================================%
%  guessA:  Find guess for zeros of function g(A)                             %
%                                                                             %
%  USAGE:  A = guessA( phi0, phi1 );                                          %
%                                                                             %
%  On input:                                                                  %
%       phi0, phi1 = Angles used in the clothoid fitting problem.             %
%                                                                             %
%  On output:                                                                 %
%       A = an approximate zero of function g(A).                             %
%                                                                             %
%=============================================================================%
function A = guessA( phi0, phi1 )
  CF = [ 2.989696028701907, ...
         0.716228953608281, ...
        -0.458969738821509, ...
        -0.502821153340377, ...
         0.261062141752652, ...
        -0.045854475238709 ] ;
  X  = phi0/pi ;
  Y  = phi1/pi ;
  xy = X*Y ;
  A  = (phi0+phi1) * ( CF(1) + xy * ( CF(2) + xy * CF(3) ) + ...
                       (CF(4) + xy * CF(5)) * (X^2+Y^2) + ...
                        CF(6) * (X^4+Y^4) ) ;
end
