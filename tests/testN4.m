%=============================================================================%
%                                                                             %
%  Autors: Enrico Bertolazzi and Marco Frego                                  %
%          Department of Industrial Engineering                               %
%          University of Trento                                               %
%          enrico.bertolazzi@unitn.it                                         %
%          m.fregox@gmail.com                                                 %
%                                                                             %
%=============================================================================%
% Driver test program to check clothoid computation                           %
%=============================================================================%

function testN4

  addpath('../G1fitting') ;

  mind   = 100 ;
  maxINT = 0 ;
  NN     = 128 ; % 1024 ;

  for phi1=[-pi*0.9999:pi/NN:pi*0.9999]
    for phi0=[-pi*0.9999:pi/NN:pi*0.9999]
      delta = phi1 - phi0 ;
      % initial point
      Aguess = guessA( phi0, phi1 ) ;

      % Newton iteration
      [A,iter] = findA( guessA( phi0, phi1 ), delta, phi0, 1e-10 ) ;
      
      [intC,intS] = GeneralizedFresnelCS( 3, 2*A, delta-A, phi0 ) ;
      g  = intS(1) ;
      dg = intC(3)-intC(2) ;
      mind   = min([mind abs(dg)]) ;
      maxINT = max([maxINT abs(Aguess-A)]) ;
    end
  end
  fprintf(1,'Minium value of derivative of g at solution %g\n',mind) ;
  fprintf(1,'Maximum distance from guess and solution of interpolation problem %g\n',maxINT) ;
  
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
