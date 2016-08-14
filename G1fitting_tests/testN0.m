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

addpath('../G1fitting');

close all ;

%tol  = 1e-12 ;
npts = 400 ;

% initial point with angle direction
x0     = 0.5 ;
y0     = 0   ;
theta0 = 1   ;

% final point with angle direction
x1     = 0.5 ;
y1     = 0.5 ;
theta1 = 1   ;

fprintf('Testing G1 Clothoid interpolation\n') ;
fprintf('initial point (%g,%g) initial angle = %g\n', x0, y0, theta0) ;
fprintf('final point (%g,%g) final angle = %g\n', x1, y1, theta1) ;

% compute clothoid parameters
[k,dk,L] = buildClothoid( x0, y0, theta0, x1, y1, theta1 ) ;

fprintf('Computed parameters: k = %g, k'' = %g, L = %g\n', k, dk, L ) ;

% compute points on clothoid
XY = pointsOnClothoid( x0, y0, theta0, k, dk, L, npts ) ;

plot( XY(1,:), XY(2,:), '-r' ) ;

axis equal
