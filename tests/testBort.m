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

npts = 400 ;

%=============================================================================%
P0 = [0.000, 0.000] ; P1 = [0.100, 0.000] ;
X = [ P0(1) P1(1) ] ;
Y = [ P0(2) P1(2) ] ;
plot( X, Y, '-r' ) ;
axis equal ;
hold on ;

%=============================================================================%
P0 = [0.100, 0.000] ; P1 = [0.085, 0.035] ; ANG = [1.571, 2.356] ;
[k,dk,L] = buildClothoid( P0(1), P0(2), ANG(1), P1(1), P1(2), ANG(2) ) ;
fprintf('Computed parameters: k = %g, k'' = %g, L = %g\n', k, dk, L ) ;
XY = pointsOnClothoid( P0(1), P0(2), ANG(1), k, dk, L, npts ) ;
plot( XY(1,:), XY(2,:), '-b' ) ;

%=============================================================================%
P0 = [0.085, 0.035] ; P1 = [0.007, 0.025] ; ANG = [2.356, 4.189] ;
[k,dk,L] = buildClothoid( P0(1), P0(2), ANG(1), P1(1), P1(2), ANG(2) ) ;
fprintf('Computed parameters: k = %g, k'' = %g, L = %g\n', k, dk, L ) ;
XY = pointsOnClothoid( P0(1), P0(2), ANG(1), k, dk, L, npts ) ;
plot( XY(1,:), XY(2,:), '-r' ) ;

%=============================================================================%
P0 = [0.007, 0.025] ; P1 = [0.100, 0.000] ; ANG = [4.189, 1.571] ;
[k,dk,L] = buildClothoid( P0(1), P0(2), ANG(1), P1(1), P1(2), ANG(2) ) ;
fprintf('Computed parameters: k = %g, k'' = %g, L = %g\n', k, dk, L ) ;
XY = pointsOnClothoid( P0(1), P0(2), ANG(1), k, dk, L, npts ) ;
plot( XY(1,:), XY(2,:), '-b' ) ;

%=============================================================================%
P0 = [0.100, 0.000] ; P1 = [0.200, 0.100] ;
X = [ P0(1) P1(1) ] ;
Y = [ P0(2) P1(2) ] ;
plot( X, Y, '-r' ) ;

%=============================================================================%
P0 = [0.200, 0.100] ; P1 = [0.240, 0.120] ; ANG = [-4.249, -7.390] ;
[k,dk,L] = buildClothoid( P0(1), P0(2), ANG(1), P1(1), P1(2), ANG(2) ) ;
fprintf('Computed parameters: k = %g, k'' = %g, L = %g\n', k, dk, L ) ;
XY = pointsOnClothoid( P0(1), P0(2), ANG(1), k, dk, L, npts ) ;
plot( XY(1,:), XY(2,:), '-b' ) ;

%=============================================================================%
P0 = [0.240, 0.120] ; P1 = [0.200, 0.100] ; ANG = [-7.390, -4.249] ;
[k,dk,L] = buildClothoid( P0(1), P0(2), ANG(1), P1(1), P1(2), ANG(2) ) ;
fprintf('Computed parameters: k = %g, k'' = %g, L = %g\n', k, dk, L ) ;
XY = pointsOnClothoid( P0(1), P0(2), ANG(1), k, dk, L, npts ) ;
plot( XY(1,:), XY(2,:), '-r' ) ;

%=============================================================================%
P0 = [0.200, 0.100] ; P1 = [0.000, 0.000] ;
X = [ P0(1) P1(1) ] ;
Y = [ P0(2) P1(2) ] ;
plot( X, Y, '-b' ) ;
