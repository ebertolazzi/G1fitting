%
%  Autors: Enrico Bertolazzi and Marco Frego
%          Universita' degli Studi di Trento
%          enrico.bertolazzi@unitn.it
%          m.fregox@gmail.com
%
%******************************************************************************
% Driver test program to check clothoid computation                           %
%******************************************************************************

addpath('../G1fitting');
addpath('../meekAndWalton');
close all ;

tol     = 1E-14 ;
maxiter = 1000 ;

% test 1
tests(1).P1 = [5;4] ;
tests(1).a1 = pi/3 ;
tests(1).P2 = [5;6] ;
tests(1).a2 = 7*pi/6 ;

% test 2
tests(2).P1 = [3;5] ;
tests(2).a1 = 2.14676 ;
tests(2).P2 = [6;5] ;
tests(2).a2 = 2.86234 ;

% test 3
tests(3).P1 = [3;6] ;
tests(3).a1 = 3.05433 ;
tests(3).P2 = [6;6] ;
tests(3).a2 = 3.14159 ;

% test 4
tests(4).P1 = [3;6] ;
tests(4).a1 = 0.08727 ;
tests(4).P2 = [6;6] ;
tests(4).a2 = 3.05433 ;

% test 5
tests(5).P1 = [5;4] ;
tests(5).a1 = 0.34907 ;
tests(5).P2 = [4;5] ;
tests(5).a2 = 4.48550 ;

% test 6
tests(6).P1 = [5;5] ;
tests(6).a1 = 0.52360 ;
tests(6).P2 = [4;4] ;
tests(6).a2 = 4.66003 ;

% test 7
tests(7).P1 = [0;0] ;
tests(7).a1 = 0 ;
tests(7).P2 = [1;0] ;
tests(7).a2 = pi/10000 ;

% test 8
tests(8).P1 = [0;0] ;
tests(8).a1 = 0 ;
tests(8).P2 = [sqrt(2)/2;sqrt(2)/2] ;
tests(8).a2 = pi/2-1e-8 ;

npts = 100 ;

disp('-----------------------------------') ;
disp('  Test N1,...,N6 of the manuscript ') ;
disp('-----------------------------------') ;

for kk=1:6

  P1 = tests(kk).P1 ;
  a1 = tests(kk).a1 ;
  P2 = tests(kk).P2 ;
  a2 = tests(kk).a2 ;

  [P0,T0,N0,a,t1,t2,iter,failFlag,reflectFlag,reverseFlag,kind] = ...
    completeShape( P1, [cos(a1);sin(a1)], P2, [cos(a2);sin(a2)], tol, maxiter ) ;

  subplot(3,2,kk) ;
  XY = pointsOnShape( P0, T0, N0, a, t1-0.1*(t2-t1), t2+0.1*(t2-t1), npts ) ;
  plot( XY(1,:), XY(2,:), '-r' ) ;
  title(sprintf('Test N.%d',kk));
  hold on ;
  
  % save data on files
  % fileID = fopen(sprintf('TestN%d_MeekWalton.dat',kk),'w');
  % fprintf(fileID,'%12s\t%12s\n','x','y');
  % fprintf(fileID,'%12.8f\t%12.8f\n',XY(:,1:10:end));
  % fclose(fileID);

  % 
  XY = pointsOnShape( P0, T0, N0, a, t1, t2, npts ) ;
  
  disp('-------------------') ;
  
  if reverseFlag
    e1 = XY(:,1)-P2 ;
    e2 = XY(:,end)-P1 ;
  else
    e1 = XY(:,1)-P1 ;
    e2 = XY(:,end)-P2 ;
  end
  
  fprintf('Meek & Walton,      iter %2d error %10g %10g\n', iter, norm(e1,1), norm(e2,1) ) ; 
  
  [k,dk,Lsol,iter] = buildClothoid( P1(1), P1(2), a1, P2(1), P2(2), a2 ) ;
  XY               = pointsOnClothoid( P1(1), P1(2), a1, k, dk, Lsol, npts ) ;
  plot( XY(1,:), XY(2,:), '-b', 'LineWidth', 2 ) ;

  axis equal
  
  ee1 = XY(:,1)-P1 ;
  ee2 = XY(:,end)-P2 ;
  fprintf('Bertolazzi & Frego, iter %2d error %10g %10g\n', iter, norm(ee1,1), norm(ee2,1) ) ; 

  % save data ion files
  % fileID = fopen(sprintf('TestN%d.dat',kk),'w');
  % fprintf(fileID,'%12s\t%12s\n','x','y');
  % fprintf(fileID,'%12.8f\t%12.8f\n',XY(:,1:10:end));
  % fclose(fileID);
  
end
