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

tol     = 1E-12 ;
maxIter = 100 ;
npts    = 400 ;

R  = 100 ;
P1 = [0;0] ;
a1 = 0 ;
P2 = [100;0] ;

disp('----------------------------') ;
disp('  Test N7 of the manuscript ') ;
disp('----------------------------') ;

for kk=1:10
    
  fprintf('\nCase k = %d\n',kk);

  a1 =  0.01/2^kk ;
  a2 = -0.02/2^kk ;

  [P0,T0,N0,a,t1,t2,iter,failFlag,reflectFlag,reverseFlag,kind] = ...
    completeShape( P1, [cos(a1);sin(a1)], P2, [cos(a2);sin(a2)], tol, maxIter ) ;

  XY = pointsOnShape( P0, T0, N0, a, t1, t2, npts ) ;
  
  if kind == 0
    XY = pointsOnShape( P0, T0, N0, a, t1, t2, npts ) ;
    if reverseFlag
      e1 = XY(:,1)-P2 ;
      e2 = XY(:,end)-P1 ;
    else
      e1 = XY(:,1)-P1 ;
      e2 = XY(:,end)-P2 ;
    end
  elseif kind == 1 % circle
    dk = 0 ;
    d  = sqrt(sum((P1-P2).^2)) ;
    k  = 2*sin(a)/d ;
    L  = a*d/sin(a) ;
    if reverseFlag
      k = -k ;
    end
    if reflectFlag
      k = -k ;
    end
    e1 = [0;0] ;
    [Xe, Ye] = evalClothoid( P1(1), P1(2), a1, k, dk, Lsol ) ;
    e2 = [Xe; Ye]-P2 ;

  else % straight line
    dk = 0 ;
    k  = 0 ;
    L  = sqrt(sum((P1-P2).^2)) ;
    e1 = [0;0] ;
    [Xe, Ye] = evalClothoid( P1(1), P1(2), a1, k, dk, Lsol ) ;
    e2 = [Xe; Ye]-P2 ;
  end
  
  fprintf('Meek & Walton,      iter %3d error %10g %10g [kind %d]\n', iter, norm(e1,2), norm(e2,2),kind ) ; 
  
  [k,dk,Lsol,iter] = buildClothoid( P1(1), P1(2), a1, P2(1), P2(2), a2 ) ;
  [Xe,Ye]          = evalClothoid( P1(1), P1(2), a1, k, dk, Lsol ) ;
  
  ee2 = [Xe;Ye]-P2 ;
  fprintf('Bertolazzi & Frego, iter %3d error %10g (only end point)\n', iter, norm(ee2,2) ) ; 
  
end
