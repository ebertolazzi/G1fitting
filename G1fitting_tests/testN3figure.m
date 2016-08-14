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
addpath('../G1fitting') ;

tol = 1e-10 ;

xx1 = [ -10, -7, -4, -3, -2, -1,   0, 1, 2, 3, 4, 7, 10 ] ; %--------------------
yy1 = [   0,  0,  0,  0,  0,  0, 0.5, 1, 1, 1, 1, 1, 1  ] ;

tt2 = pi:pi/16:(3/2)*pi ;
xx2 = [ cos(tt2), 0.1:0.1:0.9, 1-sin(tt2) ] ;
yy2 = [ sin(tt2), -ones(1,9), cos(tt2) ] ;
  
xx3 = [ cos(0:pi/4:2*pi)+1e-7*cos(0:pi/8:pi)] ; % MW crisi
yy3 = [ sin(0:pi/4:2*pi) ] ;


xx4 = [ 0:0.05:2*pi ] ;
yy4 = 1e-5*sin(xx4) ;

close all ;

PPP{1} = [xx1',yy1'] ;
PPP{2} = [xx2',yy2'] ;
PPP{3} = [xx3',yy3'] ;
PPP{4} = [xx4',yy4'] ;

mw = {'','MW'} ;

our_or_MW = {'use Bertolazzi and Frego algorithm','use Meek and Walton algorithm'} ;

opt = {'Display', 'off', 'TolX', tol, 'TolFun', tol, 'Algorithm', {'levenberg-marquardt',.01} } ;

for kk=3:3
  for jj=0:1

    fprintf('\n\n\nTest N.%d [%s]\n',kk,our_or_MW{jj+1} ) ;    
    figure('WindowStyle','docked');

    PNTS = PPP{kk} ;
    if jj==0
      [theta,k,dk,L,nevalG1,nevalF,iter,Fvalue,Fgradnorm] = G1spline( PNTS, opt{:} ) ;
    else
      [theta,k,dk,L,nevalG1,nevalF,iter,Fvalue,Fgradnorm] = G1splineMW( PNTS, opt{:} ) ;
    end

    fprintf('G1 clothoid evaluations = %d\n', nevalG1) ;
    fprintf('F evaluations           = %d\n', nevalF) ;
    fprintf('iterations              = %d\n', iter) ;
    fprintf('last F value            = %d\n', Fvalue) ;
    fprintf('last grad F norm        = %d\n', Fgradnorm) ;

    N      = length(theta) ;
    colors = {'-r', '-g'} ;

    subplot(2,1,1) ;
    hold on
    ll = 1 ;
    % save data on files
    fileID = fopen(sprintf('AppTestN%s%d.dat',mw{jj+1},kk),'w');
    fprintf(fileID,'%15s\t%15s\n','x','y');
    for j=1:N-1
      XY = pointsOnClothoid( PNTS(j,1), PNTS(j,2), theta(j), k(j), dk(j), L(j), 100 ) ;
      plot( XY(1,:), XY(2,:), colors{ll}, 'LineWidth', 2 ) ;
      fprintf(fileID,'%20.10g\t%20.10g\n',XY(:,1:10:end));
      ll = 3-ll ;
    end
    fclose(fileID);
    axis equal ;
 
    subplot(2,1,2) ;
    hold on ;
    ll = 1 ;
    xx = 0 ;
    % save data on files
    fileID = fopen(sprintf('AppTestN%s%d_curv.dat',mw{jj+1},kk),'w');
    fprintf(fileID,'%15s\t%15s\n','x','k');
    for j=1:N-1
      plot( [xx, xx+L(j)], [k(j), k(j)+dk(j)*L(j)], colors{ll}, 'LineWidth', 2 ) ;
      fprintf(fileID,'%20.10g\t%20.10g\n',xx,k(j));
      fprintf(fileID,'%20.10g\t%20.10g\n\n',xx+L(j),k(j)+dk(j)*L(j));
      l = 3-ll ;
      xx = xx+L(j) ;
    end
    fclose(fileID);
  end
end
