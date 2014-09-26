%=============================================================================%
%  intXY:  Compute G1 interpolating curve                                     %
%                                                                             %
%  USAGE: [X,Y] = G1spline( PNTS, varargin ) ;                                %
%                                                                             %
%   compute clotoids segments which passes for points PNTS and joins          %
%   with G1 continuity. A minimization process is performed in order          % 
%   to minimize the jump of curvature                                         %
%                                                                             %
%  On input:                                                                  %
%                                                                             %
%       PNTS = interpolation points                                           %
%                                                                             %
%  On output:                                                                 %
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
function [theta,k,dk,L,nevalG1,nevalF,iter,Fvalue,Fgradnorm] = G1spline( PNTS, varargin )
  global G1SPLINE_PNTS G1SPLINE_NBUILD_CLOTHOID G1SPLINE_NBUILD_CLOTHOID ;
  G1SPLINE_PNTS            = PNTS ;
  G1SPLINE_NBUILD_CLOTHOID = 0 ;
  %
  % Compute guess angles
  %  
  theta          = zeros(size(PNTS,1),1) ;
  theta(1)       = atan2( PNTS(2,2)-PNTS(1,2), PNTS(2,1)-PNTS(1,1) ) ;
  theta(2:end-1) = atan2( PNTS(3:end,2)-PNTS(1:end-2,2), PNTS(3:end,1)-PNTS(1:end-2,1) ) ;
  theta(end)     = atan2( PNTS(end,2)-PNTS(end-1,2), PNTS(end,1)-PNTS(end-1,1) ) ;
  %
  % Compute solution using lsqnonlin
  %  
  options = optimset(varargin{:});
  [theta,resnorm,residual,exitflag,output,lambda,jacobian] = lsqnonlin(@target, theta, [], [], options ) ;
  %
  % Compute spline parameters
  %
  N  = length(theta) ;
  k  = zeros(N-1,1) ;
  dk = zeros(N-1,1) ;
  L  = zeros(N-1,1) ;
  for j=1:N-1
    xL = PNTS(j,1)   ; yL = PNTS(j,2)   ; tL = theta(j) ;
    xR = PNTS(j+1,1) ; yR = PNTS(j+1,2) ; tR = theta(j+1) ;
    [ k(j), dk(j), L(j) ] = buildClothoid( xL, yL, tL, xR, yR, tR ) ;
  end
  nevalG1   = G1SPLINE_NBUILD_CLOTHOID ;
  nevalF    = output.funcCount ;  
  iter      = output.iterations ;
  Fvalue    = sqrt(resnorm/N) ;
  Fgradnorm = output.firstorderopt ;
end
%
%
%
function RES = target( THETA )
  global G1SPLINE_PNTS G1SPLINE_NBUILD_CLOTHOID G1SPLINE_NBUILD_CLOTHOID ;
  N   = length(THETA) ;
  k   = zeros(N-1,1) ;
  dk  = zeros(N-1,1) ;
  L   = zeros(N-1,1) ;
  for j=1:N-1
    xL = G1SPLINE_PNTS(j,1)   ; yL = G1SPLINE_PNTS(j,2)   ; tL = THETA(j) ;
    xR = G1SPLINE_PNTS(j+1,1) ; yR = G1SPLINE_PNTS(j+1,2) ; tR = THETA(j+1) ;
    [ k(j), dk(j), L(j), iter ] = buildClothoid( xL, yL, tL, xR, yR, tR ) ;
  end
  G1SPLINE_NBUILD_CLOTHOID = G1SPLINE_NBUILD_CLOTHOID + N-1 ;
  kL  = k+dk.*L ;
  RES = [ dk(1) ; k(2:end)-kL(1:end-1) ; dk(end) ] ;
end
