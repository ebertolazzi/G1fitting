%=============================================================================%
%  intXY:  Compute piecewise clothoid spline                                  %
%                                                                             %
%  USAGE: [X,Y] = GeneralizedFresnelCS( nk, a, b, c ) ;                       %
%                                                                             %
%  Integrals are defined as:                                                  %
%                                                                             %
%  X_k(a,b,c) = \int_0^1 t^k * cos( (a/2)*t^2 + b*t + c ) dt                  %
%  Y_k(a,b,c) = \int_0^1 t^k * sin( (a/2)*t^2 + b*t + c ) dt                  %
%                                                                             %
%  On input:                                                                  %
%                                                                             %
%       nk      = number of momentae to be computed                           %
%       a, b, c = the parameters of the integrals                             %
%                                                                             %
%  On output:                                                                 %
%                                                                             %
%       X = vector with Fresnel cosine momenta [X_0,X_1,...,X_{nk-1}]         %
%       Y = vector with Fresnel sine momenta   [Y_0,Y_1,...,Y_{nk-1}]         %
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
function [theta,k,dk,L,nevalG1,nevalG1fail,nevalG1deg,nevalF,iter,Fvalue,Fgradnorm] = G1splineMW( PNTS, varargin )
  global G1SPLINE_MW_PNTS G1SPLINE_MW_TOL MW_fails MW_calls MW_calls_deg ;
  G1SPLINE_MW_PNTS = PNTS ;
  G1SPLINE_MW_TOL  = 1e-10 ;
  MW_fails         = 0 ;
  MW_calls         = 0 ;
  MW_calls_deg     = 0 ;
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
    [ k(j), dk(j), L(j) ] = buildClothoidMW( xL, yL, tL, xR, yR, tR, G1SPLINE_MW_TOL ) ;
  end
  nevalG1      = MW_calls ;
  nevalG1fail  = MW_fails ;
  nevalG1deg   = MW_calls_deg ;
  nevalF       = output.funcCount ;  
  iter         = output.iterations ;
  Fvalue       = sqrt(resnorm/N)  ;
  Fgradnorm    = output.firstorderopt ;
end
%
%
%
function RES = target( THETA )
  global G1SPLINE_MW_PNTS G1SPLINE_MW_TOL G1SPLINE_MW_NBUILD_CLOTHOID G1SPLINE_MW_NBUILD_CLOTHOID ;
  N   = length(THETA) ;
  k   = zeros(N-1,1) ;
  dk  = zeros(N-1,1) ;
  L   = zeros(N-1,1) ;
  for j=1:N-1
    xL = G1SPLINE_MW_PNTS(j,1)   ; yL = G1SPLINE_MW_PNTS(j,2)   ; tL = THETA(j) ;
    xR = G1SPLINE_MW_PNTS(j+1,1) ; yR = G1SPLINE_MW_PNTS(j+1,2) ; tR = THETA(j+1) ;
    [ k(j), dk(j), L(j), iter ] = buildClothoidMW( xL, yL, tL, xR, yR, tR, G1SPLINE_MW_TOL ) ;
  end
  G1SPLINE_MW_NBUILD_CLOTHOID = G1SPLINE_MW_NBUILD_CLOTHOID + N-1 ;
  kL  = k+dk.*L ;
  RES = [ dk(1) ; k(2:end)-kL(1:end-1) ; dk(end) ] ;
end
