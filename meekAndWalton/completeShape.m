%
% Implementation of algorithm described in
%
%  An improved Euler spiral algorithm for shape completion
%  by D.J.Walton and D.S.Meek
%  Canadian Conference on Computer and Robot Vision  
%  2008 IEEE
%  DOI 10.1109/CRV.2008.11
%
%
%  Autors: Enrico Bertolazzi
%          Universita' degli Studi di Trento
%          enrico.bertolazzi@unitn.it
%          m.fregox@gmail.com
%
%  Input:
%    P1         initial point
%    T1         initial direction
%    P2         end point
%    T2         end direction
%    tol        tolerance for Newton iteration
%    iterLim    maximum number of iteration allowed
%
%  Output:
%    P0
%    T0
%    N0
%    a
%    t1
%    t2
%    iter          iteration used
%    failFlag      true if computation failed
%    reflectFlag   true if the computed curve must be reflected
%    reverseFlag   true if the computed curve must be reversed
%    kind          0 = clothoid, 1 = circle, 2 = straight line
%
function [P0,T0,N0,a,t1,t2,iter,failFlag,reflectFlag,reverseFlag,kind] = completeShape( P1, T1, P2, T2, tol, iterLim )
  kind = 0 ;
  D    = P2 - P1 ;
  d    = norm(D,2) ;
  if d < tol
    error('degenerate case\n') ;
  end
  if tol < 1e-7
    tol = 1e-7 ;
  end
  phi1 = atan2( cross2(T1,D), dot2(T1,D) ) ;
  phi2 = atan2( cross2(D,T2), dot2(D,T2) ) ;

  P0   = [0;0] ;
  T0   = [0;0] ;
  N0   = [0;0] ;
  t1   = 0 ;
  t2   = 0 ;
  iter = 0 ;
  failFlag = false ;

  reverseFlag = false ;
  if abs(phi1) > abs(phi2)
    [ phi1, phi2, P1, D ] = reverse( phi1, phi2, P1, D ) ;
    reverseFlag = true ;
  end
  reflectFlag = false ;
  if phi2 < 0
    [ phi1, phi2 ] = reflect( phi1, phi2 ) ;
    reflectFlag = true ;
  end
  if ( (phi1 == 0 ) & (phi2 == pi) ) | ...
     ( (phi1 == pi) & (phi2 == 0 ) ) | ...
     ( (phi1 == pi) & (phi2 == pi) )
    error('ambiguous case %f %f\n', phi1, phi2) ;
  elseif abs(phi1) <= tol & abs(phi2) <= tol
    a    = 0 ;
    kind = 2 ; % straight line
  elseif abs(phi1-phi2) <= tol
    a    = phi1 ;
    kind = 1 ; % circle
  else
    [P0,T0,N0,a,t1,t2,iter,failFlag] = fitEuler( P1, D/d, d, phi1, phi2, tol, iterLim, reflectFlag ) ;
  end
end

function RES = cross2( A, B )
  RES = A(1) * B(2) - A(2) * B(1) ;
end

function RES = dot2( A, B )
  RES = A(1) * B(1) + A(2) * B(2) ;
end

function [phi1_out,phi2_out,P1_out,D_out] = reverse( phi1, phi2, P1, D )
  D_out    = -D ;
  P1_out   = P1 + D ;
  phi1_out = -phi2 ;
  phi2_out = -phi1 ;
end

function [phi1_out,phi2_out] = reflect( phi1, phi2 )
  phi1_out = -phi1 ;
  phi2_out = -phi2 ;
end

function [W] = rotate(V,alpha)
  W    = zeros(2,1) ;
  W(1) = V(1) * cos(alpha) - V(2) * sin(alpha) ;
  W(2) = V(1) * sin(alpha) + V(2) * cos(alpha) ;
end

function [P0,T0,N0,a,t1,t2,iter,failFlag] = fitEuler(P1,T,d,phi1,phi2,tol,iterLim,reflectFlag)
  theta = 0 ;
  segno = 1 ;
  t1    = 0 ;
  t2    = sqrt( 2*(phi1+phi2)/pi ) ;
  [C,S] = FresnelCS(t2) ;
  h     = S*cos(phi1)-C*sin(phi1) ;
  if (phi1 > 0) & (h<=0)
    % C shaped
    if h > tol
      failFlag = false ; % solution theta = 0
    else
      lambda = (1-cos(phi1))/(1-cos(phi2)) ;
      theta0 = lambda^2/(1-lambda^2)*(phi1+phi2) ;
      [theta,iter,failFlag] = solve(0,theta0,phi1,phi2,segno,tol,iterLim) ;
    end
  else
    segno  = -1 ;
    theta0 = max(0,-phi1) ;
    theta1 = pi/2-phi1 ;
    [theta,iter,failFlag] = solve(theta0,theta1,phi1,phi2,segno,tol,iterLim) ;    
  end
  
  t1      = segno*sqrt(2*theta/pi) ;
  t2      = sqrt(2*(theta+phi1+phi2)/pi) ;
  [C1,S1] = FresnelCS(t1) ;
  [C2,S2] = FresnelCS(t2) ;
  phi     = phi1 + theta ;
  a       = d/((S2-S1)*sin(phi)+(C2-C1)*cos(phi)) ;
  if reflectFlag
    T0 = rotate(T,phi) ;
    N0 = rotate(T0,-pi/2) ;
  else
    T0 = rotate(T,-phi) ;
    N0 = rotate(T0,pi/2) ;
  end
  P0 = P1 - a*(C1*T0+S1*N0) ;
end

function [C,S] = Fresnel2(theta)
  [C,S] = FresnelCS(sqrt(2*abs(theta)/pi)) ;
  C = C*sign(theta) ;
  S = S*sign(theta) ;
end

function [f,df] = evaluate(theta,phi1,phi2,segno)
  [C1,S1] = Fresnel2(theta) ;
  [C2,S2] = Fresnel2(theta+phi1+phi2) ;
  c  = cos(theta+phi1) ;
  s  = sin(theta+phi1) ;
  f  = sqrt(2*pi)*(c*(S2-segno*S1)-s*(C2-segno*C1)) ;
  df = sin(phi2)/sqrt(theta+phi1+phi2) + segno*sin(phi1)/sqrt(theta) ...
     - sqrt(2*pi)*(s*(S2-segno*S1)+c*(C2-segno*C1)) ; 
end

function [theta,iter,failFlag] = solve( a, b, phi1, phi2, segno, tol, iterLim )
  [fa,dfa] = evaluate(a,phi1,phi2,segno) ;
  [fb,dfb] = evaluate(b,phi1,phi2,segno) ;
  theta  = (a+b)/2 ;
  [f,df] = evaluate(theta,phi1,phi2,segno) ;
  err    = b-a ;
  iter   = 0 ;
  while (err > tol) & (iter < iterLim)
    iter = iter + 1 ;
    NewtonFail = true ;
    if abs(df) > tol
      thetaiter = theta - f/df ;
      delta     = abs(theta-thetaiter) ;
      if (thetaiter > a) & (thetaiter < b) & ( delta < 0.5 * err )
        NewtonFail = false ;
        theta      = thetaiter ;
        err        = delta ;
      end
    end
    if NewtonFail
      if fa * f < 0
        b  = theta ;
        fb = f ;
      else
        a  = theta ;
        fa = f ;
      end
      theta = (a+b)/2 ;
      err   = b-a ;
    end
    [f,df] = evaluate(theta, phi1, phi2, segno) ;
  end
  failFlag = iter >= iterLim ;
end
