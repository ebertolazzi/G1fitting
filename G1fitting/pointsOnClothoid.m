%=============================================================================%
% pointsOnClothoid:  Compute points on a clothoid curve.                      %
%                    Used for plotting purpose.                               %
%                                                                             %
% USAGE: XY = pointsOnClothoid( x0, y0, theta0, kappa, dkappa, L, npts ) ;    %
%                                                                             %
% On input:                                                                   %
%                                                                             %
%      x0, y0  = coodinate of initial point                                   %
%      theta0  = orientation (angle) of the clothoid at initial point         %
%      kappa   = curvature at initial point                                   %
%      dkappa  = derivative of curvature respect to arclength                 %
%      L       = the lenght of the clothoid curve from initial to final point %
%      npts    = number of points along the clothoid                          %
%                                                                             %
% On output:                                                                  %
%                                                                             %
%      XY      = matrix 2 x NPTS whose column are the points of the clothoid  %
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
function varargout = pointsOnClothoid( x0, y0, theta0, kappa, dkappa, L, npts )
  X = [] ;
  Y = [] ;
  tvec = [0:L/npts:L] ;
  for t=tvec
    [C,S] = GeneralizedFresnelCS( 1, dkappa*t^2, kappa*t, theta0 ) ;
    X = [ X x0 + t*C ] ;
    Y = [ Y y0 + t*S ] ;
  end
  if nargout > 1
    varargout{1} = X ;
    varargout{2} = Y ;
  else
    varargout{1} = [X ; Y] ;
  end
end
