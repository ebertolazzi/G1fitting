%======================================================================%
% pointsOnClothoid:  Compute points on a clothoid curve.               %
%                    Used for plotting purpose.                        %
%                                                                      %
% USAGE: XY    = pointsOnClothoid( x0, y0, theta0, k, dk, L, npts ) ;  %
% USAGE: [X,Y] = pointsOnClothoid( x0, y0, theta0, k, dk, L, npts ) ;  %
% USAGE: XY    = pointsOnClothoid( x0, y0, theta0, k, dk, L ) ;        %
% USAGE: [X,Y] = pointsOnClothoid( x0, y0, theta0, k, dk, L ) ;        %
%                                                                      %
% On input:                                                            %
%                                                                      %
%  x0, y0  = coodinate of initial point                                %
%  theta0  = orientation (angle) of the clothoid at initial point      %
%  k       = curvature at initial point                                %
%  dk      = derivative of curvature respect to arclength              %
%  L       = the lenght of the clothoid curve or a vector of length    %
%            where to compute the clothoid values                      %
%  npts    = number of points along the clothoid                       %
%                                                                      %
% On output: (1 argument)                                              %
%                                                                      %
%  XY = matrix 2 x NPTS whose column are the points of the clothoid    %
%                                                                      %
% On output: (2 argument)                                              %
%                                                                      %
%   X  = matrix 1 x NPTS X coordinate of points of the clothoid        %
%   Y  = matrix 1 x NPTS Y coordinate of points of the clothoid        %
%                                                                      %
%======================================================================%
%                                                                      %
%  Autor: Enrico Bertolazzi                                            %
%         Department of Industrial Engineering                         %
%         University of Trento                                         %
%         enrico.bertolazzi@unitn.it                                   %
%                                                                      %
%======================================================================%
function varargout = pointsOnClothoid( x0, y0, theta0, kappa, dkappa, L, npts )
  X = [] ;
  Y = [] ;
  tvec = [0:L/(npts-1):L] ;
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
