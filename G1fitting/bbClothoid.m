%======================================================================%
%  bbClothoid:  Compute a series of bounding triangles for             %
%               a clothoid curve.                                      %
%                                                                      %
%  USAGE:                                                              %
%    TT = bbClothoid( x0, y0, theta0, k, dk, L, angle, size ) ;        %
%    TT = bbClothoid( x0, y0, theta0, k, dk, L, angle, size, offs ) ;  %
%                                                                      %
%  On input:                                                           %
%                                                                      %
%  x0, y0 = coodinate of initial point                                 %
%  theta0 = orientation (angle) of the clothoid at initial point       %
%  k      = curvature at initial point                                 %
%  dk     = derivative of curvature respect to arclength               %
%  L      = the lenght of the clothoid curve or a vector of length     %
%           where to compute the clothoid values                       %
%  angle  = maximum variation of angle in the bounding box triangle    %
%  size   = maximum height of the bounding box triangle                %
%                                                                      %
%  On output:                                                          %
%                                                                      %
%  TT = matrix 6 x N whose column are the coordinates of bounding      %
%       triangles.                                                     %
%       TT(:,i) = [ x0, y0, x1, y1, x2, y2 ].'                         %
%                                                                      %
%======================================================================%
%                                                                      %
%  Autor: Enrico Bertolazzi                                            %
%         Department of Industrial Engineering                         %
%         University of Trento                                         %
%         enrico.bertolazzi@unitn.it                                   %
%                                                                      %
%======================================================================%
function TT = bbClothoid( x0, y0, theta0, k, dk, L, angle, size, offs )
  error('bbClothoid is availble only as compiled mex file') ;
end
