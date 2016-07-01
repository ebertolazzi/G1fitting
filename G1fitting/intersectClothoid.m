%======================================================================%
%  intersectClothoid:  Compute intersections betweed clothoids         %
%                                                                      %
%  USAGE:                                                              %
%    s1, s2 = intersectClothoid( clot1, clot2 ) ;                      %
%                                                                      %
%  On input:                                                           %
%                                                                      %
%    clot1 = structure with the definition of the first clothoid       %
%            must contain the field                                    %
%      theta0 = orientation (angle) of the clothoid at initial point   %
%      kappa  = curvature at initial point                             %
%      dkappa = derivative of curvature respect to arclength           %
%      L      = the lenght of the clothoid curve                       %
%    clot2 = structure with the definition of the second clothoid      %
%                                                                      %
%  On output:                                                          %
%                                                                      %
%  s1 = curvilinear coordinates of intersections on clot1              %
%  s2 = curvilinear coordinates of intersections on clot2              %
%                                                                      %
%======================================================================%
%                                                                      %
%  Autor: Enrico Bertolazzi                                            %
%         Department of Industrial Engineering                         %
%         University of Trento                                         %
%         enrico.bertolazzi@unitn.it                                   %
%                                                                      %
%======================================================================%
function [s1,s2] = intersectClothoid( clot1, clot2 )
  error('intersectClothoid is availble only as compiled mex file') ;
end
