%
%  Autors: Enrico Bertolazzi
%          Universita' degli Studi di Trento
%          enrico.bertolazzi@unitn.it
%          m.fregox@gmail.com
%
function XY = pointsOnShape( P0, T0, N0, a, t1, t2, npts )
  XY = [] ;
  for t=[t1:(t2-t1)/npts:t2]
    [C,S] = FresnelCS(t) ;
    Q     = P0 + a*(C*T0+S*N0) ;
    XY    = [XY Q] ;
  end
end
