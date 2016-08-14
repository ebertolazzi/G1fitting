
#include "Clothoid.hh"
#include <iostream>

Clothoid::valueType m_pi = 3.14159265358979323846264338328 ;

using namespace std ;

int
main() {
  Clothoid::ClothoidCurve c1(0,0,m_pi*0.7,-0.1,0.15,20) ;
  Clothoid::ClothoidCurve c2(-10,0,m_pi/4,0.1,-0.05,20) ;
  Clothoid::valueType max_angle = m_pi/2 ;
  Clothoid::valueType max_size  = 0.5 ;

  std::vector<Clothoid::ClothoidCurve> c ;
  std::vector<Clothoid::Triangle2D>    t ;
  c.clear() ;
  t.clear() ;
  c1.bbSplit( max_angle, max_size, 0, c, t ) ;
  cout << "n = " << c.size() << '\n' ;
  cout << "n = " << t.size() << '\n' ;

  //c.clear() ;
  //t.clear() ;
  //c2.bbSplit( max_angle, max_size, 0, c, t ) ;
  //cout  << "n = " << c.size() << '\n' ;
  //cout  << "n = " << t.size() << '\n' ;
  /*
  for ( Clothoid::indexType i = 0 ; i < c.size() ; ++i )
    cout << c[i] << '\n' ;
  for ( Clothoid::indexType i = 0 ; i < t.size() ; ++i )
    cout << t[i] << '\n' ;
  */
  return 0 ;
}
