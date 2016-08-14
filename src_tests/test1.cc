
#include "Clothoid.hh"
#include <iostream>

Clothoid::valueType m_pi = 3.14159265358979323846264338328 ;

using namespace std ;

int
main() {
  Clothoid::ClothoidCurve c1(0,0,m_pi*0.7,-0.1,0.05,20) ;
  Clothoid::ClothoidCurve c2(-10,0,m_pi/4,0.1,-0.05,20) ;
  std::vector<Clothoid::valueType> s1, s2 ;
  Clothoid::indexType max_iter  = 10 ;
  Clothoid::valueType tolerance = 1e-8 ;
  try {
    c1.intersect( 0, c2, 0, s1, s2, max_iter, tolerance ) ;
    cout  << "ok = TRUE\n" ;
  }
  catch (...) {
    cout  << "ok = FALSE\n" ;
  }
  for ( Clothoid::indexType i = 0 ; i < s1.size() ; ++i )
    cout << "S1[" << i << "] = " << s1[i] << '\n' ;
  for ( Clothoid::indexType i = 0 ; i < s2.size() ; ++i )
    cout << "S2[" << i << "] = " << s2[i] << '\n' ;
  return 0 ;
}
