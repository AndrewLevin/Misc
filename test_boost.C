#include "TLorentzVector.h"
#include <iostream>

void test_boost ()
{

  TLorentzVector v1(TVector3(1.,2.,3.),100.);
  TLorentzVector v2(TVector3(5.,-9.,3.),100.);

  TVector3 b(v1.BoostVector());

  std::cout << v2.Px() << std::endl;

  std::cout << v1.Angle(v2.Vect()) << std::endl;
  
}
