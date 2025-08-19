/*
C++ functions to use in the package
*/
#include <iostream> 
#include <cmath>
#include "functions.h"

// hpge efficiency function
double polynomial_efficiency(double e, double a0,double a1,double a2,double a3){

  //declare terms of the equation
  double term0 = a0;
  double term1 = a1*std::pow(std::log(e),1);
  double term2 = a2*std::pow(std::log(e),2);
  double term3 = a3*std::pow(std::log(e),3);

  //declare the final value
  double poly_term = std::exp(term0+term1+term2+term3);

  std::cout << "testing C++ eff fn \n";
  //std::cout << "test result = "<< numexp <<"\n";
  return poly_term;
}

// main function
//int main(){
//  std::cout << "testing C++ functionality\n";
//  return 0;
//}