/*
C++ functions to use in the package
*/
#include <iostream> 
#include <cmath>
#include "functions.h"

// exponential function
double exponential(double num, double exponent){

  //declare numexp variable
  double numexp=pow(num,exponent);

  std::cout << "testing C++ functionality\n";
  //std::cout << "test result = "<< numexp <<"\n";
  return numexp;
}

// main function
//int main(){
//  std::cout << "testing C++ functionality\n";
//  return 0;
//}