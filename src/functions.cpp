/*
C++ functions to use in the package
*/

#include <iostream> 
#include <cmath>
#include "functions.h"

using namespace std;

// exponential function
double exponential(double num, double exponent){

  //declare numexp variable
  double numexp=pow(num,exponent);

  cout << "testing C++ functionality\n";
  cout << "test result = "<< numexp <<"\n";
  return numexp;
}

// main function
int main(){
  cout << "testing C++ functionality\n";
  return 0;
}