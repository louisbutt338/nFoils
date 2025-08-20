/*
Header file for functions
include smart pointers and extern C to stop data C++ wranglings
*/

// to import in a python module via ctypes:
// first do "cd nFoils" "cmake -S src -B build" "cmake --build build" in cmd line
// import ctypes
// cpp_functions = ctypes.CDLL("/Users/ljb841@student.bham.ac.uk/nFoils/build/libfunctions.so")
// cpp_poly_eff = cpp_functions.polynomial_efficiency
// cpp_poly_eff.argtypes = [ctypes.c_double,ctypes.c_double,ctypes.c_double, ctypes.c_double,ctypes.c_double]
// cpp_poly_eff.restype = ctypes.c_double

// to stop C++ wrangling
extern "C" {

// put the functions you want in here
double polynomial_efficiency(double e, double a0,double a1,double a2,double a3);

} // extern "C"