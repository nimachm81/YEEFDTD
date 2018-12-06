

#ifndef FDTD_NUMBERTYPES_H_
#define FDTD_NUMBERTYPES_H_

// Defines floating point data type used in the FDTD algorithm
// and the integer types that describe indices for large arrays.
// Since the size of the arrays may pass INT_MAX, int64_t is provided as the
// int type.
// "FPNumber = double" provides higher accuracy, while "FPNumber = float"
// provides memory efficiency and it can be vectorized more efficiently.


#include <cstdint>    // std::int64_t
#include <complex>

using FPNumber = std::complex<double>;  // floating point number. could be float, double, complex<float>, complex<double>
using BigIntNumber = std::int64_t;


#endif  // FDTD_NUMBERTYPES_H_



