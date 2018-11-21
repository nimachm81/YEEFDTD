

#ifndef FDTD_NUMBERTYPES_H_
#define FDTD_NUMBERTYPES_H_

// Defines floating point data type used in the FDTD algorithm
// and the integer types that describe indices for large arrays.
// Since the size of the arrays may pass INT_MAX, int64_t is provided as the
// int type.
// "RealNumber = double" provides higher accuracy, while "RealNumber = float"
// provides memory efficiency and it can be vectorized more efficiently.

#include <cstdint>    // std::int64_t
#include <complex>

using RealNumber = double;
using BigIntNumber = std::int64_t;


#endif  // FDTD_NUMBERTYPES_H_



