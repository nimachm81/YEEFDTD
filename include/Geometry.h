#ifndef FDTD_GEOMETRY_H_
#define FDTD_GEOMETRY_H_

#include <array>
#include <vector>

#include "NumberTypes.h"
#include "MultiDimArray.hpp"

class Geometry {
    public:
    virtual ~Geometry() { };

    virtual bool IsPointInsideOrOn(std::array<FPNumber, 3> point) = 0;
    virtual void AreGridPointsInsideOrOn(const NumberArray3D<FPNumber>& gridArray,
                                         const std::array<FPNumber, 3>& r0,
                                         const std::array<FPNumber, 3>& dr,
                                         std::int8_t*** areInside
                                         ) = 0;
    virtual void SubdevideSurface2D(FPNumber x_cut,     // sbdeviides the  y-z cross section at x = x_cut to small devisions
                            FPNumber maxArcLength,      // with maximum length maxArcLength
                            std::vector<std::array<FPNumber, 3>>& centerPoints,     // output: the point at the middle of each patch
                            std::vector<std::array<FPNumber, 3>>& normalVecs,   // output: normal unit vector pointing outwards
                            std::vector<FPNumber>& arcLenghts      // output: patch area
                            ) = 0;

};

#endif // FDTD_GEOMETRY_H_


