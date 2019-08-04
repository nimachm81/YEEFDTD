

#ifndef _FDTD_CYLINDER_H_
#define _FDTD_CYLINDER_H_

#include <vector>
#include <array>
#include <cstddef>

#include "NumberTypes.h"
#include "Geometry.h"

class Cylinder : public Geometry {
    public:
    virtual ~Cylinder() { };
    void SetRadius(const FPNumber r);
    void SetTopPlane(const std::array<FPNumber, 3> center);
    void SetHeight(const FPNumber h);
    void SetEvenAlignment(bool alignEven);  // aligns the edges on even (%2==0) grid points

    bool IsPointInsideOrOn(std::array<FPNumber, 3> point);
    void ArePointsInsideOrOn(std::vector<std::array<FPNumber, 3>>& points,
                             std::vector<bool>& areInside);
    void AreGridPointsInsideOrOn(const NumberArray3D<FPNumber>& gridArray,
                                 const std::array<FPNumber, 3>& r0,
                                 const std::array<FPNumber, 3>& dr,
                                 std::int8_t*** areInside
                                 );

    void SubdevideSurface2D(FPNumber x_cut,     // sbdeviides the yz cross section at x = x_cut to
                            FPNumber maxArcLength,      // maximum length of each section
                            std::vector<std::array<FPNumber, 3>>& centerPoints,     // the point at the middle of each patch
                            std::vector<std::array<FPNumber, 3>>& normalVecs,   // normal unit vector pointing outwards
                            std::vector<FPNumber>& arcLenghts,      // patch area
                            std::vector<std::vector<std::array<FPNumber, 3>>>* arcSubdivisionPoints = nullptr,    // if defined, provides
                            std::size_t numSubSubdivisionPoints = 1
                            );
    void GetBoundingBox2D(FPNumber x_cut,
                          std::array<FPNumber, 3>& lowerLeftCorner,
                          std::array<FPNumber, 3>& upperRightCorner
                          );

    virtual void SubdevideSurface(
                            FPNumber maxElementSurfaceArea,      // with maximum length maxArcLength
                            std::vector<std::array<FPNumber, 3>>& centerPoints,     // output: the point at the middle of each patch
                            std::vector<std::array<FPNumber, 3>>& normalVecs,   // output: normal unit vector pointing outwards
                            std::vector<FPNumber>& elementSurfaceAreas,      // output: patch area
                            std::vector<std::vector<std::array<FPNumber, 3>>>* surfaceSubdivisionPoints = nullptr,    // if defined, provides
                            std::size_t numSubSubdivisionPoints = 1     // this many equaly spaced points on each subsection
                            );

    virtual void GetBoundingBox(
                            std::array<FPNumber, 3>& lowerLeftCorner,
                            std::array<FPNumber, 3>& upperRightCorner
                            );


    private:
    FPNumber radius = 0.0;   // for tipRadius > 0 the tip is rounded
    std::array<FPNumber, 3> topPlaneCenter;    // position of the unrounded sharp tip
    FPNumber height;    // position of the rounded tip (apex)
    bool align_even = false;
};

#endif // FDTD_CONE_H_
