#ifndef FDTD_HYPERBOLOID_H_
#define FDTD_HYPERBOLOID_H_

#include <vector>
#include <array>
#include <cstddef>

#include "NumberTypes.h"
#include "Geometry.h"

class Hyperboloid : public Geometry {
    public:
    virtual ~Hyperboloid() { };

    void SetConeAngle(const FPNumber angle);
    void SetConeAngleInDegrees(const FPNumber angle);
    void SetApexRadius(const FPNumber radius);
    void SetApexPosition(const std::array<FPNumber, 3> pos);
    void SetHeight(const FPNumber h);

    inline void GetCanonicalScaleFators(FPNumber& a, FPNumber& b);

    virtual bool IsPointInsideOrOn(std::array<FPNumber, 3> point);
    virtual void ArePointsInsideOrOn(std::vector<std::array<FPNumber, 3>>& points,
                             std::vector<bool>& areInside);
    virtual void AreGridPointsInsideOrOn(const NumberArray3D<FPNumber>& gridArray,
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
    FPNumber coneAngle = M_PI/4.0;  // opening angle of the asymtotic cone  (full angle)
    FPNumber apexRadius = 1.0;    // radius of curvature at the tip
    std::array<FPNumber, 3> coneTipPosition;   // position of the tip of the asymptotic cone
    std::array<FPNumber, 3> apexPosition;   // position of the tip of the asymptotic cone
    FPNumber height = 1.0;

};

#endif  // FDTD_HYPERBOLOID_H_
