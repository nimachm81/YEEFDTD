
#ifndef FDTD_CONE_H_
#define FDTD_CONE_H_

#include <vector>
#include <array>
#include <cstddef>

#include "NumberTypes.h"
#include "Geometry.h"

class Cone : public Geometry {
    public:
    virtual ~Cone() { };
    void SetConeAngle(const FPNumber angle);
    void SetConeAngleInDegrees(const FPNumber angle_degree);
    void SetTipRadius(const FPNumber radius);
    void SetApexToBaseDistance(const FPNumber apexToBaseDistance);
    void SetApexPosition(const std::array<FPNumber, 3> pos);
    void CloseBase(bool close);

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
    FPNumber coneAngle = M_PI/4.0;         // wedge angle in radians (full angle)
    FPNumber apexRadius = 0.0;   // for tipRadius > 0 the tip is rounded
    std::array<FPNumber, 3> tipPosition;    // position of the unrounded sharp tip
    std::array<FPNumber, 3> apexPosition;    // position of the rounded tip (apex)
    FPNumber apexToBaseDistance = 1.0;
    FPNumber sharpConeHeight = 1.0;         // distance from unrounded tip to base of the cone
    bool closeBase = true;

};

#endif // FDTD_CONE_H_
