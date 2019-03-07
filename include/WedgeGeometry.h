#ifndef FDTD_WEDGEGEOMETRY_H_
#define FDTD_WEDGEGEOMETRY_H_

#include <vector>
#include <array>
#include <cstddef>

#include "NumberTypes.h"
#include "Geometry.h"

class WedgeGeometry : public Geometry{
    public:
    virtual ~WedgeGeometry() { };
    void SetWedgeAngle(const FPNumber angle);
    void SetWedgeAngleInDegrees(const FPNumber angle_degree);
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
    private:
    int uniformAxis = 0;    // along this axis there is no varation (0:x, 1:y, 2:z)
    int wedgeDirection = 1;      // the wedge is pointing along this direction (0:x, 1:y, 2:z)
    FPNumber wedgeAngle = M_PI/4.0;         // wedge angle in radians
    FPNumber apexRadius = 0.0;   // for tipRadius > 0 the tip is rounded
    std::array<FPNumber, 3> tipPosition;    // position of the unrounded tip
    std::array<FPNumber, 3> apexPosition;    // position of the rounded tip (apex)
    FPNumber apexToBaseDistance = 1.0;
    FPNumber wedgeHeight = 1.0;         // distance from unrounded tip to base of the wedge
    bool closeBase = true;

};

#endif // FDTD_WEDGEGEOMETRY_H_
