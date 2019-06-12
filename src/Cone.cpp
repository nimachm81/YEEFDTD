
#include <iostream>
#include <cassert>

#include "Cone.h"

void Cone::SetConeAngle(const FPNumber angle) {
    coneAngle = angle;
    assert(coneAngle > 0.0 && coneAngle < M_PI);
    // tip position and height depend on angle and should be recalculated
    SetApexPosition(apexPosition);
    SetApexToBaseDistance(apexToBaseDistance);
}

void Cone::SetConeAngleInDegrees(const FPNumber angle_degree) {
    coneAngle = angle_degree * M_PI / 180.0;
    assert(coneAngle > 0.0 && coneAngle < M_PI);
    // tip position and height should be recalculated
    SetApexPosition(apexPosition);
    SetApexToBaseDistance(apexToBaseDistance);
}

void Cone::SetTipRadius(const FPNumber radius) {
    apexRadius = radius;
    // tip position and height depend on radius and should be reset
    SetApexPosition(apexPosition);
    SetApexToBaseDistance(apexToBaseDistance);
}

void Cone::SetApexToBaseDistance(const FPNumber distance) {
    apexToBaseDistance = distance;
    assert(apexToBaseDistance > 0.0);

    FPNumber tip_to_flattop = apexRadius * std::cos(coneAngle/2.0) / std::tan(coneAngle/2.0);
    FPNumber tip_to_circleCenter = tip_to_flattop + apexRadius*std::sin(coneAngle/2.0);    // flattop : rounded cap removed
    FPNumber tip_to_roundedtop = tip_to_circleCenter - apexRadius;

    sharpConeHeight = tip_to_roundedtop + apexToBaseDistance;
    assert(sharpConeHeight*std::tan(coneAngle/2.0) > apexRadius*std::cos(coneAngle/2.0)); // base width > cap width
}

void Cone::SetApexPosition(const std::array<FPNumber, 3> pos) {
    apexPosition = pos;

    FPNumber tip_to_flattop = apexRadius * std::cos(coneAngle/2.0) / std::tan(coneAngle/2.0);
    FPNumber tip_to_circleCenter = tip_to_flattop + apexRadius*std::sin(coneAngle/2.0);
    FPNumber tip_to_roundedtop = tip_to_circleCenter - apexRadius;

    tipPosition = std::array<FPNumber, 3>{pos[0],
                                          pos[1] + tip_to_roundedtop,
                                          pos[2]};
}

void Cone::CloseBase(bool close) {
    closeBase = close;
}

bool Cone::IsPointInsideOrOn(std::array<FPNumber, 3> point) {
    FPNumber x_tip = tipPosition[0];
    FPNumber y_tip = tipPosition[1];
    FPNumber z_tip = tipPosition[2];

    FPNumber base_radius = sharpConeHeight * std::tan(coneAngle/2.0);
    FPNumber tip_to_top = apexRadius * std::cos(coneAngle/2.0) / std::tan(coneAngle/2.0);
    FPNumber tip_to_circleCenter = tip_to_top + apexRadius*std::sin(coneAngle/2.0);

    FPNumber x = point[0] - x_tip;
    FPNumber y = point[1] - y_tip;
    FPNumber z = point[2] - z_tip;
    FPNumber rho = std::sqrt(x*x + z*z);

    FPNumber slope_p = std::tan(M_PI/2.0 + coneAngle/2.0);
    FPNumber slope_m = std::tan(M_PI/2.0 - coneAngle/2.0);

    if(y > 0.0 || y < -sharpConeHeight) {
        return false;
    } else if(y <= -tip_to_top) {
        if(std::abs(rho) > base_radius) {
            return false;
        } else {

            if(y <= rho*slope_p && y <= rho*slope_m) {
                return true;
            } else {
                return false;
            }
        }
    } else {    // check whether inside cap circle
        if((y + tip_to_circleCenter)*(y + tip_to_circleCenter) + rho*rho <= apexRadius*apexRadius) {
            return true;
        } else {
            return false;
        }
    }
}


void Cone::ArePointsInsideOrOn(std::vector<std::array<FPNumber, 3>>& points,
                                        std::vector<bool>& areInside) {
    if(areInside.size() < points.size()) {
        areInside.resize(points.size());
    }

    FPNumber x_tip = tipPosition[0];
    FPNumber y_tip = tipPosition[1];
    FPNumber z_tip = tipPosition[2];

    FPNumber base_radius = sharpConeHeight * std::tan(coneAngle/2.0);
    FPNumber tip_to_top = apexRadius * std::cos(coneAngle/2.0) / std::tan(coneAngle/2.0);
    FPNumber tip_to_circleCenter = tip_to_top + apexRadius*std::sin(coneAngle/2.0);

    FPNumber slope_p = std::tan(M_PI/2.0 + coneAngle/2.0);
    FPNumber slope_m = std::tan(M_PI/2.0 - coneAngle/2.0);

    for(std::size_t i = 0; i < points.size(); ++i) {
        auto& point = points[i];
        FPNumber x = point[0] - x_tip;
        FPNumber y = point[1] - y_tip;
        FPNumber z = point[2] - z_tip;
        FPNumber rho = std::sqrt(x*x + z*z);

        if(y > 0.0 || y < -sharpConeHeight) {
            areInside[i] = false;
            continue;
        } else if(y <= -tip_to_top) {
            if(std::abs(rho) > base_radius) {
                areInside[i] = false;
                continue;
            } else {
                if(y <= rho*slope_p && y <= rho*slope_m) {
                    areInside[i] = true;
                    continue;
                } else {
                    areInside[i] = false;
                    continue;
                }
            }
        } else {    // check whether inside cap circle
            if((y + tip_to_circleCenter)*(y + tip_to_circleCenter) + rho*rho <= apexRadius*apexRadius) {
                areInside[i] = true;
                continue;
            } else {
                areInside[i] = false;
                continue;
            }
        }
    }
}


void Cone::AreGridPointsInsideOrOn(const NumberArray3D<FPNumber>& gridArray,
                                    const std::array<FPNumber, 3>& r0,
                                    const std::array<FPNumber, 3>& dr,
                                    std::int8_t*** areInside
                                    ) {
    FPNumber*** arrayData = gridArray.GetArrayData();
    const std::array<std::size_t, 3>& arrayShape = gridArray.GetShape();
    const std::array<std::size_t, 3>& arrayIndStart = gridArray.GetIndStart();

    std::size_t n0 = arrayShape[0];
    std::size_t n1 = arrayShape[1];
    std::size_t n2 = arrayShape[2];
    std::size_t ind0 = arrayIndStart[0];
    std::size_t ind1 = arrayIndStart[1];
    std::size_t ind2 = arrayIndStart[2];

    FPNumber x0 = r0[0];
    FPNumber y0 = r0[1];
    FPNumber z0 = r0[2];

    FPNumber dx = dr[0];
    FPNumber dy = dr[1];
    FPNumber dz = dr[2];

    FPNumber x_tip = tipPosition[0];
    FPNumber y_tip = tipPosition[1];
    FPNumber z_tip = tipPosition[2];

    FPNumber base_radius = sharpConeHeight * std::tan(coneAngle/2.0);
    FPNumber tip_to_top = apexRadius * std::cos(coneAngle/2.0) / std::tan(coneAngle/2.0);
    FPNumber tip_to_circleCenter = tip_to_top + apexRadius*std::sin(coneAngle/2.0);

    FPNumber x, y, z, rho, rho_sq;
    FPNumber apexRadius_sq = apexRadius*apexRadius;
    FPNumber base_radius_sq = base_radius*base_radius;

    FPNumber slope_p = std::tan(M_PI/2.0 + coneAngle/2.0);
    FPNumber slope_m = std::tan(M_PI/2.0 - coneAngle/2.0);

    for(std::size_t i0 = 0; i0 < n0; ++i0) {
        x = x0 + (FPNumber)i0*dx - x_tip;
        for(std::size_t i1 = 0; i1 < n1; ++i1) {
            y = y0 + (FPNumber)i1*dy - y_tip;   // y coordinate relative to tip
            for(std::size_t i2 = 0; i2 < n2; ++i2) {
                if(y > 0.0 || y < -sharpConeHeight) {
                    areInside[i0][i1][i2] = 0;
                } else if(y <= -tip_to_top) {
                    z = z0 + (FPNumber)i2*dz - z_tip;   // z coordinate relative to tip
                    rho_sq = x*x + z*z;
                    if(rho_sq > base_radius_sq) {
                        areInside[i0][i1][i2] = 0;
                    } else {
                        rho = std::sqrt(rho_sq);
                        if(y <= rho*slope_p && y <= rho*slope_m) {
                            areInside[i0][i1][i2] = 1;
                        } else {
                            areInside[i0][i1][i2] = 0;
                        }
                    }
                } else {    // check whether inside cap circle
                    z = z0 + (FPNumber)i2*dz - z_tip;   // z coordinate relative to tip
                    rho_sq = x*x + z*z;
                    if((y + tip_to_circleCenter)*(y + tip_to_circleCenter) + rho_sq <= apexRadius_sq) {
                        areInside[i0][i1][i2] = 1;
                    } else {
                        areInside[i0][i1][i2] = 0;
                    }
                }
            }
        }
    }
}


void Cone::SubdevideSurface2D(FPNumber x_cut,
                                       FPNumber maxArcLength,
                                       std::vector<std::array<FPNumber, 3>>& centerPoints,
                                       std::vector<std::array<FPNumber, 3>>& normalVecs,
                                       std::vector<FPNumber>& arcLenghts,
                                       std::vector<std::vector<std::array<FPNumber, 3>>>* arcSubdivisionPoints,
                                       std::size_t numSubSubdivisionPoints
                                       ) {
    assert(false);
}


void Cone::GetBoundingBox2D(FPNumber x_cut,
                                     std::array<FPNumber, 3>& lowerLeftCorner,
                                     std::array<FPNumber, 3>& upperRightCorner
                                    ) {
    assert(false);
}

void Cone::SubdevideSurface(
        FPNumber maxElementSurfaceArea,      // with maximum length maxArcLength
        std::vector<std::array<FPNumber, 3>>& centerPoints,     // output: the point at the middle of each patch
        std::vector<std::array<FPNumber, 3>>& normalVecs,   // output: normal unit vector pointing outwards
        std::vector<FPNumber>& elementSurfaceAreas,      // output: patch area
        std::vector<std::vector<std::array<FPNumber, 3>>>* surfaceSubdivisionPoints,    // if defined, provides
        std::size_t numSubSubdivisionPoints     // this many equaly spaced points on each subsection
        ) {
    assert(false);
}

void Cone::GetBoundingBox(
        std::array<FPNumber, 3>& lowerLeftCorner,
        std::array<FPNumber, 3>& upperRightCorner
        ) {
    assert(false);
}


