
#include <iostream>
#include <cassert>

#include "WedgeGeometry.h"

void WedgeGeometry::SetWedgeAngle(const FPNumber angle) {
    wedgeAngle = angle;
    assert(wedgeAngle > 0.0 && wedgeAngle < M_PI);
    // tip position and height depend on angle and should be recalculated
    SetApexPosition(apexPosition);
    SetApexToBaseDistance(apexToBaseDistance);
}

void WedgeGeometry::SetWedgeAngleInDegrees(const FPNumber angle_degree) {
    wedgeAngle = angle_degree * M_PI / 180.0;
    assert(wedgeAngle > 0.0 && wedgeAngle < M_PI);
    // tip position and height should be recalculated
    SetApexPosition(apexPosition);
    SetApexToBaseDistance(apexToBaseDistance);
}

void WedgeGeometry::SetTipRadius(const FPNumber radius) {
    apexRadius = radius;
    // tip position and height depend on radius and should be reset
    SetApexPosition(apexPosition);
    SetApexToBaseDistance(apexToBaseDistance);
}

void WedgeGeometry::SetApexToBaseDistance(const FPNumber distance) {
    apexToBaseDistance = distance;
    assert(apexToBaseDistance > 0.0);

    FPNumber tip_to_flattop = apexRadius * std::cos(wedgeAngle/2.0) / std::tan(wedgeAngle/2.0);
    FPNumber tip_to_circleCenter = tip_to_flattop + apexRadius*std::sin(wedgeAngle/2.0);    // flattop : rounded cap removed
    FPNumber tip_to_roundedtop = tip_to_circleCenter - apexRadius;

    wedgeHeight = tip_to_roundedtop + apexToBaseDistance;
    assert(wedgeHeight*std::tan(wedgeAngle/2.0) > apexRadius*std::cos(wedgeAngle/2.0)); // base width > cap width
}

void WedgeGeometry::SetApexPosition(const std::array<FPNumber, 3> pos) {
    apexPosition = pos;

    FPNumber tip_to_flattop = apexRadius * std::cos(wedgeAngle/2.0) / std::tan(wedgeAngle/2.0);
    FPNumber tip_to_circleCenter = tip_to_flattop + apexRadius*std::sin(wedgeAngle/2.0);
    FPNumber tip_to_roundedtop = tip_to_circleCenter - apexRadius;

    tipPosition = std::array<FPNumber, 3>{pos[0],
                                          pos[1] + tip_to_roundedtop,
                                          pos[2]};
}

void WedgeGeometry::CloseBase(bool close) {
    closeBase = close;
}

bool WedgeGeometry::IsPointInsideOrOn(std::array<FPNumber, 3> point) {
    FPNumber x_tip = tipPosition[0];
    FPNumber y_tip = tipPosition[1];
    FPNumber z_tip = tipPosition[2];

    FPNumber base_halfWidth = wedgeHeight * std::tan(wedgeAngle/2.0);
    FPNumber tip_to_top = apexRadius * std::cos(wedgeAngle/2.0) / std::tan(wedgeAngle/2.0);
    FPNumber tip_to_circleCenter = tip_to_top + apexRadius*std::sin(wedgeAngle/2.0);

    FPNumber x = point[0] - x_tip;
    FPNumber y = point[1] - y_tip;
    FPNumber z = point[2] - z_tip;

    FPNumber slope_p = std::tan(M_PI/2.0 + wedgeAngle/2.0);
    FPNumber slope_m = std::tan(M_PI/2.0 - wedgeAngle/2.0);

    if(uniformAxis == 0 && wedgeDirection == 1) {
        if(y > 0.0 || y < -wedgeHeight) {
            return false;
        } else if(y <= -tip_to_top) {
            if(std::abs(z) > base_halfWidth) {
                return false;
            } else {

                if(y <= z*slope_p && y <= z*slope_m) {
                    return true;
                } else {
                    return false;
                }
            }
        } else {    // check whether inside cap circle
            if((y + tip_to_circleCenter)*(y + tip_to_circleCenter) + z*z <= apexRadius*apexRadius) {
                return true;
            } else {
                return false;
            }
        }
    } else {
        std::cout << "error: other wedge directions not implemented." << std::endl;
        assert(false);
    }
}

void WedgeGeometry::ArePointsInsideOrOn(std::vector<std::array<FPNumber, 3>>& points,
                                        std::vector<bool>& areInside) {
    if(areInside.size() < points.size()) {
        areInside.resize(points.size());
    }

    FPNumber x_tip = tipPosition[0];
    FPNumber y_tip = tipPosition[1];
    FPNumber z_tip = tipPosition[2];

    FPNumber base_halfWidth = wedgeHeight * std::tan(wedgeAngle/2.0);
    FPNumber tip_to_top = apexRadius * std::cos(wedgeAngle/2.0) / std::tan(wedgeAngle/2.0);
    FPNumber tip_to_circleCenter = tip_to_top + apexRadius*std::sin(wedgeAngle/2.0);

    FPNumber slope_p = std::tan(M_PI/2.0 + wedgeAngle/2.0);
    FPNumber slope_m = std::tan(M_PI/2.0 - wedgeAngle/2.0);

    for(std::size_t i = 0; i < points.size(); ++i) {
        auto& point = points[i];
        FPNumber x = point[0] - x_tip;
        FPNumber y = point[1] - y_tip;
        FPNumber z = point[2] - z_tip;

        if(uniformAxis == 0 && wedgeDirection == 1) {
            if(y > 0.0 || y < -wedgeHeight) {
                areInside[i] = false;
                continue;
            } else if(y <= -tip_to_top) {
                if(std::abs(z) > base_halfWidth) {
                    areInside[i] = false;
                    continue;
                } else {
                    if(y <= z*slope_p && y <= z*slope_m) {
                        areInside[i] = true;
                        continue;
                    } else {
                        areInside[i] = false;
                        continue;
                    }
                }
            } else {    // check whether inside cap circle
                if((y + tip_to_circleCenter)*(y + tip_to_circleCenter) + z*z <= apexRadius*apexRadius) {
                    areInside[i] = true;
                    continue;
                } else {
                    areInside[i] = false;
                    continue;
                }
            }
        } else {
            std::cout << "error: other wedge directions not implemented." << std::endl;
            assert(false);
        }
    }
}

void WedgeGeometry::AreGridPointsInsideOrOn(const NumberArray3D<FPNumber>& gridArray,
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

    FPNumber base_halfWidth = wedgeHeight * std::tan(wedgeAngle/2.0);
    FPNumber tip_to_top = apexRadius * std::cos(wedgeAngle/2.0) / std::tan(wedgeAngle/2.0);
    FPNumber tip_to_circleCenter = tip_to_top + apexRadius*std::sin(wedgeAngle/2.0);

    FPNumber x, y, z;

    FPNumber slope_p = std::tan(M_PI/2.0 + wedgeAngle/2.0);
    FPNumber slope_m = std::tan(M_PI/2.0 - wedgeAngle/2.0);

    if(uniformAxis == 0 && wedgeDirection == 1) {
        for(std::size_t i0 = 0; i0 < n0; ++i0) {
            //x = x0 + (FPNumber)i0*dx;
            for(std::size_t i1 = 0; i1 < n1; ++i1) {
                y = y0 + (FPNumber)i1*dy - y_tip;   // y coordinate relative to tip
                for(std::size_t i2 = 0; i2 < n2; ++i2) {
                    if(y > 0.0 || y < -wedgeHeight) {
                        areInside[i0][i1][i2] = 0;
                    } else if(y <= -tip_to_top) {
                        z = z0 + (FPNumber)i2*dz - z_tip;   // z coordinate relative to tip
                        if(std::abs(z) > base_halfWidth) {
                            areInside[i0][i1][i2] = 0;
                        } else {
                            if(y <= z*slope_p && y <= z*slope_m) {
                                areInside[i0][i1][i2] = 1;
                            } else {
                                areInside[i0][i1][i2] = 0;
                            }
                        }
                    } else {    // check whether inside cap circle
                        z = z0 + (FPNumber)i2*dz - z_tip;   // z coordinate relative to tip
                        if((y + tip_to_circleCenter)*(y + tip_to_circleCenter) + z*z <= apexRadius*apexRadius) {
                            areInside[i0][i1][i2] = 1;
                        } else {
                            areInside[i0][i1][i2] = 0;
                        }
                    }
                }
            }
        }
    } else {
        std::cout << "error: other wedge directions not implemented." << std::endl;
        assert(false);
    }
}

void WedgeGeometry::SubdevideSurface2D(FPNumber x_cut,
                                       FPNumber maxArcLength,
                                       std::vector<std::array<FPNumber, 3>>& centerPoints,
                                       std::vector<std::array<FPNumber, 3>>& normalVecs,
                                       std::vector<FPNumber>& arcLenghts,
                                       std::vector<std::vector<std::array<FPNumber, 3>>>* arcSubdivisionPoints,
                                       std::size_t numSubSubdivisionPoints
                                       ) {

    const FPNumber sin_t2 = std::sin(wedgeAngle/2.0);
    const FPNumber cos_t2 = std::cos(wedgeAngle/2.0);
    const FPNumber tan_t2 = sin_t2/cos_t2;


    const FPNumber x_tip = tipPosition[0];
    const FPNumber y_tip = tipPosition[1];
    const FPNumber z_tip = tipPosition[2];

    const FPNumber base_halfWidth = wedgeHeight * tan_t2;
    const FPNumber tip_to_top = apexRadius * cos_t2 / tan_t2;
    const FPNumber tip_to_circleCenter = tip_to_top + apexRadius*sin_t2;

    FPNumber slope_p = std::tan(M_PI/2.0 + wedgeAngle/2.0);
    FPNumber slope_m = std::tan(M_PI/2.0 - wedgeAngle/2.0);

    FPNumber side_length = (wedgeHeight - tip_to_top) / cos_t2;     // length of the sides excluding the cap
    FPNumber capArc_length = apexRadius*(M_PI - wedgeAngle);     // size of the cap arc (apex)
    FPNumber base_length = 2*base_halfWidth;

    std::size_t n_cap = (std::size_t)(capArc_length / maxArcLength) + 1;
    if(apexRadius == 0.0) {
        n_cap = 0;
    }
    std::size_t n_side = (std::size_t)(side_length / maxArcLength) + 1;
    if(side_length == 0.0) {
        n_side = 0;
    }
    std::size_t n_base = (std::size_t)(base_length / maxArcLength) + 1;

    std::size_t n_total = 2*n_side + n_cap;
    if(n_base > 0 && closeBase) {
        n_total += n_base;
    }

    centerPoints.resize(n_total);
    normalVecs.resize(n_total);
    arcLenghts.resize(n_total);
    if(arcSubdivisionPoints != nullptr) {
        assert(numSubSubdivisionPoints > 0);
        arcSubdivisionPoints->resize(n_total);
        for(std::size_t i = 0; i < n_total; ++i){
            arcSubdivisionPoints->operator[](i).resize(numSubSubdivisionPoints);
        }
    }

    FPNumber y, z;
    const FPNumber x = x_cut;
    FPNumber vy, vz;
    const FPNumber vx = 0.0;

    std::size_t ind_start = 0;

    if(n_cap > 0) {
        FPNumber y_circleCenter = y_tip - tip_to_circleCenter;
        FPNumber z_circleCenter = z_tip;
        FPNumber theta_0 = wedgeAngle/2.0;
        FPNumber d_theta = (M_PI - wedgeAngle)/n_cap;
        FPNumber theta_i;
        FPNumber th0_i, th1_i, theta_ij, d_th_ij, y_ij, z_ij;    // subsubdivision variables (each subdevided arc element is subdevided again to numSubSubdivisionPoints
        FPNumber dA = apexRadius*d_theta;   // arc length

        for(std::size_t i = 0; i < n_cap; ++i) {
            theta_i = theta_0 + ((FPNumber)i + 0.5)*d_theta;
            y = y_circleCenter + apexRadius*std::sin(theta_i);
            z = z_circleCenter + apexRadius*std::cos(theta_i);

            vy = std::sin(theta_i);
            vz = std::cos(theta_i);

            centerPoints[i + ind_start] = std::array<FPNumber, 3>{x, y, z};
            normalVecs[i + ind_start] = std::array<FPNumber, 3>{vx, vy, vz};
            arcLenghts[i + ind_start] = dA;

            if(arcSubdivisionPoints != nullptr) {
                auto& arcSubdivisionPoints_i = arcSubdivisionPoints->operator[](i + ind_start);

                th0_i = theta_0 + i*d_theta;
                d_th_ij = d_theta / numSubSubdivisionPoints;

                for(std::size_t j = 0; j < numSubSubdivisionPoints; ++j) {
                    theta_ij = th0_i + ((FPNumber)j + 0.5)*d_th_ij;
                    y_ij = y_circleCenter + apexRadius*std::sin(theta_ij);
                    z_ij = z_circleCenter + apexRadius*std::cos(theta_ij);
                    arcSubdivisionPoints_i[j] = std::array<FPNumber, 3>{x, y_ij, z_ij};
                }
            }
        }
        ind_start += n_cap;
    }

    if(n_side > 0) {
        FPNumber y0_side = y_tip - tip_to_top;      // y coordinate of the top point of the sides
        FPNumber z0_side_p = z_tip + tip_to_top*tan_t2;   // z coordinate of top point of the right side
        FPNumber z0_side_m = z_tip - tip_to_top*tan_t2;   // z coordinate of top point of the right side
        FPNumber y1_side = y_tip - wedgeHeight;
        FPNumber z1_side_p = z_tip + base_halfWidth;
        FPNumber z1_side_m = z_tip - base_halfWidth;

        FPNumber dy = (y1_side - y0_side)/n_side;
        FPNumber dz_p = (z1_side_p - z0_side_p)/n_side;
        FPNumber dz_m = (z1_side_m - z0_side_m)/n_side;
        FPNumber dA = side_length/n_side;

        vy = std::sin(wedgeAngle/2.0);
        FPNumber vz_p = std::cos(wedgeAngle/2.0);
        FPNumber vz_m = -std::cos(wedgeAngle/2.0);

        FPNumber z_p, z_m;
        FPNumber y_ij, zp_ij, zm_ij;

        for(std::size_t i = 0; i < n_side; ++i) {
            y = y0_side + ((FPNumber)i + 0.5)*dy;
            z_p = z0_side_p + ((FPNumber)i + 0.5)*dz_p;
            z_m = z0_side_m + ((FPNumber)i + 0.5)*dz_m;

            centerPoints[2*i + ind_start] = std::array<FPNumber, 3>{x, y, z_p};
            normalVecs[2*i + ind_start] = std::array<FPNumber, 3>{vx, vy, vz_p};
            arcLenghts[2*i + ind_start] = dA;
            centerPoints[2*i + 1 + ind_start] = std::array<FPNumber, 3>{x, y, z_m};
            normalVecs[2*i + 1 + ind_start] = std::array<FPNumber, 3>{vx, vy, vz_m};
            arcLenghts[2*i + 1 + ind_start] = dA;

            if(arcSubdivisionPoints != nullptr) {
                auto& arcSubdivisionPoints_2i = arcSubdivisionPoints->operator[](2*i + ind_start);
                auto& arcSubdivisionPoints_2ip1 = arcSubdivisionPoints->operator[](2*i + 1 + ind_start);

                FPNumber y0_i = y0_side + (FPNumber)i*dy;
                FPNumber dy_ij = dy / numSubSubdivisionPoints;

                FPNumber z0_p_i = z0_side_p + (FPNumber)i*dz_p;
                FPNumber z0_m_i = z0_side_m + (FPNumber)i*dz_m;
                FPNumber dz_p_ij = dz_p / numSubSubdivisionPoints;
                FPNumber dz_m_ij = dz_m / numSubSubdivisionPoints;

                for(std::size_t j = 0; j < numSubSubdivisionPoints; ++j) {
                    y_ij = y0_i + ((FPNumber)j + 0.5)*dy_ij;
                    zp_ij = z0_p_i + ((FPNumber)j + 0.5)*dz_p_ij;
                    zm_ij = z0_m_i + ((FPNumber)j + 0.5)*dz_m_ij;
                    arcSubdivisionPoints_2i[j] = std::array<FPNumber, 3>{x, y_ij, zp_ij};
                    arcSubdivisionPoints_2ip1[j] = std::array<FPNumber, 3>{x, y_ij, zm_ij};
                }
            }

        }
        ind_start += 2*n_side;
    }

    if(n_base > 0 && closeBase) {
        FPNumber y0_base = y_tip - wedgeHeight;      // y coordinate of the lower left point of the base
        FPNumber z0_base = z_tip - base_halfWidth;   // z coordinate of the lower left point of the base

        FPNumber dz = base_length/n_base;
        y = y0_base;

        vy = -1.0;
        vz = 0.0;

        FPNumber dA = dz;

        FPNumber z_ij;

        for(std::size_t i = 0; i < n_base; ++i) {
            z = z0_base + ((FPNumber)i + 0.5)*dz;

            centerPoints[i + ind_start] = std::array<FPNumber, 3>{x, y, z};
            normalVecs[i + ind_start] = std::array<FPNumber, 3>{vx, vy, vz};
            arcLenghts[i + ind_start] = dA;

            if(arcSubdivisionPoints != nullptr) {
                auto& arcSubdivisionPoints_i = arcSubdivisionPoints->operator[](i + ind_start);

                FPNumber z0_i = z0_base + (FPNumber)i*dz;
                FPNumber dz_ij = dz / numSubSubdivisionPoints;
                for(std::size_t j = 0; j < numSubSubdivisionPoints; ++j) {
                    z_ij = z0_i + ((FPNumber)j + 0.5)*dz_ij;
                    arcSubdivisionPoints_i[j] = std::array<FPNumber, 3>{x, y, z_ij};
                }
            }
        }
        ind_start += n_base;
    }

}


void WedgeGeometry::GetBoundingBox2D(FPNumber x_cut,
                                     std::array<FPNumber, 3>& lowerLeftCorner,
                                     std::array<FPNumber, 3>& upperRightCorner
                                    ) {
    FPNumber base_halfWidth = wedgeHeight * std::tan(wedgeAngle/2.0);
    FPNumber y0 = apexPosition[1] - apexToBaseDistance;
    FPNumber y1 = apexPosition[1];
    FPNumber z0 = apexPosition[2] - base_halfWidth;
    FPNumber z1 = apexPosition[2] + base_halfWidth;

    lowerLeftCorner[0] = x_cut;
    lowerLeftCorner[1] = y0;
    lowerLeftCorner[2] = z0;
    upperRightCorner[0] = x_cut;
    upperRightCorner[1] = y1;
    upperRightCorner[2] = z1;
}

void WedgeGeometry::SubdevideSurface(
        FPNumber maxElementSurfaceArea,      // with maximum length maxArcLength
        std::vector<std::array<FPNumber, 3>>& centerPoints,     // output: the point at the middle of each patch
        std::vector<std::array<FPNumber, 3>>& normalVecs,   // output: normal unit vector pointing outwards
        std::vector<FPNumber>& elementSurfaceAreas,      // output: patch area
        std::vector<std::vector<std::array<FPNumber, 3>>>* surfaceSubdivisionPoints,    // if defined, provides
        std::size_t numSubSubdivisionPoints     // this many equaly spaced points on each subsection
        ) {
    assert(false);
}

void WedgeGeometry::GetBoundingBox(
        std::array<FPNumber, 3>& lowerLeftCorner,
        std::array<FPNumber, 3>& upperRightCorner
        ) {
    assert(false);
}
