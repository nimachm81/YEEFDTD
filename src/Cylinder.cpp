
#include "Cylinder.h"


void Cylinder::SetRadius(const FPNumber r) {
    radius = r;
}

void Cylinder::SetTopPlane(const std::array<FPNumber, 3> center) {
    topPlaneCenter = center;
}

void Cylinder::SetHeight(const FPNumber h) {
    height = h;
}

void Cylinder::SetEvenAlignment(bool alignEven) {
    align_even = alignEven;
}

bool Cylinder::IsPointInsideOrOn(std::array<FPNumber, 3> point) {
    FPNumber x_tc = topPlaneCenter[0];
    FPNumber y_tc = topPlaneCenter[1];
    FPNumber z_tc = topPlaneCenter[2];

    FPNumber y_bc = y_tc - height;      // bottom center

    FPNumber x = point[0] - x_tc;
    FPNumber y = point[1] - y_tc;
    FPNumber z = point[2] - z_tc;
    FPNumber rho = std::sqrt(x*x + z*z);

    if(y > 0.0 || y < -height) {
        return false;
    } else {
        if( rho <= radius) {
            return true;
        } else {
            return false;
        }
    }
}



void Cylinder::ArePointsInsideOrOn(std::vector<std::array<FPNumber, 3>>& points,
                                        std::vector<bool>& areInside) {
    if(areInside.size() < points.size()) {
        areInside.resize(points.size());
    }

    FPNumber x_tc = topPlaneCenter[0];  // top center
    FPNumber y_tc = topPlaneCenter[1];
    FPNumber z_tc = topPlaneCenter[2];

    FPNumber r_sq = radius * radius;

    for(std::size_t i = 0; i < points.size(); ++i) {
        auto& point = points[i];
        FPNumber x = point[0] - x_tc;
        FPNumber y = point[1] - y_tc;
        FPNumber z = point[2] - z_tc;
        FPNumber rho_sq = x*x + z*z;

        if(y > 0.0 || y < -height) {
            areInside[i] = false;
        } else {
            if( rho_sq <= r_sq) {
                areInside[i] = true;
            } else {
                areInside[i] = false;
            }
        }
    }
}


void Cylinder::AreGridPointsInsideOrOn(const NumberArray3D<FPNumber>& gridArray,
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

    FPNumber x_tc = topPlaneCenter[0];  // top center
    FPNumber y_tc = topPlaneCenter[1];
    FPNumber z_tc = topPlaneCenter[2];

    FPNumber r_sq = radius * radius;

    FPNumber x, y, z, rho_sq;

    auto shape = gridArray.GetShape();
    shape[1] = 1;
    NumberArray3D<std::int8_t> maskArray(shape, 0);
    std::int8_t*** maskArr = maskArray.GetArrayData();

    if(!align_even) {
        for(std::size_t i0 = 0; i0 < n0; ++i0) {
            x = x0 + (FPNumber)i0*dx - x_tc;
            for(std::size_t i2 = 0; i2 < n2; ++i2) {
                z = z0 + (FPNumber)i2*dz - z_tc;   // z coordinate relative to tip
                rho_sq = x*x + z*z;
                if( rho_sq <= r_sq) {
                    maskArr[i0][0][i2] = 1;
                }
            }
        }
    } else {
        for(std::size_t i0 = 0; i0 < n0; i0 += 2) {
            x = x0 + (FPNumber)i0*dx - x_tc;
            for(std::size_t i2 = 0; i2 < n2; i2 += 2) {
                z = z0 + (FPNumber)i2*dz - z_tc;
                rho_sq = x*x + z*z;
                if( rho_sq <= r_sq) {
                    maskArr[i0][0][i2] = 1;
                }
            }
        }
        /*
        // fill the gaps
        for(std::size_t i0 = 1; i0 < n0 - 1; i0 += 2) {
            for(std::size_t i2 = 1; i2 < n2 - 1; i2 += 2) {
                if( maskArr[i0 - 1][0][i2 - 1] == 1 && maskArr[i0 + 1][0][i2 - 1] == 1 &&
                    maskArr[i0 - 1][0][i2 + 1] == 1 && maskArr[i0 + 1][0][i2 + 1] == 1) {
                    maskArr[i0][0][i2] = 1;
                }
            }
        }
        for(std::size_t i0 = 1; i0 < n0 - 1; i0 += 2) {
            for(std::size_t i2 = 0; i2 < n2; i2 += 2) {
                if( maskArr[i0 - 1][0][i2] == 1 && maskArr[i0 + 1][0][i2] == 1 ) {
                    maskArr[i0][0][i2] = 1;
                }
            }
        }
        for(std::size_t i0 = 0; i0 < n0; i0 += 2) {
            for(std::size_t i2 = 1; i2 < n2 - 1; i2 += 2) {
                if( maskArr[i0][0][i2 - 1] == 1 && maskArr[i0][0][i2 + 1] == 1 ) {
                    maskArr[i0][0][i2] = 1;
                }
            }
        }
        */
    }

    /*
    for(std::size_t i0 = 0; i0 < n0; ++i0) {
        int n_set = 0;
        for(std::size_t i2 = 0; i2 < n2; ++i2) {
            if( maskArr[i0][0][i2] == 1 ) {
                n_set++;
            }
        }
        std::cout << align_even << " : " << n0 << " " << n1 << " " << n2 << " : " <<  n_set << std::endl;
    }*/

    /*
    if(align_even) {
        for(std::size_t i0 = 1; i0 < n0 - 1; ++i0) {
            for(std::size_t i2 = 1; i2 < n2 - 1; ++i2) {
                if(maskArr[i0][0][i2] == 1 && (maskArr[i0 - 1][0][i2] == 0 || maskArr[i0 + 1][0][i2] == 0)) {
                    assert(i0 % 2 == 0);
                }
                if(maskArr[i0][0][i2] == 1 && (maskArr[i0][0][i2 - 1] == 0 || maskArr[i0][0][i2 + 1] == 0)) {
                    assert(i2 % 2 == 0);
                }
            }
        }
    }*/

    for(std::size_t i0 = 0; i0 < n0; ++i0) {
        for(std::size_t i1 = 0; i1 < n1; ++i1) {
            y = y0 + (FPNumber)i1*dy - y_tc;   // y coordinate relative to tip
            for(std::size_t i2 = 0; i2 < n2; ++i2) {
                if(y > 0.0 || y < -height) {
                    areInside[i0][i1][i2] = 0;
                } else {
                    areInside[i0][i1][i2] = maskArr[i0][0][i2];
                }
            }
        }
    }
}


void Cylinder::SubdevideSurface2D(FPNumber x_cut,
                                       FPNumber maxArcLength,
                                       std::vector<std::array<FPNumber, 3>>& centerPoints,
                                       std::vector<std::array<FPNumber, 3>>& normalVecs,
                                       std::vector<FPNumber>& arcLenghts,
                                       std::vector<std::vector<std::array<FPNumber, 3>>>* arcSubdivisionPoints,
                                       std::size_t numSubSubdivisionPoints
                                       ) {
    assert(false);
}


void Cylinder::GetBoundingBox2D(FPNumber x_cut,
                                     std::array<FPNumber, 3>& lowerLeftCorner,
                                     std::array<FPNumber, 3>& upperRightCorner
                                    ) {
    assert(false);
}

void Cylinder::SubdevideSurface(
        FPNumber maxElementSurfaceArea,      // with maximum length maxArcLength
        std::vector<std::array<FPNumber, 3>>& centerPoints,     // output: the point at the middle of each patch
        std::vector<std::array<FPNumber, 3>>& normalVecs,   // output: normal unit vector pointing outwards
        std::vector<FPNumber>& elementSurfaceAreas,      // output: patch area
        std::vector<std::vector<std::array<FPNumber, 3>>>* surfaceSubdivisionPoints,    // if defined, provides
        std::size_t numSubSubdivisionPoints     // this many equaly spaced points on each subsection
        ) {
    assert(false);
}

void Cylinder::GetBoundingBox(
        std::array<FPNumber, 3>& lowerLeftCorner,
        std::array<FPNumber, 3>& upperRightCorner
        ) {
    assert(false);
}

