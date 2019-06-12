
#include "Hyperboloid.h"

void Hyperboloid::SetConeAngle(const FPNumber angle) {
    coneAngle = angle;
    assert(coneAngle < M_PI);
    SetApexPosition(apexPosition);
}

void Hyperboloid::SetConeAngleInDegrees(const FPNumber angle) {
    coneAngle = angle * (M_PI / 180.0);
    SetApexPosition(apexPosition);
}

void Hyperboloid::SetApexRadius(const FPNumber radius) {
    apexRadius = radius;
    SetApexPosition(apexPosition);
}

void Hyperboloid::SetApexPosition(const std::array<FPNumber, 3> pos) {
    apexPosition = pos;
    FPNumber a, b;
    GetCanonicalScaleFators(a, b);

    coneTipPosition = std::array<FPNumber, 3>{apexPosition[0],
                                              apexPosition[1] + a,
                                              apexPosition[2]
                                              };
}

void Hyperboloid::SetHeight(const FPNumber h) {
    height = h;
}

void Hyperboloid::GetCanonicalScaleFators(FPNumber& a, FPNumber& b) {
    // (y/a)^2 - (x/b)^2 = 1
    // hyperbola:  y = -a/b * sqrt(x^2 + b^2)
    FPNumber b_a = std::atan(coneAngle/2.0);    // b/a
    b = apexRadius / b_a;  // apexRadius = b^2 / a
    a = b / b_a;

    //std::cout << "a: " << a << " b: " << b << std::endl;
}

bool Hyperboloid::IsPointInsideOrOn(std::array<FPNumber, 3> point) {
    // point position relative to cone tip
    std::array<FPNumber, 3> r{point[0] - coneTipPosition[0],
                              point[1] - coneTipPosition[1],
                              point[2] - coneTipPosition[2]
                              };

    FPNumber a, b;
    GetCanonicalScaleFators(a, b);
    FPNumber b_sq = b*b;

    FPNumber rho_sq = r[0]*r[0] + r[2]*r[2];
    FPNumber y = -a/b * std::sqrt(rho_sq + b_sq);

    FPNumber y_base_rel = -a - height;      // relative to the cone tip

    if(r[1] <= y && r[1] >= y_base_rel) {
        return true;
    }
    return false;
}

void Hyperboloid::ArePointsInsideOrOn(
        std::vector<std::array<FPNumber, 3>>& points,
        std::vector<bool>& areInside) {

    FPNumber a, b;
    GetCanonicalScaleFators(a, b);
    FPNumber b_sq = b*b;
    FPNumber a_b = a/b;
    FPNumber y_base_rel = -a - height;      // relative to the cone tip

    for(std::size_t i = 0; i < points.size(); ++i) {
        std::array<FPNumber, 3>& point = points[i];

        std::array<FPNumber, 3> r{point[0] - coneTipPosition[0],
                                  point[1] - coneTipPosition[1],
                                  point[2] - coneTipPosition[2]
                                  };

        FPNumber rho_sq = r[0]*r[0] + r[2]*r[2];
        FPNumber y = -a_b * std::sqrt(rho_sq + b_sq);

        if(r[1] <= y && r[1] >= y_base_rel) {
            areInside[i] = true;
        } else {
            areInside[i] = false;
        }
    }
}

void Hyperboloid::AreGridPointsInsideOrOn(
        const NumberArray3D<FPNumber>& gridArray,
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

    FPNumber x, y, z;

    FPNumber a, b;
    GetCanonicalScaleFators(a, b);
    FPNumber b_sq = b*b;
    FPNumber a_b = a/b;
    FPNumber y_base_rel = -a - height;      // relative to the cone tip

    for(std::size_t i0 = 0; i0 < n0; ++i0) {
        x = x0 + (FPNumber)i0*dx - coneTipPosition[0];  // coordinates relative to cone tip
        for(std::size_t i1 = 0; i1 < n1; ++i1) {
            y = y0 + (FPNumber)i1*dy - coneTipPosition[1];
            for(std::size_t i2 = 0; i2 < n2; ++i2) {
                z = z0 + (FPNumber)i2*dz - coneTipPosition[2];

                FPNumber rho_sq = x*x + z*z;
                FPNumber y_proj = -a_b * std::sqrt(rho_sq + b_sq);  // (x, z) projected on the hyperbola

                if(y <= y_proj && y >= y_base_rel) {
                    areInside[i0][i1][i2] = 1;  // inside
                } else {
                    areInside[i0][i1][i2] = 0;
                }
            }
        }
    }
}


void Hyperboloid::SubdevideSurface2D(FPNumber x_cut,     // sbdeviides the yz cross section at x = x_cut to
        FPNumber maxArcLength,      // maximum length of each section
        std::vector<std::array<FPNumber, 3>>& centerPoints,     // the point at the middle of each patch
        std::vector<std::array<FPNumber, 3>>& normalVecs,   // normal unit vector pointing outwards
        std::vector<FPNumber>& arcLenghts,      // patch area
        std::vector<std::vector<std::array<FPNumber, 3>>>* arcSubdivisionPoints,    // if defined, provides
        std::size_t numSubSubdivisionPoints
        ) {
    assert(false);
}

void Hyperboloid::GetBoundingBox2D(FPNumber x_cut,
        std::array<FPNumber, 3>& lowerLeftCorner,
        std::array<FPNumber, 3>& upperRightCorner
        ) {
    assert(false);
}

void Hyperboloid::SubdevideSurface(
        FPNumber maxElementSurfaceArea,      // with maximum length maxArcLength
        std::vector<std::array<FPNumber, 3>>& centerPoints,     // output: the point at the middle of each patch
        std::vector<std::array<FPNumber, 3>>& normalVecs,   // output: normal unit vector pointing outwards
        std::vector<FPNumber>& elementSurfaceAreas,      // output: patch area
        std::vector<std::vector<std::array<FPNumber, 3>>>* surfaceSubdivisionPoints ,    // if defined, provides
        std::size_t numSubSubdivisionPoints    // this many equaly spaced points on each subsection
        ) {
    assert(false);
}

void Hyperboloid::GetBoundingBox(
        std::array<FPNumber, 3>& lowerLeftCorner,
        std::array<FPNumber, 3>& upperRightCorner
        ) {
    assert(false);
}

