
#include <cassert>
#include <cmath>

#include "SpatialCubeGridArrayManipulator.h"


void SpatialCubeGridArrayManipulator::SetCubeCorners(std::array<RealNumber, 3>& r0, std::array<RealNumber, 3>& r1) {
    assert(r0[0] <= r1[0] && r0[1] <= r1[1] && r0[2] <= r1[2]);
    cubeR0 = r0;
    cubeR1 = r1;
}

void SpatialCubeGridArrayManipulator::SetEdgeThickness(std::array<RealNumber, 3>& thickness) {
    assert(thickness[0] >= 0.0 && thickness[1] >= 0.0 && thickness[2] >= 0.0);
    smoothEdgeThickness = thickness;
}

void SpatialCubeGridArrayManipulator::SetInsideValue(RealNumber value) {
    insideValue = value;
}

void SpatialCubeGridArrayManipulator::SetOutsideValue(RealNumber value) {
    outsideValue = value;
}

std::pair<std::array<std::size_t, 3>, std::array<std::size_t, 3>>
    SpatialCubeGridArrayManipulator::GetArrayIndicesBetweenTwoCorners(std::array<RealNumber, 3>& cornerR0,
                                                                      std::array<RealNumber, 3>& cornerR1) {
    std::array<std::size_t, 3> shape = gridArray.GetShape();
    RealNumber dx = (r1[0] - r0[0])/shape[0];
    RealNumber dy = (r1[1] - r0[1])/shape[1];
    RealNumber dz = (r1[2] - r0[2])/shape[2];

    std::array<std::size_t, 3> indStart = {(std::size_t)std::ceil((cornerR0[0] - r0[0])/dx),
                                           (std::size_t)std::ceil((cornerR0[1] - r0[1])/dy),
                                           (std::size_t)std::ceil((cornerR0[2] - r0[2])/dz)};
    if(cornerR0[0] <= r0[0]) {indStart[0] = 0;}
    if(cornerR0[1] <= r0[1]) {indStart[1] = 0;}
    if(cornerR0[2] <= r0[2]) {indStart[2] = 0;}

    std::array<std::size_t, 3> indEnd = {(std::size_t)std::floor((cornerR1[0] - r0[0])/dx),
                                         (std::size_t)std::floor((cornerR1[1] - r0[1])/dy),
                                         (std::size_t)std::floor((cornerR1[2] - r0[2])/dz)};

    if(indEnd[0] > shape[0]) {indEnd[0] = shape[0];}
    if(indEnd[1] > shape[1]) {indEnd[1] = shape[1];}
    if(indEnd[2] > shape[2]) {indEnd[2] = shape[2];}
    if(cornerR1[0] <= r0[0]) {indEnd[0] = 0;}
    if(cornerR1[1] <= r0[1]) {indEnd[1] = 0;}
    if(cornerR1[2] <= r0[2]) {indEnd[2] = 0;}

    std::pair<std::array<std::size_t, 3>, std::array<std::size_t, 3>> indexPair(indStart, indEnd);
    return indexPair;
}

void SpatialCubeGridArrayManipulator::UpdateArray(const RealNumber t) {
    gridArray.SetToNumber(0.0);

    std::array<std::size_t, 3> shape = gridArray.GetShape();
    RealNumber dx = (r1[0] - r0[0])/shape[0];
    RealNumber dy = (r1[1] - r0[1])/shape[1];
    RealNumber dz = (r1[2] - r0[2])/shape[2];

    // get an slice on the cube and set its value to insideValue
    std::array<RealNumber, 3> cubeR0_p{cubeR0[0] + smoothEdgeThickness[0]/2,
                                       cubeR0[1] + smoothEdgeThickness[1]/2,
                                       cubeR0[2] + smoothEdgeThickness[2]/2};
    std::array<RealNumber, 3> cubeR0_m{cubeR0[0] - smoothEdgeThickness[0]/2,
                                       cubeR0[1] - smoothEdgeThickness[1]/2,
                                       cubeR0[2] - smoothEdgeThickness[2]/2};
    std::array<RealNumber, 3> cubeR1_p{cubeR1[0] + smoothEdgeThickness[0]/2,
                                       cubeR1[1] + smoothEdgeThickness[1]/2,
                                       cubeR1[2] + smoothEdgeThickness[2]/2};
    std::array<RealNumber, 3> cubeR1_m{cubeR1[0] - smoothEdgeThickness[0]/2,
                                       cubeR1[1] - smoothEdgeThickness[1]/2,
                                       cubeR1[2] - smoothEdgeThickness[2]/2};

    std::pair<std::array<std::size_t, 3>, std::array<std::size_t, 3>> indexPair_mp =
            GetArrayIndicesBetweenTwoCorners(cubeR0_m, cubeR1_p);
    std::pair<std::array<std::size_t, 3>, std::array<std::size_t, 3>> indexPair_pm =
            GetArrayIndicesBetweenTwoCorners(cubeR0_p, cubeR1_m);

    std::array<std::size_t, 3>& indStart_p = indexPair_pm.first;
    std::array<std::size_t, 3>& indStart_m = indexPair_mp.first;
    std::array<std::size_t, 3>& indEnd_p = indexPair_mp.second;
    std::array<std::size_t, 3>& indEnd_m = indexPair_pm.second;

    assert(indStart_m[0] <= indEnd_p[0] && indStart_m[1] <= indEnd_p[1] && indStart_m[2] <= indEnd_p[2]);
    assert(indStart_m[0] < indEnd_p[0] || indStart_m[1] < indEnd_p[1] || indStart_m[2] < indEnd_p[2]);

    NumberArray3D slice = gridArray.GetSlice(indStart_p, indEnd_m);
    slice.SetToNumber(1.0);

    if(smoothEdgeThickness[2] > 0.0) {
    // get an slice for the cube face with normal along +z or -z

        std::array<std::size_t, 3> indStart_zm{indStart_m[0], indStart_m[1], indStart_m[2]};
        std::array<std::size_t, 3> indEnd_zm{indEnd_p[0], indEnd_p[1], indStart_p[2]};
        std::array<std::size_t, 3> indStart_zp{indStart_m[0], indStart_m[1], indEnd_m[2]};
        std::array<std::size_t, 3> indEnd_zp{indEnd_p[0], indEnd_p[1], indEnd_p[2]};

        NumberArray3D slice_zm = gridArray.GetSlice(indStart_zm, indEnd_zm);
        slice_zm.SetToNumber(1.0);
        NumberArray3D slice_zp = gridArray.GetSlice(indStart_zp, indEnd_zp);
        slice_zp.SetToNumber(1.0);

        // create meshgrid
        std::array<RealNumber, 3> r0_zm{indStart_zm[0]*dx, indStart_zm[1]*dy, indStart_zm[2]*dz};
        std::array<RealNumber, 3> r1_zm{indEnd_zm[0]*dx, indEnd_zm[1]*dy, indEnd_zm[2]*dz};
        NumberArray3D<RealNumber> zm = NumberArray3D<RealNumber>::GetMeshGrid(slice_zm.GetShape(), r0_zm, r1_zm, 2);
        slice_zm *= 0.5*NumberArray3D<RealNumber>::sin((zm - r0_zm[2])*M_PI/(r1_zm[2] - r0_zm[2]) - M_PI/2) + 0.5;

        std::array<RealNumber, 3> r0_zp{r0[0] + indStart_zp[0]*dx, r0[1] + indStart_zp[1]*dy, r0[2] + indStart_zp[2]*dz};
        std::array<RealNumber, 3> r1_zp{r0[0] + indEnd_zp[0]*dx, r0[1] + indEnd_zp[1]*dy, r0[2] + indEnd_zp[2]*dz};
        NumberArray3D<RealNumber> zp = NumberArray3D<RealNumber>::GetMeshGrid(slice_zp.GetShape(), r0_zp, r1_zp, 2);
        slice_zp *= 0.5*NumberArray3D<RealNumber>::cos((zp - r0_zp[2])*M_PI/(r1_zp[2] - r0_zp[2])) + 0.5;
    }

    gridArray *= (insideValue - outsideValue);
    gridArray += outsideValue;
}

RealNumber SpatialCubeGridArrayManipulator::CalculateTime(const RealNumber dt, const std::size_t timeIndex) {
    return timeIndex*dt;
}


