
#include <cassert>
#include <cmath>

#include "SpatialCubeGridArrayManipulator.h"


void SpatialCubeGridArrayManipulator::SetCubeCorners(std::array<RealNumber, 3>& cube_r0, std::array<RealNumber, 3>& cube_r1) {
    assert(std::real(cube_r0[0]) <= std::real(cube_r1[0]) &&
           std::real(cube_r0[1]) <= std::real(cube_r1[1]) &&
           std::real(cube_r0[2]) <= std::real(cube_r1[2]));
    cubeR0 = cube_r0;
    cubeR1 = cube_r1;
}

void SpatialCubeGridArrayManipulator::SetEdgeThickness(std::array<RealNumber, 3>& thickness) {
    assert(std::real(thickness[0]) >= 0.0 &&
           std::real(thickness[1]) >= 0.0 &&
           std::real(thickness[2]) >= 0.0);
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
    RealNumber dx = dr[0];
    RealNumber dy = dr[1];
    RealNumber dz = dr[2];

    std::array<std::size_t, 3> indStart = {(std::size_t)std::ceil(std::real((cornerR0[0] - r0[0])/dx)),
                                           (std::size_t)std::ceil(std::real((cornerR0[1] - r0[1])/dy)),
                                           (std::size_t)std::ceil(std::real((cornerR0[2] - r0[2])/dz))};
    if(std::real(cornerR0[0]) <= std::real(r0[0])) {indStart[0] = 0;}
    if(std::real(cornerR0[1]) <= std::real(r0[1])) {indStart[1] = 0;}
    if(std::real(cornerR0[2]) <= std::real(r0[2])) {indStart[2] = 0;}

    std::array<std::size_t, 3> indEnd = {(std::size_t)std::floor(std::real((cornerR1[0] - r0[0])/dx)),
                                         (std::size_t)std::floor(std::real((cornerR1[1] - r0[1])/dy)),
                                         (std::size_t)std::floor(std::real((cornerR1[2] - r0[2])/dz))};

    if(indEnd[0] > shape[0]) {indEnd[0] = shape[0];}
    if(indEnd[1] > shape[1]) {indEnd[1] = shape[1];}
    if(indEnd[2] > shape[2]) {indEnd[2] = shape[2];}
    if(std::real(cornerR1[0]) <= std::real(r0[0])) {indEnd[0] = 0;}
    if(std::real(cornerR1[1]) <= std::real(r0[1])) {indEnd[1] = 0;}
    if(std::real(cornerR1[2]) <= std::real(r0[2])) {indEnd[2] = 0;}

    std::pair<std::array<std::size_t, 3>, std::array<std::size_t, 3>> indexPair(indStart, indEnd);
    return indexPair;
}

void SpatialCubeGridArrayManipulator::UpdateArray(const RealNumber t) {
    gridArray.SetToNumber(0.0);

    std::array<std::size_t, 3> shape = gridArray.GetShape();
    RealNumber dx = dr[0];
    RealNumber dy = dr[1];
    RealNumber dz = dr[2];

    // get an slice on the cube and set its value to insideValue
    std::array<RealNumber, 3> cubeR0_p{cubeR0[0] + smoothEdgeThickness[0]/(RealNumber)2,
                                       cubeR0[1] + smoothEdgeThickness[1]/(RealNumber)2,
                                       cubeR0[2] + smoothEdgeThickness[2]/(RealNumber)2};
    std::array<RealNumber, 3> cubeR0_m{cubeR0[0] - smoothEdgeThickness[0]/(RealNumber)2,
                                       cubeR0[1] - smoothEdgeThickness[1]/(RealNumber)2,
                                       cubeR0[2] - smoothEdgeThickness[2]/(RealNumber)2};
    std::array<RealNumber, 3> cubeR1_p{cubeR1[0] + smoothEdgeThickness[0]/(RealNumber)2,
                                       cubeR1[1] + smoothEdgeThickness[1]/(RealNumber)2,
                                       cubeR1[2] + smoothEdgeThickness[2]/(RealNumber)2};
    std::array<RealNumber, 3> cubeR1_m{cubeR1[0] - smoothEdgeThickness[0]/(RealNumber)2,
                                       cubeR1[1] - smoothEdgeThickness[1]/(RealNumber)2,
                                       cubeR1[2] - smoothEdgeThickness[2]/(RealNumber)2};

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

    // set the edge frames to 1.0
    if(std::real(smoothEdgeThickness[0]) > 0.0) {
    // get an slice for the cube face with normal along +x or -x

        std::array<std::size_t, 3> indStart_xm{indStart_m[0], indStart_m[1], indStart_m[2]};
        std::array<std::size_t, 3> indEnd_xm{indStart_p[0], indEnd_p[1], indEnd_p[2]};
        std::array<std::size_t, 3> indStart_xp{indEnd_m[0], indStart_m[1], indStart_m[2]};
        std::array<std::size_t, 3> indEnd_xp{indEnd_p[0], indEnd_p[1], indEnd_p[2]};

        NumberArray3D slice_xm = gridArray.GetSlice(indStart_xm, indEnd_xm);
        slice_xm.SetToNumber(1.0);
        NumberArray3D slice_xp = gridArray.GetSlice(indStart_xp, indEnd_xp);
        slice_xp.SetToNumber(1.0);
    }
    if(std::real(smoothEdgeThickness[1]) > 0.0) {
    // get an slice for the cube face with normal along +y or -y

        std::array<std::size_t, 3> indStart_ym{indStart_m[0], indStart_m[1], indStart_m[2]};
        std::array<std::size_t, 3> indEnd_ym{indEnd_p[0], indStart_p[1], indEnd_p[2]};
        std::array<std::size_t, 3> indStart_yp{indStart_m[0], indEnd_m[1], indStart_m[2]};
        std::array<std::size_t, 3> indEnd_yp{indEnd_p[0], indEnd_p[1], indEnd_p[2]};

        NumberArray3D slice_ym = gridArray.GetSlice(indStart_ym, indEnd_ym);
        slice_ym.SetToNumber(1.0);
        NumberArray3D slice_yp = gridArray.GetSlice(indStart_yp, indEnd_yp);
        slice_yp.SetToNumber(1.0);
    }
    if(std::real(smoothEdgeThickness[2]) > 0.0) {
    // get an slice for the cube face with normal along +z or -z

        std::array<std::size_t, 3> indStart_zm{indStart_m[0], indStart_m[1], indStart_m[2]};
        std::array<std::size_t, 3> indEnd_zm{indEnd_p[0], indEnd_p[1], indStart_p[2]};
        std::array<std::size_t, 3> indStart_zp{indStart_m[0], indStart_m[1], indEnd_m[2]};
        std::array<std::size_t, 3> indEnd_zp{indEnd_p[0], indEnd_p[1], indEnd_p[2]};

        NumberArray3D slice_zm = gridArray.GetSlice(indStart_zm, indEnd_zm);
        slice_zm.SetToNumber(1.0);
        NumberArray3D slice_zp = gridArray.GetSlice(indStart_zp, indEnd_zp);
        slice_zp.SetToNumber(1.0);
    }

    // smoothen the edges
    if(std::real(smoothEdgeThickness[0]) > 0.0) {
    // get an slice for the cube face with normal along +x or -x

        std::array<std::size_t, 3> indStart_xm{indStart_m[0], indStart_m[1], indStart_m[2]};
        std::array<std::size_t, 3> indEnd_xm{indStart_p[0], indEnd_p[1], indEnd_p[2]};
        std::array<std::size_t, 3> indStart_xp{indEnd_m[0], indStart_m[1], indStart_m[2]};
        std::array<std::size_t, 3> indEnd_xp{indEnd_p[0], indEnd_p[1], indEnd_p[2]};

        NumberArray3D slice_xm = gridArray.GetSlice(indStart_xm, indEnd_xm);
        NumberArray3D slice_xp = gridArray.GetSlice(indStart_xp, indEnd_xp);

        // create meshgrid
        std::array<RealNumber, 3> r0_xm{(RealNumber)(indStart_xm[0])*dx,
                                        (RealNumber)(indStart_xm[1])*dy,
                                        (RealNumber)(indStart_xm[2])*dz};
        std::array<RealNumber, 3> r1_xm{(RealNumber)(indEnd_xm[0])*dx,
                                        (RealNumber)(indEnd_xm[1])*dy,
                                        (RealNumber)(indEnd_xm[2])*dz};
        NumberArray3D<RealNumber> xm = NumberArray3D<RealNumber>::GetMeshGrid(slice_xm.GetShape(), r0_xm, r1_xm, 0);
        slice_xm *= 0.5*NumberArray3D<RealNumber>::sin((xm - r0_xm[0])*M_PI/(r1_xm[0] - r0_xm[0]) - M_PI/2) + 0.5;

        std::array<RealNumber, 3> r0_xp{r0[0] + (RealNumber)(indStart_xp[0])*dx,
                                        r0[1] + (RealNumber)(indStart_xp[1])*dy,
                                        r0[2] + (RealNumber)(indStart_xp[2])*dz};

        std::array<RealNumber, 3> r1_xp{r0[0] + (RealNumber)(indEnd_xp[0])*dx,
                                        r0[1] + (RealNumber)(indEnd_xp[1])*dy,
                                        r0[2] + (RealNumber)(indEnd_xp[2])*dz};

        NumberArray3D<RealNumber> xp = NumberArray3D<RealNumber>::GetMeshGrid(slice_xp.GetShape(), r0_xp, r1_xp, 0);
        slice_xp *= 0.5*NumberArray3D<RealNumber>::cos((xp - r0_xp[0])*M_PI/(r1_xp[0] - r0_xp[0])) + 0.5;
    }
    if(std::real(smoothEdgeThickness[1]) > 0.0) {
    // get an slice for the cube face with normal along +y or -y

        std::array<std::size_t, 3> indStart_ym{indStart_m[0], indStart_m[1], indStart_m[2]};
        std::array<std::size_t, 3> indEnd_ym{indEnd_p[0], indStart_p[1], indEnd_p[2]};
        std::array<std::size_t, 3> indStart_yp{indStart_m[0], indEnd_m[1], indStart_m[2]};
        std::array<std::size_t, 3> indEnd_yp{indEnd_p[0], indEnd_p[1], indEnd_p[2]};

        NumberArray3D slice_ym = gridArray.GetSlice(indStart_ym, indEnd_ym);
        NumberArray3D slice_yp = gridArray.GetSlice(indStart_yp, indEnd_yp);
        // create meshgrid
        std::array<RealNumber, 3> r0_ym{(RealNumber)(indStart_ym[0])*dx,
                                        (RealNumber)(indStart_ym[1])*dy,
                                        (RealNumber)(indStart_ym[2])*dz};

        std::array<RealNumber, 3> r1_ym{(RealNumber)(indEnd_ym[0])*dx,
                                        (RealNumber)(indEnd_ym[1])*dy,
                                        (RealNumber)(indEnd_ym[2])*dz};

        NumberArray3D<RealNumber> ym = NumberArray3D<RealNumber>::GetMeshGrid(slice_ym.GetShape(), r0_ym, r1_ym, 1);
        slice_ym *= 0.5*NumberArray3D<RealNumber>::sin((ym - r0_ym[1])*M_PI/(r1_ym[1] - r0_ym[1]) - M_PI/2) + 0.5;

        std::array<RealNumber, 3> r0_yp{r0[0] + (RealNumber)(indStart_yp[0])*dx,
                                        r0[1] + (RealNumber)(indStart_yp[1])*dy,
                                        r0[2] + (RealNumber)(indStart_yp[2])*dz};

        std::array<RealNumber, 3> r1_yp{r0[0] + (RealNumber)(indEnd_yp[0])*dx,
                                        r0[1] + (RealNumber)(indEnd_yp[1])*dy,
                                        r0[2] + (RealNumber)(indEnd_yp[2])*dz};

        NumberArray3D<RealNumber> yp = NumberArray3D<RealNumber>::GetMeshGrid(slice_yp.GetShape(), r0_yp, r1_yp, 1);
        slice_yp *= 0.5*NumberArray3D<RealNumber>::cos((yp - r0_yp[1])*M_PI/(r1_yp[1] - r0_yp[1])) + 0.5;
    }
    if(std::real(smoothEdgeThickness[2]) > 0.0) {
    // get an slice for the cube face with normal along +z or -z

        std::array<std::size_t, 3> indStart_zm{indStart_m[0], indStart_m[1], indStart_m[2]};
        std::array<std::size_t, 3> indEnd_zm{indEnd_p[0], indEnd_p[1], indStart_p[2]};
        std::array<std::size_t, 3> indStart_zp{indStart_m[0], indStart_m[1], indEnd_m[2]};
        std::array<std::size_t, 3> indEnd_zp{indEnd_p[0], indEnd_p[1], indEnd_p[2]};

        NumberArray3D slice_zm = gridArray.GetSlice(indStart_zm, indEnd_zm);
        NumberArray3D slice_zp = gridArray.GetSlice(indStart_zp, indEnd_zp);
        // create meshgrid
        std::array<RealNumber, 3> r0_zm{(RealNumber)(indStart_zm[0])*dx,
                                        (RealNumber)(indStart_zm[1])*dy,
                                        (RealNumber)(indStart_zm[2])*dz};

        std::array<RealNumber, 3> r1_zm{(RealNumber)(indEnd_zm[0])*dx,
                                        (RealNumber)(indEnd_zm[1])*dy,
                                        (RealNumber)(indEnd_zm[2])*dz};
        NumberArray3D<RealNumber> zm = NumberArray3D<RealNumber>::GetMeshGrid(slice_zm.GetShape(), r0_zm, r1_zm, 2);
        slice_zm *= 0.5*NumberArray3D<RealNumber>::sin((zm - r0_zm[2])*M_PI/(r1_zm[2] - r0_zm[2]) - M_PI/2) + 0.5;

        std::array<RealNumber, 3> r0_zp{r0[0] + (RealNumber)(indStart_zp[0])*dx,
                                        r0[1] + (RealNumber)(indStart_zp[1])*dy,
                                        r0[2] + (RealNumber)(indStart_zp[2])*dz};
        std::array<RealNumber, 3> r1_zp{r0[0] + (RealNumber)(indEnd_zp[0])*dx,
                                        r0[1] + (RealNumber)(indEnd_zp[1])*dy,
                                        r0[2] + (RealNumber)(indEnd_zp[2])*dz};
        NumberArray3D<RealNumber> zp = NumberArray3D<RealNumber>::GetMeshGrid(slice_zp.GetShape(), r0_zp, r1_zp, 2);
        slice_zp *= 0.5*NumberArray3D<RealNumber>::cos((zp - r0_zp[2])*M_PI/(r1_zp[2] - r0_zp[2])) + 0.5;
    }

    gridArray *= (insideValue - outsideValue);
    gridArray += outsideValue;
}

RealNumber SpatialCubeGridArrayManipulator::CalculateTime(const RealNumber dt, const std::size_t timeIndex) {
    return (RealNumber)timeIndex*dt;
}


