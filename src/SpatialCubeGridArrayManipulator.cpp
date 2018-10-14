
#include <cassert>
#include <cmath>

#include "SpatialCubeGridArrayManipulator.h"


void SpatialCubeGridArrayManipulator::SetCubeCorners(std::array<RealNumber, 3>& r0, std::array<RealNumber, 3>& r1) {
    assert(r0[0] <= r1[0] && r0[1] <= r1[1] && r0[2] <= r1[2]);
    cubeR0 = r0;
    cubeR1 = r1;
}

void SpatialCubeGridArrayManipulator::SetInsideValue(RealNumber value) {
    insideValue = value;
}

void SpatialCubeGridArrayManipulator::SetOutsideValue(RealNumber value) {
    outsideValue = value;
}

void SpatialCubeGridArrayManipulator::UpdateArray(const RealNumber t) {
    gridArray.SetToNumber(outsideValue);

    // get an slice on the cube and set its value to insideValue
    std::array<std::size_t, 3> shape = gridArray.GetShape();
    RealNumber dx = (r1[0] - r0[0])/shape[0];
    RealNumber dy = (r1[1] - r0[1])/shape[1];
    RealNumber dz = (r1[2] - r0[2])/shape[2];

    std::array<std::size_t, 3> indStart = {(std::size_t)std::ceil((cubeR0[0] - r0[0])/dx),
                                           (std::size_t)std::ceil((cubeR0[1] - r0[1])/dy),
                                           (std::size_t)std::ceil((cubeR0[2] - r0[2])/dz)};
    if(cubeR0[0] <= r0[0]) {indStart[0] = 0;}
    if(cubeR0[1] <= r0[1]) {indStart[1] = 0;}
    if(cubeR0[2] <= r0[2]) {indStart[2] = 0;}

    std::array<std::size_t, 3> indEnd = {(std::size_t)std::floor((cubeR1[0] - r0[0])/dx),
                                         (std::size_t)std::floor((cubeR1[1] - r0[1])/dy),
                                         (std::size_t)std::floor((cubeR1[2] - r0[2])/dz)};

    if(indEnd[0] > shape[0]) {indEnd[0] = shape[0];}
    if(indEnd[1] > shape[1]) {indEnd[1] = shape[1];}
    if(indEnd[2] > shape[2]) {indEnd[2] = shape[2];}
    if(cubeR1[0] <= r0[0]) {indEnd[0] = 0;}
    if(cubeR1[1] <= r0[1]) {indEnd[1] = 0;}
    if(cubeR1[2] <= r0[2]) {indEnd[2] = 0;}

    assert(indStart[0] <= indEnd[0] && indStart[1] <= indEnd[1] && indStart[2] <= indEnd[2]);
    assert(indStart[0] < indEnd[0] || indStart[1] < indEnd[1] || indStart[2] < indEnd[2]);

    NumberArray3D slice = gridArray.GetSlice(indStart, indEnd);
    slice.SetToNumber(insideValue);
}

RealNumber SpatialCubeGridArrayManipulator::CalculateTime(const RealNumber dt, const std::size_t timeIndex) {
    return timeIndex*dt;
}


