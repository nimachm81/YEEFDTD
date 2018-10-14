#ifndef FDTD_SPATIALCUBEGRIDARRAYMANIPULATOR_H_
#define FDTD_SPATIALCUBEGRIDARRAYMANIPULATOR_H_

#include <array>

#include "GridArrayManipulator.h"

class SpatialCubeGridArrayManipulator : public GridArrayManipulator{
    public:
    virtual ~SpatialCubeGridArrayManipulator() { };
    void SetCubeCorners(std::array<RealNumber, 3>& r0, std::array<RealNumber, 3>& r1);
    void SetInsideValue(RealNumber value);
    void SetOutsideValue(RealNumber value);

    void UpdateArray(const RealNumber t);
    RealNumber CalculateTime(const RealNumber dt, const std::size_t timeIndex);

    private:
    std::array<RealNumber, 3> cubeR0;       // lower left corner of the cube
    std::array<RealNumber, 3> cubeR1;       // upper right corner
    RealNumber insideValue;
    RealNumber outsideValue;

};


#endif // FDTD_SPATIALCUBEGRIDARRAYMANIPULATOR_H_


