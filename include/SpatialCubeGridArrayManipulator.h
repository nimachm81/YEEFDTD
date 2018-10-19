#ifndef FDTD_SPATIALCUBEGRIDARRAYMANIPULATOR_H_
#define FDTD_SPATIALCUBEGRIDARRAYMANIPULATOR_H_

#include <array>

#include "GridArrayManipulator.h"

class SpatialCubeGridArrayManipulator : public GridArrayManipulator{
    public:
    virtual ~SpatialCubeGridArrayManipulator() { };
    void SetCubeCorners(std::array<RealNumber, 3>& r0, std::array<RealNumber, 3>& r1);
    void SetEdgeThickness(std::array<RealNumber, 3>& thickness);
    void SetInsideValue(RealNumber value);
    void SetOutsideValue(RealNumber value);

    std::pair<std::array<std::size_t, 3>, std::array<std::size_t, 3>>
        GetArrayIndicesBetweenTwoCorners(std::array<RealNumber, 3>& r0, std::array<RealNumber, 3>& r1);

    void UpdateArray(const RealNumber t);
    RealNumber CalculateTime(const RealNumber dt, const std::size_t timeIndex);

    private:
    std::array<RealNumber, 3> cubeR0;       // lower left corner of the cube
    std::array<RealNumber, 3> cubeR1;       // upper right corner
    std::array<RealNumber, 3> smoothEdgeThickness;       // if non-zero the edges are smooth
    RealNumber insideValue;
    RealNumber outsideValue;

};


#endif // FDTD_SPATIALCUBEGRIDARRAYMANIPULATOR_H_


