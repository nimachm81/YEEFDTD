#ifndef FDTD_SPATIALCUBEGRIDARRAYMANIPULATOR_H_
#define FDTD_SPATIALCUBEGRIDARRAYMANIPULATOR_H_

#include <array>

#include "GridArrayManipulator.h"

class SpatialCubeGridArrayManipulator : public GridArrayManipulator{
    public:
    virtual ~SpatialCubeGridArrayManipulator() { };
    void SetCubeCorners(std::array<FPNumber, 3>& r0, std::array<FPNumber, 3>& r1);
    void SetEdgeThickness(std::array<FPNumber, 3>& thickness);
    void SetInsideValue(FPNumber value);
    void SetOutsideValue(FPNumber value);

    std::pair<std::array<std::size_t, 3>, std::array<std::size_t, 3>>
        GetArrayIndicesBetweenTwoCorners(std::array<FPNumber, 3>& r0, std::array<FPNumber, 3>& r1);

    void UpdateArray(const FPNumber t, GAManipulatorInstructionCode instruction);
    FPNumber CalculateTime(const FPNumber dt, const std::size_t timeIndex);

    private:
    std::array<FPNumber, 3> cubeR0;       // lower left corner of the cube
    std::array<FPNumber, 3> cubeR1;       // upper right corner
    std::array<FPNumber, 3> smoothEdgeThickness;       // if non-zero the edges are smooth
    FPNumber insideValue;
    FPNumber outsideValue;

};


#endif // FDTD_SPATIALCUBEGRIDARRAYMANIPULATOR_H_


