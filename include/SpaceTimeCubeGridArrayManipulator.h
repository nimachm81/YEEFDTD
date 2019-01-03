
#ifndef FDTD_SPACETIMECUBEGRIDARRAYMANIPULATOR_H_
#define FDTD_SPACETIMECUBEGRIDARRAYMANIPULATOR_H_

#include <array>

#include "GridArrayManipulator.h"

class SpaceTimeCubeGridArrayManipulator : public GridArrayManipulator{
    public:
    virtual ~SpaceTimeCubeGridArrayManipulator() { };
    void SetCubeCorners(std::array<FPNumber, 4>& r0, std::array<FPNumber, 4>& r1);
    void SetEdgeThickness(std::array<FPNumber, 4>& thickness);
    void SetInsideValue(FPNumber value);
    void SetOutsideValue(FPNumber value);
    void SetTimeOffsetFraction(const FPNumber offsetFraction);

    std::pair<std::array<std::size_t, 3>, std::array<std::size_t, 3>>
        GetArrayIndicesBetweenTwoCorners(std::array<FPNumber, 3>& r0, std::array<FPNumber, 3>& r1);

    void UpdateArray(const FPNumber t, GAManipulatorInstructionCode instruction);
    FPNumber CalculateTime(const FPNumber dt, const std::size_t timeIndex);

    private:
    std::array<FPNumber, 4> cubeR0;       // lower left corner of the hypercube
    std::array<FPNumber, 4> cubeR1;       // upper right corner
    std::array<FPNumber, 4> smoothEdgeThickness;       // if non-zero the edges are smooth
    FPNumber insideValue;
    FPNumber outsideValue;
    FPNumber timeOffsetFraction = 0.0;

    int temporalState = -1; // it determines whether the array has to be recalculated. Only at the rising and falling states (1, 3)
                            // the array has to be recalculated at each time step
                            // < 0   : undefined
                            // 0     : the array's initial state has been set (t < cubeR0[3])
                            // 1     : rising state
                            // 2     : middle plateau
                            // 3     : falling state
                            // 4     : final state has been set (t > cubeR1[3])
};


#endif // FDTD_SPACETIMECUBEGRIDARRAYMANIPULATOR_H_


