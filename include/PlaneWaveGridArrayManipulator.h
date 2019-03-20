
#ifndef FDTD_PLANEWAVEGRIDARRAYMANIPULATOR_H_
#define FDTD_PLANEWAVEGRIDARRAYMANIPULATOR_H_

#include <cassert>
#include <cmath>

#include "GridArrayManipulator.h"

// f(t - k.r/v) where k is the propagation direction
class PlaneWaveGridArrayManipulator : public GridArrayManipulator {
    public:
    PlaneWaveGridArrayManipulator();
    virtual ~PlaneWaveGridArrayManipulator() { };

    void SetAmplitude(FPNumber amp);
    void SetCenterTime(FPNumber t_c);   void SetPropagationVelocity(const FPNumber v);
    void SetPropagationDirection(const std::array<FPNumber, 3> direction);
    void SetTimeOffsetFraction(const FPNumber offsetFraction);

    virtual FPNumber GetNormalizedTemporalWaveform(FPNumber t) = 0;

    void UpdateArray(const FPNumber t, GAManipulatorInstructionCode instruction);
    FPNumber CalculateTime(const FPNumber dt, const std::size_t timeIndex);

    private:
    std::array<FPNumber, 3> propagationDirection;
    FPNumber velocity;

    FPNumber amplitude = 1.0;
    FPNumber t_center = 0.0;

    FPNumber timeOffsetFraction = 0.0;
};


#endif // FDTD_PLANEWAVEGRIDARRAYMANIPULATOR_H_


