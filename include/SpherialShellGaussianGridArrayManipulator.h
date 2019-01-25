
#ifndef FDTD_SPHERICALSHELLGAUSSIANGRIDARRAYMANIPULATOR_H_
#define FDTD_SPHERICALSHELLGAUSSIANGRIDARRAYMANIPULATOR_H_

#include <cassert>
#include <cmath>

#include "GridArrayManipulator.h"

class SpherialShellGaussianGridArrayManipulator : public GridArrayManipulator {
    public:
    virtual ~SpherialShellGaussianGridArrayManipulator() { };
    void SetAmplitude(const FPNumber amplitude);
    void SetCenterPoint(const std::array<FPNumber, 3>& r_center);
    void SetRadius(const FPNumber radius);
    void SetDecayRate(const FPNumber r_decay_rate);
    void SetModulationFrequency(const FPNumber r_modulationFrequency);
    void SetModulationPhase(const FPNumber r_modulationPhase);
    void SetTimeOffsetFraction(const FPNumber offsetFraction);

    void UpdateArray(const FPNumber t, GAManipulatorInstructionCode instruction);
    FPNumber CalculateTime(const FPNumber dt, const std::size_t timeIndex);

    private:
    FPNumber amplitude;
    std::array<FPNumber, 3> centerPoint;
    FPNumber radius;
    FPNumber r_decay_rate;
    FPNumber r_modulationFrequency;
    FPNumber r_modulationPhase;
    FPNumber timeOffsetFraction = 0.0;

};


#endif // FDTD_SPHERICALSHELLGAUSSIANGRIDARRAYMANIPULATOR_H_


