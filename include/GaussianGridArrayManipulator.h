
#ifndef FDTD_GAUSSIANGRIDARRAYMANIPULATOR_H_
#define FDTD_GAUSSIANGRIDARRAYMANIPULATOR_H_

#include <cassert>
#include <cmath>

#include "GridArrayManipulator.h"

class GaussianGridArrayManipulator : public GridArrayManipulator {
    public:
    virtual ~GaussianGridArrayManipulator() { };
    void SetAmplitude(const FPNumber amplitude);
    void SetCenterTime(const FPNumber t_center);
    void SetDecayTime(const FPNumber t_decay);
    void SetModulationFrequency(const FPNumber modulationFrequency);
    void SetModulationPhase(const FPNumber modulationPhase);
    void SetTimeOffsetFraction(const FPNumber offsetFraction);

    void UpdateArray(const FPNumber t, GAManipulatorInstructionCode instruction);
    FPNumber CalculateTime(const FPNumber dt, const std::size_t timeIndex);

    private:
    FPNumber amplitude;
    FPNumber t_center;
    FPNumber t_decay;
    FPNumber modulationFrequency;
    FPNumber modulationPhase;
    FPNumber timeOffsetFraction = 0.0;

};


#endif // FDTD_GAUSSIANGRIDARRAYMANIPULATOR_H_

