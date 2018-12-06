
#ifndef FDTD_GAUSSIANSPACETMEGRIDARRAYMANIPULATOR_H_
#define FDTD_GAUSSIANSPACETMEGRIDARRAYMANIPULATOR_H_

#include <cassert>
#include <cmath>

#include "GridArrayManipulator.h"

class GaussianSaceTimeGridArrayManipulator : public GridArrayManipulator {
    public:
    GaussianSaceTimeGridArrayManipulator();
    virtual ~GaussianSaceTimeGridArrayManipulator() { };
    void SetAmplitude(const FPNumber amplitude);
    void SetSpaceTimeCenterPoint(const std::array<FPNumber, 4> st_center);
    void SetSpeceTimeDecayRate(const std::array<FPNumber, 4> st_decay_rate);
    void SetSpaceTimeModulationFrequency(const std::array<FPNumber, 4> modulationFrequency);
    void SetSpaceTimeModulationPhase(const std::array<FPNumber, 4> modulationPhase);
    void SetTimeOffsetFraction(const FPNumber offsetFraction);

    void UpdateArray(const FPNumber t);
    FPNumber CalculateTime(const FPNumber dt, const std::size_t timeIndex);

    private:
    FPNumber amplitude;
    std::array<FPNumber, 4> st_center;
    std::array<FPNumber, 4> st_decay_rate;
    std::array<FPNumber, 4> st_modulationFrequency;
    std::array<FPNumber, 4> st_modulationPhase;
    FPNumber timeOffsetFraction = 0.0;

};


#endif // FDTD_GAUSSIANSPACETMEGRIDARRAYMANIPULATOR_H_

