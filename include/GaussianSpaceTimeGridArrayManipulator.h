
#ifndef FDTD_GAUSSIANSPACETMEGRIDARRAYMANIPULATOR_H_
#define FDTD_GAUSSIANSPACETMEGRIDARRAYMANIPULATOR_H_

#include <cassert>
#include <cmath>

#include "GridArrayManipulator.h"

class GaussianSaceTimeGridArrayManipulator : public GridArrayManipulator {
    public:
    GaussianSaceTimeGridArrayManipulator();
    virtual ~GaussianSaceTimeGridArrayManipulator() { };
    void SetAmplitude(const RealNumber amplitude);
    void SetSpaceTimeCenterPoint(const std::array<RealNumber, 4> st_center);
    void SetSpeceTimeDecayRate(const std::array<RealNumber, 4> st_decay_rate);
    void SetSpaceTimeModulationFrequency(const std::array<RealNumber, 4> modulationFrequency);
    void SetSpaceTimeModulationPhase(const std::array<RealNumber, 4> modulationPhase);
    void SetTimeOffsetFraction(const RealNumber offsetFraction);

    void UpdateArray(const RealNumber t);
    RealNumber CalculateTime(const RealNumber dt, const std::size_t timeIndex);

    private:
    RealNumber amplitude;
    std::array<RealNumber, 4> st_center;
    std::array<RealNumber, 4> st_decay_rate;
    std::array<RealNumber, 4> st_modulationFrequency;
    std::array<RealNumber, 4> st_modulationPhase;
    RealNumber timeOffsetFraction = 0.0;

};


#endif // FDTD_GAUSSIANSPACETMEGRIDARRAYMANIPULATOR_H_

