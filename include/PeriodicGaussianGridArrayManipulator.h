
#ifndef FDTD_PERIODICGAUSSIANGRIDARRAYMANIPULATOR_H_
#define FDTD_PERIODICGAUSSIANGRIDARRAYMANIPULATOR_H_

#include <cassert>
#include <cmath>

#include "PeriodicGridArrayManipulator.h"

class PeriodicGaussianGridArrayManipulator : public PeriodicGridArrayManipulator {
    public:
    virtual ~PeriodicGaussianGridArrayManipulator() { };
    void SetGaussianAmplitude(RealNumber amplitude);
    void SetGaussianCenter(std::array<RealNumber, 3>& center);
    void SetGaussianDecayRate(std::array<RealNumber, 3>& decayRate);

    RealNumber Func(const std::array<RealNumber, 3> r, const RealNumber t);
    RealNumber CalculateTime(const RealNumber dt, const std::size_t timeIndex);

    private:
    RealNumber gaussianAmplitude;
    std::array<RealNumber, 3> gaussianCenter;
    std::array<RealNumber, 3> gaussianDecayRate;
};


#endif  // FDTD_PERIODICGAUSSIANGRIDARRAYMANIPULATOR_H_
