
#ifndef FDTD_PERIODICGAUSSIANGRIDARRAYMANIPULATOR_H_
#define FDTD_PERIODICGAUSSIANGRIDARRAYMANIPULATOR_H_

#include <cassert>
#include <cmath>

#include "PeriodicGridArrayManipulator.h"

class PeriodicGaussianGridArrayManipulator : public PeriodicGridArrayManipulator {
    public:
    virtual ~PeriodicGaussianGridArrayManipulator() { };
    void SetGaussianAmplitude(FPNumber amplitude);
    void SetGaussianCenter(std::array<FPNumber, 3>& center);
    void SetGaussianDecayRate(std::array<FPNumber, 3>& decayRate);

    FPNumber Func(const std::array<FPNumber, 3> r, const FPNumber t);
    FPNumber CalculateTime(const FPNumber dt, const std::size_t timeIndex);

    private:
    FPNumber gaussianAmplitude;
    std::array<FPNumber, 3> gaussianCenter;
    std::array<FPNumber, 3> gaussianDecayRate;
};


#endif  // FDTD_PERIODICGAUSSIANGRIDARRAYMANIPULATOR_H_
