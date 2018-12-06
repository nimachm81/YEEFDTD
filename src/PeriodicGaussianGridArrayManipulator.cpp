

#include "PeriodicGaussianGridArrayManipulator.h"

void PeriodicGaussianGridArrayManipulator::SetGaussianAmplitude(FPNumber amplitude) {
    gaussianAmplitude = amplitude;
}

void PeriodicGaussianGridArrayManipulator::SetGaussianCenter(std::array<FPNumber, 3>& center) {
    gaussianCenter = center;
}

void PeriodicGaussianGridArrayManipulator::SetGaussianDecayRate(std::array<FPNumber, 3>& decayRate) {
    gaussianDecayRate = decayRate;
}

FPNumber PeriodicGaussianGridArrayManipulator::Func(const std::array<FPNumber, 3> r, const FPNumber t) {
    return gaussianAmplitude * std::exp( -(
        (r[0] - gaussianCenter[0])*(r[0] - gaussianCenter[0])*(gaussianDecayRate[0]*gaussianDecayRate[0]) +
        (r[1] - gaussianCenter[1])*(r[1] - gaussianCenter[1])*(gaussianDecayRate[1]*gaussianDecayRate[1]) +
        (r[2] - gaussianCenter[2])*(r[2] - gaussianCenter[2])*(gaussianDecayRate[2]*gaussianDecayRate[2])
                                          )
                                       );
}

FPNumber PeriodicGaussianGridArrayManipulator::CalculateTime(const FPNumber dt, const std::size_t timeIndex) {
    return dt*(FPNumber)timeIndex;
}



