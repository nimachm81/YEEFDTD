

#include "PeriodicGaussianGridArrayManipulator.h"

void PeriodicGaussianGridArrayManipulator::SetGaussianAmplitude(RealNumber amplitude) {
    gaussianAmplitude = amplitude;
}

void PeriodicGaussianGridArrayManipulator::SetGaussianCenter(std::array<RealNumber, 3>& center) {
    gaussianCenter = center;
}

void PeriodicGaussianGridArrayManipulator::SetGaussianDecayRate(std::array<RealNumber, 3>& decayRate) {
    gaussianDecayRate = decayRate;
}

RealNumber PeriodicGaussianGridArrayManipulator::Func(const std::array<RealNumber, 3> r, const RealNumber t) {
    return gaussianAmplitude * std::exp( -(
        (r[0] - gaussianCenter[0])*(r[0] - gaussianCenter[0])*(gaussianDecayRate[0]*gaussianDecayRate[0]) +
        (r[1] - gaussianCenter[1])*(r[1] - gaussianCenter[1])*(gaussianDecayRate[1]*gaussianDecayRate[1]) +
        (r[2] - gaussianCenter[2])*(r[2] - gaussianCenter[2])*(gaussianDecayRate[2]*gaussianDecayRate[2])
                                          )
                                       );
}

RealNumber PeriodicGaussianGridArrayManipulator::CalculateTime(const RealNumber dt, const std::size_t timeIndex) {
    return dt*timeIndex;
}



