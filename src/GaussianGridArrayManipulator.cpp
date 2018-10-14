
//#define _USE_MATH_DEFINES


#include <cmath>

#include "GaussianGridArrayManipulator.h"

void GaussianGridArrayManipulator::SetAmplitude(const RealNumber amplitude) {
    GaussianGridArrayManipulator::amplitude = amplitude;
}

void GaussianGridArrayManipulator::SetCenterTime(const RealNumber t_center) {
    GaussianGridArrayManipulator::t_center = t_center;
}

void GaussianGridArrayManipulator::SetDecayTime(const RealNumber t_decay) {
    GaussianGridArrayManipulator::t_decay = t_decay;
}

void GaussianGridArrayManipulator::SetModulationFrequency(const RealNumber modulationFrequency) {
    GaussianGridArrayManipulator::modulationFrequency = modulationFrequency;
}

void GaussianGridArrayManipulator::SetModulationPhase(const RealNumber modulationPhase) {
    GaussianGridArrayManipulator::modulationPhase = modulationPhase;
}

void GaussianGridArrayManipulator::SetTimeOffsetFraction(const RealNumber offsetFraction) {
    timeOffsetFraction = offsetFraction;
}

RealNumber GaussianGridArrayManipulator::CalculateTime(const RealNumber dt, const std::size_t timeIndex) {
    return (timeIndex + timeOffsetFraction) * dt;
}

void GaussianGridArrayManipulator::UpdateArray(const RealNumber t) {
    RealNumber gaussianValue = amplitude * std::exp(-(t - t_center)*(t - t_center) / (t_decay*t_decay)) *
                               std::cos(2.0*M_PI*modulationFrequency*t + modulationPhase);
    gridArray.SetToNumber(gaussianValue);
}

