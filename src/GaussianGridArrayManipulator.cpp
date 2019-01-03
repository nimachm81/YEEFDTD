
//#define _USE_MATH_DEFINES


#include <cmath>

#include "GaussianGridArrayManipulator.h"

void GaussianGridArrayManipulator::SetAmplitude(const FPNumber amplitude) {
    GaussianGridArrayManipulator::amplitude = amplitude;
}

void GaussianGridArrayManipulator::SetCenterTime(const FPNumber t_center) {
    GaussianGridArrayManipulator::t_center = t_center;
}

void GaussianGridArrayManipulator::SetDecayTime(const FPNumber t_decay) {
    GaussianGridArrayManipulator::t_decay = t_decay;
}

void GaussianGridArrayManipulator::SetModulationFrequency(const FPNumber modulationFrequency) {
    GaussianGridArrayManipulator::modulationFrequency = modulationFrequency;
}

void GaussianGridArrayManipulator::SetModulationPhase(const FPNumber modulationPhase) {
    GaussianGridArrayManipulator::modulationPhase = modulationPhase;
}

void GaussianGridArrayManipulator::SetTimeOffsetFraction(const FPNumber offsetFraction) {
    timeOffsetFraction = offsetFraction;
}

FPNumber GaussianGridArrayManipulator::CalculateTime(const FPNumber dt, const std::size_t timeIndex) {
    return ((FPNumber)timeIndex + timeOffsetFraction) * dt;
}

void GaussianGridArrayManipulator::UpdateArray(const FPNumber t, GAManipulatorInstructionCode instruction) {
    FPNumber gaussianValue = amplitude * std::exp(-(t - t_center)*(t - t_center) / (t_decay*t_decay)) *
                               std::cos((FPNumber)(2.0*M_PI)*modulationFrequency*(t - t_center) + modulationPhase);
    if(instruction == GAManipulatorInstructionCode::Equal) {
        gridArray = gaussianValue;
    } else if(instruction == GAManipulatorInstructionCode::PlusEqual) {
        gridArray += gaussianValue;
    } else if(instruction == GAManipulatorInstructionCode::MultiplyEqual) {
        gridArray *= gaussianValue;
    } else {
        std::cout << "Not implemented!!!" << std::endl;
        assert(false);
    }
}

