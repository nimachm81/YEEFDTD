
#include "GaussianPlaneWaveGridArrayManipulator.h"

GaussianPlaneWaveGridArrayManipulator::GaussianPlaneWaveGridArrayManipulator() :
        PlaneWaveGridArrayManipulator() {

}

void GaussianPlaneWaveGridArrayManipulator::SetTimeDecayRate(FPNumber rate) {
    t_decayRate = rate;
}

void GaussianPlaneWaveGridArrayManipulator::SetModulationFrequency(FPNumber freq) {
    t_modulationFrequency = freq;
}

void GaussianPlaneWaveGridArrayManipulator::SetModulationPhase(FPNumber phase) {
    t_modulationPhase = phase;
}

FPNumber GaussianPlaneWaveGridArrayManipulator::GetNormalizedTemporalWaveform(FPNumber t) {
    return std::exp(-t*t*t_decayRate*t_decayRate) * std::cos(2.0*M_PI*t_modulationFrequency*t + t_modulationPhase);
}



