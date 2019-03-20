

#include "GaussianPlaneWaveVectorField.h"

void GaussianPlaneWaveVectorField::SetTimeDecayRate(FPNumber rate) {
    t_decayRate = rate;
}

void GaussianPlaneWaveVectorField::SetModulationFrequency(FPNumber freq) {
    t_modulationFrequency = freq;
}

void GaussianPlaneWaveVectorField::SetModulationPhase(FPNumber phase) {
    t_modulationPhase = phase;
}

FPNumber GaussianPlaneWaveVectorField::GetNormalizedWaveform(const FPNumber& t) {
    return std::exp(-t*t*t_decayRate*t_decayRate) * std::cos(2.0*M_PI*t_modulationFrequency*t + t_modulationPhase);
}


