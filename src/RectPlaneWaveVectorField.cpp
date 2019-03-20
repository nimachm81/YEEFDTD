
#include "RectPlaneWaveVectorField.h"

void RectPlaneWaveVectorField::SetRectWidth(FPNumber width) {
    t_rectWidth = width;
}
void RectPlaneWaveVectorField::SetRectEdgeWidth(FPNumber width) {
    t_edgeWidth = width;
}

void RectPlaneWaveVectorField::SetModulationFrequency(FPNumber freq) {
    t_modulationFrequency = freq;
}

void RectPlaneWaveVectorField::SetModulationPhase(FPNumber phase) {
    t_modulationPhase = phase;
}

FPNumber RectPlaneWaveVectorField::GetNormalizedWaveform(const FPNumber& t) {
    FPNumber amp = 0.0;
    if(std::abs(t) <= (t_rectWidth - t_edgeWidth)/2.0) {
        amp = 1.0;
    } else if(std::abs(t) <= (t_rectWidth + t_edgeWidth)/2.0) {
        if(t >= 0.0) {
            amp = 0.5*std::cos( (t - (t_rectWidth - t_edgeWidth)/2.0)*M_PI/t_edgeWidth ) + 0.5;
        } else {
            amp = 0.5*std::sin( (t + (t_rectWidth + t_edgeWidth)/2.0)*M_PI/t_edgeWidth - M_PI/2) + 0.5;
        }
    }
    return amp * std::cos(2.0*M_PI*t_modulationFrequency*t + t_modulationPhase);
}

