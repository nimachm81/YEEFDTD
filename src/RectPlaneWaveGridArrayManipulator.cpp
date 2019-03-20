
#include "RectPlaneWaveGridArrayManipulator.h"

RectPlaneWaveGridArrayManipulator::RectPlaneWaveGridArrayManipulator() {

}

void RectPlaneWaveGridArrayManipulator::SetRectWidth(FPNumber width) {
    t_rectWidth = width;
}
void RectPlaneWaveGridArrayManipulator::SetRectEdgeWidth(FPNumber width) {
    t_edgeWidth = width;
}

void RectPlaneWaveGridArrayManipulator::SetModulationFrequency(FPNumber freq) {
    t_modulationFrequency = freq;
}

void RectPlaneWaveGridArrayManipulator::SetModulationPhase(FPNumber phase) {
    t_modulationPhase = phase;
}


FPNumber RectPlaneWaveGridArrayManipulator::GetNormalizedTemporalWaveform(FPNumber t) {
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
