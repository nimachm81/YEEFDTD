
#include "PlaneWaveVectorField.h"


void PlaneWaveVectorField::SetPropagationVelocity(const FPNumber v) {
    velocity = v;
}

void PlaneWaveVectorField::SetPropagationDirection(const std::array<FPNumber, 3> direction) {
    FPNumber magnitude = std::sqrt(direction[0]*direction[0] + direction[1]*direction[1] + direction[2]*direction[2]);
    propagationDirection = std::array<FPNumber, 3>{direction[0]/magnitude,
                                                   direction[1]/magnitude,
                                                   direction[2]/magnitude
                                                   };
}

void PlaneWaveVectorField::SetAmplitude(std::array<FPNumber, 3>& amp) {
    amplitude = amp;
}

void PlaneWaveVectorField::SetCenterTime(FPNumber t_c) {
    t_center = t_c;
}

void PlaneWaveVectorField::GetFieldValueAtPoint(FPNumber time,
                      std::array<FPNumber, 3>& position,
                      std::array<FPNumber, 3>& fieldValue
                      ) {
    FPNumber u_dot_r = position[0]*propagationDirection[0] +
                       position[1]*propagationDirection[1] +
                       position[2]*propagationDirection[2];

    FPNumber t = (time - t_center) - u_dot_r / velocity;
    FPNumber waveformValue = GetNormalizedWaveform(t);

    fieldValue[0] = amplitude[0] * waveformValue;
    fieldValue[1] = amplitude[1] * waveformValue;
    fieldValue[2] = amplitude[2] * waveformValue;
}

void PlaneWaveVectorField::GetFieldValuesAtPoints(FPNumber time,
                      std::vector<std::array<FPNumber, 3>>& positions,
                      std::vector<std::array<FPNumber, 3>>& fieldValues
                      ) {
    std::size_t numPoints = positions.size();
    if(fieldValues.size() != numPoints) {
        fieldValues.resize(numPoints);
    }

    for(std::size_t i = 0; i < numPoints; ++i) {
        std::array<FPNumber, 3>& position = positions[i];
        std::array<FPNumber, 3>& fieldValue = fieldValues[i];


        FPNumber u_dot_r = position[0]*propagationDirection[0] +
                           position[1]*propagationDirection[1] +
                           position[2]*propagationDirection[2];

        FPNumber t = (time - t_center) - u_dot_r / velocity;
        FPNumber waveformValue = GetNormalizedWaveform(t);


        fieldValue[0] = amplitude[0] * waveformValue;
        fieldValue[1] = amplitude[1] * waveformValue;
        fieldValue[2] = amplitude[2] * waveformValue;
    }
}



