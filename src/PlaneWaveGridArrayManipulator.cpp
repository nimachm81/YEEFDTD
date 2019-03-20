
#include "PlaneWaveGridArrayManipulator.h"

PlaneWaveGridArrayManipulator::PlaneWaveGridArrayManipulator() {
    propagationDirection = std::array<FPNumber, 3>{0.0, 0.0, 1.0};
    velocity = 1.0;
}

void PlaneWaveGridArrayManipulator::SetPropagationVelocity(const FPNumber v) {
    velocity = v;
}

void PlaneWaveGridArrayManipulator::SetAmplitude(FPNumber amp) {
    amplitude = amp;
}

void PlaneWaveGridArrayManipulator::SetCenterTime(FPNumber t_c) {
    t_center = t_c;
}


void PlaneWaveGridArrayManipulator::SetPropagationDirection(const std::array<FPNumber, 3> direction) {
    FPNumber magnitude = std::sqrt(direction[0]*direction[0] + direction[1]*direction[1] + direction[2]*direction[2]);
    propagationDirection = std::array<FPNumber, 3>{direction[0]/magnitude,
                                                   direction[1]/magnitude,
                                                   direction[2]/magnitude
                                                   };
}

void PlaneWaveGridArrayManipulator::SetTimeOffsetFraction(const FPNumber offsetFraction) {
    timeOffsetFraction = offsetFraction;
}

FPNumber PlaneWaveGridArrayManipulator::CalculateTime(const FPNumber dt, const std::size_t timeIndex) {
    return ((FPNumber)timeIndex + timeOffsetFraction) * dt;
}

void PlaneWaveGridArrayManipulator::UpdateArray(const FPNumber t, GAManipulatorInstructionCode instruction) {

    FPNumber*** arrayData = gridArray.GetArrayData();
    std::array<std::size_t, 3>& arrayShape = gridArray.GetShape();
    std::array<std::size_t, 3>& arrayIndStart = gridArray.GetIndStart();

    std::size_t n0 = arrayShape[0];
    std::size_t n1 = arrayShape[1];
    std::size_t n2 = arrayShape[2];
    std::size_t ind0 = arrayIndStart[0];
    std::size_t ind1 = arrayIndStart[1];
    std::size_t ind2 = arrayIndStart[2];

    FPNumber x0 = r0[0];
    FPNumber y0 = r0[1];
    FPNumber z0 = r0[2];

    FPNumber dx = dr[0];
    FPNumber dy = dr[1];
    FPNumber dz = dr[2];

    FPNumber x, y, z;
    FPNumber x_ax, y_ay, z_az;

    FPNumber ax = propagationDirection[0];
    FPNumber ay = propagationDirection[1];
    FPNumber az = propagationDirection[2];

    for(std::size_t i0 = 0; i0 < n0; ++i0) {
        x = x0 + (FPNumber)i0*dx;
        x_ax = x * ax;
        for(std::size_t i1 = 0; i1 < n1; ++i1) {
            y = y0 + (FPNumber)i1*dy;
            y_ay = y * ay;
            for(std::size_t i2 = 0; i2 < n2; ++i2) {
                z = z0 + (FPNumber)i2*dz;
                z_az = z * az;

                if(instruction == GAManipulatorInstructionCode::Equal) {
                    arrayData[ind0 + i0][ind1 + i1][ind2 + i2] = amplitude * GetNormalizedTemporalWaveform(t - t_center - (x_ax + y_ay + z_az)/velocity);
                } else if(instruction == GAManipulatorInstructionCode::PlusEqual) {
                    arrayData[ind0 + i0][ind1 + i1][ind2 + i2] += amplitude * GetNormalizedTemporalWaveform(t - t_center - (x_ax + y_ay + z_az)/velocity);
                } else {
                    std::cout << "error: Not implemented!!!" << std::endl;
                    assert(false);
                }

            }
        }
    }
}


