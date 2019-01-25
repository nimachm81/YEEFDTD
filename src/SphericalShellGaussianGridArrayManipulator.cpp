
#include "SpherialShellGaussianGridArrayManipulator.h"


void SpherialShellGaussianGridArrayManipulator::SetAmplitude(const FPNumber amplitude) {
    SpherialShellGaussianGridArrayManipulator::amplitude = amplitude;
}

void SpherialShellGaussianGridArrayManipulator::SetCenterPoint(const std::array<FPNumber, 3>& r_center) {
    SpherialShellGaussianGridArrayManipulator::centerPoint = r_center;
}

void SpherialShellGaussianGridArrayManipulator::SetRadius(const FPNumber radius) {
    SpherialShellGaussianGridArrayManipulator::radius = radius;
}


void SpherialShellGaussianGridArrayManipulator::SetDecayRate(const FPNumber r_decay_rate) {
    SpherialShellGaussianGridArrayManipulator::r_decay_rate = r_decay_rate;
}

void SpherialShellGaussianGridArrayManipulator::SetModulationFrequency(const FPNumber r_modulationFrequency) {
    SpherialShellGaussianGridArrayManipulator::r_modulationFrequency = r_modulationFrequency;
}

void SpherialShellGaussianGridArrayManipulator::SetModulationPhase(const FPNumber r_modulationPhase) {
    SpherialShellGaussianGridArrayManipulator::r_modulationPhase = r_modulationPhase;
}

void SpherialShellGaussianGridArrayManipulator::SetTimeOffsetFraction(const FPNumber offsetFraction) {
    SpherialShellGaussianGridArrayManipulator::timeOffsetFraction = offsetFraction;
}

FPNumber SpherialShellGaussianGridArrayManipulator::CalculateTime(const FPNumber dt, const std::size_t timeIndex) {
    return ((FPNumber)timeIndex + timeOffsetFraction) * dt;
}

void SpherialShellGaussianGridArrayManipulator::UpdateArray(const FPNumber t, GAManipulatorInstructionCode instruction) {

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

    FPNumber x, y, z, r;

    if(instruction != GAManipulatorInstructionCode::Equal) {
        std::cout << "Not implemented!!!" << std::endl;
        assert(false);
    }

    for(std::size_t i0 = 0; i0 < n0; ++i0) {
        x = x0 + (FPNumber)i0*dx - centerPoint[0];
        for(std::size_t i1 = 0; i1 < n1; ++i1) {
            y = y0 + (FPNumber)i1*dy - centerPoint[1];
            for(std::size_t i2 = 0; i2 < n2; ++i2) {
                z = z0 + (FPNumber)i2*dz - centerPoint[2];
                r = std::sqrt(x*x + y*y + z*z);
                FPNumber gaussianValue_r = std::exp(-(r - radius)*(r - radius) *
                                                    (r_decay_rate*r_decay_rate))
                * std::cos((FPNumber)(2.0*M_PI)*r_modulationFrequency*(r - radius) + r_modulationPhase);

                arrayData[ind0 + i0][ind1 + i1][ind2 + i2] = gaussianValue_r;
            }
        }
    }
}




