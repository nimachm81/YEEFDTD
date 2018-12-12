
//#define _USE_MATH_DEFINES


#include <cmath>

#include "GaussianSpaceTimeGridArrayManipulator.h"

GaussianSaceTimeGridArrayManipulator::GaussianSaceTimeGridArrayManipulator() :
        amplitude(0.0),
        st_center{0.0, 0.0, 0.0, 0.0},
        st_decay_rate{0.0, 0.0, 0.0, 0.0},
        st_modulationFrequency{0.0, 0.0, 0.0, 0.0},
        st_modulationPhase{0.0, 0.0, 0.0, 0.0} {
}

void GaussianSaceTimeGridArrayManipulator::SetAmplitude(const FPNumber amplitude) {
    GaussianSaceTimeGridArrayManipulator::amplitude = amplitude;
}

void GaussianSaceTimeGridArrayManipulator::SetSpaceTimeCenterPoint(const std::array<FPNumber, 4> st_center) {
    GaussianSaceTimeGridArrayManipulator::st_center = st_center;
}

void GaussianSaceTimeGridArrayManipulator::SetSpeceTimeDecayRate(const std::array<FPNumber, 4> st_decay_rate) {
    GaussianSaceTimeGridArrayManipulator::st_decay_rate = st_decay_rate;
}

void GaussianSaceTimeGridArrayManipulator::SetSpaceTimeModulationFrequency(
        const std::array<FPNumber, 4> st_modulationFrequency) {
    GaussianSaceTimeGridArrayManipulator::st_modulationFrequency = st_modulationFrequency;
}

void GaussianSaceTimeGridArrayManipulator::SetSpaceTimeModulationPhase(
        const std::array<FPNumber, 4> st_modulationPhase) {
    GaussianSaceTimeGridArrayManipulator::st_modulationPhase = st_modulationPhase;
}

void GaussianSaceTimeGridArrayManipulator::SetTimeOffsetFraction(const FPNumber offsetFraction) {
    timeOffsetFraction = offsetFraction;
}

FPNumber GaussianSaceTimeGridArrayManipulator::CalculateTime(const FPNumber dt, const std::size_t timeIndex) {
    return ((FPNumber)timeIndex + timeOffsetFraction) * dt;
}

void GaussianSaceTimeGridArrayManipulator::UpdateArray(const FPNumber t, GAManipulatorInstructionCode instruction) {
    FPNumber gaussianValue_t = amplitude * std::exp(-(t - st_center[3])*(t - st_center[3]) *
                                                     (st_decay_rate[3]*st_decay_rate[3]))
                                           * std::cos((FPNumber)(2.0*M_PI)*st_modulationFrequency[3]*t + st_modulationPhase[3]);

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

    if(instruction != GAManipulatorInstructionCode::Equal) {
        std::cout << "Not implemented!!!" << std::endl;
        assert(false);
    }

    for(std::size_t i0 = 0; i0 < n0; ++i0) {
        x = x0 + (FPNumber)i0*dx;
        FPNumber gaussianValue_x = std::exp(-(x - st_center[0])*(x - st_center[0]) *
                                            (st_decay_rate[0]*st_decay_rate[0]))
                                     * std::cos((FPNumber)(2.0*M_PI)*st_modulationFrequency[0]*x + st_modulationPhase[0]);
        for(std::size_t i1 = 0; i1 < n1; ++i1) {
            y = y0 + (FPNumber)i1*dy;
            FPNumber gaussianValue_y = std::exp(-(y - st_center[1])*(y - st_center[1]) *
                                                (st_decay_rate[1]*st_decay_rate[1]))
                                         * std::cos((FPNumber)(2.0*M_PI)*st_modulationFrequency[1]*y + st_modulationPhase[1]);
            for(std::size_t i2 = 0; i2 < n2; ++i2) {
                z = z0 + (FPNumber)i2*dz;
                FPNumber gaussianValue_z = std::exp(-(z - st_center[2])*(z - st_center[2]) *
                                                    (st_decay_rate[2]*st_decay_rate[2]))
                                             * std::cos((FPNumber)(2.0*M_PI)*st_modulationFrequency[2]*z + st_modulationPhase[2]);

                arrayData[ind0 + i0][ind1 + i1][ind2 + i2] = gaussianValue_t*
                                                             gaussianValue_x*gaussianValue_y*gaussianValue_z;
            }
        }
    }

}


