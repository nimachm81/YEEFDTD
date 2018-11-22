
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

void GaussianSaceTimeGridArrayManipulator::SetAmplitude(const RealNumber amplitude) {
    GaussianSaceTimeGridArrayManipulator::amplitude = amplitude;
}

void GaussianSaceTimeGridArrayManipulator::SetSpaceTimeCenterPoint(const std::array<RealNumber, 4> st_center) {
    GaussianSaceTimeGridArrayManipulator::st_center = st_center;
}

void GaussianSaceTimeGridArrayManipulator::SetSpeceTimeDecayRate(const std::array<RealNumber, 4> st_decay_rate) {
    GaussianSaceTimeGridArrayManipulator::st_decay_rate = st_decay_rate;
}

void GaussianSaceTimeGridArrayManipulator::SetSpaceTimeModulationFrequency(
        const std::array<RealNumber, 4> st_modulationFrequency) {
    GaussianSaceTimeGridArrayManipulator::st_modulationFrequency = st_modulationFrequency;
}

void GaussianSaceTimeGridArrayManipulator::SetSpaceTimeModulationPhase(
        const std::array<RealNumber, 4> st_modulationPhase) {
    GaussianSaceTimeGridArrayManipulator::st_modulationPhase = st_modulationPhase;
}

void GaussianSaceTimeGridArrayManipulator::SetTimeOffsetFraction(const RealNumber offsetFraction) {
    timeOffsetFraction = offsetFraction;
}

RealNumber GaussianSaceTimeGridArrayManipulator::CalculateTime(const RealNumber dt, const std::size_t timeIndex) {
    return ((RealNumber)timeIndex + timeOffsetFraction) * dt;
}

void GaussianSaceTimeGridArrayManipulator::UpdateArray(const RealNumber t) {
    RealNumber gaussianValue_t = amplitude * std::exp(-(t - st_center[3])*(t - st_center[3]) *
                                                     (st_decay_rate[3]*st_decay_rate[3]))
                                           * std::cos(2.0*M_PI*st_modulationFrequency[3]*t + st_modulationPhase[3]);

    RealNumber*** arrayData = gridArray.GetArrayData();
    std::array<std::size_t, 3>& arrayShape = gridArray.GetShape();
    std::array<std::size_t, 3>& arrayIndStart = gridArray.GetIndStart();

    std::size_t n0 = arrayShape[0];
    std::size_t n1 = arrayShape[1];
    std::size_t n2 = arrayShape[2];
    std::size_t ind0 = arrayIndStart[0];
    std::size_t ind1 = arrayIndStart[1];
    std::size_t ind2 = arrayIndStart[2];

    RealNumber x0 = r0[0];
    RealNumber y0 = r0[1];
    RealNumber z0 = r0[2];

    RealNumber dx = dr[0];
    RealNumber dy = dr[1];
    RealNumber dz = dr[2];

    RealNumber x, y, z;

    for(std::size_t i0 = 0; i0 < n0; ++i0) {
        x = x0 + (RealNumber)i0*dx;
        RealNumber gaussianValue_x = std::exp(-(x - st_center[0])*(x - st_center[0]) *
                                            (st_decay_rate[0]*st_decay_rate[0]))
                                     * std::cos(2.0*M_PI*st_modulationFrequency[0]*x + st_modulationPhase[0]);
        for(std::size_t i1 = 0; i1 < n1; ++i1) {
            y = y0 + (RealNumber)i1*dy;
            RealNumber gaussianValue_y = std::exp(-(y - st_center[1])*(y - st_center[1]) *
                                                (st_decay_rate[1]*st_decay_rate[1]))
                                         * std::cos(2.0*M_PI*st_modulationFrequency[1]*y + st_modulationPhase[1]);
            for(std::size_t i2 = 0; i2 < n2; ++i2) {
                z = z0 + (RealNumber)i2*dz;
                RealNumber gaussianValue_z = std::exp(-(z - st_center[2])*(z - st_center[2]) *
                                                    (st_decay_rate[2]*st_decay_rate[2]))
                                             * std::cos(2.0*M_PI*st_modulationFrequency[2]*z + st_modulationPhase[2]);

                arrayData[ind0 + i0][ind1 + i1][ind2 + i2] = gaussianValue_t*
                                                             gaussianValue_x*gaussianValue_y*gaussianValue_z;
            }
        }
    }

}


