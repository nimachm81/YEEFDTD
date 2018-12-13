


#include "PeriodicGridArrayManipulator.h"

void PeriodicGridArrayManipulator::SetUnitCellOrigin(const std::array<FPNumber, 3> origin) {
    unitCellOrigin = origin;
}

void PeriodicGridArrayManipulator::SetPrimitiveVectors(const std::array<FPNumber, 3> v0,
                             const std::array<FPNumber, 3> v1,
                             const std::array<FPNumber, 3> v2) {

    primitiveVectors[0] = v0;
    primitiveVectors[1] = v1;
    primitiveVectors[2] = v2;
}


void PeriodicGridArrayManipulator::UpdateArray(const FPNumber t, GAManipulatorInstructionCode instruction) {
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
        for(std::size_t i1 = 0; i1 < n1; ++i1) {
            y = y0 + (FPNumber)i1*dy;
            for(std::size_t i2 = 0; i2 < n2; ++i2) {
                z = z0 + (FPNumber)i2*dz;
                std::array<FPNumber, 3> r{x, y, z};
                std::array<FPNumber, 3> r_in = BringPointInsideUnitCell(r);
                arrayData[ind0 + i0][ind1 + i1][ind2 + i2] = Func(r_in, t);
            }
        }
    }
}

std::array<FPNumber, 3> PeriodicGridArrayManipulator::BringPointInsideUnitCell(const std::array<FPNumber, 3>& r) {
    std::array<FPNumber, 3>& v0 = primitiveVectors[0];
    FPNumber r_dot_v0 = (r[0])*v0[0] +
                        (r[1])*v0[1] +
                        (r[2])*v0[2];
    FPNumber v0_len = std::sqrt(v0[0]*v0[0] + v0[1]*v0[1] + v0[2]*v0[2]);
    FPNumber rem_0 = std::fmod(std::real(r_dot_v0), std::real(v0_len));
    if(std::real(rem_0) < -0.5*v0_len) {
        rem_0 += v0_len;
    }
    if(std::real(rem_0) > 0.5*v0_len) {
        rem_0 -= v0_len;
    }

    std::array<FPNumber, 3>& v1 = primitiveVectors[1];
    FPNumber r_dot_v1 = (r[0])*v1[0] +
                        (r[1])*v1[1] +
                        (r[2])*v1[2];
    FPNumber v1_len = std::sqrt(v1[0]*v1[0] + v1[1]*v1[1] + v1[2]*v1[2]);
    FPNumber rem_1 = std::fmod(std::real(r_dot_v1), std::real(v1_len));
    if(std::real(rem_1) < -0.5*v1_len) {
        rem_1 += v1_len;
    }
    if(std::real(rem_1) > 0.5*v1_len) {
        rem_1 -= v1_len;
    }

    std::array<FPNumber, 3>& v2 = primitiveVectors[2];
    FPNumber r_dot_v2 = (r[0])*v2[0] +
                        (r[1])*v2[1] +
                        (r[2])*v2[2];
    FPNumber v2_len = std::sqrt(v2[0]*v2[0] + v2[1]*v2[1] + v2[2]*v2[2]);
    FPNumber rem_2 = std::fmod(std::real(r_dot_v2), std::real(v2_len));
    if(std::real(rem_2) < -0.5*v2_len) {
        rem_2 += v2_len;
    }
    if(std::real(rem_2) > 0.5*v2_len) {
        rem_2 -= v2_len;
    }

    return std::array<FPNumber, 3>{rem_0*v0[0] + rem_1*v1[0] + rem_2*v2[0],
                                   rem_0*v0[1] + rem_1*v1[1] + rem_2*v2[1],
                                   rem_0*v0[2] + rem_1*v1[2] + rem_2*v2[2]};
}

