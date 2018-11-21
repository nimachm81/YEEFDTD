


#include "PeriodicGridArrayManipulator.h"

void PeriodicGridArrayManipulator::SetUnitCellOrigin(const std::array<RealNumber, 3> origin) {
    unitCellOrigin = origin;
}

void PeriodicGridArrayManipulator::SetPrimitiveVectors(const std::array<RealNumber, 3> v0,
                             const std::array<RealNumber, 3> v1,
                             const std::array<RealNumber, 3> v2) {

    primitiveVectors[0] = v0;
    primitiveVectors[1] = v1;
    primitiveVectors[2] = v2;
}


void PeriodicGridArrayManipulator::UpdateArray(const RealNumber t) {
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
        x = x0 + i0*dx;
        for(std::size_t i1 = 0; i1 < n1; ++i1) {
            y = y0 + i1*dy;
            for(std::size_t i2 = 0; i2 < n2; ++i2) {
                z = z0 + i2*dz;
                std::array<RealNumber, 3> r{x, y, z};
                std::array<RealNumber, 3> r_in = BringPointInsideUnitCell(r);
                arrayData[ind0 + i0][ind1 + i1][ind2 + i2] = Func(r_in, t);
            }
        }
    }
}

std::array<RealNumber, 3> PeriodicGridArrayManipulator::BringPointInsideUnitCell(const std::array<RealNumber, 3>& r) {
    std::array<RealNumber, 3>& v0 = primitiveVectors[0];
    RealNumber r_dot_v0 = (r[0] - unitCellOrigin[0])*v0[0] +
                          (r[1] - unitCellOrigin[1])*v0[1] +
                          (r[2] - unitCellOrigin[2])*v0[2];
    RealNumber v0_len = std::sqrt(v0[0]*v0[0] + v0[1]*v0[1] + v0[2]*v0[2]);
    RealNumber rem_0 = std::fmod(r_dot_v0, v0_len);
    if(rem_0 < 0.0) {
        rem_0 += v0_len;
    }

    std::array<RealNumber, 3>& v1 = primitiveVectors[1];
    RealNumber r_dot_v1 = (r[0] - unitCellOrigin[0])*v1[0] +
                          (r[1] - unitCellOrigin[1])*v1[1] +
                          (r[2] - unitCellOrigin[2])*v1[2];
    RealNumber v1_len = std::sqrt(v1[0]*v1[0] + v1[1]*v1[1] + v1[2]*v1[2]);
    RealNumber rem_1 = std::fmod(r_dot_v1, v1_len);
    if(rem_1 < 0.0) {
        rem_1 += v1_len;
    }

    std::array<RealNumber, 3>& v2 = primitiveVectors[2];
    RealNumber r_dot_v2 = (r[0] - unitCellOrigin[0])*v2[0] +
                          (r[1] - unitCellOrigin[1])*v2[1] +
                          (r[2] - unitCellOrigin[2])*v2[2];
    RealNumber v2_len = std::sqrt(v2[0]*v2[0] + v2[1]*v2[1] + v2[2]*v2[2]);
    RealNumber rem_2 = std::fmod(r_dot_v2, v2_len);
    if(rem_2 < 0.0) {
        rem_2 += v2_len;
    }

    return std::array<RealNumber, 3>{rem_0*v0[0] + rem_1*v1[0] + rem_2*v2[0],
                                     rem_0*v0[1] + rem_1*v1[1] + rem_2*v2[1],
                                     rem_0*v0[2] + rem_1*v1[2] + rem_2*v2[2]};
}

