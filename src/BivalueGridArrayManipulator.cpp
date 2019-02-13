

#include "BivalueGridArrayManipulator.h"

void BivalueGridArrayManipulator::SetupConditionArray() {
    conditionArray.ReInitialize(gridArray.GetShape(), 0);
}

void BivalueGridArrayManipulator::SetGeometry(std::shared_ptr<Geometry> geometryPtr) {
    geometry = geometryPtr;
}

void BivalueGridArrayManipulator::SetInsideValue(FPNumber value) {
    valueInside = value;
}

void BivalueGridArrayManipulator::SetOutsideValue(FPNumber value) {
    valueOutside = value;
}

FPNumber BivalueGridArrayManipulator::CalculateTime(const FPNumber dt, const std::size_t timeIndex) {
    return ((FPNumber)timeIndex + timeOffsetFraction) * dt;
}

/*
void BivalueGridArrayManipulator::UpdateArray(const FPNumber t, GAManipulatorInstructionCode instruction) {
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

    if(instruction == GAManipulatorInstructionCode::Equal) {
        for(std::size_t i0 = 0; i0 < n0; ++i0) {
            x = x0 + (FPNumber)i0*dx;
            for(std::size_t i1 = 0; i1 < n1; ++i1) {
                y = y0 + (FPNumber)i1*dy;
                for(std::size_t i2 = 0; i2 < n2; ++i2) {
                    z = z0 + (FPNumber)i2*dz;
                    if(geometry->IsPointInsideOrOn({x, y, z})) {
                        arrayData[ind0 + i0][ind1 + i1][ind2 + i2] = valueInside;
                    } else {
                        arrayData[ind0 + i0][ind1 + i1][ind2 + i2] = valueOutside;
                    }
                }
            }
        }
    } else if(instruction == GAManipulatorInstructionCode::MultiplyEqual) {
        for(std::size_t i0 = 0; i0 < n0; ++i0) {
            x = x0 + (FPNumber)i0*dx;
            for(std::size_t i1 = 0; i1 < n1; ++i1) {
                y = y0 + (FPNumber)i1*dy;
                for(std::size_t i2 = 0; i2 < n2; ++i2) {
                    z = z0 + (FPNumber)i2*dz;
                    if(geometry->IsPointInsideOrOn({x, y, z})) {
                        arrayData[ind0 + i0][ind1 + i1][ind2 + i2] *= valueInside;
                    } else {
                        arrayData[ind0 + i0][ind1 + i1][ind2 + i2] *= valueOutside;
                    }
                }
            }
        }
    } else {
            std::cout << "error: GAManipulatorInstructionCode operation not implemented." << std::endl;
            assert(false);
    }
}
*/

void BivalueGridArrayManipulator::UpdateArray(const FPNumber t, GAManipulatorInstructionCode instruction) {
    FPNumber*** arrayData = gridArray.GetArrayData();
    std::array<std::size_t, 3>& arrayShape = gridArray.GetShape();
    std::array<std::size_t, 3>& arrayIndStart = gridArray.GetIndStart();
    std::int8_t*** conditionArrayData = conditionArray.GetArrayData();

    std::size_t n0 = arrayShape[0];
    std::size_t n1 = arrayShape[1];
    std::size_t n2 = arrayShape[2];
    std::size_t ind0 = arrayIndStart[0];
    std::size_t ind1 = arrayIndStart[1];
    std::size_t ind2 = arrayIndStart[2];

    geometry->AreGridPointsInsideOrOn(gridArray, r0, dr, conditionArrayData);

    if(instruction == GAManipulatorInstructionCode::Equal) {
        for(std::size_t i0 = 0; i0 < n0; ++i0) {
            for(std::size_t i1 = 0; i1 < n1; ++i1) {
                for(std::size_t i2 = 0; i2 < n2; ++i2) {
                    if(conditionArrayData[i0][i1][i2] == 1) {
                        arrayData[ind0 + i0][ind1 + i1][ind2 + i2] = valueInside;
                    } else {
                        arrayData[ind0 + i0][ind1 + i1][ind2 + i2] = valueOutside;
                    }
                }
            }
        }
    } else if(instruction == GAManipulatorInstructionCode::MultiplyEqual) {
        for(std::size_t i0 = 0; i0 < n0; ++i0) {
            for(std::size_t i1 = 0; i1 < n1; ++i1) {
                for(std::size_t i2 = 0; i2 < n2; ++i2) {
                    if(conditionArrayData[i0][i1][i2] == 1) {
                        arrayData[ind0 + i0][ind1 + i1][ind2 + i2] *= valueInside;
                    } else {
                        arrayData[ind0 + i0][ind1 + i1][ind2 + i2] *= valueOutside;
                    }
                }
            }
        }
    } else {
            std::cout << "error: GAManipulatorInstructionCode operation not implemented." << std::endl;
            assert(false);
    }
}


