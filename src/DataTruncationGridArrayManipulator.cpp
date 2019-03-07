
#include "DataTruncationGridArrayManipulator.h"

void DataTruncationGridArrayManipulator::TruncateUp(bool truncate) {
    truncateUp = truncate;
}

void DataTruncationGridArrayManipulator::TruncateDown(bool truncate) {
    truncateDown = truncate;
}

void DataTruncationGridArrayManipulator::SetMaxValue(FPNumber value) {
    maxValue = value;
}

void DataTruncationGridArrayManipulator::SetMinValue(FPNumber value) {
    minValue = value;
}

FPNumber DataTruncationGridArrayManipulator::CalculateTime(const FPNumber dt, const std::size_t timeIndex) {
    return ((FPNumber)timeIndex + timeOffsetFraction) * dt;
}

void DataTruncationGridArrayManipulator::UpdateArray(const FPNumber t, GAManipulatorInstructionCode instruction) {
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

    if(truncateUp && truncateDown) {
        for(std::size_t i0 = 0; i0 < n0; ++i0) {
            for(std::size_t i1 = 0; i1 < n1; ++i1) {
                for(std::size_t i2 = 0; i2 < n2; ++i2) {
                    if(arrayData[ind0 + i0][ind1 + i1][ind2 + i2] > maxValue) {
                        arrayData[ind0 + i0][ind1 + i1][ind2 + i2] = maxValue;
                    }
                    if(arrayData[ind0 + i0][ind1 + i1][ind2 + i2] < minValue) {
                        arrayData[ind0 + i0][ind1 + i1][ind2 + i2] = minValue;
                    }
                }
            }
        }
    } else if(truncateUp) {
        for(std::size_t i0 = 0; i0 < n0; ++i0) {
            for(std::size_t i1 = 0; i1 < n1; ++i1) {
                for(std::size_t i2 = 0; i2 < n2; ++i2) {
                    if(arrayData[ind0 + i0][ind1 + i1][ind2 + i2] > maxValue) {
                        arrayData[ind0 + i0][ind1 + i1][ind2 + i2] = maxValue;
                    }
                }
            }
        }
    } else if(truncateDown) {
        for(std::size_t i0 = 0; i0 < n0; ++i0) {
            for(std::size_t i1 = 0; i1 < n1; ++i1) {
                for(std::size_t i2 = 0; i2 < n2; ++i2) {
                    if(arrayData[ind0 + i0][ind1 + i1][ind2 + i2] < minValue) {
                        arrayData[ind0 + i0][ind1 + i1][ind2 + i2] = minValue;
                    }
                }
            }
        }
    }
}




