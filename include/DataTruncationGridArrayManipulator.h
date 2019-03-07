#ifndef FDTD_DATATRUNCATIONGRIDARRAYMANIPULATOR_H_
#define FDTD_DATATRUNCATIONGRIDARRAYMANIPULATOR_H_


#include <cassert>
#include <cmath>

#include "GridArrayManipulator.h"

class DataTruncationGridArrayManipulator : public GridArrayManipulator {
    public:
    virtual ~DataTruncationGridArrayManipulator() { };

    void TruncateUp(bool truncate);
    void TruncateDown(bool truncate);
    void SetMaxValue(FPNumber value);
    void SetMinValue(FPNumber value);

    void UpdateArray(const FPNumber t, GAManipulatorInstructionCode instruction);
    FPNumber CalculateTime(const FPNumber dt, const std::size_t timeIndex);

    private:
    bool truncateUp = true;     // truncate from above
    bool truncateDown = true;       // truncate from below

    FPNumber maxValue = 1.0;
    FPNumber minValue = 0.0;
    FPNumber timeOffsetFraction = 0.0;

    NumberArray3D<std::int8_t> conditionArray;     // 0: array point inside 1: array point outside
};

#endif // FDTD_DATATRUNCATIONGRIDARRAYMANIPULATOR_H_
