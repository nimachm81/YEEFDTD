
#ifndef FDTD_DISCRETEPOINTSGRIDARRAYMANIPULATOR_H_
#define FDTD_DISCRETEPOINTSGRIDARRAYMANIPULATOR_H_

#include <cassert>
#include <cmath>
#include <vector>

#include "GridArrayManipulator.h"

class DiscretePointsGridArrayManipulator : public GridArrayManipulator {
    public:
    virtual ~DiscretePointsGridArrayManipulator() { };

    void SetPositions(std::vector<std::array<FPNumber, 3>>* pos);
    void SetValues(std::vector<FPNumber>* vals);

    void UpdateArray(const FPNumber t, GAManipulatorInstructionCode instruction);
    FPNumber CalculateTime(const FPNumber dt, const std::size_t timeIndex);

    private:
    std::vector<std::array<FPNumber, 3>>* positions;
    std::vector<FPNumber>* values;
    FPNumber timeOffsetFraction = 0.0;

    int interpolationType;      // 0: closest,  1: linear

};


#endif // FDTD_DISCRETEPOINTSGRIDARRAYMANIPULATOR_H_


