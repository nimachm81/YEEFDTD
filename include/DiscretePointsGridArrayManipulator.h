
#ifndef FDTD_DISCRETEPOINTSGRIDARRAYMANIPULATOR_H_
#define FDTD_DISCRETEPOINTSGRIDARRAYMANIPULATOR_H_

#include <cassert>
#include <cmath>
#include <vector>
#include <string>


#include "GridArrayManipulator.h"
#include "DiscretePointsGAMDataUpdater.h"

class DiscretePointsGridArrayManipulator : public GridArrayManipulator {
    public:
    virtual ~DiscretePointsGridArrayManipulator() { };

    void AddDataUpdater(DiscretePointsGAMDataUpdater* updater,
                        std::string dataName,   // name of data to attach
                        int direction           // which component (x, y or z) of the data to attach
                        );

    void UpdateArray(const FPNumber t, GAManipulatorInstructionCode instruction);
    FPNumber CalculateTime(const FPNumber dt, const std::size_t timeIndex);

    protected:
    std::vector<std::array<FPNumber, 3>>* positions;
    std::vector<FPNumber>* values;
    FPNumber timeOffsetFraction = 0.0;

    int interpolationType = 1;      // 0: closest,  1: linear

    DiscretePointsGAMDataUpdater* dataUpdater;

};


#endif // FDTD_DISCRETEPOINTSGRIDARRAYMANIPULATOR_H_


