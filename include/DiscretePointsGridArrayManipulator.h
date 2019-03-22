
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

    void SetInterpolationType(int intpolType);
    void SetSameCellTreatmentType(int treatmentType);

    virtual void AddDataUpdater(DiscretePointsGAMDataUpdater* updater,
                        std::string dataName,   // name of data to attach
                        int direction           // which component (x, y or z) of the data to attach, in case data is defined as
                        ) = 0;                      // an array of vectors. For scalar data it will be ignored.

    virtual FPNumber GetDataValueByInddex(std::size_t ind) = 0;

    void UpdateArray(const FPNumber t, GAManipulatorInstructionCode instruction);
    FPNumber CalculateTime(const FPNumber dt, const std::size_t timeIndex);

    protected:
    std::vector<std::array<FPNumber, 3>>* positions;

    FPNumber timeOffsetFraction = 0.0;

    int interpolationType = 1;      // 0: closest,  1: linear
    int sameCellTreatmentType = 0;  // if two points fall in the same cell     0 : add them together
                                    //                                         1 : take the one with maximum absolute value

    DiscretePointsGAMDataUpdater* dataUpdater;
};


#endif // FDTD_DISCRETEPOINTSGRIDARRAYMANIPULATOR_H_



