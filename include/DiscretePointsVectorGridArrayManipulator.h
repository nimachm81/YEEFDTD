
#ifndef FDTD_DISCRETEPOINTSVECTORGRIDARRAYMANIPULATOR_H_
#define FDTD_DISCRETEPOINTSVECTORGRIDARRAYMANIPULATOR_H_

#include <cassert>
#include <cmath>
#include <vector>
#include <string>


#include "DiscretePointsGridArrayManipulator.h"

class DiscretePointsVectorGridArrayManipulator : public DiscretePointsGridArrayManipulator {
    public:
    virtual ~DiscretePointsVectorGridArrayManipulator() { };

    virtual void AddDataUpdater(DiscretePointsGAMDataUpdater* updater,
                        std::string dataName,   // name of data to attach
                        int direction           // which component (x, y or z) of the data to attach, in case data is defined as
                        );                      // an array of vectors. For scalar data it will be ignored.

    virtual FPNumber GetDataValueByInddex(std::size_t ind);

    protected:
    std::vector<std::array<FPNumber, 3>>* values;
    int direction = 0;      // refers to x, y or z component of the values. Only this component will be used to update the grid..
};


#endif // FDTD_DISCRETEPOINTSVECTORGRIDARRAYMANIPULATOR_H_



