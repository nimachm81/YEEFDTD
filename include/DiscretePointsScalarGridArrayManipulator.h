
#ifndef FDTD_DISCRETEPOINTSSCALARGRIDARRAYMANIPULATOR_H_
#define FDTD_DISCRETEPOINTSSCALARGRIDARRAYMANIPULATOR_H_

#include <cassert>
#include <cmath>
#include <vector>
#include <string>


#include "DiscretePointsGridArrayManipulator.h"

class DiscretePointsScalarGridArrayManipulator : public DiscretePointsGridArrayManipulator {
    public:
    virtual ~DiscretePointsScalarGridArrayManipulator() { };

    virtual void AddDataUpdater(DiscretePointsGAMDataUpdater* updater,
                        std::string dataName,   // name of data to attach
                        int direction           // which component (x, y or z) of the data to attach, in case data is defined as
                        );                      // an array of vectors. For scalar data it will be ignored.

    virtual FPNumber GetDataValueByInddex(std::size_t ind);

    protected:
    std::vector<FPNumber>* values;
};


#endif // FDTD_DISCRETEPOINTSSCALARGRIDARRAYMANIPULATOR_H_


