
#ifndef FDTD_PERIODICGRIDARRAYMANIPULATOR_H_
#define FDTD_PERIODICGRIDARRAYMANIPULATOR_H_

#include <cassert>
#include <cmath>

#include "GridArrayManipulator.h"

class PeriodicGridArrayManipulator : public GridArrayManipulator {
    public:
    virtual ~PeriodicGridArrayManipulator() { };
    void SetUnitCellOrigin(const std::array<FPNumber, 3> origin);
    void SetPrimitiveVectors(const std::array<FPNumber, 3> v0,
                             const std::array<FPNumber, 3> v1,
                             const std::array<FPNumber, 3> v2);

    std::array<FPNumber, 3> BringPointInsideUnitCell(const std::array<FPNumber, 3>& r);

    void UpdateArray(const FPNumber t);

    virtual FPNumber Func(const std::array<FPNumber, 3> r, const FPNumber t) = 0;

    protected:
    std::array<FPNumber, 3> unitCellOrigin;       // the unit cell starts at this point and ends at the end of
    std::array<std::array<FPNumber, 3>, 3> primitiveVectors;     // the primitive vectors
};


#endif // FDTD_PERIODICGRIDARRAYMANIPULATOR_H_


