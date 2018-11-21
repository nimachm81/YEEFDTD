
#ifndef FDTD_PERIODICGRIDARRAYMANIPULATOR_H_
#define FDTD_PERIODICGRIDARRAYMANIPULATOR_H_

#include <cassert>
#include <cmath>

#include "GridArrayManipulator.h"

class PeriodicGridArrayManipulator : public GridArrayManipulator {
    public:
    virtual ~PeriodicGridArrayManipulator() { };
    void SetUnitCellOrigin(const std::array<RealNumber, 3> origin);
    void SetPrimitiveVectors(const std::array<RealNumber, 3> v0,
                             const std::array<RealNumber, 3> v1,
                             const std::array<RealNumber, 3> v2);

    std::array<RealNumber, 3> BringPointInsideUnitCell(const std::array<RealNumber, 3>& r);

    void UpdateArray(const RealNumber t);

    virtual RealNumber Func(const std::array<RealNumber, 3> r, const RealNumber t) = 0;

    protected:
    std::array<RealNumber, 3> unitCellOrigin;       // the unit cell starts at this point and ends at the end of
    std::array<std::array<RealNumber, 3>, 3> primitiveVectors;     // the primitive vectors
};


#endif // FDTD_PERIODICGRIDARRAYMANIPULATOR_H_


