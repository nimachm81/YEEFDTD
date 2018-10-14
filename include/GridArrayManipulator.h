#ifndef FDTD_GRIDARRAYMANIPULATOR_H_
#define FDTD_GRIDARRAYMANIPULATOR_H_

#include <array>
#include <memory>

#include "NumberTypes.h"
#include "YeeGridDataTypes.h"


class GridArrayManipulator {
    public:
    virtual ~GridArrayManipulator() { };
    void SetCornerCoordinates(std::array<RealNumber, 3>& r0, std::array<RealNumber, 3>& r1);
    void SetGridArrayTo(NumberArray3D<RealNumber>& gridData);

    virtual void UpdateArray(const RealNumber t) = 0;
    virtual RealNumber CalculateTime(const RealNumber dt, const std::size_t timeIndex) = 0;

    protected:
    std::array<RealNumber, 3> r0;  // coordinates of the first element [0,0,0] of gridArray
    std::array<RealNumber, 3> r1;  // coordinates of the last element of gridArray
    NumberArray3D<RealNumber> gridArray;     // array slice

};

#endif // FDTD_GRIDARRAYMANIPULATOR_H_
