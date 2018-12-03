#ifndef FDTD_GRIDARRAYMANIPULATOR_H_
#define FDTD_GRIDARRAYMANIPULATOR_H_

#include <array>
#include <memory>

#include "NumberTypes.h"
#include "YeeGridDataTypes.h"


class GridArrayManipulator {
    public:
    virtual ~GridArrayManipulator() { };
    void SetCornerCoordinate(std::array<RealNumber, 3>& r0);
    void SetGridSpacing(std::array<RealNumber, 3>& dr);
    void SetGridArrayTo(NumberArray3D<RealNumber>& gridData);
    void SetGridArrayTo(NumberArray3D<RealNumber>&& gridData);

    virtual void UpdateArray(const RealNumber t) = 0;
    virtual RealNumber CalculateTime(const RealNumber dt, const std::size_t timeIndex) = 0;

    protected:
    std::array<RealNumber, 3> r0;  // coordinates of the first element [0,0,0] of gridArray
    std::array<RealNumber, 3> dr;  // distance between elements. If the array has only one component along a given direction
                                   // this distance represents the distance between 2 elemnts of the background grid in
                                   // the same direction.
    NumberArray3D<RealNumber> gridArray;     // array slice

};

#endif // FDTD_GRIDARRAYMANIPULATOR_H_
