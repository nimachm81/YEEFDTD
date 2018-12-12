#ifndef FDTD_GRIDARRAYMANIPULATOR_H_
#define FDTD_GRIDARRAYMANIPULATOR_H_

#include <array>
#include <memory>

#include "NumberTypes.h"
#include "YeeGridDataTypes.h"


enum class GAManipulatorInstructionCode {
    Equal,
    PlusEqual,
    MinusEqual,
    MultiplyEqual,
    DevideEqual
};

class GridArrayManipulator {
    public:
    virtual ~GridArrayManipulator() { };
    void SetCornerCoordinate(std::array<FPNumber, 3>& r0);
    void SetGridSpacing(std::array<FPNumber, 3>& dr);
    void SetGridArrayTo(NumberArray3D<FPNumber>& gridData);
    void SetGridArrayTo(NumberArray3D<FPNumber>&& gridData);

    virtual void UpdateArray(const FPNumber t, GAManipulatorInstructionCode instruction) = 0;
    virtual FPNumber CalculateTime(const FPNumber dt, const std::size_t timeIndex) = 0;

    protected:
    std::array<FPNumber, 3> r0;  // coordinates of the first element [0,0,0] of gridArray
    std::array<FPNumber, 3> dr;  // distance between elements. If the array has only one component along a given direction
                                   // this distance represents the distance between 2 elemnts of the background grid in
                                   // the same direction.
    NumberArray3D<FPNumber> gridArray;     // array slice

};

#endif // FDTD_GRIDARRAYMANIPULATOR_H_
