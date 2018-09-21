#ifndef FDTD_GRIDARRAYMANIPULATOR_H_
#define FDTD_GRIDARRAYMANIPULATOR_H_

#include <array>

#include "NumberTypes.h"
#include "YeeGridDataTypes.h"


class GridArrayManipulator {
    public:
    void SetNumberOfCells(std::array<std::size_t, 3>& nCells);
    void SetCornerCoordinates(std::array<RealNumber, 3>& r_0, std::array<RealNumber, 3>& r_1);
    virtual void UpdateArray(YeeGridData3D& gridData, RealNumber t) = 0;

    protected:
    std::array<RealNumber, 3> r_0;  // coordinates of the lower left corner
    std::array<RealNumber, 3> r_1;  // coordinates of the upper right corner
    std::array<std::size_t, 3> nCells;      // number of Yee cells

};

#endif // FDTD_GRIDARRAYMANIPULATOR_H_
