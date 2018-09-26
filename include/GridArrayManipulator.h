#ifndef FDTD_GRIDARRAYMANIPULATOR_H_
#define FDTD_GRIDARRAYMANIPULATOR_H_

#include <array>
#include <memory>

#include "NumberTypes.h"
#include "YeeGridDataTypes.h"


class GridArrayManipulator {
    public:
    virtual ~GridArrayManipulator() { };
    void SetCornerCoordinates(std::array<RealNumber, 3>& r_0, std::array<RealNumber, 3>& r_1);
    void SetGridData(std::shared_ptr<YeeGridData3D>& gridData);

    virtual void UpdateArray(const RealNumber t) = 0;

    protected:
    std::array<RealNumber, 3> r_0;  // coordinates of the lower left corner
    std::array<RealNumber, 3> r_1;  // coordinates of the upper right corner
    std::shared_ptr<YeeGridData3D> gridData;

};

#endif // FDTD_GRIDARRAYMANIPULATOR_H_
