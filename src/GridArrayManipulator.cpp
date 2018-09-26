

#include "GridArrayManipulator.h"



void GridArrayManipulator::SetCornerCoordinates(std::array<RealNumber, 3>& r_0, std::array<RealNumber, 3>& r_1) {
    GridArrayManipulator::r_0 = r_0;
    GridArrayManipulator::r_1 = r_1;
}

void GridArrayManipulator::SetGridData(std::shared_ptr<YeeGridData3D>& gridData) {
    GridArrayManipulator::gridData = gridData;
}

