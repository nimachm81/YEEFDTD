

#include "GridArrayManipulator.h"



void GridArrayManipulator::SetCornerCoordinates(std::array<RealNumber, 3>& r0, std::array<RealNumber, 3>& r1) {
    GridArrayManipulator::r0 = r0;
    GridArrayManipulator::r1 = r1;
}

void GridArrayManipulator::SetGridArrayTo(NumberArray3D<RealNumber>& gridArray) {
    GridArrayManipulator::gridArray.MakeThisASliceOf(gridArray);
}

