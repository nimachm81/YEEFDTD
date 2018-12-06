

#include "GridArrayManipulator.h"



void GridArrayManipulator::SetCornerCoordinate(std::array<FPNumber, 3>& r0) {
    GridArrayManipulator::r0 = r0;
}

void GridArrayManipulator::SetGridSpacing(std::array<FPNumber, 3>& dr) {
    GridArrayManipulator::dr = dr;
}

void GridArrayManipulator::SetGridArrayTo(NumberArray3D<FPNumber>& gridArray) {
    GridArrayManipulator::gridArray.MakeThisASliceOf(gridArray);
}

void GridArrayManipulator::SetGridArrayTo(NumberArray3D<FPNumber>&& gridArray) {
    GridArrayManipulator::gridArray.MakeThisASliceOf(gridArray);
}

