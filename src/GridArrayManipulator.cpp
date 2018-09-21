

#include "GridArrayManipulator.h"

void GridArrayManipulator::SetNumberOfCells(std::array<std::size_t, 3>& nCells) :
         nCells(yeeGrid.GetNumberOfCells()) { }


void GridArrayManipulator::SetCornerCoordinates(std::array<RealNumber, 3>& r_0, std::array<RealNumber, 3>& r_1) :
         r_0(yeeGrid.GetCornerR0())
        ,r_1(yeeGrid.GetCornerR1()) { }



