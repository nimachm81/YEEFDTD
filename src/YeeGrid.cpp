

#include <cstddef>      //std::size_t
#include <vector>
#include <memory>

#include "NumberTypes.h"
#include "YeeGrid.h"



YeeGrid3D::YeeGrid3D(std::array<std::size_t, 3>& nCells) : nCells(nCells) { }

void YeeGrid3D::AddGridElement(const std::string name, ElementType elemType) {
    gridElements[name] = std::make_unique<YeeGridData3D>(elemType, nCells);
}

YeeGridData3D& YeeGrid3D::GetGridElement(const std::string name) {
    return *gridElements[name];
}




