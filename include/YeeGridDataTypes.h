
#ifndef FDTD_YEEGRIDDATATYPES_H_
#define FDTD_YEEGRIDDATATYPES_H_

#include <cstddef>      //std::size_t, nullptr
#include  <array>       //std::array
#include <iostream>

#include "NumberTypes.h"
#include "MultiDimArray.hpp"
#include "ElementType.h"


class YeeGridData3D {
    public:
    YeeGridData3D(ElementType elemType, std::array<std::size_t, 3>& nCells,
                                        std::array<std::size_t, 3> indOrigin={0,0,0});

    ElementType& GetElemType();
    NumberArray3D<FPNumber>& GetNumArray(int i);
    std::array<std::size_t, 3>& GetIndexOfOrigin();

    private:
    ElementType elemType;
    NumberArray3D<FPNumber> numArray[3];
    std::array<std::size_t, 3> numOfCells;
    std::array<std::size_t, 3> indexOfOrigin;   // index of the first cell
                                                // {0,0,0} starts on the lower left corner of the background grid
                                                // otherwise starts in the middle of the background grid
};

std::ostream& operator<<(std::ostream& out, YeeGridData3D& gridData);

#endif  // FDTD_YEEGRIDDATATYPES_H_


