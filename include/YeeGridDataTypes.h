
#ifndef FDTD_YEEGRIDDATATYPES_H_
#define FDTD_YEEGRIDDATATYPES_H_

#include <cstddef>      //std::size_t, nullptr
#include  <array>       //std::array

#include "NumberTypes.h"
#include "MultiDimArray.hpp"


enum class ElementType {
    ConstantScalar,
    EdgeE,
    EdgeH,
    NodeE,
    NodeH
};


class YeeGridData3D {
    public:
    //YeeGridData3D();
    //YeeGridData3D(const YeeGridData3D& gridData);
    //YeeGridData3D(const YeeGridData3D&& gridData);
    YeeGridData3D(ElementType elemType, std::array<std::size_t, 3>& nCells,
                                        std::array<std::size_t, 3> indOrigin={0,0,0});

    ElementType& GetElemType();
    NumberArray3D<RealNumber>& GetNumArray(int i);

    private:
    ElementType elemType;
    NumberArray3D<RealNumber> numArray[3];
    std::array<std::size_t, 3> indexOfOrigin;   // {0,0,0} starts on the lower left corner of the background grid
                                                // otherwise starts in the middle of the background grid
};

#endif  // FDTD_YEEGRIDDATATYPES_H_


