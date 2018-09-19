

#include <cassert>
#include <cstddef>      //std::size_t, nullptr
#include <array>        //std::array

#include "NumberTypes.h"
#include "YeeGridDataTypes.h"

/*
YeeGridData3D::YeeGridData3D() {}

YeeGridData3D::YeeGridData3D(const YeeGridData3D& gridData) :
        elemType(gridData.GetElemType()),
        numArray{gridData.GetNumArray(0), gridData.GetNumArray(1), gridData.GetNumArray(2)} { }

YeeGridData3D::YeeGridData3D(const YeeGridData3D&& gridData) :
        elemType(std::move(gridData.GetElemType())),
        numArray{std::move(gridData.GetNumArray(0)),
                 std::move(gridData.GetNumArray(1)),
                 std::move(gridData.GetNumArray(2))} { }

*/

YeeGridData3D::YeeGridData3D(ElementType elemType, std::array<std::size_t, 3>& nCells) : elemType(elemType) {
    if(elemType==ElementType::EdgeE) {
        numArray[0].ReInitialize({nCells[0]    , nCells[1] + 1, nCells[2] + 1}, 0.0);
        numArray[1].ReInitialize({nCells[0] + 1, nCells[1]    , nCells[2] + 1}, 0.0);
        numArray[2].ReInitialize({nCells[0] + 1, nCells[1] + 1, nCells[2]    }, 0.0);
    }else if(elemType==ElementType::EdgeH) {
        numArray[0].ReInitialize({nCells[0] + 1, nCells[1]    , nCells[2]    }, 0.0);
        numArray[0].ReInitialize({nCells[0]    , nCells[1] + 1, nCells[2]    }, 0.0);
        numArray[0].ReInitialize({nCells[0]    , nCells[1]    , nCells[2] + 1}, 0.0);
    }
}

ElementType& YeeGridData3D::GetElemType() {
    return elemType;
}

NumberArray3D<RealNumber>& YeeGridData3D::GetNumArray(int i) {
    return numArray[i];
}

