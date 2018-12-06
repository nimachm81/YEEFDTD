

#include <cassert>
#include <cstddef>      //std::size_t, nullptr
#include <array>        //std::array

#include "NumberTypes.h"
#include "YeeGridDataTypes.h"


YeeGridData3D::YeeGridData3D(ElementType elemType, std::array<std::size_t, 3>& nCells,
                                                   std::array<std::size_t, 3> indOrigin) :
                                           elemType(elemType), indexOfOrigin(indOrigin), numOfCells(nCells) {
    if(elemType==ElementType::EdgeE) {
        numArray[0].ReInitialize({nCells[0]    , nCells[1] + 1, nCells[2] + 1}, 0.0);
        numArray[1].ReInitialize({nCells[0] + 1, nCells[1]    , nCells[2] + 1}, 0.0);
        numArray[2].ReInitialize({nCells[0] + 1, nCells[1] + 1, nCells[2]    }, 0.0);
    }else if(elemType==ElementType::EdgeH) {
        numArray[0].ReInitialize({nCells[0] + 1, nCells[1]    , nCells[2]    }, 0.0);
        numArray[1].ReInitialize({nCells[0]    , nCells[1] + 1, nCells[2]    }, 0.0);
        numArray[2].ReInitialize({nCells[0]    , nCells[1]    , nCells[2] + 1}, 0.0);
    }
}

ElementType& YeeGridData3D::GetElemType() {
    return elemType;
}

NumberArray3D<FPNumber>& YeeGridData3D::GetNumArray(int i) {
    assert(i>=0 && i<3);
    return numArray[i];
}

std::array<std::size_t, 3>& YeeGridData3D::GetIndexOfOrigin() {
    return indexOfOrigin;
}

std::ostream& operator<<(std::ostream& out, YeeGridData3D& gridData) {
    out << "x component : " << std::endl << gridData.GetNumArray(0);
    out << "y component : " << std::endl << gridData.GetNumArray(1);
    out << "z component : " << std::endl << gridData.GetNumArray(2);
    return out;
}
