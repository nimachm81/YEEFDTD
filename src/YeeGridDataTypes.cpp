

#include <cassert>
#include <cstddef>      //std::size_t, nullptr
#include <array>        //std::array

#include "NumberTypes.h"
#include "YeeGridDataTypes.h"

YeeGridData3D::YeeGridData3D(const std::string name, std::array<std::size_t, 3>& nCells) :
        name(name), nCells(nCells) {
    assert(nCells[0] > 0 && nCells[1] > 0 && nCells[2] > 0);
}

YeeGridScalar3D::YeeGridScalar3D(const std::string name, std::array<std::size_t, 3>& nCells) :
        YeeGridData3D(name, nCells) { }


void YeeGridScalar3D::SetValue(RealNumber value) {
    YeeGridScalar3D::value = value;
}

RealNumber YeeGridScalar3D::GetValue() const {
    return value;
}

YeeGridEtypeEdgeVectorArray3D::YeeGridEtypeEdgeVectorArray3D(
                    const std::string name, std::array<std::size_t, 3>& nCells) :
        YeeGridData3D(name, nCells),
        numArray{NumberArray3D<RealNumber>({nCells[0]    , nCells[1] + 1, nCells[2] + 1}, 0.0),
                 NumberArray3D<RealNumber>({nCells[0] + 1, nCells[1]    , nCells[2] + 1}, 0.0),
                 NumberArray3D<RealNumber>({nCells[0] + 1, nCells[1] + 1, nCells[2]    }, 0.0)} { }



YeeGridHtypeEdgeVectorArray3D::YeeGridHtypeEdgeVectorArray3D(
                    const std::string name, std::array<std::size_t, 3>& nCells) :
        YeeGridData3D(name, nCells),
        numArray{NumberArray3D<RealNumber>({nCells[0] + 1, nCells[1]    , nCells[2]    }, 0.0),
                 NumberArray3D<RealNumber>({nCells[0]    , nCells[1] + 1, nCells[2]    }, 0.0),
                 NumberArray3D<RealNumber>({nCells[0]    , nCells[1]    , nCells[2] + 1}, 0.0)} { }





