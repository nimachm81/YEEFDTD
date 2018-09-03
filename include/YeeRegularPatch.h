
#ifndef FDTD_YEEREGULARPATCH_H_
#define FDTD_YEEREGULARPATCH_H_

#include <iostream>
#include <array>

#include "NumTypes.h"

namespace fdtd {

class YeeRegularPatch {
    public:
    void setStartIndices(std::array<std::size_t, NumDimensions>& indsStart);
    void setEndIndices(std::array<std::size_t, NumDimensions>& indsEnd);

    private:
    std::array<std::size_t, NumDimensions> indsStart;
    std::array<std::size_t, NumDimensions> indsEnd;
}

}

#endif  // FDTD_YEEREGULARPATCH_H_

