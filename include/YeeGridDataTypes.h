
#ifndef FDTD_YEEGRIDDATATYPES_H_
#define FDTD_YEEGRIDDATATYPES_H_

#include <cstddef>      //std::size_t, nullptr
#include  <array>       //std::array

#include "NumberTypes.h"
#include "MultiDimArray.hpp"


class YeeGridData3D {
    public:
    YeeGridData3D(const std::string name, std::array<std::size_t, 3>& nCells);
    
    protected:
    const std::string name;
    std::array<std::size_t, 3> nCells;
};

class YeeGridScalar3D : public YeeGridData3D {
    public:
    YeeGridScalar3D(const std::string name, std::array<std::size_t, 3>& nCells);
    void SetValue(RealNumber value);
    RealNumber GetValue() const;
    
    private:
    RealNumber value;

};

class YeeGridEVectorArray3D : public YeeGridData3D {
    public:
    YeeGridEVectorArray3D(const std::string name, std::array<std::size_t, 3>& nCells);
    
    private:
    NumberArray3D<RealNumber> numArray[3];
};

class YeeGridHVectorArray3D : public YeeGridData3D {
    public:
    YeeGridHVectorArray3D(const std::string name, std::array<std::size_t, 3>& nCells);
        
    private:
    NumberArray3D<RealNumber> numArray[3];
};


#endif  // FDTD_YEEGRIDDATATYPES_H_


