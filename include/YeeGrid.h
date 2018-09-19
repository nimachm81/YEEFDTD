
#ifndef FDTD_YEEGRID_H_
#define FDTD_YEEGRID_H_


#include <cstddef>      //std::size_t
#include <string>
#include <unordered_map>
#include <memory>

#include "NumberTypes.h"
#include "YeeGridDataTypes.h"



class YeeGrid3D {
    public:
    YeeGrid3D(std::array<std::size_t, 3>& nCells);
    void AddGridElement(const std::string name, ElementType elemType);
    YeeGridData3D& GetGridElement(const std::string name);


    private:
    std::array<std::size_t, 3> nCells;
    std::unordered_map<std::string, std::unique_ptr<YeeGridData3D>> gridElements;

};


#endif  // FDTD_YEEGRID_H_


