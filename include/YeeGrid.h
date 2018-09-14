
#ifndef FDTD_YEEGRID_H_
#define FDTD_YEEGRID_H_


#include <cstddef>      //std::size_t

#include "NumberTypes.h"
#include "YeeGridDataTypes.h"

class YeeGrid3D {
    public:
    YeeGrid3D();
    
    
    private:
    std::array<std::size_t, 3> nCells;
    
    std::vector<YeeGridScalar3D*> scalars;

    std::vector<YeeGridEVectorArray3D*> eVectorArrays; 
    std::vector<YeeGridHVectorArray3D*> hVectorArrays; 
    
        
};


#endif  // FDTD_YEEGRID_H_


