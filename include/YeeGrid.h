
#ifndef FDTD_YEEGRID_H_
#define FDTD_YEEGRID_H_


#include "NumberTypes.h"

class YeeGrid {
    
    YeeGrid(int nDim);
    
    int nDim;
    std::vector<std::size_t> nPoints;
    
    std::vector<YeeGridScalar&> scalars;

    std::vector<YeeGridEScalarArray&> eScalarArrays; 
    std::vector<YeeGridHScalarArray&> hScalarArrays; 

    std::vector<YeeGridEVectorArray&> eVectorArrays; 
    std::vector<YeeGridHVectorArray&> hVectorArrays; 
    
    std::vector<YeeGridEMatrixArray&> eMatrixArrays; 
    std::vector<YeeGridHMatrixArray&> hMatrixArrays; 
    
}




#endif  // FDTD_YEEGRID_H_


