
#ifndef FDTD_YEEGRID_H_
#define FDTD_YEEGRID_H_


#include <cstddef>      //std::size_t

#include "NumberTypes.h"
#include "YeeGridDataTypes.h"


enum class ElementType {
  EdgeE = 1,
  EdgeH = 2,
  NodeE = 3,
  NodeH = 4
};


class YeeGrid3D {
    public:
    YeeGrid3D();
    //UpdateElement(elemType, updateInstructions);


    private:
    std::array<std::size_t, 3> nCells;

    std::vector<YeeGridScalar3D*> scalars;

    std::vector<YeeGridEtypeEdgeVectorArray3D*> eVectorArrays;
    std::vector<YeeGridHtypeEdgeVectorArray3D*> hVectorArrays;


};


#endif  // FDTD_YEEGRID_H_


