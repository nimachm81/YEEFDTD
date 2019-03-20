
#ifndef FDTD_VECTORFIELD_H_
#define FDTD_VECTORFIELD_H_

#include <array>
#include <vector>

#include "NumberTypes.h"

class VectorField {
    public:
    virtual ~VectorField() { };
    virtual void GetFieldValueAtPoint(FPNumber time,
                          std::array<FPNumber, 3>& position,
                          std::array<FPNumber, 3>& fieldValue
                          ) = 0;
    virtual void GetFieldValuesAtPoints(FPNumber time,
                          std::vector<std::array<FPNumber, 3>>& positions,
                          std::vector<std::array<FPNumber, 3>>& fieldValues
                          ) = 0;

};

#endif // FDTD_VECTORFIELD_H_

