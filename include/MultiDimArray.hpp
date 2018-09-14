
#ifndef FDTD_MULTIDIMARRAY_H_
#define FDTD_MULTIDIMARRAY_H_


#include <cassert>
#include <cstddef>      //std::size_t, nullptr
#include <array>       //std::array
#include <iostream>

#include "MultiDimArrayAllocator.hpp"
#include "MultiDimArrayPrinting.hpp"

template <typename T>
class NumberArray3D { 
    public:
    NumberArray3D(std::array<std::size_t, 3> shape, T initValue) : isSlice(false), shape(shape), indStart{0, 0, 0} {
        arrayData = Create3DNumberArray(shape, initValue);
    }
    
    NumberArray3D(T*** arraySlice, std::array<std::size_t, 3> shape, std::array<std::size_t, 3> indStart) : 
            shape(shape), indStart(indStart) {
        arrayData = arraySlice;
    }

    ~NumberArray3D() {
        if(!isSlice) {
            Delete3DNumberArray(arrayData, shape);
        }
    }
    
    T*** GetArrayData() const {
        return arrayData;
    }
    
    const std::array<std::size_t, 3>& GetIndStart() const {
        return indStart;
    }
    
    const std::array<std::size_t, 3>& GetShape() const {
        return shape;
    }

    NumberArray3D GetSlice(std::array<std::size_t, 3> indStart, std::array<std::size_t, 3> indEnd) {
        std::array<std::size_t, 3> shape;
        for(std::size_t i = 0; i < 3; ++i) {
            shape[i] = indEnd[i] - indStart[i];
        }
        return NumberArray3D(arrayData, shape, indStart);
    }
    
    //------------------------------ Operators ---------------------------------------
    
    NumberArray3D& operator+=(const NumberArray3D& rhs) {
        assert( rhs.GetShape() == shape ); 
        std::size_t n0 = shape[0];
        std::size_t n1 = shape[1];
        std::size_t n2 = shape[2];
        std::size_t ind0 = indStart[0];
        std::size_t ind1 = indStart[1];
        std::size_t ind2 = indStart[2];
        
        T*** rhs_arrayData = rhs.GetArrayData();
        const std::array<std::size_t, 3>& rhs_indStart = rhs.GetIndStart();
        std::size_t rhs_ind0 = rhs_indStart[0];
        std::size_t rhs_ind1 = rhs_indStart[1];
        std::size_t rhs_ind2 = rhs_indStart[2];
        
        for(std::size_t i0 = 0; i0 < n0; ++i0) {
            for(std::size_t i1 = 0; i1 < n1; ++i1) {
                for(std::size_t i2 = 0; i2 < n2; ++i2) {
                    arrayData[ind0 + i0][ind1 + i1][ind2 + i2] += 
                            rhs_arrayData[rhs_ind0 + i0][rhs_ind1 + i1][rhs_ind2 + i2];
                }
            }
        }
        return *this;
    }

    NumberArray3D& operator+=(const T rhs) {
        std::size_t n0 = shape[0];
        std::size_t n1 = shape[1];
        std::size_t n2 = shape[2];
        std::size_t ind0 = indStart[0];
        std::size_t ind1 = indStart[1];
        std::size_t ind2 = indStart[2];
        
        for(std::size_t i0 = 0; i0 < n0; ++i0) {
            for(std::size_t i1 = 0; i1 < n1; ++i1) {
                for(std::size_t i2 = 0; i2 < n2; ++i2) {
                    arrayData[ind0 + i0][ind1 + i1][ind2 + i2] += rhs;
                }
            }
        }
        return *this;
    }


    NumberArray3D& operator*=(const NumberArray3D& rhs) {
        assert( rhs.GetShape() == shape ); 
        std::size_t n0 = shape[0];
        std::size_t n1 = shape[1];
        std::size_t n2 = shape[2];
        std::size_t ind0 = indStart[0];
        std::size_t ind1 = indStart[1];
        std::size_t ind2 = indStart[2];
        
        T*** rhs_arrayData = rhs.GetArrayData();
        const std::array<std::size_t, 3>& rhs_indStart = rhs.GetIndStart();
        std::size_t rhs_ind0 = rhs_indStart[0];
        std::size_t rhs_ind1 = rhs_indStart[1];
        std::size_t rhs_ind2 = rhs_indStart[2];
        
        for(std::size_t i0 = 0; i0 < n0; ++i0) {
            for(std::size_t i1 = 0; i1 < n1; ++i1) {
                for(std::size_t i2 = 0; i2 < n2; ++i2) {
                    arrayData[ind0 + i0][ind1 + i1][ind2 + i2] *= 
                            rhs_arrayData[rhs_ind0 + i0][rhs_ind1 + i1][rhs_ind2 + i2];
                }
            }
        }
        return *this;
    }

    NumberArray3D& operator*=(const T rhs) {
        std::size_t n0 = shape[0];
        std::size_t n1 = shape[1];
        std::size_t n2 = shape[2];
        std::size_t ind0 = indStart[0];
        std::size_t ind1 = indStart[1];
        std::size_t ind2 = indStart[2];
        
        for(std::size_t i0 = 0; i0 < n0; ++i0) {
            for(std::size_t i1 = 0; i1 < n1; ++i1) {
                for(std::size_t i2 = 0; i2 < n2; ++i2) {
                    arrayData[ind0 + i0][ind1 + i1][ind2 + i2] *= rhs;
                }
            }
        }
        return *this;
    }
    
    NumberArray3D& operator/=(const NumberArray3D& rhs) {
        assert( rhs.GetShape() == shape ); 
        std::size_t n0 = shape[0];
        std::size_t n1 = shape[1];
        std::size_t n2 = shape[2];
        std::size_t ind0 = indStart[0];
        std::size_t ind1 = indStart[1];
        std::size_t ind2 = indStart[2];
        
        T*** rhs_arrayData = rhs.GetArrayData();
        const std::array<std::size_t, 3>& rhs_indStart = rhs.GetIndStart();
        std::size_t rhs_ind0 = rhs_indStart[0];
        std::size_t rhs_ind1 = rhs_indStart[1];
        std::size_t rhs_ind2 = rhs_indStart[2];
        
        for(std::size_t i0 = 0; i0 < n0; ++i0) {
            for(std::size_t i1 = 0; i1 < n1; ++i1) {
                for(std::size_t i2 = 0; i2 < n2; ++i2) {
                    arrayData[ind0 + i0][ind1 + i1][ind2 + i2] /= 
                            rhs_arrayData[rhs_ind0 + i0][rhs_ind1 + i1][rhs_ind2 + i2];
                }
            }
        }
        return *this;
    }

    NumberArray3D& operator/=(const T rhs) {
        std::size_t n0 = shape[0];
        std::size_t n1 = shape[1];
        std::size_t n2 = shape[2];
        std::size_t ind0 = indStart[0];
        std::size_t ind1 = indStart[1];
        std::size_t ind2 = indStart[2];
        
        for(std::size_t i0 = 0; i0 < n0; ++i0) {
            for(std::size_t i1 = 0; i1 < n1; ++i1) {
                for(std::size_t i2 = 0; i2 < n2; ++i2) {
                    arrayData[ind0 + i0][ind1 + i1][ind2 + i2] /= rhs;
                }
            }
        }
        return *this;
    }
    
    
    //------------------------------------------------------------------------------------
    
    void Print() {
        Print3DNumberArray(shape, arrayData);
    }
    
    private:
    const std::array<std::size_t, 3> indStart;
    const std::array<std::size_t, 3> shape;
    T*** arrayData = nullptr;
    const bool isSlice = true;
};


#endif  // FDTD_MULTIDIMARRAY_H_
