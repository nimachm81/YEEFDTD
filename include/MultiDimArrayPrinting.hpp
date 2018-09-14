#ifndef FDTD_MULTIDIMARRAYPRINTING_H_
#define FDTD_MULTIDIMARRAYPRINTING_H_


#include <cstddef>      //std::size_t, nullptr
#include <array>       //std::array
#include <iostream>


template <typename T>
void Print3DNumberArray(std::array<std::size_t, 3> arrayShape, T*** arrayData) {
    std::size_t n0 = arrayShape[0];
    std::size_t n1 = arrayShape[1];
    std::size_t n2 = arrayShape[2];
    std::cout << "[" ;
    for(std::size_t i0 = 0; i0 < n0; ++i0) {
        if(i0 > 0) std::cout << " ";
        std::cout << "[";
        for(std::size_t i1 = 0; i1 < n1; ++i1) {
            if(i1 > 0) std::cout << "  ";
            for(std::size_t i2 = 0; i2 < n2; ++i2) {
                if(i2 > 0) std::cout << " ";
                std::cout << arrayData[i0][i1][i2]; 
            }
            if(i1 < n1 - 1) std::cout << std::endl;
        }
        std::cout << "]";
        if(i0 < n0 - 1) std::cout << std::endl;
    }
    std::cout << "]" << std::endl;
};



#endif  // FDTD_MULTIDIMARRAYPRINTING_H_
