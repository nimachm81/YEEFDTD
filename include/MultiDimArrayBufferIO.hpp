#ifndef FDTD_MULTIDIMARRAYBUFFERIO_H_
#define FDTD_MULTIDIMARRAYBUFFERIO_H_

#include <cassert>
#include <fstream>

namespace ndarray {
namespace buffer {

template <typename T>
void Write3D(T* buffer, std::array<std::size_t, 3>& arrayShape, std::array<std::size_t, 3>& indStart,
             T*** arrayData) {
    std::size_t n0 = arrayShape[0];
    std::size_t n1 = arrayShape[1];
    std::size_t n2 = arrayShape[2];

    std::size_t ind0 = indStart[0];
    std::size_t ind1 = indStart[1];
    std::size_t ind2 = indStart[2];

    std::size_t indBuffer = 0;
    for(std::size_t i0 = 0; i0 < n0; ++i0) {
        for(std::size_t i1 = 0; i1 < n1; ++i1) {
            for(std::size_t i2 = 0; i2 < n2; ++i2) {
                buffer[indBuffer] = arrayData[ind0 + i0][ind1 + i1][ind2 + i2];
                ++indBuffer;
            }
        }
    }
};


template <typename T>
void Read3D(T* buffer, std::array<std::size_t, 3>& arrayShape, std::array<std::size_t, 3>& indStart,
            T*** arrayData) {
    std::size_t n0 = arrayShape[0];
    std::size_t n1 = arrayShape[1];
    std::size_t n2 = arrayShape[2];

    std::size_t ind0 = indStart[0];
    std::size_t ind1 = indStart[1];
    std::size_t ind2 = indStart[2];

    std::size_t indBuffer = 0;
    for(std::size_t i0 = 0; i0 < n0; ++i0) {
        for(std::size_t i1 = 0; i1 < n1; ++i1) {
            for(std::size_t i2 = 0; i2 < n2; ++i2) {
                arrayData[ind0 + i0][ind1 + i1][ind2 + i2] = buffer[indBuffer];
                ++indBuffer;
            }
        }
    }
};


} // namespace buffer
} // namespace ndarray


#endif // FDTD_MULTIDIMARRAYBUFFERIO_H_
