
#ifndef FDTD_MULTIDIMARRAYALLOCATOR_H_
#define FDTD_MULTIDIMARRAYALLOCATOR_H_


#include <cstddef>      //std::size_t, nullptr
#include  <array>       //std::array

template <typename T>
T* Create1DNumberArray(std::size_t arrayShape, T initValue) {
    std::size_t n0 = arrayShape;
    T* tArray = new T[n0];
    for(std::size_t i0 = 0; i0 < n0; ++i0) {
        tArray[i0] = initValue;
    }
    return tArray;
};


template <typename T>
T* Create1DNumberArray(std::array<std::size_t, 1> arrayShape, T initValue) {
    std::size_t n0 = arrayShape[0];
    T* tArray = new T[n0];
    for(std::size_t i0 = 0; i0 < n0; ++i0) {
        tArray[i0] = initValue;
    }
    return tArray;
};

template <typename T>
void Delete1DNumberArray(T* tArray) {
    if(tArray != nullptr) {
        delete[] tArray;
        tArray = nullptr;
    }
};


template <typename T>
T** Create2DNumberArray(std::array<std::size_t, 2> arrayShape, T initValue) {
    std::size_t n0 = arrayShape[0];
    std::size_t n1 = arrayShape[1];
    T** tArray = new T*[n0];
    for(std::size_t i0 = 0; i0 < n0; ++i0) {
        tArray[i0] = new T[n1];
        for(std::size_t i1 = 0; i1 < n1; ++i1) {
            tArray[i0][i1] = initValue;
        }
    }
    return tArray;
};

template <typename T>
void Delete2DNumberArray(T** tArray, std::array<std::size_t, 3>& arrayShape) {
    if(tArray != nullptr) {
        std::size_t n0 = arrayShape[0];
        std::size_t n1 = arrayShape[1];
        for(std::size_t i0 = 0; i0 < n0; ++i0) {
            delete[] tArray[i0];
        }
        delete[] tArray;
        tArray = nullptr;
    }
};


template <typename T>
T*** Create3DNumberArray(std::array<std::size_t, 3>& arrayShape, T initValue) {
    std::size_t n0 = arrayShape[0];
    std::size_t n1 = arrayShape[1];
    std::size_t n2 = arrayShape[2];
    T*** tArray = new T**[n0];
    for(std::size_t i0 = 0; i0 < n0; ++i0) {
        tArray[i0] = new T*[n1];
        for(std::size_t i1 = 0; i1 < n1; ++i1) {
            tArray[i0][i1] = new T[n2];
            for(std::size_t i2 = 0; i2 < n2; ++i2) {
                tArray[i0][i1][i2] = initValue;
            }
        }
    }
    return tArray;
};


template <typename T>
void Delete3DNumberArray(T*** tArray, std::array<std::size_t, 3>& arrayShape) {
    if(tArray != nullptr) {
        std::size_t n0 = arrayShape[0];
        std::size_t n1 = arrayShape[1];
        std::size_t n2 = arrayShape[2];
        for(std::size_t i0 = 0; i0 < n0; ++i0) {
            for(std::size_t i1 = 0; i1 < n1; ++i1) {
                delete[] tArray[i0][i1];
            }
            delete[] tArray[i0];
        }
        delete[] tArray;
        tArray = nullptr;
    }
};


template <typename T>
T**** Create4DNumberArray(std::array<std::size_t, 4> arrayShape, T initValue) {
    std::size_t n0 = arrayShape[0];
    std::size_t n1 = arrayShape[1];
    std::size_t n2 = arrayShape[2];
    std::size_t n3 = arrayShape[3];
    T**** tArray = new T***[n0];
    for(std::size_t i0 = 0; i0 < n0; ++i0) {
        tArray[i0] = new T**[n1];
        for(std::size_t i1 = 0; i1 < n1; ++i1) {
            tArray[i0][i1] = new T*[n2];
            for(std::size_t i2 = 0; i2 < n2; ++i2) {
                tArray[i0][i1][i2] = new T*[n3];
                for(std::size_t i3 = 0; i3 < n3; ++i3) {
                    tArray[i0][i1][i2][i3] = initValue;
                }
            }
        }
    }
    return tArray;
};


template <typename T>
void Delete4DNumberArray(T**** tArray, std::array<std::size_t, 4> arrayShape) {
    if(tArray != nullptr) {
        std::size_t n0 = arrayShape[0];
        std::size_t n1 = arrayShape[1];
        std::size_t n2 = arrayShape[2];
        std::size_t n3 = arrayShape[3];
        for(std::size_t i0 = 0; i0 < n0; ++i0) {
            for(std::size_t i1 = 0; i1 < n1; ++i1) {
                for(std::size_t i2 = 0; i2 < n2; ++i2) {
                    delete[] tArray[i0][i1][i2];
                }
                delete[] tArray[i0][i1];
            }
            delete[] tArray[i0];
        }
        delete[] tArray;
        tArray = nullptr;
    }
};

#endif   // FDTD_MULTIDIMARRAYALLOCATOR_H_



