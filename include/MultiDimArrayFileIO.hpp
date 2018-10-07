#ifndef FDTD_MULTIDIMARRAYFILEIO_H_
#define FDTD_MULTIDIMARRAYFILEIO_H_

#include <cassert>
#include <fstream>


template <typename T>
void Write3DNumberArrayData(std::ofstream* fileOut, std::array<std::size_t, 3>& arrayShape,
                                      T*** arrayData,
                                      bool writeShape = false,
                                      bool writeDataTypeSize = false) {
    std::size_t n0 = arrayShape[0];
    std::size_t n1 = arrayShape[1];
    std::size_t n2 = arrayShape[2];

    std::size_t sizeOfT = sizeof(T);

    if(writeDataTypeSize) {
        fileOut->write((char*)&sizeOfT, sizeof(std::size_t));
    }
    if(writeShape) {
        fileOut->write((char*)&n0, sizeof(std::size_t));
        fileOut->write((char*)&n1, sizeof(std::size_t));
        fileOut->write((char*)&n2, sizeof(std::size_t));
    }

    std::streampos sizeOfLastDimension = n2*sizeof(T);
    for(std::size_t i0 = 0; i0 < n0; ++i0) {
        for(std::size_t i1 = 0; i1 < n1; ++i1) {
            fileOut->write((char*)(arrayData[i0][i1]), sizeOfLastDimension);
        }
    }
};

template <typename T>
void Write3DNumberArrayData(std::ofstream* fileOut, std::array<std::size_t, 3>& arrayShape,
                                      std::array<std::size_t, 3>& indStart,
                                      T*** arrayData,
                                      bool writeShape = false,
                                      bool writeDataTypeSize = false) {
    std::size_t n0 = arrayShape[0];
    std::size_t n1 = arrayShape[1];
    std::size_t n2 = arrayShape[2];

    std::size_t ind0 = indStart[0];
    std::size_t ind1 = indStart[1];
    std::size_t ind2 = indStart[2];

    std::size_t sizeOfT = sizeof(T);

    if(writeDataTypeSize) {
        fileOut->write((char*)&sizeOfT, sizeof(std::size_t));
    }
    if(writeShape) {
        fileOut->write((char*)&n0, sizeof(std::size_t));
        fileOut->write((char*)&n1, sizeof(std::size_t));
        fileOut->write((char*)&n2, sizeof(std::size_t));
    }

    std::streampos sizeOfLastDimension = n2*sizeof(T);
    for(std::size_t i0 = 0; i0 < n0; ++i0) {
        for(std::size_t i1 = 0; i1 < n1; ++i1) {
            fileOut->write((char*)(arrayData[ind0 + i0][ind1 + i1]), sizeOfLastDimension);
        }
    }
};

template <typename T>
void Read3DNumberArrayData(std::ifstream* fileIn, std::array<std::size_t, 3>& arrayShape,
                                                    T*** arrayData) {
    std::size_t n0 = arrayShape[0];
    std::size_t n1 = arrayShape[1];
    std::size_t n2 = arrayShape[2];

    std::streampos sizeOfLastDimension = n2*sizeof(T);

    for(std::size_t i0 = 0; i0 < n0; ++i0) {
        for(std::size_t i1 = 0; i1 < n1; ++i1) {
            fileIn->read((char*)(arrayData[i0][i1]), sizeOfLastDimension);
        }
    }
};

template <typename T>
void Read3DNumberArrayData(std::ifstream* fileIn, std::array<std::size_t, 3>& arrayShape,
                                     std::array<std::size_t, 3>& indStart,
                                     T*** arrayData) {
    std::size_t n0 = arrayShape[0];
    std::size_t n1 = arrayShape[1];
    std::size_t n2 = arrayShape[2];

    std::size_t ind0 = indStart[0];
    std::size_t ind1 = indStart[1];
    std::size_t ind2 = indStart[2];

    std::streampos sizeOfLastDimension = n2*sizeof(T);

    for(std::size_t i0 = 0; i0 < n0; ++i0) {
        for(std::size_t i1 = 0; i1 < n1; ++i1) {
            fileIn->read((char*)(arrayData[ind0 + i0][ind1 + i1]), sizeOfLastDimension);
        }
    }
};

#endif // FDTD_MULTIDIMARRAYFILEIO_H_

