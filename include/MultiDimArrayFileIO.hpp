#ifndef FDTD_MULTIDIMARRAYFILEIO_H_
#define FDTD_MULTIDIMARRAYFILEIO_H_

#include <cassert>
#include <fstream>


template <typename T>
void Write3DNumberArrayData(std::ofstream* fileOut,
                            std::array<std::size_t, 3>& arrayShape,
                            T*** arrayData,
                            int dataypeCode = -1,       // if code >=0 it will be written to file
                            bool writeShape = false,
                            bool writeDataTypeSize = false) {
    std::size_t n0 = arrayShape[0];
    std::size_t n1 = arrayShape[1];
    std::size_t n2 = arrayShape[2];

    std::size_t sizeOfT = sizeof(T);

    if(dataypeCode >= 0) {
        fileOut->write((char*)&dataypeCode, sizeof(int));
    }

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
void Write3DNumberArrayData(std::ofstream* fileOut,
                            std::array<std::size_t, 3>& arrayShape,
                            std::array<std::size_t, 3>& indStart,
                            T*** arrayData,
                            int dataypeCode = -1,       // if code >=0 it will be written to file
                            bool writeShape = false,
                            bool writeDataTypeSize = false) {
    std::size_t n0 = arrayShape[0];
    std::size_t n1 = arrayShape[1];
    std::size_t n2 = arrayShape[2];

    std::size_t ind0 = indStart[0];
    std::size_t ind1 = indStart[1];
    std::size_t ind2 = indStart[2];

    std::size_t sizeOfT = sizeof(T);

    if(dataypeCode >= 0) {
        fileOut->write((char*)&dataypeCode, sizeof(int));
    }

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
            fileOut->write((char*)(&arrayData[ind0 + i0][ind1 + i1][ind2]), sizeOfLastDimension);
        }
    }
};

inline void WriteToBuffer(char* buffer, std::size_t& bufferInd, char* charData, std::size_t charDataSize) {
    for(std::size_t i = 0; i < charDataSize; ++i) {
        buffer[bufferInd] = charData[i];
        ++bufferInd;
    }
}


//Warning MultiDimArray::GetMaxDataSize is based on the maximum size allocated for each array (including preamble)
// by the functions defined in MultiDimArrayFileIO
template <typename T>
void Write3DNumberArrayDataToMemory(char* buffer,
                            std::size_t& bufferInd,  //start writing to buffer from this index
                            std::array<std::size_t, 3>& arrayShape,
                            std::array<std::size_t, 3>& indStart,
                            T*** arrayData,
                            int dataypeCode = -1,       // if code >=0 it will be written to file
                            bool writeShape = false,
                            bool writeDataTypeSize = false) {
    std::size_t n0 = arrayShape[0];
    std::size_t n1 = arrayShape[1];
    std::size_t n2 = arrayShape[2];

    std::size_t ind0 = indStart[0];
    std::size_t ind1 = indStart[1];
    std::size_t ind2 = indStart[2];

    std::size_t sizeOfT = sizeof(T);

    if(dataypeCode >= 0) {
        WriteToBuffer(buffer, bufferInd, (char*)&dataypeCode, sizeof(int));
    }

    if(writeDataTypeSize) {
        WriteToBuffer(buffer, bufferInd, (char*)&sizeOfT, sizeof(std::size_t));
    }
    if(writeShape) {
        WriteToBuffer(buffer, bufferInd, (char*)&n0, sizeof(std::size_t));
        WriteToBuffer(buffer, bufferInd, (char*)&n1, sizeof(std::size_t));
        WriteToBuffer(buffer, bufferInd, (char*)&n2, sizeof(std::size_t));
    }

    std::size_t sizeOfLastDimension = n2*sizeof(T);
    for(std::size_t i0 = 0; i0 < n0; ++i0) {
        for(std::size_t i1 = 0; i1 < n1; ++i1) {
            WriteToBuffer(buffer, bufferInd, (char*)(&arrayData[ind0 + i0][ind1 + i1][ind2]), sizeOfLastDimension);
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
            fileIn->read((char*)(&arrayData[ind0 + i0][ind1 + i1][ind2]), sizeOfLastDimension);
        }
    }
};

#endif // FDTD_MULTIDIMARRAYFILEIO_H_

