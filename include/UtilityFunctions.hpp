#ifndef FDTD_UTILITYFUNCTIONS
#define FDTD_UTILITYFUNCTIONS

#include <iostream>
#include <fstream>

#include "boost/lexical_cast.hpp"
#include "DemotableComplex.hpp"
#include "NumberTypes.h"


class  UtilityFunctions {
    public:
    static FPNumber FWHMtoDecayRate(FPNumber FWHM) {
        return (FPNumber)(2.0*std::sqrt(std::log(2)))/FWHM;
    };


    static std::string CastComplexToJSONString(FPNumber complxNum) {
        DemotableComplex<double> c(complxNum);
        std::string stringCast = std::string("[") +
                                boost::lexical_cast<std::string>(c.real()) +
                                std::string(",") +
                                boost::lexical_cast<std::string>(c.imag()) +
                                std::string("]");
        return stringCast;
    };

    template<typename T>
    static int GetDatatypeNumericalCode() {
        int dataTypeCode = -1;
        if(typeid(T) == typeid(float)) {
            dataTypeCode = 1;
        } else if(typeid(T) == typeid(double)) {
            dataTypeCode = 2;
        } else if(typeid(T) == typeid(std::complex<float>)) {
            dataTypeCode = 3;
        } else if(typeid(T) == typeid(std::complex<double>)) {
            dataTypeCode = 4;
        } else {
            std::cout << "error: unexpected data type." << std::endl;
            assert(false);
        }

        return dataTypeCode;
    };

    template<typename T>
    static char GetDatatypeCharCode() {
        char dtype;
        if(typeid(T)==typeid(float)) {
            dtype = 'f';
        } else if(typeid(T)==typeid(double)) {
            dtype = 'd';
        } else if(typeid(T)==typeid(std::size_t)) {
            dtype = 'u';
        } else if(typeid(T)==typeid(int)) {
            dtype = 'i';
        } else {
            std::cout << "error: unexpected data type." << std::endl;
            assert(false);
        }

        return dtype;
    };

    template<typename T>
    static int WriteParamToFile(std::ofstream& paramFileOut, T param, std::string paramName) {
        char dtype = GetDatatypeCharCode<T>();

        paramFileOut.write(&dtype, sizeof(dtype));
        paramFileOut.write((char*)&param, sizeof(T));
        const std::size_t nameLength = paramName.size();
        paramFileOut.write((char*)&nameLength, sizeof(std::size_t));
        const char* paramNameChars = paramName.c_str();
        paramFileOut.write((char*)paramNameChars, nameLength*sizeof(char));
    };

    template<typename T>
    static void Read1DArrayFromFile(std::ifstream& fileIn, std::vector<T>& dataOut) {
        fileIn.seekg(0, fileIn.end);
        std::size_t fileSize = fileIn.tellg();
        assert(fileSize % sizeof(T) == 0);

        std::size_t arraySize = fileSize / sizeof(T);
        dataOut.resize(arraySize);

        fileIn.seekg(0);
        fileIn.read((char*)(dataOut.data()), fileSize);
    };

    template<typename T>
    static void Read1DArrayFromFile(std::string fileName, std::vector<T>& dataOut) {
        std::ifstream file;
        file.open(fileName.c_str(), std::ios::in | std::ios::binary);
        assert(file.is_open());

        Read1DArrayFromFile(file, dataOut);
        file.close();
    }

    template<typename T>
    static void Write1DArrayToFile(std::ofstream& fileOut, std::vector<T>& dataOut) {
        std::size_t dataSize = dataOut.size() * sizeof(T);
        fileOut.write((char*)(dataOut.data()), dataSize);
    };

    template<typename T>
    static void Write1DArrayToFile(std::string fileName, std::vector<T>& dataOut) {
        std::ofstream file;
        file.open(fileName.c_str(), std::ios::out | std::ios::binary);
        assert(file.is_open());

        Write1DArrayToFile(file, dataOut);
        file.close();
    };

    static int map_processor_to_parameter(int processRank, const int numOfParams,
                    std::size_t* maxNumOfEachParameter,
                    std::size_t* mappedParameterIndices,
                    int processRankOffset = 0       // it is added to process rank
                    ) {
        std::size_t cumulativeMaxNumOfParameters[numOfParams];
        cumulativeMaxNumOfParameters[0] = maxNumOfEachParameter[0];
        for(int i = 1; i < numOfParams; ++i) {
            cumulativeMaxNumOfParameters[i] = cumulativeMaxNumOfParameters[i - 1]*maxNumOfEachParameter[i];
        }

        for(int i = 0; i < numOfParams; ++i) {
            mappedParameterIndices[i] = 0;
        }

        int statusFlag = 0;
        int offsetedProcessRank = processRank + processRankOffset;

        if(offsetedProcessRank >= cumulativeMaxNumOfParameters[numOfParams - 1]) {
            offsetedProcessRank -= (offsetedProcessRank/cumulativeMaxNumOfParameters[numOfParams - 1])*
                                        cumulativeMaxNumOfParameters[numOfParams - 1];
            statusFlag = 1;   // processRank is outside expected range
        }

        for(int i = numOfParams - 1; i > 0; --i) {
            if(offsetedProcessRank >= cumulativeMaxNumOfParameters[i - 1]) {
                mappedParameterIndices[i] = offsetedProcessRank/cumulativeMaxNumOfParameters[i - 1];
                offsetedProcessRank -= mappedParameterIndices[i]*cumulativeMaxNumOfParameters[i - 1];
            }
        }
        mappedParameterIndices[0] = offsetedProcessRank;

        return statusFlag;
    };

    // warning: no bound checks
    static void WriteToBuffer(char* buffer, std::size_t& bufferInd, char* charData, std::size_t charDataSize) {
        for(std::size_t i = 0; i < charDataSize; ++i) {
            buffer[bufferInd] = charData[i];
            ++bufferInd;
        }
    }

};

#endif // FDTD_UTILITYFUNCTIONS

