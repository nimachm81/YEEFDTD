#ifndef FDTD_UTILITYFUNCTIONS
#define FDTD_UTILITYFUNCTIONS


#include "boost/lexical_cast.hpp"
#include "DemotableComplex.hpp"
#include "NumberTypes.h"

FPNumber FWHMtoDecayRate(FPNumber FWHM) {
    return (FPNumber)(2.0*std::sqrt(std::log(2)))/FWHM;
};


std::string CastComplexToJSONString(FPNumber complxNum) {
    DemotableComplex<double> c(complxNum);
    std::string stringCast = std::string("[") +
                            boost::lexical_cast<std::string>(c.real()) +
                            std::string(",") +
                            boost::lexical_cast<std::string>(c.imag()) +
                            std::string("]");
    return stringCast;
};

template<typename T>
int WriteParamToFile(std::ofstream& paramFileOut, T& param, std::string paramName) {
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
    }
    paramFileOut.write(&dtype, sizeof(dtype));
    paramFileOut.write((char*)&param, sizeof(T));
    const std::size_t nameLength = paramName.size();
    paramFileOut.write((char*)&nameLength, sizeof(std::size_t));
    const char* paramNameChars = paramName.c_str();
    paramFileOut.write((char*)paramNameChars, nameLength*sizeof(char));
};

int map_processor_to_parameter(int processRank, const int numOfParams,
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


#endif // FDTD_UTILITYFUNCTIONS

