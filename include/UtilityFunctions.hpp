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


#endif // FDTD_UTILITYFUNCTIONS

