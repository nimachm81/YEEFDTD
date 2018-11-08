
#ifndef FDTD_SINGLEGRIDPARAMETEREXTRACTOR_H_
#define FDTD_SINGLEGRIDPARAMETEREXTRACTOR_H_

#include <string>
#include <array>
#include <vector>
#include <tuple>

#include <boost/property_tree/ptree.hpp>
#include "NumberTypes.h"
#include "ParameterExtractor.h"

class SingleGridParameterExtractor : ParameterExtractor {
    SingleGridParameterExtractor(boost::property_tree::ptree ptree);

    std::vector<std::string> GetEntireGridArrayNames(std::string path);
    std::vector<std::tuple<std::string, std::array<std::size_t, 3>, std::array<std::size_t, 3>>>
            GetPartialGridArrayParameters(std::string path);
    std::vector<std::tuple<std::string, boost::property_tree::ptree>>
            GetTypeStringAndParameterSubtree(std::string path);
};

#endif  // FDTD_SINGLEGRIDPARAMETEREXTRACTOR_H_
