
#ifndef FDTD_GRIDCOLLECTIONPARAMETEREXTRACTOR_H_
#define FDTD_GRIDCOLLECTIONPARAMETEREXTRACTOR_H_

#include <string>
#include <array>
#include <vector>
#include <tuple>

#include <boost/property_tree/ptree.hpp>
#include "NumberTypes.h"
#include "ParameterExtractor.h"

class GridCollectionParameterExtractor : public ParameterExtractor {
    public:
    GridCollectionParameterExtractor(const std::string filename);
    GridCollectionParameterExtractor(boost::property_tree::ptree ptree);

    std::vector<std::tuple<std::size_t, std::size_t, std::vector<std::pair<std::string, std::string>>>>
    GetRunSequence(const std::string path);

};

#endif  // FDTD_GRIDCOLLECTIONPARAMETEREXTRACTOR_H_
