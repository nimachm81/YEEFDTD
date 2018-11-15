
#ifndef FDTD_SINGLEGRIDPARAMETEREXTRACTOR_H_
#define FDTD_SINGLEGRIDPARAMETEREXTRACTOR_H_

#include <string>
#include <array>
#include <vector>
#include <tuple>

#include <boost/property_tree/ptree.hpp>
#include "NumberTypes.h"
#include "ParameterExtractor.h"

class SingleGridParameterExtractor : public ParameterExtractor {
    public:
    SingleGridParameterExtractor(const std::string filename);
    SingleGridParameterExtractor(boost::property_tree::ptree ptree);

    std::vector<std::array<std::string, 2>> GetEntireGridArrayParameters(const std::string path); // it takes a path to the entireArrayGrid names
                                                                // in the .json file and returns their names

    std::vector<std::tuple<std::string, std::string, std::array<std::size_t, 3>, std::array<std::size_t, 3>>>
            GetPartialGridArrayParameters(const std::string path);    // it takes a path to an array of partialGridArrays and returns
                                                                // a vector where each element contains
                                                                // arrayName, indStart, nCells

    std::vector<std::tuple<std::string, boost::property_tree::ptree>>   // it takes a path to a node with the following
            GetTypeStringAndParameterSubtree(const std::string path);   // structure
                                                                // [{"type":"some_type", "parameters":{some parameters}},
                                                                //  {"type":"some_other_type", "parameters":{some parameters}}
                                                                //  ... ]
                                                                // and returns a vector where each element contains the type
                                                                // and the ptree to its parameters

    std::vector<std::tuple<std::string, std::vector<std::string>>> GetUpdateSequences(const std::string path);

    std::vector<std::tuple<std::string, std::size_t, std::size_t>> GetRunSequence(const std::string path);
};

#endif  // FDTD_SINGLEGRIDPARAMETEREXTRACTOR_H_
