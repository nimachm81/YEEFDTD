
#ifndef FDTD_PARAMETEREXTRACTOR_H_
#define FDTD_PARAMETEREXTRACTOR_H_

#include <string>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include "NumberTypes.h"

class ParameterExtractor {
    public:
    ParameterExtractor(const std::string filename);
    std::string GetStringProperty(const std::string path);
    std::size_t GetUintProperty(const std::string path);
    RealNumber GetRealProperty(const std::string path);
    std::array<std::size_t, 3> Get3VecUintProperty(const std::string path);
    std::array<RealNumber, 3> Get3VecRealProperty(const std::string path);

    private:
    boost::property_tree::ptree treeRoot;
};


#endif // FDTD_PARAMETEREXTRACTOR_H_
