
#ifndef FDTD_PARAMETEREXTRACTOR_H_
#define FDTD_PARAMETEREXTRACTOR_H_

#include <cstddef>
#include <string>
#include <iostream>
#include <unordered_map>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include "NumberTypes.h"

class ParameterExtractor {
    public:
    ParameterExtractor(const std::string filename);
    ParameterExtractor(boost::property_tree::ptree ptree);

    std::size_t GetSize();
    std::size_t GetCount(const std::string key);

    boost::property_tree::ptree GetSubTreeRootNode(const std::string path);
    boost::property_tree::ptree GetSubTreeByIndex(const std::size_t index);

    std::string GetStringProperty(const std::string path);
    std::size_t GetUintProperty(const std::string path);
    RealNumber GetRealProperty(const std::string path);
    std::array<std::size_t, 3> Get3VecUintProperty(const std::string path);
    std::array<RealNumber, 3> Get3VecRealProperty(const std::string path);

    void static ReplaceStringsInFile(const std::string filename, const std::string filename_replaced,
            std::unordered_map<std::string, std::string>& str_replacewith);

    protected:
    boost::property_tree::ptree treeRoot;
};


#endif // FDTD_PARAMETEREXTRACTOR_H_
