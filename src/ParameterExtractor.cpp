

#include "ParameterExtractor.h"


ParameterExtractor::ParameterExtractor(const std::string filename) {
    boost::property_tree::read_json(filename, treeRoot);
}

std::string ParameterExtractor::GetStringProperty(const std::string path) {
    return treeRoot.get<std::string>(path);
}

std::size_t ParameterExtractor::GetUintProperty(const std::string path) {
    return treeRoot.get<std::size_t>(path);
}

RealNumber ParameterExtractor::GetRealProperty(const std::string path) {
    return treeRoot.get<RealNumber>(path);
}

std::array<std::size_t, 3> ParameterExtractor::Get3VecUintProperty(const std::string path) {
    std::array<std::size_t, 3> vec;
    int i = 0;
    for(boost::property_tree::ptree::value_type &vec_i : treeRoot.get_child(path)) {
        vec[i] = vec_i.second.get_value<std::size_t>();
        ++i;
    }
    return vec;
}

std::array<RealNumber, 3> ParameterExtractor::Get3VecRealProperty(const std::string path) {
    std::array<RealNumber, 3> vec;
    int i = 0;
    for(boost::property_tree::ptree::value_type &vec_i : treeRoot.get_child(path)) {
        vec[i] = vec_i.second.get_value<RealNumber>();
        ++i;
    }
    return vec;
}



