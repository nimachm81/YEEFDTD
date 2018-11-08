
#include <cassert>

#include "ParameterExtractor.h"


ParameterExtractor::ParameterExtractor(const std::string filename) {
    boost::property_tree::read_json(filename, treeRoot);
}

ParameterExtractor::ParameterExtractor(boost::property_tree::ptree ptree) {
    treeRoot = ptree;
}

boost::property_tree::ptree ParameterExtractor::GetSubTreeRootNode(const std::string path) {
    return treeRoot.get_child(path);
}

boost::property_tree::ptree ParameterExtractor::GetSubTreeByIndex(const std::size_t index) {
    //boost::property_tree::ptree::value_type& it = treeRoot.front();
    //return it.second;//.get_value<boost::property_tree::ptree>();
    assert(index < GetSize());
    auto it = treeRoot.begin();
    std::advance(it, index);
    return it->second;
}

std::size_t ParameterExtractor::GetSize() {
    return treeRoot.size();
}

std::size_t ParameterExtractor::GetCount(const std::string key) {
    return treeRoot.count(key);
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


void ParameterExtractor::ReplaceStringsInFile(const std::string filename, const std::string filename_replaced,
        std::unordered_map<std::string, std::string>& str_replacewith) {
    std::ifstream infile;
    infile.open(filename.c_str(), std::ios::in);
    assert(infile.is_open());

    std::ofstream outfile;
    outfile.open(filename_replaced.c_str(), std::ios::out);
    assert(outfile.is_open());

    std::string line;
    while (std::getline(infile, line)) {
        for(auto& it : str_replacewith) {
            std::size_t pos = 0;
            while(true) {
                pos = line.find(it.first, pos);
                if (pos != std::string::npos) {
                    line.replace(pos, it.first.size(), it.second);
                } else {
                    break;
                }
            }
        }
        outfile << line << std::endl;
    }
    infile.close();
    outfile.close();
}

