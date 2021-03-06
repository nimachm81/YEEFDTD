
#ifndef FDTD_PARAMETEREXTRACTOR_H_
#define FDTD_PARAMETEREXTRACTOR_H_

#include <cstddef>
#include <string>
#include <iostream>
#include <unordered_map>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include "NumberTypes.h"
#include "DemotableComplex.hpp"

class ParameterExtractor {
    public:
    ParameterExtractor(const std::string filename);
    ParameterExtractor(boost::property_tree::ptree ptree);

    boost::property_tree::ptree& GetPropertyTreeRoot();
    std::size_t GetSize();
    std::size_t GetCount(const std::string key);

    boost::property_tree::ptree GetSubTreeRootNode(const std::string path);     // given a path to a branch node in
                                                                // a property tree, returns the branch node

    boost::property_tree::ptree GetSubTreeByIndex(const std::size_t index);     // given the array
                                                                // [subtree_0, subtree_1, ...]
                                                                // returns the subtree corresponding to the given index

    std::string GetStringProperty(const std::string path);      // gets a path to a key whose value is a string
                                                                // and returns the value as an string
                                                                // for exampe given a path to the element "key":"value"
                                                                // returns the "value"

    std::size_t GetUintProperty(const std::string path);        // given a path to "key":1162 retuens 1162 as size_t

    FPNumber GetRealProperty(const std::string path);         // given a path to "key":1.2 returns 1.2 the return type
                                                                // is FPNumber (which could have been defined to be of
                                                                // complex typme in the header NumberTypes.h)

    double GetDoubleProperty(const std::string path);           // given a path to "key":1.2 returns 1.2

    DemotableComplex<double> GetComplexProperty(const std::string path); // given a path to "key":[1.2, 3.1] returns
                                                                         // the complex number 1.2 + 3.1j

    std::array<std::size_t, 3> Get3VecUintProperty(const std::string path);     // given a path to
                                                                // "key":[1, 2, 3] returns the array [1, 2, 3]

    std::array<FPNumber, 3> Get3VecRealProperty(const std::string path);      // given a path to
                                                                // "key":[1.2, 3.1, 2.1] returns the array [1.2, 3.1, 2.1]

    std::array<FPNumber, 4> Get4VecRealProperty(const std::string path);      // given a path to
                                                                // "key":[1.2, 3.1, 2.1, 7.2] returns the array [1.2, 3.1, 2.1, 7.2]

    std::vector<std::string> GetStringArray(const std::string path);  // gets a path to a string array such as ["a", "few", "words"]
                                                                // and returns the strings in a vector

    std::pair<std::string, std::string> GetStringPairProperty(const std::string path); // gets a path to a string pair such
                                                                // as ["first", "second"] and returns the strings as a pair

    std::vector<std::pair<std::string, std::string>>
    GetStringPairArray(const std::string path);     // gets a path to a string pair such
                                            // as [["a", "b"], ["aa", "bb"], ...]  and returns them as vector<pair>

    std::vector<FPNumber> GetRealArray(const std::string path);     // gets a path to a real array such as [1.2, 2.3, 9.81 ...]
                                                                // and returns the numbers in a vector

    std::vector<std::array<std::size_t, 3>> Get3VecUintArray(const std::string path); // gets a path to an array
                                                                // such as [[1 , 2, 3], [2, 5, 1], [5, 6, 9], [6, 3, 0] ...]
                                                                // and returns the 3size_t-arrays in a vector

    std::vector<std::array<FPNumber, 3>> Get3VecRealArray(const std::string path);  // gets a path to an array
                                                                // such as [[1.1 , 2.2, 3.1], [2.4, 5.7, 1.5], [5.4, 6, 9.2], [6, 3.2, 0] ...]
                                                                // and returns the 3Real-arrays in a vector

    void static ReplaceStringsInFile(const std::string filename, const std::string filename_replaced,
            std::unordered_map<std::string, std::string>& str_replacewith);     // replaces the words as speified in
                                                                                // str_replacewith

    protected:
    boost::property_tree::ptree treeRoot;
};


#endif // FDTD_PARAMETEREXTRACTOR_H_
