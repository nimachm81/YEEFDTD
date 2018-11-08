
#include "SingleGridParameterExtractor.h"

SingleGridParameterExtractor::SingleGridParameterExtractor(boost::property_tree::ptree ptree) :
        ParameterExtractor(ptree) {
}

std::vector<std::string> SingleGridParameterExtractor::GetEntireGridArrayNames(std::string path) {
    std::vector<std::string> gridArrayNames;
    ParameterExtractor entireGridArrayExtractor(GetSubTreeRootNode(path));
    for(std::size_t i = 0; i < entireGridArrayExtractor.GetSize(); ++i){
        ParameterExtractor entireGridArrayExtractor_i(entireGridArrayExtractor.GetSubTreeByIndex(i));
        gridArrayNames.emplace_back(entireGridArrayExtractor_i.GetStringProperty("name"));
    }
    return gridArrayNames;
}

std::vector<std::tuple<std::string, std::array<std::size_t, 3>, std::array<std::size_t, 3>>>
        SingleGridParameterExtractor::GetPartialGridArrayParameters(std::string path) {
    std::vector<std::tuple<std::string, std::array<std::size_t, 3>, std::array<std::size_t, 3>>> gridArrayParams;
    ParameterExtractor partialGridArrayExtractor(GetSubTreeRootNode(path));
    for(std::size_t i = 0; i < partialGridArrayExtractor.GetSize(); ++i){
        ParameterExtractor partialGridArrayExtractor_i(partialGridArrayExtractor.GetSubTreeByIndex(i));
        std::tuple<std::string, std::array<std::size_t, 3>, std::array<std::size_t, 3>> params(
                partialGridArrayExtractor_i.GetStringProperty("name"),
                partialGridArrayExtractor_i.Get3VecUintProperty("indStart"),
                partialGridArrayExtractor_i.Get3VecUintProperty("nCells"));
        gridArrayParams.push_back(params);
    }
    return gridArrayParams;
}

std::vector<std::tuple<std::string, boost::property_tree::ptree>>
        SingleGridParameterExtractor::GetTypeStringAndParameterSubtree(std::string path) {
    std::vector<std::tuple<std::string, boost::property_tree::ptree>> gridArrayManipulatorTypeAndParamNode;
    ParameterExtractor gridArrayManipulatorExtractor(GetSubTreeRootNode(path));
    for(std::size_t i = 0; i < gridArrayManipulatorExtractor.GetSize(); ++i){
        ParameterExtractor gridArrayManipulatorExtractor_i(gridArrayManipulatorExtractor.GetSubTreeByIndex(i));
        std::tuple<std::string, boost::property_tree::ptree> typeAndParamNode(
            gridArrayManipulatorExtractor_i.GetStringProperty("type"),
            gridArrayManipulatorExtractor.GetSubTreeRootNode("parameters"));
        gridArrayManipulatorTypeAndParamNode.push_back(typeAndParamNode);
    }
    return gridArrayManipulatorTypeAndParamNode;
}


