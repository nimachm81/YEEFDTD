
#include "SingleGridParameterExtractor.h"

SingleGridParameterExtractor::SingleGridParameterExtractor(const std::string filename) :
        ParameterExtractor(filename) {
}

SingleGridParameterExtractor::SingleGridParameterExtractor(boost::property_tree::ptree ptree) :
        ParameterExtractor(ptree) {
}

std::vector<std::array<std::string, 2>>
            SingleGridParameterExtractor::GetEntireGridArrayParameters(const std::string path) {

    std::vector<std::array<std::string, 2>> gridArrayParams;
    ParameterExtractor entireGridArrayExtractor(GetSubTreeRootNode(path));
    for(std::size_t i = 0; i < entireGridArrayExtractor.GetSize(); ++i){
        ParameterExtractor entireGridArrayExtractor_i(entireGridArrayExtractor.GetSubTreeByIndex(i));
        std::array<std::string, 2> nameType{entireGridArrayExtractor_i.GetStringProperty("name"),
                                            entireGridArrayExtractor_i.GetStringProperty("type")};
        gridArrayParams.push_back(nameType);
    }
    return gridArrayParams;
}

std::vector<std::tuple<std::string, std::string, std::array<std::size_t, 3>, std::array<std::size_t, 3>>>
        SingleGridParameterExtractor::GetPartialGridArrayParameters(const std::string path) {

    std::vector<std::tuple<std::string, std::string, std::array<std::size_t, 3>, std::array<std::size_t, 3>>> gridArrayParams;
    ParameterExtractor partialGridArrayExtractor(GetSubTreeRootNode(path));
    for(std::size_t i = 0; i < partialGridArrayExtractor.GetSize(); ++i){
        ParameterExtractor partialGridArrayExtractor_i(partialGridArrayExtractor.GetSubTreeByIndex(i));
        std::tuple<std::string, std::string, std::array<std::size_t, 3>, std::array<std::size_t, 3>> params(
                partialGridArrayExtractor_i.GetStringProperty("name"),
                partialGridArrayExtractor_i.GetStringProperty("type"),
                partialGridArrayExtractor_i.Get3VecUintProperty("indStart"),
                partialGridArrayExtractor_i.Get3VecUintProperty("nCells"));
        gridArrayParams.push_back(params);
    }
    return gridArrayParams;
}

std::vector<std::tuple<std::string, boost::property_tree::ptree>>
        SingleGridParameterExtractor::GetTypeStringAndParameterSubtree(const std::string path) {
    std::vector<std::tuple<std::string, boost::property_tree::ptree>> gridArrayManipulatorTypeAndParamNode;
    ParameterExtractor gridArrayManipulatorExtractor(GetSubTreeRootNode(path));
    for(std::size_t i = 0; i < gridArrayManipulatorExtractor.GetSize(); ++i){
        ParameterExtractor gridArrayManipulatorExtractor_i(gridArrayManipulatorExtractor.GetSubTreeByIndex(i));
        std::tuple<std::string, boost::property_tree::ptree> typeAndParamNode(
            gridArrayManipulatorExtractor_i.GetStringProperty("type"),
            gridArrayManipulatorExtractor_i.GetSubTreeRootNode("parameters"));
        gridArrayManipulatorTypeAndParamNode.push_back(typeAndParamNode);
    }
    return gridArrayManipulatorTypeAndParamNode;
}

std::vector<std::tuple<std::string, std::vector<std::string>>>
                        SingleGridParameterExtractor::GetUpdateSequences(const std::string path) {
    std::vector<std::tuple<std::string, std::vector<std::string>>> updateSequences;
    ParameterExtractor updateSequenceExtractor(GetSubTreeRootNode(path));
    for(std::size_t i = 0; i < updateSequenceExtractor.GetSize(); ++i){
        ParameterExtractor updateSequenceExtractor_i(updateSequenceExtractor.GetSubTreeByIndex(i));
        std::tuple<std::string, std::vector<std::string>> updateSequence_i(
                updateSequenceExtractor_i.GetStringProperty("name"),
                updateSequenceExtractor_i.GetStringArray("sequence"));
        updateSequences.push_back(updateSequence_i);
    }
    return updateSequences;
}

std::vector<std::tuple<std::string, std::size_t, std::size_t>>
                        SingleGridParameterExtractor::GetRunSequence(const std::string path) {
    std::vector<std::tuple<std::string, std::size_t, std::size_t>> runSequence;
    ParameterExtractor runSequenceExtractor(GetSubTreeRootNode(path));
    for(std::size_t i = 0; i < runSequenceExtractor.GetSize(); ++i){
        ParameterExtractor runSequenceExtractor_i(runSequenceExtractor.GetSubTreeByIndex(i));
        std::tuple<std::string, std::size_t, std::size_t> name_indStart_indEnd(
                runSequenceExtractor_i.GetStringProperty("name"),
                runSequenceExtractor_i.GetUintProperty("timeIndStart"),
                runSequenceExtractor_i.GetUintProperty("timeIndEnd")
                );
        runSequence.push_back(name_indStart_indEnd);
    }
    return runSequence;
}


