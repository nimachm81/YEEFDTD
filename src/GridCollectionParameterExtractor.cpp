

#include "GridCollectionParameterExtractor.h"


GridCollectionParameterExtractor::GridCollectionParameterExtractor(const std::string filename) :
        ParameterExtractor(filename) {
}

GridCollectionParameterExtractor::GridCollectionParameterExtractor(boost::property_tree::ptree ptree)  :
        ParameterExtractor(ptree) {
}

std::vector<std::tuple<std::size_t, std::size_t, std::vector<std::pair<std::string, std::string>>>>
GridCollectionParameterExtractor::GetRunSequence(const std::string path) {
    std::vector<std::tuple<std::size_t, std::size_t, std::vector<std::pair<std::string, std::string>>>> runSequence;
    ParameterExtractor runSequenceExtractor(GetSubTreeRootNode(path));
    for(std::size_t i = 0; i < runSequenceExtractor.GetSize(); ++i){
        ParameterExtractor runSequenceExtractor_i(runSequenceExtractor.GetSubTreeByIndex(i));
        std::tuple<std::size_t, std::size_t, std::vector<std::pair<std::string, std::string>>> indStart_indEnd_sequences(
                runSequenceExtractor_i.GetUintProperty("timeIndStart"),
                runSequenceExtractor_i.GetUintProperty("timeIndEnd"),
                runSequenceExtractor_i.GetStringPairArray("sequence")
                );
        runSequence.push_back(indStart_indEnd_sequences);
    }
    return runSequence;
}


