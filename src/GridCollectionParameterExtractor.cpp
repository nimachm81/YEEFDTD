

#include "GridCollectionParameterExtractor.h"


GridCollectionParameterExtractor::GridCollectionParameterExtractor(const std::string filename) :
        ParameterExtractor(filename) {
}

GridCollectionParameterExtractor::GridCollectionParameterExtractor(boost::property_tree::ptree ptree)  :
        ParameterExtractor(ptree) {
}

std::vector<std::tuple<std::size_t, std::size_t, std::vector<std::pair<std::string, std::string>>, bool>>
GridCollectionParameterExtractor::GetRunSequence(const std::string path) {
    std::vector<std::tuple<std::size_t, std::size_t, std::vector<std::pair<std::string, std::string>>, bool>> runSequence;
    ParameterExtractor runSequenceExtractor(GetSubTreeRootNode(path));
    for(std::size_t i = 0; i < runSequenceExtractor.GetSize(); ++i){
        ParameterExtractor runSequenceExtractor_i(runSequenceExtractor.GetSubTreeByIndex(i));
        bool manualUpdate = false;
        if(runSequenceExtractor_i.GetCount("manualTimeUpdate") > 0) {
            if(runSequenceExtractor_i.GetStringProperty("manualTimeUpdate") == "yes") {
                manualUpdate = true;
            } else {
                assert(runSequenceExtractor_i.GetStringProperty("manualTimeUpdate") == "no");
            }
        }
        std::tuple<std::size_t, std::size_t, std::vector<std::pair<std::string, std::string>>, bool> indStart_indEnd_sequences_manual(
                runSequenceExtractor_i.GetUintProperty("timeIndStart"),
                runSequenceExtractor_i.GetUintProperty("timeIndEnd"),
                runSequenceExtractor_i.GetStringPairArray("sequence"),
                manualUpdate
                );
        runSequence.push_back(indStart_indEnd_sequences_manual);
    }
    return runSequence;
}


