

#include "ParameterExtractor.h"


void test_read_json() {
    ParameterExtractor paramExtractor("instructions/MaxwellYee1D_processed.json");

    auto simulationType = paramExtractor.GetStringProperty("simulationType");

    ParameterExtractor dimensionsExtractor(paramExtractor.GetSubTreeRootNode("simulationParameters.dimensions"));
    auto r0 = dimensionsExtractor.Get3VecRealProperty("r0");
    auto r1 = dimensionsExtractor.Get3VecRealProperty("r1");
    auto nCells = dimensionsExtractor.Get3VecUintProperty("nCells");

    std::cout << "simulationType : " << simulationType << std::endl;
    std::cout << "r0 : " << r0[0] << " " << r0[1] << " " << r0[2] << std::endl;
    std::cout << "r1 : " << r1[0] << " " << r1[1] << " " << r1[2] << std::endl;
    std::cout << "nCells : " << nCells[0] << " " << nCells[1] << " " << nCells[2] << std::endl;

    ParameterExtractor entireGridArrayExtractor(paramExtractor.GetSubTreeRootNode("simulationParameters.entireGridArrays"));
    ParameterExtractor entireGridArrayExtractorFirst(entireGridArrayExtractor.GetSubTreeByIndex(1));
    std::string gridName = entireGridArrayExtractorFirst.GetStringProperty("name");
    std::cout << "GridArrayName : " << gridName << std::endl;

    ParameterExtractor updateSequenceExtractor(paramExtractor.GetSubTreeRootNode("simulationParameters.updateSequences"));
    for(std::size_t i = 0; i < updateSequenceExtractor.GetSize(); ++i) {
        ParameterExtractor updateSequenceExtractor_i(updateSequenceExtractor.GetSubTreeByIndex(i));
        auto update_sequence_name = updateSequenceExtractor_i.GetStringProperty("name");
        std::cout << "update sequence : " << update_sequence_name <<  std::endl;
        ParameterExtractor updateSequence_sequenceExtractor(updateSequenceExtractor_i.GetSubTreeRootNode("sequence"));
        for(std::size_t j = 0; j < updateSequence_sequenceExtractor.GetSize(); ++j) {
            ParameterExtractor updateSequence_sequenceExtractor_j(updateSequence_sequenceExtractor.GetSubTreeByIndex(j));
            auto update_name_j = updateSequence_sequenceExtractor_j.GetStringProperty("");
            std::cout << update_name_j <<  std::endl;
        }
    }
}


