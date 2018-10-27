
#include <cassert>

#include "YeeGridCollection.h"


std::size_t YeeGridCollection::AddGrid() {
    grids.emplace_back();
    return grids.size() - 1;
}


YeeGrid3D& YeeGridCollection::GetGrid(std::size_t index) {
    assert(index < grids.size());
    return grids[index];
}


void YeeGridCollection::RunInstructionsPeriodically(std::size_t timeIndStart, std::size_t timeIndEnd,
                                                    std::vector<std::string> instructions) {
    for(std::size_t timeInd = timeIndStart; timeInd < timeIndEnd; ++timeInd) {
        for(auto& instruction : instructions) {
            for(auto& grid : grids) {
                grid.SetTimeIndex(timeInd);
                grid.ApplyInstructions(instruction);
                grid.WriteAllGridElemViewsToFile();
            }
        }
    }
    for(auto& grid : grids) {
        grid.CloseGridViewFiles();
    }
}


