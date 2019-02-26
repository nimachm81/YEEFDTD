
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
                                                    std::vector<std::pair<std::size_t, std::string>> instructions) {
    for(std::size_t timeInd = timeIndStart; timeInd < timeIndEnd; ++timeInd) {
        for(auto& instruction : instructions) {
            YeeGrid3D& grid = grids[instruction.first];
            grid.SetTimeIndex(timeInd);
            grid.ApplyInstructions(instruction.second, timeInd, timeInd + 1, /*writeToFile = */ false);
        }
        for(auto& grid : grids) {
            grid.WriteAllGridElemViewsToFile();             // TODO: if a grid is not involved in the iteration its time index remains zero, writing out the gread views at each time step
            //std::cout << "wrote to file." << std::endl;
        }
    }
}


