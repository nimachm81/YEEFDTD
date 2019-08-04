
#include <cassert>

#include "YeeGridCollection.h"

YeeGridCollection::YeeGridCollection(std::size_t n_threads) : numOfThreads(n_threads)
                                                            {
    if(n_threads > 1) {
        useThreads = true;
        threadedAlgera = std::make_unique<ThreadedAlgebra>(n_threads);
    } else {
        useThreads = false;
    }
}

std::size_t YeeGridCollection::AddGrid() {
    grids.emplace_back();
    grids.back().SetThreadOptions(useThreads, threadedAlgera.get());
    return grids.size() - 1;
}

YeeGrid3D& YeeGridCollection::GetGrid(std::size_t index) {
    assert(index < grids.size());
    return grids[index];
}

void YeeGridCollection::RunInstructionsPeriodically(std::size_t timeIndStart, std::size_t timeIndEnd,
                                                    std::vector<std::pair<std::size_t, std::string>> instructions,
                                                    bool updateTimeIndexAutomatically
                                                    ) {
    for(std::size_t timeInd = timeIndStart; timeInd < timeIndEnd; ++timeInd) {
        for(auto& instruction : instructions) {
            YeeGrid3D& grid = grids[instruction.first];
            if(updateTimeIndexAutomatically) {
                grid.SetTimeIndex(timeInd);
                grid.ApplyInstructions(instruction.second, timeInd, timeInd + 1, /*writeToFile = */ false);
            } else {
                grid.ApplyInstructionsOnce(instruction.second);
            }
        }
        for(auto& grid : grids) {
            grid.WriteAllGridElemViewsToFile();             // TODO: if a grid is not involved in the iteration its time
                                                            // index remains zero, writing out the grid views at each time step
            //std::cout << "wrote to file." << std::endl;
        }
    }
}


