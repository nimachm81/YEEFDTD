
#ifndef FDTD_YEEGRIDCOLLECTION_H_
#define FDTD_YEEGRIDCOLLECTION_H_

#include <cstddef>
#include <string>
#include <vector>

#include "NumberTypes.h"
#include "YeeGrid.h"

class YeeGridCollection {
    public:
    std::size_t AddGrid();
    YeeGrid3D& GetGrid(std::size_t index);

    void RunInstructionsPeriodically(std::size_t timeIndStart, std::size_t timeIndEnd,
            std::vector<std::pair<std::size_t, std::string>> instructions); // each instruction is composed of a grid index
                                                                            // and the name of the instruction to be invoked
                                                                            // by the grid
    protected:
    std::vector<YeeGrid3D> grids;
};


#endif // FDTD_YEEGRIDCOLLECTION_H_


