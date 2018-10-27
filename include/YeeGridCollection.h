
#ifndef FDTD_YEEGRIDCOLLECTION_H_
#define FDTD_YEEGRIDCOLLECTION_H_

#include <cstddef>
#include <vector>

#include "YeeGrid.h"

class YeeGridCollection {
    public:
    std::size_t AddGrid();
    YeeGrid3D& GetGrid(std::size_t index);
    void RunInstructionsPeriodically(std::size_t timeIndStart, std::size_t timeIndEnd,
            std::vector<std::string> instructions   // apply these instructions to all grids periodically
            );

    private:
    std::vector<YeeGrid3D> grids;

};


#endif // FDTD_YEEGRIDCOLLECTION_H_


