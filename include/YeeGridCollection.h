
#ifndef FDTD_YEEGRIDCOLLECTION_H_
#define FDTD_YEEGRIDCOLLECTION_H_

#include <cstddef>
#include <string>
#include <vector>

#include "NumberTypes.h"
#include "YeeGrid.h"
#include "ThreadedAlgebra.h"


class YeeGridCollection {
    public:
    YeeGridCollection(std::size_t n_threads = 1);
    std::size_t AddGrid();  // Warning: AddGrid() may resize the grids vector, making pointers to their elements invalid
                            // After taking any pointers or references to grids elements stop calling this function

    YeeGrid3D& GetGrid(std::size_t index);

    void RunInstructionsPeriodically(std::size_t timeIndStart, std::size_t timeIndEnd,
            std::vector<std::pair<std::size_t, std::string>> instructions,  // each instruction is composed of a grid index
                                                                            // and the name of the instruction to be invoked
                                                                            // by the grid
            bool updateTimeIndexAutomatically
            );


    protected:
    std::vector<YeeGrid3D> grids;

    bool useThreads = false;
    std::size_t numOfThreads = 1;
    std::unique_ptr<ThreadedAlgebra> threadedAlgera = nullptr;
};


#endif // FDTD_YEEGRIDCOLLECTION_H_


