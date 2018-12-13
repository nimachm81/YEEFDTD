
#include <mpi.h>

#include "Test_FDTD_2D_LargePeriodicLattice_GaussianPlasma_Json.hpp"


void test_run_fdtd_large_periodic_gaussian_plasma_2d_from_json_parallel() {

    int numOfProcesses;;
    int processRank;
    char processorName[MPI_MAX_PROCESSOR_NAME];
    int processorNameLen;

    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &numOfProcesses);
    MPI_Comm_rank(MPI_COMM_WORLD, &processRank);
    MPI_Get_processor_name(processorName, &processorNameLen);

    FPNumber theta_deg = 45.0*processRank/numOfProcesses;

    std::cout << "Running process from processor " << processorName << ", rank: " << processRank
              << "/" << numOfProcesses << " theta: " << theta_deg << std::endl;

    test_run_fdtd_large_periodic_gaussian_plasma_2d_from_json(theta_deg);

    // Finalize the MPI environment.
    MPI_Finalize();
}


