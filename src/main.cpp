

#include "Test_FDTD_1D_Unidirectional_Source_Json.hpp"
#include "Test_FDTD_2D_LargePeriodicLattice_GaussianPlasma_Json.hpp"
#include "Test_FDTD_2D_LargePeriodicLattice_GaussianPlasma_MPI_Json.hpp"

int main(int argc, char** argv) {
    //test_run_fdtd_large_periodic_gaussian_plasma_2d_from_json();
    //test_run_fdtd_1d_unidirectional_source_from_json();

    test_run_fdtd_large_periodic_gaussian_plasma_2d_from_json_parallel();
}


