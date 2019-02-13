

#include "Test_FDTD_1D_Unidirectional_Source_Json.hpp"
#include "Test_FDTD_2D_LargePeriodicLattice_GaussianPlasma_Json.hpp"
#include "Test_FDTD_2D_LargePeriodicLattice_GaussianPlasma_MPI_Json.hpp"
#include "Test_FDTD_2D_LargePeriodicLattice_GaussianPlasma_TimeSwitch_Json.hpp"
#include "Test_FDTD_2D_LargePeriodicLattice_GaussianPlasma_TimeSwitch_MPI_Json.hpp"
#include "Test_FDTD_2D_LargePeriodicLattice_GaussianPlasma_TimeSwitch_GridCollection_Json.hpp"
#include "Test_FDTD_2D_SwitchedPlasma_CurvedLines_GridCollection_Json.hpp"
#include "Test_FDTD_2D_MonopoleCharge_Json.hpp"
#include "Test_FDTD_2D_MovingCharge_Json.hpp"
#include "Test_FDTD_2D_TM_Json.hpp"
#include "Test_FDTD_2D_PEC_Wedge_Json.hpp"
#include "Test_FDTD_2D_PEC_Wedge_electron_emitter_Json.hpp"

int main(int argc, char** argv) {
    //test_run_fdtd_2d_TM_from_json();
    //test_run_fdtd_curved_lines_gaussian_plasma_time_switch_2d_from_json();
    //test_run_fdtd_large_periodic_gaussian_plasma_time_switch_grid_collection_2d_from_json();
    //test_run_fdtd_large_periodic_gaussian_plasma_time_switch_2d_from_json();
    //test_run_fdtd_large_periodic_gaussian_plasma_time_switch_2d_from_json_parallel(0);

    //test_run_fdtd_large_periodic_gaussian_plasma_2d_from_json();
    //test_run_fdtd_1d_unidirectional_source_from_json();

    //test_run_fdtd_large_periodic_gaussian_plasma_2d_from_json_parallel();

    //test_run_fdtd_2d_moving_charge_from_json();
    //test_run_fdtd_2d_monopole_charge_from_json();

    //test_run_fdtd_2d_pec_wedge_from_json();
    test_run_fdtd_2d_pec_wedge_electron_emitter_from_json();
}


