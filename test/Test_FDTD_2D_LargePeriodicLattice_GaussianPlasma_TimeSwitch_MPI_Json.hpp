
#include <mpi.h>

#include "Test_FDTD_2D_LargePeriodicLattice_GaussianPlasma_TimeSwitch_Json.hpp"


int map_processor_to_parameter(int processRank, const int numOfParams,
                std::size_t* maxNumOfEachParameter,
                std::size_t* mappedParameterIndices,
                int processRankOffset = 0       // it is added to process rank
                ) {
    std::size_t cumulativeMaxNumOfParameters[numOfParams];
    cumulativeMaxNumOfParameters[0] = maxNumOfEachParameter[0];
    for(int i = 1; i < numOfParams; ++i) {
        cumulativeMaxNumOfParameters[i] = cumulativeMaxNumOfParameters[i - 1]*maxNumOfEachParameter[i];
    }

    for(int i = 0; i < numOfParams; ++i) {
        mappedParameterIndices[i] = 0;
    }

    int statusFlag = 0;
    int offsetedProcessRank = processRank + processRankOffset;

    if(offsetedProcessRank >= cumulativeMaxNumOfParameters[numOfParams - 1]) {
        offsetedProcessRank -= (offsetedProcessRank/cumulativeMaxNumOfParameters[numOfParams - 1])*
                                    cumulativeMaxNumOfParameters[numOfParams - 1];
        statusFlag = 1;   // processRank is outside expected range
    }

    for(int i = numOfParams - 1; i > 0; --i) {
        if(offsetedProcessRank >= cumulativeMaxNumOfParameters[i - 1]) {
            mappedParameterIndices[i] = offsetedProcessRank/cumulativeMaxNumOfParameters[i - 1];
            offsetedProcessRank -= mappedParameterIndices[i]*cumulativeMaxNumOfParameters[i - 1];
        }
    }
    mappedParameterIndices[0] = offsetedProcessRank;

    return statusFlag;
}

void test_run_fdtd_large_periodic_gaussian_plasma_time_switch_2d_from_json_parallel(int processRankOffset = 0) {

    int numOfProcesses;;
    int processRank;
    char processorName[MPI_MAX_PROCESSOR_NAME];
    int processorNameLen;

    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &numOfProcesses);
    MPI_Comm_rank(MPI_COMM_WORLD, &processRank);
    MPI_Get_processor_name(processorName, &processorNameLen);

    const std::size_t n_theta = 3;
    FPNumber theta_deg[n_theta] = {0.0, 22.5, 45.0};
    const std::size_t n_pitch_to_unitlength = 3;
    FPNumber pitch_to_unitlength[n_pitch_to_unitlength] = {0.2, 1.0, 2.0};
    const std::size_t n_fwhm_to_pitch = 2;
    FPNumber fwhm_to_pitch[n_fwhm_to_pitch] = {0.3, 0.6};
    const std::size_t n_wp_2p_thz = 2;
    FPNumber wp_2p_thz[n_wp_2p_thz] = {1.0, 1.5};
    const std::size_t n_gamma_thz = 1;
    FPNumber gamma_thz[n_gamma_thz] = {1.0};
    const std::size_t n_wp_switch_dt = 2;
    FPNumber wp_switch_dt[n_wp_switch_dt] = {0.1, 2.0};

    const int numOfParams = 6;
    std::size_t maxNumOfEachParameter[numOfParams] = {n_theta, n_pitch_to_unitlength, n_fwhm_to_pitch,
                                                      n_wp_2p_thz, n_gamma_thz, n_wp_switch_dt};
    std::size_t mappedParameterIndices[numOfParams];

    int statusFlag = map_processor_to_parameter(processRank, numOfParams,
                maxNumOfEachParameter,
                mappedParameterIndices,
                processRankOffset
                );

    if(statusFlag == 0) {
        FPNumber theta_deg_i = theta_deg[mappedParameterIndices[0]];
        FPNumber pitch_to_unitlength_i = pitch_to_unitlength[mappedParameterIndices[1]];
        FPNumber fwhm_to_pitch_i = fwhm_to_pitch[mappedParameterIndices[2]];
        FPNumber wp_2p_thz_i = wp_2p_thz[mappedParameterIndices[3]];
        FPNumber gamma_thz_i = gamma_thz[mappedParameterIndices[4]];
        FPNumber wp_switch_dt_i = wp_switch_dt[mappedParameterIndices[5]];

        std::cout << "Running process from processor " << processorName << ", rank: " << processRank
                  << "/" << numOfProcesses << " theta[" << mappedParameterIndices[0] << "]: " << theta_deg_i
                                           << " pitch_to_unitlength[" << mappedParameterIndices[1] << "]: " << pitch_to_unitlength_i
                                           << " fwhm_to_pitch[" << mappedParameterIndices[2] << "]: " << fwhm_to_pitch_i
                                           << " wp_2p_thz[" << mappedParameterIndices[3] << "]: " << wp_2p_thz_i
                                           << " gamma_thz[" << mappedParameterIndices[4] << "]: " << gamma_thz_i
                                           << " wp_switch_dt[" << mappedParameterIndices[5] << "]: " << wp_switch_dt_i
                                           << std::endl;

        test_run_fdtd_large_periodic_gaussian_plasma_time_switch_2d_from_json(
                        theta_deg_i,
                        pitch_to_unitlength_i,
                        fwhm_to_pitch_i,
                        wp_2p_thz_i,
                        gamma_thz_i,
                        wp_switch_dt_i
                        );
    } else {
        std::cout << "Error message from processor " << processorName << ", rank: " << processRank
                  << "/" << numOfProcesses << " status flag: " << statusFlag
                  << std::endl;

    }

    // Finalize the MPI environment.
    MPI_Finalize();
}


