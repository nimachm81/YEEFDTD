
#ifndef TEST_FDTD_2D_LARGEPERIODICLATTICE_GAUSSIANPLASMA_TIMESWITCH_JSON
#define TEST_FDTD_2D_LARGEPERIODICLATTICE_GAUSSIANPLASMA_TIMESWITCH_JSON

#include "boost/lexical_cast.hpp"
#include "NumberTypes.h"
#include "ParameterExtractor.h"
#include "ParamFileTranslator.h"
#include "UtilityFunctions.hpp"

void test_run_fdtd_large_periodic_gaussian_plasma_time_switch_2d_from_json(FPNumber theta_deg = 0.0) {

    FPNumber pitch = 124.0;
    FPNumber FWHM = 2000.0;//54.0;
    FPNumber eps_r = 11.7;

    std::size_t numOfSamplesPerUnitLength = 100;

    FPNumber y0 = -5.0;
    FPNumber y1 = 5.0;
    FPNumber z0 = -6.0;
    FPNumber z1 = 6.0;
    std::size_t ny = static_cast<std::size_t>(std::real(y1 - y0) * numOfSamplesPerUnitLength);;
    std::size_t nz = static_cast<std::size_t>(std::real(z1 - z0) * numOfSamplesPerUnitLength);
    FPNumber dy = (y1 - y0)/(FPNumber)(ny);
    FPNumber dz = (z1 - z0)/(FPNumber)(nz);
    FPNumber stabilityFactor = 0.95;
    FPNumber dt = (FPNumber)1.0/std::sqrt((FPNumber)1.0/(dy*dy) + (FPNumber)1.0/(dz*dz))*stabilityFactor;
    FPNumber z_j = -2.0;
    std::size_t indzJ = std::round(std::real((z_j - z0)/dz));
    FPNumber j_center_y = (y0 + y1)/(FPNumber)2.0;
    FPNumber j_decay_rate_y = 0.5;
    FPNumber j_center_t = 2.5;
    FPNumber j_decay_rate_t = 0.8;
    FPNumber jm_center_t = j_center_t + dt/(FPNumber)2.0*std::sqrt(eps_r);
    FPNumber jm_amplitude = std::sqrt(eps_r);
    FPNumber j_mod_freq = 0.1;
    FPNumber j_mod_phase = M_PI/2.0;

    FPNumber gamma = 1.0e12/(3.0e8/(pitch*1.0e-6));
    FPNumber wp = 1.0e12*(2.0*M_PI)/(3.0e8/(pitch*1.0e-6));
    std::cout << "wp : " << wp << std::endl;

    FPNumber wp2_decayrate_y = FWHMtoDecayRate(FWHM/pitch);
    FPNumber wp2_decayrate_z = FWHMtoDecayRate(FWHM/pitch);
    FPNumber wp2_center_y = 0.0;
    FPNumber wp2_center_z = 0.0;

    FPNumber wp2_mask_z0 = -6.0;
    FPNumber wp2_mask_t0 = -z_j*std::sqrt(eps_r) + j_center_t;
    FPNumber wp2_mask_t1 = 100.0;
    FPNumber wp2_mask_dt = 0.1;

    std::size_t numOfTimeSamplesBeforeSwitch = static_cast<std::size_t>((wp2_mask_t0 - wp2_mask_dt)/dt);
    std::size_t numOfTimeSamplesDuringSwitch = static_cast<std::size_t>(2.0*wp2_mask_dt/dt);
    std::size_t numOfTimeSamplesAfterSwitch = 2000;

    std::size_t nt_0 = numOfTimeSamplesBeforeSwitch;
    std::size_t nt_1 = nt_0 + numOfTimeSamplesDuringSwitch;
    std::size_t nt_2 = nt_1 + numOfTimeSamplesAfterSwitch;

    std::cout << "Nt : " << nt_2 << std::endl;

    FPNumber ei_z_record = -1.5;
    FPNumber er_z_record = -2.5;
    FPNumber et_z_record = +2.5;
    FPNumber wp2_z_record = 0;
    std::size_t ei_indz_record = std::round(std::real((ei_z_record - z0)/dz));
    std::size_t er_indz_record = std::round(std::real((er_z_record - z0)/dz));
    std::size_t et_indz_record = std::round(std::real((et_z_record - z0)/dz));
    std::size_t wp2_indz_record = std::round(std::real((wp2_z_record - z0)/dz));

    FPNumber theta = theta_deg/180.0*M_PI;
    FPNumber a1_y = std::cos(theta);
    FPNumber a1_z = -std::sin(theta);
    FPNumber a2_y = std::sin(theta);
    FPNumber a2_z = std::cos(theta);

    std::string e_slice_i_name = std::string("\"") + "2D/Ei-x-slice-" + boost::lexical_cast<std::string>(std::real(theta_deg)) +
                               std::string("\"");
    std::string e_slice_r_name = std::string("\"") + "2D/Er-x-slice-" + boost::lexical_cast<std::string>(std::real(theta_deg)) +
                               std::string("\"");
    std::string e_slice_t_name = std::string("\"") + "2D/Et-x-slice-" + boost::lexical_cast<std::string>(std::real(theta_deg)) +
                               std::string("\"");
    std::string e_name = std::string("\"") + "2D/E-x-" + boost::lexical_cast<std::string>(std::real(theta_deg)) +
                               std::string("\"");
    std::string wp2_name = std::string("\"") + "2D/Wp2-x-" + boost::lexical_cast<std::string>(std::real(theta_deg)) +
                               std::string("\"");
    std::string wp2_slice_name = std::string("\"") + "2D/Wp2-slice-x-" + boost::lexical_cast<std::string>(std::real(theta_deg)) +
                               std::string("\"");


    std::unordered_map<std::string, std::string> str_replacewith{
            {"\"_y0_\"", boost::lexical_cast<std::string>(std::real(y0))},
            {"\"_y1_\"", boost::lexical_cast<std::string>(std::real(y1))},
            {"\"_z0_\"", boost::lexical_cast<std::string>(std::real(z0))},
            {"\"_z1_\"", boost::lexical_cast<std::string>(std::real(z1))},
            {"\"_ny_\"", boost::lexical_cast<std::string>(ny)},
            {"\"_nz_\"", boost::lexical_cast<std::string>(nz)},
            {"\"_ny_p1_\"", boost::lexical_cast<std::string>(ny + 1)},
            {"\"_nz_p1_\"", boost::lexical_cast<std::string>(nz + 1)},
            {"\"_ny_m1_\"", boost::lexical_cast<std::string>(ny - 1)},
            {"\"_nz_m1_\"", boost::lexical_cast<std::string>(nz - 1)},
            {"\"_dy_\"", boost::lexical_cast<std::string>(std::real(dy))},
            {"\"_dz_\"", boost::lexical_cast<std::string>(std::real(dz))},
            {"\"_dt_\"", boost::lexical_cast<std::string>(std::real(dt))},
            {"\"_m_dt_\"", boost::lexical_cast<std::string>(std::real(-dt))},
            {"\"_dt_dy_\"", boost::lexical_cast<std::string>(std::real(dt/dy))},
            {"\"_m_dt_dy_\"", boost::lexical_cast<std::string>(std::real(-dt/dy))},
            {"\"_dt_dz_\"", boost::lexical_cast<std::string>(std::real(dt/dz))},
            {"\"_m_dt_dz_\"", boost::lexical_cast<std::string>(std::real(-dt/dz))},
            {"\"_dt_dy_eps_\"", boost::lexical_cast<std::string>(std::real(dt/dy/eps_r))},
            {"\"_m_dt_dy_eps_\"", boost::lexical_cast<std::string>(std::real(-dt/dy/eps_r))},
            {"\"_dt_dz_eps_\"", boost::lexical_cast<std::string>(std::real(dt/dz/eps_r))},
            {"\"_m_dt_dz_eps_\"", boost::lexical_cast<std::string>(std::real(-dt/dz/eps_r))},
            {"\"_indzJ_\"", boost::lexical_cast<std::string>(indzJ)},
            {"\"_indzJ_p1_\"", boost::lexical_cast<std::string>(indzJ + 1)},
            {"\"_indzJ_m1_\"", boost::lexical_cast<std::string>(indzJ - 1)},
            {"\"_J_center_y_\"", boost::lexical_cast<std::string>(std::real(j_center_y))},
            {"\"_J_decayrate_y_\"", boost::lexical_cast<std::string>(std::real(j_decay_rate_y))},
            {"\"_J_center_t_\"", boost::lexical_cast<std::string>(std::real(j_center_t))},
            {"\"_Jm_center_t_\"", boost::lexical_cast<std::string>(std::real(jm_center_t))},
            {"\"_J_decayrate_t_\"", boost::lexical_cast<std::string>(std::real(j_decay_rate_t))},
            {"\"_Jm_amp_\"", boost::lexical_cast<std::string>(std::real(jm_amplitude))},
            {"\"_J_mod_freq_\"", boost::lexical_cast<std::string>(std::real(j_mod_freq))},
            {"\"_J_mod_phase_\"", boost::lexical_cast<std::string>(std::real(j_mod_phase))},
            {"\"_m_dt_tau_\"", boost::lexical_cast<std::string>(std::real(-dt*gamma))},
            {"\"_wp_sq_\"", boost::lexical_cast<std::string>(std::real(wp*wp))},
            {"\"_wp_sq_center_y_\"", boost::lexical_cast<std::string>(std::real(wp2_center_y))},
            {"\"_wp_sq_center_z_\"", boost::lexical_cast<std::string>(std::real(wp2_center_z))},
            {"\"_wp_sq_decayrate_y_\"", boost::lexical_cast<std::string>(std::real(wp2_decayrate_y))},
            {"\"_wp_sq_decayrate_z_\"", boost::lexical_cast<std::string>(std::real(wp2_decayrate_z))},
            {"\"_cube_y0_\"", boost::lexical_cast<std::string>(std::real(y0))},
            {"\"_cube_y1_\"", boost::lexical_cast<std::string>(std::real(y1))},
            {"\"_cube_z0_\"", boost::lexical_cast<std::string>(std::real(wp2_mask_z0))},
            {"\"_cube_z1_\"", boost::lexical_cast<std::string>(std::real(z1))},
            {"\"_cube_t0_\"", boost::lexical_cast<std::string>(std::real(wp2_mask_t0))},
            {"\"_cube_t1_\"", boost::lexical_cast<std::string>(std::real(wp2_mask_t1))},
            {"\"_cube_dy_\"", boost::lexical_cast<std::string>(std::real(0.0))},
            {"\"_cube_dz_\"", boost::lexical_cast<std::string>(std::real(0.1))},
            {"\"_cube_dt_\"", boost::lexical_cast<std::string>(std::real(wp2_mask_dt))},
            {"\"_e_ind_z_record_i_\"", boost::lexical_cast<std::string>(ei_indz_record)},
            {"\"_e_ind_z_record_i_p1_\"", boost::lexical_cast<std::string>(ei_indz_record + 1)},
            {"\"_e_ind_z_record_r_\"", boost::lexical_cast<std::string>(er_indz_record)},
            {"\"_e_ind_z_record_r_p1_\"", boost::lexical_cast<std::string>(er_indz_record + 1)},
            {"\"_e_ind_z_record_t_\"", boost::lexical_cast<std::string>(et_indz_record)},
            {"\"_e_ind_z_record_t_p1_\"", boost::lexical_cast<std::string>(et_indz_record + 1)},
            {"\"_wp2_ind_z_record_\"", boost::lexical_cast<std::string>(wp2_indz_record)},
            {"\"_wp2_ind_z_record_p1_\"", boost::lexical_cast<std::string>(wp2_indz_record + 1)},
            {"\"_nt_0_\"", boost::lexical_cast<std::string>(nt_0)},
            {"\"_nt_1_\"", boost::lexical_cast<std::string>(nt_1)},
            {"\"_nt_2_\"", boost::lexical_cast<std::string>(nt_2)},
            {"\"_a1_y_\"", boost::lexical_cast<std::string>(std::real(a1_y))},
            {"\"_a1_z_\"", boost::lexical_cast<std::string>(std::real(a1_z))},
            {"\"_a2_y_\"", boost::lexical_cast<std::string>(std::real(a2_y))},
            {"\"_a2_z_\"", boost::lexical_cast<std::string>(std::real(a2_z))},
            {"\"_Eslice_name_i_\"", e_slice_i_name},
            {"\"_Eslice_name_r_\"", e_slice_r_name},
            {"\"_Eslice_name_t_\"", e_slice_t_name},
            {"\"_E_name_\"", e_name},
            {"\"_Wp2Slice_name_\"", wp2_slice_name},
            {"\"_Wp2_name_\"", wp2_name}
            };

    // for parallel processing this file's name should be unique to the processor. Here the suffix theta is supposes to
    // make it unique (each process is assumed treating a different theta)
    std::string processedFileName =
            std::string("instructions/processed/MaxwellYee2D_LargePeriodicGaussianPlasma_TimeSwitch_processed_") +
            boost::lexical_cast<std::string>(std::real(theta_deg)) + ".json";

    ParameterExtractor::ReplaceStringsInFile("instructions/MaxwellYee2D_LargePeriodicGaussianPlasma_TimeSwitch.json",
                processedFileName, str_replacewith);

    ParamFileTranslator fileTranslator(processedFileName);
    fileTranslator.Translate();
};

#endif // TEST_FDTD_2D_LARGEPERIODICLATTICE_GAUSSIANPLASMA_TIMESWITCH_JSON

