
#ifndef TEST_FDTD_2D_LARGEPERIODICLATTICE_GAUSSIANPLASMA_TIMESWITCH_JSON
#define TEST_FDTD_2D_LARGEPERIODICLATTICE_GAUSSIANPLASMA_TIMESWITCH_JSON

#include<typeinfo>

#include "boost/lexical_cast.hpp"
#include "NumberTypes.h"
#include "ParameterExtractor.h"
#include "ParamFileTranslator.h"
#include "UtilityFunctions.hpp"


void test_run_fdtd_large_periodic_gaussian_plasma_time_switch_2d_from_json(
                        FPNumber theta_deg = 0.0,   // roration angle of the periodic plasma
                        FPNumber pitch_to_unitlength = 1.0,  // periodicity in natural units
                        FPNumber fwhm_to_pitch = 0.3,    // fwhm/pitch
                        FPNumber wp_2p_thz = 1.0,   // wp/2pi
                        FPNumber gamma_thz = 1.0,   // scattering rate
                        FPNumber wp_switch_dt = 0.1,          // switch time in natural units
                        FPNumber celldisplacement_to_pitch = 0.5     // the structure is shifted by celldisplacement_to_pitch*pitch_to_unitlength
                        ) {
    FPNumber unit_length_si = 124.0e-6;     // unit length in SI
    FPNumber pitch_um = pitch_to_unitlength*unit_length_si;

    FPNumber eps_r = 11.7;      // silicon

    std::size_t numOfSamplesPerUnitLength = 50;

    FPNumber y0 = -10.0;
    FPNumber y1 = 10.0;
    FPNumber z0 = -12.0;
    FPNumber z1 = 12.0;
    std::size_t ny = static_cast<std::size_t>(std::real(y1 - y0) * numOfSamplesPerUnitLength);
    std::size_t nz = static_cast<std::size_t>(std::real(z1 - z0) * numOfSamplesPerUnitLength);
    FPNumber dy = (y1 - y0)/(FPNumber)(ny);
    FPNumber dz = (z1 - z0)/(FPNumber)(nz);
    FPNumber stabilityFactor = 0.95;
    FPNumber dt = (FPNumber)1.0/std::sqrt((FPNumber)1.0/(dy*dy) + (FPNumber)1.0/(dz*dz))*stabilityFactor;
    FPNumber z_j = -2.0;
    std::size_t indzJ = std::round(std::real((z_j - z0)/dz));
    FPNumber j_center_y = (y0 + y1)/(FPNumber)2.0;
    FPNumber j_decay_rate_y = 0.5;
    FPNumber j_center_t = 3.5;
    FPNumber j_decay_rate_t = 0.8;
    FPNumber jm_center_t = j_center_t + dt/(FPNumber)2.0*std::sqrt(eps_r);
    FPNumber jm_amplitude = std::sqrt(eps_r);
    FPNumber j_mod_freq = 0.1;
    FPNumber j_mod_phase = M_PI/2.0;

    FPNumber gamma = gamma_thz*(FPNumber)1.0e12/((FPNumber)3.0e8/unit_length_si);
    FPNumber wp = wp_2p_thz*(FPNumber)1.0e12
                  *((FPNumber)2.0*(FPNumber)M_PI)
                  /((FPNumber)3.0e8/unit_length_si);
    std::cout << "wp : " << wp << std::endl;

    FPNumber wp2_decayrate_y = FWHMtoDecayRate(fwhm_to_pitch*pitch_to_unitlength);
    FPNumber wp2_decayrate_z = FWHMtoDecayRate(fwhm_to_pitch*pitch_to_unitlength);
    FPNumber wp2_center_y = 0.0;
    FPNumber wp2_center_z = celldisplacement_to_pitch*pitch_to_unitlength;

    FPNumber wp2_mask_z0 = -5.0;
    FPNumber wp2_mask_z1 = +5.0;
    FPNumber wp2_mask_t0 = -z_j*std::sqrt(eps_r) + j_center_t;
    FPNumber wp2_mask_t1 = 1000.0;
    FPNumber wp2_mask_dt = wp_switch_dt;

    std::size_t numOfTimeSamplesBeforeSwitch = static_cast<std::size_t>(std::real((wp2_mask_t0 - wp2_mask_dt)/dt));
    std::size_t numOfTimeSamplesDuringSwitch = static_cast<std::size_t>(std::real((FPNumber)2.0*wp2_mask_dt/dt));
    std::size_t numOfTimeSamplesAfterSwitch = static_cast<std::size_t>(std::real(
            numOfTimeSamplesBeforeSwitch*(4.0 + std::real(wp)*0.2)
            ));
    if(std::real(fwhm_to_pitch) <= 1.0) {
        numOfTimeSamplesAfterSwitch = static_cast<std::size_t>(std::real(
            (FPNumber)numOfTimeSamplesBeforeSwitch*((FPNumber)4.0 + std::real(wp)*(FPNumber)0.2*fwhm_to_pitch)
            ));
    }

    std::size_t nt_0 = numOfTimeSamplesBeforeSwitch;
    std::size_t nt_1 = nt_0 + numOfTimeSamplesDuringSwitch;
    std::size_t nt_2 = nt_1 + numOfTimeSamplesAfterSwitch;

    std::cout << "Nt : " << nt_2 << std::endl;

    FPNumber ei_z_record = -1.5;
    FPNumber er_z_record = -5.5;
    FPNumber et_z_record = +5.5;
    FPNumber wp2_z_record = 0;
    std::size_t ei_indz_record = std::round(std::real((ei_z_record - z0)/dz));
    std::size_t er_indz_record = std::round(std::real((er_z_record - z0)/dz));
    std::size_t et_indz_record = std::round(std::real((et_z_record - z0)/dz));
    std::size_t wp2_indz_record = std::round(std::real((wp2_z_record - z0)/dz));

    FPNumber field_entire_z_m_corner_record = -6.0;
    FPNumber field_entire_z_p_corner_record = +6.0;
    FPNumber field_entire_y_m_corner_record = -6.0;
    FPNumber field_entire_y_p_corner_record = +6.0;
    std::size_t field_entire_ind_z_m_record = std::round(std::real((field_entire_z_m_corner_record - z0)/dz));
    std::size_t field_entire_ind_z_p_record = std::round(std::real((field_entire_z_p_corner_record - z0)/dz));
    std::size_t field_entire_ind_y_m_record = std::round(std::real((field_entire_y_m_corner_record - y0)/dy));
    std::size_t field_entire_ind_y_p_record = std::round(std::real((field_entire_y_p_corner_record - y0)/dy));


    FPNumber theta = theta_deg/(FPNumber)180.0*(FPNumber)M_PI;
    FPNumber a1_y = pitch_to_unitlength*std::cos(theta);
    FPNumber a1_z = -pitch_to_unitlength*std::sin(theta);
    FPNumber a2_y = pitch_to_unitlength*std::sin(theta);
    FPNumber a2_z = pitch_to_unitlength*std::cos(theta);

    FPNumber wp_switch_dt_si = wp_switch_dt/((FPNumber)3.0e8/unit_length_si);

    std::string outputFolder = "LargePeriodicLattice-GaussianPlasma-TimeSwitched/";

    std::string file_suffix = std::string("-rot=") + boost::lexical_cast<std::string>(std::real(theta_deg)) +
                 "-fp=" + boost::lexical_cast<std::string>(std::real(wp_2p_thz)) +
                 "-gamma=" + boost::lexical_cast<std::string>(std::real(gamma_thz)) +
                 "-pitch=" + boost::lexical_cast<std::string>(std::real(pitch_um*(FPNumber)1.0e6)) +
                 "-fwhmToPitch=" + boost::lexical_cast<std::string>(std::real(fwhm_to_pitch)) +
                 "-swithTime=" + boost::lexical_cast<std::string>(std::real(wp_switch_dt)) +
                 "-displacementFactor=" + boost::lexical_cast<std::string>(std::real(celldisplacement_to_pitch)) +
                 "-res=" + boost::lexical_cast<std::string>(numOfSamplesPerUnitLength)
                 ;

    std::string e_slice_i_name = std::string("\"") + outputFolder + "Ei-x-slice" +
                                 file_suffix +
                                 std::string("\"");
    std::string e_slice_r_name = std::string("\"") + outputFolder + "Er-x-slice" +
                                 file_suffix +
                                 std::string("\"");
    std::string e_slice_t_name = std::string("\"") + outputFolder + "Et-x-slice" +
                                 file_suffix +
                                 std::string("\"");
    std::string e_name = std::string("\"") + outputFolder + "E-x" +
                                 file_suffix +
                                 std::string("\"");
    std::string wp2_name = std::string("\"") + outputFolder + "Wp2-x" +
                                 file_suffix +
                                 std::string("\"");
    std::string wp2_slice_name = std::string("\"") + outputFolder + "Wp2-slice-x" +
                                 file_suffix +
                                 std::string("\"");

    std::size_t EorWp_entire_sample_rate = 50;

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
            {"\"_cube_z1_\"", boost::lexical_cast<std::string>(std::real(wp2_mask_z1))},
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
            {"\"_Wp2_name_\"", wp2_name},
            {"\"_sample_rate_\"", boost::lexical_cast<std::string>(EorWp_entire_sample_rate)},
            {"\"_e_entire_ind_ym_record_\"", boost::lexical_cast<std::string>(field_entire_ind_y_m_record)},
            {"\"_e_entire_ind_yp_record_\"", boost::lexical_cast<std::string>(field_entire_ind_y_p_record)},
            {"\"_e_entire_ind_zm_record_\"", boost::lexical_cast<std::string>(field_entire_ind_z_m_record)},
            {"\"_e_entire_ind_zp_record_\"", boost::lexical_cast<std::string>(field_entire_ind_z_p_record)}
            };

    // for parallel processing this file's name should be unique to the processor. Here the suffix theta is supposes to
    // make it unique (each process is assumed treating a different theta)
    std::string processedFileName =
            std::string("instructions/processed/MaxwellYee2D_LargePeriodicGaussianPlasma_TimeSwitch_processed_") +
            file_suffix + ".json";

    ParameterExtractor::ReplaceStringsInFile("instructions/MaxwellYee2D_LargePeriodicGaussianPlasma_TimeSwitch.json",
                processedFileName, str_replacewith);

    ParamFileTranslator fileTranslator(processedFileName);
    fileTranslator.Translate();

    std::string parametersFileName =
            std::string("data/LargePeriodicLattice-GaussianPlasma-TimeSwitched/params") +
            file_suffix +
            ".param";
    std::ofstream paramFileOut(parametersFileName.c_str(), std::ios::out | std::ios::binary);
    WriteParamToFile<FPNumber>(paramFileOut, dt, "dt");   // 0
    WriteParamToFile<FPNumber>(paramFileOut, dy, "dy");   // 1
    WriteParamToFile<FPNumber>(paramFileOut, dz, "dz");   // 2
    WriteParamToFile<FPNumber>(paramFileOut, unit_length_si, "unit_length_si");   // 3
    WriteParamToFile<FPNumber>(paramFileOut, pitch_to_unitlength, "pitch_to_unitlength");   // 4
    WriteParamToFile<FPNumber>(paramFileOut, fwhm_to_pitch, "fwhm_to_pitch");   // 5
    WriteParamToFile<FPNumber>(paramFileOut, wp_2p_thz, "wp_2p_thz");   // 6
    WriteParamToFile<FPNumber>(paramFileOut, gamma_thz, "gamma_thz");   // 7
    WriteParamToFile<FPNumber>(paramFileOut, wp_switch_dt, "wp_switch_dt");   // 8
    WriteParamToFile<std::size_t>(paramFileOut, EorWp_entire_sample_rate, "EorWp_entire_sample_rate");   // 9
    WriteParamToFile<FPNumber>(paramFileOut, er_z_record, "er_z_record");   // 10
    WriteParamToFile<FPNumber>(paramFileOut, et_z_record, "et_z_record");   // 11
    WriteParamToFile<FPNumber>(paramFileOut, theta_deg, "theta_deg");   // 12
    WriteParamToFile<FPNumber>(paramFileOut, celldisplacement_to_pitch, "celldisplacement_to_pitch");   // 13
    paramFileOut.close();
};

#endif // TEST_FDTD_2D_LARGEPERIODICLATTICE_GAUSSIANPLASMA_TIMESWITCH_JSON

