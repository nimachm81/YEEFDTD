

#ifndef TEST_FDTD_2D_SWITCHEDPLASMA_CURVEDLINES_GRIDCOLLECTION_JSON
#define TEST_FDTD_2D_SWITCHEDPLASMA_CURVEDLINES_GRIDCOLLECTION_JSON

#include<typeinfo>

#include "boost/lexical_cast.hpp"
#include "NumberTypes.h"
#include "ParameterExtractor.h"
#include "ParamFileTranslator.h"
#include "UtilityFunctions.hpp"


void test_run_fdtd_curved_lines_gaussian_plasma_time_switch_2d_from_json(
                        FPNumber lineDistance_in_unitlength = 1.2e-3/124.0e-6,  // line distance in natural units
                        FPNumber radius_right_in_unitlength = 128.0e-3/124.0e-6,    // radius of curvature of the right line
                        FPNumber radius_left_in_unitlength = 10000.0,      // radius of curvature of the left line
                        FPNumber fwhm_left_in_unitlength = 42.0/124.0,    // fwhm in natural units
                        FPNumber fwhm_right_in_unitlength = 64.0/124.0,    // fwhm in natural units
                        FPNumber wp_2p_thz = 12.75,   // wp/2pi
                        FPNumber gamma_thz = 1.0,   // scattering rate
                        FPNumber wp_switch_dt = 0.1          // switch time in natural units
                        ) {
    FPNumber unit_length_si = 124.0e-6;     // unit length in SI
    FPNumber lineDistance_si = lineDistance_in_unitlength*unit_length_si;
    FPNumber fwhm_left_si = fwhm_left_in_unitlength*unit_length_si;
    FPNumber fwhm_right_si = fwhm_right_in_unitlength*unit_length_si;
    FPNumber radius_right_si = radius_right_in_unitlength*unit_length_si;
    FPNumber radius_left_si = radius_left_in_unitlength*unit_length_si;

    std::cout << "lineDistance_si : " << lineDistance_si << std::endl;
    std::cout << "fwhm_left_si : " << fwhm_left_si << std::endl;
    std::cout << "fwhm_right_si : " << fwhm_right_si << std::endl;
    std::cout << "radius_right_in_unitlength : " << radius_right_si << std::endl;
    std::cout << "radius_left_in_unitlength : " << radius_left_si << std::endl;

    FPNumber eps_r = 11.7;      // silicon

    std::size_t numOfSamplesPerUnitLength = 50;

    FPNumber y0 = -40.0;
    FPNumber y1 = 40.0;
    FPNumber z0 = -7.0;
    FPNumber z1 = 7.0;
    FPNumber pml_r_z0 = z1;
    FPNumber pml_r_z1 = z1 + 2.0;
    FPNumber pml_l_z0 = z0 - 2.0;
    FPNumber pml_l_z1 = z0;

    std::size_t ny = static_cast<std::size_t>(std::real(y1 - y0) * numOfSamplesPerUnitLength);
    std::size_t nz = static_cast<std::size_t>(std::real(z1 - z0) * numOfSamplesPerUnitLength);
    std::size_t pml_r_nz = static_cast<std::size_t>(std::real(pml_r_z1 - pml_r_z0) * numOfSamplesPerUnitLength);
    std::size_t pml_l_nz = static_cast<std::size_t>(std::real(pml_l_z1 - pml_l_z0) * numOfSamplesPerUnitLength);

    FPNumber dy = (y1 - y0)/(FPNumber)(ny);
    FPNumber dz = (z1 - z0)/(FPNumber)(nz);
    FPNumber pml_r_dz = (pml_r_z1 - pml_r_z0)/(FPNumber)(pml_r_nz);
    FPNumber pml_l_dz = (pml_l_z1 - pml_l_z0)/(FPNumber)(pml_l_nz);

    FPNumber pml_r_cube_dz = 1.0;
    FPNumber pml_r_cube_sige_z0 = z0;
    FPNumber pml_r_cube_sige_z1 = z1 + pml_r_cube_dz/2.0 + 0.1;
    FPNumber pml_r_cube_sigh_z0 = pml_r_cube_sige_z0;
    FPNumber pml_r_cube_sigh_z1 = pml_r_cube_sige_z1 + pml_r_dz/2.0;

    FPNumber pml_l_cube_dz = 1.0;
    FPNumber pml_l_cube_sige_z0 = z0 - pml_l_cube_dz/2.0 - 0.1;
    FPNumber pml_l_cube_sige_z1 = z1;
    FPNumber pml_l_cube_sigh_z0 = pml_l_cube_sige_z0 - pml_l_dz/2.0;
    FPNumber pml_l_cube_sigh_z1 = pml_l_cube_sige_z1;

    FPNumber pml_r_sig_E = 2.0;
    FPNumber pml_r_sig_H = 2.0;
    FPNumber pml_l_sig_E = 2.0;
    FPNumber pml_l_sig_H = 2.0;

    FPNumber stabilityFactor = 0.95;
    FPNumber dt = (FPNumber)1.0/std::sqrt((FPNumber)1.0/(dy*dy) + (FPNumber)1.0/(dz*dz))*stabilityFactor;
    FPNumber z_j = -3.0;
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

    FPNumber wp2_decayrate_r_left = UtilityFunctions::FWHMtoDecayRate(fwhm_left_in_unitlength);
    FPNumber wp2_decayrate_r_right = UtilityFunctions::FWHMtoDecayRate(fwhm_right_in_unitlength);
    FPNumber wp2_center_z_right = -radius_right_in_unitlength + lineDistance_in_unitlength/2.0;
    FPNumber wp2_center_z_left = -radius_left_in_unitlength - lineDistance_in_unitlength/2.0;

    FPNumber wp2_mask_z0 = z0 + 0.5;
    FPNumber wp2_mask_z1 = z1 - 0.5;
    FPNumber wp2_mask_t0 = -z_j*std::sqrt(eps_r) + j_center_t;
    FPNumber wp2_mask_t1 = 1000.0;
    FPNumber wp2_mask_dt = wp_switch_dt;

    std::size_t numOfTimeSamplesBeforeSwitch = static_cast<std::size_t>(std::real((wp2_mask_t0 - wp2_mask_dt)/dt));
    std::size_t numOfTimeSamplesDuringSwitch = static_cast<std::size_t>(std::real((FPNumber)2.0*wp2_mask_dt/dt));
    std::size_t numOfTimeSamplesAfterSwitch = 15000;

    std::size_t nt_0 = numOfTimeSamplesBeforeSwitch;
    std::size_t nt_1 = nt_0 + numOfTimeSamplesDuringSwitch;
    std::size_t nt_2 = nt_1 + numOfTimeSamplesAfterSwitch;

    std::cout << "Nt : " << nt_2 << std::endl;

    FPNumber ei_z_record = z_j + 0.2;
    FPNumber er_z_record = z0 + 0.2;
    FPNumber et_z_record = z1 - 0.2;
    FPNumber ec_z_record = 0.0;
    FPNumber wp2_z_record = lineDistance_in_unitlength/2.0;
    std::size_t ei_indz_record = std::round(std::real((ei_z_record - z0)/dz));
    std::size_t er_indz_record = std::round(std::real((er_z_record - z0)/dz));
    std::size_t et_indz_record = std::round(std::real((et_z_record - z0)/dz));
    std::size_t ec_indz_record = std::round(std::real((ec_z_record - z0)/dz));
    std::size_t wp2_indz_record = std::round(std::real((wp2_z_record - z0)/dz));

    FPNumber field_entire_z_m_corner_record = z0 + 0.1;
    FPNumber field_entire_z_p_corner_record = z1 - 0.1;
    FPNumber field_entire_y_m_corner_record = -6.0;
    FPNumber field_entire_y_p_corner_record = +6.0;
    std::size_t field_entire_ind_z_m_record = std::round(std::real((field_entire_z_m_corner_record - z0)/dz));
    std::size_t field_entire_ind_z_p_record = std::round(std::real((field_entire_z_p_corner_record - z0)/dz));
    std::size_t field_entire_ind_y_m_record = std::round(std::real((field_entire_y_m_corner_record - y0)/dy));
    std::size_t field_entire_ind_y_p_record = std::round(std::real((field_entire_y_p_corner_record - y0)/dy));


    FPNumber wp_switch_dt_si = wp_switch_dt/((FPNumber)3.0e8/unit_length_si);

    std::string outputFolder = "GaussianPlasmaCurvedLines-TimeSwitched/";

    std::string file_suffix =
                 "-fp=" + boost::lexical_cast<std::string>(std::real(wp_2p_thz)) +
                 "-gamma=" + boost::lexical_cast<std::string>(std::real(gamma_thz)) +
                 "-lineDistance=" + boost::lexical_cast<std::string>(std::real(lineDistance_si*(FPNumber)1.0e6)) +
                 "-fwhmLeft=" + boost::lexical_cast<std::string>(std::real(fwhm_left_si*(FPNumber)1.0e6)) +
                 "-fwhmRight=" + boost::lexical_cast<std::string>(std::real(fwhm_right_si*(FPNumber)1.0e6)) +
                 "-swithTime=" + boost::lexical_cast<std::string>(std::real(wp_switch_dt)) +
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
    std::string e_slice_c_name = std::string("\"") + outputFolder + "Ec-x-slice" +
                                 file_suffix +
                                 std::string("\"");
    std::string h_slice_c_name = std::string("\"") + outputFolder + "Hc-y-slice" +
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
    std::string pml_r_e_name = std::string("\"") + outputFolder + "pml-r-E-x" +
                                 file_suffix +
                                 std::string("\"");
    std::string pml_l_e_name = std::string("\"") + outputFolder + "pml-l-E-x" +
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
            {"\"_m_dt_eps_\"", boost::lexical_cast<std::string>(std::real(-dt/eps_r))},
            {"\"_1_dy_\"", boost::lexical_cast<std::string>(std::real(1.0/dy))},
            {"\"_m_1_dy_\"", boost::lexical_cast<std::string>(std::real(-1.0/dy))},
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
            {"\"_wp_sq_left_center_z_\"", boost::lexical_cast<std::string>(std::real(wp2_center_z_left))},
            {"\"_wp_sq_right_center_z_\"", boost::lexical_cast<std::string>(std::real(wp2_center_z_right))},
            {"\"_wp_sq_left_radius_\"", boost::lexical_cast<std::string>(std::real(radius_left_in_unitlength))},
            {"\"_wp_sq_right_radius_\"", boost::lexical_cast<std::string>(std::real(radius_right_in_unitlength))},
            {"\"_wp_sq_left_decayrate_r_\"", boost::lexical_cast<std::string>(std::real(wp2_decayrate_r_left))},
            {"\"_wp_sq_right_decayrate_r_\"", boost::lexical_cast<std::string>(std::real(wp2_decayrate_r_right))},
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
            {"\"_e_ind_z_record_c_\"", boost::lexical_cast<std::string>(ec_indz_record)},
            {"\"_e_ind_z_record_c_p1_\"", boost::lexical_cast<std::string>(ec_indz_record + 1)},
            {"\"_wp2_ind_z_record_\"", boost::lexical_cast<std::string>(wp2_indz_record)},
            {"\"_wp2_ind_z_record_p1_\"", boost::lexical_cast<std::string>(wp2_indz_record + 1)},
            {"\"_nt_0_\"", boost::lexical_cast<std::string>(nt_0)},
            {"\"_nt_1_\"", boost::lexical_cast<std::string>(nt_1)},
            {"\"_nt_2_\"", boost::lexical_cast<std::string>(nt_2)},
            {"\"_Eslice_name_i_\"", e_slice_i_name},
            {"\"_Eslice_name_r_\"", e_slice_r_name},
            {"\"_Eslice_name_t_\"", e_slice_t_name},
            {"\"_Eslice_name_c_\"", e_slice_c_name},
            {"\"_Hslice_name_c_\"", h_slice_c_name},
            {"\"_E_name_\"", e_name},
            {"\"_Wp2Slice_name_\"", wp2_slice_name},
            {"\"_Wp2_name_\"", wp2_name},
            {"\"_sample_rate_\"", boost::lexical_cast<std::string>(EorWp_entire_sample_rate)},
            {"\"_e_entire_ind_ym_record_\"", boost::lexical_cast<std::string>(field_entire_ind_y_m_record)},
            {"\"_e_entire_ind_yp_record_\"", boost::lexical_cast<std::string>(field_entire_ind_y_p_record)},
            {"\"_e_entire_ind_zm_record_\"", boost::lexical_cast<std::string>(field_entire_ind_z_m_record)},
            {"\"_e_entire_ind_zp_record_\"", boost::lexical_cast<std::string>(field_entire_ind_z_p_record)},
            {"\"_gr_z0_\"", boost::lexical_cast<std::string>(std::real(pml_r_z0))},
            {"\"_gr_z1_\"", boost::lexical_cast<std::string>(std::real(pml_r_z1))},
            {"\"_gr_nz_\"", boost::lexical_cast<std::string>(pml_r_nz)},
            {"\"_gr_nz_p1_\"", boost::lexical_cast<std::string>(pml_r_nz + 1)},
            {"\"_gr_cube_sige_z0_\"", boost::lexical_cast<std::string>(std::real(pml_r_cube_sige_z0))},
            {"\"_gr_cube_sige_z1_\"", boost::lexical_cast<std::string>(std::real(pml_r_cube_sige_z1))},
            {"\"_gr_cube_sigh_z0_\"", boost::lexical_cast<std::string>(std::real(pml_r_cube_sigh_z0))},
            {"\"_gr_cube_sigh_z1_\"", boost::lexical_cast<std::string>(std::real(pml_r_cube_sigh_z1))},
            {"\"_gr_cube_dz_\"", boost::lexical_cast<std::string>(std::real(pml_r_cube_dz))},
            {"\"_gr_dt_dz_\"", boost::lexical_cast<std::string>(std::real(dt/pml_r_dz))},
            {"\"_gr_m_dt_dz_\"", boost::lexical_cast<std::string>(std::real(-dt/pml_r_dz))},
            {"\"_gr_dt_dz_eps_\"", boost::lexical_cast<std::string>(std::real(dt/pml_r_dz/eps_r))},
            {"\"_gr_m_dt_dz_eps_\"", boost::lexical_cast<std::string>(std::real(-dt/pml_r_dz/eps_r))},
            {"\"_gr_sig_e_\"", boost::lexical_cast<std::string>(std::real(pml_r_sig_E))},
            {"\"_gr_sig_h_\"", boost::lexical_cast<std::string>(std::real(pml_r_sig_H))},
            {"\"_gr_E_name_\"", pml_r_e_name},
            {"\"_gl_z0_\"", boost::lexical_cast<std::string>(std::real(pml_l_z0))},
            {"\"_gl_z1_\"", boost::lexical_cast<std::string>(std::real(pml_l_z1))},
            {"\"_gl_nz_\"", boost::lexical_cast<std::string>(pml_l_nz)},
            {"\"_gl_nz_p1_\"", boost::lexical_cast<std::string>(pml_l_nz + 1)},
            {"\"_gl_nz_m1_\"", boost::lexical_cast<std::string>(pml_l_nz - 1)},
            {"\"_gl_cube_sige_z0_\"", boost::lexical_cast<std::string>(std::real(pml_l_cube_sige_z0))},
            {"\"_gl_cube_sige_z1_\"", boost::lexical_cast<std::string>(std::real(pml_l_cube_sige_z1))},
            {"\"_gl_cube_sigh_z0_\"", boost::lexical_cast<std::string>(std::real(pml_l_cube_sigh_z0))},
            {"\"_gl_cube_sigh_z1_\"", boost::lexical_cast<std::string>(std::real(pml_l_cube_sigh_z1))},
            {"\"_gl_cube_dz_\"", boost::lexical_cast<std::string>(std::real(pml_l_cube_dz))},
            {"\"_gl_dt_dz_\"", boost::lexical_cast<std::string>(std::real(dt/pml_l_dz))},
            {"\"_gl_m_dt_dz_\"", boost::lexical_cast<std::string>(std::real(-dt/pml_l_dz))},
            {"\"_gl_dt_dz_eps_\"", boost::lexical_cast<std::string>(std::real(dt/pml_l_dz/eps_r))},
            {"\"_gl_m_dt_dz_eps_\"", boost::lexical_cast<std::string>(std::real(-dt/pml_l_dz/eps_r))},
            {"\"_gl_sig_e_\"", boost::lexical_cast<std::string>(std::real(pml_l_sig_E))},
            {"\"_gl_sig_h_\"", boost::lexical_cast<std::string>(std::real(pml_l_sig_H))},
            {"\"_gl_E_name_\"", pml_l_e_name}
            };

    // for parallel processing this file's name should be unique to the processor. Here the suffix theta is supposes to
    // make it unique (each process is assumed treating a different theta)
    std::string processedFileName =
            std::string("instructions/processed/MaxwellYee2D_CurvedLinesGaussianPlasma_TimeSwitch_GridCollection_processed_") +
            file_suffix + ".json";

    ParameterExtractor::ReplaceStringsInFile("instructions/MaxwellYee2D_CurvedLinesGaussianPlasma_TimeSwitch_GridCollection.json",
                processedFileName, str_replacewith);

    ParamFileTranslator fileTranslator(processedFileName);
    fileTranslator.Translate();

    std::string parametersFileName =
            std::string("data/") + outputFolder + "params" + file_suffix + ".param";

    std::ofstream paramFileOut(parametersFileName.c_str(), std::ios::out | std::ios::binary);
    UtilityFunctions::WriteParamToFile<FPNumber>(paramFileOut, dt, "dt");   // 0
    UtilityFunctions::WriteParamToFile<FPNumber>(paramFileOut, dy, "dy");   // 1
    UtilityFunctions::WriteParamToFile<FPNumber>(paramFileOut, dz, "dz");   // 2
    UtilityFunctions::WriteParamToFile<FPNumber>(paramFileOut, unit_length_si, "unit_length_si");   // 3
    UtilityFunctions::WriteParamToFile<FPNumber>(paramFileOut, lineDistance_in_unitlength, "lineDistance_in_unitlength");   // 4
    UtilityFunctions::WriteParamToFile<FPNumber>(paramFileOut, radius_left_in_unitlength, "radius_left_in_unitlength");   // 4
    UtilityFunctions::WriteParamToFile<FPNumber>(paramFileOut, radius_right_in_unitlength, "radius_right_in_unitlength");   // 4
    UtilityFunctions::WriteParamToFile<FPNumber>(paramFileOut, fwhm_left_in_unitlength, "fwhm_left_in_unitlength");   // 5
    UtilityFunctions::WriteParamToFile<FPNumber>(paramFileOut, fwhm_right_in_unitlength, "fwhm_right_in_unitlength");   // 5
    UtilityFunctions::WriteParamToFile<FPNumber>(paramFileOut, wp_2p_thz, "wp_2p_thz");   // 6
    UtilityFunctions::WriteParamToFile<FPNumber>(paramFileOut, gamma_thz, "gamma_thz");   // 7
    UtilityFunctions::WriteParamToFile<FPNumber>(paramFileOut, wp_switch_dt, "wp_switch_dt");   // 8
    UtilityFunctions::WriteParamToFile<std::size_t>(paramFileOut, EorWp_entire_sample_rate, "EorWp_entire_sample_rate");   // 9
    UtilityFunctions::WriteParamToFile<FPNumber>(paramFileOut, er_z_record, "er_z_record");   // 10
    UtilityFunctions::WriteParamToFile<FPNumber>(paramFileOut, et_z_record, "et_z_record");   // 11
    paramFileOut.close();
};

#endif // TEST_FDTD_2D_SWITCHEDPLASMA_CURVEDLINES_GRIDCOLLECTION_JSON

