


#include "NumberTypes.h"
#include "PhysicalUnits.hpp"
#include "UtilityFunctions.hpp"

void test_run_fdtd_2d_metal_wedge_electron_emitter_pureScatteredFDTD_gridCollection_from_json() {
    FPNumber fdtd_unit_length = 10.0e-6;
    PhysicalUnits units(fdtd_unit_length);

    FPNumber y0 = -0.5;
    FPNumber y1 = 0.5;
    FPNumber z0 = -0.5;
    FPNumber z1 = 0.5;
    FPNumber pml_thickness = 0.3;
    FPNumber pml_r_z0 = z1;
    FPNumber pml_r_z1 = z1 + pml_thickness;
    FPNumber pml_l_z0 = z0 - pml_thickness;
    FPNumber pml_l_z1 = z0;
    FPNumber pml_t_y0 = y1;
    FPNumber pml_t_y1 = y1 + pml_thickness;

    std::size_t numOfSamplesPerUnitLength = 300;

    std::size_t ny = static_cast<std::size_t>(std::real(y1 - y0) * numOfSamplesPerUnitLength);
    std::size_t nz = static_cast<std::size_t>(std::real(z1 - z0) * numOfSamplesPerUnitLength);
    std::size_t pml_r_nz = static_cast<std::size_t>(std::real(pml_r_z1 - pml_r_z0) * numOfSamplesPerUnitLength);
    std::size_t pml_l_nz = static_cast<std::size_t>(std::real(pml_l_z1 - pml_l_z0) * numOfSamplesPerUnitLength);
    std::size_t pml_t_ny = static_cast<std::size_t>(std::real(pml_t_y1 - pml_t_y0) * numOfSamplesPerUnitLength);

    FPNumber dy = (y1 - y0)/(FPNumber)(ny);
    FPNumber dz = (z1 - z0)/(FPNumber)(nz);
    FPNumber pml_r_dz = (pml_r_z1 - pml_r_z0)/(FPNumber)(pml_r_nz);
    FPNumber pml_l_dz = (pml_l_z1 - pml_l_z0)/(FPNumber)(pml_l_nz);
    FPNumber pml_t_dy = (pml_t_y1 - pml_t_y0)/(FPNumber)(pml_t_ny);

    FPNumber pml_r_cube_dz = (pml_r_z1 - pml_r_z0) / 2.0;
    FPNumber pml_r_cube_sige_z0 = z0;
    FPNumber pml_r_cube_sige_z1 = z1 + pml_r_cube_dz/2.0 + 2.0*pml_r_dz;
    FPNumber pml_r_cube_sigh_z0 = pml_r_cube_sige_z0;
    FPNumber pml_r_cube_sigh_z1 = pml_r_cube_sige_z1 + pml_r_dz/2.0;

    FPNumber pml_l_cube_dz = (pml_l_z1 - pml_l_z0) / 2.0;
    FPNumber pml_l_cube_sige_z0 = z0 - pml_l_cube_dz/2.0 - 2.0*pml_l_dz;
    FPNumber pml_l_cube_sige_z1 = z1;
    FPNumber pml_l_cube_sigh_z0 = pml_l_cube_sige_z0 + pml_r_dz/2.0;
    FPNumber pml_l_cube_sigh_z1 = pml_l_cube_sige_z1;

    FPNumber pml_t_cube_dy = (pml_t_y1 - pml_t_y0) / 2.0;
    FPNumber pml_t_cube_sige_y0 = y0;
    FPNumber pml_t_cube_sige_y1 = y1 + pml_t_cube_dy/2.0 + 2.0*pml_t_dy;
    FPNumber pml_t_cube_sigh_y0 = pml_t_cube_sige_y0;
    FPNumber pml_t_cube_sigh_y1 = pml_t_cube_sige_y1 + pml_t_dy/2.0;

    FPNumber sig_eh = 15.0;
    FPNumber pml_r_sig_E = sig_eh;
    FPNumber pml_r_sig_H = sig_eh;
    FPNumber pml_l_sig_E = sig_eh;
    FPNumber pml_l_sig_H = sig_eh;
    FPNumber pml_t_sig_E = sig_eh;
    FPNumber pml_t_sig_H = sig_eh;


    FPNumber stabilityFactor = 0.95;
    FPNumber dt = (FPNumber)1.0/std::sqrt((FPNumber)1.0/(dy*dy) + (FPNumber)1.0/(dz*dz))*stabilityFactor;

    const FPNumber eps_r = 1.0;     // only to set jm_amplitude... for eps_r != 1 json file should be updated

    FPNumber eFieldMax_SI = 2.0e8;     // V/m
    FPNumber eFieldMax_FD = units.ConvertSIElectricFieldToFDUnits(eFieldMax_SI);
    std::cout << "eFieldMax_SI: " << eFieldMax_SI << " ,eFieldMax_FD: " << eFieldMax_FD << std::endl;

    FPNumber eField_FD_convertto_SI = units.ConvertFDElectricFieldToSIUnits(1.0);

    //plane wave (pw) parameters
    std::array<FPNumber, 3> pw_direction{0.0, 0.0, 1.0};
    FPNumber pw_velocity = 1.0;
    FPNumber pw_mod_freq = units.ConvertSIFrequencyToFDUnits(3.0e12);
    FPNumber pw_mod_phase = M_PI/2.0;
    FPNumber pw_decay_rate_t = pw_mod_freq / 0.5;
    FPNumber pw_center_t = 1.0 / pw_mod_freq;
    FPNumber pw_amplitude = -eFieldMax_FD;
    FPNumber pw_ey_amplitude = pw_amplitude;
    FPNumber pw_hx_amplitude = -pw_amplitude;

    std::cout << "pw_mod_freq : " << pw_mod_freq << std::endl;
    std::cout << "pw_center_t : " << pw_center_t << " , nt_center : " << (pw_center_t/dt) << std::endl;


    // wedge parameters
    FPNumber wedgeAngle = 4.0/180.0*M_PI;
    FPNumber wedgeTipRadius = units.ConvertSILengthToFDUnits(100.0e-9);
    FPNumber wedgeHeight = 1.0;
    std::array<FPNumber, 3> wedgeTipPosition{0.0, 0.0, 0.0};
    FPNumber wedgeTopHeight = 2.0*wedgeTipRadius;  // only this part of the wedge is discretized for the emitter

    FPNumber maxSurfElemSize = 0.01*wedgeTipRadius;
    std::cout << "wedgeTipRadius: " << wedgeTipRadius << " , maxSurfElemSize: " << maxSurfElemSize << std::endl;

    // particles parameters
    FPNumber electronCharge = -units.GetElectronChargeInFDUnits();
    FPNumber electronMass = units.GetElectronMassInFDUnits();
    FPNumber holeCharge = -electronCharge;

    std::cout << "electronCharge: " << electronCharge << " , electronMass: " << electronMass << std::endl;

    FPNumber plasmaFrequency = units.ConvertSIFrequencyToFDUnits(300.0e12);
    FPNumber gamma = units.ConvertSIFrequencyToFDUnits(5.0e12);;   // scattering rate
    std::cout << "plasmaFrequency: " << plasmaFrequency << " , gamma: " << gamma << std::endl;


    FPNumber t_max = units.ConvertSITimeToFDUnits(1.0e-12);
    FPNumber dt_data_save = units.ConvertSITimeToFDUnits(0.01e-12);

    std::size_t data_save_rate = (std::size_t)(dt_data_save/dt);
    std::size_t numOfTimeSamples = (std::size_t)(t_max/dt);

    std::cout << "data_save_rate: " << data_save_rate << std::endl;
    std::cout << "numOfTimeSamples: " << numOfTimeSamples << std::endl;

    std::unordered_map<std::string, std::string> str_replacewith{
            {"\"_y0_\"", boost::lexical_cast<std::string>(std::real(y0))},
            {"\"_y1_\"", boost::lexical_cast<std::string>(std::real(y1))},
            {"\"_z0_\"", boost::lexical_cast<std::string>(std::real(z0))},
            {"\"_z1_\"", boost::lexical_cast<std::string>(std::real(z1))},
            {"\"_ny_\"", boost::lexical_cast<std::string>(ny)},
            {"\"_nz_\"", boost::lexical_cast<std::string>(nz)},
            {"\"_ny_p1_\"", boost::lexical_cast<std::string>(ny + 1)},
            {"\"_ny_m1_\"", boost::lexical_cast<std::string>(ny - 1)},
            {"\"_nz_p1_\"", boost::lexical_cast<std::string>(nz + 1)},
            {"\"_nz_m1_\"", boost::lexical_cast<std::string>(nz - 1)},
            {"\"_dy_\"", boost::lexical_cast<std::string>(std::real(dy))},
            {"\"_dz_\"", boost::lexical_cast<std::string>(std::real(dz))},
            {"\"_dt_\"", boost::lexical_cast<std::string>(std::real(dt))},
            {"\"_m_dt_\"", boost::lexical_cast<std::string>(std::real(-dt))},
            {"\"_dt_2_\"", boost::lexical_cast<std::string>(std::real(dt/2.0))},
            {"\"_m_dt_2_\"", boost::lexical_cast<std::string>(std::real(-dt/2.0))},
            {"\"_dt_dy_\"", boost::lexical_cast<std::string>(std::real(dt/dy))},
            {"\"_m_dt_dy_\"", boost::lexical_cast<std::string>(std::real(-dt/dy))},
            {"\"_m_dtdz_dy_\"", boost::lexical_cast<std::string>(std::real(-dt*dz/dy))},
            {"\"_dt_dy_dz_\"", boost::lexical_cast<std::string>(std::real(dt/dy/(dz)))},
            {"\"_m_dt_dy_dz_\"", boost::lexical_cast<std::string>(std::real(-dt/dy/(dz)))},
            {"\"_1_dy_\"", boost::lexical_cast<std::string>(std::real(1.0/dy))},
            {"\"_m_1_dy_\"", boost::lexical_cast<std::string>(std::real(-1.0/dy))},
            {"\"_dt_dz_\"", boost::lexical_cast<std::string>(std::real(dt/dz))},
            {"\"_m_dt_dz_\"", boost::lexical_cast<std::string>(std::real(-dt/dz))},
            {"\"_m_dtdy_dz_\"", boost::lexical_cast<std::string>(std::real(-dt*dy/dz))},
            {"\"_dt_dz_dy_\"", boost::lexical_cast<std::string>(std::real(dt/dz/(dy)))},
            {"\"_m_dt_dz_dy_\"", boost::lexical_cast<std::string>(std::real(-dt/dz/(dy)))},
            {"\"_1_dz_\"", boost::lexical_cast<std::string>(std::real(1.0/dz))},
            {"\"_m_1_dz_\"", boost::lexical_cast<std::string>(std::real(-1.0/dz))},
            {"\"_gr_z0_\"", boost::lexical_cast<std::string>(std::real(pml_r_z0))},
            {"\"_gr_z1_\"", boost::lexical_cast<std::string>(std::real(pml_r_z1))},
            {"\"_gr_nz_\"", boost::lexical_cast<std::string>(pml_r_nz)},
            {"\"_gr_nz_p1_\"", boost::lexical_cast<std::string>(pml_r_nz + 1)},
            {"\"_gr_nz_m1_\"", boost::lexical_cast<std::string>(pml_r_nz - 1)},
            {"\"_gr_cube_sige_z0_\"", boost::lexical_cast<std::string>(std::real(pml_r_cube_sige_z0))},
            {"\"_gr_cube_sige_z1_\"", boost::lexical_cast<std::string>(std::real(pml_r_cube_sige_z1))},
            {"\"_gr_cube_sigh_z0_\"", boost::lexical_cast<std::string>(std::real(pml_r_cube_sigh_z0))},
            {"\"_gr_cube_sigh_z1_\"", boost::lexical_cast<std::string>(std::real(pml_r_cube_sigh_z1))},
            {"\"_gr_cube_dz_\"", boost::lexical_cast<std::string>(std::real(pml_r_cube_dz))},
            {"\"_gr_dt_dz_\"", boost::lexical_cast<std::string>(std::real(dt/pml_r_dz))},
            {"\"_gr_m_dt_dz_\"", boost::lexical_cast<std::string>(std::real(-dt/pml_r_dz))},
            {"\"_gr_1_dz_\"", boost::lexical_cast<std::string>(std::real(1.0/pml_r_dz))},
            {"\"_gr_m_1_dz_\"", boost::lexical_cast<std::string>(std::real(-1.0/pml_r_dz))},
            {"\"_gr_sig_e_\"", boost::lexical_cast<std::string>(std::real(pml_r_sig_E))},
            {"\"_gr_sig_h_\"", boost::lexical_cast<std::string>(std::real(pml_r_sig_H))},
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
            {"\"_gl_1_dz_\"", boost::lexical_cast<std::string>(std::real(1.0/pml_l_dz))},
            {"\"_gl_m_1_dz_\"", boost::lexical_cast<std::string>(std::real(-1.0/pml_l_dz))},
            {"\"_gl_sig_e_\"", boost::lexical_cast<std::string>(std::real(pml_l_sig_E))},
            {"\"_gl_sig_h_\"", boost::lexical_cast<std::string>(std::real(pml_l_sig_H))},
            {"\"_gt_y0_\"", boost::lexical_cast<std::string>(std::real(pml_t_y0))},
            {"\"_gt_y1_\"", boost::lexical_cast<std::string>(std::real(pml_t_y1))},
            {"\"_gt_ny_\"", boost::lexical_cast<std::string>(pml_t_ny)},
            {"\"_gt_ny_p1_\"", boost::lexical_cast<std::string>(pml_t_ny + 1)},
            {"\"_gt_ny_m1_\"", boost::lexical_cast<std::string>(pml_t_ny - 1)},
            {"\"_gt_cube_sige_y0_\"", boost::lexical_cast<std::string>(std::real(pml_t_cube_sige_y0))},
            {"\"_gt_cube_sige_y1_\"", boost::lexical_cast<std::string>(std::real(pml_t_cube_sige_y1))},
            {"\"_gt_cube_sigh_y0_\"", boost::lexical_cast<std::string>(std::real(pml_t_cube_sigh_y0))},
            {"\"_gt_cube_sigh_y1_\"", boost::lexical_cast<std::string>(std::real(pml_t_cube_sigh_y1))},
            {"\"_gt_cube_dy_\"", boost::lexical_cast<std::string>(std::real(pml_t_cube_dy))},
            {"\"_gt_dt_dy_\"", boost::lexical_cast<std::string>(std::real(dt/pml_t_dy))},
            {"\"_gt_m_dt_dy_\"", boost::lexical_cast<std::string>(std::real(-dt/pml_t_dy))},
            {"\"_gt_1_dy_\"", boost::lexical_cast<std::string>(std::real(1.0/pml_t_dy))},
            {"\"_gt_m_1_dy_\"", boost::lexical_cast<std::string>(std::real(-1.0/pml_t_dy))},
            {"\"_gt_sig_e_\"", boost::lexical_cast<std::string>(std::real(pml_t_sig_E))},
            {"\"_gt_sig_h_\"", boost::lexical_cast<std::string>(std::real(pml_t_sig_H))},
            {"\"_echarge_\"", boost::lexical_cast<std::string>(std::real(electronCharge))},
            {"\"_hcharge_\"", boost::lexical_cast<std::string>(std::real(holeCharge))},
            {"\"_mass_\"", boost::lexical_cast<std::string>(std::real(electronMass))},
            {"\"_echarge_mass_\"", boost::lexical_cast<std::string>(std::real(electronCharge/electronMass))},
            {"\"_surf_dl_\"", boost::lexical_cast<std::string>(std::real(maxSurfElemSize))},
            {"\"_Eyinc_amp_\"", boost::lexical_cast<std::string>(std::real(pw_ey_amplitude))},
            {"\"_Hxinc_amp_\"", boost::lexical_cast<std::string>(std::real(pw_hx_amplitude))},
            {"\"_planewave_direction_x_\"", boost::lexical_cast<std::string>(std::real(pw_direction[0]))},
            {"\"_planewave_direction_y_\"", boost::lexical_cast<std::string>(std::real(pw_direction[1]))},
            {"\"_planewave_direction_z_\"", boost::lexical_cast<std::string>(std::real(pw_direction[2]))},
            {"\"_planewave_velocity_\"", boost::lexical_cast<std::string>(std::real(pw_velocity))},
            {"\"_gaussianplanewave_decay_rate_\"", boost::lexical_cast<std::string>(std::real(pw_decay_rate_t))},
            {"\"_planewave_t_center_\"", boost::lexical_cast<std::string>(std::real(pw_center_t))},
            {"\"_planewave_freq_\"", boost::lexical_cast<std::string>(std::real(pw_mod_freq))},
            {"\"_planewave_phase_\"", boost::lexical_cast<std::string>(std::real(pw_mod_phase))},
            {"\"_wedge_angle_\"", boost::lexical_cast<std::string>(std::real(wedgeAngle))},
            {"\"_wedge_tip_radius_\"", boost::lexical_cast<std::string>(std::real(wedgeTipRadius))},
            {"\"_wedge_height_\"", boost::lexical_cast<std::string>(std::real(wedgeHeight))},
            {"\"_wedgetop_height_\"", boost::lexical_cast<std::string>(std::real(wedgeTopHeight))},
            {"\"_wedge_tip_y_\"", boost::lexical_cast<std::string>(std::real(wedgeTipPosition[1]))},
            {"\"_wedge_tip_z\"", boost::lexical_cast<std::string>(std::real(wedgeTipPosition[2]))},
            {"\"_fdtd_unit_length_\"", boost::lexical_cast<std::string>(std::real(fdtd_unit_length))},
            {"\"_gamma_\"", boost::lexical_cast<std::string>(std::real(gamma))},
            {"\"_wp_sq_\"", boost::lexical_cast<std::string>(std::real(plasmaFrequency*plasmaFrequency))},
            {"\"_save_rate_\"", boost::lexical_cast<std::string>(data_save_rate)},
            {"\"_nt_\"", boost::lexical_cast<std::string>(numOfTimeSamples)}
            };
    std::string fileName = "instructions/MaxwellYee2D_Metal_Wedge_electron_emitter_pureScatteredFDTD_GridCollection.json";

    ParameterExtractor::ReplaceStringsInFile(fileName,
                "instructions/processed/MaxwellYee2D_Metal_Wedge_electron_emitter_GridCollection_processed.json", str_replacewith);

    ParamFileTranslator fileTranslator("instructions/processed/MaxwellYee2D_Metal_Wedge_electron_emitter_GridCollection_processed.json");
    fileTranslator.Translate();


    std::string parametersFileName =
            std::string("data/2D/") + "params" + ".param";

    std::ofstream paramFileOut(parametersFileName.c_str(), std::ios::out | std::ios::binary);
    UtilityFunctions::WriteParamToFile<FPNumber>(paramFileOut, (FPNumber)units.ConvertFDTimeToSIUnits(dt), "dt");   // 0
    UtilityFunctions::WriteParamToFile<FPNumber>(paramFileOut, (FPNumber)units.ConvertFDLengthToSIUnits(dy), "dy");   // 1
    UtilityFunctions::WriteParamToFile<FPNumber>(paramFileOut, (FPNumber)units.ConvertFDLengthToSIUnits(dz), "dz");   // 2
    UtilityFunctions::WriteParamToFile<FPNumber>(paramFileOut, fdtd_unit_length, "fdtd_unit_length");
    UtilityFunctions::WriteParamToFile<FPNumber>(paramFileOut, eField_FD_convertto_SI, "eField_FD_convertto_SI");
    UtilityFunctions::WriteParamToFile<FPNumber>(paramFileOut, data_save_rate, "data_save_rate");
    UtilityFunctions::WriteParamToFile<FPNumber>(paramFileOut, (FPNumber)units.ConvertFDFrequencyToSIUnits(1.0), "frequency_conversion_factor");   // 0

}




