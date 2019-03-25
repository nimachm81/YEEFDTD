



#include "NumberTypes.h"
#include "PhysicalUnits.hpp"
#include "UtilityFunctions.hpp"
#include "WedgeGeometry.h"

void test_run_fdtd_2d_metal_wedge_electron_emitter_pureScatteredFDTD_partialGrid_gridCollection_from_json() {
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
    FPNumber pml_b_y0 = y0 - pml_thickness;
    FPNumber pml_b_y1 = y0;

    std::size_t numOfSamplesPerUnitLength = 300;

    std::size_t ny = static_cast<std::size_t>(std::real(y1 - y0) * numOfSamplesPerUnitLength);
    std::size_t nz = static_cast<std::size_t>(std::real(z1 - z0) * numOfSamplesPerUnitLength);
    std::size_t pml_r_nz = static_cast<std::size_t>(std::real(pml_r_z1 - pml_r_z0) * numOfSamplesPerUnitLength);
    std::size_t pml_l_nz = static_cast<std::size_t>(std::real(pml_l_z1 - pml_l_z0) * numOfSamplesPerUnitLength);
    std::size_t pml_t_ny = static_cast<std::size_t>(std::real(pml_t_y1 - pml_t_y0) * numOfSamplesPerUnitLength);
    std::size_t pml_b_ny = static_cast<std::size_t>(std::real(pml_b_y1 - pml_b_y0) * numOfSamplesPerUnitLength);

    FPNumber dy = (y1 - y0)/(FPNumber)(ny);
    FPNumber dz = (z1 - z0)/(FPNumber)(nz);
    FPNumber pml_r_dz = (pml_r_z1 - pml_r_z0)/(FPNumber)(pml_r_nz);
    FPNumber pml_l_dz = (pml_l_z1 - pml_l_z0)/(FPNumber)(pml_l_nz);
    FPNumber pml_t_dy = (pml_t_y1 - pml_t_y0)/(FPNumber)(pml_t_ny);
    FPNumber pml_b_dy = (pml_b_y1 - pml_b_y0)/(FPNumber)(pml_b_ny);

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

    FPNumber pml_b_cube_dy = (pml_b_y1 - pml_b_y0) / 2.0;
    FPNumber pml_b_cube_sige_y0 = y0 - pml_b_cube_dy/2.0 - 2.0*pml_b_dy;
    FPNumber pml_b_cube_sige_y1 = y1;
    FPNumber pml_b_cube_sigh_y0 = pml_b_cube_sige_y0 + pml_b_dy/2.0;
    FPNumber pml_b_cube_sigh_y1 = pml_b_cube_sige_y1;

    FPNumber sig_eh = 10.0;
    FPNumber pml_r_sig_E = sig_eh;
    FPNumber pml_r_sig_H = sig_eh;
    FPNumber pml_l_sig_E = sig_eh;
    FPNumber pml_l_sig_H = sig_eh;
    FPNumber pml_t_sig_E = sig_eh;
    FPNumber pml_t_sig_H = sig_eh;
    FPNumber pml_b_sig_E = sig_eh;
    FPNumber pml_b_sig_H = sig_eh;


    FPNumber stabilityFactor = 0.95;
    FPNumber dt = (FPNumber)1.0/std::sqrt((FPNumber)1.0/(dy*dy) + (FPNumber)1.0/(dz*dz))*stabilityFactor;

    const FPNumber eps_r = 1.0;     // only to set jm_amplitude... for eps_r != 1 json file should be updated

    FPNumber eFieldMax_SI = 5.0e7;     // V/m
    FPNumber eFieldMax_FD = units.ConvertSIElectricFieldToFDUnits(eFieldMax_SI);
    std::cout << "eFieldMax_SI: " << eFieldMax_SI << " ,eFieldMax_FD: " << eFieldMax_FD << std::endl;

    FPNumber eField_FD_convertto_SI = units.ConvertFDElectricFieldToSIUnits(1.0);

    //plane wave (pw) parameters
    std::array<FPNumber, 3> pw_direction{0.0, 0.0, 1.0};
    FPNumber pw_velocity = 1.0;
    FPNumber pw_mod_freq = units.ConvertSIFrequencyToFDUnits(1.0e12);
    FPNumber pw_mod_phase = M_PI/2.0;

    FPNumber pw_center_t = 0.5 / pw_mod_freq;

    FPNumber pw_decay_rate_t = 3.0 / pw_center_t;   // for Gaussian plane wave

    FPNumber pw_rect_width = 1.0 / pw_mod_freq;     // for rect plane wave
    FPNumber pw_rect_edge_width = 0.1 * pw_rect_width;  // for rect plane wave

    FPNumber pw_amplitude = eFieldMax_FD;
    FPNumber pw_ey_amplitude = +pw_amplitude;
    FPNumber pw_bx_amplitude = -pw_amplitude;

    std::cout << "pw_mod_freq : " << pw_mod_freq << std::endl;
    std::cout << "pw_center_t : " << pw_center_t << " , nt_center : " << (pw_center_t/dt) << std::endl;


    // wedge parameters
    FPNumber emitterWedgeAngle = 1.0/180.0*M_PI;
    FPNumber emitterWedgeTipRadius = units.ConvertSILengthToFDUnits(300.0e-9);
    std::array<FPNumber, 3> emitterWedgeTipPosition{0.0, -0.0, 0.0};
    FPNumber emitterWedgeHeight = emitterWedgeTipPosition[1] - y0;
    FPNumber emitterWedgeTopHeight = 2.0*emitterWedgeTipRadius;  // only this part of the wedge is discretized for the emitter

    FPNumber metalWedgeReductionSize = std::sqrt(dy*dy + dz*dz);    // about one cell smaller than the emitter wedge so that
                                                        // the emitter interacts with the normal electric field outside the metal
    FPNumber metalWedgeAngle = emitterWedgeAngle;
    FPNumber metalWedgeTipRadius = emitterWedgeTipRadius - metalWedgeReductionSize;
    std::array<FPNumber, 3> metalWedgeTipPosition = emitterWedgeTipPosition;
    metalWedgeTipPosition[1] -= metalWedgeReductionSize;
    FPNumber metalWedgeHeight = emitterWedgeHeight - metalWedgeReductionSize;


    FPNumber maxSurfElemSize = 0.01*emitterWedgeTipRadius;
    std::cout << "emitterWedgeTipRadius: " << emitterWedgeTipRadius << " , maxSurfElemSize: " << maxSurfElemSize << std::endl;
    std::cout << "emitterWedgeHeight: " << emitterWedgeHeight << std::endl;

    WedgeGeometry wedge;
    wedge.SetWedgeAngle(emitterWedgeAngle);
    wedge.SetTipRadius(emitterWedgeTipRadius);
    wedge.SetApexToBaseDistance(emitterWedgeHeight);
    wedge.SetApexPosition(emitterWedgeTipPosition);
    std::array<FPNumber, 3> wedge_r0;
    std::array<FPNumber, 3> wedge_r1;
    wedge.GetBoundingBox2D(0.0, wedge_r0, wedge_r1);
    std::intmax_t ind_wedge_y0 = (wedge_r0[1] - y0) / dy - 2;
    if(ind_wedge_y0 < 0) {
        ind_wedge_y0 = 0;
    } else if(ind_wedge_y0 > ny) {
        ind_wedge_y0 = ny;
    }
    std::intmax_t ind_wedge_y1 = (wedge_r1[1] - y0) / dy + 2;
    if(ind_wedge_y1 < 0) {
        ind_wedge_y1 = 1;
    } else if(ind_wedge_y1 > ny) {
        ind_wedge_y1 = ny;
    }

    std::intmax_t ind_wedge_z0 = (wedge_r0[2] - z0) / dz - 2;
    if(ind_wedge_z0 < 0) {
        ind_wedge_z0 = 0;
    } else if(ind_wedge_z0 > nz) {
        ind_wedge_z0 = nz;
    }
    std::intmax_t ind_wedge_z1 = (wedge_r1[2] - z0) / dz + 2;
    if(ind_wedge_z1 < 0) {
        ind_wedge_z1 = 1;
    } else if(ind_wedge_z1 > nz) {
        ind_wedge_z1 = nz;
    }

    std::size_t ny_wedge = ind_wedge_y1 - ind_wedge_y0;
    std::size_t nz_wedge = ind_wedge_z1 - ind_wedge_z0;

   std::intmax_t ind_wedge_base_z0 = std::ceil((wedge_r0[2] - z0) / dz) ;
    if(ind_wedge_base_z0 < 0) {
        ind_wedge_base_z0 = 0;
    } else if(ind_wedge_base_z0 > nz) {
        ind_wedge_base_z0 = nz;
    }
    std::intmax_t ind_wedge_base_z1 = std::floor((wedge_r1[2] - z0) / dz);
    if(ind_wedge_base_z1 < 0) {
        ind_wedge_base_z1 = 0;
    } else if(ind_wedge_base_z1 > nz) {
        ind_wedge_base_z1 = nz;
    }

    // particles parameters
    FPNumber electronCharge = -units.GetElectronChargeInFDUnits();
    FPNumber electronMass = units.GetElectronMassInFDUnits();
    FPNumber holeCharge = -electronCharge;

    std::cout << "electronCharge: " << electronCharge << " , electronMass: " << electronMass << std::endl;

    FPNumber plasmaFrequency = units.ConvertSIFrequencyToFDUnits(1000.0e12);
    FPNumber gamma = units.ConvertSIFrequencyToFDUnits(20.0e12);;   // scattering rate
    std::cout << "plasmaFrequency: " << plasmaFrequency << " , gamma: " << gamma << std::endl;


    FPNumber t_max = units.ConvertSITimeToFDUnits(1.5e-12);
    FPNumber dt_data_save = units.ConvertSITimeToFDUnits(0.01e-12);

    std::size_t particle_bunch_size = 1;

    std::size_t data_save_rate = (std::size_t)(dt_data_save/dt);
    std::size_t data_save_rate_2 = data_save_rate;
    std::size_t data_save_rate_pml = 10000; //data_save_rate;
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
            {"\"_gb_y0_\"", boost::lexical_cast<std::string>(std::real(pml_b_y0))},
            {"\"_gb_y1_\"", boost::lexical_cast<std::string>(std::real(pml_b_y1))},
            {"\"_gb_ny_\"", boost::lexical_cast<std::string>(pml_b_ny)},
            {"\"_gb_ny_p1_\"", boost::lexical_cast<std::string>(pml_b_ny + 1)},
            {"\"_gb_ny_m1_\"", boost::lexical_cast<std::string>(pml_b_ny - 1)},
            {"\"_gb_cube_sige_y0_\"", boost::lexical_cast<std::string>(std::real(pml_b_cube_sige_y0))},
            {"\"_gb_cube_sige_y1_\"", boost::lexical_cast<std::string>(std::real(pml_b_cube_sige_y1))},
            {"\"_gb_cube_sigh_y0_\"", boost::lexical_cast<std::string>(std::real(pml_b_cube_sigh_y0))},
            {"\"_gb_cube_sigh_y1_\"", boost::lexical_cast<std::string>(std::real(pml_b_cube_sigh_y1))},
            {"\"_gb_cube_dy_\"", boost::lexical_cast<std::string>(std::real(pml_b_cube_dy))},
            {"\"_gb_dt_dy_\"", boost::lexical_cast<std::string>(std::real(dt/pml_b_dy))},
            {"\"_gb_m_dt_dy_\"", boost::lexical_cast<std::string>(std::real(-dt/pml_b_dy))},
            {"\"_gb_1_dy_\"", boost::lexical_cast<std::string>(std::real(1.0/pml_b_dy))},
            {"\"_gb_m_1_dy_\"", boost::lexical_cast<std::string>(std::real(-1.0/pml_b_dy))},
            {"\"_gb_sig_e_\"", boost::lexical_cast<std::string>(std::real(pml_b_sig_E))},
            {"\"_gb_sig_h_\"", boost::lexical_cast<std::string>(std::real(pml_b_sig_H))},
            {"\"_echarge_\"", boost::lexical_cast<std::string>(std::real(electronCharge))},
            {"\"_hcharge_\"", boost::lexical_cast<std::string>(std::real(holeCharge))},
            {"\"_mass_\"", boost::lexical_cast<std::string>(std::real(electronMass))},
            {"\"_echarge_mass_\"", boost::lexical_cast<std::string>(std::real(electronCharge/electronMass))},
            {"\"_surf_dl_\"", boost::lexical_cast<std::string>(std::real(maxSurfElemSize))},
            {"\"_Eyinc_amp_\"", boost::lexical_cast<std::string>(std::real(pw_ey_amplitude))},
            {"\"_Bxinc_amp_\"", boost::lexical_cast<std::string>(std::real(pw_bx_amplitude))},
            {"\"_planewave_direction_x_\"", boost::lexical_cast<std::string>(std::real(pw_direction[0]))},
            {"\"_planewave_direction_y_\"", boost::lexical_cast<std::string>(std::real(pw_direction[1]))},
            {"\"_planewave_direction_z_\"", boost::lexical_cast<std::string>(std::real(pw_direction[2]))},
            {"\"_planewave_velocity_\"", boost::lexical_cast<std::string>(std::real(pw_velocity))},
            {"\"_gaussianplanewave_decay_rate_\"", boost::lexical_cast<std::string>(std::real(pw_decay_rate_t))},
            {"\"_planewave_t_center_\"", boost::lexical_cast<std::string>(std::real(pw_center_t))},
            {"\"_planewave_freq_\"", boost::lexical_cast<std::string>(std::real(pw_mod_freq))},
            {"\"_planewave_phase_\"", boost::lexical_cast<std::string>(std::real(pw_mod_phase))},
            {"\"_rect_width_\"", boost::lexical_cast<std::string>(std::real(pw_rect_width))},
            {"\"_rect_edge_width_\"", boost::lexical_cast<std::string>(std::real(pw_rect_edge_width))},
            {"\"_emitterwedge_angle_\"", boost::lexical_cast<std::string>(std::real(emitterWedgeAngle))},
            {"\"_emitterwedge_tip_radius_\"", boost::lexical_cast<std::string>(std::real(emitterWedgeTipRadius))},
            {"\"_emitterwedge_height_\"", boost::lexical_cast<std::string>(std::real(emitterWedgeHeight))},
            {"\"_emitterwedgetop_height_\"", boost::lexical_cast<std::string>(std::real(emitterWedgeTopHeight))},
            {"\"_emitterwedge_tip_y_\"", boost::lexical_cast<std::string>(std::real(emitterWedgeTipPosition[1]))},
            {"\"_emitterwedge_tip_z\"", boost::lexical_cast<std::string>(std::real(emitterWedgeTipPosition[2]))},
            {"\"_metalwedge_angle_\"", boost::lexical_cast<std::string>(std::real(metalWedgeAngle))},
            {"\"_metalwedge_tip_radius_\"", boost::lexical_cast<std::string>(std::real(metalWedgeTipRadius))},
            {"\"_metalwedge_height_\"", boost::lexical_cast<std::string>(std::real(metalWedgeHeight))},
            {"\"_metalwedge_tip_y_\"", boost::lexical_cast<std::string>(std::real(metalWedgeTipPosition[1]))},
            {"\"_metalwedge_tip_z\"", boost::lexical_cast<std::string>(std::real(metalWedgeTipPosition[2]))},
            {"\"_indy0Wedge_\"", boost::lexical_cast<std::string>(ind_wedge_y0)},
            {"\"_indy1Wedge_\"", boost::lexical_cast<std::string>(ind_wedge_y1)},
            {"\"_indy1Wedge_p1_\"", boost::lexical_cast<std::string>(ind_wedge_y1 + 1)},
            {"\"_indz0Wedge_\"", boost::lexical_cast<std::string>(ind_wedge_z0)},
            {"\"_indz1Wedge_\"", boost::lexical_cast<std::string>(ind_wedge_z1)},
            {"\"_indz1Wedge_p1_\"", boost::lexical_cast<std::string>(ind_wedge_z1 + 1)},
            {"\"_ny_wedge_\"", boost::lexical_cast<std::string>(ny_wedge)},
            {"\"_ny_wedge_p1_\"", boost::lexical_cast<std::string>(ny_wedge + 1)},
            {"\"_nz_wedge_\"", boost::lexical_cast<std::string>(nz_wedge)},
            {"\"_nz_wedge_p1_\"", boost::lexical_cast<std::string>(nz_wedge + 1)},
            {"\"_indz0Wedge_base_\"", boost::lexical_cast<std::string>(ind_wedge_base_z0)},
            {"\"_indz1Wedge_base_\"", boost::lexical_cast<std::string>(ind_wedge_base_z1)},
            {"\"_fdtd_unit_length_\"", boost::lexical_cast<std::string>(std::real(fdtd_unit_length))},
            {"\"_gamma_\"", boost::lexical_cast<std::string>(std::real(gamma))},
            {"\"_wp_sq_\"", boost::lexical_cast<std::string>(std::real(plasmaFrequency*plasmaFrequency))},
            {"\"_bunch_size_\"", boost::lexical_cast<std::string>(particle_bunch_size)},
            {"\"_save_rate_\"", boost::lexical_cast<std::string>(data_save_rate)},
            {"\"_save_rate_2_\"", boost::lexical_cast<std::string>(data_save_rate_2)},
            {"\"_save_rate_pml_\"", boost::lexical_cast<std::string>(data_save_rate_pml)},
            {"\"_nt_\"", boost::lexical_cast<std::string>(numOfTimeSamples)}
            };
    std::string fileName = "instructions/MaxwellYee2D_Metal_Wedge_electron_emitter_pureScatteredFDTD_partialGrid_GridCollection.json";

    ParameterExtractor::ReplaceStringsInFile(fileName,
                "instructions/processed/MaxwellYee2D_Metal_Wedge_electron_emitter_GridCollection_processed.json", str_replacewith);

    ParamFileTranslator fileTranslator("instructions/processed/MaxwellYee2D_Metal_Wedge_electron_emitter_GridCollection_processed.json");
    fileTranslator.Translate();


    std::string parametersFileName =
            std::string("data/2D/") + "params" + ".param";

    std::ofstream paramFileOut(parametersFileName.c_str(), std::ios::out | std::ios::binary);
    UtilityFunctions::WriteParamToFile<FPNumber>(paramFileOut, (FPNumber)units.ConvertFDTimeToSIUnits(dt), "dt_si");   // 0
    UtilityFunctions::WriteParamToFile<FPNumber>(paramFileOut, (FPNumber)units.ConvertFDLengthToSIUnits(dy), "dy_si");   // 1
    UtilityFunctions::WriteParamToFile<FPNumber>(paramFileOut, (FPNumber)units.ConvertFDLengthToSIUnits(dz), "dz_si");   // 2
    UtilityFunctions::WriteParamToFile<FPNumber>(paramFileOut, y0, "y0");
    UtilityFunctions::WriteParamToFile<FPNumber>(paramFileOut, y1, "y1");
    UtilityFunctions::WriteParamToFile<FPNumber>(paramFileOut, z0, "z0");
    UtilityFunctions::WriteParamToFile<FPNumber>(paramFileOut, z1, "z1");
    UtilityFunctions::WriteParamToFile<FPNumber>(paramFileOut, dt, "dt");
    UtilityFunctions::WriteParamToFile<FPNumber>(paramFileOut, dy, "dy");
    UtilityFunctions::WriteParamToFile<FPNumber>(paramFileOut, dz, "dz");
    UtilityFunctions::WriteParamToFile<FPNumber>(paramFileOut, fdtd_unit_length, "fdtd_unit_length");
    UtilityFunctions::WriteParamToFile<FPNumber>(paramFileOut, eField_FD_convertto_SI, "eField_FD_convertto_SI");
    UtilityFunctions::WriteParamToFile<FPNumber>(paramFileOut, data_save_rate, "data_save_rate");
    UtilityFunctions::WriteParamToFile<FPNumber>(paramFileOut, (FPNumber)units.ConvertFDFrequencyToSIUnits(1.0), "frequency_conversion_factor");   // 0
    UtilityFunctions::WriteParamToFile<FPNumber>(paramFileOut, ind_wedge_y0, "ind_wedge_y0");
    UtilityFunctions::WriteParamToFile<FPNumber>(paramFileOut, ind_wedge_y1, "ind_wedge_y1");
    UtilityFunctions::WriteParamToFile<FPNumber>(paramFileOut, ind_wedge_z0, "ind_wedge_z0");
    UtilityFunctions::WriteParamToFile<FPNumber>(paramFileOut, ind_wedge_z1, "ind_wedge_z1");
    UtilityFunctions::WriteParamToFile<FPNumber>(paramFileOut, pw_amplitude, "pw_amplitude");
    UtilityFunctions::WriteParamToFile<FPNumber>(paramFileOut, pw_center_t, "pw_center_t");
    UtilityFunctions::WriteParamToFile<FPNumber>(paramFileOut, pw_decay_rate_t, "pw_decay_rate_t");
    UtilityFunctions::WriteParamToFile<FPNumber>(paramFileOut, pw_mod_freq, "pw_mod_freq");
    UtilityFunctions::WriteParamToFile<FPNumber>(paramFileOut, pw_mod_phase, "pw_mod_phase");
    UtilityFunctions::WriteParamToFile<FPNumber>(paramFileOut, pw_velocity, "pw_velocity");
    UtilityFunctions::WriteParamToFile<FPNumber>(paramFileOut, pw_ey_amplitude, "pw_ey_amplitude");
    UtilityFunctions::WriteParamToFile<FPNumber>(paramFileOut, pw_bx_amplitude, "pw_bx_amplitude");

}




