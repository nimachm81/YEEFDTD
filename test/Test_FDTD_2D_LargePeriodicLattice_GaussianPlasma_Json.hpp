

#include "boost/lexical_cast.hpp"
#include "NumberTypes.h"
#include "ParameterExtractor.h"
#include "ParamFileTranslator.h"
#include "UtilityFunctions.hpp"

void test_run_fdtd_large_periodic_gaussian_plasma_2d_from_json() {

    FPNumber pitch = 127.0;
    FPNumber FWHM = 62.0;
    FPNumber eps_r = 11.7;


    FPNumber y0 = -5.0;
    FPNumber y1 = 5.0;
    FPNumber z0 = -8.0;
    FPNumber z1 = 5.0;
    std::size_t ny = 1000;
    std::size_t nz = 1300;
    FPNumber dy = (y1 - y0)/(FPNumber)(ny);
    FPNumber dz = (z1 - z0)/(FPNumber)(nz);
    FPNumber stabilityFactor = 0.95;
    FPNumber dt = (FPNumber)1.0/std::sqrt((FPNumber)1.0/(dy*dy) + (FPNumber)1.0/(dz*dz))*stabilityFactor;
    FPNumber z_j = -5.5;
    std::size_t indzJ = std::round(std::real((z_j - z0)/dz));
    FPNumber j_center_y = (y0 + y1)/(FPNumber)2.0;
    FPNumber j_decay_rate_y = 0.5;
    FPNumber j_center_t = 1.5;
    FPNumber jm_center_t = j_center_t + dt/(FPNumber)2.0*std::sqrt(eps_r);
    FPNumber jm_amplitude = std::sqrt(eps_r);
    FPNumber j_mod_freq = 1.0;
    FPNumber j_mod_phase = M_PI/2.0;


    std::size_t numOfTimeSamples = 2400;

    FPNumber gamma = 1.0e12/(3.0e8/(pitch*1.0e-6));
    FPNumber wp = 5.0e12/(3.0e8/(pitch*1.0e-6));
    std::cout << "wp : " << wp << std::endl;

    FPNumber wp2_decayrate_y = FWHMtoDecayRate(FWHM/pitch);
    FPNumber wp2_decayrate_z = FWHMtoDecayRate(FWHM/pitch);
    FPNumber wp2_center_y = 0.5;
    FPNumber wp2_center_z = 0.5;

    FPNumber wp2_mask_z0 = -5.0;

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
            {"\"_cube_dy_\"", boost::lexical_cast<std::string>(std::real(0.0))},
            {"\"_cube_dz_\"", boost::lexical_cast<std::string>(std::real(0.1))},
            {"\"_nt_\"", boost::lexical_cast<std::string>(numOfTimeSamples)}
            };
    ParameterExtractor::ReplaceStringsInFile("instructions/MaxwellYee2D_LargePeriodicGaussianPlasma.json",
                "instructions/MaxwellYee2D_LargePeriodicGaussianPlasma_processed.json", str_replacewith);

    ParamFileTranslator fileTranslator("instructions/MaxwellYee2D_LargePeriodicGaussianPlasma_processed.json");
    fileTranslator.Translate();
}


