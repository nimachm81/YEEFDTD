
#include "boost/lexical_cast.hpp"
#include "NumberTypes.h"
#include "ParameterExtractor.h"
#include "ParamFileTranslator.h"
#include "UtilityFunctions.hpp"

void test_run_fdtd_gaussian_plasma_2d_periodic_square_lattice_from_json(
            double ky_kymax = 0.0,      // ky/ky_max
            double kz_kzmax = 0.5,     // kz_kz_max
            std::string fileNamePostfix = "-GX-",
            int fileInd = 0
            ) {
    assert(typeid(FPNumber) == typeid(std::complex<double>) || typeid(FPNumber) == typeid(std::complex<float>));

    FPNumber pitch = 127.0;
    FPNumber FWHM = 62.0;
    FPNumber eps_r = 11.7;

    FPNumber y0 = -0.5;
    FPNumber y1 = 0.5;
    FPNumber z0 = -0.5;
    FPNumber z1 = 0.5;
    std::size_t ny = 110;
    std::size_t nz = 110;
    FPNumber dy = (y1 - y0)/(FPNumber)(ny);
    FPNumber dz = (z1 - z0)/(FPNumber)(nz);
    FPNumber stabilityFactor = 0.99;
    FPNumber dt = (FPNumber)1.0/std::sqrt((FPNumber)1.0/(dy*dy) + (FPNumber)1.0/(dz*dz))*stabilityFactor;
    FPNumber y_j = (FPNumber)0.5*y0 + (FPNumber)0.5*y1;
    FPNumber z_j = (FPNumber)0.9*z0 + (FPNumber)0.1*z1;
    std::size_t indyJ = std::round(std::real((y_j - y0)/dy));
    std::size_t indzJ = std::round(std::real((z_j - z0)/dz));
    std::size_t numOfTimeSamples = 1500;

    FPNumber gamma = 0.0;
    FPNumber wp = 5.0*1.0e12/(3.0e8/(pitch*1.0e-6));
    std::cout << "wp : " << wp << std::endl;

    FPNumber wp2_decayrate_y = FWHMtoDecayRate(FWHM/pitch);
    FPNumber wp2_decayrate_z = FWHMtoDecayRate(FWHM/pitch);
    FPNumber wp2_center_y = 0.0;
    FPNumber wp2_center_z = 0.0;
    std::cout << "wp2_decayrate_y : " << wp2_decayrate_y << std::endl;
    std::cout << "wp2_decayrate_z : " << wp2_decayrate_z << std::endl;

    FPNumber ky_max = (FPNumber)(2.0*M_PI)/(y1 - y0);
    FPNumber kz_max = (FPNumber)(2.0*M_PI)/(z1 - z0);
    FPNumber k_y = (FPNumber)ky_kymax*ky_max;
    FPNumber k_z = (FPNumber)kz_kzmax*kz_max;


    FPNumber _j = (FPNumber)(DemotableComplex<double>(0.0, 1.0));

    FPNumber _m_dt_dy_exp_jky_y10_eps_ = -dt/dy*std::exp(_j*k_y*(y1 - y0))/eps_r;
    FPNumber _dt_dz_exp_jkz_z10_eps_  = dt/dz*std::exp(_j*k_z*(z1 - z0))/eps_r;
    FPNumber _exp_mjky_y10_ = std::exp(-_j*k_y*(y1 - y0));
    FPNumber _exp_mjkz_z10_ = std::exp(-_j*k_z*(z1 - z0));

    std::string ExFilename = std::string("\"") +
                             "GaussianPlasmaPBC/E-x" +
                             fileNamePostfix +
                             boost::lexical_cast<std::string>(fileInd) +
                             "\"";
    std::string HyFilename = std::string("\"") +
                             "GaussianPlasmaPBC/H-y" +
                             fileNamePostfix +
                             boost::lexical_cast<std::string>(fileInd) +
                             "\"";
    std::string Wp2Filename = std::string("\"") +
                              "GaussianPlasmaPBC/Wp2-x" +
                              fileNamePostfix +
                              boost::lexical_cast<std::string>(fileInd) +
                              "\"";

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
            {"\"_dt_dy_\"", boost::lexical_cast<std::string>(std::real(dt/dy))},
            {"\"_m_dt_dy_\"", boost::lexical_cast<std::string>(std::real(-dt/dy))},
            {"\"_dt_dz_\"", boost::lexical_cast<std::string>(std::real(dt/dz))},
            {"\"_m_dt_dz_\"", boost::lexical_cast<std::string>(std::real(-dt/dz))},
            {"\"_dt_dy_eps_\"", boost::lexical_cast<std::string>(std::real(dt/dy/eps_r))},
            {"\"_m_dt_dy_eps_\"", boost::lexical_cast<std::string>(std::real(-dt/dy/eps_r))},
            {"\"_dt_dz_eps_\"", boost::lexical_cast<std::string>(std::real(dt/dz/eps_r))},
            {"\"_m_dt_dz_eps_\"", boost::lexical_cast<std::string>(std::real(-dt/dz/eps_r))},
            {"\"_m_dt_dydz_\"", boost::lexical_cast<std::string>(std::real(-dt/(dy*dz)))},
            {"\"_y_j_\"", boost::lexical_cast<std::string>(std::real(y_j))},
            {"\"_z_j_\"", boost::lexical_cast<std::string>(std::real(z_j))},
            {"\"_indyJ_\"", boost::lexical_cast<std::string>(indyJ)},
            {"\"_indyJ_p1_\"", boost::lexical_cast<std::string>(indyJ + 1)},
            {"\"_indzJ_\"", boost::lexical_cast<std::string>(indzJ)},
            {"\"_indzJ_p1_\"", boost::lexical_cast<std::string>(indzJ + 1)},
            {"\"_m_dt_tau_\"", boost::lexical_cast<std::string>(std::real(-dt*gamma))},
            {"\"_wp_sq_\"", boost::lexical_cast<std::string>(std::real(wp*wp))},
            {"\"_wp_sq_center_y_\"", boost::lexical_cast<std::string>(std::real(wp2_center_y))},
            {"\"_wp_sq_center_z_\"", boost::lexical_cast<std::string>(std::real(wp2_center_z))},
            {"\"_wp_sq_decayrate_y_\"", boost::lexical_cast<std::string>(std::real(wp2_decayrate_y))},
            {"\"_wp_sq_decayrate_z_\"", boost::lexical_cast<std::string>(std::real(wp2_decayrate_z))},
            {"\"_m_dt_dy_exp_jky_y10_eps_\"", CastComplexToJSONString(_m_dt_dy_exp_jky_y10_eps_)},
            {"\"_dt_dz_exp_jkz_z10_eps_\"", CastComplexToJSONString(_dt_dz_exp_jkz_z10_eps_)},
            {"\"_exp_mjky_y10_\"", CastComplexToJSONString(_exp_mjky_y10_)},
            {"\"_exp_mjkz_z10_\"", CastComplexToJSONString(_exp_mjkz_z10_)},
            {"\"_nt_\"", boost::lexical_cast<std::string>(numOfTimeSamples)},
            {"\"_Ex_filename_\"", boost::lexical_cast<std::string>(ExFilename)},
            {"\"_Hy_filename_\"", boost::lexical_cast<std::string>(HyFilename)},
            {"\"_Wp2_filename_\"", boost::lexical_cast<std::string>(Wp2Filename)}
            };
    ParameterExtractor::ReplaceStringsInFile("instructions/MaxwellYee2D_GaussianPlasma_PBC.json",
                "instructions/MaxwellYee2D_GaussianPlasma_PBC_processed.json", str_replacewith);

    ParamFileTranslator fileTranslator("instructions/MaxwellYee2D_GaussianPlasma_PBC_processed.json");
    fileTranslator.Translate();
}

void test_run_fdtd_gaussian_plasma_2d_periodic_square_lattice_sweep_over_BZ(int n_pts = 100) {
    for(int i = 0; i < n_pts; ++i) {
        double ky_kymax = 0.0;
        double kz_kzmax = (double)i/n_pts;
        test_run_fdtd_gaussian_plasma_2d_periodic_square_lattice_from_json(ky_kymax, kz_kzmax, "-GX-", i);
    }
}

