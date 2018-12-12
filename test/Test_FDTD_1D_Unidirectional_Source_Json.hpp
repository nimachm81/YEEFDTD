

#include "boost/lexical_cast.hpp"
#include "NumberTypes.h"
#include "ParameterExtractor.h"
#include "ParamFileTranslator.h"


void test_run_fdtd_1d_unidirectional_source_from_json() {

    FPNumber eps_r = 11.7;

    FPNumber z0 = 0.0;
    FPNumber z1 = 10.0;
    std::size_t nz = 800;
    FPNumber dz = (z1 - z0)/(FPNumber)(nz);
    FPNumber stabilityFactor = 0.99;
    FPNumber dt = dz*stabilityFactor;
    FPNumber z_j = 5.0;
    std::size_t indJ = std::round(std::real((z_j - z0)/dz));
    std::size_t numOfTimeSamples = 200;

    FPNumber je_t_center = 1.5;
    FPNumber jm_t_center = je_t_center + dt/(FPNumber)2.0*std::sqrt(eps_r);
    FPNumber j_mod_freq = 1.0;
    FPNumber j_mod_phase = M_PI/2.0;

    FPNumber jm_amplitude = std::sqrt(eps_r);

    std::unordered_map<std::string, std::string> str_replacewith{
            {"\"_z0_\"", boost::lexical_cast<std::string>(std::real(z0))},
            {"\"_z1_\"", boost::lexical_cast<std::string>(std::real(z1))},
            {"\"_nz_\"", boost::lexical_cast<std::string>(nz)},
            {"\"_dz_\"", boost::lexical_cast<std::string>(std::real(dz))},
            {"\"_dt_\"", boost::lexical_cast<std::string>(std::real(dt))},
            {"\"_dt_dz_\"", boost::lexical_cast<std::string>(std::real(dt/dz))},
            {"\"_m_dt_dz_\"", boost::lexical_cast<std::string>(std::real(-dt/dz))},
            {"\"_dt_dz_eps_\"", boost::lexical_cast<std::string>(std::real(dt/dz/eps_r))},
            {"\"_m_dt_dz_eps_\"", boost::lexical_cast<std::string>(std::real(-dt/dz/eps_r))},
            {"\"_Je_t_center_\"", boost::lexical_cast<std::string>(std::real(je_t_center))},
            {"\"_Jm_t_center_\"", boost::lexical_cast<std::string>(std::real(jm_t_center))},
            {"\"_Jm_amp_\"", boost::lexical_cast<std::string>(std::real(jm_amplitude))},
            {"\"_J_mod_freq_\"", boost::lexical_cast<std::string>(std::real(j_mod_freq))},
            {"\"_J_mod_phase_\"", boost::lexical_cast<std::string>(std::real(j_mod_phase))},
            {"\"_indJ_\"", boost::lexical_cast<std::string>(indJ)},
            {"\"_indJ_p1_\"", boost::lexical_cast<std::string>(indJ + 1)},
            {"\"_indJ_m1_\"", boost::lexical_cast<std::string>(indJ - 1)},
            {"\"_nt_\"", boost::lexical_cast<std::string>(numOfTimeSamples)}
            };
    ParameterExtractor::ReplaceStringsInFile("instructions/MaxwellYee1D_UnidirectionalSource.json",
                "instructions/MaxwellYee1D_UnidirectionalSource_processed.json", str_replacewith);

    ParamFileTranslator fileTranslator("instructions/MaxwellYee1D_UnidirectionalSource_processed.json");
    fileTranslator.Translate();
}




