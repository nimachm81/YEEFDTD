

#include "boost/lexical_cast.hpp"
#include "NumberTypes.h"
#include "ParameterExtractor.h"
#include "ParamFileTranslator.h"


void test_run_fdtd_1d_nonuniform_collection_from_json() {
    FPNumber grid0_z0 = 0.0;
    FPNumber grid0_z1 = 10.0;
    std::size_t grid0_nz = 400;
    FPNumber grid0_dz = (grid0_z1 - grid0_z0)/(FPNumber)(grid0_nz);
    FPNumber stabilityFactor = 0.99;
    FPNumber grid0_dt = grid0_dz*stabilityFactor;
    FPNumber grid0_z_j = 5.0;
    std::size_t grid0_indJ = std::round(std::real((grid0_z_j - grid0_z0)/grid0_dz));

    FPNumber grid1_z0 = 10.0;
    FPNumber grid1_z1 = 20.0;
    std::size_t grid1_nz = grid0_nz/2;
    FPNumber grid1_dz = (grid1_z1 - grid1_z0)/(FPNumber)(grid1_nz);
    FPNumber grid1_dt = grid1_dz*stabilityFactor;

    std::size_t numOfTimeSamples = 1000;

    std::size_t grid1_saveRate = 1;
    std::size_t grid0_saveRate = 2*grid1_saveRate;

    std::unordered_map<std::string, std::string> str_replacewith{
            {"\"_g0_z0_\"", boost::lexical_cast<std::string>(std::real(grid0_z0))},
            {"\"_g0_z1_\"", boost::lexical_cast<std::string>(std::real(grid0_z1))},
            {"\"_g0_nz_\"", boost::lexical_cast<std::string>(grid0_nz)},
            {"\"_g0_nz_m1_\"", boost::lexical_cast<std::string>(grid0_nz - 1)},
            {"\"_g0_nz_m2_\"", boost::lexical_cast<std::string>(grid0_nz - 2)},
            {"\"_g0_nz_p1_\"", boost::lexical_cast<std::string>(grid0_nz + 1)},
            {"\"_g0_dz_\"", boost::lexical_cast<std::string>(std::real(grid0_dz))},
            {"\"_g0_dt_\"", boost::lexical_cast<std::string>(std::real(grid0_dt))},
            {"\"_g0_dt_dz_\"", boost::lexical_cast<std::string>(std::real(grid0_dt/grid0_dz))},
            {"\"_g0_m_dt_dz_\"", boost::lexical_cast<std::string>(std::real(-grid0_dt/grid0_dz))},
            {"\"_g0_indJ_\"", boost::lexical_cast<std::string>(grid0_indJ)},
            {"\"_g0_indJ_p1_\"", boost::lexical_cast<std::string>(grid0_indJ + 1)},
            {"\"_g1_z0_\"", boost::lexical_cast<std::string>(std::real(grid1_z0))},
            {"\"_g1_z1_\"", boost::lexical_cast<std::string>(std::real(grid1_z1))},
            {"\"_g1_nz_\"", boost::lexical_cast<std::string>(grid1_nz)},
            {"\"_g1_nz_m1_\"", boost::lexical_cast<std::string>(grid1_nz - 1)},
            {"\"_g1_dz_\"", boost::lexical_cast<std::string>(std::real(grid1_dz))},
            {"\"_g1_dt_\"", boost::lexical_cast<std::string>(std::real(grid1_dt))},
            {"\"_g1_dt_dz_\"", boost::lexical_cast<std::string>(std::real(grid1_dt/grid1_dz))},
            {"\"_g1_m_dt_dz_\"", boost::lexical_cast<std::string>(std::real(-grid1_dt/grid1_dz))},
            {"\"_g1_dt_dz_2_\"", boost::lexical_cast<std::string>(std::real(grid1_dt/grid1_dz/2.0))},
            {"\"_g1_m_dt_dz_2_\"", boost::lexical_cast<std::string>(std::real(-grid1_dt/grid1_dz/2.0))},
            {"\"_gall_nt_\"", boost::lexical_cast<std::string>(numOfTimeSamples)},
            {"\"_g0_saverate_\"", boost::lexical_cast<std::string>(grid0_saveRate)},
            {"\"_g1_saverate_\"", boost::lexical_cast<std::string>(grid1_saveRate)}
            };
    ParameterExtractor::ReplaceStringsInFile("instructions/MaxwellYee1D_nonuniform_GridCollection.json",
                "instructions/processed/MaxwellYee1D_nonuniform_GridCollection_processed.json", str_replacewith);

    ParamFileTranslator fileTranslator("instructions/processed/MaxwellYee1D_nonuniform_GridCollection_processed.json");
    fileTranslator.Translate();
}


