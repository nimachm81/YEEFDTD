


void test_run_fdtd_1d_from_json() {
    FPNumber z0 = 0.0;
    FPNumber z1 = 10.0;
    std::size_t nz = 300;
    FPNumber dz = (z1 - z0)/(FPNumber)(nz);
    FPNumber stabilityFactor = 0.99;
    FPNumber dt = dz*stabilityFactor;
    FPNumber z_j = 5.0;
    std::size_t indJ = std::round(std::real((z_j - z0)/dz));
    std::size_t numOfTimeSamples = 600;

    std::unordered_map<std::string, std::string> str_replacewith{
            {"\"_z0_\"", boost::lexical_cast<std::string>(std::real(z0))},
            {"\"_z1_\"", boost::lexical_cast<std::string>(std::real(z1))},
            {"\"_nz_\"", boost::lexical_cast<std::string>(nz)},
            {"\"_dz_\"", boost::lexical_cast<std::string>(std::real(dz))},
            {"\"_dt_\"", boost::lexical_cast<std::string>(std::real(dt))},
            {"\"_dt_dz_\"", boost::lexical_cast<std::string>(std::real(dt/dz))},
            {"\"_m_dt_dz_\"", boost::lexical_cast<std::string>(std::real(-dt/dz))},
            {"\"_z_j_\"", boost::lexical_cast<std::string>(std::real(z_j))},
            {"\"_indJ_\"", boost::lexical_cast<std::string>(indJ)},
            {"\"_indJ_p1_\"", boost::lexical_cast<std::string>(indJ + 1)},
            {"\"_nt_\"", boost::lexical_cast<std::string>(numOfTimeSamples)}
            };
    ParameterExtractor::ReplaceStringsInFile("instructions/MaxwellYee1D.json",
                "instructions/MaxwellYee1D_processed.json", str_replacewith);

    ParamFileTranslator fileTranslator("instructions/MaxwellYee1D_processed.json");
    fileTranslator.Translate();
}




