

void test_run_fdtd_plasma_1d_from_json() {
    FPNumber z0 = 0.0;
    FPNumber z1 = 10.0;
    std::size_t nz = 1000;
    FPNumber dz = (z1 - z0)/(FPNumber)(nz);
    FPNumber stabilityFactor = 0.99;
    FPNumber dt = dz*stabilityFactor;
    FPNumber z_j = 5.0;
    std::size_t indJ = std::round(std::real((z_j - z0)/dz));
    std::size_t numOfTimeSamples = 3000;

    FPNumber cube_z0 = 7.0;
    FPNumber cube_z1 = 10.0;
    FPNumber cube_dz = 0.0;

    FPNumber gamma = 0.0;
    FPNumber wp = 2.0*M_PI*3.0;

    std::unordered_map<std::string, std::string> str_replacewith{
            {"\"_z0_\"", boost::lexical_cast<std::string>(std::real(z0))},
            {"\"_z1_\"", boost::lexical_cast<std::string>(std::real(z1))},
            {"\"_nz_\"", boost::lexical_cast<std::string>(nz)},
            {"\"_dz_\"", boost::lexical_cast<std::string>(std::real(dz))},
            {"\"_dt_\"", boost::lexical_cast<std::string>(std::real(dt))},
            {"\"_m_dt_\"", boost::lexical_cast<std::string>(std::real(-dt))},
            {"\"_dt_dz_\"", boost::lexical_cast<std::string>(std::real(dt/dz))},
            {"\"_m_dt_dz_\"", boost::lexical_cast<std::string>(std::real(-dt/dz))},
            {"\"_z_j_\"", boost::lexical_cast<std::string>(std::real(z_j))},
            {"\"_indJ_\"", boost::lexical_cast<std::string>(indJ)},
            {"\"_indJ_p1_\"", boost::lexical_cast<std::string>(indJ + 1)},
            {"\"_cube_z0_\"", boost::lexical_cast<std::string>(std::real(cube_z0))},
            {"\"_cube_z1_\"", boost::lexical_cast<std::string>(std::real(cube_z1))},
            {"\"_cube_dz_\"", boost::lexical_cast<std::string>(std::real(cube_dz))},
            {"\"_wp_sq_\"", boost::lexical_cast<std::string>(std::real(wp*wp))},
            {"\"_m_dt_tau_\"", boost::lexical_cast<std::string>(std::real(-dt*gamma))},
            {"\"_nt_\"", boost::lexical_cast<std::string>(numOfTimeSamples)}
            };
    ParameterExtractor::ReplaceStringsInFile("instructions/MaxwellYee1D_plasma.json",
                "instructions/MaxwellYee1D_plasma_processed.json", str_replacewith);

    ParamFileTranslator fileTranslator("instructions/MaxwellYee1D_plasma_processed.json");
    fileTranslator.Translate();
}



