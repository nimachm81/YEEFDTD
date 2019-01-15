


void test_run_fdtd_2d_monopole_charge_from_json() {
    FPNumber y0 = -3.0;
    FPNumber y1 = 3.0;
    FPNumber z0 = -3.0;
    FPNumber z1 = 3.0;
    std::size_t ny = 1000;
    std::size_t nz = 1000;
    FPNumber dy = (y1 - y0)/(FPNumber)(ny);
    FPNumber dz = (z1 - z0)/(FPNumber)(nz);
    FPNumber stabilityFactor = 0.99;
    FPNumber dt = (FPNumber)1.0/std::sqrt((FPNumber)1.0/(dy*dy) + (FPNumber)1.0/(dz*dz))*stabilityFactor;
    FPNumber y0_j = (FPNumber)1.0*y0 + (FPNumber)0.0*y1;
    FPNumber y1_j = (y0 + y1)/(FPNumber)2.0;
    FPNumber z_j = (z0 + z1)/(FPNumber)2.0;
    std::size_t indy0J = std::round(std::real((y0_j - y0)/dy));
    std::size_t indy1J = std::round(std::real((y1_j - y0)/dy));
    std::size_t indzJ = std::round(std::real((z_j - z0)/dz));
    std::size_t numOfTimeSamples = 401;

    std::unordered_map<std::string, std::string> str_replacewith{
            {"\"_y0_\"", boost::lexical_cast<std::string>(std::real(y0))},
            {"\"_y1_\"", boost::lexical_cast<std::string>(std::real(y1))},
            {"\"_z0_\"", boost::lexical_cast<std::string>(std::real(z0))},
            {"\"_z1_\"", boost::lexical_cast<std::string>(std::real(z1))},
            {"\"_ny_\"", boost::lexical_cast<std::string>(ny)},
            {"\"_nz_\"", boost::lexical_cast<std::string>(nz)},
            {"\"_ny_p1_\"", boost::lexical_cast<std::string>(ny + 1)},
            {"\"_nz_p1_\"", boost::lexical_cast<std::string>(nz + 1)},
            {"\"_dy_\"", boost::lexical_cast<std::string>(std::real(dy))},
            {"\"_dz_\"", boost::lexical_cast<std::string>(std::real(dz))},
            {"\"_dt_\"", boost::lexical_cast<std::string>(std::real(dt))},
            {"\"_dt_dy_\"", boost::lexical_cast<std::string>(std::real(dt/dy))},
            {"\"_m_dt_dy_\"", boost::lexical_cast<std::string>(std::real(-dt/dy))},
            {"\"_dt_dz_\"", boost::lexical_cast<std::string>(std::real(dt/dz))},
            {"\"_m_dt_dz_\"", boost::lexical_cast<std::string>(std::real(-dt/dz))},
            {"\"_indy0J_\"", boost::lexical_cast<std::string>(indy0J)},
            {"\"_indy1J_\"", boost::lexical_cast<std::string>(indy1J)},
            {"\"_ny_J_\"", boost::lexical_cast<std::string>(indy1J - indy0J)},
            {"\"_indzJ_\"", boost::lexical_cast<std::string>(indzJ)},
            {"\"_indzJ_p1_\"", boost::lexical_cast<std::string>(indzJ + 1)},
            {"\"_nt_\"", boost::lexical_cast<std::string>(numOfTimeSamples)}
            };
    ParameterExtractor::ReplaceStringsInFile("instructions/MaxwellYee2D_monopoleCharge.json",
                "instructions/processed/MaxwellYee2D_monopoleCharge_processed.json", str_replacewith);

    ParamFileTranslator fileTranslator("instructions/processed/MaxwellYee2D_monopoleCharge_processed.json");
    fileTranslator.Translate();
}


