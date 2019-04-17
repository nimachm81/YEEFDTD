


void test_run_fdtd_2d_nonuniform_GridCollection_from_json() {
    FPNumber y0 = -5.0;
    FPNumber y1 = 5.0;
    FPNumber z0 = -5.0;
    FPNumber z1 = 5.0;
    std::size_t ny = 400;
    std::size_t nz = 400;
    FPNumber dy = (y1 - y0)/(FPNumber)(ny);
    FPNumber dz = (z1 - z0)/(FPNumber)(nz);
    FPNumber stabilityFactor = 0.99;
    FPNumber dt = (FPNumber)1.0/std::sqrt((FPNumber)1.0/(dy*dy) + (FPNumber)1.0/(dz*dz))*stabilityFactor;
    FPNumber y_j = (y0 + y1)/(FPNumber)2.0;
    FPNumber z_j = (z0 + z1)/(FPNumber)2.0;
    std::size_t indyJ = std::round(std::real((y_j - y0)/dy));
    std::size_t indzJ = std::round(std::real((z_j - z0)/dz));

    FPNumber gr_y0 = -5.0;
    FPNumber gr_y1 = 5.0;
    FPNumber gr_z0 = 5.0;
    FPNumber gr_z1 = 15.0;
    std::size_t gr_ny = ny / 2;
    std::size_t gr_nz = nz / 2;
    FPNumber gr_dy = (gr_y1 - gr_y0)/(FPNumber)(gr_ny);
    FPNumber gr_dz = (gr_z1 - gr_z0)/(FPNumber)(gr_nz);
    FPNumber gr_dt = 2.0*dt;


    std::size_t numOfTimeSamples = 201;

    std::size_t gr_save_rate = 2;
    std::size_t save_rate = 2*gr_save_rate;

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
            {"\"_nz_m2_\"", boost::lexical_cast<std::string>(nz - 2)},
            {"\"_dy_\"", boost::lexical_cast<std::string>(std::real(dy))},
            {"\"_dz_\"", boost::lexical_cast<std::string>(std::real(dz))},
            {"\"_dt_\"", boost::lexical_cast<std::string>(std::real(dt))},
            {"\"_dt_dy_\"", boost::lexical_cast<std::string>(std::real(dt/dy))},
            {"\"_m_dt_dy_\"", boost::lexical_cast<std::string>(std::real(-dt/dy))},
            {"\"_dt_dz_\"", boost::lexical_cast<std::string>(std::real(dt/dz))},
            {"\"_m_dt_dz_\"", boost::lexical_cast<std::string>(std::real(-dt/dz))},
            {"\"_m_dt_dydz_\"", boost::lexical_cast<std::string>(std::real(-dt/(dy*dz)))},
            {"\"_indyJ_\"", boost::lexical_cast<std::string>(indyJ)},
            {"\"_indyJ_p1_\"", boost::lexical_cast<std::string>(indyJ + 1)},
            {"\"_indzJ_\"", boost::lexical_cast<std::string>(indzJ)},
            {"\"_indzJ_p1_\"", boost::lexical_cast<std::string>(indzJ + 1)},
            {"\"_gr_y0_\"", boost::lexical_cast<std::string>(std::real(gr_y0))},
            {"\"_gr_y1_\"", boost::lexical_cast<std::string>(std::real(gr_y1))},
            {"\"_gr_z0_\"", boost::lexical_cast<std::string>(std::real(gr_z0))},
            {"\"_gr_z1_\"", boost::lexical_cast<std::string>(std::real(gr_z1))},
            {"\"_gr_ny_\"", boost::lexical_cast<std::string>(gr_ny)},
            {"\"_gr_nz_\"", boost::lexical_cast<std::string>(gr_nz)},
            {"\"_gr_ny_p1_\"", boost::lexical_cast<std::string>(gr_ny + 1)},
            {"\"_gr_nz_p1_\"", boost::lexical_cast<std::string>(gr_nz + 1)},
            {"\"_gr_dy_\"", boost::lexical_cast<std::string>(std::real(gr_dy))},
            {"\"_gr_dz_\"", boost::lexical_cast<std::string>(std::real(gr_dz))},
            {"\"_gr_dt_\"", boost::lexical_cast<std::string>(std::real(gr_dt))},
            {"\"_gr_dt_dy_\"", boost::lexical_cast<std::string>(std::real(gr_dt/gr_dy))},
            {"\"_gr_m_dt_dy_\"", boost::lexical_cast<std::string>(std::real(-gr_dt/gr_dy))},
            {"\"_gr_dt_dz_\"", boost::lexical_cast<std::string>(std::real(gr_dt/gr_dz))},
            {"\"_gr_m_dt_dz_\"", boost::lexical_cast<std::string>(std::real(-gr_dt/gr_dz))},
            {"\"_gr_dt_dy_2_\"", boost::lexical_cast<std::string>(std::real(gr_dt/gr_dy/2.0))},
            {"\"_gr_m_dt_dy_2_\"", boost::lexical_cast<std::string>(std::real(-gr_dt/gr_dy/2.0))},
            {"\"_gr_dt_dz_2_\"", boost::lexical_cast<std::string>(std::real(gr_dt/gr_dz/2.0))},
            {"\"_gr_m_dt_dz_2_\"", boost::lexical_cast<std::string>(std::real(-gr_dt/gr_dz/2.0))},
            {"\"_nt_coarse_\"", boost::lexical_cast<std::string>(numOfTimeSamples)},
            {"\"_save_rate_\"", boost::lexical_cast<std::string>(save_rate)},
            {"\"_gr_save_rate_\"", boost::lexical_cast<std::string>(gr_save_rate)}
            };
    ParameterExtractor::ReplaceStringsInFile("instructions/MaxwellYee2D_nonuniform_GridCollection.json",
                "instructions/processed/MaxwellYee2D_nonuniform_GridCollection_processed.json", str_replacewith);

    ParamFileTranslator fileTranslator("instructions/processed/MaxwellYee2D_nonuniform_GridCollection_processed.json");
    fileTranslator.Translate();

    std::string parametersFileName =
            std::string("data/2D/") + "params" + ".param";

    std::ofstream paramFileOut(parametersFileName.c_str(), std::ios::out | std::ios::binary);
    UtilityFunctions::WriteParamToFile<FPNumber>(paramFileOut, y0, "y0");
    UtilityFunctions::WriteParamToFile<FPNumber>(paramFileOut, y1, "y1");
    UtilityFunctions::WriteParamToFile<FPNumber>(paramFileOut, z0, "z0");
    UtilityFunctions::WriteParamToFile<FPNumber>(paramFileOut, z1, "z1");
    UtilityFunctions::WriteParamToFile<FPNumber>(paramFileOut, gr_y0, "gr_y0");
    UtilityFunctions::WriteParamToFile<FPNumber>(paramFileOut, gr_y1, "gr_y1");
    UtilityFunctions::WriteParamToFile<FPNumber>(paramFileOut, gr_z0, "gr_z0");
    UtilityFunctions::WriteParamToFile<FPNumber>(paramFileOut, gr_z1, "gr_z1");

}



