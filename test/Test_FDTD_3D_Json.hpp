

void test_run_fdtd_3d_from_json() {
    FPNumber x0 = 0;
    FPNumber x1 = 10.0;
    FPNumber y0 = 0;
    FPNumber y1 = 10.0;
    FPNumber z0 = 0.0;
    FPNumber z1 = 10.0;
    std::size_t nx = 200;
    std::size_t ny = 200;
    std::size_t nz = 200;
    FPNumber dx = (x1 - x0)/(FPNumber)(nx);
    FPNumber dy = (y1 - y0)/(FPNumber)(ny);
    FPNumber dz = (z1 - z0)/(FPNumber)(nz);
    FPNumber stabilityFactor = 0.99;
    FPNumber dt = (FPNumber)1.0/std::sqrt((FPNumber)1.0/(dx*dx) + (FPNumber)1.0/(dy*dy) + (FPNumber)1.0/(dz*dz))*stabilityFactor;
    FPNumber x_j = (x0 + x1)/(FPNumber)2.0;
    FPNumber y_j = (y0 + y1)/(FPNumber)2.0;
    FPNumber z_j = (z0 + z1)/(FPNumber)2.0;
    std::size_t indxJ = std::round(std::real((x_j - x0)/dx));
    std::size_t indyJ = std::round(std::real((y_j - y0)/dy));
    std::size_t indzJ = std::round(std::real((z_j - z0)/dz));
    std::size_t numOfTimeSamples = 150;

    std::unordered_map<std::string, std::string> str_replacewith{
            {"\"_x0_\"", boost::lexical_cast<std::string>(std::real(x0))},
            {"\"_x1_\"", boost::lexical_cast<std::string>(std::real(x1))},
            {"\"_y0_\"", boost::lexical_cast<std::string>(std::real(y0))},
            {"\"_y1_\"", boost::lexical_cast<std::string>(std::real(y1))},
            {"\"_z0_\"", boost::lexical_cast<std::string>(std::real(z0))},
            {"\"_z1_\"", boost::lexical_cast<std::string>(std::real(z1))},
            {"\"_nx_\"", boost::lexical_cast<std::string>(nx)},
            {"\"_ny_\"", boost::lexical_cast<std::string>(ny)},
            {"\"_nz_\"", boost::lexical_cast<std::string>(nz)},
            {"\"_nx_p1_\"", boost::lexical_cast<std::string>(nx + 1)},
            {"\"_ny_p1_\"", boost::lexical_cast<std::string>(ny + 1)},
            {"\"_nz_p1_\"", boost::lexical_cast<std::string>(nz + 1)},
            {"\"_dx_\"", boost::lexical_cast<std::string>(std::real(dx))},
            {"\"_dy_\"", boost::lexical_cast<std::string>(std::real(dy))},
            {"\"_dz_\"", boost::lexical_cast<std::string>(std::real(dz))},
            {"\"_dt_\"", boost::lexical_cast<std::string>(std::real(dt))},
            {"\"_dt_dx_\"", boost::lexical_cast<std::string>(std::real(dt/dx))},
            {"\"_m_dt_dx_\"", boost::lexical_cast<std::string>(std::real(-dt/dx))},
            {"\"_dt_dy_\"", boost::lexical_cast<std::string>(std::real(dt/dy))},
            {"\"_m_dt_dy_\"", boost::lexical_cast<std::string>(std::real(-dt/dy))},
            {"\"_dt_dz_\"", boost::lexical_cast<std::string>(std::real(dt/dz))},
            {"\"_m_dt_dz_\"", boost::lexical_cast<std::string>(std::real(-dt/dz))},
            {"\"_m_dt_dxdydz_\"", boost::lexical_cast<std::string>(std::real(-dt/(dx*dy*dz)))},
            {"\"_x_j_\"", boost::lexical_cast<std::string>(std::real(x_j))},
            {"\"_y_j_\"", boost::lexical_cast<std::string>(std::real(y_j))},
            {"\"_z_j_\"", boost::lexical_cast<std::string>(std::real(z_j))},
            {"\"_indxJ_\"", boost::lexical_cast<std::string>(indxJ)},
            {"\"_indxJ_p1_\"", boost::lexical_cast<std::string>(indxJ + 1)},
            {"\"_indyJ_\"", boost::lexical_cast<std::string>(indyJ)},
            {"\"_indyJ_p1_\"", boost::lexical_cast<std::string>(indyJ + 1)},
            {"\"_indzJ_\"", boost::lexical_cast<std::string>(indzJ)},
            {"\"_indzJ_p1_\"", boost::lexical_cast<std::string>(indzJ + 1)},
            {"\"_nt_\"", boost::lexical_cast<std::string>(numOfTimeSamples)},
            {"\"_mod_phase_\"", boost::lexical_cast<std::string>(M_PI/2.0)}
            };
    ParameterExtractor::ReplaceStringsInFile("instructions/MaxwellYee3D.json",
                "instructions/MaxwellYee3D_processed.json", str_replacewith);

    ParamFileTranslator fileTranslator("instructions/MaxwellYee3D_processed.json");
    fileTranslator.Translate();
}



