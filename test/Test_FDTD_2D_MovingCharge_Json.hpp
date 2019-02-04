


void test_run_fdtd_2d_moving_charge_from_json() {
    FPNumber y0 = -3.0;
    FPNumber y1 = 3.0;
    FPNumber z0 = -3.0;
    FPNumber z1 = 3.0;
    std::size_t ny = 100;
    std::size_t nz = 100;
    FPNumber dy = (y1 - y0)/(FPNumber)(ny);
    FPNumber dz = (z1 - z0)/(FPNumber)(nz);
    FPNumber stabilityFactor = 0.99;
    FPNumber dt = (FPNumber)1.0/std::sqrt((FPNumber)1.0/(dy*dy) + (FPNumber)1.0/(dz*dz))*stabilityFactor;

    FPNumber epsilon_r = 9.0;

    FPNumber q = 1.0;//1.602e-19;
    FPNumber m = 1.0e4;
    FPNumber y0_q = (FPNumber)0.5*y0 + (FPNumber)0.5*y1;
    FPNumber z0_q = z0 - 2.0*dz;
    FPNumber vy0_q = 0.0;
    FPNumber vz0_q = 0.5;

    std::size_t numOfTimeSamples = 501;

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
            {"\"_dt_dy_eps_\"", boost::lexical_cast<std::string>(std::real(dt/dy/epsilon_r))},
            {"\"_m_dt_dy_eps_\"", boost::lexical_cast<std::string>(std::real(-dt/dy/epsilon_r))},
            {"\"_dt_dz_eps_\"", boost::lexical_cast<std::string>(std::real(dt/dz/epsilon_r))},
            {"\"_m_dt_dz_eps_\"", boost::lexical_cast<std::string>(std::real(-dt/dz/epsilon_r))},
            {"\"_p0_q_\"", boost::lexical_cast<std::string>(std::real(q))},
            {"\"_p0_m_\"", boost::lexical_cast<std::string>(std::real(m))},
            {"\"_p0_r0_y_\"", boost::lexical_cast<std::string>(std::real(y0_q))},
            {"\"_p0_r0_z_\"", boost::lexical_cast<std::string>(std::real(z0_q))},
            {"\"_p0_v0_y_\"", boost::lexical_cast<std::string>(std::real(vy0_q))},
            {"\"_p0_v0_z_\"", boost::lexical_cast<std::string>(std::real(vz0_q))},
            {"\"_nt_\"", boost::lexical_cast<std::string>(numOfTimeSamples)}
            };
    ParameterExtractor::ReplaceStringsInFile("instructions/MaxwellYee2D_MovingCharge.json",
                "instructions/processed/MaxwellYee2D_MovingCharge_processed.json", str_replacewith);

    ParamFileTranslator fileTranslator("instructions/processed/MaxwellYee2D_MovingCharge_processed.json");
    fileTranslator.Translate();
}


