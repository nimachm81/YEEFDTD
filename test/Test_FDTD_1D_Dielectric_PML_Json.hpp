

void test_run_fdtd_dielectric_pml_1d_from_json() {
    FPNumber z0 = 0.0;
    FPNumber z1 = 15.0;
    std::size_t nz = 600;
    FPNumber dz = (z1 - z0)/(FPNumber)(nz);
    FPNumber stabilityFactor = 0.95;
    FPNumber dt = dz*stabilityFactor;
    FPNumber z_j = 5.0;
    std::size_t indJ = std::round(std::real((z_j - z0)/dz));
    std::size_t numOfTimeSamples = 800;

    FPNumber cube_z0 = 9.0;
    FPNumber cube_z1 = 15.0;
    FPNumber cube_dz = 0.0;

    FPNumber epsilon = 4.0;

    FPNumber pml_z0 = 3.0;
    FPNumber pml_z1 = 12.0;
    FPNumber pml_dz = 1.0;

    FPNumber sig_e = 1.0;
    FPNumber sig_h = 1.0;


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
            {"\"_eps_inv_\"", boost::lexical_cast<std::string>(std::real((FPNumber)1.0/epsilon))},
            {"\"_pml_z0_\"", boost::lexical_cast<std::string>(std::real(pml_z0))},
            {"\"_pml_z1_\"", boost::lexical_cast<std::string>(std::real(pml_z1))},
            {"\"_pml_dz_\"", boost::lexical_cast<std::string>(std::real(pml_dz))},
            {"\"_sig_e_\"", boost::lexical_cast<std::string>(std::real(sig_e))},
            {"\"_sig_h_\"", boost::lexical_cast<std::string>(std::real(sig_h))},
            {"\"_nt_\"", boost::lexical_cast<std::string>(numOfTimeSamples)}
            };
    ParameterExtractor::ReplaceStringsInFile("instructions/MaxwellYee1D_dielectric_pml.json",
                "instructions/MaxwellYee1D_dielectric_pml_processed.json", str_replacewith);

    ParamFileTranslator fileTranslator("instructions/MaxwellYee1D_dielectric_pml_processed.json");
    fileTranslator.Translate();
}



