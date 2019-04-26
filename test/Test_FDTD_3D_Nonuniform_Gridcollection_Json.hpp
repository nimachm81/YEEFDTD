


void test_run_fdtd_3d_nonuniform_gridcollection_from_json() {

    std::array<FPNumber, 3> box_in_r0 = {-1.0, -1.0, -1.0};
    std::array<FPNumber, 3> box_in_r1 = {+1.0, +1.0, +1.0};

    std::array<FPNumber, 3> box_out_r0 = {-3.0, -3.0, -3.0};
    std::array<FPNumber, 3> box_out_r1 = {+3.0, +3.0, +3.0};

    std::size_t n_xyz_in = 100;

    FPNumber x0 = box_in_r0[0];
    FPNumber x1 = box_in_r1[0];
    FPNumber y0 = box_in_r0[1];
    FPNumber y1 = box_in_r1[1];
    FPNumber z0 = box_in_r0[2];
    FPNumber z1 = box_in_r1[2];
    std::size_t nx = n_xyz_in;
    std::size_t ny = n_xyz_in;
    std::size_t nz = n_xyz_in;
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
    if(indxJ % 2 == 1) {
        indxJ += 1;
    }

    FPNumber gr_x0 = box_in_r0[0];
    FPNumber gr_x1 = box_in_r1[0];
    FPNumber gr_y0 = box_in_r0[1];
    FPNumber gr_y1 = box_in_r1[1];
    FPNumber gr_z0 = box_in_r1[2];
    FPNumber gr_z1 = box_out_r1[2];
    std::size_t gr_nx = n_xyz_in/2;
    std::size_t gr_ny = n_xyz_in/2;
    std::size_t gr_nz = n_xyz_in/2;
    FPNumber gr_dx = (gr_x1 - gr_x0)/(FPNumber)(gr_nx);
    FPNumber gr_dy = (gr_y1 - gr_y0)/(FPNumber)(gr_ny);
    FPNumber gr_dz = (gr_z1 - gr_z0)/(FPNumber)(gr_nz);
    FPNumber gr_dt = 2.0*dt;

    FPNumber gl_x0 = box_in_r0[0];
    FPNumber gl_x1 = box_in_r1[0];
    FPNumber gl_y0 = box_in_r0[1];
    FPNumber gl_y1 = box_in_r1[1];
    FPNumber gl_z0 = box_out_r0[2];
    FPNumber gl_z1 = box_in_r0[2];
    std::size_t gl_nx = n_xyz_in/2;
    std::size_t gl_ny = n_xyz_in/2;
    std::size_t gl_nz = n_xyz_in/2;
    FPNumber gl_dx = (gl_x1 - gl_x0)/(FPNumber)(gl_nx);
    FPNumber gl_dy = (gl_y1 - gl_y0)/(FPNumber)(gl_ny);
    FPNumber gl_dz = (gl_z1 - gl_z0)/(FPNumber)(gl_nz);
    FPNumber gl_dt = 2.0*dt;

    FPNumber gu_x0 = box_in_r0[0];
    FPNumber gu_x1 = box_in_r1[0];
    FPNumber gu_y0 = box_in_r1[1];
    FPNumber gu_y1 = box_out_r1[1];
    FPNumber gu_z0 = box_out_r0[2];
    FPNumber gu_z1 = box_out_r1[2];
    std::size_t gu_nx = n_xyz_in/2;
    std::size_t gu_ny = n_xyz_in/2;
    std::size_t gu_nz = n_xyz_in/2*3;
    FPNumber gu_dx = (gu_x1 - gu_x0)/(FPNumber)(gu_nx);
    FPNumber gu_dy = (gu_y1 - gu_y0)/(FPNumber)(gu_ny);
    FPNumber gu_dz = (gu_z1 - gu_z0)/(FPNumber)(gu_nz);
    FPNumber gu_dt = 2.0*dt;

    std::size_t gr_indxSave = indxJ / 2;
    std::size_t gl_indxSave = gr_indxSave;
    std::size_t gu_indxSave = gr_indxSave;

    std::size_t gr_save_rate = 1;
    std::size_t gl_save_rate = gr_save_rate;
    std::size_t gu_save_rate = gr_save_rate;
    std::size_t save_rate = 2*gr_save_rate;

    std::size_t numOfTimeSamples = 150;
    int interpolationOrder = 1;
    std::size_t nt_ip0, nt_ip1;
    if(interpolationOrder == 0) {
        nt_ip0 = numOfTimeSamples;
        nt_ip1 = 0;
    } else if(interpolationOrder == 1) {
        nt_ip1 = numOfTimeSamples;
        nt_ip0 = 0;
    } else {
        assert(false);
    }


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
            {"\"_nx_m1_\"", boost::lexical_cast<std::string>(nx - 1)},
            {"\"_nx_m2_\"", boost::lexical_cast<std::string>(nx - 2)},
            {"\"_ny_p1_\"", boost::lexical_cast<std::string>(ny + 1)},
            {"\"_ny_m1_\"", boost::lexical_cast<std::string>(ny - 1)},
            {"\"_ny_m2_\"", boost::lexical_cast<std::string>(ny - 2)},
            {"\"_nz_p1_\"", boost::lexical_cast<std::string>(nz + 1)},
            {"\"_nz_m1_\"", boost::lexical_cast<std::string>(nz - 1)},
            {"\"_nz_m2_\"", boost::lexical_cast<std::string>(nz - 2)},
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
            {"\"_indxSave_\"", boost::lexical_cast<std::string>(indxJ)},
            {"\"_indxSave_p1_\"", boost::lexical_cast<std::string>(indxJ + 1)},
            {"\"_save_rate_\"", boost::lexical_cast<std::string>(save_rate)},
            {"\"_gr_x0_\"", boost::lexical_cast<std::string>(std::real(gr_x0))},
            {"\"_gr_x1_\"", boost::lexical_cast<std::string>(std::real(gr_x1))},
            {"\"_gr_y0_\"", boost::lexical_cast<std::string>(std::real(gr_y0))},
            {"\"_gr_y1_\"", boost::lexical_cast<std::string>(std::real(gr_y1))},
            {"\"_gr_z0_\"", boost::lexical_cast<std::string>(std::real(gr_z0))},
            {"\"_gr_z1_\"", boost::lexical_cast<std::string>(std::real(gr_z1))},
            {"\"_gr_nx_\"", boost::lexical_cast<std::string>(gr_nx)},
            {"\"_gr_ny_\"", boost::lexical_cast<std::string>(gr_ny)},
            {"\"_gr_nz_\"", boost::lexical_cast<std::string>(gr_nz)},
            {"\"_gr_nx_p1_\"", boost::lexical_cast<std::string>(gr_nx + 1)},
            {"\"_gr_ny_p1_\"", boost::lexical_cast<std::string>(gr_ny + 1)},
            {"\"_gr_nz_p1_\"", boost::lexical_cast<std::string>(gr_nz + 1)},
            {"\"_gr_nx_m1_\"", boost::lexical_cast<std::string>(gr_nx - 1)},
            {"\"_gr_ny_m1_\"", boost::lexical_cast<std::string>(gr_ny - 1)},
            {"\"_gr_nz_m1_\"", boost::lexical_cast<std::string>(gr_nz - 1)},
            {"\"_gr_dx_\"", boost::lexical_cast<std::string>(std::real(gr_dx))},
            {"\"_gr_dy_\"", boost::lexical_cast<std::string>(std::real(gr_dy))},
            {"\"_gr_dz_\"", boost::lexical_cast<std::string>(std::real(gr_dz))},
            {"\"_gr_dt_\"", boost::lexical_cast<std::string>(std::real(gr_dt))},
            {"\"_gr_dt_dx_\"", boost::lexical_cast<std::string>(std::real(gr_dt/gr_dx))},
            {"\"_gr_m_dt_dx_\"", boost::lexical_cast<std::string>(std::real(-gr_dt/gr_dx))},
            {"\"_gr_dt_dy_\"", boost::lexical_cast<std::string>(std::real(gr_dt/gr_dy))},
            {"\"_gr_m_dt_dy_\"", boost::lexical_cast<std::string>(std::real(-gr_dt/gr_dy))},
            {"\"_gr_dt_dz_\"", boost::lexical_cast<std::string>(std::real(gr_dt/gr_dz))},
            {"\"_gr_m_dt_dz_\"", boost::lexical_cast<std::string>(std::real(-gr_dt/gr_dz))},
            {"\"_gr_dt_dx_2_\"", boost::lexical_cast<std::string>(std::real(gr_dt/gr_dx/2.0))},
            {"\"_gr_m_dt_dx_2_\"", boost::lexical_cast<std::string>(std::real(-gr_dt/gr_dx/2.0))},
            {"\"_gr_dt_dy_2_\"", boost::lexical_cast<std::string>(std::real(gr_dt/gr_dy/2.0))},
            {"\"_gr_m_dt_dy_2_\"", boost::lexical_cast<std::string>(std::real(-gr_dt/gr_dy/2.0))},
            {"\"_gr_dt_dz_2_\"", boost::lexical_cast<std::string>(std::real(gr_dt/gr_dz/2.0))},
            {"\"_gr_m_dt_dz_2_\"", boost::lexical_cast<std::string>(std::real(-gr_dt/gr_dz/2.0))},
            {"\"_gr_indxSave_\"", boost::lexical_cast<std::string>(gr_indxSave)},
            {"\"_gr_indxSave_p1_\"", boost::lexical_cast<std::string>(gr_indxSave + 1)},
            {"\"_gr_save_rate_\"", boost::lexical_cast<std::string>(gr_save_rate)},
            {"\"_gl_x0_\"", boost::lexical_cast<std::string>(std::real(gl_x0))},
            {"\"_gl_x1_\"", boost::lexical_cast<std::string>(std::real(gl_x1))},
            {"\"_gl_y0_\"", boost::lexical_cast<std::string>(std::real(gl_y0))},
            {"\"_gl_y1_\"", boost::lexical_cast<std::string>(std::real(gl_y1))},
            {"\"_gl_z0_\"", boost::lexical_cast<std::string>(std::real(gl_z0))},
            {"\"_gl_z1_\"", boost::lexical_cast<std::string>(std::real(gl_z1))},
            {"\"_gl_nx_\"", boost::lexical_cast<std::string>(gl_nx)},
            {"\"_gl_ny_\"", boost::lexical_cast<std::string>(gl_ny)},
            {"\"_gl_nz_\"", boost::lexical_cast<std::string>(gl_nz)},
            {"\"_gl_nx_p1_\"", boost::lexical_cast<std::string>(gl_nx + 1)},
            {"\"_gl_ny_p1_\"", boost::lexical_cast<std::string>(gl_ny + 1)},
            {"\"_gl_nz_p1_\"", boost::lexical_cast<std::string>(gl_nz + 1)},
            {"\"_gl_nx_m1_\"", boost::lexical_cast<std::string>(gl_nx - 1)},
            {"\"_gl_ny_m1_\"", boost::lexical_cast<std::string>(gl_ny - 1)},
            {"\"_gl_nz_m1_\"", boost::lexical_cast<std::string>(gl_nz - 1)},
            {"\"_gl_nz_p_nz2_\"", boost::lexical_cast<std::string>(gl_nz + nz/2)},
            {"\"_gl_nz_p_nz2_m1_\"", boost::lexical_cast<std::string>(gl_nz + nz/2 - 1)},
            {"\"_gl_dx_\"", boost::lexical_cast<std::string>(std::real(gl_dx))},
            {"\"_gl_dy_\"", boost::lexical_cast<std::string>(std::real(gl_dy))},
            {"\"_gl_dz_\"", boost::lexical_cast<std::string>(std::real(gl_dz))},
            {"\"_gl_dt_\"", boost::lexical_cast<std::string>(std::real(gl_dt))},
            {"\"_gl_dt_dx_\"", boost::lexical_cast<std::string>(std::real(gl_dt/gl_dx))},
            {"\"_gl_m_dt_dx_\"", boost::lexical_cast<std::string>(std::real(-gl_dt/gl_dx))},
            {"\"_gl_dt_dy_\"", boost::lexical_cast<std::string>(std::real(gl_dt/gl_dy))},
            {"\"_gl_m_dt_dy_\"", boost::lexical_cast<std::string>(std::real(-gl_dt/gl_dy))},
            {"\"_gl_dt_dz_\"", boost::lexical_cast<std::string>(std::real(gl_dt/gl_dz))},
            {"\"_gl_m_dt_dz_\"", boost::lexical_cast<std::string>(std::real(-gl_dt/gl_dz))},
            {"\"_gl_dt_dx_2_\"", boost::lexical_cast<std::string>(std::real(gl_dt/gl_dx/2.0))},
            {"\"_gl_m_dt_dx_2_\"", boost::lexical_cast<std::string>(std::real(-gl_dt/gl_dx/2.0))},
            {"\"_gl_dt_dy_2_\"", boost::lexical_cast<std::string>(std::real(gl_dt/gl_dy/2.0))},
            {"\"_gl_m_dt_dy_2_\"", boost::lexical_cast<std::string>(std::real(-gl_dt/gl_dy/2.0))},
            {"\"_gl_dt_dz_2_\"", boost::lexical_cast<std::string>(std::real(gl_dt/gl_dz/2.0))},
            {"\"_gl_m_dt_dz_2_\"", boost::lexical_cast<std::string>(std::real(-gl_dt/gl_dz/2.0))},
            {"\"_gl_indxSave_\"", boost::lexical_cast<std::string>(gl_indxSave)},
            {"\"_gl_indxSave_p1_\"", boost::lexical_cast<std::string>(gl_indxSave + 1)},
            {"\"_gl_save_rate_\"", boost::lexical_cast<std::string>(gl_save_rate)},
            {"\"_gu_x0_\"", boost::lexical_cast<std::string>(std::real(gu_x0))},
            {"\"_gu_x1_\"", boost::lexical_cast<std::string>(std::real(gu_x1))},
            {"\"_gu_y0_\"", boost::lexical_cast<std::string>(std::real(gu_y0))},
            {"\"_gu_y1_\"", boost::lexical_cast<std::string>(std::real(gu_y1))},
            {"\"_gu_z0_\"", boost::lexical_cast<std::string>(std::real(gu_z0))},
            {"\"_gu_z1_\"", boost::lexical_cast<std::string>(std::real(gu_z1))},
            {"\"_gu_nx_\"", boost::lexical_cast<std::string>(gu_nx)},
            {"\"_gu_ny_\"", boost::lexical_cast<std::string>(gu_ny)},
            {"\"_gu_nz_\"", boost::lexical_cast<std::string>(gu_nz)},
            {"\"_gu_nx_p1_\"", boost::lexical_cast<std::string>(gu_nx + 1)},
            {"\"_gu_ny_p1_\"", boost::lexical_cast<std::string>(gu_ny + 1)},
            {"\"_gu_nz_p1_\"", boost::lexical_cast<std::string>(gu_nz + 1)},
            {"\"_gu_nx_m1_\"", boost::lexical_cast<std::string>(gu_nx - 1)},
            {"\"_gu_ny_m1_\"", boost::lexical_cast<std::string>(gu_ny - 1)},
            {"\"_gu_nz_m1_\"", boost::lexical_cast<std::string>(gu_nz - 1)},
            {"\"_gu_dx_\"", boost::lexical_cast<std::string>(std::real(gu_dx))},
            {"\"_gu_dy_\"", boost::lexical_cast<std::string>(std::real(gu_dy))},
            {"\"_gu_dz_\"", boost::lexical_cast<std::string>(std::real(gu_dz))},
            {"\"_gu_dt_\"", boost::lexical_cast<std::string>(std::real(gu_dt))},
            {"\"_gu_dt_dx_\"", boost::lexical_cast<std::string>(std::real(gu_dt/gu_dx))},
            {"\"_gu_m_dt_dx_\"", boost::lexical_cast<std::string>(std::real(-gu_dt/gu_dx))},
            {"\"_gu_dt_dy_\"", boost::lexical_cast<std::string>(std::real(gu_dt/gu_dy))},
            {"\"_gu_m_dt_dy_\"", boost::lexical_cast<std::string>(std::real(-gu_dt/gu_dy))},
            {"\"_gu_dt_dz_\"", boost::lexical_cast<std::string>(std::real(gu_dt/gu_dz))},
            {"\"_gu_m_dt_dz_\"", boost::lexical_cast<std::string>(std::real(-gu_dt/gu_dz))},
            {"\"_gu_dt_dx_2_\"", boost::lexical_cast<std::string>(std::real(gu_dt/gu_dx/2.0))},
            {"\"_gu_m_dt_dx_2_\"", boost::lexical_cast<std::string>(std::real(-gu_dt/gu_dx/2.0))},
            {"\"_gu_dt_dy_2_\"", boost::lexical_cast<std::string>(std::real(gu_dt/gu_dy/2.0))},
            {"\"_gu_m_dt_dy_2_\"", boost::lexical_cast<std::string>(std::real(-gu_dt/gu_dy/2.0))},
            {"\"_gu_dt_dz_2_\"", boost::lexical_cast<std::string>(std::real(gu_dt/gu_dz/2.0))},
            {"\"_gu_m_dt_dz_2_\"", boost::lexical_cast<std::string>(std::real(-gu_dt/gu_dz/2.0))},
            {"\"_gu_indxSave_\"", boost::lexical_cast<std::string>(gu_indxSave)},
            {"\"_gu_indxSave_p1_\"", boost::lexical_cast<std::string>(gu_indxSave + 1)},
            {"\"_gu_save_rate_\"", boost::lexical_cast<std::string>(gu_save_rate)},
            {"\"_nt_coarse_IP0_\"", boost::lexical_cast<std::string>(nt_ip0)},
            {"\"_nt_coarse_IP1_\"", boost::lexical_cast<std::string>(nt_ip1)},
            {"\"_mod_phase_\"", boost::lexical_cast<std::string>(M_PI/2.0)}
            };
    ParameterExtractor::ReplaceStringsInFile("instructions/MaxwellYee3D_Nonuniform_GridCollection.json",
                "instructions/processed/MaxwellYee3D_Nonuniform_GridCollection_processed.json", str_replacewith);

    ParamFileTranslator fileTranslator("instructions/processed/MaxwellYee3D_Nonuniform_GridCollection_processed.json");
    fileTranslator.Translate();

   std::string parametersFileName =
            std::string("data/2D/") + "params" + ".param";

    std::ofstream paramFileOut(parametersFileName.c_str(), std::ios::out | std::ios::binary);
    UtilityFunctions::WriteParamToFile<FPNumber>(paramFileOut, x0, "x0");
    UtilityFunctions::WriteParamToFile<FPNumber>(paramFileOut, x1, "x1");
    UtilityFunctions::WriteParamToFile<FPNumber>(paramFileOut, y0, "y0");
    UtilityFunctions::WriteParamToFile<FPNumber>(paramFileOut, y1, "y1");
    UtilityFunctions::WriteParamToFile<FPNumber>(paramFileOut, z0, "z0");
    UtilityFunctions::WriteParamToFile<FPNumber>(paramFileOut, z1, "z1");
    UtilityFunctions::WriteParamToFile<FPNumber>(paramFileOut, gr_x0, "gr_x0");
    UtilityFunctions::WriteParamToFile<FPNumber>(paramFileOut, gr_x1, "gr_x1");
    UtilityFunctions::WriteParamToFile<FPNumber>(paramFileOut, gr_y0, "gr_y0");
    UtilityFunctions::WriteParamToFile<FPNumber>(paramFileOut, gr_y1, "gr_y1");
    UtilityFunctions::WriteParamToFile<FPNumber>(paramFileOut, gr_z0, "gr_z0");
    UtilityFunctions::WriteParamToFile<FPNumber>(paramFileOut, gr_z1, "gr_z1");
    UtilityFunctions::WriteParamToFile<FPNumber>(paramFileOut, gl_x0, "gl_x0");
    UtilityFunctions::WriteParamToFile<FPNumber>(paramFileOut, gl_x1, "gl_x1");
    UtilityFunctions::WriteParamToFile<FPNumber>(paramFileOut, gl_y0, "gl_y0");
    UtilityFunctions::WriteParamToFile<FPNumber>(paramFileOut, gl_y1, "gl_y1");
    UtilityFunctions::WriteParamToFile<FPNumber>(paramFileOut, gl_z0, "gl_z0");
    UtilityFunctions::WriteParamToFile<FPNumber>(paramFileOut, gl_z1, "gl_z1");
    UtilityFunctions::WriteParamToFile<FPNumber>(paramFileOut, gu_x0, "gu_x0");
    UtilityFunctions::WriteParamToFile<FPNumber>(paramFileOut, gu_x1, "gu_x1");
    UtilityFunctions::WriteParamToFile<FPNumber>(paramFileOut, gu_y0, "gu_y0");
    UtilityFunctions::WriteParamToFile<FPNumber>(paramFileOut, gu_y1, "gu_y1");
    UtilityFunctions::WriteParamToFile<FPNumber>(paramFileOut, gu_z0, "gu_z0");
    UtilityFunctions::WriteParamToFile<FPNumber>(paramFileOut, gu_z1, "gu_z1");

}



