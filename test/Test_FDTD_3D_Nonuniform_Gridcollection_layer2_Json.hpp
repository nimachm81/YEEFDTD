


void test_run_fdtd_3d_nonuniform_gridcollection_layer2_from_json() {

    std::array<FPNumber, 3> box_lev0_r0 = {-0.5, -0.5, -0.5};
    std::array<FPNumber, 3> box_lev0_r1 = {+0.5, +0.5, +0.5};

    std::size_t n_xyz_lev0 = 100;

    FPNumber x0 = box_lev0_r0[0];
    FPNumber x1 = box_lev0_r1[0];
    FPNumber y0 = box_lev0_r0[1];
    FPNumber y1 = box_lev0_r1[1];
    FPNumber z0 = box_lev0_r0[2];
    FPNumber z1 = box_lev0_r1[2];
    std::size_t nx = n_xyz_lev0;
    std::size_t ny = n_xyz_lev0;
    std::size_t nz = n_xyz_lev0;
    FPNumber dx = (x1 - x0)/(FPNumber)(nx);
    FPNumber dy = (y1 - y0)/(FPNumber)(ny);
    FPNumber dz = (z1 - z0)/(FPNumber)(nz);
    FPNumber stabilityFactor = 0.99;
    FPNumber dt = (FPNumber)1.0/std::sqrt((FPNumber)1.0/(dx*dx) + (FPNumber)1.0/(dy*dy) + (FPNumber)1.0/(dz*dz))*stabilityFactor;

    FPNumber dx_lev1 = 2.0*dx;
    FPNumber dy_lev1 = 2.0*dy;
    FPNumber dz_lev1 = 2.0*dz;
    FPNumber dt_lev1 = 2.0*dt;
    std::size_t dn_lev1 = 20;
    std::array<FPNumber, 3> box_lev1_r0 = {box_lev0_r0[0] - dn_lev1*dx_lev1,
                                           box_lev0_r0[1] - dn_lev1*dy_lev1,
                                           box_lev0_r0[2] - dn_lev1*dz_lev1};
    std::array<FPNumber, 3> box_lev1_r1 = {box_lev0_r1[0] + dn_lev1*dx_lev1,
                                           box_lev0_r1[1] + dn_lev1*dy_lev1,
                                           box_lev0_r1[2] + dn_lev1*dz_lev1};

    FPNumber dx_lev2 = 2.0*dx_lev1;
    FPNumber dy_lev2 = 2.0*dy_lev1;
    FPNumber dz_lev2 = 2.0*dz_lev1;
    FPNumber dt_lev2 = 2.0*dt_lev1;
    std::size_t dn_lev2 = 20;
    std::array<FPNumber, 3> box_lev2_r0 = {box_lev1_r0[0] - dn_lev2*dx_lev2,
                                           box_lev1_r0[1] - dn_lev2*dy_lev2,
                                           box_lev1_r0[2] - dn_lev2*dz_lev2};
    std::array<FPNumber, 3> box_lev2_r1 = {box_lev1_r1[0] + dn_lev2*dx_lev2,
                                           box_lev1_r1[1] + dn_lev2*dy_lev2,
                                           box_lev1_r1[2] + dn_lev2*dz_lev2};


    FPNumber x_j = (x0 + x1)/(FPNumber)2.0;
    FPNumber y_j = (y0 + y1)/(FPNumber)2.0;
    FPNumber z_j = (z0 + z1)/(FPNumber)2.0;
    std::size_t indxJ = std::round(std::real((x_j - x0)/dx));
    std::size_t indyJ = std::round(std::real((y_j - y0)/dy));
    std::size_t indzJ = std::round(std::real((z_j - z0)/dz));

    std::size_t nx0_lev1 = nx/2;
    std::size_t ny0_lev1 = ny/2;
    std::size_t nz0_lev1 = nz/2;

    FPNumber gr_x0 = box_lev0_r0[0];
    FPNumber gr_x1 = box_lev0_r1[0];
    FPNumber gr_y0 = box_lev0_r0[1];
    FPNumber gr_y1 = box_lev0_r1[1];
    FPNumber gr_z0 = box_lev0_r1[2];
    FPNumber gr_z1 = box_lev1_r1[2];
    std::size_t gr_nx = nx0_lev1;
    std::size_t gr_ny = ny0_lev1;
    std::size_t gr_nz = dn_lev1;
    FPNumber gr_dx = dx_lev1;
    FPNumber gr_dy = dy_lev1;
    FPNumber gr_dz = dz_lev1;
    FPNumber gr_dt = dt_lev1;

    FPNumber gl_x0 = box_lev0_r0[0];
    FPNumber gl_x1 = box_lev0_r1[0];
    FPNumber gl_y0 = box_lev0_r0[1];
    FPNumber gl_y1 = box_lev0_r1[1];
    FPNumber gl_z0 = box_lev1_r0[2];
    FPNumber gl_z1 = box_lev0_r0[2];
    std::size_t gl_nx = nx0_lev1;
    std::size_t gl_ny = ny0_lev1;
    std::size_t gl_nz = dn_lev1;
    FPNumber gl_dx = dx_lev1;
    FPNumber gl_dy = dy_lev1;
    FPNumber gl_dz = dz_lev1;
    FPNumber gl_dt = dt_lev1;

    FPNumber gu_x0 = box_lev0_r0[0];
    FPNumber gu_x1 = box_lev0_r1[0];
    FPNumber gu_y0 = box_lev0_r1[1];
    FPNumber gu_y1 = box_lev1_r1[1];
    FPNumber gu_z0 = box_lev1_r0[2];
    FPNumber gu_z1 = box_lev1_r1[2];
    std::size_t gu_nx = nx0_lev1;
    std::size_t gu_ny = dn_lev1;
    std::size_t gu_nz = nz0_lev1 + 2*dn_lev1;
    FPNumber gu_dx = dx_lev1;
    FPNumber gu_dy = dy_lev1;
    FPNumber gu_dz = dz_lev1;
    FPNumber gu_dt = dt_lev1;

    FPNumber gd_x0 = box_lev0_r0[0];
    FPNumber gd_x1 = box_lev0_r1[0];
    FPNumber gd_y0 = box_lev1_r0[1];
    FPNumber gd_y1 = box_lev0_r0[1];
    FPNumber gd_z0 = box_lev1_r0[2];
    FPNumber gd_z1 = box_lev1_r1[2];
    std::size_t gd_nx = nx0_lev1;
    std::size_t gd_ny = dn_lev1;
    std::size_t gd_nz = nz0_lev1 + 2*dn_lev1;
    FPNumber gd_dx = dx_lev1;
    FPNumber gd_dy = dy_lev1;
    FPNumber gd_dz = dz_lev1;
    FPNumber gd_dt = dt_lev1;

    FPNumber gf_x0 = box_lev0_r1[0];
    FPNumber gf_x1 = box_lev1_r1[0];
    FPNumber gf_y0 = box_lev1_r0[1];
    FPNumber gf_y1 = box_lev1_r1[1];
    FPNumber gf_z0 = box_lev1_r0[2];
    FPNumber gf_z1 = box_lev1_r1[2];
    std::size_t gf_nx = dn_lev1;
    std::size_t gf_ny = ny0_lev1 + 2*dn_lev1;
    std::size_t gf_nz = nz0_lev1 + 2*dn_lev1;
    FPNumber gf_dx = dx_lev1;
    FPNumber gf_dy = dy_lev1;
    FPNumber gf_dz = dz_lev1;
    FPNumber gf_dt = dt_lev1;

    FPNumber gb_x0 = box_lev1_r0[0];
    FPNumber gb_x1 = box_lev0_r0[0];
    FPNumber gb_y0 = box_lev1_r0[1];
    FPNumber gb_y1 = box_lev1_r1[1];
    FPNumber gb_z0 = box_lev1_r0[2];
    FPNumber gb_z1 = box_lev1_r1[2];
    std::size_t gb_nx = dn_lev1;
    std::size_t gb_ny = ny0_lev1 + 2*dn_lev1;
    std::size_t gb_nz = nz0_lev1 + 2*dn_lev1;
    FPNumber gb_dx = dx_lev1;
    FPNumber gb_dy = dy_lev1;
    FPNumber gb_dz = dz_lev1;
    FPNumber gb_dt = dt_lev1;

    std::size_t nx0_lev2 = gr_nx/2 + gb_nx/2 + gf_nx/2;
    std::size_t ny0_lev2 = gb_ny/2;
    std::size_t nz0_lev2 = gb_nz/2;

    FPNumber grr_x0 = box_lev1_r0[0];
    FPNumber grr_x1 = box_lev1_r1[0];
    FPNumber grr_y0 = box_lev1_r0[1];
    FPNumber grr_y1 = box_lev1_r1[1];
    FPNumber grr_z0 = box_lev1_r1[2];
    FPNumber grr_z1 = box_lev2_r1[2];
    std::size_t grr_nx = nx0_lev2;
    std::size_t grr_ny = ny0_lev2;
    std::size_t grr_nz = dn_lev2;
    FPNumber grr_dx = dx_lev2;
    FPNumber grr_dy = dy_lev2;
    FPNumber grr_dz = dz_lev2;
    FPNumber grr_dt = dt_lev2;

    FPNumber gll_x0 = box_lev1_r0[0];
    FPNumber gll_x1 = box_lev1_r1[0];
    FPNumber gll_y0 = box_lev1_r0[1];
    FPNumber gll_y1 = box_lev1_r1[1];
    FPNumber gll_z0 = box_lev2_r0[2];
    FPNumber gll_z1 = box_lev1_r0[2];
    std::size_t gll_nx = nx0_lev2;
    std::size_t gll_ny = ny0_lev2;
    std::size_t gll_nz = dn_lev2;
    FPNumber gll_dx = dx_lev2;
    FPNumber gll_dy = dy_lev2;
    FPNumber gll_dz = dz_lev2;
    FPNumber gll_dt = dt_lev2;


    std::string j_polarization = "\"x\"";

    FPNumber x_save = 0.0;
    std::size_t gm_indxSave = std::round((x_save - x0 - dx/2) / dx);
    std::size_t gr_indxSave = std::round((x_save - gr_x0 - gr_dx/2) / gr_dx);
    std::size_t gl_indxSave = gr_indxSave;
    std::size_t gu_indxSave = gr_indxSave;
    std::size_t gd_indxSave = gr_indxSave;
    std::size_t grr_indxSave = std::round((x_save - grr_x0 - grr_dx/2) / grr_dx);
    std::size_t gll_indxSave = grr_indxSave;

    FPNumber y_save = 0.0;
    std::size_t gm_indySave = std::round((y_save - y0) / dy);
    std::size_t gr_indySave = std::round((y_save - gr_y0) / gr_dy);
    std::size_t gl_indySave = gr_indySave;
    std::size_t gf_indySave = std::round((y_save - gf_y0) / gf_dy);
    std::size_t gb_indySave = gf_indySave;
    std::size_t grr_indySave = std::round((y_save - grr_y0) / grr_dy);
    std::size_t gll_indySave = grr_indySave;


    std::size_t grr_save_rate = 1;
    std::size_t gll_save_rate = grr_save_rate;
    std::size_t gr_save_rate = 2*grr_save_rate;
    std::size_t gl_save_rate = gr_save_rate;
    std::size_t gu_save_rate = gr_save_rate;
    std::size_t gd_save_rate = gr_save_rate;
    std::size_t gf_save_rate = gr_save_rate;
    std::size_t gb_save_rate = gr_save_rate;
    std::size_t save_rate = 2*gr_save_rate;

    std::size_t numOfTimeSamples = 150;
    std::size_t nt_ip0, nt_ip1;
    nt_ip1 = numOfTimeSamples;

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
            {"\"_indxSave_\"", boost::lexical_cast<std::string>(gm_indxSave)},
            {"\"_indxSave_p1_\"", boost::lexical_cast<std::string>(gm_indxSave + 1)},
            {"\"_indySave_\"", boost::lexical_cast<std::string>(gm_indySave)},
            {"\"_indySave_p1_\"", boost::lexical_cast<std::string>(gm_indySave + 1)},
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
            {"\"_gr_nz_m2_\"", boost::lexical_cast<std::string>(gr_nz - 2)},
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
            {"\"_gr_indySave_\"", boost::lexical_cast<std::string>(gr_indySave)},
            {"\"_gr_indySave_p1_\"", boost::lexical_cast<std::string>(gr_indySave + 1)},
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
            {"\"_gl_indySave_\"", boost::lexical_cast<std::string>(gl_indySave)},
            {"\"_gl_indySave_p1_\"", boost::lexical_cast<std::string>(gl_indySave + 1)},
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
            {"\"_gu_ny_m2_\"", boost::lexical_cast<std::string>(gu_ny - 2)},
            {"\"_gu_nz_m2_\"", boost::lexical_cast<std::string>(gu_nz - 2)},
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
            {"\"_gd_x0_\"", boost::lexical_cast<std::string>(std::real(gd_x0))},
            {"\"_gd_x1_\"", boost::lexical_cast<std::string>(std::real(gd_x1))},
            {"\"_gd_y0_\"", boost::lexical_cast<std::string>(std::real(gd_y0))},
            {"\"_gd_y1_\"", boost::lexical_cast<std::string>(std::real(gd_y1))},
            {"\"_gd_z0_\"", boost::lexical_cast<std::string>(std::real(gd_z0))},
            {"\"_gd_z1_\"", boost::lexical_cast<std::string>(std::real(gd_z1))},
            {"\"_gd_nx_\"", boost::lexical_cast<std::string>(gd_nx)},
            {"\"_gd_ny_\"", boost::lexical_cast<std::string>(gd_ny)},
            {"\"_gd_nz_\"", boost::lexical_cast<std::string>(gd_nz)},
            {"\"_gd_nx_p1_\"", boost::lexical_cast<std::string>(gd_nx + 1)},
            {"\"_gd_ny_p1_\"", boost::lexical_cast<std::string>(gd_ny + 1)},
            {"\"_gd_nz_p1_\"", boost::lexical_cast<std::string>(gd_nz + 1)},
            {"\"_gd_nx_m1_\"", boost::lexical_cast<std::string>(gd_nx - 1)},
            {"\"_gd_ny_m1_\"", boost::lexical_cast<std::string>(gd_ny - 1)},
            {"\"_gd_nz_m1_\"", boost::lexical_cast<std::string>(gd_nz - 1)},
            {"\"_gd_nz_m2_\"", boost::lexical_cast<std::string>(gd_nz - 2)},
            {"\"_gd_ny_p_ny2_\"", boost::lexical_cast<std::string>(gd_ny + ny/2)},
            {"\"_gd_ny_p_ny2_m1_\"", boost::lexical_cast<std::string>(gd_ny + ny/2 - 1)},
            {"\"_gd_ny2_\"", boost::lexical_cast<std::string>(gd_ny/2)},
            {"\"_gd_ny2_p1_\"", boost::lexical_cast<std::string>(gd_ny/2 + 1)},
            {"\"_gd_ny2_m1_\"", boost::lexical_cast<std::string>(gd_ny/2 - 1)},
            {"\"_gd_ny2_p_gr_ny2_\"", boost::lexical_cast<std::string>(gd_ny/2 + gr_ny/2)},
            {"\"_gd_ny2_p_gr_ny2_p1_\"", boost::lexical_cast<std::string>(gd_ny/2 + gr_ny/2 + 1)},
            {"\"_gd_ny2_p_gr_ny2_m1_\"", boost::lexical_cast<std::string>(gd_ny/2 + gr_ny/2 - 1)},
            {"\"_gd_dx_\"", boost::lexical_cast<std::string>(std::real(gd_dx))},
            {"\"_gd_dy_\"", boost::lexical_cast<std::string>(std::real(gd_dy))},
            {"\"_gd_dz_\"", boost::lexical_cast<std::string>(std::real(gd_dz))},
            {"\"_gd_dt_\"", boost::lexical_cast<std::string>(std::real(gd_dt))},
            {"\"_gd_dt_dx_\"", boost::lexical_cast<std::string>(std::real(gd_dt/gd_dx))},
            {"\"_gd_m_dt_dx_\"", boost::lexical_cast<std::string>(std::real(-gd_dt/gd_dx))},
            {"\"_gd_dt_dy_\"", boost::lexical_cast<std::string>(std::real(gd_dt/gd_dy))},
            {"\"_gd_m_dt_dy_\"", boost::lexical_cast<std::string>(std::real(-gd_dt/gd_dy))},
            {"\"_gd_dt_dz_\"", boost::lexical_cast<std::string>(std::real(gd_dt/gd_dz))},
            {"\"_gd_m_dt_dz_\"", boost::lexical_cast<std::string>(std::real(-gd_dt/gd_dz))},
            {"\"_gd_dt_dx_2_\"", boost::lexical_cast<std::string>(std::real(gd_dt/gd_dx/2.0))},
            {"\"_gd_m_dt_dx_2_\"", boost::lexical_cast<std::string>(std::real(-gd_dt/gd_dx/2.0))},
            {"\"_gd_dt_dy_2_\"", boost::lexical_cast<std::string>(std::real(gd_dt/gd_dy/2.0))},
            {"\"_gd_m_dt_dy_2_\"", boost::lexical_cast<std::string>(std::real(-gd_dt/gd_dy/2.0))},
            {"\"_gd_dt_dz_2_\"", boost::lexical_cast<std::string>(std::real(gd_dt/gd_dz/2.0))},
            {"\"_gd_m_dt_dz_2_\"", boost::lexical_cast<std::string>(std::real(-gd_dt/gd_dz/2.0))},
            {"\"_gd_indxSave_\"", boost::lexical_cast<std::string>(gd_indxSave)},
            {"\"_gd_indxSave_p1_\"", boost::lexical_cast<std::string>(gd_indxSave + 1)},
            {"\"_gd_save_rate_\"", boost::lexical_cast<std::string>(gd_save_rate)},
            {"\"_gf_x0_\"", boost::lexical_cast<std::string>(std::real(gf_x0))},
            {"\"_gf_x1_\"", boost::lexical_cast<std::string>(std::real(gf_x1))},
            {"\"_gf_y0_\"", boost::lexical_cast<std::string>(std::real(gf_y0))},
            {"\"_gf_y1_\"", boost::lexical_cast<std::string>(std::real(gf_y1))},
            {"\"_gf_z0_\"", boost::lexical_cast<std::string>(std::real(gf_z0))},
            {"\"_gf_z1_\"", boost::lexical_cast<std::string>(std::real(gf_z1))},
            {"\"_gf_nx_\"", boost::lexical_cast<std::string>(gf_nx)},
            {"\"_gf_ny_\"", boost::lexical_cast<std::string>(gf_ny)},
            {"\"_gf_nz_\"", boost::lexical_cast<std::string>(gf_nz)},
            {"\"_gf_nx_p1_\"", boost::lexical_cast<std::string>(gf_nx + 1)},
            {"\"_gf_ny_p1_\"", boost::lexical_cast<std::string>(gf_ny + 1)},
            {"\"_gf_nz_p1_\"", boost::lexical_cast<std::string>(gf_nz + 1)},
            {"\"_gf_nx_m1_\"", boost::lexical_cast<std::string>(gf_nx - 1)},
            {"\"_gf_ny_m1_\"", boost::lexical_cast<std::string>(gf_ny - 1)},
            {"\"_gf_nz_m1_\"", boost::lexical_cast<std::string>(gf_nz - 1)},
            {"\"_gf_nx_m2_\"", boost::lexical_cast<std::string>(gf_nx - 2)},
            {"\"_gf_ny_m2_\"", boost::lexical_cast<std::string>(gf_ny - 2)},
            {"\"_gf_nz_m2_\"", boost::lexical_cast<std::string>(gf_nz - 2)},
            {"\"_gf_dx_\"", boost::lexical_cast<std::string>(std::real(gf_dx))},
            {"\"_gf_dy_\"", boost::lexical_cast<std::string>(std::real(gf_dy))},
            {"\"_gf_dz_\"", boost::lexical_cast<std::string>(std::real(gf_dz))},
            {"\"_gf_dt_\"", boost::lexical_cast<std::string>(std::real(gf_dt))},
            {"\"_gf_dt_dx_\"", boost::lexical_cast<std::string>(std::real(gf_dt/gf_dx))},
            {"\"_gf_m_dt_dx_\"", boost::lexical_cast<std::string>(std::real(-gf_dt/gf_dx))},
            {"\"_gf_dt_dy_\"", boost::lexical_cast<std::string>(std::real(gf_dt/gf_dy))},
            {"\"_gf_m_dt_dy_\"", boost::lexical_cast<std::string>(std::real(-gf_dt/gf_dy))},
            {"\"_gf_dt_dz_\"", boost::lexical_cast<std::string>(std::real(gf_dt/gf_dz))},
            {"\"_gf_m_dt_dz_\"", boost::lexical_cast<std::string>(std::real(-gf_dt/gf_dz))},
            {"\"_gf_dt_dx_2_\"", boost::lexical_cast<std::string>(std::real(gf_dt/gf_dx/2.0))},
            {"\"_gf_m_dt_dx_2_\"", boost::lexical_cast<std::string>(std::real(-gf_dt/gf_dx/2.0))},
            {"\"_gf_dt_dy_2_\"", boost::lexical_cast<std::string>(std::real(gf_dt/gf_dy/2.0))},
            {"\"_gf_m_dt_dy_2_\"", boost::lexical_cast<std::string>(std::real(-gf_dt/gf_dy/2.0))},
            {"\"_gf_dt_dz_2_\"", boost::lexical_cast<std::string>(std::real(gf_dt/gf_dz/2.0))},
            {"\"_gf_m_dt_dz_2_\"", boost::lexical_cast<std::string>(std::real(-gf_dt/gf_dz/2.0))},
            {"\"_gf_indySave_\"", boost::lexical_cast<std::string>(gf_indySave)},
            {"\"_gf_indySave_p1_\"", boost::lexical_cast<std::string>(gf_indySave + 1)},
            {"\"_gf_save_rate_\"", boost::lexical_cast<std::string>(gf_save_rate)},
            {"\"_gb_x0_\"", boost::lexical_cast<std::string>(std::real(gb_x0))},
            {"\"_gb_x1_\"", boost::lexical_cast<std::string>(std::real(gb_x1))},
            {"\"_gb_y0_\"", boost::lexical_cast<std::string>(std::real(gb_y0))},
            {"\"_gb_y1_\"", boost::lexical_cast<std::string>(std::real(gb_y1))},
            {"\"_gb_z0_\"", boost::lexical_cast<std::string>(std::real(gb_z0))},
            {"\"_gb_z1_\"", boost::lexical_cast<std::string>(std::real(gb_z1))},
            {"\"_gb_nx_\"", boost::lexical_cast<std::string>(gb_nx)},
            {"\"_gb_ny_\"", boost::lexical_cast<std::string>(gb_ny)},
            {"\"_gb_nz_\"", boost::lexical_cast<std::string>(gb_nz)},
            {"\"_gb_nx_p1_\"", boost::lexical_cast<std::string>(gb_nx + 1)},
            {"\"_gb_ny_p1_\"", boost::lexical_cast<std::string>(gb_ny + 1)},
            {"\"_gb_nz_p1_\"", boost::lexical_cast<std::string>(gb_nz + 1)},
            {"\"_gb_nx_m1_\"", boost::lexical_cast<std::string>(gb_nx - 1)},
            {"\"_gb_ny_m1_\"", boost::lexical_cast<std::string>(gb_ny - 1)},
            {"\"_gb_nz_m1_\"", boost::lexical_cast<std::string>(gb_nz - 1)},
            {"\"_gb_ny_m2_\"", boost::lexical_cast<std::string>(gb_ny - 2)},
            {"\"_gb_nz_m2_\"", boost::lexical_cast<std::string>(gb_nz - 2)},
            {"\"_gb_nx2_\"", boost::lexical_cast<std::string>(gb_nx/2)},
            {"\"_gb_nx2_p1_\"", boost::lexical_cast<std::string>(gb_nx/2 + 1)},
            {"\"_gb_nx2_m1_\"", boost::lexical_cast<std::string>(gb_nx/2 - 1)},
            {"\"_gb_nx2_p_gr_nx2_\"", boost::lexical_cast<std::string>(gb_nx/2 + gr_nx/2)},
            {"\"_gb_nx2_p_gr_nx2_p1_\"", boost::lexical_cast<std::string>(gb_nx/2 + gr_nx/2 + 1)},
            {"\"_gb_nx2_p_gr_nx2_m1_\"", boost::lexical_cast<std::string>(gb_nx/2 + gr_nx/2 - 1)},
            {"\"_gb_dx_\"", boost::lexical_cast<std::string>(std::real(gb_dx))},
            {"\"_gb_dy_\"", boost::lexical_cast<std::string>(std::real(gb_dy))},
            {"\"_gb_dz_\"", boost::lexical_cast<std::string>(std::real(gb_dz))},
            {"\"_gb_dt_\"", boost::lexical_cast<std::string>(std::real(gb_dt))},
            {"\"_gb_dt_dx_\"", boost::lexical_cast<std::string>(std::real(gb_dt/gb_dx))},
            {"\"_gb_m_dt_dx_\"", boost::lexical_cast<std::string>(std::real(-gb_dt/gb_dx))},
            {"\"_gb_dt_dy_\"", boost::lexical_cast<std::string>(std::real(gb_dt/gb_dy))},
            {"\"_gb_m_dt_dy_\"", boost::lexical_cast<std::string>(std::real(-gb_dt/gb_dy))},
            {"\"_gb_dt_dz_\"", boost::lexical_cast<std::string>(std::real(gb_dt/gb_dz))},
            {"\"_gb_m_dt_dz_\"", boost::lexical_cast<std::string>(std::real(-gb_dt/gb_dz))},
            {"\"_gb_dt_dx_2_\"", boost::lexical_cast<std::string>(std::real(gb_dt/gb_dx/2.0))},
            {"\"_gb_m_dt_dx_2_\"", boost::lexical_cast<std::string>(std::real(-gb_dt/gb_dx/2.0))},
            {"\"_gb_dt_dy_2_\"", boost::lexical_cast<std::string>(std::real(gb_dt/gb_dy/2.0))},
            {"\"_gb_m_dt_dy_2_\"", boost::lexical_cast<std::string>(std::real(-gb_dt/gb_dy/2.0))},
            {"\"_gb_dt_dz_2_\"", boost::lexical_cast<std::string>(std::real(gb_dt/gb_dz/2.0))},
            {"\"_gb_m_dt_dz_2_\"", boost::lexical_cast<std::string>(std::real(-gb_dt/gb_dz/2.0))},
            {"\"_gb_indySave_\"", boost::lexical_cast<std::string>(gb_indySave)},
            {"\"_gb_indySave_p1_\"", boost::lexical_cast<std::string>(gb_indySave + 1)},
            {"\"_gb_save_rate_\"", boost::lexical_cast<std::string>(gb_save_rate)},
            {"\"_grr_x0_\"", boost::lexical_cast<std::string>(std::real(grr_x0))},
            {"\"_grr_x1_\"", boost::lexical_cast<std::string>(std::real(grr_x1))},
            {"\"_grr_y0_\"", boost::lexical_cast<std::string>(std::real(grr_y0))},
            {"\"_grr_y1_\"", boost::lexical_cast<std::string>(std::real(grr_y1))},
            {"\"_grr_z0_\"", boost::lexical_cast<std::string>(std::real(grr_z0))},
            {"\"_grr_z1_\"", boost::lexical_cast<std::string>(std::real(grr_z1))},
            {"\"_grr_nx_\"", boost::lexical_cast<std::string>(grr_nx)},
            {"\"_grr_ny_\"", boost::lexical_cast<std::string>(grr_ny)},
            {"\"_grr_nz_\"", boost::lexical_cast<std::string>(grr_nz)},
            {"\"_grr_nx_p1_\"", boost::lexical_cast<std::string>(grr_nx + 1)},
            {"\"_grr_ny_p1_\"", boost::lexical_cast<std::string>(grr_ny + 1)},
            {"\"_grr_nz_p1_\"", boost::lexical_cast<std::string>(grr_nz + 1)},
            {"\"_grr_nx_m1_\"", boost::lexical_cast<std::string>(grr_nx - 1)},
            {"\"_grr_ny_m1_\"", boost::lexical_cast<std::string>(grr_ny - 1)},
            {"\"_grr_nz_m1_\"", boost::lexical_cast<std::string>(grr_nz - 1)},
            {"\"_grr_dx_\"", boost::lexical_cast<std::string>(std::real(grr_dx))},
            {"\"_grr_dy_\"", boost::lexical_cast<std::string>(std::real(grr_dy))},
            {"\"_grr_dz_\"", boost::lexical_cast<std::string>(std::real(grr_dz))},
            {"\"_grr_dt_\"", boost::lexical_cast<std::string>(std::real(grr_dt))},
            {"\"_grr_dt_dx_\"", boost::lexical_cast<std::string>(std::real(grr_dt/grr_dx))},
            {"\"_grr_m_dt_dx_\"", boost::lexical_cast<std::string>(std::real(-grr_dt/grr_dx))},
            {"\"_grr_dt_dy_\"", boost::lexical_cast<std::string>(std::real(grr_dt/grr_dy))},
            {"\"_grr_m_dt_dy_\"", boost::lexical_cast<std::string>(std::real(-grr_dt/grr_dy))},
            {"\"_grr_dt_dz_\"", boost::lexical_cast<std::string>(std::real(grr_dt/grr_dz))},
            {"\"_grr_m_dt_dz_\"", boost::lexical_cast<std::string>(std::real(-grr_dt/grr_dz))},
            {"\"_grr_dt_dx_2_\"", boost::lexical_cast<std::string>(std::real(grr_dt/grr_dx/2.0))},
            {"\"_grr_m_dt_dx_2_\"", boost::lexical_cast<std::string>(std::real(-grr_dt/grr_dx/2.0))},
            {"\"_grr_dt_dy_2_\"", boost::lexical_cast<std::string>(std::real(grr_dt/grr_dy/2.0))},
            {"\"_grr_m_dt_dy_2_\"", boost::lexical_cast<std::string>(std::real(-grr_dt/grr_dy/2.0))},
            {"\"_grr_dt_dz_2_\"", boost::lexical_cast<std::string>(std::real(grr_dt/grr_dz/2.0))},
            {"\"_grr_m_dt_dz_2_\"", boost::lexical_cast<std::string>(std::real(-grr_dt/grr_dz/2.0))},
            {"\"_grr_indxSave_\"", boost::lexical_cast<std::string>(grr_indxSave)},
            {"\"_grr_indxSave_p1_\"", boost::lexical_cast<std::string>(grr_indxSave + 1)},
            {"\"_grr_indySave_\"", boost::lexical_cast<std::string>(grr_indySave)},
            {"\"_grr_indySave_p1_\"", boost::lexical_cast<std::string>(grr_indySave + 1)},
            {"\"_grr_save_rate_\"", boost::lexical_cast<std::string>(grr_save_rate)},
            {"\"_gll_x0_\"", boost::lexical_cast<std::string>(std::real(gll_x0))},
            {"\"_gll_x1_\"", boost::lexical_cast<std::string>(std::real(gll_x1))},
            {"\"_gll_y0_\"", boost::lexical_cast<std::string>(std::real(gll_y0))},
            {"\"_gll_y1_\"", boost::lexical_cast<std::string>(std::real(gll_y1))},
            {"\"_gll_z0_\"", boost::lexical_cast<std::string>(std::real(gll_z0))},
            {"\"_gll_z1_\"", boost::lexical_cast<std::string>(std::real(gll_z1))},
            {"\"_gll_nx_\"", boost::lexical_cast<std::string>(gll_nx)},
            {"\"_gll_ny_\"", boost::lexical_cast<std::string>(gll_ny)},
            {"\"_gll_nz_\"", boost::lexical_cast<std::string>(gll_nz)},
            {"\"_gll_nx_p1_\"", boost::lexical_cast<std::string>(gll_nx + 1)},
            {"\"_gll_ny_p1_\"", boost::lexical_cast<std::string>(gll_ny + 1)},
            {"\"_gll_nz_p1_\"", boost::lexical_cast<std::string>(gll_nz + 1)},
            {"\"_gll_nx_m1_\"", boost::lexical_cast<std::string>(gll_nx - 1)},
            {"\"_gll_ny_m1_\"", boost::lexical_cast<std::string>(gll_ny - 1)},
            {"\"_gll_nz_m1_\"", boost::lexical_cast<std::string>(gll_nz - 1)},
            {"\"_gll_dx_\"", boost::lexical_cast<std::string>(std::real(gll_dx))},
            {"\"_gll_dy_\"", boost::lexical_cast<std::string>(std::real(gll_dy))},
            {"\"_gll_dz_\"", boost::lexical_cast<std::string>(std::real(gll_dz))},
            {"\"_gll_dt_\"", boost::lexical_cast<std::string>(std::real(gll_dt))},
            {"\"_gll_dt_dx_\"", boost::lexical_cast<std::string>(std::real(gll_dt/gll_dx))},
            {"\"_gll_m_dt_dx_\"", boost::lexical_cast<std::string>(std::real(-gll_dt/gll_dx))},
            {"\"_gll_dt_dy_\"", boost::lexical_cast<std::string>(std::real(gll_dt/gll_dy))},
            {"\"_gll_m_dt_dy_\"", boost::lexical_cast<std::string>(std::real(-gll_dt/gll_dy))},
            {"\"_gll_dt_dz_\"", boost::lexical_cast<std::string>(std::real(gll_dt/gll_dz))},
            {"\"_gll_m_dt_dz_\"", boost::lexical_cast<std::string>(std::real(-gll_dt/gll_dz))},
            {"\"_gll_dt_dx_2_\"", boost::lexical_cast<std::string>(std::real(gll_dt/gll_dx/2.0))},
            {"\"_gll_m_dt_dx_2_\"", boost::lexical_cast<std::string>(std::real(-gll_dt/gll_dx/2.0))},
            {"\"_gll_dt_dy_2_\"", boost::lexical_cast<std::string>(std::real(gll_dt/gll_dy/2.0))},
            {"\"_gll_m_dt_dy_2_\"", boost::lexical_cast<std::string>(std::real(-gll_dt/gll_dy/2.0))},
            {"\"_gll_dt_dz_2_\"", boost::lexical_cast<std::string>(std::real(gll_dt/gll_dz/2.0))},
            {"\"_gll_m_dt_dz_2_\"", boost::lexical_cast<std::string>(std::real(-gll_dt/gll_dz/2.0))},
            {"\"_gll_indxSave_\"", boost::lexical_cast<std::string>(gll_indxSave)},
            {"\"_gll_indxSave_p1_\"", boost::lexical_cast<std::string>(gll_indxSave + 1)},
            {"\"_gll_indySave_\"", boost::lexical_cast<std::string>(gll_indySave)},
            {"\"_gll_indySave_p1_\"", boost::lexical_cast<std::string>(gll_indySave + 1)},
            {"\"_gll_save_rate_\"", boost::lexical_cast<std::string>(gll_save_rate)},
            {"\"_nt_coarse_IP1_\"", boost::lexical_cast<std::string>(nt_ip1)},
            {"\"_mod_phase_\"", boost::lexical_cast<std::string>(M_PI/2.0)},
            {"\"_j_polarization_\"", j_polarization}
            };
    ParameterExtractor::ReplaceStringsInFile("instructions/MaxwellYee3D_Nonuniform_layer2.json",
                "instructions/processed/MaxwellYee3D_Nonuniform_layer2_processed.json", str_replacewith);

    ParamFileTranslator fileTranslator("instructions/processed/MaxwellYee3D_Nonuniform_layer2_processed.json");
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
    UtilityFunctions::WriteParamToFile<FPNumber>(paramFileOut, gd_x0, "gd_x0");
    UtilityFunctions::WriteParamToFile<FPNumber>(paramFileOut, gd_x1, "gd_x1");
    UtilityFunctions::WriteParamToFile<FPNumber>(paramFileOut, gd_y0, "gd_y0");
    UtilityFunctions::WriteParamToFile<FPNumber>(paramFileOut, gd_y1, "gd_y1");
    UtilityFunctions::WriteParamToFile<FPNumber>(paramFileOut, gd_z0, "gd_z0");
    UtilityFunctions::WriteParamToFile<FPNumber>(paramFileOut, gd_z1, "gd_z1");
    UtilityFunctions::WriteParamToFile<FPNumber>(paramFileOut, gf_x0, "gf_x0");
    UtilityFunctions::WriteParamToFile<FPNumber>(paramFileOut, gf_x1, "gf_x1");
    UtilityFunctions::WriteParamToFile<FPNumber>(paramFileOut, gf_y0, "gf_y0");
    UtilityFunctions::WriteParamToFile<FPNumber>(paramFileOut, gf_y1, "gf_y1");
    UtilityFunctions::WriteParamToFile<FPNumber>(paramFileOut, gf_z0, "gf_z0");
    UtilityFunctions::WriteParamToFile<FPNumber>(paramFileOut, gf_z1, "gf_z1");
    UtilityFunctions::WriteParamToFile<FPNumber>(paramFileOut, gb_x0, "gb_x0");
    UtilityFunctions::WriteParamToFile<FPNumber>(paramFileOut, gb_x1, "gb_x1");
    UtilityFunctions::WriteParamToFile<FPNumber>(paramFileOut, gb_y0, "gb_y0");
    UtilityFunctions::WriteParamToFile<FPNumber>(paramFileOut, gb_y1, "gb_y1");
    UtilityFunctions::WriteParamToFile<FPNumber>(paramFileOut, gb_z0, "gb_z0");
    UtilityFunctions::WriteParamToFile<FPNumber>(paramFileOut, gb_z1, "gb_z1");
    UtilityFunctions::WriteParamToFile<FPNumber>(paramFileOut, grr_x0, "grr_x0");
    UtilityFunctions::WriteParamToFile<FPNumber>(paramFileOut, grr_x1, "grr_x1");
    UtilityFunctions::WriteParamToFile<FPNumber>(paramFileOut, grr_y0, "grr_y0");
    UtilityFunctions::WriteParamToFile<FPNumber>(paramFileOut, grr_y1, "grr_y1");
    UtilityFunctions::WriteParamToFile<FPNumber>(paramFileOut, grr_z0, "grr_z0");
    UtilityFunctions::WriteParamToFile<FPNumber>(paramFileOut, grr_z1, "grr_z1");
    UtilityFunctions::WriteParamToFile<FPNumber>(paramFileOut, gll_x0, "gll_x0");
    UtilityFunctions::WriteParamToFile<FPNumber>(paramFileOut, gll_x1, "gll_x1");
    UtilityFunctions::WriteParamToFile<FPNumber>(paramFileOut, gll_y0, "gll_y0");
    UtilityFunctions::WriteParamToFile<FPNumber>(paramFileOut, gll_y1, "gll_y1");
    UtilityFunctions::WriteParamToFile<FPNumber>(paramFileOut, gll_z0, "gll_z0");
    UtilityFunctions::WriteParamToFile<FPNumber>(paramFileOut, gll_z1, "gll_z1");

}

