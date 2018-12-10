


void test_yeegrid_3d() {
    std::size_t nx = 100;
    std::size_t nz = 100;
    std::size_t ny = 100;
    std::size_t indxJx = nx/2;
    std::size_t indzJx = nz/2;
    std::size_t indyJx = ny/2;
    std::array<FPNumber, 3> r0{0.0, 0.0, 0.0};
    std::array<FPNumber, 3> r1{10.0, 10.0, 10.0};
    std::array<std::size_t, 3> nCells{nx, ny, nz};
    FPNumber dx = (r1[0] - r0[0])/(FPNumber)(nx);
    FPNumber dy = (r1[1] - r0[1])/(FPNumber)(ny);
    FPNumber dz = (r1[2] - r0[2])/(FPNumber)(nz);
    FPNumber dt = dz/std::sqrt((FPNumber)3.0)*(FPNumber)0.99;
    YeeGrid3D yee;
    yee.SetCornerCoordinates(r0, r1);
    yee.SetNumOfCells(nCells);
    yee.SetTimeResolution(dt);
    yee.AddEntireGridElement("E", ElementType::EdgeE);
    yee.AddEntireGridElement("H", ElementType::EdgeH);
    yee.AddPartialGridElement("J", ElementType::EdgeE, {indxJx, indyJx, indzJx}, {1, 0, 0});
    yee.AddGaussianGridArrayManipulator("JUpdater", "J", 0 /*x*/, 1.0 /*amp*/, 1.0 /*t_center*/, 0.2 /*t_decay*/,
            1.0 /*modulation frequency*/, M_PI/2.0 /*modulation phase*/, -0.5 /*time offset fraction*/);

    void* Ex_Hy_update_params = FDInstructionFactory::Get_AxEdgeE_plusEqual_b_Curl_CyEdgeH(
                yee,
                {0, 1, 1},        // indStart
                {nx, ny, nz},      // indEnd
                "E",              // name of A
                +dt,             // b value
                "H"                  // name of C
                );
    void* Ex_Hz_update_params = FDInstructionFactory::Get_AxEdgeE_plusEqual_b_Curl_CzEdgeH(
                yee,
                {0, 1, 1},        // indStart
                {nx, ny, nz},      // indEnd
                "E",              // name of A
                +dt,             // b value
                "H"                  // name of C
                );
    void* Ey_Hx_update_params = FDInstructionFactory::Get_AyEdgeE_plusEqual_b_Curl_CxEdgeH(
                yee,
                {1, 0, 1},        // indStart
                {nx, ny, nz},      // indEnd
                "E",              // name of A
                +dt,             // b value
                "H"                  // name of C
                );
    void* Ey_Hz_update_params = FDInstructionFactory::Get_AyEdgeE_plusEqual_b_Curl_CzEdgeH(
                yee,
                {1, 0, 1},        // indStart
                {nx, ny, nz},      // indEnd
                "E",              // name of A
                +dt,             // b value
                "H"                  // name of C
                );
    void* Ez_Hx_update_params = FDInstructionFactory::Get_AzEdgeE_plusEqual_b_Curl_CxEdgeH(
                yee,
                {1, 1, 0},        // indStart
                {nx, ny, nz},      // indEnd
                "E",              // name of A
                +dt,             // b value
                "H"                  // name of C
                );
    void* Ez_Hy_update_params = FDInstructionFactory::Get_AzEdgeE_plusEqual_b_Curl_CyEdgeH(
                yee,
                {1, 1, 0},        // indStart
                {nx, ny, nz},      // indEnd
                "E",              // name of A
                +dt,             // b value
                "H"                  // name of C
                );
    void* E_J_update_params = yee.ConstructParams_A_plusequal_sum_b_C(
        {indxJx, indyJx, indzJx},
        {indxJx + 1, indyJx + 1, indzJx + 1},
        "E",
        0,
        {-dt/(dx*dy*dz)},
        {"J"},
        {0},
        {{0, 0, 0}}
    );
    void* Hx_Ey_update_params = FDInstructionFactory::Get_AxEdgeH_plusEqual_b_Curl_CyEdgeE(
            yee,
            {0, 0, 0},
            {nx + 1, ny, nz},
            "H",
            -dt,
            "E");
    void* Hx_Ez_update_params = FDInstructionFactory::Get_AxEdgeH_plusEqual_b_Curl_CzEdgeE(
            yee,
            {0, 0, 0},
            {nx + 1, ny, nz},
            "H",
            -dt,
            "E");
    void* Hy_Ex_update_params = FDInstructionFactory::Get_AyEdgeH_plusEqual_b_Curl_CxEdgeE(
            yee,
            {0, 0, 0},
            {nx, ny + 1, nz},
            "H",
            -dt,
            "E");
    void* Hy_Ez_update_params = FDInstructionFactory::Get_AyEdgeH_plusEqual_b_Curl_CzEdgeE(
            yee,
            {0, 0, 0},
            {nx, ny + 1, nz},
            "H",
            -dt,
            "E");
    void* Hz_Ex_update_params = FDInstructionFactory::Get_AzEdgeH_plusEqual_b_Curl_CxEdgeE(
            yee,
            {0, 0, 0},
            {nx, ny, nz + 1},
            "H",
            -dt,
            "E");
    void* Hz_Ey_update_params = FDInstructionFactory::Get_AzEdgeH_plusEqual_b_Curl_CyEdgeE(
            yee,
            {0, 0, 0},
            {nx, ny, nz + 1},
            "H",
            -dt,
            "E");
    void* J_update_params = yee.ConstructParams_A_equal_func_r_t(
        "JUpdater"
    );
    yee.AddUpdateInstruction("J-update", FDInstructionCode::A_equal_func_r_t, J_update_params);
    yee.AddUpdateInstruction("Ex-Hy-update", FDInstructionCode::A_plusequal_sum_b_C, Ex_Hy_update_params);
    yee.AddUpdateInstruction("Ex-Hz-update", FDInstructionCode::A_plusequal_sum_b_C, Ex_Hz_update_params);
    yee.AddUpdateInstruction("Ey-Hx-update", FDInstructionCode::A_plusequal_sum_b_C, Ey_Hx_update_params);
    yee.AddUpdateInstruction("Ey-Hz-update", FDInstructionCode::A_plusequal_sum_b_C, Ey_Hz_update_params);
    yee.AddUpdateInstruction("Ez-Hx-update", FDInstructionCode::A_plusequal_sum_b_C, Ez_Hx_update_params);
    yee.AddUpdateInstruction("Ez-Hy-update", FDInstructionCode::A_plusequal_sum_b_C, Ez_Hy_update_params);
    yee.AddUpdateInstruction("E-J-update", FDInstructionCode::A_plusequal_sum_b_C, E_J_update_params);
    yee.AddUpdateInstruction("Hx-Ey-update", FDInstructionCode::A_plusequal_sum_b_C, Hx_Ey_update_params);
    yee.AddUpdateInstruction("Hx-Ez-update", FDInstructionCode::A_plusequal_sum_b_C, Hx_Ez_update_params);
    yee.AddUpdateInstruction("Hy-Ex-update", FDInstructionCode::A_plusequal_sum_b_C, Hy_Ex_update_params);
    yee.AddUpdateInstruction("Hy-Ez-update", FDInstructionCode::A_plusequal_sum_b_C, Hy_Ez_update_params);
    yee.AddUpdateInstruction("Hz-Ex-update", FDInstructionCode::A_plusequal_sum_b_C, Hz_Ex_update_params);
    yee.AddUpdateInstruction("Hz-Ey-update", FDInstructionCode::A_plusequal_sum_b_C, Hz_Ey_update_params);
    yee.SetIterativeSequence({"J-update", "Ex-Hy-update", "Ex-Hz-update",
                                          "Ey-Hx-update", "Ey-Hz-update",
                                          "Ez-Hx-update", "Ez-Hy-update",
                                          "E-J-update",
                                          "Hx-Ey-update", "Hx-Ez-update",
                                          "Hy-Ex-update", "Hy-Ez-update",
                                          "Hz-Ex-update", "Hz-Ey-update"});
    yee.AddGridElementView("E-x", "E", 0, {indxJx, 0, 0}, {indxJx + 1, ny + 1, nz + 1});
    yee.AddGridElementView("H-y", "H", 1, {indxJx + int(nx/4), 0, 0}, {indxJx + int(nx/4) + 1, ny + 1, nz});
    yee.DeleteOlderViewFiles();
    yee.SetDataStoreRate("E-x", 1);
    yee.SetDataStoreRate("H-y", 1);
    yee.ApplyIterativeInstructions(150);
}



