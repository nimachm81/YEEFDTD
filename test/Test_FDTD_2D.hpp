


void test_yeegrid_2d() {
    std::size_t nz = 200;
    std::size_t ny = 200;
    std::size_t indzJx = nz/2;
    std::size_t indyJx = ny/2;
    std::array<FPNumber, 3> r0{0.0, 0.0, 0.0};
    std::array<FPNumber, 3> r1{0.1, 10.0, 10.0};
    std::array<std::size_t, 3> nCells{1, ny, nz};
    FPNumber dz = (r1[2] - r0[2])/(FPNumber)(nz);
    FPNumber dy = (r1[1] - r0[1])/(FPNumber)(ny);
    FPNumber dt = dz/std::sqrt((FPNumber)2.0)*(FPNumber)0.99;
    YeeGrid3D yee;
    yee.SetCornerCoordinates(r0, r1);
    yee.SetNumOfCells(nCells);
    yee.SetTimeResolution(dt);
    yee.AddEntireGridElement("E", ElementType::EdgeE);
    yee.AddEntireGridElement("H", ElementType::EdgeH);
    yee.AddPartialGridElement("J", ElementType::EdgeE, {0, indyJx, indzJx}, {1, 0, 0});
    yee.AddGaussianGridArrayManipulator("JUpdater", "J", 0 /*x*/, 1.0 /*amp*/, 1.0 /*t_center*/, 0.2 /*t_decay*/,
            0.0 /*modulation frequency*/, 0.0 /*modulation phase*/, -0.5 /*time offset fraction*/);
    void* Ex_Hy_update_params = FDInstructionFactory::Get_AxEdgeE_plusEqual_b_Curl_CyEdgeH(  // A_x += b*curl C
                yee,
                {0, 1, 1},        // indStart
                {1, ny, nz},      // indEnd
                "E",              // name of A
                +dt,             // b value
                "H"                  // name of C
                );
    void* Ex_Hz_update_params = FDInstructionFactory::Get_AxEdgeE_plusEqual_b_Curl_CzEdgeH(  // A_x += b*curl C
                yee,
                {0, 1, 1},        // indStart
                {1, ny, nz},      // indEnd
                "E",              // name of A
                +dt,             // b value
                "H"                  // name of C
                );
    void* E_J_update_params = yee.ConstructParams_A_plusequal_sum_b_C(
        {0, indyJx, indzJx},
        {1, indyJx + 1, indzJx + 1},
        "E",
        0,
        {-(FPNumber)1.0*dt/(dy*dz)},
        {"J"},
        {0},
        {{0, 0, 0}}
    );
    void* Hy_update_params = FDInstructionFactory::Get_AyEdgeH_plusEqual_b_Curl_CxEdgeE(
            yee,
            {0, 0, 0},
            {1, ny + 1, nz},
            "H",
            -(FPNumber)1.0*dt,
            "E");
    void* Hz_update_params = FDInstructionFactory::Get_AzEdgeH_plusEqual_b_Curl_CxEdgeE(
            yee,
            {0, 0, 0},
            {1, ny, nz + 1},
            "H",
            -(FPNumber)1.0*dt,
            "E");
    void* J_update_params = yee.ConstructParams_A_equal_func_r_t(
        "JUpdater"
    );
    yee.AddUpdateInstruction("J-update", FDInstructionCode::A_equal_func_r_t, J_update_params);
    yee.AddUpdateInstruction("Ex-Hy-update", FDInstructionCode::A_plusequal_sum_b_C, Ex_Hy_update_params);
    yee.AddUpdateInstruction("Ex-Hz-update", FDInstructionCode::A_plusequal_sum_b_C, Ex_Hz_update_params);
    yee.AddUpdateInstruction("E-J-update", FDInstructionCode::A_plusequal_sum_b_C, E_J_update_params);
    yee.AddUpdateInstruction("Hy-update", FDInstructionCode::A_plusequal_sum_b_C, Hy_update_params);
    yee.AddUpdateInstruction("Hz-update", FDInstructionCode::A_plusequal_sum_b_C, Hz_update_params);
    yee.SetIterativeSequence({"J-update", "Ex-Hy-update", "Ex-Hz-update", "E-J-update", "Hy-update", "Hz-update"});
    yee.AddFullGridElementView("E-x", "E", 0);
    yee.AddFullGridElementView("H-y", "H", 1);
    yee.DeleteOlderViewFiles();
    yee.SetDataStoreRate("E-x", 1);
    yee.SetDataStoreRate("H-y", 1);
    yee.ApplyIterativeInstructions(401);
}



