


void test_yeegrid_dielectric_pml_1d() {
    std::size_t nz = 600;
    std::size_t indJ = nz/2;
    std::array<FPNumber, 3> r0{0.0, 0.0, 0.0};
    std::array<FPNumber, 3> r1{0.1, 0.1, 15.0};
    std::array<std::size_t, 3> nCells{1, 1, nz};
    FPNumber dz = (r1[2] - r0[2])/(FPNumber)(nz);
    FPNumber dt = dz*(FPNumber)0.95;
    YeeGrid3D yee;
    yee.SetCornerCoordinates(r0, r1);
    yee.SetNumOfCells(nCells);
    yee.SetTimeResolution(dt);
    yee.AddEntireGridElement("E", ElementType::EdgeE);
    yee.AddEntireGridElement("D", ElementType::EdgeE);
    yee.AddEntireGridElement("sig-E", ElementType::EdgeE);
    yee.AddEntireGridElement("EpsilonInv", ElementType::EdgeE);
    yee.AddEntireGridElement("H", ElementType::EdgeH);
    yee.AddEntireGridElement("sig-H", ElementType::EdgeH);
    yee.AddPartialGridElement("J", ElementType::EdgeE, {0, 0, indJ}, {1, 1, 0});
    yee.AddGaussianGridArrayManipulator("JUpdater", "J", 0 /*x*/, (FPNumber)1.0 /*amp*/, (FPNumber)1.0 /*t_center*/,
            (FPNumber)0.2 /*t_decay*/, (FPNumber)0.0 /*modulation frequency*/,
            (FPNumber)0.0 /*modulation phase*/, -(FPNumber)0.5 /*time offset fraction*/);
    yee.AddSpatialCubeGridArrayManipulator("EpsilonUpdater", "EpsilonInv", 0 /*x*/,
            {r0[0] - (FPNumber)0.1, r0[1] - (FPNumber)0.1, (FPNumber)9.0} /*box corner 0*/,
            {r1[0] + (FPNumber)0.1, r1[1] + (FPNumber)0.1, (FPNumber)15.0} /*box corner 1*/,
            {0.0, 0.0, 0.0} /*edge thickness*/,
            (FPNumber)1.0/(FPNumber)4.0 /*insideValue*/, (FPNumber)1.0 /*outsideValue*/);
    yee.AddSpatialCubeGridArrayManipulator("SigEUpdater", "sig-E", 2 /*z*/,
            {r0[0] - (FPNumber)0.1, r0[1] - (FPNumber)0.1, (FPNumber)4.01} /*box corner 0*/,
            {r1[0] + (FPNumber)0.1, r1[1] + (FPNumber)0.1, (FPNumber)11.01} /*box corner 1*/,
            {0.0, 0.0, (FPNumber)1.0} /*edge thickness*/,
            0.0 /*insideValue*/, +(FPNumber)1.0 /*outsideValue*/);
    yee.AddSpatialCubeGridArrayManipulator("SigHUpdater", "sig-H", 2 /*z*/,
            {r0[0] - (FPNumber)0.1, r0[1] - (FPNumber)0.1, (FPNumber)4.01} /*box corner 0*/,
            {r1[0] + (FPNumber)0.1, r1[1] + (FPNumber)0.1, (FPNumber)11.01} /*box corner 1*/,
            {0.0, 0.0, (FPNumber)1.0} /*edge thickness*/,
            0.0 /*insideValue*/, +(FPNumber)1.0 /*outsideValue*/);

    void* D_curlH_update_params = FDInstructionFactory::Get_AxEdgeE_plusEqual_b_Curl_CyEdgeH(
            yee,    // grid
            {0, 0, 1},  // indStart
            {1, 1, nz}, // indEnd
            "D",        // A name
            dt,         // b value
            "H"         // C name
            );
    void* D_sigz_update_params = yee.ConstructParams_A_plusequal_sum_bB_C(
        {0, 0, 1},
        {1, 1, nz},
        "D",
        0,
        {-dt},          // bValues
        {"sig-E"},
        {2},
        {{0, 0, 1}},
        {"D"},
        {0},
        {{0, 0, 1}}
    );
    void* E_D_update_params = yee.ConstructParams_A_plusequal_sum_bB_C(
        {0, 0, 1},
        {1, 1, nz},
        "E",
        0,
        {(FPNumber)1.0},          // bValues
        {"EpsilonInv"},
        {0},
        {{0, 0, 1}},
        {"D"},
        {0},
        {{0, 0, 1}}
    );
    void* D_J_update_params = yee.ConstructParams_A_plusequal_sum_b_C(
        {0, 0, indJ},
        {1, 1, indJ+1},
        "D",
        0,
        {-dt/dz},
        {"J"},
        {0},
        {{0, 0, 0}}
    );
    void* H_curlE_update_params = FDInstructionFactory::Get_AyEdgeH_plusEqual_b_Curl_CxEdgeE(
            yee,
            {0, 0, 0},
            {1, 1, nz},
            "H",
            -dt,
            "E");
    void* H_sigz_update_params = yee.ConstructParams_A_plusequal_sum_bB_C(
        {0, 0, 0},
        {1, 1, nz},
        "H",
        1,
        {-dt},          // bValues
        {"sig-H"},
        {2},
        {{0, 0, 0}},
        {"H"},
        {1},
        {{0, 0, 0}}
    );
    void* J_update_params = yee.ConstructParams_A_equal_func_r_t(
        "JUpdater"
    );
    void* Epsilon_update_params = yee.ConstructParams_A_equal_func_r_t(
        "EpsilonUpdater"
    );
    void* SigE_update_params = yee.ConstructParams_A_equal_func_r_t(
        "SigEUpdater"
    );
    void* SigH_update_params = yee.ConstructParams_A_equal_func_r_t(
        "SigHUpdater"
    );
    yee.AddUpdateInstruction("Epsilon-update", FDInstructionCode::A_equal_func_r_t, Epsilon_update_params);
    yee.AddUpdateInstruction("SigE-update", FDInstructionCode::A_equal_func_r_t, SigE_update_params);
    yee.AddUpdateInstruction("SigH-update", FDInstructionCode::A_equal_func_r_t, SigH_update_params);
    yee.AddUpdateInstruction("J-update", FDInstructionCode::A_equal_func_r_t, J_update_params);
    yee.AddUpdateInstruction("D-curlH-update", FDInstructionCode::A_plusequal_sum_b_C, D_curlH_update_params);
    yee.AddUpdateInstruction("D-sigz-update", FDInstructionCode::A_plusequal_sum_bB_C, D_sigz_update_params);
    yee.AddUpdateInstruction("D-J-update", FDInstructionCode::A_plusequal_sum_b_C, D_J_update_params);
    yee.AddUpdateInstruction("E-D-update", FDInstructionCode::A_equal_sum_bB_C, E_D_update_params);
    yee.AddUpdateInstruction("H-curlE-update", FDInstructionCode::A_plusequal_sum_b_C, H_curlE_update_params);
    yee.AddUpdateInstruction("H-sigz-update", FDInstructionCode::A_plusequal_sum_bB_C, H_sigz_update_params);
    yee.SetSingleRunSequence({"Epsilon-update", "SigE-update", "SigH-update"});
    yee.SetIterativeSequence({"J-update",
                              "D-sigz-update", "D-curlH-update", "D-J-update", "E-D-update",
                              "H-sigz-update", "H-curlE-update"});
    yee.AddFullGridElementView("E-x", "E", 0);
    yee.AddFullGridElementView("H-y", "H", 1);
    yee.AddFullGridElementView("Eps", "EpsilonInv", 0);
    yee.SetDataStoreRate("E-x", 1);
    yee.SetDataStoreRate("H-y", 1);
    yee.SetDataStoreRate("Eps", 1);
    yee.DeleteOlderViewFiles();
    yee.ApplySingleRunInstructions();
    yee.ApplyIterativeInstructions(600);
}




