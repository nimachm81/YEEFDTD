



void test_yeegrid_2d_pml() {
    std::size_t nz = 200;
    std::size_t ny = 200;
    std::size_t indzJx = nz/2;
    std::size_t indyJx = ny/2;
    std::array<FPNumber, 3> r0{0.0, 0.0, 0.0};
    std::array<FPNumber, 3> r1{(FPNumber)0.1, (FPNumber)10.0, (FPNumber)10.0};
    std::array<std::size_t, 3> nCells{1, ny, nz};
    FPNumber dz = (r1[2] - r0[2])/(FPNumber)(nz);
    FPNumber dy = (r1[1] - r0[1])/(FPNumber)(ny);
    FPNumber dt = dz/std::sqrt((FPNumber)2.0)*(FPNumber)0.95;
    YeeGrid3D yee;
    yee.SetCornerCoordinates(r0, r1);
    yee.SetNumOfCells(nCells);
    yee.SetTimeResolution(dt);
    yee.AddEntireGridElement("E", ElementType::EdgeE);
    yee.AddEntireGridElement("D", ElementType::EdgeE);
    yee.AddEntireGridElement("F", ElementType::EdgeE);
    yee.AddEntireGridElement("dF", ElementType::EdgeE);
    yee.AddEntireGridElement("sig-E", ElementType::EdgeE);

    yee.AddEntireGridElement("H", ElementType::EdgeH);
    yee.AddEntireGridElement("B", ElementType::EdgeH);
    yee.AddEntireGridElement("G", ElementType::EdgeH);
    yee.AddEntireGridElement("dG", ElementType::EdgeH);
    yee.AddEntireGridElement("sig-H", ElementType::EdgeH);

    yee.AddPartialGridElement("J", ElementType::EdgeE, {0, indyJx, indzJx}, {1, 0, 0});

    yee.AddGaussianGridArrayManipulator("JUpdater", "J", 0 /*x*/, (FPNumber)1.0 /*amp*/, (FPNumber)1.0 /*t_center*/,
            (FPNumber)0.2 /*t_decay*/, 0.0 /*modulation frequency*/, 0.0 /*modulation phase*/,
            -(FPNumber)0.5 /*time offset fraction*/);
    yee.AddSpatialCubeGridArrayManipulator("SigEzUpdater", "sig-E", 2 /*z*/,
            {r0[0] - (FPNumber)0.1, r0[1] - (FPNumber)0.1, (FPNumber)2.01} /*box corner 0*/,
            {r1[0] + (FPNumber)0.1, r1[1] + (FPNumber)0.1, (FPNumber)8.01} /*box corner 1*/,
            {0.0, 0.0, (FPNumber)1.0} /*edge thickness*/,
            0.0 /*insideValue*/, +(FPNumber)2.0 /*outsideValue*/);
    yee.AddSpatialCubeGridArrayManipulator("SigEyUpdater", "sig-E", 1 /*y*/,
            {r0[0] - (FPNumber)0.1, (FPNumber)2.01, r0[2] - (FPNumber)0.1} /*box corner 0*/,
            {r1[0] + (FPNumber)0.1, (FPNumber)8.01, r1[2] + (FPNumber)0.1} /*box corner 1*/,
            {0.0, 1.0, 0.0} /*edge thickness*/,
            0.0 /*insideValue*/, +(FPNumber)2.0 /*outsideValue*/);
    yee.AddSpatialCubeGridArrayManipulator("SigHzUpdater", "sig-H", 2 /*z*/,
            {r0[0] - (FPNumber)0.1, r0[1] - (FPNumber)0.1, (FPNumber)2.01} /*box corner 0*/,
            {r1[0] + (FPNumber)0.1, r1[1] + (FPNumber)0.1, (FPNumber)8.01} /*box corner 1*/,
            {0.0, 0.0, (FPNumber)1.0} /*edge thickness*/,
            0.0 /*insideValue*/, +2.0 /*outsideValue*/);
    yee.AddSpatialCubeGridArrayManipulator("SigHyUpdater", "sig-H", 1 /*y*/,
            {r0[0] - (FPNumber)0.1, (FPNumber)2.01, r0[2] - (FPNumber)0.1} /*box corner 0*/,
            {r1[0] + (FPNumber)0.1, (FPNumber)8.01, r1[2] + (FPNumber)0.1} /*box corner 1*/,
            {0.0, (FPNumber)1.0, 0.0} /*edge thickness*/,
            0.0 /*insideValue*/, +(FPNumber)2.0 /*outsideValue*/);
    void* dFx_sigzFx_update_params = yee.ConstructParams_A_plusequal_sum_bB_C(
        {0, 1, 1},
        {1, ny, nz},
        "dF",
        0,
        {-(FPNumber)1.0},          // bValues
        {"sig-E"},
        {2},
        {{0, 1, 1}},
        {"F"},
        {0},
        {{0, 1, 1}}
    );
    void* dFx_Hy_update_params = FDInstructionFactory::Get_AxEdgeE_plusEqual_b_Curl_CyEdgeH(  // A_x += b*curl C
                yee,
                {0, 1, 1},        // indStart
                {1, ny, nz},      // indEnd
                "dF",              // name of A
                +1.0,             // b value
                "H"                  // name of C
                );
    void* dFx_Hz_update_params = FDInstructionFactory::Get_AxEdgeE_plusEqual_b_Curl_CzEdgeH(  // A_x += b*curl C
                yee,
                {0, 1, 1},        // indStart
                {1, ny, nz},      // indEnd
                "dF",              // name of A
                +(FPNumber)1.0,             // b value
                "H"                  // name of C
                );
    void* dFx_Jx_update_params = yee.ConstructParams_A_plusequal_sum_b_C(
        {0, indyJx, indzJx},
        {1, indyJx + 1, indzJx + 1},
        "dF",    //
        0,
        {-(FPNumber)1.0/(dy*dz)},
        {"J"},
        {0},
        {{0, 0, 0}}
    );
    void* Fx_dFx_update_params = yee.ConstructParams_A_plusequal_sum_b_C(
        {0, 1, 1},
        {1, ny, nz},
        "F",
        0,
        {+dt},          // bValues
        {"dF"},
        {0},
        {{0, 1, 1}}
    );
    void* Dx_dFx_update_params = yee.ConstructParams_A_plusequal_sum_b_C(
        {0, 1, 1},
        {1, ny, nz},
        "D",
        0,
        {+dt},          // bValues
        {"dF"},
        {0},
        {{0, 1, 1}}
    );
    void* Dx_sigyDx_update_params = yee.ConstructParams_A_plusequal_sum_bB_C(
        {0, 1, 1},
        {1, ny, nz},
        "D",
        0,
        {-dt},          // bValues
        {"sig-E"},
        {1},            // y component
        {{0, 1, 1}},
        {"D"},
        {0},
        {{0, 1, 1}}
    );
    void* Ex_Dx_update_params = yee.ConstructParams_A_plusequal_sum_b_C(
        {0, 1, 1},
        {1, ny, nz},
        "E",
        0,
        {+(FPNumber)1.0},          // bValues
        {"D"},
        {0},
        {{0, 1, 1}}
    );
    void* dGy_Ex_update_params = FDInstructionFactory::Get_AyEdgeH_plusEqual_b_Curl_CxEdgeE(
            yee,
            {0, 0, 0},
            {1, ny + 1, nz},
            "dG",
            -(FPNumber)1.0,
            "E");
    void* dGz_Ex_update_params = FDInstructionFactory::Get_AzEdgeH_plusEqual_b_Curl_CxEdgeE(
            yee,
            {0, 0, 0},
            {1, ny, nz + 1},
            "dG",
            (FPNumber)(-1.0),
            "E");
    void* dGz_sigyGz_update_params = yee.ConstructParams_A_plusequal_sum_bB_C(
        {0, 0, 0},
        {1, ny, nz + 1},
        "dG",
        2,
        {(FPNumber)(-1.0)},          // bValues
        {"sig-H"},
        {1},
        {{0, 0, 0}},
        {"G"},
        {2},
        {{0, 0, 0}}
    );
    void* Bz_sigzGz_update_params = yee.ConstructParams_A_plusequal_sum_bB_C(
        {0, 0, 0},
        {1, ny, nz + 1},
        "B",
        2,
        {+dt},          // bValues
        {"sig-H"},
        {2},
        {{0, 0, 0}},
        {"G"},
        {2},
        {{0, 0, 0}}
    );
    void* Gz_dGz_update_params = yee.ConstructParams_A_plusequal_sum_b_C(
        {0, 0, 0},
        {1, ny, nz + 1},
        "G",
        2,
        {+dt},          // bValues
        {"dG"},
        {2},
        {{0, 0, 0}}
    );
    void* Gy_dGy_update_params = yee.ConstructParams_A_plusequal_sum_b_C(
        {0, 0, 0},
        {1, ny + 1, nz},
        "G",
        1,
        {+dt},          // bValues
        {"dG"},
        {1},
        {{0, 0, 0}}
    );
    void* Bz_dGz_update_params = yee.ConstructParams_A_plusequal_sum_b_C(
        {0, 0, 0},
        {1, ny, nz + 1},
        "B",
        2,
        {+dt},          // bValues
        {"dG"},
        {2},
        {{0, 0, 0}}
    );
    void* By_sigzBy_update_params = yee.ConstructParams_A_plusequal_sum_bB_C(
        {0, 0, 0},
        {1, ny, nz},    // TODO: ny+1 should be done separately
        "B",
        1,
        {-dt},          // bValues
        {"sig-H"},
        {2},
        {{0, 0, 0}},
        {"B"},
        {1},
        {{0, 0, 0}}
    );
    void* By_dGy_update_params = yee.ConstructParams_A_plusequal_sum_b_C(
        {0, 0, 0},
        {1, ny + 1, nz},
        "B",
        1,
        {+dt},          // bValues
        {"dG"},
        {1},
        {{0, 0, 0}}
    );
    void* Hy_By_update_params = yee.ConstructParams_A_plusequal_sum_b_C(
        {0, 0, 0},
        {1, ny + 1, nz},
        "H",
        1,
        {(FPNumber)1.0},          // bValues
        {"B"},
        {1},
        {{0, 0, 0}}
    );
    void* Hz_Bz_update_params = yee.ConstructParams_A_plusequal_sum_b_C(
        {0, 0, 0},
        {1, ny, nz + 1},
        "H",
        2,
        {(FPNumber)1.0},          // bValues
        {"B"},
        {2},
        {{0, 0, 0}}
    );
    void* J_update_params = yee.ConstructParams_A_equal_func_r_t(
        "JUpdater"
    );
    void* SigEz_update_params = yee.ConstructParams_A_equal_func_r_t(
        "SigEzUpdater"
    );
    void* SigEy_update_params = yee.ConstructParams_A_equal_func_r_t(
        "SigEyUpdater"
    );
    void* SigHz_update_params = yee.ConstructParams_A_equal_func_r_t(
        "SigHzUpdater"
    );
    void* SigHy_update_params = yee.ConstructParams_A_equal_func_r_t(
        "SigHyUpdater"
    );
    void* Ex_Fx_update_params = yee.ConstructParams_A_plusequal_sum_b_C(
        {0, 1, 1},
        {1, ny, nz},
        "E",
        0,
        {(FPNumber)1.0},          // bValues
        {"F"},
        {0},
        {{0, 1, 1}}
    );
    void* Hy_Gy_update_params = yee.ConstructParams_A_plusequal_sum_b_C(
        {0, 0, 0},
        {1, ny + 1, nz},
        "H",
        1,
        {(FPNumber)1.0},          // bValues
        {"G"},
        {1},
        {{0, 0, 0}}
    );
    void* Hz_Gz_update_params = yee.ConstructParams_A_plusequal_sum_b_C(
        {0, 0, 0},
        {1, ny, nz + 1},
        "H",
        2,
        {(FPNumber)1.0},          // bValues
        {"G"},
        {2},
        {{0, 0, 0}}
    );
    yee.AddUpdateInstruction("J-update", FDInstructionCode::A_equal_func_r_t, J_update_params);
    yee.AddUpdateInstruction("sigEz-update", FDInstructionCode::A_equal_func_r_t, SigEz_update_params);
    yee.AddUpdateInstruction("sigEy-update", FDInstructionCode::A_equal_func_r_t, SigEy_update_params);
    yee.AddUpdateInstruction("sigHz-update", FDInstructionCode::A_equal_func_r_t, SigHz_update_params);
    yee.AddUpdateInstruction("sigHy-update", FDInstructionCode::A_equal_func_r_t, SigHy_update_params);

    yee.AddUpdateInstruction("dFx-Jx-update", FDInstructionCode::A_plusequal_sum_b_C, dFx_Jx_update_params);
    yee.AddUpdateInstruction("dFx-Hy-update", FDInstructionCode::A_equal_sum_b_C, dFx_Hy_update_params);
    yee.AddUpdateInstruction("dFx-Hz-update", FDInstructionCode::A_plusequal_sum_b_C, dFx_Hz_update_params);
    yee.AddUpdateInstruction("dFx-sigz-Fx-update", FDInstructionCode::A_plusequal_sum_bB_C, dFx_sigzFx_update_params);
    yee.AddUpdateInstruction("Fx-dFx-update", FDInstructionCode::A_plusequal_sum_b_C, Fx_dFx_update_params);
    yee.AddUpdateInstruction("Dx-dFx-update", FDInstructionCode::A_plusequal_sum_b_C, Dx_dFx_update_params);
    yee.AddUpdateInstruction("Dx-sigy-Dx-update", FDInstructionCode::A_plusequal_sum_bB_C, Dx_sigyDx_update_params);
    yee.AddUpdateInstruction("Ex-Dx-update", FDInstructionCode::A_equal_sum_b_C, Ex_Dx_update_params);

    yee.AddUpdateInstruction("dGy-Ex-update", FDInstructionCode::A_equal_sum_b_C, dGy_Ex_update_params);
    yee.AddUpdateInstruction("dGz-Ex-update", FDInstructionCode::A_equal_sum_b_C, dGz_Ex_update_params);
    yee.AddUpdateInstruction("dGz-sigy-Gz-update", FDInstructionCode::A_plusequal_sum_bB_C, dGz_sigyGz_update_params);
    yee.AddUpdateInstruction("Gz-dGz-update", FDInstructionCode::A_plusequal_sum_b_C, Gz_dGz_update_params);
    yee.AddUpdateInstruction("Gy-dGy-update", FDInstructionCode::A_plusequal_sum_b_C, Gy_dGy_update_params);
    yee.AddUpdateInstruction("Bz-dGz-update", FDInstructionCode::A_plusequal_sum_b_C, Bz_dGz_update_params);
    yee.AddUpdateInstruction("Bz-sigz-Gz-update", FDInstructionCode::A_plusequal_sum_bB_C, Bz_sigzGz_update_params);
    yee.AddUpdateInstruction("By-dGy-update", FDInstructionCode::A_plusequal_sum_b_C, By_dGy_update_params);
    yee.AddUpdateInstruction("By-sigz-By-update", FDInstructionCode::A_plusequal_sum_bB_C, By_sigzBy_update_params);
    yee.AddUpdateInstruction("Hy-By-update", FDInstructionCode::A_equal_sum_b_C, Hy_By_update_params);
    yee.AddUpdateInstruction("Hz-Bz-update", FDInstructionCode::A_equal_sum_b_C, Hz_Bz_update_params);

    yee.AddUpdateInstruction("Ex-Fx-update", FDInstructionCode::A_equal_sum_b_C, Ex_Fx_update_params);
    yee.AddUpdateInstruction("Hy-Gy-update", FDInstructionCode::A_equal_sum_b_C, Hy_Gy_update_params);
    yee.AddUpdateInstruction("Hz-Gz-update", FDInstructionCode::A_equal_sum_b_C, Hz_Gz_update_params);

    yee.SetSingleRunSequence({"sigEz-update", "sigEy-update", "sigHy-update", "sigHz-update"});
    yee.SetIterativeSequence({"J-update",
                              "dFx-Hy-update", "dFx-Jx-update", "dFx-Hz-update", "dFx-sigz-Fx-update",
                              "Fx-dFx-update",
                              "Dx-dFx-update", "Dx-sigy-Dx-update",
                              "Ex-Dx-update",
                              "dGy-Ex-update", "dGz-Ex-update", "dGz-sigy-Gz-update",
                              "Gz-dGz-update", "Gy-dGy-update",
                              "Bz-dGz-update", "Bz-sigz-Gz-update",
                              "By-dGy-update", "By-sigz-By-update",
                              "Hy-By-update", "Hz-Bz-update"
                              });
    /*yee.SetIterativeSequence({"J-update",
                              "dFx-Hy-update", "dFx-Jx-update", "dFx-Hz-update", "dFx-sigz-Fx-update",
                              "Fx-dFx-update",
                              "Ex-Fx-update",
                              "dGy-Ex-update", "dGz-Ex-update", "dGz-sigy-Gz-update",
                              "Gz-dGz-update", "Gy-dGy-update",
                              "Hy-Gy-update", "Hz-Gz-update"
                              });*/
    yee.AddFullGridElementView("E-x", "E", 0);
    yee.AddFullGridElementView("H-y", "H", 1);
    yee.AddFullGridElementView("sigEz", "sig-E", 2);
    yee.AddFullGridElementView("sigEy", "sig-E", 1);
    yee.DeleteOlderViewFiles();
    yee.SetDataStoreRate("E-x", 1);
    yee.SetDataStoreRate("H-y", 1);
    yee.SetDataStoreRate("sigEz", 10);
    yee.SetDataStoreRate("sigEy", 10);
    yee.ApplySingleRunInstructions();
    yee.ApplyIterativeInstructions(401);
}


