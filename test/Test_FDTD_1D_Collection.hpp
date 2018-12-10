

#include "YeeGridCollection.h"


void test_yeegrid_1d_collection() {
    YeeGridCollection yeeCollection;
    std::size_t yeeLeftInd = yeeCollection.AddGrid();
    std::size_t yeeRightInd = yeeCollection.AddGrid();

    std::size_t nz = 300/2;
    std::size_t indJ = nz/2;
    std::array<FPNumber, 3> r0{0.0, 0.0, 0.0};
    std::array<FPNumber, 3> r1{0.1, 0.0, (FPNumber)5.0};
    std::array<std::size_t, 3> nCells{1, 1, nz};
    FPNumber dz = (r1[2] - r0[2])/(FPNumber)(nz);
    FPNumber dt = dz*(FPNumber)0.99;
    //------ left side
    YeeGrid3D& yee_l = yeeCollection.GetGrid(yeeLeftInd);
    YeeGrid3D& yee_r = yeeCollection.GetGrid(yeeRightInd);

    yee_l.SetCornerCoordinates(r0, r1);
    yee_l.SetNumOfCells(nCells);
    yee_l.SetTimeResolution(dt);
    yee_l.AddEntireGridElement("E", ElementType::EdgeE);
    yee_l.AddEntireGridElement("H", ElementType::EdgeH);
    yee_l.AddPartialGridElement("J", ElementType::EdgeE, {0, 0, indJ}, {1, 1, 0});
    yee_l.AddGaussianGridArrayManipulator("JUpdater", "J", 0 /*x*/, 1.0 /*amp*/, 1.0 /*t_center*/, 0.2 /*t_decay*/,
            0.0 /*modulation frequency*/, 0.0 /*modulation phase*/, -0.5 /*time offset fraction*/);

    void* E_update_params_l = FDInstructionFactory::Get_AxEdgeE_plusEqual_b_Curl_CyEdgeH(
            yee_l,    // grid
            {0, 0, 1},  // indStart
            {1, 1, nz}, // indEnd
            "E",        // A name
            dt,     // b value
            "H"         // C name
            );
    void* E_J_update_params_l = yee_l.ConstructParams_A_plusequal_sum_b_C(
        {0, 0, indJ},
        {1, 1, indJ+1},
        "E",
        0,
        {-dt/dz},
        {"J"},
        {0},
        {{0, 0, 0}}
    );

    void* H_update_params_l = FDInstructionFactory::Get_AyEdgeH_plusEqual_b_Curl_CxEdgeE(
            yee_l,
            {0, 0, 0},
            {1, 1, nz - 1},     // except the last H element
            "H",
            -dt,
            "E");
    void* H_update_params_nz_inside_l = yee_l.ConstructParams_A_plusequal_sum_b_C(
            {0, 0, nz - 1},  // indStart
            {1, 1, nz},  // indEnd
            "H",        // A name
            1,
            {+dt/dz},     // b values
            {"E"},           // C names
            {0},
            {{0, 0, nz - 1}}
            );
    void* H_update_params_nz_outside_l = yee_l.ConstructParams_A_plusequal_sum_b_C_neighbor(
            &yee_r,     // neighbor grid
            {0, 0, nz - 1},  // indStart
            {1, 1, nz},  // indEnd
            "H",        // A name
            1,
            {-dt/dz},     // b values
            {"E"},           // C names
            {0},
            {{0, 0, 0}}
            );
    void* J_update_params_l = yee_l.ConstructParams_A_equal_func_r_t(
        "JUpdater"
    );
    yee_l.AddUpdateInstruction("J-update", FDInstructionCode::A_equal_func_r_t, J_update_params_l);
    yee_l.AddUpdateInstruction("E-update", FDInstructionCode::A_plusequal_sum_b_C, E_update_params_l);
    yee_l.AddUpdateInstruction("E-J-update", FDInstructionCode::A_plusequal_sum_b_C, E_J_update_params_l);
    yee_l.AddUpdateInstruction("H-update", FDInstructionCode::A_plusequal_sum_b_C, H_update_params_l);
    yee_l.AddUpdateInstruction("H-update-last-inside", FDInstructionCode::A_plusequal_sum_b_C,
                                H_update_params_nz_inside_l);
    yee_l.AddUpdateInstruction("H-update-last-outside", FDInstructionCode::A_plusequal_sum_b_C_neighbor,
                                H_update_params_nz_outside_l);

    //-------- right side
    r0[2] += (FPNumber)5.0;
    r1[2] += (FPNumber)5.0;
    yee_r.SetCornerCoordinates(r0, r1);
    yee_r.SetNumOfCells(nCells);
    yee_r.SetTimeResolution(dt);
    yee_r.AddEntireGridElement("E", ElementType::EdgeE);
    yee_r.AddEntireGridElement("H", ElementType::EdgeH);

    void* E_update_params_r = FDInstructionFactory::Get_AxEdgeE_plusEqual_b_Curl_CyEdgeH(
            yee_l,    // grid
            {0, 0, 1},  // indStart
            {1, 1, nz}, // indEnd
            "E",        // A name
            dt,         // b value
            "H"         // C name
            );
    void* E_update_params_0_inside_r = yee_r.ConstructParams_A_plusequal_sum_b_C(
            {0, 0, 0},  // indStart
            {1, 1, 1},  // indEnd
            "E",        // A name
            0,
            {-dt/dz},        // b values
            {"H"},           // C names
            {1},
            {{0, 0, 0}}
            );
    void* E_update_params_0_outside_r = yee_r.ConstructParams_A_plusequal_sum_b_C_neighbor(
            &yee_l,     // neighbor grid
            {0, 0, 0},  // indStart
            {1, 1, 1},  // indEnd
            "E",        // A name
            0,
            {+dt/dz},     // b values
            {"H"},           // C names
            {1},
            {{0, 0, nz - 1}}
            );

    void* H_update_params_r = FDInstructionFactory::Get_AyEdgeH_plusEqual_b_Curl_CxEdgeE(
            yee_l,
            {0, 0, 0},
            {1, 1, nz},
            "H",
            -dt,
            "E");

    yee_r.AddUpdateInstruction("E-update", FDInstructionCode::A_plusequal_sum_b_C, E_update_params_r);
    yee_r.AddUpdateInstruction("H-update", FDInstructionCode::A_plusequal_sum_b_C, H_update_params_r);
    yee_r.AddUpdateInstruction("E-update-first-inside", FDInstructionCode::A_plusequal_sum_b_C,
                                E_update_params_0_inside_r);
    yee_r.AddUpdateInstruction("E-update-first-outside", FDInstructionCode::A_plusequal_sum_b_C_neighbor,
                                E_update_params_0_outside_r);



    yee_l.AddInstructionSequence("iterative-E-update", {"J-update", "E-update", "E-J-update"});
    yee_l.AddInstructionSequence("iterative-H-update", {"H-update", "H-update-last-inside", "H-update-last-outside"});

    yee_r.AddInstructionSequence("iterative-E-update", {"E-update", "E-update-first-inside", "E-update-first-outside"});
    yee_r.AddInstructionSequence("iterative-H-update", {"H-update"});


    yee_l.AddFullGridElementView("E-x-l", "E", 0);
    yee_l.AddFullGridElementView("H-y-l", "H", 1);
    yee_l.SetDataStoreRate("E-x-l", 1);
    yee_l.SetDataStoreRate("H-y-l", 1);
    yee_l.DeleteOlderViewFiles();

    yee_r.AddFullGridElementView("E-x-r", "E", 0);
    yee_r.AddFullGridElementView("H-y-r", "H", 1);
    yee_r.SetDataStoreRate("E-x-r", 1);
    yee_r.SetDataStoreRate("H-y-r", 1);
    yee_r.DeleteOlderViewFiles();

    yeeCollection.RunInstructionsPeriodically(0, 600, {std::pair<std::size_t, std::string>{0, "iterative-E-update"},
                                                       std::pair<std::size_t, std::string>{1, "iterative-E-update"},
                                                       std::pair<std::size_t, std::string>{0, "iterative-H-update"},
                                                       std::pair<std::size_t, std::string>{1, "iterative-H-update"}}
                                                       );
}


