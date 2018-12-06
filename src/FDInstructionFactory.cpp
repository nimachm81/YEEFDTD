

#include "NumberTypes.h"
#include "FDInstructionCode.h"
#include "FDInstructionFactory.h"


void* FDInstructionFactory::Get_AxEdgeE_plusEqual_b_Curl_CyEdgeH(
            const YeeGrid3D& yee,
            std::array<std::size_t, 3> indStart,
            std::array<std::size_t, 3> indEnd,
            std::string aArrayName,
            FPNumber bNumber,
            std::string cArrayName) {

    FPNumber dz = yee.GetSpaceResolution(2);
    void* aArray_update_params = yee.ConstructParams_A_plusequal_sum_b_C(
        indStart,       // ind_start_A
        indEnd,         // ind_end_A
        aArrayName,     // arrayA_name
        0,              // arrayA_component
        {-(FPNumber)1.0*bNumber/dz,    (FPNumber)1.0*bNumber/dz},                      // bValues
        {cArrayName,         cArrayName},                                  // arrayC_names
        {1,                  1},                                           // arraC_components
        {indStart,           {indStart[0], indStart[1], indStart[2] - 1}}  // arrayC_indsStart
    );
    return aArray_update_params;
}

void* FDInstructionFactory::Get_AxEdgeE_plusEqual_b_Curl_CzEdgeH(
            const YeeGrid3D& yee,
            std::array<std::size_t, 3> indStart,
            std::array<std::size_t, 3> indEnd,
            std::string aArrayName,
            FPNumber bNumber,
            std::string cArrayName
            ) {
    FPNumber dy = yee.GetSpaceResolution(1);
    void* aArray_update_params = yee.ConstructParams_A_plusequal_sum_b_C(
        indStart,      // ind_start_A
        indEnd,    // ind_end_A
        aArrayName,            // arrayA_name
        0,              // arrayA_component
        {+(FPNumber)1.0*bNumber/dy, -(FPNumber)1.0*bNumber/dy},   // bValues
        {cArrayName,      cArrayName},           // arrayC_names
        {2,               2},                    // arraC_components
        {indStart,        {indStart[0], indStart[1] - 1, indStart[2]}}     // arrayC_indsStart
    );
    return aArray_update_params;
}

void* FDInstructionFactory::Get_AyEdgeE_plusEqual_b_Curl_CxEdgeH(
            const YeeGrid3D& yee,
            std::array<std::size_t, 3> indStart,
            std::array<std::size_t, 3> indEnd,
            std::string aArrayName,
            FPNumber bNumber,
            std::string cArrayName) {

    FPNumber dz = yee.GetSpaceResolution(2);
    void* aArray_update_params = yee.ConstructParams_A_plusequal_sum_b_C(
        indStart,      // ind_start_A
        indEnd,    // ind_end_A
        aArrayName,            // arrayA_name
        1,              // arrayA_component
        {+(FPNumber)1.0*bNumber/dz, -(FPNumber)1.0*bNumber/dz},    // bValues
        {cArrayName,       cArrayName},     // arrayC_names
        {0,                0},         // arraC_components
        {indStart,         {indStart[0], indStart[1], indStart[2] - 1}}     // arrayC_indsStart
    );
    return aArray_update_params;
}

void* FDInstructionFactory::Get_AyEdgeE_plusEqual_b_Curl_CzEdgeH(
            const YeeGrid3D& yee,
            std::array<std::size_t, 3> indStart,
            std::array<std::size_t, 3> indEnd,
            std::string aArrayName,
            FPNumber bNumber,
            std::string cArrayName) {

    FPNumber dx = yee.GetSpaceResolution(0);
    void* aArray_update_params = yee.ConstructParams_A_plusequal_sum_b_C(
        indStart,      // ind_start_A
        indEnd,    // ind_end_A
        aArrayName,            // arrayA_name
        1,              // arrayA_component
        {-(FPNumber)1.0*bNumber/dx, +(FPNumber)1.0*bNumber/dx},    // bValues
        {cArrayName,       cArrayName},     // arrayC_names
        {2,         2},         // arraC_components
        {indStart, {indStart[0] - 1, indStart[1], indStart[2]}}     // arrayC_indsStart
    );
    return aArray_update_params;
}

void* FDInstructionFactory::Get_AzEdgeE_plusEqual_b_Curl_CxEdgeH(
            const YeeGrid3D& yee,
            std::array<std::size_t, 3> indStart,
            std::array<std::size_t, 3> indEnd,
            std::string aArrayName,
            FPNumber bNumber,
            std::string cArrayName) {

    FPNumber dy = yee.GetSpaceResolution(1);
    void* aArray_update_params = yee.ConstructParams_A_plusequal_sum_b_C(
        indStart,      // ind_start_A
        indEnd,    // ind_end_A
        aArrayName,            // arrayA_name
        2,              // arrayA_component
        {-(FPNumber)1.0*bNumber/dy, +(FPNumber)1.0*bNumber/dy},    // bValues
        {cArrayName,        cArrayName},     // arrayC_names
        {0,           0},         // arraC_components
        {indStart,   {indStart[0], indStart[1] - 1, indStart[2]}}     // arrayC_indsStart
    );
    return aArray_update_params;
}

void* FDInstructionFactory::Get_AzEdgeE_plusEqual_b_Curl_CyEdgeH(
            const YeeGrid3D& yee,
            std::array<std::size_t, 3> indStart,
            std::array<std::size_t, 3> indEnd,
            std::string aArrayName,
            FPNumber bNumber,
            std::string cArrayName) {

    FPNumber dx = yee.GetSpaceResolution(0);
    void* aArray_update_params = yee.ConstructParams_A_plusequal_sum_b_C(
        indStart,      // ind_start_A
        indEnd,    // ind_end_A
        aArrayName,            // arrayA_name
        2,              // arrayA_component
        {+(FPNumber)1.0*bNumber/dx, -(FPNumber)1.0*bNumber/dx},    // bValues
        {cArrayName,        cArrayName},     // arrayC_names
        {1,          1},         // arraC_components
        {indStart,  {indStart[0] - 1, indStart[1], indStart[2]}}     // arrayC_indsStart
    );
    return aArray_update_params;
}

void* FDInstructionFactory::Get_AxEdgeH_plusEqual_b_Curl_CyEdgeE(
            const YeeGrid3D& yee,
            std::array<std::size_t, 3> indStart,
            std::array<std::size_t, 3> indEnd,
            std::string aArrayName,
            FPNumber bNumber,
            std::string cArrayName
            ) {

    FPNumber dz = yee.GetSpaceResolution(2);
    void* aArray_update_params = yee.ConstructParams_A_plusequal_sum_b_C(
        indStart,
        indEnd,
        aArrayName,
        0,
        {-(FPNumber)1.0*bNumber/dz,                 +(FPNumber)1.0*bNumber/dz},
        {cArrayName,                                  cArrayName},
        {1,                                           1},
        {{indStart[0], indStart[1], indStart[2] + 1}, indStart}
    );
    return aArray_update_params;
}

void* FDInstructionFactory::Get_AxEdgeH_plusEqual_b_Curl_CzEdgeE(
            const YeeGrid3D& yee,
            std::array<std::size_t, 3> indStart,
            std::array<std::size_t, 3> indEnd,
            std::string aArrayName,
            FPNumber bNumber,
            std::string cArrayName
            ) {

    FPNumber dy = yee.GetSpaceResolution(1);
    void* aArray_update_params = yee.ConstructParams_A_plusequal_sum_b_C(
        indStart,
        indEnd,
        aArrayName,
        0,
        {+(FPNumber)1.0*bNumber/dy,                  -(FPNumber)1.0*bNumber/dy},
        {cArrayName,                                   cArrayName},
        {2,                                            2},
        {{indStart[0], indStart[1] + 1, indStart[2]},  indStart}
    );
    return aArray_update_params;
}

void* FDInstructionFactory::Get_AyEdgeH_plusEqual_b_Curl_CxEdgeE(
            const YeeGrid3D& yee,
            std::array<std::size_t, 3> indStart,
            std::array<std::size_t, 3> indEnd,
            std::string aArrayName,
            FPNumber bNumber,
            std::string cArrayName
            ) {

    FPNumber dz = yee.GetSpaceResolution(2);
    void* aArray_update_params = yee.ConstructParams_A_plusequal_sum_b_C(
        indStart,       // ind_start_A
        indEnd,         // ind_end_A
        aArrayName,     // arrayA_name
        1,              // arrayA_component
        {+(FPNumber)1.0*bNumber/dz,                 -(FPNumber)1.0*bNumber/dz},      // bValues
        {cArrayName,                                  cArrayName},              // arrayC_names
        {0,                                           0},                       // arraC_components
        {{indStart[0], indStart[1], indStart[2] + 1}, indStart}                 // arrayC_indsStart
    );

    return aArray_update_params;
}

void* FDInstructionFactory::Get_AyEdgeH_plusEqual_b_Curl_CzEdgeE(
            const YeeGrid3D& yee,
            std::array<std::size_t, 3> indStart,
            std::array<std::size_t, 3> indEnd,
            std::string aArrayName,
            FPNumber bNumber,
            std::string cArrayName
            ) {
    FPNumber dx = yee.GetSpaceResolution(0);
    void* aArray_update_params = yee.ConstructParams_A_plusequal_sum_b_C(
        indStart,
        indEnd,
        aArrayName,
        1,
        {-(FPNumber)1.0*bNumber/dx,                  +(FPNumber)1.0*bNumber/dx},
        {cArrayName,                                   cArrayName},
        {2,                                            2},
        {{indStart[0] + 1, indStart[1], indStart[2]},  indStart}
    );
    return aArray_update_params;
}

void* FDInstructionFactory::Get_AzEdgeH_plusEqual_b_Curl_CxEdgeE(
            const YeeGrid3D& yee,
            std::array<std::size_t, 3> indStart,
            std::array<std::size_t, 3> indEnd,
            std::string aArrayName,
            FPNumber bNumber,
            std::string cArrayName
            ) {
    FPNumber dy = yee.GetSpaceResolution(1);
    void* aArray_update_params = yee.ConstructParams_A_plusequal_sum_b_C(
        indStart,
        indEnd,
        aArrayName,
        2,
        {-(FPNumber)1.0*bNumber/dy, +(FPNumber)1.0*bNumber/dy},
        {cArrayName, cArrayName},
        {0, 0},
        {{indStart[0], indStart[1] + 1, indStart[2]}, indStart}
    );
    return aArray_update_params;
}

void* FDInstructionFactory::Get_AzEdgeH_plusEqual_b_Curl_CyEdgeE(
            const YeeGrid3D& yee,
            std::array<std::size_t, 3> indStart,
            std::array<std::size_t, 3> indEnd,
            std::string aArrayName,
            FPNumber bNumber,
            std::string cArrayName
            ) {
    FPNumber dx = yee.GetSpaceResolution(0);
    void* aArray_update_params = yee.ConstructParams_A_plusequal_sum_b_C(
        indStart,
        indEnd,
        aArrayName,
        2,
        {+(FPNumber)1.0*bNumber/dx, -(FPNumber)1.0*bNumber/dx},
        {cArrayName,        cArrayName},
        {1,           1},
        {{indStart[0] + 1, indStart[1], indStart[2]},   indStart}
    );
    return aArray_update_params;
}



