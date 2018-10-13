
#ifndef FDTD_INSTRUCTIONFACTORY_H_
#define FDTD_INSTRUCTIONFACTORY_H_

#include <array>
#include <string>

#include "YeeGrid.h"

class FDInstructionFactory {
    public:
    // Ax type: EdgeE
    static void* Get_AxEdgeE_plusEqual_b_Curl_CyEdgeH(  // A_x += b*curl C_y     A_x type: EdgeE, C_y type: EdgeH
                const YeeGrid3D& yee,
                std::array<std::size_t, 3> indStart,    // A_x starts here
                std::array<std::size_t, 3> indEnd,      // A_x ends here (indEnd excluded)
                std::string aArrayName,
                RealNumber bNumber,
                std::string cArrayName
                );
    static void* Get_AxEdgeE_plusEqual_b_Curl_CzEdgeH(  // A_x += b*curl C     A_x type: EdgeE, C type: EdgeH
                const YeeGrid3D& yee,
                std::array<std::size_t, 3> indStart,
                std::array<std::size_t, 3> indEnd,
                std::string aArrayName,
                RealNumber bNumber,
                std::string cArrayName
                );
    // Ay type: EdgeE
    static void* Get_AyEdgeE_plusEqual_b_Curl_CxEdgeH(  // A_y += b*curl C_x     A_y type: EdgeE, C_x type: EdgeH
                const YeeGrid3D& yee,
                std::array<std::size_t, 3> indStart,
                std::array<std::size_t, 3> indEnd,
                std::string aArrayName,
                RealNumber bNumber,
                std::string cArrayName
                );
    static void* Get_AyEdgeE_plusEqual_b_Curl_CzEdgeH(  // A_y += b*curl C_z     A_y type: EdgeE, C_z type: EdgeH
                const YeeGrid3D& yee,
                std::array<std::size_t, 3> indStart,
                std::array<std::size_t, 3> indEnd,
                std::string aArrayName,
                RealNumber bNumber,
                std::string cArrayName
                );
    // Az type: EdgeE
    static void* Get_AzEdgeE_plusEqual_b_Curl_CxEdgeH(  // A_z += b*curl C_x     A_z type: EdgeE, C_x type: EdgeH
                const YeeGrid3D& yee,
                std::array<std::size_t, 3> indStart,
                std::array<std::size_t, 3> indEnd,
                std::string aArrayName,
                RealNumber bNumber,
                std::string cArrayName
                );
    static void* Get_AzEdgeE_plusEqual_b_Curl_CyEdgeH(  // A_z += b*curl C_y     A_z type: EdgeE, C_z type: EdgeH
                const YeeGrid3D& yee,
                std::array<std::size_t, 3> indStart,
                std::array<std::size_t, 3> indEnd,
                std::string aArrayName,
                RealNumber bNumber,
                std::string cArrayName
                );

    // Ax type: EdgeH
    static void* Get_AxEdgeH_plusEqual_b_Curl_CyEdgeE(  // A_x += b*curl C_y     A_x type: EdgeH, C_y type: EdgeE
                const YeeGrid3D& yee,
                std::array<std::size_t, 3> indStart,
                std::array<std::size_t, 3> indEnd,
                std::string aArrayName,
                RealNumber bNumber,
                std::string cArrayName
                );
    static void* Get_AxEdgeH_plusEqual_b_Curl_CzEdgeE(  // A_x += b*curl C_z     A_x type: EdgeH, C_z type: EdgeE
                const YeeGrid3D& yee,
                std::array<std::size_t, 3> indStart,
                std::array<std::size_t, 3> indEnd,
                std::string aArrayName,
                RealNumber bNumber,
                std::string cArrayName
                );
    // Ay type: EdgeH
    static void* Get_AyEdgeH_plusEqual_b_Curl_CxEdgeE(  // A_y += b*curl C_x     A_y type: EdgeH, C_x type: EdgeE
                const YeeGrid3D& yee,
                std::array<std::size_t, 3> indStart,
                std::array<std::size_t, 3> indEnd,
                std::string aArrayName,
                RealNumber bNumber,
                std::string cArrayName
                );
    static void* Get_AyEdgeH_plusEqual_b_Curl_CzEdgeE(  // A_y += b*curl C_z     A_y type: EdgeH, C_z type: EdgeE
                const YeeGrid3D& yee,
                std::array<std::size_t, 3> indStart,
                std::array<std::size_t, 3> indEnd,
                std::string aArrayName,
                RealNumber bNumber,
                std::string cArrayName
                );
    // Az type: EdgeH
    static void* Get_AzEdgeH_plusEqual_b_Curl_CxEdgeE(  // A_z += b*curl C_x     A_z type: EdgeH, C_x type: EdgeE
                const YeeGrid3D& yee,
                std::array<std::size_t, 3> indStart,
                std::array<std::size_t, 3> indEnd,
                std::string aArrayName,
                RealNumber bNumber,
                std::string cArrayName
                );
    static void* Get_AzEdgeH_plusEqual_b_Curl_CyEdgeE(  // A_z += b*curl C_y     A_z type: EdgeH, C_y type: EdgeE
                const YeeGrid3D& yee,
                std::array<std::size_t, 3> indStart,
                std::array<std::size_t, 3> indEnd,
                std::string aArrayName,
                RealNumber bNumber,
                std::string cArrayName
                );

};


#endif // FDTD_INSTRUCTIONFACTORY_H_
