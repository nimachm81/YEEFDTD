
#ifndef FDTD_YEEGRID_H_
#define FDTD_YEEGRID_H_


#include <cstddef>      //std::size_t
#include <string>
#include <unordered_map>
#include <utility>
#include <memory>

#include "NumberTypes.h"
#include "YeeGridDataTypes.h"


enum class FDInstructionCode {
    A_pe_b_Cijk   // A[i][j][k] += b0*C[i+i0][j+j0][k+k0] + b1*C[i+i1][j+j1][k+k1] + ... for ijk in the range
                  // ind_start...ind_end
                  // parameters : tuple{pair(ind_start, ind_end), A, ind_xyz_A, vector<b>, vector<C>, vector<ind_xyz_C>,
                  //              vector<[i0,j0,k0]>} with types
                  // std::tuple<std::pair<std::array<std::size_t, 3>, std::array<std::size_t, 3>>,      // summation range
                  //            std::string,                    // A array names
                  //            int,                            // 0:A_x, 1:A_y, 2:A_z
                  //            std::vector<RealNumber>,        // b scalars
                  //            std::vector<std::string>,       // C arrays names
                  //            std::vector<int>,               // C arrays components : 0:C_x, 1:C_y, 2:C_z
                  //            std::vector<std::array<std::size_t, 3>>>    // [i0,j0,k0], ...
                  //
};

class YeeGrid3D {
    public:
    YeeGrid3D(std::array<std::size_t, 3>& nCells);
    void AddGridElement(const std::string name, ElementType elemType);
    YeeGridData3D& GetGridElement(const std::string name);

    void AddUpdateInstruction(const std::string name, FDInstructionCode instructionCode, void* params);
    void ApplyUpdateInstruction(FDInstructionCode instructionCode, void* params);

    private:
    std::array<std::size_t, 3> nCells;
    std::unordered_map<std::string, std::unique_ptr<YeeGridData3D>> gridElements;
    std::unordered_map<std::string, std::pair<FDInstructionCode, void*>> instructions;

};


#endif  // FDTD_YEEGRID_H_


