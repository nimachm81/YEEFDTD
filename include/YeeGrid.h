
#ifndef FDTD_YEEGRID_H_
#define FDTD_YEEGRID_H_


#include <cstddef>      //std::size_t
#include <string>
#include <unordered_map>
#include <memory>

#include "NumberTypes.h"
#include "YeeGridDataTypes.h"


enum class FDInstructionCode {
    A_pe_b_Cijk   // A[i][j][k] += b*(C[i+i0][j+j0][k+k0] + C[i+i1][j+j1][k+k1] + ... ) fir ijk in the range
                  // ind_start...ind_end
                  // parameters : tuple{pair(ind_start, ind_end), b, C, vector<[i0,j0,k0]>} with types
                  // std::tuple<std::pair<std::array<int, 3>>, RealNumber, std::string, std::vector<std::array<RealNumber, 3>>>
                  //
};

class YeeGrid3D {
    public:
    YeeGrid3D(std::array<std::size_t, 3>& nCells);
    void AddGridElement(const std::string name, ElementType elemType);
    YeeGridData3D& GetGridElement(const std::string name);

    //void AddUpdateInstruction(instructionCode, parameters);


    private:
    std::array<std::size_t, 3> nCells;
    std::unordered_map<std::string, std::unique_ptr<YeeGridData3D>> gridElements;

};


#endif  // FDTD_YEEGRID_H_


