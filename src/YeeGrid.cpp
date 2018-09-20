

#include <cstddef>      //std::size_t
#include <vector>
#include <memory>

#include "NumberTypes.h"
#include "YeeGrid.h"



YeeGrid3D::YeeGrid3D(std::array<std::size_t, 3>& nCells) : nCells(nCells) { }

void YeeGrid3D::AddGridElement(const std::string name, ElementType elemType) {
    gridElements[name] = std::make_unique<YeeGridData3D>(elemType, nCells);
}

YeeGridData3D& YeeGrid3D::GetGridElement(const std::string name) {
    return *gridElements[name];
}

void YeeGrid3D::AddUpdateInstruction(const std::string name, FDInstructionCode instructionCode, void* params) {
    instructions[name] = std::pair<FDInstructionCode, void*>(instructionCode, params);
}

void YeeGrid3D::ApplyUpdateInstruction(FDInstructionCode instructionCode, void* params) {
    if(instructionCode == FDInstructionCode::A_pe_b_Cijk) {
        auto& params_tuple =
                *static_cast<std::tuple<std::pair<std::array<std::size_t, 3>, std::array<std::size_t, 3>>,   // 0
                    std::string,    // 1
                    int,    // 2
                    std::vector<RealNumber>,    // 3
                    std::vector<std::string>,   // 4
                    std::vector<int>,           // 5
                    std::vector<std::array<std::size_t, 3>>     // 6
                    >*>(params);
    }
}


