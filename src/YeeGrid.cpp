

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

void* YeeGrid3D::ConstructParams_A_plusequal_sum_b_C(
                                  std::array<std::size_t, 3> ind_start,
                                  std::array<std::size_t, 3> ind_end,
                                  std::string arrayA_name,
                                  int arrayA_component,
                                  std::vector<RealNumber> bValues,
                                  std::vector<std::string> arrayC_names,
                                  std::vector<int> arrayC_components,
                                  std::vector<std::array<std::size_t, 3>> arrayC_indexShifts
                                  ) {
    auto* params_tuple =
        new std::tuple<
            std::pair<std::array<std::size_t, 3>, std::array<std::size_t, 3>>,   // 0
            std::string,    // 1
            int,    // 2
            std::vector<RealNumber>,    // 3
            std::vector<std::string>,   // 4
            std::vector<int>,           // 5
            std::vector<std::array<std::size_t, 3>>     // 6
        >(
            std::pair<std::array<std::size_t, 3>, std::array<std::size_t, 3>>(ind_start, ind_end),  // 0
            arrayA_name,        // 1
            arrayA_component,   // 2
            bValues,            // 3
            arrayC_names,       // 4
            arrayC_components,  // 5
            arrayC_indexShifts  // 6
        );
     return static_cast<void*>(params_tuple);
}


void YeeGrid3D::ApplyUpdateInstruction(FDInstructionCode instructionCode, void* params) {
    if(instructionCode == FDInstructionCode::A_plusequal_sum_b_C) {
        auto& params_tuple =
            *static_cast<
                std::tuple<
                    std::pair<std::array<std::size_t, 3>, std::array<std::size_t, 3>>,   // 0
                    std::string,    // 1
                    int,    // 2
                    std::vector<RealNumber>,    // 3
                    std::vector<std::string>,   // 4
                    std::vector<int>,           // 5
                    std::vector<std::array<std::size_t, 3>>     // 6
                >*
            >(params);
        std::array<std::size_t, 3>& ind_start = std::get<0>(params_tuple).first;
        std::array<std::size_t, 3>& ind_end = std::get<0>(params_tuple).second;
        std::string& arrayA_name = std::get<1>(params_tuple);
        int arrayA_component = std::get<2>(params_tuple);
        std::vector<RealNumber>& bValue = std::get<3>(params_tuple);
        std::vector<std::string>& arrayC_names = std::get<4>(params_tuple);
        std::vector<int>& arrayC_components = std::get<5>(params_tuple);
        std::vector<std::array<std::size_t, 3>>& arrayC_indexShifts = std::get<6>(params_tuple);
    }
}


