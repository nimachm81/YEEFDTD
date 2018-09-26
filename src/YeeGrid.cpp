

#include <cstddef>      //std::size_t
#include <vector>
#include <memory>

#include "NumberTypes.h"
#include "YeeGrid.h"
#include "GaussianGridArrayManipulator.h"


YeeGrid3D::YeeGrid3D(std::array<std::size_t, 3>& nCells) : nCells(nCells) { }

YeeGrid3D::~YeeGrid3D() {
    for(std::string updateName : iterationSequence) {
        auto& instructCode_param_pair = instructions[updateName];
        auto& instructionCode = instructCode_param_pair.first;
        void* params = instructCode_param_pair.second;
        if(instructionCode == FDInstructionCode::A_plusequal_sum_b_C) {
            auto* params_tuple =
                static_cast<
                    std::tuple<
                        std::pair<std::array<std::size_t, 3>, std::array<std::size_t, 3>>,   // 0
                        std::string,    // 1
                        int,    // 2
                        std::vector<RealNumber>,    // 3
                        std::vector<std::string>,   // 4
                        std::vector<int>,           // 5
                        std::vector<std::array<int, 3>>     // 6
                    >*
                >(params);
            delete params_tuple;
        }
    }
}

void YeeGrid3D::SetCornerCoordinates(std::array<RealNumber, 3> r_0, std::array<RealNumber, 3> r_1) {
    YeeGrid3D::r_0 = r_0;
    YeeGrid3D::r_1 = r_1;
}

const std::array<std::size_t, 3>& YeeGrid3D::GetNumberOfCells() const {
    return nCells;
}

const std::array<RealNumber, 3>& YeeGrid3D::GetCornerR0() const {
    return r_0;
}

const std::array<RealNumber, 3>& YeeGrid3D::GetCornerR1() const {
    return r_1;
}


void YeeGrid3D::AddEntireGridElement(const std::string name, ElementType elemType) {
    gridElements[name] = std::make_shared<YeeGridData3D>(elemType, nCells);
}

void YeeGrid3D::AddPartialGridElement(const std::string name, ElementType elemType
        ,std::array<std::size_t, 3> startCell ,std::array<std::size_t, 3> numCells) {
    gridElements[name] = std::make_shared<YeeGridData3D>(elemType, numCells, startCell);
}


YeeGridData3D& YeeGrid3D::GetGridElement(const std::string name) {
    return *gridElements[name];
}

void YeeGrid3D::AddUpdateInstruction(const std::string name, FDInstructionCode instructionCode, void* params) {
    instructions[name] = std::pair<FDInstructionCode, void*>(instructionCode, params);
}

void YeeGrid3D::SetIterationSequence(std::vector<std::string> sequence) {
    iterationSequence = sequence;
}

void* YeeGrid3D::ConstructParams_A_plusequal_sum_b_C(
                                  std::array<std::size_t, 3> ind_start,
                                  std::array<std::size_t, 3> ind_end,
                                  std::string arrayA_name,
                                  int arrayA_component,
                                  std::vector<RealNumber> bValues,
                                  std::vector<std::string> arrayC_names,
                                  std::vector<int> arrayC_components,
                                  std::vector<std::array<int, 3>> arrayC_indexShifts
                                  ) {
    auto* params_tuple =
        new std::tuple<
            std::pair<std::array<std::size_t, 3>, std::array<std::size_t, 3>>,   // 0
            std::string,    // 1
            int,    // 2
            std::vector<RealNumber>,    // 3
            std::vector<std::string>,   // 4
            std::vector<int>,           // 5
            std::vector<std::array<int, 3>>     // 6
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
                    std::vector<std::array<int, 3>>     // 6
                >*
            >(params);
        std::array<std::size_t, 3>& ind_start = std::get<0>(params_tuple).first;
        std::array<std::size_t, 3>& ind_end = std::get<0>(params_tuple).second;
        std::string& arrayA_name = std::get<1>(params_tuple);
        int arrayA_component = std::get<2>(params_tuple);
        std::vector<RealNumber>& bValue = std::get<3>(params_tuple);
        std::vector<std::string>& arrayC_names = std::get<4>(params_tuple);
        std::vector<int>& arrayC_components = std::get<5>(params_tuple);
        std::vector<std::array<int, 3>>& arrayC_indexShifts = std::get<6>(params_tuple);

        std::size_t numRhs = bValue.size();
        assert(arrayC_names.size() == numRhs &&
               arrayC_components.size() == numRhs && arrayC_indexShifts.size() == numRhs);

        NumberArray3D<RealNumber>& arrayA = (*(gridElements[arrayA_name])).GetNumArray(arrayA_component);
        NumberArray3D<RealNumber> arrayASlice = arrayA.GetSlice(ind_start, ind_end);

        for(std::size_t i = 0; i < numRhs; ++i) {
            RealNumber b = bValue[i];
            NumberArray3D<RealNumber>& arrayC = (*(gridElements[arrayC_names[i]])).GetNumArray(arrayC_components[i]);
            std::array<int, 3>& arrayC_shift = arrayC_indexShifts[i];
            std::array<std::size_t, 3> ind_start_C{ind_start[0] + arrayC_shift[0],
                                                   ind_start[1] + arrayC_shift[1],
                                                   ind_start[2] + arrayC_shift[2]};
            std::array<std::size_t, 3> ind_end_C{ind_end[0] + arrayC_shift[0],
                                                 ind_end[1] + arrayC_shift[1],
                                                 ind_end[2] + arrayC_shift[2]};
            NumberArray3D<RealNumber> arrayCSlice = arrayC.GetSlice(ind_start_C, ind_end_C);

            arrayASlice += b*arrayCSlice;   // TODO : Do A += b*C in place, without creating a temp rhs
        }
    }
}


void YeeGrid3D::ApplyUpdateInstructions(std::size_t numIterations) {
    for(std::size_t i = 0; i < numIterations; ++i) {
        for(std::string updateName : iterationSequence) {
            auto& instructCode_param_pair = instructions[updateName];
            ApplyUpdateInstruction(instructCode_param_pair.first, instructCode_param_pair.second);
        }
    }
}


void YeeGrid3D::AddGaussianPointSource(const std::string name, const std::string gridDataName,
        int direction, std::array<std::size_t, 3> index, RealNumber amplitude,
        RealNumber t_center, RealNumber t_decay, RealNumber modulationFrequecy,
        RealNumber modulatioPhase) {
    std::shared_ptr<GaussianGridArrayManipulator> source(new GaussianGridArrayManipulator);
    source->SetDirection(direction);
    source->SetAmplitude(amplitude);
    source->SetCenterTime(t_center);
    source->SetDecayTime(t_decay);
    source->SetModulationFrequency(modulationFrequecy);
    source->SetModulationPhase(modulatioPhase);
    source->SetGridData(gridElements[gridDataName]);
    gridArrayManipulators[name] = source;
}


