

#include <cstddef>      //std::size_t
#include <vector>
#include <memory>

#include "NumberTypes.h"
#include "YeeGrid.h"
#include "GaussianGridArrayManipulator.h"


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
                        std::vector<std::array<std::size_t, 3>>     // 6
                    >*
                >(params);
            delete params_tuple;
        } else if(instructionCode == FDInstructionCode::A_equal_func_r_t) {
            auto* params_tuple =
                static_cast<
                    std::tuple<
                        std::string
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

void YeeGrid3D::SetNumOfCells(std::array<std::size_t, 3>& nCells) {
    YeeGrid3D::nCells = nCells;
    for(int i = 0; i < 3; ++i) {
        assert(nCells[i] >= 0);
        if(nCells[i] > 0) {
            dr[i] = (r_1[i] - r_0[i]) / nCells[i];
        } else {
            assert(r_1[i]==r_0[i]);
            dr[i] = 0.0;
        }
    }
}

void YeeGrid3D::SetTimeResolution(const RealNumber dt) {
    YeeGrid3D::dt = dt;
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
                                  std::array<std::size_t, 3> ind_start_A,
                                  std::array<std::size_t, 3> ind_end_A,
                                  std::string arrayA_name,
                                  int arrayA_component,
                                  std::vector<RealNumber> bValues,
                                  std::vector<std::string> arrayC_names,
                                  std::vector<int> arrayC_components,
                                  std::vector<std::array<std::size_t, 3>> arrayC_indsStart
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
            std::pair<std::array<std::size_t, 3>, std::array<std::size_t, 3>>(ind_start_A, ind_end_A),  // 0
            arrayA_name,        // 1
            arrayA_component,   // 2
            bValues,            // 3
            arrayC_names,       // 4
            arrayC_components,  // 5
            arrayC_indsStart    // 6
        );
     return static_cast<void*>(params_tuple);
}


void* YeeGrid3D::ConstructParams_A_equal_func_r_t(std::string gridManipulator_name) {
    auto* params_tuple = new std::tuple<std::string>(
            gridManipulator_name    // 0
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
        std::array<std::size_t, 3>& ind_start_A = std::get<0>(params_tuple).first;
        std::array<std::size_t, 3>& ind_end_A = std::get<0>(params_tuple).second;
        std::string& arrayA_name = std::get<1>(params_tuple);
        int arrayA_component = std::get<2>(params_tuple);
        std::vector<RealNumber>& bValue = std::get<3>(params_tuple);
        std::vector<std::string>& arrayC_names = std::get<4>(params_tuple);
        std::vector<int>& arrayC_components = std::get<5>(params_tuple);
        std::vector<std::array<std::size_t, 3>>& arrayC_indsStart = std::get<6>(params_tuple);

        std::size_t numRhs = bValue.size();
        assert(arrayC_names.size() == numRhs &&
               arrayC_components.size() == numRhs && arrayC_indsStart.size() == numRhs);

        NumberArray3D<RealNumber>& arrayA = (*(gridElements[arrayA_name])).GetNumArray(arrayA_component);
        NumberArray3D<RealNumber> arrayASlice = arrayA.GetSlice(ind_start_A, ind_end_A);

        for(std::size_t i = 0; i < numRhs; ++i) {
            RealNumber b = bValue[i];
            NumberArray3D<RealNumber>& arrayC = (*(gridElements[arrayC_names[i]])).GetNumArray(arrayC_components[i]);
            std::array<std::size_t, 3>& ind_start_C = arrayC_indsStart[i];
            std::array<std::size_t, 3> ind_end_C{ind_start_C[0] + ind_end_A[0] - ind_start_A[0],
                                                 ind_start_C[1] + ind_end_A[1] - ind_start_A[1],
                                                 ind_start_C[2] + ind_end_A[2] - ind_start_A[2]};
            NumberArray3D<RealNumber> arrayCSlice = arrayC.GetSlice(ind_start_C, ind_end_C);

            arrayASlice += b*arrayCSlice;   // TODO : Do A += b*C in place, without creating a temp rhs
        }
    } else if(instructionCode == FDInstructionCode::A_equal_func_r_t) {
        auto& params_tuple =
            *static_cast<
                std::tuple<
                    std::string    // 0
                >*
            >(params);
        std::string& gridManipulator_name = std::get<0>(params_tuple);
        GridArrayManipulator& gridManipulator = *gridArrayManipulators[gridManipulator_name];
        RealNumber t = gridManipulator.CalculateTime(dt, timeIndex);
        gridManipulator.UpdateArray(t);
    }
    PrintAllGridData();
}


void YeeGrid3D::ApplyUpdateInstructions(std::size_t numIterations) {
    for(std::size_t i = 0; i < numIterations; ++i) {
        timeIndex = i;
        for(std::string updateName : iterationSequence) {
            auto& instructCode_param_pair = instructions[updateName];
            ApplyUpdateInstruction(instructCode_param_pair.first, instructCode_param_pair.second);
        }
    }
}


void YeeGrid3D::AddGaussianPointSource(const std::string name, const std::string gridDataName,
        int direction, RealNumber amplitude,
        RealNumber t_center, RealNumber t_decay, RealNumber modulationFrequecy,
        RealNumber modulatioPhase, RealNumber timeOffsetFraction) {
    std::shared_ptr<GaussianGridArrayManipulator> source(new GaussianGridArrayManipulator);
    source->SetDirection(direction);
    source->SetAmplitude(amplitude);
    source->SetCenterTime(t_center);
    source->SetDecayTime(t_decay);
    source->SetModulationFrequency(modulationFrequecy);
    source->SetModulationPhase(modulatioPhase);
    source->SetGridData(gridElements[gridDataName]);
    source->SetTimeOffsetFraction(timeOffsetFraction);
    gridArrayManipulators[name] = source;
}

void YeeGrid3D::PrintAllGridData() {
    for(auto it = gridElements.begin(); it != gridElements.end(); ++it) {
        std::cout << it->first << std::endl << *(it->second);
    }
}
