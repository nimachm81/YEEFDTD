

#include <cstddef>
#include <cstdio>
#include <cmath>
#include <vector>
#include <memory>

#include "NumberTypes.h"
#include "YeeGrid.h"
#include "GaussianGridArrayManipulator.h"
#include "SpatialCubeGridArrayManipulator.h"
#include "GaussianSpaceTimeGridArrayManipulator.h"
#include "PeriodicGaussianGridArrayManipulator.h"
#include "SpaceTimeCubeGridArrayManipulator.h"
#include "SpherialShellGaussianGridArrayManipulator.h"
#include "ChargedParticlesTracer.h"
#include "DiscretePointsGridArrayManipulator.h"
#include "BivalueGridArrayManipulator.h"
#include "DataTruncationGridArrayManipulator.h"
#include "WedgeGeometry.h"
#include "ChargedParticleEmitter.h"
#include "ManualChargedParticleEmitter.h"
#include "GaussianPlaneWaveGridArrayManipulator.h"
#include "RectPlaneWaveGridArrayManipulator.h"
#include "GaussianPlaneWaveVectorField.h"
#include "RectPlaneWaveVectorField.h"


YeeGrid3D::~YeeGrid3D() {
    CloseGridViewFiles();
    for(std::string updateName : iterativeSequence) {
        auto& instructCode_param_pair = updateInstructions[updateName];
        auto& instructionCode = instructCode_param_pair.first;
        void* params = instructCode_param_pair.second;
        if(instructionCode == FDInstructionCode::A_plusequal_sum_b_C ||
           instructionCode == FDInstructionCode::A_equal_sum_b_C ||
           instructionCode == FDInstructionCode::A_plusequal_sum_b_C_shifted) {
            auto* params_tuple =
                static_cast<
                    std::tuple<
                        std::pair<std::array<std::size_t, 3>, std::array<std::size_t, 3>>,   // 0
                        std::string,    // 1
                        int,    // 2
                        std::vector<FPNumber>,    // 3
                        std::vector<std::string>,   // 4
                        std::vector<int>,           // 5
                        std::vector<std::array<std::size_t, 3>>     // 6
                    >*
                >(params);
            delete params_tuple;
        } else if(instructionCode == FDInstructionCode::A_plusequal_sum_bB_C ||
           instructionCode == FDInstructionCode::A_equal_sum_bB_C) {
            auto* params_tuple =
                static_cast<
                    std::tuple<
                    std::pair<std::array<std::size_t, 3>, std::array<std::size_t, 3>>,   // 0
                    std::string,    // 1
                    int,    // 2
                    std::vector<FPNumber>,    // 3
                    std::vector<std::string>,   // 4
                    std::vector<int>,           // 5
                    std::vector<std::array<std::size_t, 3>>,    // 6
                    std::vector<std::string>,   // 7
                    std::vector<int>,           // 8
                    std::vector<std::array<std::size_t, 3>>     // 9
                    >*
                >(params);
            delete params_tuple;
        } else if(instructionCode == FDInstructionCode::A_equal_func_r_t ||
                 instructionCode == FDInstructionCode::A_plusequal_func_r_t ||
                 instructionCode == FDInstructionCode::A_multequal_func_r_t) {
            auto* params_tuple =
                static_cast<
                    std::tuple<
                        std::string
                    >*
                >(params);
            delete params_tuple;
        } else if(instructionCode == FDInstructionCode::A_plusequal_sum_b_C_neighbor ||
                  instructionCode == FDInstructionCode::A_equal_sum_b_C_neighbor) {
            auto* params_tuple =
                static_cast<
                    std::tuple<
                        YeeGrid3D*,     // 0
                        std::pair<std::array<std::size_t, 3>, std::array<std::size_t, 3>>,   // 1
                        std::string,    // 2
                        int,    // 3
                        std::vector<FPNumber>,    // 4
                        std::vector<std::string>,   // 5
                        std::vector<int>,           // 6
                        std::vector<std::array<std::size_t, 3>>     // 7
                    >*
                >(params);
            delete params_tuple;
        }
    }
}

void YeeGrid3D::SetCornerCoordinates(std::array<FPNumber, 3> r_0, std::array<FPNumber, 3> r_1) {
    YeeGrid3D::r_0 = r_0;
    YeeGrid3D::r_1 = r_1;
}

void YeeGrid3D::SetNumOfCells(std::array<std::size_t, 3>& nCells) {
    YeeGrid3D::nCells = nCells;
    for(int i = 0; i < 3; ++i) {
        assert(nCells[i] >= 0);
        if(nCells[i] > 0) {
            dr[i] = (r_1[i] - r_0[i]) / (FPNumber)(nCells[i]);
        } else {
            assert(r_1[i]==r_0[i]);
            dr[i] = 0.0;
        }
    }
}

void YeeGrid3D::SetTimeResolution(const FPNumber dt) {
    YeeGrid3D::dt = dt;
}

void YeeGrid3D::SetTimeIndex(const std::size_t ind) {
    timeIndex = ind;
}

FPNumber YeeGrid3D::GetTimeResolution() const {
    return dt;
}

FPNumber YeeGrid3D::GetSpaceResolution(int i) const {
    return dr[i];
}

const std::array<FPNumber, 3>& YeeGrid3D::GetSpaceResolution() const {
    return dr;
}

const std::array<std::size_t, 3>& YeeGrid3D::GetNumberOfCells() const {
    return nCells;
}

const std::array<FPNumber, 3>& YeeGrid3D::GetCornerR0() const {
    return r_0;
}

const std::array<FPNumber, 3>& YeeGrid3D::GetCornerR1() const {
    return r_1;
}


void YeeGrid3D::AddEntireGridElement(const std::string name, ElementType elemType) {
    auto found = gridElements.find(name);
    assert(found == gridElements.end()); // make sure name does not already exist.
    gridElements[name] = std::make_shared<YeeGridData3D>(elemType, nCells);
}

void YeeGrid3D::AddPartialGridElement(const std::string name, ElementType elemType
        ,std::array<std::size_t, 3> startCell ,std::array<std::size_t, 3> numCells) {
    auto found = gridElements.find(name);
    assert(found == gridElements.end()); // make sure name does not already exist.
    gridElements[name] = std::make_shared<YeeGridData3D>(elemType, numCells, startCell);
}

YeeGridData3D& YeeGrid3D::GetGridElement(const std::string name) {
    auto found = gridElements.find(name);
    assert(found != gridElements.end()); // make sure name is valid.
    return *gridElements[name];
}

void YeeGrid3D::AddUpdateInstruction(const std::string name, FDInstructionCode instructionCode, void* params) {
    auto found = updateInstructions.find(name);
    assert(found == updateInstructions.end()); // make sure name does not already exist.
    updateInstructions[name] = std::pair<FDInstructionCode, void*>(instructionCode, params);
}

void YeeGrid3D::SetIterativeSequence(std::vector<std::string> sequence) {
    for(auto instructionName : sequence) {
        // make sure sequence instruction names are valid.
        auto found = updateInstructions.find(instructionName);
        assert(found != updateInstructions.end());
    }
    iterativeSequence = sequence;
}

void YeeGrid3D::SetSingleRunSequence(std::vector<std::string> sequence) {
    for(auto instructionName : sequence) {
        // make sure sequence instruction names are valid.
        auto found = updateInstructions.find(instructionName);
        assert(found != updateInstructions.end());
    }
    singleRunSequence = sequence;
}

void YeeGrid3D::AddInstructionSequence(std::string name, std::vector<std::string> sequence) {
    auto found = instructionSequences.find(name);
    assert(found == instructionSequences.end());  // name is valid (not already taken)

    for(auto instructionName : sequence) {
        // make sure sequence instruction names are valid.
        auto found = updateInstructions.find(instructionName);
        if(found == updateInstructions.end()) {
            std::cout << "Instruction not defined: " << instructionName << std::endl;
            assert(false);
        }
    }
    instructionSequences[name] = sequence;
}

void* YeeGrid3D::ConstructParams_A_plusequal_sum_b_C(
                                  std::array<std::size_t, 3> ind_start_A,
                                  std::array<std::size_t, 3> ind_end_A,
                                  std::string arrayA_name,
                                  int arrayA_component,
                                  std::vector<FPNumber> bValues,
                                  std::vector<std::string> arrayC_names,
                                  std::vector<int> arrayC_components,
                                  std::vector<std::array<std::size_t, 3>> arrayC_indsStart
                                  ) const {
    auto* params_tuple =
        new std::tuple<
            std::pair<std::array<std::size_t, 3>, std::array<std::size_t, 3>>,   // 0
            std::string,    // 1
            int,    // 2
            std::vector<FPNumber>,    // 3
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

void* YeeGrid3D::ConstructParams_A_plusequal_sum_bB_C(
                                  std::array<std::size_t, 3> ind_start_A,
                                  std::array<std::size_t, 3> ind_end_A,
                                  std::string arrayA_name,
                                  int arrayA_component,
                                  std::vector<FPNumber> bValues,
                                  std::vector<std::string> arrayB_names,
                                  std::vector<int> arrayB_components,
                                  std::vector<std::array<std::size_t, 3>> arrayB_indsStart,
                                  std::vector<std::string> arrayC_names,
                                  std::vector<int> arrayC_components,
                                  std::vector<std::array<std::size_t, 3>> arrayC_indsStart
                                  ) const {
    auto* params_tuple =
        new std::tuple<
            std::pair<std::array<std::size_t, 3>, std::array<std::size_t, 3>>,   // 0
            std::string,    // 1
            int,    // 2
            std::vector<FPNumber>,    // 3
            std::vector<std::string>,   // 4
            std::vector<int>,           // 5
            std::vector<std::array<std::size_t, 3>>,    // 6
            std::vector<std::string>,   // 7
            std::vector<int>,           // 8
            std::vector<std::array<std::size_t, 3>>     // 9
        >(
            std::pair<std::array<std::size_t, 3>, std::array<std::size_t, 3>>(ind_start_A, ind_end_A),  // 0
            arrayA_name,        // 1
            arrayA_component,   // 2
            bValues,            // 3
            arrayB_names,       // 4
            arrayB_components,  // 5
            arrayB_indsStart,   // 6
            arrayC_names,       // 7
            arrayC_components,  // 8
            arrayC_indsStart    // 9
        );
     return static_cast<void*>(params_tuple);
}


void* YeeGrid3D::ConstructParams_A_equal_func_r_t(std::string gridManipulator_name) const {
    auto* params_tuple = new std::tuple<std::string>(
            gridManipulator_name    // 0
            );
    return static_cast<void*>(params_tuple);
}

void* YeeGrid3D::ConstructParams_A_plusequal_sum_b_C_neighbor(
                                  YeeGrid3D* neighborGrid,
                                  std::array<std::size_t, 3> ind_start_A,
                                  std::array<std::size_t, 3> ind_end_A,
                                  std::string arrayA_name,
                                  int arrayA_component,
                                  std::vector<FPNumber> bValues,
                                  std::vector<std::string> arrayC_names,
                                  std::vector<int> arrayC_components,
                                  std::vector<std::array<std::size_t, 3>> arrayC_indsStart
                                  ) const {
    auto* params_tuple =
        new std::tuple<
            YeeGrid3D*,             // 0
            std::pair<std::array<std::size_t, 3>, std::array<std::size_t, 3>>,   // 1
            std::string,    // 2
            int,    // 3
            std::vector<FPNumber>,    // 4
            std::vector<std::string>,   // 5
            std::vector<int>,           // 6
            std::vector<std::array<std::size_t, 3>>     // 7
        >(
            neighborGrid,       // 0
            std::pair<std::array<std::size_t, 3>, std::array<std::size_t, 3>>(ind_start_A, ind_end_A),  // 1
            arrayA_name,        // 2
            arrayA_component,   // 3
            bValues,            // 4
            arrayC_names,       // 5
            arrayC_components,  // 6
            arrayC_indsStart    // 7
        );
     return static_cast<void*>(params_tuple);
}

void YeeGrid3D::ApplyUpdateInstruction(FDInstructionCode instructionCode, void* params) {
    if(instructionCode == FDInstructionCode::A_plusequal_sum_b_C ||
       instructionCode == FDInstructionCode::A_equal_sum_b_C ||
       instructionCode == FDInstructionCode::A_multequal_sum_b_C ||
       instructionCode == FDInstructionCode::A_plusequal_sum_b_C_shifted) {
        auto& params_tuple =
            *static_cast<
                std::tuple<
                    std::pair<std::array<std::size_t, 3>, std::array<std::size_t, 3>>,   // 0
                    std::string,    // 1
                    int,    // 2
                    std::vector<FPNumber>,    // 3
                    std::vector<std::string>,   // 4
                    std::vector<int>,           // 5
                    std::vector<std::array<std::size_t, 3>>     // 6
                >*
            >(params);
        std::array<std::size_t, 3>& ind_start_A = std::get<0>(params_tuple).first;
        std::array<std::size_t, 3>& ind_end_A = std::get<0>(params_tuple).second;
        std::string& arrayA_name = std::get<1>(params_tuple);
        int arrayA_component = std::get<2>(params_tuple);
        std::vector<FPNumber>& bValue = std::get<3>(params_tuple);
        std::vector<std::string>& arrayC_names = std::get<4>(params_tuple);
        std::vector<int>& arrayC_components = std::get<5>(params_tuple);
        std::vector<std::array<std::size_t, 3>>& arrayC_indsStart = std::get<6>(params_tuple);

        std::size_t numRhs = bValue.size();
        assert(arrayC_names.size() == numRhs &&
               arrayC_components.size() == numRhs && arrayC_indsStart.size() == numRhs);

        NumberArray3D<FPNumber>& arrayA = (*(gridElements[arrayA_name])).GetNumArray(arrayA_component);
        NumberArray3D<FPNumber> arrayASlice = arrayA.GetSlice(ind_start_A, ind_end_A);

        for(std::size_t i = 0; i < numRhs; ++i) {
            FPNumber b = bValue[i];
            NumberArray3D<FPNumber>& arrayC = (*(gridElements[arrayC_names[i]])).GetNumArray(arrayC_components[i]);
            std::array<std::size_t, 3>& ind_start_C = arrayC_indsStart[i];
            std::array<std::size_t, 3> ind_end_C{ind_start_C[0] + ind_end_A[0] - ind_start_A[0],
                                                 ind_start_C[1] + ind_end_A[1] - ind_start_A[1],
                                                 ind_start_C[2] + ind_end_A[2] - ind_start_A[2]};
            NumberArray3D<FPNumber> arrayCSlice = arrayC.GetSlice(ind_start_C, ind_end_C);

            if(instructionCode == FDInstructionCode::A_plusequal_sum_b_C_shifted) {
                std::array<std::size_t, 3>& indOriginA = (*(gridElements[arrayA_name])).GetIndexOfOrigin();
                std::array<std::size_t, 3>& indOriginC = (*(gridElements[arrayC_names[i]])).GetIndexOfOrigin();
                std::array<std::size_t, 3> ind_start_A_rel{ind_start_A[0] + indOriginC[0] - indOriginA[0],
                                                           ind_start_A[1] + indOriginC[1] - indOriginA[1],
                                                           ind_start_A[2] + indOriginC[2] - indOriginA[2]};
                std::array<std::size_t, 3> ind_end_A_rel{ind_end_A[0] + indOriginC[0] - indOriginA[0],
                                                         ind_end_A[1] + indOriginC[1] - indOriginA[1],
                                                         ind_end_A[2] + indOriginC[2] - indOriginA[2]};
                arrayASlice.MakeThisASliceOf(arrayA.GetSlice(ind_start_A_rel, ind_end_A_rel));
            }

            if(instructionCode == FDInstructionCode::A_plusequal_sum_b_C ||
                    instructionCode == FDInstructionCode::A_plusequal_sum_b_C_shifted) {
                arrayASlice += b*arrayCSlice;   // TODO : Do A += b*C in place, without creating a temp rhs
            } else if(instructionCode == FDInstructionCode::A_equal_sum_b_C) {
                if(i == 0) {
                    arrayASlice = b*arrayCSlice;   // TODO : Do A = b*C in place, without creating a temp rhs
                } else {
                    arrayASlice += b*arrayCSlice;   // TODO : Do A = b*C in place, without creating a temp rhs
                }
            } else if(instructionCode == FDInstructionCode::A_multequal_sum_b_C) {
                if(i == 0 && numRhs == 1) {
                    arrayASlice *= b*arrayCSlice;   // TODO : Do A = b*C in place, without creating a temp rhs
                } else {
                    std::cout << "#### error: not implemented!" << std::endl;
                    assert(false);
                }
            }
        }
    } else if(instructionCode == FDInstructionCode::A_plusequal_sum_bB_C ||
       instructionCode == FDInstructionCode::A_equal_sum_bB_C) {
        auto& params_tuple =
            *static_cast<
                std::tuple<
                    std::pair<std::array<std::size_t, 3>, std::array<std::size_t, 3>>,   // 0
                    std::string,    // 1
                    int,    // 2
                    std::vector<FPNumber>,    // 3
                    std::vector<std::string>,   // 4
                    std::vector<int>,           // 5
                    std::vector<std::array<std::size_t, 3>>,    // 6
                    std::vector<std::string>,   // 7
                    std::vector<int>,           // 8
                    std::vector<std::array<std::size_t, 3>>     // 9
                >*
            >(params);
        std::array<std::size_t, 3>& ind_start_A = std::get<0>(params_tuple).first;
        std::array<std::size_t, 3>& ind_end_A = std::get<0>(params_tuple).second;
        std::string& arrayA_name = std::get<1>(params_tuple);
        int arrayA_component = std::get<2>(params_tuple);
        std::vector<FPNumber>& bValues = std::get<3>(params_tuple);
        std::vector<std::string>& arrayB_names = std::get<4>(params_tuple);
        std::vector<int>& arrayB_components = std::get<5>(params_tuple);
        std::vector<std::array<std::size_t, 3>>& arrayB_indsStart = std::get<6>(params_tuple);
        std::vector<std::string>& arrayC_names = std::get<7>(params_tuple);
        std::vector<int>& arrayC_components = std::get<8>(params_tuple);
        std::vector<std::array<std::size_t, 3>>& arrayC_indsStart = std::get<9>(params_tuple);

        std::size_t numRhs = bValues.size();
        assert(arrayB_names.size() == numRhs && arrayC_names.size() == numRhs &&
               arrayC_components.size() == numRhs && arrayC_indsStart.size() == numRhs &&
               arrayB_components.size() == numRhs && arrayB_indsStart.size() == numRhs);

        NumberArray3D<FPNumber>& arrayA = (*(gridElements[arrayA_name])).GetNumArray(arrayA_component);
        NumberArray3D<FPNumber> arrayASlice = arrayA.GetSlice(ind_start_A, ind_end_A);

        for(std::size_t i = 0; i < numRhs; ++i) {
            FPNumber b = bValues[i];
            NumberArray3D<FPNumber>& arrayB = (*(gridElements[arrayB_names[i]])).GetNumArray(arrayB_components[i]);
            std::array<std::size_t, 3>& ind_start_B = arrayB_indsStart[i];
            std::array<std::size_t, 3> ind_end_B{ind_start_B[0] + ind_end_A[0] - ind_start_A[0],
                                                 ind_start_B[1] + ind_end_A[1] - ind_start_A[1],
                                                 ind_start_B[2] + ind_end_A[2] - ind_start_A[2]};
            NumberArray3D<FPNumber> arrayBSlice = arrayB.GetSlice(ind_start_B, ind_end_B);

            NumberArray3D<FPNumber>& arrayC = (*(gridElements[arrayC_names[i]])).GetNumArray(arrayC_components[i]);
            std::array<std::size_t, 3>& ind_start_C = arrayC_indsStart[i];
            std::array<std::size_t, 3> ind_end_C{ind_start_C[0] + ind_end_A[0] - ind_start_A[0],
                                                 ind_start_C[1] + ind_end_A[1] - ind_start_A[1],
                                                 ind_start_C[2] + ind_end_A[2] - ind_start_A[2]};
            NumberArray3D<FPNumber> arrayCSlice = arrayC.GetSlice(ind_start_C, ind_end_C);

            if(instructionCode == FDInstructionCode::A_plusequal_sum_bB_C) {
                arrayASlice += b*arrayBSlice*arrayCSlice;   // TODO : Do A += B*C in place, without creating a temp rhs
            } else if(instructionCode == FDInstructionCode::A_equal_sum_bB_C) {
                if(i == 0) {
                    arrayASlice = b*arrayBSlice*arrayCSlice;   // TODO : Do A = B*C in place, without creating a temp rhs
                } else {
                    arrayASlice += b*arrayBSlice*arrayCSlice;   // TODO : Do A += B*C in place, without creating a temp rhs
                }
            }
        }
    } else if(instructionCode == FDInstructionCode::A_equal_func_r_t ||
              instructionCode == FDInstructionCode::A_plusequal_func_r_t ||
              instructionCode == FDInstructionCode::A_multequal_func_r_t) {
        auto& params_tuple =
            *static_cast<
                std::tuple<
                    std::string    // 0
                >*
            >(params);
        std::string& gridManipulator_name = std::get<0>(params_tuple);

        auto found = gridArrayManipulators.find(gridManipulator_name);
        assert(found != gridArrayManipulators.end());  // name is valid

        GridArrayManipulator& gridManipulator = *gridArrayManipulators[gridManipulator_name];
        FPNumber t = gridManipulator.CalculateTime(dt, timeIndex);
        if(instructionCode == FDInstructionCode::A_equal_func_r_t) {
            gridManipulator.UpdateArray(t, GAManipulatorInstructionCode::Equal);
        } else if(instructionCode == FDInstructionCode::A_plusequal_func_r_t) {
            gridManipulator.UpdateArray(t, GAManipulatorInstructionCode::PlusEqual);
        }else if(instructionCode == FDInstructionCode::A_multequal_func_r_t) {
            gridManipulator.UpdateArray(t, GAManipulatorInstructionCode::MultiplyEqual);
        }
    } else if(instructionCode == FDInstructionCode::A_plusequal_sum_b_C_neighbor ||
       instructionCode == FDInstructionCode::A_equal_sum_b_C_neighbor) {
        auto& params_tuple =
            *static_cast<
                std::tuple<
                    YeeGrid3D*,    // 0
                    std::pair<std::array<std::size_t, 3>, std::array<std::size_t, 3>>,   // 1
                    std::string,    // 2
                    int,    // 3
                    std::vector<FPNumber>,    // 4
                    std::vector<std::string>,   // 5
                    std::vector<int>,           // 6
                    std::vector<std::array<std::size_t, 3>>     // 7
                >*
            >(params);
        YeeGrid3D* neighborGrid = std::get<0>(params_tuple);
        std::array<std::size_t, 3>& ind_start_A = std::get<1>(params_tuple).first;
        std::array<std::size_t, 3>& ind_end_A = std::get<1>(params_tuple).second;
        std::string& arrayA_name = std::get<2>(params_tuple);
        int arrayA_component = std::get<3>(params_tuple);
        std::vector<FPNumber>& bValue = std::get<4>(params_tuple);
        std::vector<std::string>& arrayC_names = std::get<5>(params_tuple);
        std::vector<int>& arrayC_components = std::get<6>(params_tuple);
        std::vector<std::array<std::size_t, 3>>& arrayC_indsStart = std::get<7>(params_tuple);

        std::size_t numRhs = bValue.size();
        assert(arrayC_names.size() == numRhs &&
               arrayC_components.size() == numRhs && arrayC_indsStart.size() == numRhs);

        NumberArray3D<FPNumber>& arrayA = (*(gridElements[arrayA_name])).GetNumArray(arrayA_component);
        NumberArray3D<FPNumber> arrayASlice = arrayA.GetSlice(ind_start_A, ind_end_A);

        for(std::size_t i = 0; i < numRhs; ++i) {
            FPNumber b = bValue[i];
            NumberArray3D<FPNumber>& arrayC =
                    neighborGrid->GetGridElement(arrayC_names[i]).GetNumArray(arrayC_components[i]);
            std::array<std::size_t, 3>& ind_start_C = arrayC_indsStart[i];
            std::array<std::size_t, 3> ind_end_C{ind_start_C[0] + ind_end_A[0] - ind_start_A[0],
                                                 ind_start_C[1] + ind_end_A[1] - ind_start_A[1],
                                                 ind_start_C[2] + ind_end_A[2] - ind_start_A[2]};
            NumberArray3D<FPNumber> arrayCSlice = arrayC.GetSlice(ind_start_C, ind_end_C);

            if(instructionCode == FDInstructionCode::A_plusequal_sum_b_C_neighbor) {
                arrayASlice += b*arrayCSlice;   // TODO : Do A += b*C in place, without creating a temp rhs
            } else if(instructionCode == FDInstructionCode::A_equal_sum_b_C_neighbor) {
                if(i == 0) {
                    arrayASlice = b*arrayCSlice;   // TODO : Do A = b*C in place, without creating a temp rhs
                } else {
                    arrayASlice += b*arrayCSlice;   // TODO : Do A = b*C in place, without creating a temp rhs
                }
            }
        }
    } else {
        std::cout << "#### error: instruction code not implemented!" << std::endl;
        assert(false);
    }
}


void YeeGrid3D::ApplyIterativeInstructions(std::size_t numIterations) {
    for(std::size_t i = 0; i < numIterations; ++i) {
        timeIndex = i;
        ApplyIterativeInstructionsOnce();
        WriteAllGridElemViewsToFile();
        //PrintAllGridData();
    }
    //CloseGridViewFiles();
}

void YeeGrid3D::ApplyIterativeInstructionsOnce() {
    for(std::string& updateName : iterativeSequence) {
        auto& instructCode_param_pair = updateInstructions[updateName];
        ApplyUpdateInstruction(instructCode_param_pair.first, instructCode_param_pair.second);
    }
}

void YeeGrid3D::ApplySingleRunInstructions() {
    for(std::string& updateName : singleRunSequence) {
        auto& instructCode_param_pair = updateInstructions[updateName];
        ApplyUpdateInstruction(instructCode_param_pair.first, instructCode_param_pair.second);
    }
}

void YeeGrid3D::ApplyInstructions(std::string name, std::size_t timeIndStart, std::size_t timeIndEnd,
                                  bool writeToFile) {
    for(std::size_t i = timeIndStart; i < timeIndEnd; ++i) {
        timeIndex = i;
        ApplyInstructionsOnce(name);
        if(writeToFile) {
            WriteAllGridElemViewsToFile();
        }
    }
    //CloseGridViewFiles();
}

void YeeGrid3D::ApplyInstructionsOnce(std::string name) {
    auto found = instructionSequences.find(name);
    assert(found != instructionSequences.end());  // name is valid

    //std::cout << name << std::endl;

    std::vector<std::string>& sequence = instructionSequences[name];
    for(std::string& updateName : sequence) {
        auto& instructCode_param_pair = updateInstructions[updateName];
        //std::cout << updateName << std::endl;
        ApplyUpdateInstruction(instructCode_param_pair.first, instructCode_param_pair.second);
        //std::cout << "Done!" << std::endl;
    }
}

void YeeGrid3D::AddGaussianGridArrayManipulator(const std::string name, const std::string gridDataName,
        int direction, FPNumber amplitude,
        FPNumber t_center, FPNumber t_decay, FPNumber modulationFrequecy,
        FPNumber modulatioPhase, FPNumber timeOffsetFraction) {
    auto found = gridArrayManipulators.find(name);
    assert(found == gridArrayManipulators.end()); // make sure name does not already exist.

    std::shared_ptr<GaussianGridArrayManipulator> modifier(new GaussianGridArrayManipulator);
    modifier->SetAmplitude(amplitude);
    modifier->SetCenterTime(t_center);
    modifier->SetDecayTime(t_decay);
    modifier->SetModulationFrequency(modulationFrequecy);
    modifier->SetModulationPhase(modulatioPhase);
    modifier->SetGridArrayTo(gridElements[gridDataName]->GetNumArray(direction));
    modifier->SetTimeOffsetFraction(timeOffsetFraction);
    gridArrayManipulators[name] = modifier;
}

void YeeGrid3D::AddGaussianGridArrayManipulator(const std::string name,
        const std::string gridDataName,
        int direction,
        std::array<std::size_t, 3> indStart,
        std::array<std::size_t, 3> indEnd,
        FPNumber amplitude,
        FPNumber t_center, FPNumber t_decay, FPNumber modulationFrequecy,
        FPNumber modulatioPhase, FPNumber timeOffsetFraction
        ) {
    auto found = gridArrayManipulators.find(name);
    assert(found == gridArrayManipulators.end()); // make sure name does not already exist.

    std::shared_ptr<GaussianGridArrayManipulator> modifier(new GaussianGridArrayManipulator);
    modifier->SetAmplitude(amplitude);
    modifier->SetCenterTime(t_center);
    modifier->SetDecayTime(t_decay);
    modifier->SetModulationFrequency(modulationFrequecy);
    modifier->SetModulationPhase(modulatioPhase);
    modifier->SetGridArrayTo(gridElements[gridDataName]->GetNumArray(direction).GetSlice(indStart, indEnd));
    modifier->SetTimeOffsetFraction(timeOffsetFraction);
    gridArrayManipulators[name] = modifier;
}


void YeeGrid3D::AddSpatialCubeGridArrayManipulator(const std::string name,
        const std::string gridDataName,
        int direction,
        std::array<FPNumber, 3> boxCornerR0, std::array<FPNumber, 3> boxCornerR1,
        std::array<FPNumber, 3> edgeThickness,
        FPNumber insideValue, FPNumber outsideValue
        ) {
    auto found = gridArrayManipulators.find(name);
    assert(found == gridArrayManipulators.end()); // make sure name does not already exist.

    std::shared_ptr<SpatialCubeGridArrayManipulator> modifier(new SpatialCubeGridArrayManipulator);
    modifier->SetCubeCorners(boxCornerR0, boxCornerR1);
    modifier->SetEdgeThickness(edgeThickness);
    modifier->SetInsideValue(insideValue);
    modifier->SetOutsideValue(outsideValue);

    // find the coordinates of the first element of the array
    std::array<FPNumber, 3> arrayR0 = GetCoordinatesOfFirstElementOfGridDataArray(gridDataName, direction);

    modifier->SetCornerCoordinate(arrayR0);
    modifier->SetGridSpacing(dr);
    modifier->SetGridArrayTo(gridElements[gridDataName]->GetNumArray(direction));
    gridArrayManipulators[name] = modifier;
}

void YeeGrid3D::AddSpaceTimeCubeGridArrayManipulator(const std::string name,
        const std::string gridDataName,
        int direction,
        std::array<FPNumber, 4> boxCornerR0, std::array<FPNumber, 4> boxCornerR1,
        std::array<FPNumber, 4> edgeThickness,
        FPNumber insideValue, FPNumber outsideValue,
        FPNumber timeOffsetFraction
        ) {
    auto found = gridArrayManipulators.find(name);
    assert(found == gridArrayManipulators.end()); // make sure name does not already exist.

    std::shared_ptr<SpaceTimeCubeGridArrayManipulator> modifier(new SpaceTimeCubeGridArrayManipulator);
    modifier->SetCubeCorners(boxCornerR0, boxCornerR1);
    modifier->SetEdgeThickness(edgeThickness);
    modifier->SetInsideValue(insideValue);
    modifier->SetOutsideValue(outsideValue);
    modifier->SetTimeOffsetFraction(timeOffsetFraction);

    // find the coordinates of the first element of the array
    std::array<FPNumber, 3> arrayR0 = GetCoordinatesOfFirstElementOfGridDataArray(gridDataName, direction);

    modifier->SetCornerCoordinate(arrayR0);
    modifier->SetGridSpacing(dr);
    modifier->SetGridArrayTo(gridElements[gridDataName]->GetNumArray(direction));
    gridArrayManipulators[name] = modifier;
}

void YeeGrid3D::AddGaussianSpaceTimeGridArrayManipulator(const std::string name,
        const std::string gridDataName,
        int direction,
        FPNumber amplitude,
        std::array<FPNumber, 4> st_center,
        std::array<FPNumber, 4> st_decay_rate,
        std::array<FPNumber, 4> st_modulationFrequecy,
        std::array<FPNumber, 4> st_modulatioPhase,
        FPNumber timeOffsetFraction
        ) {
    auto found = gridArrayManipulators.find(name);
    assert(found == gridArrayManipulators.end()); // make sure name does not already exist.

    std::shared_ptr<GaussianSaceTimeGridArrayManipulator> modifier(new GaussianSaceTimeGridArrayManipulator);
    modifier->SetAmplitude(amplitude);
    modifier->SetSpaceTimeCenterPoint(st_center);
    modifier->SetSpeceTimeDecayRate(st_decay_rate);
    modifier->SetSpaceTimeModulationFrequency(st_modulationFrequecy);
    modifier->SetSpaceTimeModulationPhase(st_modulatioPhase);
    modifier->SetGridArrayTo(gridElements[gridDataName]->GetNumArray(direction));
    modifier->SetTimeOffsetFraction(timeOffsetFraction);

        // find the coordinates of the first element of the array
    std::array<FPNumber, 3> arrayR0 = GetCoordinatesOfFirstElementOfGridDataArray(gridDataName, direction);

    modifier->SetCornerCoordinate(arrayR0);
    modifier->SetGridSpacing(dr);

    gridArrayManipulators[name] = modifier;
}

void YeeGrid3D::AddPeriodicGaussianGridArrayManipulator(const std::string name,
        const std::string gridDataName,
        int direction,
        FPNumber amplitude,
        std::array<FPNumber, 3> center,
        std::array<FPNumber, 3> decay_rate,
        std::array<FPNumber, 3> unitCellOrigin,
        std::array<std::array<FPNumber, 3>, 3> primitiveVectors
        ) {
    auto found = gridArrayManipulators.find(name);
    assert(found == gridArrayManipulators.end()); // make sure name does not already exist.

    std::shared_ptr<PeriodicGaussianGridArrayManipulator> modifier(new PeriodicGaussianGridArrayManipulator);
    modifier->SetGaussianAmplitude(amplitude);
    modifier->SetGaussianCenter(center);
    modifier->SetGaussianDecayRate(decay_rate);
    modifier->SetUnitCellOrigin(unitCellOrigin);
    modifier->SetPrimitiveVectors(primitiveVectors[0], primitiveVectors[1], primitiveVectors[2]);
    modifier->SetGridArrayTo(gridElements[gridDataName]->GetNumArray(direction));

        // find the coordinates of the first element of the array
    std::array<FPNumber, 3> arrayR0 = GetCoordinatesOfFirstElementOfGridDataArray(gridDataName, direction);

    modifier->SetCornerCoordinate(arrayR0);
    modifier->SetGridSpacing(dr);

    gridArrayManipulators[name] = modifier;
}

void YeeGrid3D::AddSpherialShellGaussianGridArrayManipulator(const std::string name,
        const std::string gridDataName,
        int direction,
        FPNumber amplitude,
        std::array<FPNumber, 3> centerPoint,
        FPNumber radius,
        FPNumber r_decay_rate,
        FPNumber r_modulationFrequecy,
        FPNumber r_modulatioPhase
        ) {
    auto found = gridArrayManipulators.find(name);
    assert(found == gridArrayManipulators.end()); // make sure name does not already exist.

    std::shared_ptr<SpherialShellGaussianGridArrayManipulator> modifier(new SpherialShellGaussianGridArrayManipulator);
    modifier->SetAmplitude(amplitude);
    modifier->SetCenterPoint(centerPoint);
    modifier->SetRadius(radius);
    modifier->SetDecayRate(r_decay_rate);
    modifier->SetModulationFrequency(r_modulationFrequecy);
    modifier->SetModulationPhase(r_modulatioPhase);
    modifier->SetGridArrayTo(gridElements[gridDataName]->GetNumArray(direction));

        // find the coordinates of the first element of the array
    std::array<FPNumber, 3> arrayR0 = GetCoordinatesOfFirstElementOfGridDataArray(gridDataName, direction);

    modifier->SetCornerCoordinate(arrayR0);
    modifier->SetGridSpacing(dr);

    gridArrayManipulators[name] = modifier;
}

void YeeGrid3D::AddBivalueGridArrayManipulator(const std::string name,
        const std::string gridDataName,
        int direction,
        const std::string geometryName,
        FPNumber valueInside,
        FPNumber valueOutside
        ) {
    auto found = gridArrayManipulators.find(name);
    assert(found == gridArrayManipulators.end()); // make sure name does not already exist.

    std::shared_ptr<BivalueGridArrayManipulator> modifier(new BivalueGridArrayManipulator);

    auto foundGeometry = geometries.find(geometryName);
    assert(foundGeometry != geometries.end()); // geometry name exists

    modifier->SetGeometry(geometries[geometryName]);
    modifier->SetInsideValue(valueInside);
    modifier->SetOutsideValue(valueOutside);

    modifier->SetGridArrayTo(gridElements[gridDataName]->GetNumArray(direction));
    modifier->SetupConditionArray();

        // find the coordinates of the first element of the array
    std::array<FPNumber, 3> arrayR0 = GetCoordinatesOfFirstElementOfGridDataArray(gridDataName, direction);

    modifier->SetCornerCoordinate(arrayR0);
    modifier->SetGridSpacing(dr);

    gridArrayManipulators[name] = modifier;
}

void YeeGrid3D::AddGaussianPlaneWaveGridArrayManipulator(const std::string name,
        const std::string gridDataName,
        int direction,
        std::array<FPNumber, 3> propagationDirection,
        FPNumber velocity,
        FPNumber amplitude,
        FPNumber t_center,
        FPNumber t_decay_rate,
        FPNumber t_modulationFrequecy,
        FPNumber t_modulatioPhase,
        FPNumber timeOffsetFraction
        ) {
    auto found = gridArrayManipulators.find(name);
    assert(found == gridArrayManipulators.end()); // make sure name does not already exist.

    std::shared_ptr<GaussianPlaneWaveGridArrayManipulator> modifier(new GaussianPlaneWaveGridArrayManipulator);
    modifier->SetPropagationDirection(propagationDirection);
    modifier->SetPropagationVelocity(velocity);
    modifier->SetAmplitude(amplitude);
    modifier->SetCenterTime(t_center);
    modifier->SetTimeDecayRate(t_decay_rate);
    modifier->SetModulationFrequency(t_modulationFrequecy);
    modifier->SetModulationPhase(t_modulatioPhase);
    modifier->SetGridArrayTo(gridElements[gridDataName]->GetNumArray(direction));
    modifier->SetTimeOffsetFraction(timeOffsetFraction);

        // find the coordinates of the first element of the array
    std::array<FPNumber, 3> arrayR0 = GetCoordinatesOfFirstElementOfGridDataArray(gridDataName, direction);

    modifier->SetCornerCoordinate(arrayR0);
    modifier->SetGridSpacing(dr);

    gridArrayManipulators[name] = modifier;
}

void YeeGrid3D::AddRectPlaneWaveGridArrayManipulator(const std::string name,
            const std::string gridDataName,
            int direction,
            std::array<FPNumber, 3> propagationDirection,
            FPNumber velocity,
            FPNumber amplitude,
            FPNumber t_center,
            FPNumber t_rect_width,
            FPNumber t_edge_width,
            FPNumber t_modulationFrequecy,
            FPNumber t_modulatioPhase,
            FPNumber timeOffsetFraction
            ) {
    auto found = gridArrayManipulators.find(name);
    assert(found == gridArrayManipulators.end()); // make sure name does not already exist.

    std::shared_ptr<RectPlaneWaveGridArrayManipulator> modifier(new RectPlaneWaveGridArrayManipulator);
    modifier->SetPropagationDirection(propagationDirection);
    modifier->SetPropagationVelocity(velocity);
    modifier->SetAmplitude(amplitude);
    modifier->SetCenterTime(t_center);
    modifier->SetRectWidth(t_rect_width);
    modifier->SetRectEdgeWidth(t_edge_width);
    modifier->SetModulationFrequency(t_modulationFrequecy);
    modifier->SetModulationPhase(t_modulatioPhase);
    modifier->SetGridArrayTo(gridElements[gridDataName]->GetNumArray(direction));
    modifier->SetTimeOffsetFraction(timeOffsetFraction);

        // find the coordinates of the first element of the array
    std::array<FPNumber, 3> arrayR0 = GetCoordinatesOfFirstElementOfGridDataArray(gridDataName, direction);

    modifier->SetCornerCoordinate(arrayR0);
    modifier->SetGridSpacing(dr);

    gridArrayManipulators[name] = modifier;
}

void YeeGrid3D::AddDataTruncationGridArrayManipulator(const std::string name,
            const std::string gridDataName,
            int direction,
            FPNumber minValue,
            FPNumber maxValue,
            bool truncateMin,
            bool truncateMax
            ) {
    auto found = gridArrayManipulators.find(name);
    assert(found == gridArrayManipulators.end()); // make sure name does not already exist.

    std::shared_ptr<DataTruncationGridArrayManipulator> modifier(new DataTruncationGridArrayManipulator);

    modifier->SetMinValue(minValue);
    modifier->SetMaxValue(maxValue);
    modifier->TruncateDown(truncateMin);
    modifier->TruncateUp(truncateMax);

    modifier->SetGridArrayTo(gridElements[gridDataName]->GetNumArray(direction));

    // find the coordinates of the first element of the array
    std::array<FPNumber, 3> arrayR0 = GetCoordinatesOfFirstElementOfGridDataArray(gridDataName, direction);

    modifier->SetCornerCoordinate(arrayR0);
    modifier->SetGridSpacing(dr);

    gridArrayManipulators[name] = modifier;
}

void YeeGrid3D::AddWedgeGeometry(const std::string name,
            const FPNumber wedgeAngle,
            const FPNumber tipRadius,
            const FPNumber apexToBaseDistance,
            const std::array<FPNumber, 3> apexPosition,
            const bool closeBase
            ) {
    auto found = geometries.find(name);
    assert(found == geometries.end()); // make sure name does not already exist.

    std::shared_ptr<WedgeGeometry> geometry(new WedgeGeometry);

    geometry->SetWedgeAngle(wedgeAngle);
    geometry->SetTipRadius(tipRadius);
    geometry->SetApexToBaseDistance(apexToBaseDistance);
    geometry->SetApexPosition(apexPosition);
    geometry->CloseBase(closeBase);

    geometries[name] = geometry;
}

void YeeGrid3D::AddGaussianPlaneWaveVectorField(const std::string name,
            std::array<FPNumber, 3> propagationDirection,
            FPNumber velocity,
            std::array<FPNumber, 3> amplitude,
            FPNumber t_center,
            FPNumber t_decayRate,
            FPNumber t_modulationFrequency,
            FPNumber t_modulationPhase
            ) {
    auto found = vectorFields.find(name);
    assert(found == vectorFields.end()); // make sure name does not already exist.

    std::shared_ptr<GaussianPlaneWaveVectorField> field(new GaussianPlaneWaveVectorField);

    field->SetPropagationVelocity(velocity);
    field->SetPropagationDirection(propagationDirection);
    field->SetAmplitude(amplitude);
    field->SetCenterTime(t_center);
    field->SetTimeDecayRate(t_decayRate);
    field->SetModulationFrequency(t_modulationFrequency);
    field->SetModulationPhase(t_modulationPhase);


    vectorFields[name] = field;
}

void YeeGrid3D::AddRectPlaneWaveVectorField(const std::string name,
            std::array<FPNumber, 3> propagationDirection,
            FPNumber velocity,
            std::array<FPNumber, 3> amplitude,
            FPNumber t_center,
            FPNumber t_rectWidth,
            FPNumber t_edgeWidth,
            FPNumber t_modulationFrequency,
            FPNumber t_modulationPhase
            ) {
    auto found = vectorFields.find(name);
    assert(found == vectorFields.end()); // make sure name does not already exist.

    std::shared_ptr<RectPlaneWaveVectorField> field(new RectPlaneWaveVectorField);

    field->SetPropagationVelocity(velocity);
    field->SetPropagationDirection(propagationDirection);
    field->SetAmplitude(amplitude);
    field->SetCenterTime(t_center);
    field->SetRectWidth(t_rectWidth);
    field->SetRectEdgeWidth(t_edgeWidth);
    field->SetModulationFrequency(t_modulationFrequency);
    field->SetModulationPhase(t_modulationPhase);


    vectorFields[name] = field;
}

void YeeGrid3D::AddManualChargedParticleEmitter(const std::string name,
        FPNumber particleCharge,
        FPNumber particleMasse,
        const std::vector<FPNumber>& emissionTimes,        // particles are emitted at these times
        const std::vector<FPNumber>& emissionNumbers,      // number of particles emitted at emissionTimes[i]
        const std::vector<std::array<FPNumber, 3>>& particlesInitialPositions,
        const std::vector<std::array<FPNumber, 3>>& particlesInitialVelocities
        ) {
    auto found = particleEmitters.find(name);
    assert(found == particleEmitters.end()); // make sure name does not already exist.

    std::shared_ptr<ManualChargedParticleEmitter> emitter(new ManualChargedParticleEmitter);

    emitter->SetElementaryCharge(particleCharge);
    emitter->SetElementaryMass(particleMasse);
    emitter->SetEmissionTimes(emissionTimes);
    emitter->SetNumberOfParticlesToEmitAtEachSpecifiedTime(emissionNumbers);
    emitter->SetInitialPositions(particlesInitialPositions);
    emitter->SetInitialVelocities(particlesInitialVelocities);

    particleEmitters[name] = emitter;
}

void YeeGrid3D::AddChargedParticleEmitter(const std::string name,
            FPNumber particleCharge,
            FPNumber particleMasse,
            const std::string geometryName,
            int dimensions,
            FPNumber maxElemSize,
            const std::string eFieldName,
            const std::string analyticEFieldName,
            FPNumber unitLength,
            std::size_t numOfSubPoints
            ) {
    auto found = particleEmitters.find(name);
    assert(found == particleEmitters.end()); // make sure name does not already exist.

    std::shared_ptr<ChargedParticleEmitter> emitter(new ChargedParticleEmitter);

    emitter->SetElementaryCharge(particleCharge);
    emitter->SetElementaryMass(particleMasse);

    emitter->SetUnitLength(unitLength);

    //--- set up geometry surface
    auto foundGeometry = geometries.find(geometryName);
    assert(foundGeometry != geometries.end()); // make sure name does not already exist.

    std::vector<std::array<FPNumber, 3>> centerPoints;     // the point at the middle of each patch
    std::vector<std::array<FPNumber, 3>> normalVecs;   // normal unit vector pointing outwards
    std::vector<FPNumber> arcLenghts;

    if(dimensions == 2) {
        if(numOfSubPoints >= 1) {
            std::shared_ptr<std::vector<std::vector<std::array<FPNumber, 3>>>>
                    subPoints(new std::vector<std::vector<std::array<FPNumber, 3>>>{});
            auto* subPointsPtr = subPoints.get();
            geometries[geometryName]->SubdevideSurface2D(0.0,     // x_cut at x = 0.0
                                    maxElemSize,      // maximum length of each subdevision
                                    centerPoints,     // the point at the middle of each patch
                                    normalVecs,   // normal unit vector pointing outwards
                                    arcLenghts,      // patch area
                                    subPointsPtr,
                                    numOfSubPoints
                                    );
            emitter->SetEmissionSubPoints(subPoints);
        } else {
            geometries[geometryName]->SubdevideSurface2D(0.0,     // x_cut at x = 0.0
                                    maxElemSize,      // maximum length of each subdevision
                                    centerPoints,     // the point at the middle of each patch
                                    normalVecs,   // normal unit vector pointing outwards
                                    arcLenghts      // patch area
                                    );
        }
    } else {
        std::cout << "error: 3D subdivision not implemented" << std::endl;
        assert(false);
    }

    emitter->SetEmissionPoints(centerPoints);
    emitter->SetNormalVectors(normalVecs);
    emitter->SetSurfaceAreas(arcLenghts);

    // set up electric field
    auto foundEField = gridElements.find(eFieldName);
    assert(foundEField != gridElements.end());

    emitter->SetElectricField(gridElements[eFieldName].get());
    for(int direction = 0; direction < 3; ++direction) {
        std::array<FPNumber, 3> arrayR0 = GetCoordinatesOfFirstElementOfGridDataArray(eFieldName, direction);
        emitter->SetElectricFieldGridOrigin(direction, arrayR0);
    }
    emitter->SetGridSpacing(dr);

    if(analyticEFieldName != "") {
        auto foundAnalEField = vectorFields.find(analyticEFieldName);
        assert(foundAnalEField != vectorFields.end());

        emitter->SetAnalyticElectricField(vectorFields[analyticEFieldName].get());
    }

    particleEmitters[name] = emitter;
}

void YeeGrid3D::AddChargedParticlesTracer(const std::string name,
            const std::string eFieldName,
            const std::string bFieldName,
            const std::string analyticEFieldName,
            const std::string analyticBFieldName,
            const std::string srFieldName,
            const std::string particleEmitterName,
            const std::size_t numberOfReservedParticles,
            const std::size_t bunchSize,
            const std::string constrainingGeometryName,
            bool keepPointsInside
            ) {
    auto found = gamDataUpdaters.find(name);
    assert(found == gamDataUpdaters.end()); // make sure name does not already exist.

    std::shared_ptr<ChargedParticlesTracer> updater(new ChargedParticlesTracer);

    auto foundEField = gridElements.find(eFieldName);
    assert(foundEField != gridElements.end());
    updater->SetElectricFieldGrid(gridElements[eFieldName].get());

    auto foundBField = gridElements.find(bFieldName);
    assert(foundBField != gridElements.end());
    updater->SetMagneticFieldGrid(gridElements[bFieldName].get());

    if(analyticEFieldName != "") {
        auto foundAnalEField = vectorFields.find(analyticEFieldName);
        assert(foundAnalEField != vectorFields.end());
        updater->SetAnalyticElectricField(vectorFields[analyticEFieldName].get());
    }

    if(analyticBFieldName != "") {
        auto foundAnalBField = vectorFields.find(analyticBFieldName);
        assert(foundAnalBField != vectorFields.end());
        updater->SetAnalyticMagneticField(vectorFields[analyticBFieldName].get());
    }


    if(srFieldName != "") {
        auto foundSRField = gridElements.find(srFieldName);
        assert(foundSRField != gridElements.end());
        updater->SetScatteringRateFieldGrid(gridElements[srFieldName].get());
    }

    // add particles
    auto foundEmitter = particleEmitters.find(particleEmitterName);
    assert(foundEmitter != particleEmitters.end());
    updater->SetParticleEmitter(particleEmitters[particleEmitterName].get());

    if(numberOfReservedParticles > 0) {
        updater->ReserveMemory(numberOfReservedParticles);
    }

    if(bunchSize > 0 ) {
        updater->SetMaxChargedPartcileBunchSize(bunchSize);
    }

    // set field orgins
    for(int direction = 0; direction < 3; ++direction) {
        std::array<FPNumber, 3> arrayR0 = GetCoordinatesOfFirstElementOfGridDataArray(eFieldName, direction);
        updater->SetElectricFieldGridOrigin(direction, arrayR0);
    }
    for(int direction = 0; direction < 3; ++direction) {
        std::array<FPNumber, 3> arrayR0 = GetCoordinatesOfFirstElementOfGridDataArray(bFieldName, direction);
        updater->SetMagneticFieldGridOrigin(direction, arrayR0);
    }
    if(srFieldName != "") {
        for(int direction = 0; direction < 3; ++direction) {
            if(gridElements[srFieldName]->GetElemType() == ElementType::NodeScalar) {
                std::array<FPNumber, 3> arrayR0 = GetCoordinatesOfFirstElementOfGridDataArray(srFieldName, 0);
                updater->SetScatteringRateFieldGridOrigin(direction, arrayR0);
            } else {
                std::array<FPNumber, 3> arrayR0 = GetCoordinatesOfFirstElementOfGridDataArray(srFieldName, direction);
                updater->SetScatteringRateFieldGridOrigin(direction, arrayR0);
            }
        }
    }

    if(constrainingGeometryName != "") {
        auto foundGeometry = geometries.find(constrainingGeometryName);
        assert(foundGeometry != geometries.end());
        updater->SetConstrainingGeometry(geometries[constrainingGeometryName].get(), keepPointsInside);
    }

    // set grid spacing
    updater->SetGridSpacing(dr);

    gamDataUpdaters[name] = updater;
}


void YeeGrid3D::AddDiscretePointsGridArrayManipulator(const std::string name,
            const std::string gridDataName,
            int direction,
            const std::string dataUpdaterName,
            const std::string dataUpdaterDataName,      // name of the array inside the dataUpdater to associate with gridData
            int dataUpdaterDataDirection                // direction of the array inside dataUpdater
    ) {
    auto found = gridArrayManipulators.find(name);
    assert(found == gridArrayManipulators.end()); // make sure name does not already exist.

    std::shared_ptr<DiscretePointsGridArrayManipulator> modifier(new DiscretePointsGridArrayManipulator);

    // attach updater
    modifier->AddDataUpdater(gamDataUpdaters[dataUpdaterName].get(),
                             dataUpdaterDataName,
                             dataUpdaterDataDirection
                             );

    modifier->SetGridArrayTo(gridElements[gridDataName]->GetNumArray(direction));

    // Set origin
    std::array<FPNumber, 3> arrayR0 = GetCoordinatesOfFirstElementOfGridDataArray(gridDataName, direction);

    modifier->SetCornerCoordinate(arrayR0);
    modifier->SetGridSpacing(dr);

    gridArrayManipulators[name] = modifier;
}

std::array<FPNumber, 3> YeeGrid3D::GetCoordinatesOfFirstElementOfGridDataArray(
                                        const std::string& gridDataName,
                                        int direction) {
    std::array<std::size_t, 3>& indexOfOrigin = gridElements[gridDataName]->GetIndexOfOrigin();
    // find the coordinates of the first element of the array
    std::array<FPNumber, 3> arrayR0{r_0[0] + (FPNumber)(indexOfOrigin[0])*dr[0],
                                      r_0[1] + (FPNumber)(indexOfOrigin[1])*dr[1],
                                      r_0[2] + (FPNumber)(indexOfOrigin[2])*dr[2]};
    ElementType& elemType = gridElements[gridDataName]->GetElemType();
    if(elemType == ElementType::EdgeE) {
        arrayR0[direction] += dr[direction]/(FPNumber)2.0;
    } else if(elemType == ElementType::EdgeH) {
        for(int i = 0; i < 3; ++i) {
            if(i != direction) {
                arrayR0[i] += dr[i]/(FPNumber)2.0;
            }
        }
    } else if(elemType == ElementType::NodeScalar) {
        assert(direction == 0);
        // do nothing
    } else if(elemType == ElementType::NodeVector) {
        // do nothing
    } else{
        assert(false);
    }

    return arrayR0;
}




void YeeGrid3D::PrintAllGridData() {
    for(auto it = gridElements.begin(); it != gridElements.end(); ++it) {
        std::cout << it->first << std::endl << *(it->second);
    }
}

void YeeGrid3D::AddGridElementView(std::string gridElemViewName,   // name of the gridView
            std::string gridElemName , int gridElemComponent,   // name of the gridElement and its x,y,z component
            std::array<std::size_t, 3> indStart, std::array<std::size_t, 3> indEnd // slice start and end
            ) {
    gridElementViews.emplace(gridElemViewName, GridElementView());
    GridElementView& geView = gridElementViews[gridElemViewName];
    geView.SetName(gridElemViewName);
    geView.SetNumArray(gridElements[gridElemName]->GetNumArray(gridElemComponent).GetSlice(indStart, indEnd));
}

void YeeGrid3D::AddFullGridElementView(std::string gridElemViewName,   // name of the gridView
            std::string gridElemName , int gridElemComponent   // name of the gridElement and its x,y,z component
            ) {
    gridElementViews.emplace(gridElemViewName, GridElementView());
    GridElementView& geView = gridElementViews[gridElemViewName];
    geView.SetName(gridElemViewName);
    geView.SetNumArray(gridElements[gridElemName]->GetNumArray(gridElemComponent));
}

void YeeGrid3D::SetDataStoreRate(std::string gridElemViewName, std::size_t saveEveryNSammples) {
    gridElementViews[gridElemViewName].SetSaveOnDiskFrequency(saveEveryNSammples);
}

void YeeGrid3D::WriteAllGridElemViewsToFile() {
    for(auto it = gridElementViews.begin(); it != gridElementViews.end(); ++it) {
        it->second.StoreData(timeIndex);
    }
}

void YeeGrid3D::DeleteOlderViewFiles() {
    for(auto it = gridElementViews.begin(); it != gridElementViews.end(); ++it) {
        it->second.DeleteOlderFiles();
    }
}

void YeeGrid3D::CloseGridViewFiles() {
    for(auto it = gridElementViews.begin(); it != gridElementViews.end(); ++it) {
        it->second.CloseFile();
    }
}
