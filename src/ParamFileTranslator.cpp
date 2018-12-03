

#include "YeeGrid.h"
#include "ParamFileTranslator.h"

ParamFileTranslator::ParamFileTranslator(const std::string filename) {
    ParamFileTranslator::filename = filename;
}

void ParamFileTranslator::Translate() {
    ParameterExtractor paramExtractor(filename);
    auto simulationType = paramExtractor.GetStringProperty("simulationType");
    if(simulationType == "singleGrid") {
        TranslateSingleGrid(paramExtractor.GetSubTreeRootNode("simulationParameters"));
    } else {
        assert(false);
    }
}

void ParamFileTranslator::TranslateSingleGrid(boost::property_tree::ptree node) {

    SingleGridParameterExtractor singleGridRoot(node);
    ParameterExtractor dimensionsExtractor(singleGridRoot.GetSubTreeRootNode("dimensions"));
    auto r0 = dimensionsExtractor.Get3VecRealProperty("r0");
    auto r1 = dimensionsExtractor.Get3VecRealProperty("r1");
    auto nCells = dimensionsExtractor.Get3VecUintProperty("nCells");

    RealNumber dt = dimensionsExtractor.GetRealProperty("dt");

    YeeGrid3D yee;
    yee.SetCornerCoordinates(r0, r1);
    yee.SetNumOfCells(nCells);
    yee.SetTimeResolution(dt);

    auto entireGridElements = singleGridRoot.GetEntireGridArrayParameters("entireGridArrays");

    for(auto& arrayNameType : entireGridElements) {
        yee.AddEntireGridElement(arrayNameType[0], stringToElementTypeMap[arrayNameType[1]]);
    }

    auto partialGridElements = singleGridRoot.GetPartialGridArrayParameters("partialGridArrays");
    for(auto& arrayParams : partialGridElements) {
        yee.AddPartialGridElement(std::get<0>(arrayParams),
                                  stringToElementTypeMap[std::get<1>(arrayParams)],
                                  std::get<2>(arrayParams),
                                  std::get<3>(arrayParams)
                                  );
    }

    auto gridArrayManipulators = singleGridRoot.GetTypeStringAndParameterSubtree("girdArrayManipulators");
    for(auto& manipulatorNameAndParams : gridArrayManipulators) {
        if(std::get<0>(manipulatorNameAndParams) == "GaussianGridArrayManipulator") {
            ParameterExtractor manipulatorParams(std::get<1>(manipulatorNameAndParams));
            yee.AddGaussianGridArrayManipulator(
                    manipulatorParams.GetStringProperty("name"),
                    manipulatorParams.GetStringProperty("array"),
                    stringDirectionToIntDirectionMap[manipulatorParams.GetStringProperty("direction")],
                    manipulatorParams.GetRealProperty("amplitude"),
                    manipulatorParams.GetRealProperty("t_center"),
                    manipulatorParams.GetRealProperty("t_decay"),
                    manipulatorParams.GetRealProperty("modulationFrequency"),
                    manipulatorParams.GetRealProperty("modulationPhase"),
                    manipulatorParams.GetRealProperty("timeOffsetFraction"));

        } else if(std::get<0>(manipulatorNameAndParams) == "GaussianSliceGridArrayManipulator") {
            ParameterExtractor manipulatorParams(std::get<1>(manipulatorNameAndParams));
            yee.AddGaussianGridArrayManipulator(
                    manipulatorParams.GetStringProperty("name"),
                    manipulatorParams.GetStringProperty("array"),
                    stringDirectionToIntDirectionMap[manipulatorParams.GetStringProperty("direction")],
                    manipulatorParams.Get3VecUintProperty("indStart"),
                    manipulatorParams.Get3VecUintProperty("indEnd"),
                    manipulatorParams.GetRealProperty("amplitude"),
                    manipulatorParams.GetRealProperty("t_center"),
                    manipulatorParams.GetRealProperty("t_decay"),
                    manipulatorParams.GetRealProperty("modulationFrequency"),
                    manipulatorParams.GetRealProperty("modulationPhase"),
                    manipulatorParams.GetRealProperty("timeOffsetFraction"));

        } else if(std::get<0>(manipulatorNameAndParams) == "GaussianSpaceTimeGridArrayManipulator") {
            ParameterExtractor manipulatorParams(std::get<1>(manipulatorNameAndParams));
            yee.AddGaussianSpaceTimeGridArrayManipulator(
                    manipulatorParams.GetStringProperty("name"),
                    manipulatorParams.GetStringProperty("array"),
                    stringDirectionToIntDirectionMap[manipulatorParams.GetStringProperty("direction")],
                    manipulatorParams.GetRealProperty("amplitude"),
                    manipulatorParams.Get4VecRealProperty("st_center"),
                    manipulatorParams.Get4VecRealProperty("st_decay_rate"),
                    manipulatorParams.Get4VecRealProperty("st_modulationFrequency"),
                    manipulatorParams.Get4VecRealProperty("st_modulationPhase"),
                    manipulatorParams.GetRealProperty("timeOffsetFraction"));

        } else if(std::get<0>(manipulatorNameAndParams) == "PeriodicGaussianGridArrayManipulator") {
            ParameterExtractor manipulatorParams(std::get<1>(manipulatorNameAndParams));

            auto primitiveVectors_vec = manipulatorParams.Get3VecRealArray("primitiveVectors");
            assert(primitiveVectors_vec.size() == 3);
            std::array<std::array<RealNumber, 3>, 3> primitiveVectors_array{primitiveVectors_vec[0],    // a 3-array
                                                                            primitiveVectors_vec[1],    // containing the
                                                                            primitiveVectors_vec[2]};   // 3 primitive vectors
            yee.AddPeriodicGaussianGridArrayManipulator(
                    manipulatorParams.GetStringProperty("name"),
                    manipulatorParams.GetStringProperty("array"),
                    stringDirectionToIntDirectionMap[manipulatorParams.GetStringProperty("direction")],
                    manipulatorParams.GetRealProperty("amplitude"),
                    manipulatorParams.Get3VecRealProperty("center"),
                    manipulatorParams.Get3VecRealProperty("decay_rate"),
                    manipulatorParams.Get3VecRealProperty("unitCellOrigin"),
                    primitiveVectors_array);

        } else if(std::get<0>(manipulatorNameAndParams) == "SpatialCubeGridArrayManipulator") {
            ParameterExtractor manipulatorParams(std::get<1>(manipulatorNameAndParams));
            yee.AddSpatialCubeGridArrayManipulator(
                    manipulatorParams.GetStringProperty("name"),
                    manipulatorParams.GetStringProperty("array"),
                    stringDirectionToIntDirectionMap[manipulatorParams.GetStringProperty("direction")],
                    manipulatorParams.Get3VecRealProperty("cornerR0"),
                    manipulatorParams.Get3VecRealProperty("cornerR1"),
                    manipulatorParams.Get3VecRealProperty("edgeThickness"),
                    manipulatorParams.GetRealProperty("valueInside"),
                    manipulatorParams.GetRealProperty("valueOutside")
                    );

        } else {
            assert(false);
        }
    }

    auto updateInstructions = singleGridRoot.GetTypeStringAndParameterSubtree("updateInstructions");
    for(auto& updateTypeandParams : updateInstructions) {
        auto& updateType = std::get<0>(updateTypeandParams);
        if(updateType == "A+=sumbC" || updateType == "A=sumbC") {
            ParameterExtractor updateParams(std::get<1>(updateTypeandParams));

            auto C_directions_str = updateParams.GetStringArray("C_direction");
            std::vector<int> C_directions;
            for(auto& direction : C_directions_str) {
                C_directions.emplace_back(stringDirectionToIntDirectionMap[direction]);
            }

            void* updateParamsPtr = yee.ConstructParams_A_plusequal_sum_b_C(
                updateParams.Get3VecUintProperty("A_indStart"),  // indStart
                updateParams.Get3VecUintProperty("A_indEnd"), // indEnd
                updateParams.GetStringProperty("A"),        // A name
                stringDirectionToIntDirectionMap[updateParams.GetStringProperty("A_direction")], // A direcction
                updateParams.GetRealArray("b"),     // b values
                updateParams.GetStringArray("C"),   // C names
                C_directions,
                updateParams.Get3VecUintArray("C_indStart")
                );

            if(updateType == "A+=sumbC") {
                yee.AddUpdateInstruction(updateParams.GetStringProperty("name"),
                        FDInstructionCode::A_plusequal_sum_b_C,
                        updateParamsPtr
                        );
            }else if(updateType == "A=sumbC") {
                yee.AddUpdateInstruction(updateParams.GetStringProperty("name"),
                        FDInstructionCode::A_equal_sum_b_C,
                        updateParamsPtr
                        );
            } else {
                assert(false);
            }

        } else if(updateType == "A+=sumbBC" || updateType == "A=sumbBC") {
            ParameterExtractor updateParams(std::get<1>(updateTypeandParams));

            auto B_directions_str = updateParams.GetStringArray("B_direction");
            std::vector<int> B_directions;
            for(auto& direction : B_directions_str) {
                B_directions.emplace_back(stringDirectionToIntDirectionMap[direction]);
            }

            auto C_directions_str = updateParams.GetStringArray("C_direction");
            std::vector<int> C_directions;
            for(auto& direction : C_directions_str) {
                C_directions.emplace_back(stringDirectionToIntDirectionMap[direction]);
            }

            void* updateParamsPtr = yee.ConstructParams_A_plusequal_sum_bB_C(
                updateParams.Get3VecUintProperty("A_indStart"),  // indStart
                updateParams.Get3VecUintProperty("A_indEnd"), // indEnd
                updateParams.GetStringProperty("A"),        // A name
                stringDirectionToIntDirectionMap[updateParams.GetStringProperty("A_direction")], // A direcction
                updateParams.GetRealArray("b"),     // b values
                updateParams.GetStringArray("B"),   // B names
                B_directions,
                updateParams.Get3VecUintArray("B_indStart"),
                updateParams.GetStringArray("C"),   // C names
                C_directions,
                updateParams.Get3VecUintArray("C_indStart")
                );

            if(updateType == "A+=sumbBC") {
                yee.AddUpdateInstruction(updateParams.GetStringProperty("name"),
                        FDInstructionCode::A_plusequal_sum_bB_C,
                        updateParamsPtr
                        );
            }else if(updateType == "A=sumbBC") {
                yee.AddUpdateInstruction(updateParams.GetStringProperty("name"),
                        FDInstructionCode::A_equal_sum_bB_C,
                        updateParamsPtr
                        );
            } else {
                assert(false);
            }

        } else if(updateType == "A=frt") {
            ParameterExtractor updateParams(std::get<1>(updateTypeandParams));

            void* updateParamsPtr = yee.ConstructParams_A_equal_func_r_t(
                    updateParams.GetStringProperty("girdArrayManipulator")
                    );

            yee.AddUpdateInstruction(updateParams.GetStringProperty("name"),
                    FDInstructionCode::A_equal_func_r_t,
                    updateParamsPtr
                    );
        } else {
            assert(false);
        }
    }

    auto updateSequences = singleGridRoot.GetUpdateSequences("updateSequences");
    for(auto& updateNameAndSequence : updateSequences) {
        yee.AddInstructionSequence(
                std::get<0>(updateNameAndSequence),
                std::get<1>(updateNameAndSequence)
                );
    };

    auto gridViews = singleGridRoot.GetTypeStringAndParameterSubtree("gridViews");
    for(auto& viewTypeandParams : gridViews) {
        auto& viewType = std::get<0>(viewTypeandParams);
        if(viewType == "entire") {
            ParameterExtractor viewParams(std::get<1>(viewTypeandParams));

            std::string viewName = viewParams.GetStringProperty("name");
            std::string arrayName = viewParams.GetStringProperty("array");
            int arrayDirection = stringDirectionToIntDirectionMap[viewParams.GetStringProperty("direction")];

            yee.AddFullGridElementView(viewName, arrayName, arrayDirection);

            std::size_t saveRate = viewParams.GetUintProperty("saveRate");
            yee.SetDataStoreRate(viewName, saveRate);
        } else if(viewType == "partial") {
            ParameterExtractor viewParams(std::get<1>(viewTypeandParams));

            std::string viewName = viewParams.GetStringProperty("name");
            std::string arrayName = viewParams.GetStringProperty("array");
            int arrayDirection = stringDirectionToIntDirectionMap[viewParams.GetStringProperty("direction")];
            auto indStart = viewParams.Get3VecUintProperty("indStart");
            auto indEnd = viewParams.Get3VecUintProperty("indEnd");

            yee.AddGridElementView(viewName, arrayName, arrayDirection, indStart, indEnd);

            std::size_t saveRate = viewParams.GetUintProperty("saveRate");
            yee.SetDataStoreRate(viewName, saveRate);
        } else {
            assert(false);
        }
    }

    yee.DeleteOlderViewFiles();

    auto runSequence = singleGridRoot.GetRunSequence("runSequence");
    for(auto& sequenceNameAndNumOfRuns : runSequence) {
        yee.ApplyInstructions(      //TODO: set time index at the beginning
                std::get<0>(sequenceNameAndNumOfRuns),
                std::get<1>(sequenceNameAndNumOfRuns),
                std::get<2>(sequenceNameAndNumOfRuns)
                );
    }
}

