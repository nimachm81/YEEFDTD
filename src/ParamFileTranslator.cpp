
#include <string>
#include <map>

#include "YeeGrid.h"
#include "YeeGridCollection.h"
#include "ParamFileTranslator.h"

ParamFileTranslator::ParamFileTranslator(const std::string filename) {
    ParamFileTranslator::filename = filename;
}

void ParamFileTranslator::Translate() {
    ParameterExtractor paramExtractor(filename);
    auto simulationType = paramExtractor.GetStringProperty("simulationType");
    if(simulationType == "singleGrid") {
        TranslateSingleGrid(paramExtractor.GetSubTreeRootNode("simulationParameters"));
    } else if(simulationType == "gridCollection") {
        TranslateGridCollection(paramExtractor.GetSubTreeRootNode("simulationParameters"));
    } else {
        assert(false);
    }
}

void ParamFileTranslator::TranslateSingleGrid(boost::property_tree::ptree node) {

    SingleGridParameterExtractor singleGridRoot(node);
    YeeGrid3D yee;
    std::map<std::string, YeeGrid3D*> grids;     // not used. but it should be defined to pass to SetSingleGridUpddateInstructions

    SetSingleGridDimensions(yee, singleGridRoot);
    SetSingleGridGridArrays(yee, singleGridRoot);
    SetSingleGridGridArrayManipulators(yee, singleGridRoot);
    SetSingleGridUpddateInstructions(yee, singleGridRoot, grids);
    SetSingleGridUpdateSequences(yee, singleGridRoot);
    SetSingleGridGridViews(yee, singleGridRoot);
    SetAndRunSingleGridRunSequencs(yee, singleGridRoot);
}

void ParamFileTranslator::TranslateGridCollection(boost::property_tree::ptree node) {
    ParameterExtractor simulationParamsRoot(node);
    ParameterExtractor gridsRoot(simulationParamsRoot.GetSubTreeRootNode("grids"));     // TODO: assert "grids" exists

    std::size_t numOfGrids = gridsRoot.GetSize();
    YeeGridCollection gridCollection;
    // initialize grids
    std::map<std::string, YeeGrid3D*> gridsMap;
    std::map<std::string, std::size_t> gridsInds;
    for(auto& it : gridsRoot.GetPropertyTreeRoot()) {
        std::size_t gridInd = gridCollection.AddGrid();
        std::string gridName = it.first;
        gridsInds[gridName] = gridInd;
    }
    for(auto& it : gridsRoot.GetPropertyTreeRoot()) {
        std::string gridName = it.first;
        std::size_t gridInd = gridsInds[gridName];
        gridsMap[gridName] = &gridCollection.GetGrid(gridInd);
    }

    for(auto& it : gridsRoot.GetPropertyTreeRoot()) {
        std::string gridName = it.first;
        SingleGridParameterExtractor singleGridRoot(it.second);

        YeeGrid3D& grid_i = *gridsMap[gridName];

        SetSingleGridDimensions(grid_i, singleGridRoot);
        SetSingleGridGridArrays(grid_i, singleGridRoot);
        SetSingleGridGridArrayManipulators(grid_i, singleGridRoot);
        SetSingleGridUpddateInstructions(grid_i, singleGridRoot, gridsMap);
        SetSingleGridUpdateSequences(grid_i, singleGridRoot);
        SetSingleGridGridViews(grid_i, singleGridRoot);
    }

    GridCollectionParameterExtractor gridCollectionRoot(node);
    SetAndRunGridCollectionRunSequencs(gridCollection, gridsMap, gridsInds, gridCollectionRoot);
}


void ParamFileTranslator::SetSingleGridDimensions(YeeGrid3D& yee, SingleGridParameterExtractor& singleGridRoot) {
    ParameterExtractor dimensionsExtractor(singleGridRoot.GetSubTreeRootNode("dimensions"));
    auto r0 = dimensionsExtractor.Get3VecRealProperty("r0");
    auto r1 = dimensionsExtractor.Get3VecRealProperty("r1");
    auto nCells = dimensionsExtractor.Get3VecUintProperty("nCells");

    FPNumber dt = dimensionsExtractor.GetRealProperty("dt");

    yee.SetCornerCoordinates(r0, r1);
    yee.SetNumOfCells(nCells);
    yee.SetTimeResolution(dt);
}

void ParamFileTranslator::SetSingleGridGridArrays(YeeGrid3D& yee, SingleGridParameterExtractor& singleGridRoot) {
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
}

void ParamFileTranslator::SetSingleGridGridArrayManipulators(YeeGrid3D& yee,
                                                             SingleGridParameterExtractor& singleGridRoot) {
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
            std::array<std::array<FPNumber, 3>, 3> primitiveVectors_array{primitiveVectors_vec[0],    // a 3-array
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
}

void ParamFileTranslator::SetSingleGridUpddateInstructions(YeeGrid3D& yee,
                                                           SingleGridParameterExtractor& singleGridRoot,
                                                           std::map<std::string, YeeGrid3D*>& gridsMap) {
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

        } else if(updateType == "A+=sumbC::NB" || updateType == "A=sumbC::NB") {
            ParameterExtractor updateParams(std::get<1>(updateTypeandParams));

            auto C_directions_str = updateParams.GetStringArray("C_direction");
            std::vector<int> C_directions;
            for(auto& direction : C_directions_str) {
                C_directions.emplace_back(stringDirectionToIntDirectionMap[direction]);
            }

            void* updateParamsPtr = yee.ConstructParams_A_plusequal_sum_b_C_neighbor(
                gridsMap[updateParams.GetStringProperty("neighborGrid")],
                updateParams.Get3VecUintProperty("A_indStart"),  // indStart
                updateParams.Get3VecUintProperty("A_indEnd"), // indEnd
                updateParams.GetStringProperty("A"),        // A name
                stringDirectionToIntDirectionMap[updateParams.GetStringProperty("A_direction")], // A direcction
                updateParams.GetRealArray("b"),     // b values
                updateParams.GetStringArray("C"),   // C names
                C_directions,
                updateParams.Get3VecUintArray("C_indStart")
                );

            if(updateType == "A+=sumbC::NB") {
                yee.AddUpdateInstruction(updateParams.GetStringProperty("name"),
                        FDInstructionCode::A_plusequal_sum_b_C_neighbor,
                        updateParamsPtr
                        );
            }else if(updateType == "A=sumbC::NB") {
                yee.AddUpdateInstruction(updateParams.GetStringProperty("name"),
                        FDInstructionCode::A_equal_sum_b_C_neighbor,
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
            std::cout << "Unknown update instruction: " << updateType << std::endl;
            assert(false);
        }
    }
}

void ParamFileTranslator::SetSingleGridUpdateSequences(YeeGrid3D& yee,
                                                             SingleGridParameterExtractor& singleGridRoot) {

    auto updateSequences = singleGridRoot.GetUpdateSequences("updateSequences");
    for(auto& updateNameAndSequence : updateSequences) {
        yee.AddInstructionSequence(
                std::get<0>(updateNameAndSequence),
                std::get<1>(updateNameAndSequence)
                );
    };
}

void ParamFileTranslator::SetSingleGridGridViews(YeeGrid3D& yee,
                                                             SingleGridParameterExtractor& singleGridRoot) {
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
}

void ParamFileTranslator::SetAndRunSingleGridRunSequencs(YeeGrid3D& yee,
                                                             SingleGridParameterExtractor& singleGridRoot) {
    auto runSequence = singleGridRoot.GetRunSequence("runSequence");
    for(auto& sequenceNameAndNumOfRuns : runSequence) {
        yee.ApplyInstructions(
                std::get<0>(sequenceNameAndNumOfRuns),
                std::get<1>(sequenceNameAndNumOfRuns),
                std::get<2>(sequenceNameAndNumOfRuns)
                );
    }
}

void ParamFileTranslator::SetAndRunGridCollectionRunSequencs(YeeGridCollection& gridCollection,
                                                             std::map<std::string, YeeGrid3D*>& gridsMap,
                                                             std::map<std::string, std::size_t>& gridsInds,
                                                             GridCollectionParameterExtractor& gridCollectionRoot) {

    auto runSequence = gridCollectionRoot.GetRunSequence("runSequence");
    for(auto& it : runSequence) {
        std::size_t indStart = std::get<0>(it);
        std::size_t indEnd = std::get<1>(it);
        auto& name_sequence_vec = std::get<2>(it);

        std::cout << "Running sequence: " << std::endl;
        std::cout << indStart << " " << indEnd << std::endl;

        // convert names in name_sequence_vec to indices and save them in the new vector index_sequence_vec
        std::vector<std::pair<std::size_t, std::string>> index_sequence_vec;

        for(auto& name_sequence : name_sequence_vec) {
            std::string& gridName = name_sequence.first;
            std::string& sequenceName = name_sequence.second;
            std::cout << gridName << " " << sequenceName << std::endl;

            index_sequence_vec.emplace_back(gridsInds[gridName], sequenceName);
        }

        gridCollection.RunInstructionsPeriodically(indStart, indEnd, index_sequence_vec);
    }
}




