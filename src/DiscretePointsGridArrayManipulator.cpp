
#include <cstddef>

#include "DiscretePointsGridArrayManipulator.h"


void DiscretePointsGridArrayManipulator::AddDataUpdater(DiscretePointsGAMDataUpdater* updater,
                                                        std::string dataName,
                                                        int direction
                                                        ) {
    dataUpdater = updater;
    updater->AttachDataToGAMPositions(positions);
    updater->AttachDataToGAMValues(values, dataName, direction);
}


void DiscretePointsGridArrayManipulator::UpdateArray(const FPNumber t, GAManipulatorInstructionCode instruction) {
    dataUpdater->UpdateGAMValues(t);
    if(instruction == GAManipulatorInstructionCode::Equal) {
        gridArray = (FPNumber)0.0;

        std::size_t numOfPoints = positions->size();
        assert(values->size() == numOfPoints);

        if(interpolationType == 0) {
            for(std::size_t i = 0; i < numOfPoints; ++i) {
                std::array<FPNumber, 3>& position = (*positions)[i];

                std::array<std::intmax_t, 3> lowerIndx
                                   {static_cast<std::intmax_t>(std::floor((position[0] - r0[0]) / dr[0])),
                                    static_cast<std::intmax_t>(std::floor((position[1] - r0[1]) / dr[1])),
                                    static_cast<std::intmax_t>(std::floor((position[2] - r0[2]) / dr[2]))};

                std::array<std::intmax_t, 3> closestIndx;

                for(int j = 0; j < 3; ++j) {
                    if((position[j] - r0[j]) - lowerIndx[j]*dr[j] <= dr[j]/(FPNumber)2.0) {
                        closestIndx[j] = lowerIndx[j];
                    } else {
                        closestIndx[j] = lowerIndx[j] + 1;
                    }
                }

                // assert inside bounds
                const std::array<std::size_t, 3>& gridArrayShape = gridArray.GetShape();
                if(closestIndx[0] >= 0 && closestIndx[0] < gridArrayShape[0] &&
                   closestIndx[1] >= 0 && closestIndx[1] < gridArrayShape[1] &&
                   closestIndx[2] >= 0 && closestIndx[2] < gridArrayShape[2]) {
                    const std::array<std::size_t, 3> indx_sizet{(std::size_t)closestIndx[0],
                                                                (std::size_t)closestIndx[1],
                                                                (std::size_t)closestIndx[2]};
                    gridArray[indx_sizet] += (*values)[i];
                }
            }
        } else if(interpolationType == 1) {
            for(std::size_t i = 0; i < numOfPoints; ++i) {
                std::array<FPNumber, 3>& position = (*positions)[i];

                std::array<std::intmax_t, 3> lowerIndx
                                   {static_cast<std::intmax_t>(std::floor((position[0] - r0[0]) / dr[0])),
                                    static_cast<std::intmax_t>(std::floor((position[1] - r0[1]) / dr[1])),
                                    static_cast<std::intmax_t>(std::floor((position[2] - r0[2]) / dr[2]))};

                std::array<FPNumber, 3> alpha;
                for(int j = 0; j < 3; ++j) {
                    alpha[j] = ((position[j] - r0[j]) - lowerIndx[j]*dr[j])/dr[j];
                }

                // assert inside bounds
                const std::array<std::size_t, 3>& gridArrayShape = gridArray.GetShape();

                std::array<std::intmax_t, 3> indx;
                std::array<FPNumber, 3> vol_ratio;
                FPNumber value = (*values)[i];
//                std::cout << "-------------------------------------" << std::endl;
//                std::cout << "r0:" << r0[0] << ", " << r0[1] << ", " << r0[2] << std::endl;
//                std::cout << "dr:" << dr[0] << ", " << dr[1] << ", " << dr[2] << std::endl;
//                std::cout << "r:" << position[0] << ", " << position[1] << ", " << position[2] << ", v: " << value << std::endl;
//                std::cout << "li:" << lowerIndx[0] << ", " << lowerIndx[1] << ", " << lowerIndx[2] << std::endl;
//                std::cout << "shape:" << gridArrayShape[0] << ", " << gridArrayShape[1] << ", " << gridArrayShape[2] << std::endl;
                for(int i_x = 0; i_x < 2; ++i_x) {
                    indx[0] = lowerIndx[0] + i_x;
                    if(indx[0] < 0 || indx[0] >= gridArrayShape[0]) {
                        continue;
                    }
                    if(i_x == 0) {
                        vol_ratio[0] = 1.0 - alpha[0];
                    } else {
                        vol_ratio[0] = alpha[0];
                    }
                    for(int i_y = 0; i_y < 2; ++i_y) {
                        indx[1] = lowerIndx[1] + i_y;
                        if(indx[1] < 0 || indx[1] >= gridArrayShape[1]) {
                            continue;
                        }
                        if(i_y == 0) {
                            vol_ratio[1] = 1.0 - alpha[1];
                        } else {
                            vol_ratio[1] = alpha[1];
                        }
                        for(int i_z = 0; i_z < 2; ++i_z) {
                            indx[2] = lowerIndx[2] + i_z;
                            if(indx[2] < 0 || indx[2] >= gridArrayShape[2]) {
                                continue;
                            }
                            if(i_z == 0) {
                                vol_ratio[2] = 1.0 - alpha[2];
                            } else {
                                vol_ratio[2] = alpha[2];
                            }

                            const std::array<std::size_t, 3> indx_sizet{(std::size_t)indx[0],
                                                                        (std::size_t)indx[1],
                                                                        (std::size_t)indx[2]};
                            gridArray[indx_sizet] += value*vol_ratio[0]*vol_ratio[1]*vol_ratio[2];
                            //std::cout << indx[0] << ", " << indx[1] << ", " << indx[2] << ", j:" << gridArray[indx_sizet] << std::endl;

                        }
                    }
                }
            }

        }
    } else {
        std::cout << "error: not implemented" << std::endl;
        assert(false);
    }
}

FPNumber DiscretePointsGridArrayManipulator::CalculateTime(const FPNumber dt, const std::size_t timeIndex) {
    return ((FPNumber)timeIndex + timeOffsetFraction) * dt;
}

