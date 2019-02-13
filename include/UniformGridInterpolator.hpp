#ifndef FDTD_UNIFORMGRIDINTERPOLATOR_H_
#define FDTD_UNIFORMGRIDINTERPOLATOR_H_

#include <vector>
#include <array>
#include <cstddef>

#include "NumberTypes.h"
#include "MultiDimArray.hpp"

class UniformGridInterpolator {
    public:
    static void InterpolateGridOnPoints(const NumberArray3D<FPNumber> gridArray,
                                        const std::array<FPNumber, 3> r0,       // gridArray origin
                                        const std::array<FPNumber, 3> dr,       // grid spacing
                                        const  std::vector<std::array<FPNumber, 3>> points, // interpolate gridArray on these points
                                        std::vector<FPNumber>& interpolatedValues  // will be updated by this function
                                        ) {
        std::size_t numOfPoints = points.size();
        if(interpolatedValues.size() != numOfPoints) {
            interpolatedValues.resize(numOfPoints);
        }

        for(std::size_t i = 0; i < numOfPoints; ++i) {
            const std::array<FPNumber, 3>& position = points[i];

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

            FPNumber valInterp = 0.0;     // interpolated value

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
                        indx[2] = lowerIndx[2] + i_x;
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
                        valInterp += gridArray[indx_sizet]*vol_ratio[0]*vol_ratio[1]*vol_ratio[2];
                    }
                }
            }

            interpolatedValues[i] = valInterp;
        }
    };


};

#endif // FDTD_UNIFORMGRIDINTERPOLATOR_H_

