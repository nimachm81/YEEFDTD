

#include "WedgeGridArrayManipulator.h"

void WedgeGridArrayManipulator::SetTipAngle(FPNumber angle) {
    wedgeAngle = angle;
}

void WedgeGridArrayManipulator::SetTipAngleInDegrees(FPNumber angle_degree) {
    wedgeAngle = angle_degree * M_PI / 180.0;
}

void WedgeGridArrayManipulator::SetTipRadius(FPNumber radius) {
    tipRadius = radius;
}

void WedgeGridArrayManipulator::SetWedgeHeight(FPNumber height) {
    wedgeHeight = height;
}

void WedgeGridArrayManipulator::SetTipPosition(std::array<FPNumber, 3> pos) {
    tipPosition = pos;
}

void WedgeGridArrayManipulator::SetInsideValue(FPNumber value) {
    valueInside = value;
}

void WedgeGridArrayManipulator::SetOutsideValue(FPNumber value) {
    valueOutside = value;
}

std::array<FPNumber, 3> WedgeGridArrayManipulator::GetTipPositionGivenRoundedTipPosition(
                                FPNumber wedgeAngle,
                                FPNumber tipRadius,
                                std::array<FPNumber, 3> roundedTipTopPosition
                                ) {
    FPNumber tip_to_flattop = tipRadius * std::cos(wedgeAngle/2.0) / std::tan(wedgeAngle/2.0);

    FPNumber tip_to_circleCenter = tip_to_flattop + tipRadius*std::sin(wedgeAngle/2.0);

    FPNumber tip_to_roundedtop = tip_to_circleCenter - tipRadius;

    return std::array<FPNumber, 3>{roundedTipTopPosition[0],
                                   roundedTipTopPosition[1] + tip_to_roundedtop,
                                   roundedTipTopPosition[2]};
}

FPNumber WedgeGridArrayManipulator::CalculateTime(const FPNumber dt, const std::size_t timeIndex) {
    return ((FPNumber)timeIndex + timeOffsetFraction) * dt;
}

void WedgeGridArrayManipulator::UpdateArray(const FPNumber t, GAManipulatorInstructionCode instruction) {
    FPNumber*** arrayData = gridArray.GetArrayData();
    std::array<std::size_t, 3>& arrayShape = gridArray.GetShape();
    std::array<std::size_t, 3>& arrayIndStart = gridArray.GetIndStart();

    std::size_t n0 = arrayShape[0];
    std::size_t n1 = arrayShape[1];
    std::size_t n2 = arrayShape[2];
    std::size_t ind0 = arrayIndStart[0];
    std::size_t ind1 = arrayIndStart[1];
    std::size_t ind2 = arrayIndStart[2];

    FPNumber x0 = r0[0];
    FPNumber y0 = r0[1];
    FPNumber z0 = r0[2];

    FPNumber dx = dr[0];
    FPNumber dy = dr[1];
    FPNumber dz = dr[2];

    FPNumber x, y, z;

    FPNumber slope_p = std::tan(M_PI/2.0 + wedgeAngle/2.0);
    FPNumber slope_m = std::tan(M_PI/2.0 - wedgeAngle/2.0);

    FPNumber x_tip = tipPosition[0];
    FPNumber y_tip = tipPosition[1];
    FPNumber z_tip = tipPosition[2];

    FPNumber base_halfWidth = wedgeHeight * std::tan(wedgeAngle/2.0);

    FPNumber tip_to_top = tipRadius * std::cos(wedgeAngle/2.0) / std::tan(wedgeAngle/2.0);

    FPNumber tip_to_circleCenter = tip_to_top + tipRadius*std::sin(wedgeAngle/2.0);

    if(instruction == GAManipulatorInstructionCode::Equal) {
        if(uniformAxis == 0 && wedgeDirection == 1) {
            for(std::size_t i0 = 0; i0 < n0; ++i0) {
                //x = x0 + (FPNumber)i0*dx;
                for(std::size_t i1 = 0; i1 < n1; ++i1) {
                    y = y0 + (FPNumber)i1*dy - y_tip;   // y coordinate relative to tip
                    for(std::size_t i2 = 0; i2 < n2; ++i2) {
                        if(y > 0.0 || y < -wedgeHeight) {
                            arrayData[ind0 + i0][ind1 + i1][ind2 + i2] = valueOutside;
                        } else if(y <= -tip_to_top) {
                            z = z0 + (FPNumber)i2*dz - z_tip;   // z coordinate relative to tip
                            if(std::abs(z) > base_halfWidth) {
                                arrayData[ind0 + i0][ind1 + i1][ind2 + i2] = valueOutside;
                            } else {

                                if(y <= z*slope_p && y <= z*slope_m) {
                                    arrayData[ind0 + i0][ind1 + i1][ind2 + i2] = valueInside;
                                } else {
                                    arrayData[ind0 + i0][ind1 + i1][ind2 + i2] = valueOutside;
                                }
                            }
                        } else {    // check whether inside cap circle
                            z = z0 + (FPNumber)i2*dz - z_tip;   // z coordinate relative to tip
                            if((y + tip_to_circleCenter)*(y + tip_to_circleCenter) + z*z <= tipRadius*tipRadius) {
                                arrayData[ind0 + i0][ind1 + i1][ind2 + i2] = valueInside;
                            } else {
                                arrayData[ind0 + i0][ind1 + i1][ind2 + i2] = valueOutside;
                            }
                        }
                    }
                }
            }
        } else {
            std::cout << "error: not implemented." << std::endl;
            assert(false);
        }
    } else if(instruction == GAManipulatorInstructionCode::MultiplyEqual) {
        if(uniformAxis == 0 && wedgeDirection == 1) {
            for(std::size_t i0 = 0; i0 < n0; ++i0) {
                //x = x0 + (FPNumber)i0*dx;
                for(std::size_t i1 = 0; i1 < n1; ++i1) {
                    y = y0 + (FPNumber)i1*dy - y_tip;   // y coordinate relative to tip
                    for(std::size_t i2 = 0; i2 < n2; ++i2) {
                        if(y > 0.0 || y < -wedgeHeight) {
                            arrayData[ind0 + i0][ind1 + i1][ind2 + i2] *= valueOutside;
                        } else if(y <= -tip_to_top) {
                            z = z0 + (FPNumber)i2*dz - z_tip;   // z coordinate relative to tip
                            if(std::abs(z) > base_halfWidth) {
                                arrayData[ind0 + i0][ind1 + i1][ind2 + i2] *= valueOutside;
                            } else {

                                if(y <= z*slope_p && y <= z*slope_m) {
                                    arrayData[ind0 + i0][ind1 + i1][ind2 + i2] *= valueInside;
                                } else {
                                    arrayData[ind0 + i0][ind1 + i1][ind2 + i2] *= valueOutside;
                                }
                            }
                        } else {    // check whether inside cap circle
                            z = z0 + (FPNumber)i2*dz - z_tip;   // z coordinate relative to tip
                            if((y + tip_to_circleCenter)*(y + tip_to_circleCenter) + z*z <= tipRadius*tipRadius) {
                                arrayData[ind0 + i0][ind1 + i1][ind2 + i2] *= valueInside;
                            } else {
                                arrayData[ind0 + i0][ind1 + i1][ind2 + i2] *= valueOutside;
                            }
                        }
                    }
                }
            }
        } else {
            std::cout << "error: not implemented." << std::endl;
            assert(false);
        }
    } else {
            std::cout << "error: GAManipulatorInstructionCode operation not implemented." << std::endl;
            assert(false);
    }
}


