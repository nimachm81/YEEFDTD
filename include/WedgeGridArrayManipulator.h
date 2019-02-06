#ifndef FDTD_WEDGEGRIDARRAYMANIPULATOR_H_
#define FDTD_WEDGEGRIDARRAYMANIPULATOR_H_

#include <cassert>
#include <cmath>

#include "GridArrayManipulator.h"

class WedgeGridArrayManipulator : public GridArrayManipulator {
    public:
    virtual ~WedgeGridArrayManipulator() { };

    void SetTipAngle(FPNumber angle);
    void SetTipAngleInDegrees(FPNumber angle_degree);
    void SetTipRadius(FPNumber radius);
    void SetWedgeHeight(FPNumber height);
    void SetTipPosition(std::array<FPNumber, 3> pos);
    void SetInsideValue(FPNumber value);
    void SetOutsideValue(FPNumber value);

    static
    std::array<FPNumber, 3> GetTipPositionGivenRoundedTipPosition(              // given wedgeAngle, tipRadius and
                                FPNumber wedgeAngle,                            // the top position of the rounded tip
                                FPNumber tipRadius,                             // returns the position of the unrounded tip
                                std::array<FPNumber, 3> rendedTipTopPosition);

    void UpdateArray(const FPNumber t, GAManipulatorInstructionCode instruction);
    FPNumber CalculateTime(const FPNumber dt, const std::size_t timeIndex);

    private:
    int uniformAxis = 0;    // along this axis there is no varation (0:x, 1:y, 2:z)
    int wedgeDirection = 1;      // the wedge is pointing along this direction (0:x, 1:y, 2:z)
    FPNumber wedgeAngle = M_PI/4.0;         // wedge angle in radians
    FPNumber tipRadius = 0.0;   // for tipRadius > 0 the tip is rounded
    std::array<FPNumber, 3> tipPosition;    // position of the unrounded tip
    FPNumber wedgeHeight = 1.0;         // distance from unrounded tip to base of the wedge

    FPNumber valueInside = 1.0;
    FPNumber valueOutside = 0.0;
    FPNumber timeOffsetFraction = 0.0;
};

#endif // FDTD_WEDGEGRIDARRAYMANIPULATOR_H_
