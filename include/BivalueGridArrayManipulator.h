#ifndef FDTD_BIVALUEGRIDARRAYMANIPULATOR_H_
#define FDTD_BIVALUEGRIDARRAYMANIPULATOR_H_

#include <cassert>
#include <cmath>

#include "GridArrayManipulator.h"
#include "Geometry.h"

class BivalueGridArrayManipulator : public GridArrayManipulator {
    public:
    virtual ~BivalueGridArrayManipulator() { };

    void SetupConditionArray();
    void SetGeometry(std::shared_ptr<Geometry> geometryPtr);
    void SetInsideValue(FPNumber value);
    void SetOutsideValue(FPNumber value);

    void UpdateArray(const FPNumber t, GAManipulatorInstructionCode instruction);
    FPNumber CalculateTime(const FPNumber dt, const std::size_t timeIndex);

    private:
    std::shared_ptr<Geometry> geometry;

    FPNumber valueInside = 1.0;
    FPNumber valueOutside = 0.0;
    FPNumber timeOffsetFraction = 0.0;

    NumberArray3D<std::int8_t> conditionArray;     // 0: array point inside 1: array point outside
};

#endif // FDTD_BIVALUEGRIDARRAYMANIPULATOR_H_
