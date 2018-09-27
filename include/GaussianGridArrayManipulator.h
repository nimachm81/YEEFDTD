
#ifndef FDTD_GAUSSIANGRIDARRAYMANIPULATOR_H_
#define FDTD_GAUSSIANGRIDARRAYMANIPULATOR_H_

#include <cassert>
#include <cmath>

#include "GridArrayManipulator.h"

class GaussianGridArrayManipulator : public GridArrayManipulator {
    public:
    virtual ~GaussianGridArrayManipulator() { };
    void SetDirection(const int direction);
    void SetAmplitude(const RealNumber amplitude);
    void SetCenterTime(const RealNumber t_center);
    void SetDecayTime(const RealNumber t_decay);
    void SetModulationFrequency(const RealNumber modulationFrequency);
    void SetModulationPhase(const RealNumber modulationPhase);
    void SetTimeOffsetFraction(const RealNumber offsetFraction);

    void UpdateArray(const RealNumber t);
    RealNumber CalculateTime(const RealNumber dt, const std::size_t timeIndex);

    private:
    int direction;
    RealNumber amplitude;
    RealNumber t_center;
    RealNumber t_decay;
    RealNumber modulationFrequency;
    RealNumber modulationPhase;
    RealNumber timeOffsetFraction = 0.0;

};


#endif // FDTD_GAUSSIANGRIDARRAYMANIPULATOR_H_

