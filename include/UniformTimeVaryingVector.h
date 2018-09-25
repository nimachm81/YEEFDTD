
#ifndef FDTD_UNIFORMTIMEVARYINGCURRENT_H_
#define FDTD_UNIFORMTIMEVARYINGCURRENT_H_

#include "GridArrayManipulator.h"

class UniformSingleDirectionGaussianTimeVaryingVector : public GridArrayManipulator {
    public:
    void SetDirection(const int direction);
    void SetAmplitude(const RealNumber amplitude);
    void SetCenterTime(const RealNumber t_center);
    void SetDecayTime(const RealNumber t_decay);
    void SetModulationFrequency(const RealNumber modulationFrequency);
    void SetModulationPhase(const RealNumber modulationPhase);


    virtual void UpdateArray(YeeGridData3D& gridData, RealNumber t);

    private:
    int direction;
    RealNumber amplitude;
    RealNumber t_center;
    RealNumber t_decay;
    RealNumber modulationFrequency;
    RealNumber modulationPhase;

}


#endif // FDTD_UNIFORMTIMEVARYINGCURRENT_H_

