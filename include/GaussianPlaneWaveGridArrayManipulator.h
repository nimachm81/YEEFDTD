
#ifndef FDTD_GAUSSIANPLANEWAVEGRIDARRAYMANIPULATOR_H_
#define FDTD_GAUSSIANPLANEWAVEGRIDARRAYMANIPULATOR_H_

#include "PlaneWaveGridArrayManipulator.h"

class GaussianPlaneWaveGridArrayManipulator : public PlaneWaveGridArrayManipulator {
    public:
    virtual ~GaussianPlaneWaveGridArrayManipulator() { };

    void SetTimeDecayRate(FPNumber rate);
    void SetModulationFrequency(FPNumber freq);
    void SetModulationPhase(FPNumber phase);

    virtual FPNumber GetNormalizedTemporalWaveform(FPNumber t);

    protected:
    FPNumber t_decayRate = 1.0;
    FPNumber t_modulationFrequency = 0.0;
    FPNumber t_modulationPhase = 0.0;

};


#endif // FDTD_GAUSSIANPLANEWAVEGRIDARRAYMANIPULATOR_H_


