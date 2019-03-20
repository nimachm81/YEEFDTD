
#ifndef FDTD_GAUSSIANPLANEWAVEVECTORFIELD_H_
#define FDTD_GAUSSIANPLANEWAVEVECTORFIELD_H_

#include "PlaneWaveVectorField.h"

class GaussianPlaneWaveVectorField : public PlaneWaveVectorField {
    public:
    virtual ~GaussianPlaneWaveVectorField() { };


    void SetTimeDecayRate(FPNumber rate);
    void SetModulationFrequency(FPNumber freq);
    void SetModulationPhase(FPNumber phase);

    virtual FPNumber GetNormalizedWaveform(const FPNumber& t);

    protected:
    FPNumber t_decayRate = 1.0;
    FPNumber t_modulationFrequency = 0.0;
    FPNumber t_modulationPhase = 0.0;

};

#endif // FDTD_GAUSSIANPLANEWAVEVECTORFIELD_H_



