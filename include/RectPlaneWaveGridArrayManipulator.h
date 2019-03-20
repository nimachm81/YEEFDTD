
#ifndef FDTD_RECTPLANEWAVEGRIDARRAYMANIPULATOR_H_
#define FDTD_RECTPLANEWAVEGRIDARRAYMANIPULATOR_H_

#include "PlaneWaveGridArrayManipulator.h"

class RectPlaneWaveGridArrayManipulator : public PlaneWaveGridArrayManipulator {
    public:
    RectPlaneWaveGridArrayManipulator();
    virtual ~RectPlaneWaveGridArrayManipulator() { };

    void SetRectWidth(FPNumber width);
    void SetRectEdgeWidth(FPNumber width);
    void SetModulationFrequency(FPNumber freq);
    void SetModulationPhase(FPNumber phase);

    virtual FPNumber GetNormalizedTemporalWaveform(FPNumber t);

    protected:
    FPNumber t_rectWidth = 1.0;
    FPNumber t_edgeWidth = 0.0;
    FPNumber t_modulationFrequency = 0.0;
    FPNumber t_modulationPhase = 0.0;

};


#endif // FDTD_RECTPLANEWAVEGRIDARRAYMANIPULATOR_H_


