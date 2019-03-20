
#ifndef FDTD_RECTPLANEWAVEVECTORFIELD_H_
#define FDTD_RECTPLANEWAVEVECTORFIELD_H_

#include "PlaneWaveVectorField.h"

class RectPlaneWaveVectorField : public PlaneWaveVectorField {
    public:
    virtual ~RectPlaneWaveVectorField() { };


    void SetRectWidth(FPNumber width);
    void SetRectEdgeWidth(FPNumber width);
    void SetModulationFrequency(FPNumber freq);
    void SetModulationPhase(FPNumber phase);

    virtual FPNumber GetNormalizedWaveform(const FPNumber& t);

    protected:
    FPNumber t_rectWidth = 1.0;
    FPNumber t_edgeWidth = 0.0;
    FPNumber t_modulationFrequency = 0.0;
    FPNumber t_modulationPhase = 0.0;

};

#endif // FDTD_RECTPLANEWAVEVECTORFIELD_H_

