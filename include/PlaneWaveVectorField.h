
#ifndef FDTD_PLANEWAVEVECTORFIELD_H_
#define FDTD_PLANEWAVEVECTORFIELD_H_

#include "VectorField.h"

class PlaneWaveVectorField : public VectorField {
    public:
    virtual ~PlaneWaveVectorField() { };

    void SetPropagationVelocity(const FPNumber v);
    void SetPropagationDirection(const std::array<FPNumber, 3> direction);
    void SetAmplitude(std::array<FPNumber, 3>& amp);
    void SetCenterTime(FPNumber t_c);

    virtual FPNumber GetNormalizedWaveform(const FPNumber& t) = 0;

    virtual void GetFieldValueAtPoint(FPNumber time,
                          std::array<FPNumber, 3>& position,
                          std::array<FPNumber, 3>& fieldValue
                          );
    virtual void GetFieldValuesAtPoints(FPNumber time,
                          std::vector<std::array<FPNumber, 3>>& positions,
                          std::vector<std::array<FPNumber, 3>>& fieldValues
                          );

    protected:
    std::array<FPNumber, 3> propagationDirection;
    FPNumber velocity;
    std::array<FPNumber, 3> amplitude;
    FPNumber t_center = 0.0;

};

#endif // FDTD_PLANEWAVEVECTORFIELD_H_



