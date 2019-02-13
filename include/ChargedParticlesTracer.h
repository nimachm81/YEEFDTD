#ifndef FDTD_CHARGEDPARTICLESTRACER_H_
#define FDTD_CHARGEDPARTICLESTRACER_H_


#include "YeeGridDataTypes.h"
#include "NumberTypes.h"
#include "ParticlesTracer.h"

class ChargedParticlesTracer : public ParticlesTracer {
    public:
    virtual ~ChargedParticlesTracer() { };

    void SetGridSpacing(std::array<FPNumber, 3>& dr);
    void AddParticle(const FPNumber charge,
                     const FPNumber mass,
                     const std::array<FPNumber, 3>& position,
                     const std::array<FPNumber, 3>& velocity,
                     const std::array<FPNumber, 3>& force);
    void AddParticlesEmittedByTheParticleEmitter(FPNumber t,  bool bunchParticlesAsOne = true);

    void SetElectricFieldGrid(YeeGridData3D* eField);
    void SetMagneticFieldGrid(YeeGridData3D* bField);

    void SetElectricFieldGridOrigin(int direction, std::array<FPNumber, 3>& origin);
    void SetMagneticFieldGridOrigin(int direction, std::array<FPNumber, 3>& origin);

    void UpdateElectricForce(int direction // x, y or z component. Directiion of the electric field
                             );
    void UpdateMagneticForce(int direction // 0=x, 1=y or 2=z component. Directiion of the magnetic field
                             );

    void UpdateParticlesCurrents();  // updates currents using particles velocities

    void AttachDataToGAMPositions(std::vector<std::array<FPNumber, 3>>*& positions);
    void AttachDataToGAMValues(std::vector<FPNumber>*& values, std::string dataName, int direction);
    void UpdateGAMValues(const FPNumber t);

    private:
    std::vector<FPNumber> charges;
    std::array<std::vector<FPNumber>, 3> currentComponents;
    YeeGridData3D* electricField;       // E
    std::array<std::array<FPNumber, 3>, 3> electricFieldConponentsOrigin;
    YeeGridData3D* magneticField;       // B
    std::array<std::array<FPNumber, 3>, 3> magneticFieldConponentsOrigin;
    std::array<FPNumber, 3> gridSpacing;

};

#endif // FDTD_CHARGEDPARTICLESTRACER_H_
