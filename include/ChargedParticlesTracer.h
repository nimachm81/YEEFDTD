#ifndef FDTD_CHARGEDPARTICLESTRACER_H_
#define FDTD_CHARGEDPARTICLESTRACER_H_


#include "YeeGridDataTypes.h"
#include "NumberTypes.h"
#include "ParticlesTracer.h"
#include "VectorField.h"

class ChargedParticlesTracer : public ParticlesTracer {
    public:
    virtual ~ChargedParticlesTracer() { };

    void ReserveMemory(std::size_t numberOfElements);
    void AddParticle(const FPNumber charge,
                     const FPNumber mass,
                     const std::array<FPNumber, 3>& position,
                     const std::array<FPNumber, 3>& velocity,
                     const std::array<FPNumber, 3>& force);
    void AddParticlesEmittedByTheParticleEmitter(FPNumber t);
    void SetMaxChargedPartcileBunchSize(std::size_t bunchSize);

    void SetElectricFieldGrid(YeeGridData3D* eField);
    void SetMagneticFieldGrid(YeeGridData3D* bField);

    void SetElectricFieldGridOrigin(int direction, std::array<FPNumber, 3>& origin);
    void SetMagneticFieldGridOrigin(int direction, std::array<FPNumber, 3>& origin);

    void SetAnalyticElectricField(VectorField* eField);
    void SetAnalyticMagneticField(VectorField* bField);

    void UpdateElectricForce(int direction // x, y or z component. Directiion of the electric field
                             );
    void UpdateMagneticForce(int direction // 0=x, 1=y or 2=z component. Directiion of the magnetic field
                             );

    void UpdateParticlesCurrents();  // updates currents using particles velocities

    virtual void PointToScalarData(std::vector<FPNumber>*& values, std::string dataName, int direction);
    virtual void PointToVectorData(std::vector<std::array<FPNumber, 3>>*& values,
                                             std::string dataName);
    void UpdateGAMValues(const FPNumber t);

    private:
    std::vector<FPNumber> charges;
    std::array<std::vector<FPNumber>, 3> currentComponents;

    YeeGridData3D* electricField = nullptr;       // E
    std::array<std::array<FPNumber, 3>, 3> electricFieldConponentsOrigin;
    YeeGridData3D* magneticField = nullptr;       // B
    std::array<std::array<FPNumber, 3>, 3> magneticFieldConponentsOrigin;

    VectorField* analyticElectricField = nullptr;       // E
    VectorField* analyticMagneticField = nullptr;       // B

    std::size_t maxChargedParticleBunchSize = 1;  // tries to bunch several particles into 1 with this number of particles
    std::vector<FPNumber> chargeParticleEmissionNumberTracker;

};

#endif // FDTD_CHARGEDPARTICLESTRACER_H_
