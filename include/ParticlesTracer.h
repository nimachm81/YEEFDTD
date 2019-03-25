
#ifndef FDTD_PARTICLESTRACER_H_
#define FDTD_PARTICLESTRACER_H_

#include <vector>
#include <array>
#include <cassert>

#include "NumberTypes.h"
#include "YeeGridDataTypes.h"
#include "DiscretePointsGAMDataUpdater.h"
#include "ParticleEmitter.h"
#include "Geometry.h"

class ParticlesTracer : public DiscretePointsGAMDataUpdater {
    public:
    virtual ~ParticlesTracer() { };

    void ReserveMemory(std::size_t numberOfElements);
    void AddParticle(const FPNumber mass,
                     const std::array<FPNumber, 3>& position,
                     const std::array<FPNumber, 3>& velocity,
                     const std::array<FPNumber, 3>& force);

    void SetGridSpacing(std::array<FPNumber, 3>& dr);
    void SetScatteringRateFieldGrid(YeeGridData3D* srField);
    void SetScatteringRateFieldGridOrigin(int direction, std::array<FPNumber, 3>& origin);

    virtual void PointToDataPositions(std::vector<std::array<FPNumber, 3>>*& positions);
    virtual void PointToScalarData(std::vector<FPNumber>*& values,
                                   std::string dataName,
                                   int direction);

    virtual void PointToVectorData(std::vector<std::array<FPNumber, 3>>*& values,
                                 std::string dataName);

    void SetParticleEmitter(ParticleEmitter* emitter);
    void AddParticlesEmittedByTheParticleEmitter(FPNumber t);   // adds the particles emitted at time t to the collection of particles

    void SetConstrainingGeometry(Geometry* geometry, bool keepInside);
    void UpdateParticlesConstraintStatus();

    void UpdateParticlesPositions(const FPNumber dt, bool applyConstraints = true);
    void UpdateParticlesMomentumsAndVelocities(const FPNumber dt);
    void UpdateScatteringForce(int direction);
    void UpdateTime(const FPNumber newTime);

    void UpdateParticlesMomentumVelocityPosition(const FPNumber newTime);

    void SetForce(std::size_t index, std::array<FPNumber, 3>& force);
    void ResetForces();

    protected:
    FPNumber time = 0.0;
    std::vector<FPNumber> masses;
    std::vector<std::array<FPNumber, 3>> positions;
    std::vector<std::array<FPNumber, 3>> velocities;
    std::vector<std::array<FPNumber, 3>> momentums;
    std::vector<std::array<FPNumber, 3>> forces;

    YeeGridData3D* scatteringRateField = nullptr;       // gamma as a field. dp_x = -gamma_x*px*dt where gamma_x is scattering rate at the position of the particle
    std::array<std::array<FPNumber, 3>, 3> scatteringRateFieldConponentsOrigin;
    std::array<FPNumber, 3> gridSpacing;    // grid spacing for the fields

    ParticleEmitter* particleEmitter = nullptr;

    Geometry* constrainingGeometry = nullptr;     // if initalized the particles are constrained inside/outside the gromtry
    bool keepInsideGeometry = true;     // true: keep particles inside, false: keep paarticles outside
    std::vector<bool> arePointsInside;     // true: constraint satisfied       false : not satisfied
};

#endif  // FDTD_PARTICLESTRACER_H_
