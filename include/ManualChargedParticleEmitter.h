#ifndef FDTD_MANUALCHARGEDPARTICLEEMITTER_H_
#define FDTD_MANUALCHARGEDPARTICLEEMITTER_H_

#include "NumberTypes.h"
#include "ChargedParticleEmitter.h"

// Emits particles with initial positions and velocities at specified times
class ManualChargedParticleEmitter : public ChargedParticleEmitter {
    public:
    virtual ~ManualChargedParticleEmitter() { };

    void SetEmissionTimes(const std::vector<FPNumber>& times);
    void SetNumberOfParticlesToEmitAtEachSpecifiedTime(const std::vector<FPNumber>& number);
    void SetInitialPositions(const std::vector<std::array<FPNumber, 3>>& positions);
    void SetInitialVelocities(const std::vector<std::array<FPNumber, 3>>& velocities);

    // overrrides ChargedParticleEmitter
    const std::vector<FPNumber>& GetEmissionNumber(FPNumber t);       // get number of emissitted particles per surface area at time t


    private:
    std::vector<FPNumber> emissionTimes;
    std::vector<bool> isEmitted;
    std::vector<FPNumber> numberOfPatriclesToEmit;
    std::vector<std::array<FPNumber, 3>> initialPositions;     // point on the surface of the object
    std::vector<std::array<FPNumber, 3>> initialVelocities;     // initial velocities of the emitted particles

};

#endif // FDTD_MANUALCHARGEDPARTICLEEMITTER_H_

