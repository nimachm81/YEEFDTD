
#ifndef FDTD_PARTICLESTRACER_H_
#define FDTD_PARTICLESTRACER_H_

#include <vector>
#include <array>

#include "NumberTypes.h"

class ParticlesTracer {
    public:
    void AddParticle(const FPNumber mass,
                     const std::array<FPNumber, 3>& position,
                     const std::array<FPNumber, 3>& velocity,
                     const std::array<FPNumber, 3>& force);
    void UpdateParticlesPositions(const FPNumber dt);
    void UpdateParticlesMomentumsAndVelocities(const FPNumber dt);
    inline void SetForce(std::size_t index, std::array<FPNumber, 3>& force);

    protected:
    std::vector<FPNumber> masses;
    std::vector<std::array<FPNumber, 3>> positions;
    std::vector<std::array<FPNumber, 3>> velocities;
    std::vector<std::array<FPNumber, 3>> momentums;
    std::vector<std::array<FPNumber, 3>> forces;
};

#endif  // FDTD_PARTICLESTRACER_H_
