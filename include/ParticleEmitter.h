#ifndef FDTD_PARTICLEEMITTER_H_
#define FDTD_PARTICLEEMITTER_H_

#include <vector>
#include <array>
#include <cstddef>
#include <unordered_map>

#include "NumberTypes.h"

class ParticleEmitter {
    public:
    virtual ~ParticleEmitter() { };

    virtual const std::vector<FPNumber>& GetEmissionNumber(FPNumber t) = 0; // number of emitted particles at time t
    virtual const std::vector<std::array<FPNumber, 3>>& GetEmissionPoints() = 0;
    virtual const std::vector<std::array<FPNumber, 3>>& GetEmissionVelocities() = 0;

    virtual std::vector<std::vector<std::array<FPNumber, 3>>>* GetEmissionSubPoints() = 0;

    virtual const std::unordered_map<std::string, FPNumber>& GetParticleParameters() = 0;

    void UpdateTime(FPNumber t) {
        time = t;
    }

    protected:
    FPNumber time = 0.0;

};

#endif // FDTD_PARTICLEEMITTER_H_

