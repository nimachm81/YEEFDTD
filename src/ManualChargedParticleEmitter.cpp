

#include "ManualChargedParticleEmitter.h"


void ManualChargedParticleEmitter::SetEmissionTimes(const std::vector<FPNumber>& times) {
    emissionTimes = times;
    isEmitted.resize(emissionTimes.size());
    for(std::size_t i = 0; i < isEmitted.size(); ++i) {
        isEmitted[i] = false;
    }
}

void ManualChargedParticleEmitter::SetNumberOfParticlesToEmitAtEachSpecifiedTime(const std::vector<FPNumber>& number) {
    numberOfPatriclesToEmit = number;
}


void ManualChargedParticleEmitter::SetInitialPositions(const std::vector<std::array<FPNumber, 3>>& positions) {
    initialPositions = positions;
}


void ManualChargedParticleEmitter::SetInitialVelocities(const std::vector<std::array<FPNumber, 3>>& velocities) {
    initialVelocities = velocities;
}

const std::vector<FPNumber>& ManualChargedParticleEmitter::GetEmissionNumber(FPNumber t) {
    std::size_t numOfParticles = emissionTimes.size();
    emissionPoints.clear();
    emissionVelocities.clear();
    emissionNumbers.clear();
    for(std::size_t i = 0; i < numOfParticles; ++i) {
        if(!isEmitted[i] && t >= emissionTimes[i]) {
            emissionPoints.push_back(initialPositions[i]);
            emissionVelocities.push_back(initialVelocities[i]);
            emissionNumbers.push_back(numberOfPatriclesToEmit[i]);
            isEmitted[i] = true;
        }
    }

    return emissionNumbers;
}

