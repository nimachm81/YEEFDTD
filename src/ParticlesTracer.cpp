

#include "UniformGridInterpolator.hpp"
#include "ParticlesTracer.h"

void ParticlesTracer::ReserveMemory(std::size_t numberOfElements) {
    masses.reserve(numberOfElements);
    positions.reserve(numberOfElements);
    velocities.reserve(numberOfElements);
    momentums.reserve(numberOfElements);
    forces.reserve(numberOfElements);
}

void ParticlesTracer::AddParticle(const FPNumber mass,
                 const std::array<FPNumber, 3>& position,
                 const std::array<FPNumber, 3>& velocity,
                 const std::array<FPNumber, 3>& force) {
    FPNumber velocitySquared = velocity[0]*velocity[0] + velocity[1]*velocity[1] + velocity[2]*velocity[2];
    FPNumber gamma = 1.0/std::sqrt(1 - velocitySquared);
    std::array<FPNumber, 3> momentum{gamma*mass*velocity[0],
                                     gamma*mass*velocity[1],
                                     gamma*mass*velocity[2]};
    masses.push_back(mass);
    positions.push_back(position);
    velocities.push_back(velocity);
    momentums.push_back(momentum);
    forces.push_back(force);
}

void ParticlesTracer::SetScatteringRateFieldGrid(YeeGridData3D* srField) {
    scatteringRateField = srField;
}

void ParticlesTracer::SetScatteringRateFieldGridOrigin(int direction, std::array<FPNumber, 3>& origin) {
    scatteringRateFieldConponentsOrigin[direction] = origin;
}


void ParticlesTracer::SetGridSpacing(std::array<FPNumber, 3>& dr) {
    gridSpacing = dr;
}

void ParticlesTracer::SetParticleEmitter(ParticleEmitter* emitter) {
    particleEmitter = emitter;
}


void ParticlesTracer::AddParticlesEmittedByTheParticleEmitter(FPNumber t) {
    if(particleEmitter == nullptr) {
        return;
    }

    std::unordered_map<std::string, FPNumber> particleParams = particleEmitter->GetParticleParameters();
    const FPNumber mass = particleParams["mass"];
    const std::vector<FPNumber>& numOfEmittedParticles = particleEmitter->GetEmissionNumber(t);
    const std::vector<std::array<FPNumber, 3>>& emissionPoints = particleEmitter->GetEmissionPoints();
    const std::vector<std::array<FPNumber, 3>>& emissionVelocities = particleEmitter->GetEmissionVelocities();

    std::size_t numEmissions = numOfEmittedParticles.size();
    assert(emissionPoints.size() == numEmissions && emissionVelocities.size() == numEmissions);
    std::array<FPNumber, 3> force{0.0, 0.0, 0.0};

    for(std::size_t i = 0; i < numEmissions; ++i) {
        AddParticle(mass*numOfEmittedParticles[i], emissionPoints[i], emissionVelocities[i], force);
    }
}

void ParticlesTracer::SetConstrainingGeometry(Geometry* geometry, bool keepInside) {
    constrainingGeometry = geometry;
    keepInsideGeometry = keepInside;
}


void ParticlesTracer::UpdateParticlesConstraintStatus() {
    if(constrainingGeometry != nullptr) {
        constrainingGeometry->ArePointsInsideOrOn(positions, arePointsInside);
    }
}

void ParticlesTracer::UpdateParticlesPositions(const FPNumber dt, bool applyConstraints) {
    for(std::size_t i = 0; i < positions.size(); ++i) {
        std::array<FPNumber, 3>& r = positions[i];
        std::array<FPNumber, 3>& v = velocities[i];
        r[0] += v[0]*dt;
        r[1] += v[1]*dt;
        r[2] += v[2]*dt;
    }
    if(applyConstraints && constrainingGeometry != nullptr) {
        UpdateParticlesConstraintStatus();
        for(std::size_t i = 0; i < positions.size(); ++i) {
            if(arePointsInside[i] != keepInsideGeometry) {
                std::array<FPNumber, 3>& r = positions[i];
                std::array<FPNumber, 3>& v = velocities[i];
                r[0] -= v[0]*dt;
                r[1] -= v[1]*dt;
                r[2] -= v[2]*dt;
                v[0] = 0.0;
                v[1] = 0.0;
                v[2] = 0.0;
            }
        }
    }
}

void ParticlesTracer::UpdateParticlesMomentumsAndVelocities(const FPNumber dt) {
    for(std::size_t i = 0; i < positions.size(); ++i) {
        std::array<FPNumber, 3>& p = momentums[i];
        std::array<FPNumber, 3>& f = forces[i];
        p[0] += f[0]*dt;
        p[1] += f[1]*dt;
        p[2] += f[2]*dt;

        FPNumber mSquared = masses[i]*masses[i];
        FPNumber pSquared = p[0]*p[0] + p[1]*p[1] + p[2]*p[2];
        std::array<FPNumber, 3>& v = velocities[i];
        FPNumber mRel = std::sqrt(mSquared + pSquared);
        v[0] = p[0]/mRel;
        v[1] = p[1]/mRel;
        v[2] = p[2]/mRel;
    }
}

void ParticlesTracer::UpdateScatteringForce(int direction) {
    if(scatteringRateField == nullptr) {
        return;
    }

    const std::array<FPNumber, 3>& r0 = scatteringRateFieldConponentsOrigin[direction];
    const std::array<FPNumber, 3>& dr = gridSpacing;

    const std::size_t numOfParticles = masses.size();
    NumberArray3D<FPNumber>* sr_direction;
    if(scatteringRateField->GetElemType() == ElementType::NodeScalar) {
        sr_direction = &(scatteringRateField->GetNumArray(0));
    } else {
        sr_direction = &(scatteringRateField->GetNumArray(direction));
    }

    std::vector<FPNumber> sr_interpolated(numOfParticles);
    UniformGridInterpolator::InterpolateGridOnPoints(*sr_direction, r0, dr, positions, sr_interpolated);

    for(std::size_t i = 0; i < numOfParticles; ++i) {
          forces[i][direction] -= momentums[i][direction]*sr_interpolated[i];
    }
}


void ParticlesTracer::UpdateTime(const FPNumber newTime) {
    time = newTime;
}

void ParticlesTracer::UpdateParticlesMomentumVelocityPosition(const FPNumber newTime) {
    const FPNumber dt = newTime - time;
    UpdateParticlesMomentumsAndVelocities(dt);
    UpdateParticlesPositions(dt);
    // time must be updated at the end by the calling object using UpdateTime..
}


void ParticlesTracer::SetForce(std::size_t index, std::array<FPNumber, 3>& force) {
    forces[index] = force;
}

void ParticlesTracer::ResetForces() {
    for(std::size_t i = 0; i < positions.size(); ++i) {
        forces[i][0] = 0.0;
        forces[i][1] = 0.0;
        forces[i][2] = 0.0;
    }
}

