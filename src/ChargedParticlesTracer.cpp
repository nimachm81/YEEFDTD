

#include "UniformGridInterpolator.hpp"
#include "ChargedParticlesTracer.h"

void ChargedParticlesTracer::ReserveMemory(std::size_t numberOfElements) {
    ParticlesTracer::ReserveMemory(numberOfElements);
    charges.reserve(numberOfElements);
    currentComponents[0].reserve(numberOfElements);
    currentComponents[1].reserve(numberOfElements);
    currentComponents[2].reserve(numberOfElements);
}

void ChargedParticlesTracer::AddParticle(const FPNumber charge,
                                         const FPNumber mass,
                                         const std::array<FPNumber, 3>& position,
                                         const std::array<FPNumber, 3>& velocity,
                                         const std::array<FPNumber, 3>& force) {
    ParticlesTracer::AddParticle(mass, position, velocity, force);
    charges.push_back(charge);
    currentComponents[0].push_back(charge*velocity[0]);
    currentComponents[1].push_back(charge*velocity[1]);
    currentComponents[2].push_back(charge*velocity[2]);
}

void ChargedParticlesTracer::AddParticlesEmittedByTheParticleEmitter(FPNumber t,
                                                                    std::size_t bunchSize
                                                                    ) {
    if(particleEmitter == nullptr) {
        return;
    }

    std::unordered_map<std::string, FPNumber> particleParams = particleEmitter->GetParticleParameters();
    const FPNumber mass = particleParams["mass"];
    const FPNumber charge = particleParams["charge"];
    const std::vector<FPNumber>& numOfEmittedParticles = particleEmitter->GetEmissionNumber(t);
    const std::vector<std::array<FPNumber, 3>>& emissionPoints = particleEmitter->GetEmissionPoints();
    const std::vector<std::array<FPNumber, 3>>& emissionVelocities = particleEmitter->GetEmissionVelocities();

    std::size_t numEmissions = numOfEmittedParticles.size();
    assert(emissionPoints.size() == numEmissions && emissionVelocities.size() == numEmissions);
    std::array<FPNumber, 3> force{0.0, 0.0, 0.0};

    for(std::size_t i = 0; i < numEmissions; ++i) {
        if(numOfEmittedParticles[i] > 1.0) {
            for(std::size_t i_p = 0; i_p < numOfEmittedParticles[i]/bunchSize; ++i_p) {
                AddParticle(charge*bunchSize, mass*bunchSize, emissionPoints[i], emissionVelocities[i], force);
            }
        }
    }
}

void ChargedParticlesTracer::SetGridSpacing(std::array<FPNumber, 3>& dr) {
    gridSpacing = dr;
}

void ChargedParticlesTracer::SetElectricFieldGrid(YeeGridData3D* eField) {
    electricField = eField;
}

void ChargedParticlesTracer::SetMagneticFieldGrid(YeeGridData3D* bField) {
    magneticField = bField;
}

void ChargedParticlesTracer::SetElectricFieldGridOrigin(int direction, std::array<FPNumber, 3>& origin) {
    electricFieldConponentsOrigin[direction] = origin;
}

void ChargedParticlesTracer::SetMagneticFieldGridOrigin(int direction, std::array<FPNumber, 3>& origin) {
    magneticFieldConponentsOrigin[direction] = origin;
}

void ChargedParticlesTracer::UpdateElectricForce(int direction) {
    const std::array<FPNumber, 3>& r0 = electricFieldConponentsOrigin[direction];
    const std::array<FPNumber, 3>& dr = gridSpacing;

    const std::size_t numOfParticles = charges.size();
    const NumberArray3D<FPNumber>& e_direction = electricField->GetNumArray(direction);

    std::vector<FPNumber> e_interpolated(numOfParticles);
    UniformGridInterpolator::InterpolateGridOnPoints(e_direction, r0, dr, positions, e_interpolated);

    for(std::size_t i = 0; i < numOfParticles; ++i) {
          forces[i][direction] += charges[i]*e_interpolated[i];
    }
}


void ChargedParticlesTracer::UpdateMagneticForce(int direction) {
    const std::array<FPNumber, 3>& r0 = magneticFieldConponentsOrigin[direction];
    const std::array<FPNumber, 3>& dr = gridSpacing;

    std::size_t numOfParticles = charges.size();
    const NumberArray3D<FPNumber>& b_direction = magneticField->GetNumArray(direction);

    std::vector<FPNumber> b_interpolated(numOfParticles);
    UniformGridInterpolator::InterpolateGridOnPoints(b_direction, r0, dr, positions, b_interpolated);


    for(std::size_t i = 0; i < numOfParticles; ++i) {
        FPNumber bInterp = b_interpolated[i];

        // F = q v x B
        FPNumber leviCivita3[3][3][3];
        leviCivita3[0][1][2] = 1.0;
        leviCivita3[1][2][0] = 1.0;
        leviCivita3[2][0][1] = 1.0;
        leviCivita3[0][2][1] = -1.0;
        leviCivita3[1][0][2] = -1.0;
        leviCivita3[2][1][0] = -1.0;

        for(int j = 0; j < 3; ++j) {
            if(j == direction) {continue;}
            for(int k = 0; k < 3; ++k) {
                if(k == j || k == direction) {continue;}
                forces[i][j] += charges[i]*leviCivita3[j][k][direction]*velocities[i][k]*bInterp;
            }
        }
    }
}

void ChargedParticlesTracer::UpdateParticlesCurrents() {
    const std::array<FPNumber, 3>& dr = gridSpacing;

    std::size_t numOfParticles = charges.size();
    std::vector<FPNumber>& Jx = currentComponents[0];
    std::vector<FPNumber>& Jy = currentComponents[1];
    std::vector<FPNumber>& Jz = currentComponents[2];
    FPNumber dA_xy = dr[0]*dr[1];
    FPNumber dA_yz = dr[1]*dr[2];
    FPNumber dA_xz = dr[0]*dr[2];
    for(std::size_t i = 0; i < numOfParticles; ++i) {
        const FPNumber& q = charges[i];
        std::array<FPNumber, 3>& v = velocities[i];
        Jx[i] = q*v[0]/dA_yz;
        Jy[i] = q*v[1]/dA_xz;
        Jz[i] = q*v[2]/dA_xy;
    }

    //std::cout << numOfParticles << " ";
}


void ChargedParticlesTracer::AttachDataToGAMPositions(std::vector<std::array<FPNumber, 3>>*& positions) {
    positions = &(this->positions);
}

void ChargedParticlesTracer::AttachDataToGAMValues(std::vector<FPNumber>*& values, std::string dataName, int direction) {
    if(dataName == "current") {
        values = &currentComponents[direction];
    } else {
        std::cout << "error: " << dataName << " is not a valid data name." << std::endl;
        assert(false);
    }
}

void ChargedParticlesTracer::UpdateGAMValues(const FPNumber t) {
    if( t > time ) {
        AddParticlesEmittedByTheParticleEmitter(t);
        ResetForces();
        for(int direction = 0; direction < 3; ++direction) {
            UpdateElectricForce(direction);
            UpdateMagneticForce(direction);
        }
        UpdateParticlesMomentumVelocityPosition(t);
        UpdateParticlesCurrents();
    }
}

