

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

void ChargedParticlesTracer::SetMaxChargedPartcileBunchSize(std::size_t bunchSize) {
    maxChargedParticleBunchSize = bunchSize;
}


void ChargedParticlesTracer::AddParticlesEmittedByTheParticleEmitter(FPNumber t) {
    if(particleEmitter == nullptr) {
        return;
    }

    std::unordered_map<std::string, FPNumber> particleParams = particleEmitter->GetParticleParameters();
    const FPNumber mass = particleParams["mass"];
    const FPNumber charge = particleParams["charge"];
    const std::vector<FPNumber>& numOfEmittedParticles = particleEmitter->GetEmissionNumber(t);
    const std::vector<std::array<FPNumber, 3>>& emissionPoints = particleEmitter->GetEmissionPoints();
    const std::vector<std::array<FPNumber, 3>>& emissionVelocities = particleEmitter->GetEmissionVelocities();
    std::vector<std::vector<std::array<FPNumber, 3>>>* emissionSubPoints = particleEmitter->GetEmissionSubPoints();

    std::size_t numEmissions = numOfEmittedParticles.size();
    assert(emissionPoints.size() == numEmissions && emissionVelocities.size() == numEmissions);
    std::array<FPNumber, 3> force{0.0, 0.0, 0.0};

    if(chargeParticleEmissionNumberTracker.size() != numEmissions) {
        chargeParticleEmissionNumberTracker.resize(numEmissions, 0.0);
    }

    std::size_t n_particle_previous = charges.size();

    for(std::size_t i = 0; i < numEmissions; ++i) {
        FPNumber numParticle_i = numOfEmittedParticles[i];
        const std::array<FPNumber, 3>& velocity_i = emissionVelocities[i];
        if(numParticle_i == 0.0) {
            chargeParticleEmissionNumberTracker[i] = 0.0;
        } else {
            chargeParticleEmissionNumberTracker[i] += numParticle_i;
            numParticle_i = chargeParticleEmissionNumberTracker[i];
        }

        if(numParticle_i >= 1.0) {
            if(emissionSubPoints != nullptr && numParticle_i > 2*maxChargedParticleBunchSize
                                            && maxChargedParticleBunchSize >= 1) {
                // bunch particles
                auto& emissionSubPoints_i = emissionSubPoints->operator[](i);
                std::size_t numOfSubPts = emissionSubPoints_i.size();
                std::size_t numOfBunches = numParticle_i / maxChargedParticleBunchSize;
                if(numOfBunches > numOfSubPts) {
                    numOfBunches = numOfSubPts;
                }

                FPNumber bunchSize = numParticle_i / numOfBunches;   // number of particles in each bunch
                std::size_t subPtsStep = std::floor(numOfSubPts / numOfBunches);
                assert(subPtsStep >= 1);
                for(std::size_t j = 0; j < numOfBunches; j += subPtsStep) {
                    AddParticle(charge*bunchSize, mass*bunchSize, emissionSubPoints_i[j], velocity_i, force);
                }
                chargeParticleEmissionNumberTracker[i] = 0.0;

            } else if(numParticle_i >= maxChargedParticleBunchSize || maxChargedParticleBunchSize < 1.0) {
                AddParticle(charge*numParticle_i, mass*numParticle_i, emissionPoints[i], velocity_i, force);
                chargeParticleEmissionNumberTracker[i] = 0.0;
            }
        }
    }

    if(charges.size() % 10 == 1 && charges.size() > n_particle_previous) {
        std::cout << "num of particles : " << charges.size() << " , bunch size: " << maxChargedParticleBunchSize << std::endl;
    }
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

void ChargedParticlesTracer::SetAnalyticElectricField(VectorField* eField) {
    analyticElectricField = eField;
}

void ChargedParticlesTracer::SetAnalyticMagneticField(VectorField* bField) {
    analyticMagneticField = bField;
}

void ChargedParticlesTracer::UpdateElectricForce(int direction) {

    const std::array<FPNumber, 3>& r0 = electricFieldConponentsOrigin[direction];
    const std::array<FPNumber, 3>& dr = gridSpacing;

    const std::size_t numOfParticles = charges.size();

    std::vector<FPNumber> e_interpolated(numOfParticles, 0.0);
    if(electricField != nullptr) {
        const NumberArray3D<FPNumber>& e_direction = electricField->GetNumArray(direction);
        UniformGridInterpolator::InterpolateGridOnPoints(e_direction, r0, dr, positions, e_interpolated);
    }
    if(analyticElectricField != nullptr) {
        std::vector<std::array<FPNumber, 3>> e_analytic;

        analyticElectricField->GetFieldValuesAtPoints(time, positions, e_analytic);

        for(std::size_t i = 0; i < numOfParticles; ++i) {
            e_interpolated[i] += e_analytic[i][direction];
        }
    }

    for(std::size_t i = 0; i < numOfParticles; ++i) {
          forces[i][direction] += charges[i]*e_interpolated[i];
    }
}


void ChargedParticlesTracer::UpdateMagneticForce(int direction) {

    const std::array<FPNumber, 3>& r0 = magneticFieldConponentsOrigin[direction];
    const std::array<FPNumber, 3>& dr = gridSpacing;

    std::size_t numOfParticles = charges.size();

    std::vector<FPNumber> b_interpolated(numOfParticles, 0.0);
    if(magneticField != nullptr) {
        const NumberArray3D<FPNumber>& b_direction = magneticField->GetNumArray(direction);
        UniformGridInterpolator::InterpolateGridOnPoints(b_direction, r0, dr, positions, b_interpolated);
    }
    if(analyticMagneticField != nullptr) {
        std::vector<std::array<FPNumber, 3>> b_analytic;

        analyticMagneticField->GetFieldValuesAtPoints(time, positions, b_analytic);

        for(std::size_t i = 0; i < numOfParticles; ++i) {
            b_interpolated[i] += b_analytic[i][direction];
        }
    }

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
    //FPNumber dA_xy = dr[0]*dr[1];
    //FPNumber dA_yz = dr[1]*dr[2];
    //FPNumber dA_xz = dr[0]*dr[2];
    for(std::size_t i = 0; i < numOfParticles; ++i) {
        const FPNumber& q = charges[i];
        std::array<FPNumber, 3>& v = velocities[i];
        Jx[i] = q*v[0];// /dA_yz;   warning: these currents should be treated as point currents as opposed to volumetric currents inside FDTD
        Jy[i] = q*v[1];// /dA_xz;
        Jz[i] = q*v[2];// /dA_xy;
    }

    //std::cout << numOfParticles << " ";
}


void ChargedParticlesTracer::AttachDataToGAMPositions(std::vector<std::array<FPNumber, 3>>*& positions) {
    positions = &(this->positions);
}

void ChargedParticlesTracer::AttachScalarDataToGAMValues(std::vector<FPNumber>*& values, std::string dataName, int direction) {
    if(dataName == "current") {
        values = &currentComponents[direction];
    } else if(dataName == "charge") {
        values = &charges;
    } else {
        std::cout << "error: " << dataName << " is not a valid data name." << std::endl;
        assert(false);
    }
}

void ChargedParticlesTracer::AttachVectorDataToGAMValues(std::vector<std::array<FPNumber, 3>>*& values, std::string dataName, int direction) {
    if(dataName == "velocity") {
        values = &velocities;
    } else if(dataName == "momentum") {
        values = &momentums;
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
            UpdateScatteringForce(direction);
            UpdateElectricForce(direction);
            UpdateMagneticForce(direction);
        }
        UpdateParticlesMomentumVelocityPosition(t);
        UpdateParticlesCurrents();

        UpdateTime(t);
        particleEmitter->UpdateTime(t);
    }
}

