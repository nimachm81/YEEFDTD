

#include "ChargedParticlesTracer.h"

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

    std::size_t numOfParticles = charges.size();
    NumberArray3D e_direction = electricField->GetNumArray(direction);
    for(std::size_t i = 0; i < numOfParticles; ++i) {
        std::array<FPNumber, 3>& position = positions[i];

        std::array<std::size_t, 3> lowerIndx
                           {static_cast<std::size_t>((position[0] - r0[0]) / dr[0]),
                            static_cast<std::size_t>((position[1] - r0[1]) / dr[1]),
                            static_cast<std::size_t>((position[2] - r0[2]) / dr[2])};

        std::array<FPNumber, 3> alpha;
        for(int j = 0; j < 3; ++j) {
            alpha[j] = ((position[j] - r0[j]) - lowerIndx[j]*dr[j])/dr[j];
        }

        // assert inside bounds
        const std::array<std::size_t, 3>& gridArrayShape = e_direction.GetShape();

        std::array<std::size_t, 3> indx;
        std::array<FPNumber, 3> vol_ratio;
        //FPNumber[2][2][2] weights;

        FPNumber eInterp = 0.0;     // interpolated electric field

        for(int i_x = 0; i_x < 2; ++i_x) {
            indx[0] = lowerIndx[0] + i_x;
            if(indx[0] < 0 || indx[0] >= gridArrayShape[0]) {
                continue;
            }
            if(i_x == 0) {
                vol_ratio[0] = 1.0 - alpha[0];
            } else {
                vol_ratio[0] = alpha[0];
            }
            for(int i_y = 0; i_y < 2; ++i_y) {
                indx[1] = lowerIndx[1] + i_y;
                if(indx[1] < 0 || indx[1] >= gridArrayShape[1]) {
                    continue;
                }
                if(i_y == 0) {
                    vol_ratio[1] = 1.0 - alpha[1];
                } else {
                    vol_ratio[1] = alpha[1];
                }
                for(int i_z = 0; i_z < 2; ++i_z) {
                    indx[2] = lowerIndx[2] + i_x;
                    if(indx[2] < 0 || indx[2] >= gridArrayShape[2]) {
                        continue;
                    }
                    if(i_z == 0) {
                        vol_ratio[2] = 1.0 - alpha[2];
                    } else {
                        vol_ratio[2] = alpha[2];
                    }

                    eInterp += e_direction[indx]*vol_ratio[0]*vol_ratio[1]*vol_ratio[2];
                }
            }
        }

        forces[i][direction] += charges[i]*eInterp;
    }
}


void ChargedParticlesTracer::UpdateMagneticForce(int direction) {
    const std::array<FPNumber, 3>& r0 = magneticFieldConponentsOrigin[direction];
    const std::array<FPNumber, 3>& dr = gridSpacing;

    std::size_t numOfParticles = charges.size();
    NumberArray3D b_direction = magneticField->GetNumArray(direction);
    for(std::size_t i = 0; i < numOfParticles; ++i) {
        std::array<FPNumber, 3>& position = positions[i];

        std::array<std::size_t, 3> lowerIndx
                           {static_cast<std::size_t>((position[0] - r0[0]) / dr[0]),
                            static_cast<std::size_t>((position[1] - r0[1]) / dr[1]),
                            static_cast<std::size_t>((position[2] - r0[2]) / dr[2])};

        std::array<FPNumber, 3> alpha;
        for(int j = 0; j < 3; ++j) {
            alpha[j] = ((position[j] - r0[j]) - lowerIndx[j]*dr[j])/dr[j];
        }

        // assert inside bounds
        const std::array<std::size_t, 3>& gridArrayShape = b_direction.GetShape();

        std::array<std::size_t, 3> indx;
        std::array<FPNumber, 3> vol_ratio;
        //FPNumber[2][2][2] wights;

        FPNumber bInterp = 0.0;

        for(int i_x = 0; i_x < 2; ++i_x) {
            indx[0] = lowerIndx[0] + i_x;
            if(indx[0] < 0 || indx[0] >= gridArrayShape[0]) {
                continue;
            }
            if(i_x == 0) {
                vol_ratio[0] = 1.0 - alpha[0];
            } else {
                vol_ratio[0] = alpha[0];
            }
            for(int i_y = 0; i_y < 2; ++i_y) {
                indx[1] = lowerIndx[1] + i_y;
                if(indx[1] < 0 || indx[1] >= gridArrayShape[1]) {
                    continue;
                }
                if(i_y == 0) {
                    vol_ratio[1] = 1.0 - alpha[1];
                } else {
                    vol_ratio[1] = alpha[1];
                }
                for(int i_z = 0; i_z < 2; ++i_z) {
                    indx[2] = lowerIndx[2] + i_x;
                    if(indx[2] < 0 || indx[2] >= gridArrayShape[2]) {
                        continue;
                    }
                    if(i_z == 0) {
                        vol_ratio[2] = 1.0 - alpha[2];
                    } else {
                        vol_ratio[2] = alpha[2];
                    }

                    bInterp += b_direction[indx]*vol_ratio[0]*vol_ratio[1]*vol_ratio[2];
                }
            }
        }

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
        std::array<FPNumber, 3>& v = velocities[i];
        Jx[i] = v[0]/dA_yz;
        Jy[i] = v[1]/dA_xz;
        Jz[i] = v[2]/dA_xy;
    }
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
    ResetForces();
    for(int direction = 0; direction < 3; ++direction) {
        UpdateElectricForce(direction);
        UpdateMagneticForce(direction);
    }
    UpdateParticlesMomentumVelocityPosition(t);
    UpdateParticlesCurrents();
}

