
#include "UniformGridInterpolator.hpp"
#include "ChargedParticleEmitter.h"

void ChargedParticleEmitter::SetElementaryCharge(FPNumber charge) {
    particleElementaryParameters["charge"] = charge;
}

void ChargedParticleEmitter::SetElementaryMass(FPNumber mass) {
    particleElementaryParameters["mass"] = mass;
}

void ChargedParticleEmitter::SetUnitLength(FPNumber length) {
    unitConverter.SetFDLengthInSIUnits(length);
}

void ChargedParticleEmitter::SetWorkFunction(FPNumber workFunction) {
    fnEmitter.SetWorkFunction(workFunction);
}

void ChargedParticleEmitter::SetEmissionPoints(std::vector<std::array<FPNumber, 3>>& emissionPts) {
    emissionPoints = emissionPts;
}

void ChargedParticleEmitter::SetNormalVectors(std::vector<std::array<FPNumber, 3>>& normalVecs) {
    normalVectors = normalVecs;
}

void ChargedParticleEmitter::SetSurfaceAreas(std::vector<FPNumber>& areas) {
    surfaceAreas = areas;
}

void ChargedParticleEmitter::SetElectricField(YeeGridData3D* eField) {
    electricField = eField;
}

void ChargedParticleEmitter::SetElectricFieldGridOrigin(int direction, std::array<FPNumber, 3>& origin) {
    electricFieldConponentsOrigin[direction] = origin;
}

void ChargedParticleEmitter::SetGridSpacing(std::array<FPNumber, 3>& dr) {
    gridSpacing = dr;
}

FPNumber ChargedParticleEmitter::GetFowlerNordheimEmissionNumber(FPNumber eFieldNormal, FPNumber surfaceArea,
                                                             FPNumber timeInterval) {
    FPNumber eFieldNormal_SI = unitConverter.ConvertFDElectricFieldToSIUnits(eFieldNormal);
    FPNumber surfaceArea_SI = unitConverter.ConvertFDLengthToSIUnits(surfaceArea);  // TODO: good for 2D.. 3D should be treated using area instead of arc length
    FPNumber timeInterval_SI = unitConverter.ConvertFDTimeToSIUnits(timeInterval);

    fnEmitter.SetEfield(eFieldNormal_SI);
    FPNumber numOfParticles = fnEmitter.GetNumberOfEmittedElectrons(surfaceArea_SI, timeInterval_SI);

    if(numOfParticles > 1000) {
        std::cout << "E_SI: " << eFieldNormal_SI << " , dA_SI: " << surfaceArea_SI << " , dt_SI: " << timeInterval_SI << " , n_e : " <<  numOfParticles << std::endl;
    }

    return numOfParticles;
}


const std::vector<FPNumber>& ChargedParticleEmitter::GetEmissionNumber(FPNumber t) {
    assert(t > time);
    FPNumber dt = t - time;

    std::size_t numOfPoints = emissionPoints.size();
    emissionNumbers.resize(numOfPoints);
    emissionVelocities.resize(numOfPoints);

    std::array<std::vector<FPNumber>, 3> e_interp;
    for(int i = 0; i < 3; ++i) {
        UniformGridInterpolator::InterpolateGridOnPoints(electricField->GetNumArray(i),
                                                         electricFieldConponentsOrigin[i],
                                                         gridSpacing,
                                                         emissionPoints,
                                                         e_interp[i]);
    }

    const std::vector<FPNumber>& ex_interp = e_interp[0];
    const std::vector<FPNumber>& ey_interp = e_interp[1];
    const std::vector<FPNumber>& ez_interp = e_interp[2];
    for(std::size_t i = 0; i < numOfPoints; ++i) {
        std::array<FPNumber, 3>& v_norm = normalVectors[i];
        FPNumber e_normal = ex_interp[i]*v_norm[0] + ey_interp[i]*v_norm[1] + ez_interp[i]*v_norm[2];

        if(e_normal < 0) {
            emissionNumbers[i] = GetFowlerNordheimEmissionNumber(std::abs(e_normal), surfaceAreas[i], dt);
            //if(emissionNumbers[i] > 1) {std::cout << "E: " << std::abs(e_normal) << " , dA: " << surfaceAreas[i] << " , dt: " << dt << " , n_e: " << emissionNumbers[i] << std::endl;}
        } else {
            emissionNumbers[i] = 0.0;
        }
        //std::cout << "em num:" << emissionNumbers[i] << ", " ;

        emissionVelocities[i] = std::array<FPNumber, 3>{0.0, 0.0, 0.0};
     }

    return emissionNumbers;
}

const std::vector<std::array<FPNumber, 3>>& ChargedParticleEmitter::GetEmissionPoints() {
    return emissionPoints;
}

const std::vector<std::array<FPNumber, 3>>& ChargedParticleEmitter::GetEmissionVelocities() {
    return emissionVelocities;
}

const std::unordered_map<std::string, FPNumber>& ChargedParticleEmitter::GetParticleParameters() {
    return particleElementaryParameters;
}


