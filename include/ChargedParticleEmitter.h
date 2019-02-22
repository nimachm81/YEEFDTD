#ifndef FDTD_CHARGEDPARTICLEEMITTER_H_
#define FDTD_CHARGEDPARTICLEEMITTER_H_

#include "NumberTypes.h"
#include "YeeGridDataTypes.h"
#include "PhysicalUnits.hpp"
#include "FowlerNordheimEmission.hpp"
#include "ParticleEmitter.h"

class ChargedParticleEmitter : public ParticleEmitter {
    public:
    virtual ~ChargedParticleEmitter() {};
    void SetElementaryCharge(FPNumber charge);
    void SetElementaryMass(FPNumber mass);
    void SetUnitLength(FPNumber length);        // used to convert between FDTD and Si units
    void SetWorkFunction(FPNumber workFunction);        // used to calculate emission rate

    void SetEmissionPoints(std::vector<std::array<FPNumber, 3>>& emissionPts);
    void SetNormalVectors(std::vector<std::array<FPNumber, 3>>& normalVecs);
    void SetSurfaceAreas(std::vector<FPNumber>& areas);

    void SetElectricField(YeeGridData3D* eField);
    void SetElectricFieldGridOrigin(int direction, std::array<FPNumber, 3>& origin);
    void SetGridSpacing(std::array<FPNumber, 3>& dr);

    FPNumber GetFowlerNordheimEmissionNumber(FPNumber eFieldNormal, // in FDTD units
                                         FPNumber surfaceArea,  // in FDTD units
                                         FPNumber timeInterval  // in FDTD units
                                         );

    // override ParticleEmitter's virtual functions
    const std::vector<FPNumber>& GetEmissionNumber(FPNumber t);       // get number of emitted particles per surface area at time t
    const std::vector<std::array<FPNumber, 3>>& GetEmissionPoints();
    const std::vector<std::array<FPNumber, 3>>& GetEmissionVelocities();

    const std::unordered_map<std::string, FPNumber>& GetParticleParameters();


    protected:
    std::vector<FPNumber> emissionNumbers;          // number of emitted particles

    std::vector<std::array<FPNumber, 3>> emissionPoints;     // point on the surface of the object
    std::vector<std::array<FPNumber, 3>> emissionVelocities;     // initial velocities of the emitted particles
    std::vector<std::array<FPNumber, 3>> normalVectors;       // normal to the surface (assumed normalized)
    std::vector<FPNumber> surfaceAreas;                       // surface area of a patch on the surface of the object

    std::unordered_map<std::string, FPNumber> particleElementaryParameters;   // elementary parameters of each particle
                                                                // such as unit charge and mass, for example:
                                                                // {"charge": q, "mass": m}

    //FPNumber temperature;
    PhysicalUnits unitConverter;      // used to convert between FDTD units and SI units
    FowlerNordheimEmission fnEmitter;       // used to calculate Fowler-Nordheim field emission
    YeeGridData3D* electricField;       // E
    std::array<std::array<FPNumber, 3>, 3> electricFieldConponentsOrigin;
    std::array<FPNumber, 3> gridSpacing;

};

#endif // FDTD_CHARGEDPARTICLEEMITTER_H_

