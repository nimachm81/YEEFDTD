#ifndef FDTD_CHARGEDPARTICLESTRACER_H_
#define FDTD_CHARGEDPARTICLESTRACER_H_


#include "YeeGridDataTypes.h"
#include "NumberTypes.h"
#include "DiscretePointsGridArrayManipulator.h"
#include "ParticlesTracer.h"

class ChargedParticlesTracer : ParticlesTracer {
    public:
    void AddParticle(const FPNumber charge,
                     const FPNumber mass,
                     const std::array<FPNumber, 3>& position,
                     const std::array<FPNumber, 3>& velocity,
                     const std::array<FPNumber, 3>& force);
    void SetElectricFieldGrid(YeeGridData3D* eField);
    void SetMagneticFieldGrid(YeeGridData3D* bField);

    void AttachCurrentToDiscreteGridArrayManipulator(DiscretePointsGridArrayManipulator& discreteGAM, int direction);
    void UpdateElectricForce(int direction, // x, y or z component. Directiion of the electric field
                             std::array<FPNumber, 3>& r0,   // coordinates of the first component of the electric field grid
                             std::array<FPNumber, 3>& dr
                             );
    void UpdateMagneticForce(int direction, // 0=x, 1=y or 2=z component. Directiion of the magnetic field
                             std::array<FPNumber, 3>& r0,   // coordinates of the first component of the electric field grid
                             std::array<FPNumber, 3>& dr
                             );
    private:
    std::vector<FPNumber> charges;
    std::array<std::vector<FPNumber>, 3> currentComponents;
    YeeGridData3D* electricField;       // E
    YeeGridData3D* magneticField;       // B
};

#endif // FDTD_CHARGEDPARTICLESTRACER_H_
