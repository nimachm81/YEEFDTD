#ifndef FDTD_PHYSICALCONSTANTS_H_
#define FDTD_PHYSICALCONSTANTS_H_


class PhysicalConstants_SI {
    public:
    constexpr static double speedOfLight = 299792458.0;     // m/s
    constexpr static double electronCharge = 1.60217662e-19; // Coulomb
    constexpr static double electronMass = 9.10938356e-31;    // Kg
    constexpr static double vacuumPermittivity = 8.854187817e-12; // F/m (Farad per meter)
    constexpr static double vacuumPermeability = 4.0*M_PI*1.0e-7; // H/m (Henry per meter)

};

#endif // FDTD_PHYSICALCONSTANTS_H_

