#ifndef FDTD_PHYSICALUNITS_H_
#define FDTD_PHYSICALUNITS_H_


#include "PhysicalConstants.hpp"


/* Base units: length (L), velocity of light (c), electron charge (e)
    c = 1.0
    e = 1.0
    L = 1.0
    epsilon_0 = 1.0
    mu_0 = 1.0
*/
class PhysicalUnits {
    public:
    PhysicalUnits(double length_si = 1.0) {
        SetFDLengthInSIUnits(length_si);
    };

    void SetFDLengthInSIUnits(double l_si) {
        length_SI = l_si;
        // set conversion factors
        length_SI_to_FD = length_SI / length_FD;
        area_SI_to_FD = length_SI_to_FD * length_SI_to_FD;
        time_SI_to_FD = (length_SI / speedOfLight_SI) / (length_FD / speedOfLight_FD);
        electricPotential_SI_to_FD = (electronCharge_SI / (vacuumPermittivity_SI * length_SI)) / (electronCharge_FD / (vacuumPermittivity_FD * length_FD));
        electricField_SI_to_FD = (electronCharge_SI / (vacuumPermittivity_SI * length_SI * length_SI)) / (electronCharge_FD / (vacuumPermittivity_FD * length_FD * length_FD));
        electricCurrent_SI_to_FD = (electronCharge_SI / electronCharge_FD) / time_SI_to_FD;
        //energy_SI_to_FD = (vacuumPermittivity_SI / vacuumPermittivity_FD) * electricField_SI_to_FD * electricField_SI_to_FD
        //                        * length_SI*length_SI*length_SI;
        energy_SI_to_FD = (electronCharge_SI / electronCharge_FD) * electricPotential_SI_to_FD;
        double c_SI_to_FD = speedOfLight_SI / speedOfLight_FD;
        mass_SI_to_FD = energy_SI_to_FD /  (c_SI_to_FD * c_SI_to_FD);
    };

    double ConvertFDLengthToSIUnits(double l_fd) {
        return l_fd * length_SI_to_FD;
    };

    double ConvertSILengthToFDUnits(double l_si) {
        return l_si / length_SI_to_FD;
    };

    double ConvertFDAreaToSIUnits(double a_fd) {
        return a_fd * area_SI_to_FD;
    };

    double ConvertSIAreaToFDUnits(double a_si) {
        return a_si / area_SI_to_FD;
    };

    double ConvertFDTimeToSIUnits(double t_fd) {
        return t_fd * time_SI_to_FD;
    };

    double ConvertSITimeToFDUnits(double t_si) {
        return t_si / time_SI_to_FD;
    };

    double ConvertFDElectricPotentialToSIUnits(double v_fd) {
        return v_fd * electricPotential_SI_to_FD;
    };

    double ConvertSIElectricPotentialToFDUnits(double v_si) {
        return v_si / electricPotential_SI_to_FD;
    };

    double ConvertFDElectricFieldToSIUnits(double eField_fd) {
        return eField_fd * electricField_SI_to_FD;
    };

    double ConvertSIElectricFieldToFDUnits(double eField_si) {
        return eField_si / electricField_SI_to_FD;
    };

    double GetElectronMassInFDUnits() {
        return electronMass_SI / mass_SI_to_FD;
    };

    double GetElectronMassInSIUnits() {
        return electronMass_SI;
    };

    double GetElectronChargeInFDUnits() {
        return electronCharge_FD;
    };

    double GetElectronChargeInSIUnits() {
        return electronCharge_SI;
    };

    double ConvertSIEnergyToFDUnits(double energy_si) {
        return energy_si / energy_SI_to_FD;
    };

    double ConvertFDEnergyToSIUnits(double energy_fd) {
        return energy_fd * energy_SI_to_FD;
    };

    private:
    // FDTD units
    const double speedOfLight_FD = 1.0;
    const double electronCharge_FD = 1.0;
    const double length_FD = 1.0;
    const double vacuumPermittivity_FD = 1.0;
    const double vacuumPermeability_FD = 1.0;

    // SI Units
    const double speedOfLight_SI = PhysicalConstants_SI::speedOfLight;
    const double electronCharge_SI = PhysicalConstants_SI::electronCharge;
    const double electronMass_SI = PhysicalConstants_SI::electronMass;
    const double vacuumPermittivity_SI = PhysicalConstants_SI::vacuumPermittivity;
    const double vacuumPermeability_SI = PhysicalConstants_SI::vacuumPermeability;
    double length_SI;     // FDTD unit length in SI units, it can be set based on a specific wavelegnth for example

    // conversion factors
    // _SI_to_FD : SI unit/ FD unit ratio
    double length_SI_to_FD;
    double area_SI_to_FD;
    double time_SI_to_FD;
    double electricPotential_SI_to_FD;
    double electricField_SI_to_FD;
    double electricCurrent_SI_to_FD;
    double energy_SI_to_FD;
    double mass_SI_to_FD;
};

#endif // FDTD_PHYSICALUNITS_H_
