
#ifndef FDTD_FOWLERNORDHEIM_H_
#define FDTD_FOWLERNORDHEIM_H_

#define __STDCPP_WANT_MATH_SPEC_FUNCS__

#include <cmath>

#include "PhysicalConstants.hpp"

class FowlerNordheimEmission {
    public:
    FowlerNordheimEmission(double eField_v_m = 0.0, double workFunction_eV = 1.0) {
        eField_voltsPerCentimeter = eField_v_m * 1.0e-2;
        workFunction_electronVolts = workFunction_eV;
    };

    void SetWorkFunction(double workFunction_eV) {
        workFunction_electronVolts = workFunction_eV;
    };

    void SetEfield(double eField_v_m) {
        eField_voltsPerCentimeter = eField_v_m * 1.0e-2;
    };

    double GetParameterVy(double y) {
        double _1_y2_sq = std::sqrt(1 - y*y);
        double k = std::sqrt(2.0*_1_y2_sq / (1.0 + _1_y2_sq));

        double ellipInt_1 = std::comp_ellint_1(k);
        double ellipInt_2 = std::comp_ellint_2(k);

        double v = std::sqrt((1.0 + _1_y2_sq) / 2.0) *
                   (ellipInt_2 - y*y*ellipInt_1/(1.0 + _1_y2_sq));
        if(y == 0.0) {
            v = std::sqrt((1.0 + _1_y2_sq) / 2.0) *
                   (ellipInt_2);
        } else if(y >= 1.0) {
            v = 0.0;
        }

        return v;
    };

    double GetEmissionCurrent() {
        double y = fn_c * std::sqrt(eField_voltsPerCentimeter) / workFunction_electronVolts;
        double v_y = GetParameterVy(y);
        double current_A_cm2 = fn_a * eField_voltsPerCentimeter*eField_voltsPerCentimeter / workFunction_electronVolts
            * std::exp(-fn_b * workFunction_electronVolts*std::sqrt(workFunction_electronVolts)
                             / eField_voltsPerCentimeter * v_y);
        double current_A_m2 = current_A_cm2 * 1.0e4;
        //std::cout << "y: " << y << " , v_y: " << v_y << std::endl;
        return current_A_m2;
    };

    double GetNumberOfEmittedElectrons(double surfaceArea_m2, double timeInterval_sec) {
        double current_A_m2 = GetEmissionCurrent();
        double charge_Coulomb = current_A_m2 * surfaceArea_m2 * timeInterval_sec;
        double numberOfElectrons = charge_Coulomb / PhysicalConstants_SI::electronCharge;
        //std::cout << "current_A_m2: " << current_A_m2 << " , charge_Coulomb: " << charge_Coulomb << " , numberOfElectrons: " << numberOfElectrons << std::endl;
        return numberOfElectrons;
    };

    private:
    const double fn_a = 1.55e-6;
    const double fn_b = 6.86e7;
    const double fn_c = 3.62e-4;

    double eField_voltsPerCentimeter;
    double workFunction_electronVolts;

};

#endif // FDTD_FOWLERNORDHEIM_H_

