

#include "PhysicalUnits.hpp"
#include "FowlerNordheimEmission.hpp"

void test_Fowler_Nordheim() {
    // TODO: separate units test from emission test
    PhysicalUnits ps(100.0e-6);
    double eField_si = 3.0e9; // v/m
    double eField_fd = ps.ConvertSIElectricFieldToFDUnits(eField_si);
    std::cout << "electron mass SI : " << ps.GetElectronMassInSIUnits() << std::endl;
    std::cout << "electron mass FD : " << ps.GetElectronMassInFDUnits() << std::endl;

    std::cout << "electric field SI : " << eField_si << std::endl;
    std::cout << "electric field FD : " << eField_fd << std::endl;

    double workFunction_eV = 4.5;
    FowlerNordheimEmission fn(eField_si, workFunction_eV);

    int ny = 20;
    double dy = 1.0/ny;
    for(int i = 0; i < ny; ++i) {
        double y = (double)i * dy;
        double vy = fn.GetParameterVy(y);

        std::cout << "y : " << y << " , v(y) : " << vy <<std::endl;
    }

    double numOfElectrons = fn.GetNumberOfEmittedElectrons(1.0e-7, 1.0e-12);
    std::cout << "numOfElectrons : " << numOfElectrons << std::endl;

};
