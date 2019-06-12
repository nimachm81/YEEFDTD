
__all__ = ["PhysicalUnits"]


from scipy import constants

""" Base units: length (L), velocity of light (c), electron charge (e)
    c = 1.0
    e = 1.0
    L = 1.0
    epsilon_0 = 1.0
    mu_0 = 1.0
"""

class PhysicalUnits:
    
    def __init__(self, length_si = 1.0):
            ## FDTD units
        self.speedOfLight_FD = 1.0
        self.electronCharge_FD = 1.0
        self.length_FD = 1.0
        self.vacuumPermittivity_FD = 1.0
        self.vacuumPermeability_FD = 1.0

        ## SI Units
        self.speedOfLight_SI = constants.c
        self.electronCharge_SI = constants.e
        self.electronMass_SI = constants.m_e
        self.vacuumPermittivity_SI = constants.epsilon_0
        self.vacuumPermeability_SI = constants.mu_0
        self.length_SI = None     # FDTD unit length in SI units, it can be set based on a specific wavelegnth for example

        ## conversion factors
        ## _SI_to_FD : SI unit/ FD unit ratio
        self.length_SI_to_FD = None    ## L_SI/L_FD
        self.area_SI_to_FD = None
        self.time_SI_to_FD = None
        self.frequency_SI_to_FD = None
        self.electricPotential_SI_to_FD = None
        self.electricField_SI_to_FD = None
        self.electricCurrent_SI_to_FD = None
        self.energy_SI_to_FD = None
        self.mass_SI_to_FD = None

        self.SetFDLengthInSIUnits(length_si)


    def SetFDLengthInSIUnits(self, l_si):
        self.length_SI = l_si
        
        ## set conversion factors
        ## length, time
        self.length_SI_to_FD = self.length_SI / self.length_FD
        self.area_SI_to_FD = self.length_SI_to_FD * self.length_SI_to_FD
        self.time_SI_to_FD = (self.length_SI / self.speedOfLight_SI) / (self.length_FD / self.speedOfLight_FD)
        self.frequency_SI_to_FD = 1.0 / self.time_SI_to_FD
        
        ## electric field, potential, current
        self.electricPotential_SI_to_FD = (self.electronCharge_SI / (self.vacuumPermittivity_SI * self.length_SI)) /   \
                                          (self.electronCharge_FD / (self.vacuumPermittivity_FD * self.length_FD))
        self.electricField_SI_to_FD = (self.electronCharge_SI / (self.vacuumPermittivity_SI * self.length_SI * self.length_SI)) /  \
                                      (self.electronCharge_FD / (self.vacuumPermittivity_FD * self.length_FD * self.length_FD))
        self.electricCurrent_SI_to_FD = (self.electronCharge_SI / self.electronCharge_FD) / self.time_SI_to_FD

        ## energy
        self.energy_SI_to_FD = (self.electronCharge_SI / self.electronCharge_FD) * self.electricPotential_SI_to_FD
        
        ## mass
        c_SI_to_FD = self.speedOfLight_SI / self.speedOfLight_FD
        self.mass_SI_to_FD = self.energy_SI_to_FD /  (c_SI_to_FD * c_SI_to_FD)
    

    def ConvertFDLengthToSIUnits(self, l_fd):
        return l_fd * self.length_SI_to_FD

    def ConvertSILengthToFDUnits(self, l_si):
        return l_si / self.length_SI_to_FD

    def ConvertFDAreaToSIUnits(self, a_fd):
        return a_fd * self.area_SI_to_FD

    def ConvertSIAreaToFDUnits(self, a_si):
        return a_si / self.area_SI_to_FD

    def ConvertFDTimeToSIUnits(self, t_fd):
        return t_fd * self.time_SI_to_FD

    def ConvertSITimeToFDUnits(self, t_si):
        return t_si / self.time_SI_to_FD

    def ConvertFDFrequencyToSIUnits(self, f_fd):
        return f_fd * self.frequency_SI_to_FD

    def ConvertSIFrequencyToFDUnits(self, f_si):
        return f_si / self.frequency_SI_to_FD

    def ConvertFDElectricPotentialToSIUnits(self, v_fd):
        return v_fd * self.electricPotential_SI_to_FD

    def ConvertSIElectricPotentialToFDUnits(self, v_si):
        return v_si / self.electricPotential_SI_to_FD

    def ConvertFDElectricFieldToSIUnits(self, eField_fd):
        return eField_fd * self.electricField_SI_to_FD

    def ConvertSIElectricFieldToFDUnits(self, eField_si):
        return eField_si / self.electricField_SI_to_FD

    def GetElectronMassInFDUnits(self):
        return self.electronMass_SI / self.mass_SI_to_FD

    def GetElectronMassInSIUnits(self):
        return self.electronMass_SI

    def GetElectronChargeInFDUnits(self):
        return self.electronCharge_FD

    def GetElectronChargeInSIUnits(self):
        return self.electronCharge_SI

    def ConvertSIEnergyToFDUnits(self, energy_si):
        return energy_si / self.energy_SI_to_FD

    def ConvertFDEnergyToSIUnits(self, energy_fd):
        return energy_fd * self.energy_SI_to_FD




