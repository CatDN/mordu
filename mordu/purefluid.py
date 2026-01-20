#27/02/2025
from dataclasses import dataclass

# class for storing pure fluid data
@dataclass
class PureFluid():
        """ Class for creating a pure fluid based on its properties

        This class stores all the properties which define a fluid. These properties are then used
        when creating an equation of state for the particular fluid.

        Attributes
        ----------
        name: str
            Name of the fluid
        formula: str
			IUPAC formula of the fluid
        M: float
			Molar mass in [kg/mol]
        P_c: float
			Critical point pressure in [Pa]
		T_c: float
            Critical point temperature in [K]
        rho_c: float
            Critical point density in [mol/m3]
        P_t: float
            Triple point pressure in [Pa]
        T_t: float
            Triple point temperature in [K]
        rho_t: float
            Triple point temperature of the vapour in [mol/m3]
        omega: float
            Accentric factor in []
        sigma: float
            Lennard-Jones characteristic length in [m]
        epsilon: float
            Lennard-Jones characteristic energy over the Boltzmann constant in [K]
        Sutherland_S: float
            Sutherland's constant for the calculation of viscosity in [K]
        Sutherland_T0: float
            Sutherland's reference temperature for the calculation of viscosity in [K]
        Sutherland_mu0: float
            Sutherland's reference viscosity for the calculation of viscosity in [Ns/m2]
        dipole: float
            Dipole moment of the fluid molecule in [D] (debyes)
        """     
        name: str
        formula: str
        M: float
        P_c: float
        T_c: float
        rho_c: float
        P_t: float
        T_t: float
        rho_t: float
        omega: float
        sigma: float
        epsilon: float
        Sutherland_S: float
        Sutherland_T0: float
        Sutherland_mu0: float
        dipole: float
        R: float = 8.31446261815324



