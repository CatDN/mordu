#27/02/2025
from dataclasses import dataclass

# class for storing pure fluid data
@dataclass
class PureFluid():
        """Class for creating pure fluids through their properties, all properties should be in SI units unless otherwise specified:

        name -> name of fluid
        formula -> molecular formula of fluid according to iupac
        M -> molar mass of fluid in [kg/m3] (SI)
        P_c -> pressure at critical point in [Pa]
        T_c -> temperature at critical point in [K]
        rho_c -> density at critical point in [kg/m3]
        P_t -> pressure at triple point in [Pa]
        T_t -> temperature at triple point [K]
        rho_t -> density at triple point in [kg/m3]
        omega -> accentric factor of pure fluid, see Wikipedia, also known as omega
        sigma -> Lennard-Jones characteristic length in [m], can be obtained from [0421] pages 779 onwards
        epsilon -> Lennard-Jones characteristic energy over boltzmann constant in [K], can be obtained from [0421] pages 779 onwards
        Sutherland_S -> Sutherland's constant used in the calculation of viscosity, see "Viscous Fluid Flow" for more information and numbers, in [K]
        Sutherland_T0 -> Sutherland's reference temperature used in the calculation of viscosity, in [K]
        Sutherland_mu0 -> Sutherland's reference viscosity used in the calculation of viscosity, in [Ns/m^2]
        dipole -> the dipole moment, measure of the polarity of a molecule, in [debyes]
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



