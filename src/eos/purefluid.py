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




#create fluids
H2 = PureFluid(
        name = "hydrogen",
        formula = "H2",
        M = 2.01588*1e-3,                 # [kg/mol], source = https://webbook.nist.gov/cgi/cbook.cgi?ID=C1333740
        P_c = 1.2964e6,                   # [Pa], source = [0313]
        T_c = 33.145,                     # [K], source = [0313]
        rho_c = 15.508*1e3,               # [mol/m3], source = [0313]
        P_t = 0.00736e6,                  # [Pa], source = [0313]
        T_t = 13.957,                     # [K], source = [0313]
        rho_t = 38.2*1e3,                 # [mol/m3], source = [0313]
        omega = -0.22,                    # [-], source = https://en.wikipedia.org/wiki/Acentric_factor
        sigma = 2.827*1e-10,              # [m], source = [0421] page  779
        epsilon = 59.7,                   # [1/K], source = [0421] page  779
        Sutherland_S = 97,                # [K], source = "Viscous Fluid Flow" page 28
        Sutherland_T0 = 273,              # [K], source = "Viscous Fluid Flow" page 28
        Sutherland_mu0 = 8.411e-6,        # [N s/m2], source = "Viscous Fluid Flow" page 28
        dipole = 0                        #it is an apolar molecule
        )

#nitrogen
N2 = PureFluid(
        name = "nitrogen",
        formula = "N2",
        M = 28.0134e-3,                   # [kg/mol], source = https://webbook.nist.gov/cgi/cbook.cgi?ID=C7727379&Mask=4
        P_c = 33.978e5,                   # [Pa], source = [0532] page 65
        T_c = 126.192,                    # [K], source = [0451] page 44
        rho_c=  11.1839*1e3,              # [mol/m3], source = [0451] page 44
        P_t =  0.012563e6,                # [Pa], source = [0532] page 62
        T_t = 63.15,                      # [K], source = [0532] page 62
        rho_t = 31.046e3,                 # [mol/m3], source = [0532] page 62
        omega = 0.040,                    # [-], source = https://en.wikipedia.org/wiki/Acentric_factor
        sigma = 3.798e-10,                # [m], source = [0421] page 780
        epsilon = 71.4,                   # [1/K], source = [0421] page 780
        Sutherland_S = 107,               # [K], source = "Viscous Fluid Flow" page 28
        Sutherland_T0 = 273,              # [K], source = "Viscous Fluid Flow" page 28
        Sutherland_mu0 = 1.663e-5,        # [N s/m2], source = "Viscous Fluid Flow" page 28
        dipole = 0                        #it is an apolar molecule
        )

NH3 = PureFluid(
        name = "ammonia",
        formula = "NH3",
        M = 17.03052*1e-3,               # [kg/mol], source = https://webbook.nist.gov/cgi/cbook.cgi?ID=C7664417
        P_c = 11.3634e6,                 # [Pa], source = [0300] table 1
        T_c = 405.56,                    # [K], source = [0300] table 1
        rho_c =  13.696e3,               # [mol/m3], source = [0300] table 1
        P_t = 6.05339e3,                 # [Pa], source = [0300] table 1
        T_t = 195.49,                    # [K], source = [0300] table 1
        rho_t = 0.0037408e3,             # [mol/m3], source = [0300] table 1 (vapor density at triple point)
        omega = 0.253,                   # [-], source = https://en.wikipedia.org/wiki/Acentric_factor
        sigma = 2.9e-10,                 # [m], source = [0421] page 780
        epsilon = 558.3,                 # [1/K], source = [0421] page 780
        Sutherland_S = 377,              # [K], source = "Viscous Fluid FLow" page 577
        Sutherland_T0 = 273,             # [K], source = "Viscous Fluid FLow" page 577
        Sutherland_mu0 = 0.96e-5,        # [N s/m2], source = "Viscous Fluid FLow" page 577
        dipole = 1.47                    # [debye], source = [0421]  page 490
        )

CH4 = PureFluid(
        name = "methane",
        formula = "CH4",
        M = 16.04246e-3,                # [kg/mol], source = [0451]
        P_c = 4.5922e6,                 # [Pa], source = [0543]
        T_c = 190.564,                  # [K], source = [0543]
        rho_c = 10.139342719e3,         # [mol/m3], source = [0451]
        P_t = 0.011696e6,               # [Pa], source = [0543]
        T_t = 90.6941,                  # [K], source = [0543]
        rho_t = None,                   # [mol/m3], source = None
        omega = 0.0110,                 # [-], source = https://www.chemeo.com/cid/27-471-9/Methane
        sigma = 3.758e-10,              # [m], source = [0421] page 779
        epsilon = 148.6,                # [1/K], source = [0421] page 779
        Sutherland_S = 198,             # [K], source = "Viscous Fluid FLow" page 577
        Sutherland_T0 = 273,            # [K], source = "Viscous Fluid FLow" page 577
        Sutherland_mu0 = 1.1996e-5,     # [N s/m2], source = "Viscous Fluid FLow" page 577
        dipole = 0                      # it is an apolar molecule
        )

#ethane
C2H6 = PureFluid(
        name = "ethane",
        formula = "C2H6",
        M =30.0690*1e-3 ,               # [kg/mol], source = https://webbook.nist.gov/cgi/cbook.cgi?ID=C74840
        P_c = 49e5,                     # [Pa], source = https://webbook.nist.gov/cgi/cbook.cgi?ID=C74840&Mask=4#Thermo-Phase
        T_c = 305.3,                    # [K], source = https://webbook.nist.gov/cgi/cbook.cgi?ID=C74840&Mask=4#Thermo-Phase
        rho_c = 6.9e3,                  # [mol/m3], source = https://webbook.nist.gov/cgi/cbook.cgi?ID=C74840&Mask=4#Thermo-Phase
        P_t = 0.000011e5,               # [Pa], source = https://webbook.nist.gov/cgi/cbook.cgi?ID=C74840&Mask=4#Thermo-Phase
        T_t = 91,                       # [K], source = https://webbook.nist.gov/cgi/cbook.cgi?ID=C74840&Mask=4#Thermo-Phase
        rho_t = None,                   # [mol/m3], source = None
        omega = 0.099,                  # [-], source = [0421] page 725
        sigma = 4.443e-10,              # [m], source = [0421] page 779
        epsilon = 215.7,                # [1/K], source = [0421] page 779
        Sutherland_S = None,            # [K], source = None
        Sutherland_T0 = None,           # [K], source = None
        Sutherland_mu0 = None,          # [N s/m2], source = None
        dipole = 0                      # it is an apolar molecule
        )

#iso octane, aka 2,2,4-trimethylpentane
C8H18 = PureFluid(
        name = "isooctane",
        formula = "C8H18",
        M = 114.2285e-3,                # [kg/mol], source = https://webbook.nist.gov/cgi/cbook.cgi?ID=C540841&Mask=FFF
        P_c = 25.7e5,                   # [Pa], source = [0533] 
        T_c = 543.8,                    # [K], source = [0533]
        rho_c =  2.14e3,                # [mol/m3], source = https://webbook.nist.gov/cgi/cbook.cgi?ID=C540841&Mask=FFF#Thermo-Phase
        P_t = 1,                        # [Pa], source = None (due to necessity for reference state object, it has been set to 1)
        T_t = 165.76,                   # [K], source = https://webbook.nist.gov/cgi/cbook.cgi?ID=C540841&Mask=FFF#Thermo-Phase
        rho_t = None,                   # [mol/m3], source = None
        omega = 0.304,                  # [-], source = [0421] page 732
        sigma = None,                   # [m], source = 
        epsilon = None,                 # [1/K], source = 
        Sutherland_S = None,            # [K], source = 
        Sutherland_T0 = None,           # [K], source = 
        Sutherland_mu0 = None,          # [N s/m2], source = 
        dipole = 0                      # it is an apolar molecule, see [0421] page 747
        )

#methanol
CH3OH = PureFluid(
        name = "methanol",
        formula = "CH3OH",
        M = 32.0419e-3,                 # [kg/mol], source = https://webbook.nist.gov/cgi/cbook.cgi?ID=C67561&Mask=4#Thermo-Phase
        P_c = 81e5,                     # [Pa], source = https://webbook.nist.gov/cgi/cbook.cgi?ID=C67561&Mask=4#Thermo-Phase
        T_c = 513,                      # [K], source = https://webbook.nist.gov/cgi/cbook.cgi?ID=C67561&Mask=4#Thermo-Phase
        rho_c = 8.51e3 ,                # [mol/m3], source = https://webbook.nist.gov/cgi/cbook.cgi?ID=C67561&Mask=4#Thermo-Phase
        P_t = 0.1835,                   # [Pa], source = [0534] table 1
        T_t = 175.5,                    # [K], source = https://webbook.nist.gov/cgi/cbook.cgi?ID=C67561&Mask=4#Thermo-Phase
        rho_t = 1.271e-7*1e3,           # [mol/m3], source = [0534] table 1, vapour
        omega = 0.565,                  # [-], source = [0421] page 723
        sigma = 3.636e-10,              # [m], source = [0421] page 779
        epsilon = 481.8,                # [1/K], source = [0421] page 779
        Sutherland_S = None,            # [K], source = None
        Sutherland_T0 = None,           # [K], source = None
        Sutherland_mu0 = None,          # [N s/m2], source = None
        dipole = 1.7                    # [debye], source = [0421] page 738
        )


#ethanol
C2H5OH = PureFluid(
        name = "ethanol",
        formula = "C2H5OH",
        M = 46.0684e-3,                 # [kg/mol], source = https://webbook.nist.gov/cgi/cbook.cgi?ID=C64175&Mask=4#Thermo-Phase
        P_c = 63e5,                     # [Pa], source = https://webbook.nist.gov/cgi/cbook.cgi?ID=C64175&Mask=4#Thermo-Phase
        T_c = 514,                      # [K], source = https://webbook.nist.gov/cgi/cbook.cgi?ID=C64175&Mask=4#Thermo-Phase
        rho_c = 6e3,                    # [mol/m3], source = https://webbook.nist.gov/cgi/cbook.cgi?ID=C64175&Mask=4#Thermo-Phase
        P_t = 7.185e-10*1e6,            # [Pa], source = [0535] table 2
        T_t = 150,                      # [K], source = https://webbook.nist.gov/cgi/cbook.cgi?ID=C64175&Mask=4#Thermo-Phase
        rho_t = None,                   # [mol/m3], source = None
        omega = 0.649,                  # [-], source = [0421] page 725
        sigma = 4.53e-10,               # [m], source = [0421] page 779
        epsilon = 362.6,                # [1/K], source = [0421] page 779
        Sutherland_S = None,            # [K], source = None
        Sutherland_T0 = None,           # [K], source = None
        Sutherland_mu0 = None,          # [N s/m2], source = None
        dipole = 1.7                    # [debye], source = [0421] page 740
        )

# n-butanol: C4H8OH
C4H9OH = PureFluid(
        name = "nbutanol",
        formula = "C4H8OH",
        M = 74.1216*1e-3,               # [kg/mol], source = https://webbook.nist.gov/cgi/cbook.cgi?ID=71-36-3&Type=IR-SPEC&Index=QUANT-IR,1
        P_c = 44.23*1e5,                # [Pa], source = [0421] page 727
        T_c = 563.05,                   # [K], source = [0421] page 727
        rho_c =  1/275*1e6,             # [mol/m3], source = [0421] page 727
        P_t = 1,                        # [Pa], source = None (due to necessity for reference state object, it has been set to 1)
        T_t = 184.54,                   # [K], source = [0537] page 58, in table 31
        rho_t = None,                   # [mol/m3], source = None
        omega = 0.59,                   # [-], source = [0421] page 727
        sigma = None,                   # [m], source = 
        epsilon = None,                 # [1/K], source = 
        Sutherland_S = None,            # [K], source = 
        Sutherland_T0 = None,           # [K], source = 
        Sutherland_mu0 = None,          # [N s/m2], source = 
        dipole = 1.8                    # [debye], source = [0421] page 742
        )

# n-heptane
C7H16 = PureFluid(
        name = "heptane",
        formula = "C7H16",
        M = 100.2019,                   # [kg/mol], source = https://webbook.nist.gov/cgi/cbook.cgi?ID=142-82-5
        P_c = 27.4*1e5,                 # [Pa], source = https://webbook.nist.gov/cgi/cbook.cgi?ID=C142825&Mask=4#Thermo-Phase
        T_c = 540,                      # [K], source = https://webbook.nist.gov/cgi/cbook.cgi?ID=C142825&Mask=4#Thermo-Phase
        rho_c =  2.35*1e3,              # [mol/m3], source = https://webbook.nist.gov/cgi/cbook.cgi?ID=C142825&Mask=4#Thermo-Phase
        P_t = 1,                        # [Pa], source = None (due to necessity for reference state object, it has been set to 1)
        T_t = 182.56,                   # [K], source = https://webbook.nist.gov/cgi/cbook.cgi?ID=C142825&Mask=4#Thermo-Phase
        rho_t = None,                   # [mol/m3], source = 
        omega = 0.35,                   # [-], source = [0421] page 731
        sigma = None,                   # [m], source = 
        epsilon = None,                 # [1/K], source = 
        Sutherland_S = None,            # [K], source = 
        Sutherland_T0 = None,           # [K], source = 
        Sutherland_mu0 = None,          # [N s/m2], source = 
        dipole = 0                      # it is an apolar molecule
        )

#water
H2O = PureFluid(
        name = "water",
        formula = "H2O",
        M = 18.0153e-3,                 # [kg/mol], source = https://webbook.nist.gov/cgi/cbook.cgi?ID=C7732185
        P_c = 220.64e5,                 # [Pa], source = https://webbook.nist.gov/cgi/cbook.cgi?ID=C7732185&Mask=4#Thermo-Phase
        T_c = 647,                      # [K], source = https://webbook.nist.gov/cgi/cbook.cgi?ID=C7732185&Mask=4#Thermo-Phase
        rho_c = 17.9e3,                 # [mol/m3], source = https://webbook.nist.gov/cgi/cbook.cgi?ID=C7732185&Mask=4#Thermo-Phase
        P_t = 0.0061e5,                 # [Pa], source = https://webbook.nist.gov/cgi/cbook.cgi?ID=C7732185&Mask=4#Thermo-Phase
        T_t = 273.16,                   # [K], source = None
        rho_t = None,                   # [mol/m3], source = None
        omega = 0.344,                  # [-], source = [0421] page 737
        sigma = 2.641e-10,              # [m], source = [0421] page 780
        epsilon = 809.1,                # [1/K], source = [0421] page 780
        Sutherland_S = 1064,            # [K], source = "Viscous Fluid Flow" page 28 (stream)
        Sutherland_T0 = 350,            # [K], source = "Viscous Fluid Flow" page 28 (steam)
        Sutherland_mu0 = 1.12e-5,       # [N s/m2], source = "Viscous Fluid Flow" page 28 (steam)
        dipole = 1.8                    # [debye], source = [0421] page 752
        )


#carbon dioxide
CO2 = PureFluid(
        name = "carbon-dioxide",
        formula = "CO2",
        M = 44.0095e-3,                 # [kg/mol], source = https://webbook.nist.gov/cgi/cbook.cgi?ID=C124389
        P_c = 73.8,                     # [Pa], source = https://webbook.nist.gov/cgi/cbook.cgi?ID=C124389&Mask=4#Thermo-Phase
        T_c = 304.18,                   # [K], source = https://webbook.nist.gov/cgi/cbook.cgi?ID=C124389&Mask=4#Thermo-Phase
        rho_c = 10.6e3,                 # [mol/m3], source = https://webbook.nist.gov/cgi/cbook.cgi?ID=C124389&Mask=4#Thermo-Phase
        P_t = 5.185e5,                  # [Pa], source = https://webbook.nist.gov/cgi/cbook.cgi?ID=C124389&Mask=4#Thermo-Phase
        T_t = 216.58,                   # [K], source = https://webbook.nist.gov/cgi/cbook.cgi?ID=C124389&Mask=4#Thermo-Phase
        rho_t = None,                   # [mol/m3], source = None
        omega = 0.225,                  # [-], source = [0421] page 724
        sigma = 3.941e-10,              # [m], source = [0421] page 779
        epsilon = 195.2,                # [1/K], source = [0421] page 779
        Sutherland_S = None,            # [K], source = None
        Sutherland_T0 = None,           # [K], source = None
        Sutherland_mu0 = None,          # [N s/m2], source = None
        dipole = 0                      # it is an apolar molecule
        )



# Uncomment the following lines to create a new fluid, and fill in the values accordingly
# fluid = PureFluid(
#         name = "",
#         formula = "",
#         M = ,                # [kg/mol], source = 
#         P_c = ,                  # [Pa], source = 
#         T_c = ,                    # [K], source = 
#         rho_c =  ,              # [mol/m3], source = 
#         P_t = ,                  # [Pa], source = 
#         T_t = ,                     # [K], source = 
#         rho_t = ,                 # [mol/m3], source = 
#         omega = ,                    # [-], source = 
#         sigma = ,              # [m], source = 
#         epsilon = ,                   # [1/K], source = 
#         Sutherland_S = ,                # [K], source = 
#         Sutherland_T0 = ,              # [K], source = 
#         Sutherland_mu0 = ,        # [N s/m2], source = 
#         dipole =                         #it is an apolar molecule
#         )