import importlib
from mordu.purefluid import PureFluid

_lazy_objects = {
    "H2": lambda: PureFluid(
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
        ),

    "N2": lambda: PureFluid(
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
        ),

    "NH3": lambda: PureFluid(
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
        ),

    "CH4": lambda: PureFluid(
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
        ),

    "C2H6": lambda: PureFluid(
        name = "ethane",
        formula = "C2H6",
        M =30.0690*1e-3 ,
        P_c = 49e5,
        T_c = 305.3,
        rho_c = 6.9e3,
        P_t = 0.000011e5,
        T_t = 91,
        rho_t = None,
        omega = 0.099,
        sigma = 4.443e-10,
        epsilon = 215.7,
        Sutherland_S = None,
        Sutherland_T0 = None,
        Sutherland_mu0 = None,
        dipole = 0
        ),

    "C8H18": lambda: PureFluid(
        name = "isooctane",
        formula = "C8H18",
        M = 114.2285e-3,
        P_c = 25.7e5,
        T_c = 543.8,
        rho_c =  2.14e3,
        P_t = 1,
        T_t = 165.76,
        rho_t = None,
        omega = 0.304,
        sigma = None,
        epsilon = None,
        Sutherland_S = None,
        Sutherland_T0 = None,
        Sutherland_mu0 = None,
        dipole = 0
        ),

    "CH3OH": lambda: PureFluid(
        name = "methanol",
        formula = "CH3OH",
        M = 32.0419e-3,
        P_c = 81e5,
        T_c = 513,
        rho_c = 8.51e3 ,
        P_t = 0.1835,
        T_t = 175.5,
        rho_t = 1.271e-7*1e3,
        omega = 0.565,
        sigma = 3.636e-10,
        epsilon = 481.8,
        Sutherland_S = None,
        Sutherland_T0 = None,
        Sutherland_mu0 = None,
        dipole = 1.7
        ),

    "C2H5OH": lambda: PureFluid(
        name = "ethanol",
        formula = "C2H5OH",
        M = 46.0684e-3,
        P_c = 63e5,
        T_c = 514,
        rho_c = 6e3,
        P_t = 7.185e-10*1e6,
        T_t = 150,
        rho_t = None,
        omega = 0.649,
        sigma = 4.53e-10,
        epsilon = 362.6,
        Sutherland_S = None,
        Sutherland_T0 = None,
        Sutherland_mu0 = None,
        dipole = 1.7
        ),

    "C4H9OH": lambda: PureFluid(
        name = "nbutanol",
        formula = "C4H8OH",
        M = 74.1216*1e-3,
        P_c = 44.23*1e5,
        T_c = 563.05,
        rho_c =  1/275*1e6,
        P_t = 1,
        T_t = 184.54,
        rho_t = None,
        omega = 0.59,
        sigma = None,
        epsilon = None,
        Sutherland_S = None,
        Sutherland_T0 = None,
        Sutherland_mu0 = None,
        dipole = 1.8
        ),

    "C7H16": lambda: PureFluid(
        name = "heptane",
        formula = "C7H16",
        M = 100.2019,
        P_c = 27.4*1e5,
        T_c = 540,
        rho_c =  2.35*1e3,
        P_t = 1,
        T_t = 182.56,
        rho_t = None,
        omega = 0.35,
        sigma = None,
        epsilon = None,
        Sutherland_S = None,
        Sutherland_T0 = None,
        Sutherland_mu0 = None,
        dipole = 0
        ),

    "H2O": lambda: PureFluid(
        name = "water",
        formula = "H2O",
        M = 18.0153e-3,
        P_c = 220.64e5,
        T_c = 647,
        rho_c = 17.9e3,
        P_t = 0.0061e5,
        T_t = 273.16,
        rho_t = None,
        omega = 0.344,
        sigma = 2.641e-10,
        epsilon = 809.1,
        Sutherland_S = 1064,
        Sutherland_T0 = 350,
        Sutherland_mu0 = 1.12e-5,
        dipole = 1.8
        ),

    "CO2": lambda: PureFluid(
        name = "carbon-dioxide",
        formula = "CO2",
        M = 44.0095e-3,
        P_c = 73.8,
        T_c = 304.18,
        rho_c = 10.6e3,
        P_t = 5.185e5,
        T_t = 216.58,
        rho_t = None,
        omega = 0.225,
        sigma = 3.941e-10,
        epsilon = 195.2,
        Sutherland_S = None,
        Sutherland_T0 = None,
        Sutherland_mu0 = None,
        dipole = 0
        ),
}

_loaded_objects = {}

def __getattr__(name: str):
    if name in _loaded_objects:
        return _loaded_objects[name]
    if name in _lazy_objects:
        obj = _lazy_objects[name]()
        _loaded_objects[name] = obj
        return obj
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")

def __dir__():
    return sorted(list(globals().keys()) + list(_lazy_objects.keys()))
