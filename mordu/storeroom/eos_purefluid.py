# 07/11/2025

# create presets for all the pure fluid EOS relevant to the paper and thesis
# classes
import importlib

from mordu.symbols import *

from mordu.eos import EOS

# object creation factories
def get_fluid(name:str):
    statement = f"importlib.import_module('.fluids', 'mordu.storeroom').{name}"
    return eval(statement)

def get_cp0(name:str):
    statement = f"importlib.import_module('.cp0s', 'mordu.storeroom').{name}"
    return eval(statement)

def get_alpha_r_func(name:str):
    statement = f"importlib.import_module('.alpha_r_func', 'mordu.storeroom').alpha_r_{name}"
    return eval(statement)


def create_ideal(fluid:str, cp0:str):
    from mordu.alpha_r_helmholtz import AlphaRHelmholtz

    return EOS("ideal", get_fluid(fluid) ,get_cp0(cp0), AlphaRHelmholtz, alpha_r_expr=sp.sympify(0))

def create_cubic(kind : str = "", fluid : str = "", cp0: str = ""):
    from mordu.alpha_r_cubic import AlphaRCubic

    if kind=="vdW":
        cubic_params = {
            "alpha_r_expr":  sp.log((1/rho/(1/rho-b)))-a/(R*T*1/rho),
            "a_c_expr": 27/64*R**2*T_c**2/P_c,
            "alpha_T_expr": sp.sympify(1),
            "b_expr": 1/8*R*T_c/P_c
        }

    if kind == "PR":
        cubic_params = {
            "alpha_r_expr": sp.log((1/rho/(1/rho-b))) + 1/(R*T)*a/(4*b)*2**0.5*sp.log(((1/rho-b*(2**0.5-1))/(1/rho+b*(2**0.5+1)))),
            "a_c_expr": 0.4572*R**2*T_c**2/P_c,
            "alpha_T_expr": (1+(0.37464+1.54226*omega - 0.2699*omega**2)*(1-(T/T_c)**0.5))**2,
            "b_expr": 0.0778*R*T_c/P_c
        }

    if kind == "RK":
        cubic_params = {
            "alpha_r_expr": sp.log((1/rho/(1/rho-b))) + 1/(R*T)*a/b*sp.log((1/rho/(1/rho+b))),
            "a_c_expr": 0.4275*R**2*T_c**2/P_c,
            "alpha_T_expr": 1/T**0.5,
            "b_expr":  0.0867*R*T_c/P_c
        }

    if kind == "SRK":
        cubic_params ={
            "alpha_r_expr": sp.log(((1/rho/(1/rho-b)))) + 1/(R*T)*a/b*sp.log((1/rho/(1/rho+b))),
            "a_c_expr": 0.42748*R**2*T_c**2/P_c,
            "alpha_T_expr": (1+ (0.48 + 1.574*omega - 0.176*omega**2)*(1-(T/T_c)**0.5))**2,
            "b_expr": 0.08664*R*T_c/P_c
        }

    if kind == "MSRK":
        cubic_params ={
            "alpha_r_expr": sp.log(((1/rho/(1/rho-b)))) + 1/(R*T)*a/b*sp.log((1/rho/(1/rho+b))),
            "a_c_expr": 0.42748*R**2*T_c**2/P_c,
            "alpha_T_expr": (1+ (0.48503 + 1.5571*omega - 0.15613*omega**2)*(1-(T/T_c)**0.5))**2,
            "b_expr": 0.08664*R*T_c/P_c
        }
    
    elif kind not in ["vdW", "PR", "RK", "SRK", "MSRK"]:
        raise ValueError(f"{kind} cubic equation of state has not been predefined, please create EOS from scratch")
    
    return EOS(kind, get_fluid(fluid), get_cp0(cp0), AlphaRCubic, **cubic_params)

def create_Helmholtz(kind: str = "", fluid: str="", cp0:str =""):
    from mordu.alpha_r_helmholtz import AlphaRHelmholtz

    alpha_r_expr = get_alpha_r_func(kind)()
    return EOS(kind, get_fluid(fluid), get_cp0(cp0), AlphaRHelmholtz, alpha_r_expr = alpha_r_expr)

def create_SAFT(kind:str="", fluid:str="", cp0:str =""):
    from mordu.alpha_r_saft import AlphaRSAFT

    if kind == "0323":
        SAFT_params ={
            "epsilon": 124.3255,      #epsilon/k, in (K)
            "sigma": 2.2334,          #sigma, in (A) angstrom
            "m": 2.7078,              #segment length for ammonia, dimensionless
            "epsilon_AB": 1115.64,    #epsilon^AB/k, in (K)
            "k_AB": 0.4595,           #k^AB
            "M": 2,                    #number of association sites per molecule
            "a": [
                [0.9105631445, -0.3084016918, -0.0906148351],
                [0.6361281449, 0.1860531159, 0.4527842806],
                [2.6861347891, -2.5030047259, 0.5962700728],
                [-26.547362491, 21.419793629, -1.7241829131],
                [97.759208784, -65.255885330, -4.1302112531],
                [-159.59154087, 83.318680481, 13.776631870],
                [91.29774084, -33.746922930, -8.6728470368]
            ],
            "b": [
                [0.7240946941, -0.5755498075, 0.0976883116],
                [2.2382791861, 0.6995095521, -0.2557574982],
                [-4.0025849485, 3.8925673390, -9.1558561530],
                [-21.003576815, -17.215471648, 20.642075974],
                [26.855641363, 192.67226447, -38.804430052],
                [206.55133841, -161.82646165, 93.626774077],
                [-355.60235612, -165.20769346, -29.666905585]
            ],
            "association_scheme": "2B",
        }

    elif kind == "0324":
        SAFT_params = {
            "epsilon": 204.63,      #epsilon/k, in (K)
            "sigma": 3.2386,          #sigma, in (A) angstrom
            "m": 1.1157,              #segment length for ammonia, dimensionless
            "epsilon_AB": 646.38,    #epsilon^AB/k, in (K)
            "k_AB": 0.00597,           #k^AB
            "M": 4,                   #number of association sites per molecule
            "association_scheme": "4B",
            "x_p": 1.3976,
            "mu": 1.469,
            "a": [
                [0.791982807, -0.623115538, -0.067775558],
                [1.071486513, 0.485734369, 0.028374114],
                [0.914746607, 1.124852696, 0.096122805],
                [-7.810606510, -2.094850156, 0.068150274],
                [25.78559770, 9.450498226, 0.059801866],
                [-56.98228765, -17.10272618, 0.286609791],
                [41.93089410, 7.776102807, -0.747016979]
            ],
            "b": [
                [0.791982807, -0.623115538, -0.067775558],
                [2.142973025, 0.971468739, 0.056748227],
                [2.744239820, 3.374558089, 0.288368414],
                [-31.24242604, -8.379400622, 0.272601097],
                [128.9279885, 47.25249113, 0.299009329],
                [-341.8937259, -102.6163571, 1.719658743],
                [ 293.5162587, 54.43271965, -5.229118852]
            ]
        }

    elif kind == "0328":
        epsilon = 283.18
        sigma = 3.3476
        m = 1.503
        tau = 0.74048
        u = epsilon*(1+1/T)                                 #in K
        d = sigma*(1-0.12*sp.exp(-3*epsilon/T))                #in angstrom

        eta = pi/6*rho*m*d**3                        #dimensionaless

        #Z_seg requires the universal constants D_ij available from [0336]
        D_ij = [
            [-8.8043, 4.1646270, -48.203555, 140.43620, -195.23339, 113.51500, 0, 0, 0],
            [2.9396, -6.0865383, 40.137956, -76.230797, -133.70055, 860.25349, -1535.3224, 1221.4261, -409.10539], 
            [-2.8225, 4.7600148, 11.257177, -66.382743, 69.248785, 0, 0, 0, 0],
            [0.3400, -3.1875014, 12.231796, -12.110681, 0, 0, 0, 0, 0]
        ]


        alpha_0_disp = sum([sum([D_ij[i][j]*(u/T)**(i+1)*(eta/tau)**(j+1) for j in range(0,9)]) for i in range(0,4)])

        alpha_disp = m*alpha_0_disp

        #fix the units
        alpha_disp = alpha_disp.subs([(rho, rho*N_av*1e-30)])

        SAFT_params = {
            "epsilon": epsilon,      #epsilon/k, in (K)
            "sigma": sigma,          #sigma, in (A) angstrom
            "m": m,              #segment length for ammonia, dimensionless
            "epsilon_AB": 893.1,    #epsilon^AB/k, in (K)
            "k_AB": 3.27e-2,           #k^AB
            "M": 3,                   #number of association sites per molecule
            "association_scheme": "3B",
            "alpha_disp": alpha_disp
        }

    elif kind == "0330":
        SAFT_params = {
            "epsilon": 75.092,      #epsilon/k, in (K)
            "sigma": 2.2677,          #sigma, in (A) angstrom
            "m": 2.5485,              #segment length for ammonia, dimensionless
            "epsilon_AB": 1041.5,    #epsilon^AB/k, in (K)
            "k_AB": 0.37213,           #k^AB
            "M": 4,                   #number of association sites per molecule
            "association_scheme": "4C",
            "a": [
                [0.9105631445, -0.3084016918, -0.0906148351],
                [0.6361281449, 0.1860531159, 0.4527842806],
                [2.6861347891, -2.5030047259, 0.5962700728],
                [-26.547362491, 21.419793629, -1.7241829131],
                [97.759208784, -65.255885330, -4.1302112531],
                [-159.59154087, 83.318680481, 13.776631870],
                [91.29774084, -33.746922930, -8.6728470368]
            ],
            "b": [
                [0.7240946941, -0.5755498075, 0.0976883116],
                [2.2382791861, 0.6995095521, -0.2557574982],
                [-4.0025849485, 3.8925673390, -9.1558561530],
                [-21.003576815, -17.215471648, 20.642075974],
                [26.855641363, 192.67226447, -38.804430052],
                [206.55133841, -161.82646165, 93.626774077],
                [-355.60235612, -165.20769346, -29.666905585]
            ],
        }
    
    elif kind == "0552":
        SAFT_params ={
            "epsilon": 30.32,
            "sigma":2.86,
            "m": 1.25,
            "a":[
                [0.9105631445, -0.3084016918, -0.0906148351],
                [0.6361281449, 0.1860531159, 0.4527842806],
                [2.6861347891, -2.5030047259, 0.5962700728],
                [-26.547362491, 21.419793629, -1.7241829131],
                [97.759208784, -65.255885330, -4.1302112531],
                [-159.59154087, 83.318680481, 13.776631870],
                [91.29774084, -33.746922930, -8.6728470368]
                ],
            "b":[
                [0.7240946941, -0.5755498075, 0.0976883116],
                [2.2382791861, 0.6995095521, -0.2557574982],
                [-4.0025849485, 3.8925673390, -9.1558561530],
                [-21.003576815, -17.215471648, 20.642075974],
                [26.855641363, 192.67226447, -38.804430052],
                [206.55133841, -161.82646165, 93.626774077],
                [-355.60235612, -165.20769346, -29.666905585]
                ]
        }

    else:
        raise ValueError(f"module create_SAFT has no attribute {kind}")
    return EOS(kind, get_fluid(fluid), get_cp0(cp0), AlphaRSAFT, **SAFT_params)

# lazy object dictionary
_lazy_objects = {
    # ideal gas law
    "NH3_ideal": lambda: create_ideal( "NH3" ,"NH3_cp0_NIST"),
    "H2_ideal": lambda: create_ideal( "H2", "H2_cp0_NIST"),
    
    # cubic EOS
    "NH3_vdW": lambda: create_cubic("vdW", "NH3", "NH3_cp0_NIST"),
    "NH3_PR": lambda: create_cubic("PR", "NH3", "NH3_cp0_NIST"),
    "NH3_RK": lambda: create_cubic("RK", "NH3", "NH3_cp0_NIST"),
    "NH3_SRK": lambda: create_cubic("SRK", "NH3", "NH3_cp0_NIST"),
    "NH3_MSRK": lambda: create_cubic("MSRK","NH3", "NH3_cp0_NIST"),

    "H2_vdW": lambda: create_cubic("vdW", "H2","H2_cp0_NIST"),
    "H2_PR": lambda: create_cubic("PR", "H2", "H2_cp0_NIST"),
    "H2_RK": lambda: create_cubic("RK", "H2", "H2_cp0_NIST"),
    "H2_SRK": lambda: create_cubic("SRK", "H2", "H2_cp0_NIST"),
    "H2_MSRK": lambda: create_cubic("MSRK", "H2", "H2_cp0_NIST"),

    # Helmholtz EOS
    "NH3_0290": lambda: create_Helmholtz("0290", "NH3", "NH3_cp0_NIST"),
    "NH3_0298": lambda: create_Helmholtz("0298", "NH3", "NH3_cp0_NIST"),
    "NH3_0300": lambda: create_Helmholtz("0300", "NH3", "NH3_cp0_NIST"),

    "H2_0313": lambda: create_Helmholtz("0313", "H2", "H2_cp0_NIST"),

    # SAFT EOS
    "NH3_0323": lambda: create_SAFT("0323", "NH3", "NH3_cp0_NIST"),
    "NH3_0324": lambda: create_SAFT("0324", "NH3", "NH3_cp0_NIST"),
    "NH3_0328": lambda: create_SAFT("0328", "NH3", "NH3_cp0_NIST"),
    "NH3_0330": lambda: create_SAFT("0330", "NH3", "NH3_cp0_NIST"),

}

_loaded_objects = {}

# handle the getting and the caching
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
