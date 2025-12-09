import importlib

from mordu.cp0 import Cp0
from mordu.symbols import *

# object creation factories
def get_fluid(name:str):
    statement = f"importlib.import_module('.fluids', 'mordu.storeroom').{name}"
    return eval(statement)

def create_from_NIST(fluid:str):
    if fluid == "H2":
        params = {
            "A1": 33.066178,
            "B1": -11.363417,
            "C1": 11.432816,
            "D1": -2.772874,
            "E1": -0.158558,
            "A2": 18.563083,
            "B2": 12.257357,
            "C2": -2.859786,
            "D2": 0.268238,
            "E2": 1.977990,
        }
        temp_list = [298, 1000, 2500]

    if fluid == "NH3":
        params = {
            "A1": 19.99563,
            "B1": 49.77119,
            "C1": -15.37599,
            "D1": 1.921168,
            "E1": 0.189174,
            "A2": 52.02427,
            "B2": 18.48801,
            "C2": -3.765128,
            "D2": 0.248541,
            "E2": -12.45799,
        }
        temp_list = [298, 1400, 6000]

    if fluid == "CH4":
        params = {
            "A1": -0.703029,
            "B1": 108.4773,
            "C1": -42.52157,
            "D1": 5.862788,
            "E1": 0.678565,
            "A2": 85.81217,
            "B2": 11.26467,
            "C2": -2.114146,
            "D2": 0.138190,
            "E2": -26.42221,
        }
        temp_list = [298,1300,6000]

    if fluid == "C2H6":
        params = {
            "A1": 4.10642981e+01,
            "B1": -7.84306228e+01,
            "C1": 5.31982109e+02,
            "D1": -4.54052919e+02,
            "E1": -2.38531897e-02,
            "A2": 3.41665199e+01,
            "B2": 1.31051616e+02,
            "C2": -4.56695012e+01,
            "D2": 5.66792081e+00,
            "E2": -2.82519136e+00,
        }
        temp_list = [0,500,3000]

    if fluid == "N2": 
        params = {
            "A1": 28.968641,
            "B1": 1.853978,
            "C1": -9.647459,
            "D1": 16.63537,
            "E1": 0.000117,
            "A2": 19.50583,
            "B2": 19.88705,
            "C2": -8.598535,
            "D2": 1.369784,
            "E2": 0.527601,
        }
        temp_list = [100, 500, 2000]

    if fluid == "CO2":
        params = {
            "A1": -24.99735,
            "B1": 55.18696,
            "C1": -33.69137,
            "D1": 7.948387,
            "E1": -0.136638,
            "A2": 58.16639,
            "B2": 2.720074,
            "C2": -0.492289,
            "D2": 0.038844,
            "E2": -6.447293,
        }
    temp_list = [298, 1200, 6000]

    if fluid == "CH3OH":
        params = {
            "A1": 4.50881244e+01,
            "B1": -9.27098010e+01,
            "C1": 3.84741126e+02,
            "D1": -2.80354544e+02,
            "E1": -1.86217182e-02,
            "A2": 3.31112285e+01,
            "B2": 8.30518895e+01,
            "C2": -2.77321246e+01,
            "D2": 3.28593850e+00,
            "E2": -2.23590669e+00,
        }
        temp_list = [0, 500, 3000]

    if fluid == "C2H5OH":
        params = {
            "A1": 3.34940603e+01,
            "B1": 5.13271408e+01,
            "C1": 2.56184420e+02,
            "D1": -2.23490374e+02,
            "E1": 1.43312738e-03,
            "A2": 7.71621564e+01,
            "B2": 1.03752357e+02,
            "C2": -3.41761165e+01,
            "D2": 4.04030850e+00,
            "E2": -8.01345560e+00,
        }
        temp_list = [0, 700, 3000]

    if fluid == "C8H18":
        params = {
            "A1": 1.00000000e+00,
            "B1": 1.00000000e+00,
            "C1": 1.00000000e+00,
            "D1": 1.00000000e+00,
            "E1": 1.00000000e+00,
            "A2": -2.93230829e+01,
            "B2": 8.52852773e+02,
            "C2": -4.75599823e+02,
            "D2": 1.06866507e+02,
            "E2": 2.89007558e-01,
        }
        temp_list =  [0, 0, 1400]

    if fluid == "C4H9OH":
        params = {
            "A1": 6.39046416e+01,
            "B1": -1.16366586e+02,
            "C1": 1.29251505e+03,
            "D1": -1.31289706e+03,
            "E1": -4.67791494e-02,
            "A2": 9.97872529e+01,
            "B2": 2.20211246e+02,
            "C2": -7.72465686e+01,
            "D2": 9.60003665e+00,
            "E2": -6.73354593e+00,
        }
        temp_list = [0, 500, 3000]
    
    if fluid == "C7H16":
        params = {
            "A1": 1.0,
            "B1": 1.0,
            "C1": 1.0,
            "D1": 1.0,
            "E1": 1.0,
            "A2": -56.19544584,
            "B2": 832.45291321,
            "C2": -521.66936176,
            "D2": 126.05381043,
            "E2": 1.47741399,
        }
        temp_list = [0, 0, 3000]
    
    if fluid == "H2O":
        params = {
            "A1": 30.09200,
            "B1": 6.832514,
            "C1": 6.793435,
            "D1": -2.534480,
            "E1": 0.082139,
            "A2": 41.96426,
            "B2": 8.622053,
            "C2": -1.499780,
            "D2": 0.098119,
            "E2": -11.15764,
        }
        temp_list = [500, 1700, 6000]
    
    if fluid not in ["H2", "NH3", "CH4", "C2H6", "N2", "CO2", "CH3OH", "C2H5OH", "C8H18", "C4H9OH", "C7H16", "H2O"]:
        raise ValueError(f"a cp0 object for {fluid} from NIST has not been predefined, please create it from scratch")

    return Cp0.from_NIST(get_fluid(fluid), temp_list, **params)

def create_from_0421(fluid:str):
    if fluid == "H2":
        params = {
            "a0": 2.883,
            "a1": 3.681,
            "a2": -0.772,
            "a3": 0.692,
            "a4": -0.213,
        }
        temp_list = [50, 1000]

    if fluid == "NH3":
        params = {
            "a0": 4.238,
            "a1": -4.215,
            "a2": 2.041,
            "a3": -2.126,
            "a4": 0.761,
        }
        temp_list = [50, 1000]

    if fluid == "CH4":
        params = {
            "a0": 4.568,
            "a1": -8.975,
            "a2": 3.631,
            "a3": -3.407,
            "a4": 1.091,
        }
        temp_list = [50,1000]

    if fluid == "C2H6":
        params = {
            "a0": 4.178,
            "a1": -4.427,
            "a2": 5.660,
            "a3": -6.651,
            "a4": 2.487,
        }
        temp_list = [50, 1000]

    if fluid == "N2": 
        params = {
            "a0": 3.539,
            "a1": -0.261,
            "a2": 0.007,
            "a3": 0.157,
            "a4": -0.099,
        }
        temp_list = [50, 1000]

    if fluid == "CO2":
        params = {
            "a0": 3.259,
            "a1": 1.356,
            "a2": 1.502,
            "a3": -2.374,
            "a4": 1.056,
        }
    temp_list = [50, 1000]

    if fluid == "CH3OH":
        params = {
            "a0": 4.714,
            "a1": -6.986,
            "a2": 4.211,
            "a3": -4.443,
            "a4": 1.535,
        }
        temp_list = [50, 1000]

    if fluid == "C2H5OH":
        params = {
            "a0": 4.396,
            "a1": 0.628,
            "a2": 5.546,
            "a3": -7.024,
            "a4": 2.685,
        }
        temp_list = [50, 1000]

    if fluid == "C8H18":
        params = {
            "a0": 0.384,
            "a1": 77.059,
            "a2": 0.665,
            "a3": -5.565,
            "a4": 2.619,
        }
        temp_list =  [50, 1000]

    if fluid == "C4H9OH":
        params = {
            "a0": 4.467,
            "a1": 16.395,
            "a2": 6.688,
            "a3": -9.69,
            "a4": 3.864,
        }
        temp_list = [50, 1000]
    
    if fluid == "C7H16":
        params = {
            "a0": 9.634,
            "a1": 4.156,
            "a2": 15.494,
            "a3": -20.066,
            "a4": 7.77,
        }
        temp_list =  [50, 1000]
    


    if fluid not in ["H2", "NH3", "CH4", "C2H6", "N2", "CO2", "CH3OH", "C2H5OH", "C8H18", "C4H9OH", "C7H16"]:
        raise ValueError(f"a cp0 object for {fluid} from NIST has not been predefined, please create it from scratch")

    return Cp0.from_0421(get_fluid(fluid), temp_list, **params)

def create_from_reference(fluid:str, reference:str):
    if fluid=="H2" and reference=="0313":
        fluid = get_fluid(fluid)
        u_k = [1.616, -0.4117, -0.792, 0.758, 1.217]
        v_k = [531, 751, 1989, 2484, 6859]
        c_0 = 2.5

        cp0 = fluid.R*(c_0 + sum([u_k[i]*(v_k[i]/T)**2*sp.exp(v_k[i]/T)/(sp.exp(v_k[i]/T)-1)**2 for i in range(0, len(u_k))]))
        cp0_int_T = fluid.R*(c_0*T + sum([u_k[i]*v_k[i]/(sp.exp(v_k[i]/T)-1) for i in range(0, len(u_k))]))
        cp0_over_T_int_T = fluid.R*(c_0*sp.log(T) + sum([u_k[i]*v_k[i]/T*(1/(sp.exp(v_k[i]/T)-1)+1) - u_k[i]*sp.log(abs(sp.exp(v_k[i]/T)-1)) for i in range(0, len(u_k)) ]))

    
    if fluid == "NH3" and reference=="0300":
        fluid = get_fluid(fluid)
        u_k = [2.224, 3.148, 0.9579]
        v_k = [1646, 3965, 7231]
        c_0 = 4.0

        cp0 = fluid.R*(c_0 + sum([u_k[i]*(v_k[i]/T)**2*sp.exp(v_k[i]/T)/(sp.exp(v_k[i]/T)-1)**2 for i in range(0, len(u_k))]))
        cp0_int_T = fluid.R*(c_0*T + sum([u_k[i]*v_k[i]/(sp.exp(v_k[i]/T)-1) for i in range(0, len(u_k))]))
        cp0_over_T_int_T = fluid.R*(c_0*sp.log(T) + sum([u_k[i]*v_k[i]/T*(1/(sp.exp(v_k[i]/T)-1)+1) - u_k[i]*sp.log(abs(sp.exp(v_k[i]/T)-1)) for i in range(0, len(u_k)) ]))

    if fluid == "H2" and reference=="0318":
        fluid = get_fluid(fluid)
        cp0 = (1.1230*T - 61.468*T**0.5 + 1259.3 - 10512*T**-0.5 + 31638*T**-1)  #J/mol K
        cp0_int_T = ((1.123*T**2/2 - 61.468*T**(3/2)*2/3 + 1259.3*T - 10512*T**0.5/0.5 +31638*sp.log(T)) )
        cp0_over_T_int_T = ((1.123*T - 61.468*T**0.5/0.5 + 1259.3*sp.log(T) -10512*T**-0.5/(-0.5) + 31638*T**-1/(-1)))


    return Cp0(fluid.formula, cp0, cp0_int_T, cp0_over_T_int_T)

# lazy objects
_lazy_objects = {
    # NIST objects
    "H2_cp0_NIST": lambda: create_from_NIST("H2"),
    "NH3_cp0_NIST": lambda: create_from_NIST("NH3"),
    "CH4_cp0_NIST": lambda: create_from_NIST("CH4"),
    "C2H6_cp0_NIST": lambda: create_from_NIST("C2H6"),
    "N2_cp0_NIST": lambda: create_from_NIST("N2"),
    "CO2_cp0_NIST": lambda: create_from_NIST("CO2"),
    "CH3OH_cp0_NIST": lambda: create_from_NIST("CH3OH"),
    "C2H5OH_cp0_NIST": lambda: create_from_NIST("C2H5OH"),
    "C8H18_cp0_NIST": lambda: create_from_NIST("C8H18"),
    "C4H9OH_cp0_NIST": lambda: create_from_NIST("C4H9OH"),
    "C7H16_cp0_NIST": lambda: create_from_NIST("C7H16"),
    "H2O_cp0_NIST": lambda: create_from_NIST("H2O"),

    # 0421 objects
    "H2_cp0_0421": lambda: create_from_0421("H2"),
    "NH3_cp0_0421": lambda: create_from_0421("NH3"),
    "CH4_cp0_0421": lambda: create_from_0421("CH4"),
    "C2H6_cp0_0421": lambda: create_from_0421("C2H6"),
    "N2_cp0_0421": lambda: create_from_0421("N2"),
    "CO2_cp0_0421": lambda: create_from_0421("CO2"),
    "CH3OH_cp0_0421": lambda: create_from_0421("CH3OH"),
    "C2H5OH_cp0_0421": lambda: create_from_0421("C2H5OH"),
    "C8H18_cp0_0421": lambda: create_from_0421("C8H18"),
    "C4H9OH_cp0_0421": lambda: create_from_0421("C4H9OH"),
    "C7H16_cp0_0421": lambda: create_from_0421("C7H16"),
    "H2O_cp0_0421": lambda: create_from_0421("H2O"),

    # from reference objects
    "H2_cp0_0300": lambda: create_from_reference("H2", "0313"),
    "NH3_cp0_0300": lambda: create_from_reference("NH3", "0300"),
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
