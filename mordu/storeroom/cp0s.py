import importlib
import sympy as sp
from mordu.cp0 import Cp0

_lazy_objects = {
    "T": lambda: importlib.import_module("mordu.symbols").T,
    "H2": lambda: importlib.import_module(".fluids", "mordu.storeroom").H2,
    "N2": lambda: importlib.import_module(".fluids", "mordu.storeroom").N2,
    "NH3": lambda: importlib.import_module(".fluids", "mordu.storeroom").NH3,
    "CH4": lambda: importlib.import_module(".fluids", "mordu.storeroom").CH4,
    "C2H6": lambda: importlib.import_module(".fluids", "mordu.storeroom").C2H6,
    "CO2": lambda: importlib.import_module(".fluids", "mordu.storeroom").CO2,
    "CH3OH": lambda: importlib.import_module(".fluids", "mordu.storeroom").CH3OH,
    "C2H5OH": lambda: importlib.import_module(".fluids", "mordu.storeroom").C2H5OH,
    "C8H18": lambda: importlib.import_module(".fluids", "mordu.storeroom").C8H18,
    "C4H9OH": lambda: importlib.import_module(".fluids", "mordu.storeroom").C4H9OH,
    "C7H16": lambda: importlib.import_module(".fluids", "mordu.storeroom").C7H16,
    "H2O": lambda: importlib.import_module(".fluids", "mordu.storeroom").H2O,

    #define all relevant cp0 objects for hydrogen or ammonia
    #HYDROGEN
    #NIST
    "NIST_H2": lambda: {"A1": 33.066178,
            "B1": -11.363417,
            "C1": 11.432816,
            "D1": -2.772874,
            "E1": -0.158558,
            "A2": 18.563083,
            "B2": 12.257357,
            "C2": -2.859786,
            "D2": 0.268238,
            "E2": 1.977990,
    },

    "H2_cp0_NIST": lambda: Cp0.from_NIST(_lazy_objects["H2"](), [298, 1000, 2500], **_lazy_objects["NIST_H2"]()),

    #0313
    "fluid_H2_0313": lambda: _lazy_objects["H2"](),
    "u_k_H2_0313": lambda: [1.616, -0.4117, -0.792, 0.758, 1.217],
    "v_k_H2_0313": lambda: [531, 751, 1989, 2484, 6859],
    "c_0_H2_0313": lambda: 2.5,

    "cp0_H2_0313_expr": lambda: _lazy_objects["fluid_H2_0313"]().R*(_lazy_objects["c_0_H2_0313"]() + sum([_lazy_objects["u_k_H2_0313"]()[i]*(_lazy_objects["v_k_H2_0313"]()[i]/_lazy_objects["T"]())**2*sp.exp(_lazy_objects["v_k_H2_0313"]()[i]/_lazy_objects["T"]())/(sp.exp(_lazy_objects["v_k_H2_0313"]()[i]/_lazy_objects["T"]())-1)**2 for i in range(0, 5)])),
    "cp0_int_T_H2_0313_expr": lambda: _lazy_objects["fluid_H2_0313"]().R*(_lazy_objects["c_0_H2_0313"]()*_lazy_objects["T"]() + sum([_lazy_objects["u_k_H2_0313"]()[i]*_lazy_objects["v_k_H2_0313"]()[i]/(sp.exp(_lazy_objects["v_k_H2_0313"]()[i]/_lazy_objects["T"]())-1) for i in range(0, len(_lazy_objects["u_k_H2_0313"]()))])),
    "cp0_over_T_int_T_H2_0313_expr": lambda: _lazy_objects["fluid_H2_0313"]().R*(_lazy_objects["c_0_H2_0313"]()*sp.log(_lazy_objects["T"]()) + sum([_lazy_objects["u_k_H2_0313"]()[i]*_lazy_objects["v_k_H2_0313"]()[i]/_lazy_objects["T"]()*(1/(sp.exp(_lazy_objects["v_k_H2_0313"]()[i]/_lazy_objects["T"]())-1)+1) - _lazy_objects["u_k_H2_0313"]()[i]*sp.log(abs(sp.exp(_lazy_objects["v_k_H2_0313"]()[i]/_lazy_objects["T"]())-1)) for i in range(0, len(_lazy_objects["u_k_H2_0313"]())) ])),
    "H2_cp0_0313": lambda: Cp0(_lazy_objects["fluid_H2_0313"]().formula, _lazy_objects["cp0_H2_0313_expr"](), _lazy_objects["cp0_int_T_H2_0313_expr"](), _lazy_objects["cp0_over_T_int_T_H2_0313_expr"]()),

    #0318
    "cp0_H2_0318_expr": lambda: (1.1230*_lazy_objects["T"]() - 61.468*_lazy_objects["T"]()**0.5 + 1259.3 - 10512*_lazy_objects["T"]()**-0.5 + 31638*_lazy_objects["T"]()**-1),  #J/mol K
    "cp0_int_T_H2_0318_expr": lambda: ((1.123*_lazy_objects["T"]()**2/2 - 61.468*_lazy_objects["T"]()**(3/2)*2/3 + 1259.3*_lazy_objects["T"]() - 10512*_lazy_objects["T"]()**0.5/0.5 +31638*sp.log(_lazy_objects["T"]())) ),
    "cp0_over_T_int_T_H2_0318_expr": lambda: ((1.123*_lazy_objects["T"]() - 61.468*_lazy_objects["T"]()**0.5/0.5 + 1259.3*sp.log(_lazy_objects["T"]()) -10512*_lazy_objects["T"]()**-0.5/(-0.5) + 31638*_lazy_objects["T"]()**-1/(-1))),
    "H2_cp0_0318": lambda: Cp0(_lazy_objects["H2"]().formula, _lazy_objects["cp0_H2_0318_expr"](), _lazy_objects["cp0_int_T_H2_0318_expr"](), _lazy_objects["cp0_over_T_int_T_H2_0318_expr"]()),

    # 0421
    "H2_0421_params": lambda: {
        "a0": 2.883,
        "a1": 3.681,
        "a2": -0.772,
        "a3": 0.692,
        "a4": -0.213,
    },

    "H2_cp0_0421": lambda: Cp0.from_0421(_lazy_objects["H2"](), [50, 1000], **_lazy_objects["H2_0421_params"]()),

    #AMMONIA
    #NIST
    "NIST_NH3": lambda: {"A1": 19.99563,
            "B1": 49.77119,
            "C1": -15.37599,
            "D1": 1.921168,
            "E1": 0.189174,
            "A2": 52.02427,
            "B2": 18.48801,
            "C2": -3.765128,
            "D2": 0.248541,
            "E2": -12.45799,
    },

    "NH3_cp0_NIST": lambda: Cp0.from_NIST(_lazy_objects["NH3"](), [298, 1400, 6000], **_lazy_objects["NIST_NH3"]()),

    #0300
    "fluid_NH3_0300": lambda: _lazy_objects["NH3"](),
    "u_k_NH3_0300": lambda: [2.224, 3.148, 0.9579],
    "v_k_NH3_0300": lambda: [1646, 3965, 7231],
    "c_0_NH3_0300": lambda: 4.0,

    "cp0_NH3_0300_expr": lambda: _lazy_objects["fluid_NH3_0300"]().R*(_lazy_objects["c_0_NH3_0300"]() + sum([_lazy_objects["u_k_NH3_0300"]()[i]*(_lazy_objects["v_k_NH3_0300"]()[i]/_lazy_objects["T"]())**2*sp.exp(_lazy_objects["v_k_NH3_0300"]()[i]/_lazy_objects["T"]())/(sp.exp(_lazy_objects["v_k_NH3_0300"]()[i]/_lazy_objects["T"]())-1)**2 for i in range(0, len(_lazy_objects["u_k_NH3_0300"]()))])),
    "cp0_int_T_NH3_0300_expr": lambda: _lazy_objects["fluid_NH3_0300"]().R*(_lazy_objects["c_0_NH3_0300"]()*_lazy_objects["T"]() + sum([_lazy_objects["u_k_NH3_0300"]()[i]*_lazy_objects["v_k_NH3_0300"]()[i]/(sp.exp(_lazy_objects["v_k_NH3_0300"]()[i]/_lazy_objects["T"]())-1) for i in range(0, len(_lazy_objects["u_k_NH3_0300"]()))])),
    "cp0_over_T_int_T_NH3_0300_expr": lambda: _lazy_objects["fluid_NH3_0300"]().R*(_lazy_objects["c_0_NH3_0300"]()*sp.log(_lazy_objects["T"]()) + sum([_lazy_objects["u_k_NH3_0300"]()[i]*_lazy_objects["v_k_NH3_0300"]()[i]/_lazy_objects["T"]()*(1/(sp.exp(_lazy_objects["v_k_NH3_0300"]()[i]/_lazy_objects["T"]())-1)+1) - _lazy_objects["u_k_NH3_0300"]()[i]*sp.log(abs(sp.exp(_lazy_objects["v_k_NH3_0300"]()[i]/_lazy_objects["T"]())-1)) for i in range(0, len(_lazy_objects["u_k_NH3_0300"]())) ])),
    "NH3_cp0_0300": lambda: Cp0(_lazy_objects["fluid_NH3_0300"]().formula, _lazy_objects["cp0_NH3_0300_expr"](), _lazy_objects["cp0_int_T_NH3_0300_expr"](), _lazy_objects["cp0_over_T_int_T_NH3_0300_expr"]()),

    # 0421
    "NH3_0421_params": lambda: {
        "a0": 4.238,
        "a1": -4.215,
        "a2": 2.041,
        "a3": -2.126,
        "a4": 0.761,
    },

    "NH3_cp0_0421": lambda: Cp0.from_0421(_lazy_objects["NH3"](), [50, 1000], **_lazy_objects["NH3_0421_params"]()),

    ####Other fluids
    #METHANE
    "NIST_CH4": lambda: {"A1": -0.703029,
            "B1": 108.4773,
            "C1": -42.52157,
            "D1": 5.862788,
            "E1": 0.678565,
            "A2": 85.81217,
            "B2": 11.26467,
            "C2": -2.114146,
            "D2": 0.138190,
            "E2": -26.42221,
    },

    "CH4_cp0_NIST": lambda: Cp0.from_NIST(_lazy_objects["CH4"](), [298,1300,6000], **_lazy_objects["NIST_CH4"]()),

    # 0421
    "CH4_0421_params": lambda: {
        "a0": 4.568,
        "a1": -8.975,
        "a2": 3.631,
        "a3": -3.407,
        "a4": 1.091,
    },

    "CH4_cp0_0421": lambda: Cp0.from_0421(_lazy_objects["CH4"](), [50, 1000], **_lazy_objects["CH4_0421_params"]()),

    #ETHANE
    "NIST_C2H6": lambda: {
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
    },

    "C2H6_cp0_NIST": lambda: Cp0.from_NIST(_lazy_objects["C2H6"](), [0,500,3000], **_lazy_objects["NIST_C2H6"]()),

    # 0421
    "C2H6_0421_params": lambda: {
        "a0": 4.178,
        "a1": -4.427,
        "a2": 5.660,
        "a3": -6.651,
        "a4": 2.487,
    },

    "C2H6_cp0_0421": lambda: Cp0.from_0421(_lazy_objects["C2H6"](), [50, 1000], **_lazy_objects["C2H6_0421_params"]()),

    #Nitrogen
    "NIST_N2": lambda: {"A1": 28.968641,
            "B1": 1.853978,
            "C1": -9.647459,
            "D1": 16.63537,
            "E1": 0.000117,
            "A2": 19.50583,
            "B2": 19.88705,
            "C2": -8.598535,
            "D2": 1.369784,
            "E2": 0.527601,
    },

    "N2_cp0_NIST": lambda: Cp0.from_NIST(_lazy_objects["N2"](), [100, 500, 2000], **_lazy_objects["NIST_N2"]()),

    # 0421
    "N2_0421_params": lambda: {
        "a0": 3.539,
        "a1": -0.261,
        "a2": 0.007,
        "a3": 0.157,
        "a4": -0.099,
    },

    "N2_cp0_0421": lambda: Cp0.from_0421(_lazy_objects["N2"](), [50, 1000], **_lazy_objects["N2_0421_params"]()),

    #Carbon dioxide
    #NIST https://webbook.nist.gov/cgi/cbook.cgi?ID=C124389&Mask=1#Thermo-Gas
    "NIST_CO2": lambda: {"A1": -24.99735,
            "B1": 55.18696,
            "C1": -33.69137,
            "D1": 7.948387,
            "E1": -0.136638,
            "A2": 58.16639,
            "B2": 2.720074,
            "C2": -0.492289,
            "D2": 0.038844,
            "E2": -6.447293,
    },

    "CO2_cp0_NIST": lambda: Cp0.from_NIST(_lazy_objects["CO2"](), [298, 1200, 6000], **_lazy_objects["NIST_CO2"]()),

    # 0421
    "CO2_0421_params": lambda: {
        "a0": 3.259,
        "a1": 1.356,
        "a2": 1.502,
        "a3": -2.374,
        "a4": 1.056,
    },

    "CO2_cp0_0421": lambda: Cp0.from_0421(_lazy_objects["CO2"](), [50, 1000], **_lazy_objects["CO2_0421_params"]()),

    #methanol
    #raw data points from https://webbook.nist.gov/cgi/cbook.cgi?ID=C67561&Mask=1#Thermo-Gas, which were then fitted to obtain the coefficients for the Shomate equation
    #see VLE.ipynb
    "NIST_CH3OH": lambda: {
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
    },
    "CH3OH_cp0_NIST": lambda: Cp0.from_NIST(_lazy_objects["CH3OH"](), [0, 500, 3000], **_lazy_objects["NIST_CH3OH"]()),

    # 0421
    "CH3OH_0421_params": lambda: {
        "a0": 4.714,
        "a1": -6.986,
        "a2": 4.211,
        "a3": -4.443,
        "a4": 1.535,
    },

    "CH3OH_cp0_0421": lambda: Cp0.from_0421(_lazy_objects["CH3OH"](), [50, 1000], **_lazy_objects["CH3OH_0421_params"]()),

    #ethanol
    "NIST_C2H5OH": lambda: {
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
    },

    "C2H5OH_cp0_NIST": lambda: Cp0.from_NIST(_lazy_objects["C2H5OH"](), [0, 700, 3000], **_lazy_objects["NIST_C2H5OH"]()),

    # 0421
    "C2H5OH_0421_params": lambda: {
        "a0": 4.396,
        "a1": 0.628,
        "a2": 5.546,
        "a3": -7.024,
        "a4": 2.685,
    },

    "C2H5OH_cp0_0421": lambda: Cp0.from_0421(_lazy_objects["C2H5OH"](), [50, 1000], **_lazy_objects["C2H5OH_0421_params"]()),

    #isooctane
    "NIST_C8H18": lambda: {
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
    },

    "C8H18_cp0_NIST": lambda: Cp0.from_NIST(_lazy_objects["C8H18"](), [0, 0, 1400], **_lazy_objects["NIST_C8H18"]()),

    # 0421
    "C8H18_0421_params": lambda: {
        "a0": 0.384,
        "a1": 77.059,
        "a2": 0.665,
        "a3": -5.565,
        "a4": 2.619,
    },

    "C8H18_cp0_0421": lambda: Cp0.from_0421(_lazy_objects["C8H18"](), [200, 1000], **_lazy_objects["C8H18_0421_params"]()),

    # n-butanol
    "NIST_C4H9OH": lambda: {
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
    },

    "C4H9OH_cp0_NIST": lambda: Cp0.from_NIST(_lazy_objects["C4H9OH"](), [0, 500, 3000], **_lazy_objects["NIST_C4H9OH"]()),

    # 0421
    "C4H9OH_0421_params": lambda: {
        "a0": 4.467,
        "a1": 16.395,
        "a2": 6.688,
        "a3": -9.69,
        "a4": 3.864,
    },

    "C4H9OH_cp0_0421": lambda: Cp0.from_0421(_lazy_objects["C4H9OH"](), [50, 1000], **_lazy_objects["C4H9OH_0421_params"]()),

    # heptane
    "NIST_C7H16": lambda: {
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
    },

    "C7H16_cp0_NIST": lambda: Cp0.from_NIST(_lazy_objects["C7H16"](), [0, 0, 3000], **_lazy_objects["NIST_C7H16"]()),

    # 0421
    "C7H16_0421_params": lambda: {
        "a0": 9.634,
        "a1": 4.156,
        "a2": 15.494,
        "a3": -20.066,
        "a4": 7.77,
    },

    "C7H16_cp0_0421": lambda: Cp0.from_0421(_lazy_objects["C7H16"](), [200, 1000], **_lazy_objects["C7H16_0421_params"]()),


    # water
    "NIST_H2O": lambda: {
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
    },

    "H2O_cp0_NIST": lambda: Cp0.from_NIST(_lazy_objects["H2O"](), [500, 1700, 6000], **_lazy_objects["NIST_H2O"]()),
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
