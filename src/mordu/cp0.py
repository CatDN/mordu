#03/03/2025
import sympy as sp
from symbols import *   #sympy symbols
from .purefluid import PureFluid

from .storeroom.fluids import H2, NH3, CH4, N2, CO2, CH3OH , C2H5OH, C2H6, C8H18, C4H9OH, C7H16, H2O

class Cp0():

    def __init__(self, fluid_formula, cp0_expression, cp0_int_T_expresion=None, cp0_over_T_intT_T_expression=None):
        self.fluid_formula = fluid_formula
        self.cp0 = cp0_expression

        #simplify and expand
        if cp0_int_T_expresion == None or cp0_over_T_intT_T_expression==None:
            cp0_terms = sp.expand(sp.simplify(self.cp0)).args

            #integral of cp0 with respect to T
            integral_T = sum([sp.integrate(term, T) for term in cp0_terms])
            self.cp0_int_T = integral_T

            #integral of cp0 over T with respect to T
            integral_T = sum([sp.integrate(term/T, T) for term in cp0_terms])
            self.cp0_over_T_int_T = integral_T

        else:
            self.cp0_int_T = cp0_int_T_expresion
            self.cp0_over_T_int_T = cp0_over_T_intT_T_expression


    @classmethod
    def from_NIST(cls, fluid: type[PureFluid] , temperature_list: list[float, float, float],**kwargs):
        """
        ...
        Args:
            fluid (PureFluid object): object from PureFluid class

            temperature_list (list[float]): list of temperature ranges for NIST expression

            **A1 (float): coefficient A for the first temperature range
            **B1 (float): coefficient B for the first temperature range
            **C1 (float): coefficient C for the first temperature range
            **D1 (float): coefficient D for the first temperature range
            **E1 (float): coefficient E for the first temperature range


            **A2 (float): coefficient A for the second temperature range
            **B2 (float): coefficient B for the second temperature range
            **C2 (float): coefficient C for the second temperature range
            **D2 (float): coefficient D for the second temperature range
            **E2 (float): coefficient E for the second temperature range

        ...
        """

        # sigmoid, to avoid usage of heaviside
        def sigmoid(x):
            return 1/(1+sp.exp(-x))

        #create the expressions from the NIST arguments
        T0, T1, T2 = temperature_list
        A, B, C, D, E, a, b, c, d, e = kwargs.values()

        t = T/1000

        #if the temperature lies beyond the accpetable bounds it will still use the equation
        cp0 = (A+B*t+C*t**2 + D*t**3 + E/t**2)*sigmoid(-T+T1) + (a+b*t+c*t**2 + d*t**3 + e/t**2)*sigmoid(T-T1)
        cp0_int_T = (A*t + B*t**2/2 + C*t**3/3 +D*t**4/4 - E/t)*1000*sigmoid(-T+T1) + (a*t + b*t**2/2 + c*t**3/3 +d*t**4/4 - e/t)*1000*sigmoid(T-T1)
        cp0_over_T_int_T = (A*sp.log(T) + B*t + C/2*t**2 + D/3*t**3 - E/2/t**2)*sigmoid(-T+T1) + (a*sp.log(T) + b*t +c/2*t**2 + d/3*t**3 - e/2/t**2)*sigmoid(T-T1)

        return cls(fluid.formula, cp0, cp0_int_T, cp0_over_T_int_T)

    @classmethod
    def from_0421(cls, fluid: type[PureFluid] , temperature_list: list[float, float],**kwargs):
        """
        ...
        Args:
            fluid (PureFluid object): object from PureFluid class

            temperature_list (list[float]): lower and upper temperature limit for [0421] expression

            **a0 (float): coefficient a0
            **a1 (float): coefficient a1*1e3
            **a2 (float): coefficient a2*1e5
            **a3 (float): coefficient a3*1e8
            **a4 (float): coefficient a4*1e11
        ...
        """
        a0, a1, a2, a3, a4 = kwargs.values()
        cp0 = (a0 + a1*1e-3*T + a2*1e-5*T**2 + a3*1e-8*T**3 + a4*1e-11*T**4)*fluid.R
        cp0_int_T = (a0*T + a1*1e-3*T**2/2 + a2*1e-5*T**3/3 + a3*1e-8*T**4/4 + a4*1e-11*T**5/5)*fluid.R
        cp0_over_T_int_T = (a0*sp.log(T) + a1*1e-3*T + a2*1e-5*T**2/2 + a3*1e-8*T**3/3 + a4*1e-11*T**4/4)*fluid.R

        return cls(fluid.formula, cp0, cp0_int_T, cp0_over_T_int_T)






#define all relevant cp0 objects for hydrogen or ammonia

#HYDROGEN
#NIST
NIST_H2 = {"A1": 33.066178,
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
H2_cp0_NIST = Cp0.from_NIST(H2, [298, 1000, 2500], **NIST_H2)

#0313
fluid = H2
u_k = [1.616, -0.4117, -0.792, 0.758, 1.217]
v_k = [531, 751, 1989, 2484, 6859]
c_0 = 2.5

cp0 = fluid.R*(c_0 + sum([u_k[i]*(v_k[i]/T)**2*sp.exp(v_k[i]/T)/(sp.exp(v_k[i]/T)-1)**2 for i in range(0, 5)]))
cp0_int_T = fluid.R*(c_0*T + sum([u_k[i]*v_k[i]/(sp.exp(v_k[i]/T)-1) for i in range(0, len(u_k))]))
cp0_over_T_int_T = fluid.R*(c_0*sp.log(T) + sum([u_k[i]*v_k[i]/T*(1/(sp.exp(v_k[i]/T)-1)+1) - u_k[i]*sp.log(abs(sp.exp(v_k[i]/T)-1)) for i in range(0, len(u_k)) ]))
H2_cp0_0313 = Cp0(fluid.formula, cp0, cp0_int_T, cp0_over_T_int_T)

#0318
cp0 = (1.1230*T - 61.468*T**0.5 + 1259.3 - 10512*T**-0.5 + 31638*T**-1)  #J/mol K
cp0_int_T = ((1.123*T**2/2 - 61.468*T**(3/2)*2/3 + 1259.3*T - 10512*T**0.5/0.5 +31638*sp.log(T)) )
cp0_over_T_int_T = ((1.123*T - 61.468*T**0.5/0.5 + 1259.3*sp.log(T) -10512*T**-0.5/(-0.5) + 31638*T**-1/(-1)))
H2_cp0_0318 = Cp0(fluid.formula, cp0, cp0_int_T, cp0_over_T_int_T)

# 0421
H2_0421 = {
    "a0": 2.883,
    "a1": 3.681,
    "a2": -0.772,
    "a3": 0.692,
    "a4": -0.213,
}

H2_cp0_0421 = Cp0.from_0421(H2, [50, 1000], **H2_0421)

#AMMONIA
#NIST
NIST_NH3 ={"A1": 19.99563,
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

NH3_cp0_NIST = Cp0.from_NIST(NH3, [298, 1400, 6000], **NIST_NH3)


#0300
fluid = NH3

u_k = [2.224, 3.148, 0.9579]
v_k = [1646, 3965, 7231]
c_0 = 4.0

cp0 = fluid.R*(c_0 + sum([u_k[i]*(v_k[i]/T)**2*sp.exp(v_k[i]/T)/(sp.exp(v_k[i]/T)-1)**2 for i in range(0, len(u_k))]))
cp0_int_T = fluid.R*(c_0*T + sum([u_k[i]*v_k[i]/(sp.exp(v_k[i]/T)-1) for i in range(0, len(u_k))]))
cp0_over_T_int_T = fluid.R*(c_0*sp.log(T) + sum([u_k[i]*v_k[i]/T*(1/(sp.exp(v_k[i]/T)-1)+1) - u_k[i]*sp.log(abs(sp.exp(v_k[i]/T)-1)) for i in range(0, len(u_k)) ]))          
NH3_cp0_0300 = Cp0(fluid.formula, cp0, cp0_int_T, cp0_over_T_int_T)

# #0495


# 0421
NH3_0421 = {
    "a0": 4.238,
    "a1": -4.215,
    "a2": 2.041,
    "a3": -2.126,
    "a4": 0.761,
}

NH3_cp0_0421 = Cp0.from_0421(NH3, [50, 1000], **NH3_0421)


####Other fluids
#METHANE
NIST_CH4 ={"A1": -0.703029,
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

CH4_cp0_NIST = Cp0.from_NIST(CH4, [298,1300,6000], **NIST_CH4)

# 0421
CH4_0421 ={
    "a0": 4.568,
    "a1": -8.975,
    "a2": 3.631,
    "a3": -3.407,
    "a4": 1.091,
}

CH4_cp0_0421 = Cp0.from_0421(CH4, [50, 1000], **CH4_0421)

#ETHANE
NIST_C2H6 = {
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

C2H6_cp0_NIST = Cp0.from_NIST(C2H6, [0,500,3000], **NIST_C2H6)

# 0421
C2H6_0421 ={
    "a0": 4.178,
    "a1": -4.427,
    "a2": 5.660,
    "a3": -6.651,
    "a4": 2.487,
}

C2H6_cp0_0421 = Cp0.from_0421(C2H6, [50, 1000], **C2H6_0421)


#Nitrogen
NIST_N2 = {"A1": 28.968641,
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

N2_cp0_NIST= Cp0.from_NIST(N2, [100, 500, 2000], **NIST_N2)

# 0421
N2_0421 ={
    "a0": 3.539,
    "a1": -0.261,
    "a2": 0.007,
    "a3": 0.157,
    "a4": -0.099,
}

N2_cp0_0421 = Cp0.from_0421(N2, [50, 1000], **N2_0421)

#Carbon dioxide
#NIST https://webbook.nist.gov/cgi/cbook.cgi?ID=C124389&Mask=1#Thermo-Gas
NIST_CO2 = {"A1": -24.99735,
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

CO2_cp0_NIST= Cp0.from_NIST(CO2, [298, 1200, 6000], **NIST_CO2)

# 0421
CO2_0421 ={
    "a0": 3.259,
    "a1": 1.356,
    "a2": 1.502,
    "a3": -2.374,
    "a4": 1.056,
}

CO2_cp0_0421 = Cp0.from_0421(CO2, [50, 1000], **CO2_0421)

#methanol
#raw data points from https://webbook.nist.gov/cgi/cbook.cgi?ID=C67561&Mask=1#Thermo-Gas, which were then fitted to obtain the coefficients for the Shomate equation
#see VLE.ipynb
NIST_CH3OH = {
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
CH3OH_cp0_NIST= Cp0.from_NIST(CH3OH, [0, 500, 3000], **NIST_CH3OH)

# 0421
CH3OH_0421 ={
    "a0": 4.714,
    "a1": -6.986,
    "a2": 4.211,
    "a3": -4.443,
    "a4": 1.535,
}

CH3OH_cp0_0421 = Cp0.from_0421(CH3OH, [50, 1000], **CH3OH_0421)

#ethanol
NIST_C2H5OH = {
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

C2H5OH_cp0_NIST= Cp0.from_NIST(C2H5OH, [0, 700, 3000], **NIST_C2H5OH)

# 0421
C2H5OH_0421 ={
    "a0": 4.396,
    "a1": 0.628,
    "a2": 5.546,
    "a3": -7.024,
    "a4": 2.685,
}

C2H5OH_cp0_0421 = Cp0.from_0421(C2H5OH, [50, 1000], **C2H5OH_0421)

#isooctane
NIST_C8H18 = {
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



C8H18_cp0_NIST= Cp0.from_NIST(C8H18, [0, 0, 1400], **NIST_C8H18)


# 0421
C8H18_0421 ={
    "a0": 0.384,
    "a1": 77.059,
    "a2": 0.665,
    "a3": -5.565,
    "a4": 2.619,
}

C8H18_cp0_0421 = Cp0.from_0421(C8H18, [200, 1000], **C8H18_0421)

# n-butanol
NIST_C4H9OH = {
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

C4H9OH_cp0_NIST= Cp0.from_NIST(C4H9OH, [0, 500, 3000], **NIST_C4H9OH)

# 0421
C4H9OH_0421 ={
    "a0": 4.467,
    "a1": 16.395,
    "a2": 6.688,
    "a3": -9.69,
    "a4": 3.864,
}

C4H9OH_cp0_0421 = Cp0.from_0421(C4H9OH, [50, 1000], **C4H9OH_0421)

# heptane
NIST_C7H16 = {
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

C7H16_cp0_NIST= Cp0.from_NIST(C7H16, [0, 0, 3000], **NIST_C7H16)

# 0421
C7H16_0421 ={
    "a0": 9.634,
    "a1": 4.156,
    "a2": 15.494,
    "a3": -20.066,
    "a4": 7.77,
}

C7H16_cp0_0421 = Cp0.from_0421(C7H16, [200, 1000], **C7H16_0421)


# water
# heptane
NIST_H20 = {
    "A1": 30.09200,
    "B1": 6.832514	,
    "C1": 6.793435,
    "D1": -2.534480,
    "E1": 0.082139	,
    "A2": 41.96426,
    "B2": 8.622053,
    "C2": -1.499780,
    "D2": 0.098119,
    "E2": -11.15764,
}

H2O_cp0_NIST= Cp0.from_NIST(H2O, [500, 1700, 6000], **NIST_H20)
