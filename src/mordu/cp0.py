#03/03/2025
from .symbols import *   #sympy symbols
from .purefluid import PureFluid

from .storeroom.fluids import H2, NH3, CH4, N2, CO2, CH3OH , C2H5OH, C2H6, C8H18, C4H9OH, C7H16, H2O

# ideal gas isobaric heat capacity
class Cp0():

    #it is possible to define it using the expressions themselves
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

    # or define it using the NIST values
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

    # or define it using the [0421] values
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




