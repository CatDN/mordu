#03/03/2025
from .symbols import *   #sympy symbols
from .purefluid import PureFluid


# ideal gas isobaric heat capacity
class Cp0():
    r"""Class defining the isobaric heat capacity objects which are necessary in the creation of the
    ideal gas reduced Helmholtz energy.

    Attributes
    ---------
    fluid_formula: str
        String which defines the fluid formula. Compatible with the formula attribute of the fluid object. Examples: "NH3", "H2".
        This attribute is used to check for compatibility between the ideal and residual reduced Helmholtz energies when 
        creating the EOS.
    cp0: sp.Expr
        Sympy expression of the ideal gas isobaric heat capacity in [J/mol K] as a function of temperature, `T`, in [K].
    cp0_int_T: sp.Expr
        Sympy expression of the integral of the ideal gas heat capacity in [J/mol K] over temperature, `T`, in [K].
    cp0_over_T_int_T: sp.Expr
        Sympy expression of the integral of the ideal gas heat capacity in [J/mol K] divided by temperature, `T`, over temperature, `T`, in [K].

    Note
    ----
    The attributes have been chosen so as to facilitate the creation of the ideal gas reduced Helmholtz energy.
    This is defined as such:

    .. math:
    \alpha^0(T, P) = \frac{1}{RT}\int_{T_0}^T c_p^0dT + \frac{h_0^0}{RT}-1-\int_{T_0}^T \frac{c_p^0}{RT} dT + \ln \left(\frac{P}{P_0} \right) - \frac{s_0^0}{R}
    """

    #it is possible to define it using the expressions themselves
    def __init__(self, fluid_formula: str, cp0_expression: sp.Expr, cp0_int_T_expresion: sp.Expr = None, cp0_over_T_intT_T_expression: sp.Expr = None):
        """Definition of the ideal gas isobaric heat capacity in [J/mol K].

        Parameters
        ----------
        fluid_formula: str
            String which defines the fluid formula. Compatible with the formula attribute of the fluid object. Examples: "NH3", "H2".
            This parameter is stored as an attribute.
        cp0_expression: sp.Expr
            Sympy expression for the ideal gas isobaric heat capacity in [J/mol K] as a function of temperature, `T`, in [K].
        cp0_int_T_expression: sp.Expr, optional
            Sympy expression of the integral of the ideal gas heat capacity in [J/mol K] over temperature, `T`, in [K].
            If not provided then the integral will be calculated using `cp0_expression` and Sympy. This is not reccomended
            as it is very time consuming. It is reccomended to always provide the expression for this parameter.
        cp0_over_T_int_T: sp.Expr, optional
            Sympy expression of the integral of the ideal gas heat capacity in [J/mol K] divided by temperature, `T`, 
            over temperature, `T`, in [K].
            If not provided then the integral will be calculated using `cp0_expression` and Sympy. This is not reccomended
            as it is very time consuming. It is reccomended to always provide the expression for this parameter.

        """
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
    def from_NIST(cls, fluid: PureFluid , temperature_list: list[float, float, float], 
                  A1: float, B1: float, C1: float, D1: float, E1: float,
                  A2: float, B2: float, C2: float, D2: float, E2: float):
        """Defining the ideal gas isobaric heat capacity in [J/mol K] using the formula provided by NIST.

        The ideal gas isboraric heat capacity defined by NIST follows the following formulation: 

        .. math:: C_p^0 = A + Bt + Ct^2 + Dt^3 + E/t^2
        
        Args
        ----
        fluid: PureFluid 
            Object from PureFluid class, corresponding to the relevant fluid.
        temperature_list: list[float, float, float]
            List of temperature ranges for NIST expression.
        A1: float 
            coefficient A for the first temperature range
        B1: float  
            coefficient B for the first temperature range
        C1: float  
            coefficient C for the first temperature range
        D1: float 
            coefficient D for the first temperature range
        E1: float 
            coefficient E for the first temperature range
        A2: float 
            coefficient A for the second temperature range
        B2: float 
            coefficient B for the second temperature range
        C2: float 
            coefficient C for the second temperature range
        D2: float 
            coefficient D for the second temperature range
        E2: float 
            coefficient E for the second temperature range

        Notes
        -----
        NIST website: `https://webbook.nist.gov/`

        """

        # sigmoid, to avoid usage of heaviside
        def sigmoid(x):
            return 1/(1+sp.exp(-x))

        #create the expressions from the NIST arguments
        T0, T1, T2 = temperature_list
        A, B, C, D, E, a, b, c, d, e = A1, B1, C1, D1, E1, A2, B2, C2, D2, E2

        t = T/1000

        #if the temperature lies beyond the accpetable bounds it will still use the equation
        cp0 = (A+B*t+C*t**2 + D*t**3 + E/t**2)*sigmoid(-T+T1) + (a+b*t+c*t**2 + d*t**3 + e/t**2)*sigmoid(T-T1)
        cp0_int_T = (A*t + B*t**2/2 + C*t**3/3 +D*t**4/4 - E/t)*1000*sigmoid(-T+T1) + (a*t + b*t**2/2 + c*t**3/3 +d*t**4/4 - e/t)*1000*sigmoid(T-T1)
        cp0_over_T_int_T = (A*sp.log(T) + B*t + C/2*t**2 + D/3*t**3 - E/2/t**2)*sigmoid(-T+T1) + (a*sp.log(T) + b*t +c/2*t**2 + d/3*t**3 - e/2/t**2)*sigmoid(T-T1)

        return cls(fluid.formula, cp0, cp0_int_T, cp0_over_T_int_T)

    # or define it using the [0421] values
    @classmethod
    def from_0421(cls, fluid: PureFluid , temperature_list: list[float, float], 
                  a0: float, a1: float, a2: float, a3: float, a4: float):
        """Defining the ideal gas isobaric heat capacity according to the formulation provided by [0421].

        Args
        ----
        fluid: PureFluid
            Object from PureFluid class of the relevant fluid.
        temperature_list: list[float, float]
            Lower and upper temperature limit for [0421] expression.
        a0: float
            Coefficient a0
        a1: float
            Coefficient a1*1e3
        a2: float
            Coefficient a2*1e5
        a3: float
            Coefficient a3*1e8
        a4: float
            Coefficient a4*1e11
        
        References
        ----------
        [0421] Poling, B. E., Prausnitz, J. M., John Paul, O. C., & Reid, R. C. (2001). The properties of gases and liquids (Vol. 5). New York: Mcgraw-hill.
        """

        cp0 = (a0 + a1*1e-3*T + a2*1e-5*T**2 + a3*1e-8*T**3 + a4*1e-11*T**4)*fluid.R
        cp0_int_T = (a0*T + a1*1e-3*T**2/2 + a2*1e-5*T**3/3 + a3*1e-8*T**4/4 + a4*1e-11*T**5/5)*fluid.R
        cp0_over_T_int_T = (a0*sp.log(T) + a1*1e-3*T + a2*1e-5*T**2/2 + a3*1e-8*T**3/3 + a4*1e-11*T**4/4)*fluid.R

        return cls(fluid.formula, cp0, cp0_int_T, cp0_over_T_int_T)




