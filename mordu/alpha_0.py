#03/03/2025
from .symbols import *
from .purefluid import PureFluid
from .cp0 import Cp0

# reference state class, stored here because its only used here
class ReferenceState():
    """Definition of the reference state pressure, temperature,  enthalpy and entropy.

    Attributes
    ----------
    P0: float
        Reference state pressure in [Pa].
    T0: float
        Reference state temperature in [K].
    h0: sp.Expr
        Ideal gas enthalpy in [J/mol].
    s0: sp.Expr
        Ideal gas entropy in [J/ mol K].
    """

    def __init__(self, fluid: PureFluid, cp0_object: Cp0, T0: float, T00: float, P0: float, P00: float):
        """Creating the reference state object.

        Parameters
        ----------
        fluid: PureFluid
            PureFluid object corresponding to the relevant fluid.
        cp0_object: Cp0
            Cp0 object corresponding to the ideal gas isobaric heat capacity of the relevant fluid.
        T0: float
            Reference temperature in [K].
        T00: float
            Temperature at which the ideal gas enthalpy is 0.
        P0: float
            Reference pressure in [Pa].
        P00: float
            Pressure at which the ideal gas entropy is 0.
        """
        cp0_int_T = cp0_object.cp0_int_T
        cp0_over_T_int_T = cp0_object.cp0_over_T_int_T

        u_00 = 0
        h_00 = u_00 + fluid.R*T00 
        h0 = h_00 + cp0_int_T.subs([(T, T0)]) - cp0_int_T.subs([(T, T00)])

        #s_0 = s_00 + int(cp0_over_T(T_0, T_00)) - Rln(T_0/T_00) - Rln(rho_0/rho_00)
        #s_0 = s_00 + int(cp0_over_T(T_0, T_00)) - Rln(P_0/P_00)
        s_00 = 0
        s0 = s_00 + cp0_over_T_int_T.subs([(T, T0)]) - cp0_over_T_int_T.subs([(T, T00)]) - fluid.R*sp.log(P0/P00)

        self.P0 = P0
        self.T0 = T0
        self.h0 = h0
        self.s0 = s0

# ideal non dimensional helmholtz energy class
# every pure fluid EOS requires an object of this class for build
class Alpha0:
    """Class for the creation of the ideal gas non-dimensional Helmholtz energy.

    Attributes
    ----------
    alpha_0: sp.Expr
        SymPy expression of the ideal gas non-dimensional Helmholtz energy as 
        a function of temperature in [K], `T`, and density in [mol/m3], `rho`.
    cp0: mordu.Cp0
        Cp0 object of the ideal gas isobaric heat capacity corresponding to the relevant fluid.
               
    """

    def __init__(self, fluid: PureFluid, cp0_object: Cp0, reference_state_values: list[float, float, float, float] = None, alpha0: sp.Expr = None):
        """Creating the ideal gas non-dimensional Helmholtz energy.

        Parameters
        ----------
        fluid: PureFluid
            PureFluid object corresponding to the relevant fluid.
        cp0_object: Cp0
            Cp0 object of the ideal gas isobaric heat capacity corresponding to the relevant fluid.
        reference_state_values: list[T0, T00, P0, P00], optional
            Values used in the definition of the ReferenceState object.
        alpha_0: sp.Expr, optional
            An expression for the ideal gas non-dimensional Helmholtz energy can be directly defined
            as a SymPy expression as a function of temperature in [K], `T`, ande density in [mol/m3],
            `rho`.
        """
        #check if the fluid used in cp0 corresponds to the fluid
        if fluid.formula != cp0_object.fluid_formula:
            raise AttributeError("The fluid used in the cp0 object does not match the desired fluid")
        #handling cp0
        cp0_int_T = cp0_object.cp0_int_T
        cp0_over_T_int_T = cp0_object.cp0_over_T_int_T

        #handling reference state
        if reference_state_values == None:
            T0 = 300
            T00 = fluid.T_t
            P0 = 1e3
            P00 = fluid.P_t

        else:
            T0, T00, P0, P00 = reference_state_values

        rf = ReferenceState(fluid, cp0_object, T0, T00, P0, P00)

        #handling alpha0
        if alpha0 == None:
            alpha0 = 1/(fluid.R*T)*cp0_int_T + rf.h0/(fluid.R*T) -1 - 1/fluid.R*cp0_over_T_int_T + sp.ln(rho*T/(rf.P0/(fluid.R*rf.T0)*rf.T0)) - rf.s0/fluid.R

        self.alpha_0 = alpha0
        self.cp0 = cp0_object
        

