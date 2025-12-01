#03/03/2025
from .symbols import *

# reference state class, stored here because its only used here
class ReferenceState():

    def __init__(self, fluid, cp0_object, T0, T00, P0, P00):
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

    def __init__(self, fluid, cp0_object, reference = None, alpha0 = None):
        #check if the fluid used in cp0 corresponds to the fluid
        if fluid.formula != cp0_object.fluid_formula:
            raise AttributeError("The fluid used in the cp0 object does not match the desired fluid")
        #handling cp0
        cp0_int_T = cp0_object.cp0_int_T
        cp0_over_T_int_T = cp0_object.cp0_over_T_int_T

        #handling reference state
        if reference == None:
            T0 = 300
            T00 = fluid.T_t
            P0 = 1e3
            P00 = fluid.P_t

        else:
            T0, T00, P0, P00 = reference

        rf = ReferenceState(fluid, cp0_object, T0, T00, P0, P00)

        #handling alpha0
        if alpha0 == None:
            alpha0 = 1/(fluid.R*T)*cp0_int_T + rf.h0/(fluid.R*T) -1 - 1/fluid.R*cp0_over_T_int_T + sp.ln(rho*T/(rf.P0/(fluid.R*rf.T0)*rf.T0)) - rf.s0/fluid.R

        self.alpha_0 = alpha0
        self.cp0 = cp0_object
        

