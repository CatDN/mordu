#06/05/2025

#sympy symbols
from .symbols import *

#pure fluid object
from purefluid import PureFluid

# cubic non dimensional residual Helmholts energy class
# necessary to create any cubic EOS for a pure fluid
class AlphaRCubic():
    """The non-dimensional residual Helmholtz energy for a cubic equation of state

    Each equation of state uses a different definition for the residual Helmholtz energy.
    This definition works for cubic equations of state only.
    This definition **does not** include an option for volume translation.
    """

    def __init__(self, alpha_r_expr: sp.Expr, a_value: float,  b_value: float):
        """Create the non-dimensional Helmholtz energy for a cubic equation of state

        Parameters
        ----------
        alpha_r_expr: sp.Expr
            Sympy expression for the non-dimensional helmholtz energy as a function 
            of the symbols T (temperature in [K]), rho (density in [mol/m3]),
            a (a parameter in []), b (b parameter in [])
        a_value: sp.Expr
            Sympy expression of the a parameter for the specific fluid for the specific EOS
            This value can be a dependent of the symbol T (temperature in [K])
        b_value: float
            Value of the b parameter for the specific fluid for the specific EOS
            This value is almost always a float
        """
        self.expression = alpha_r_expr

        self.a = a_value
        self.b = b_value

        self.a_function = sp.lambdify((T, z1, z2), self.a)
        self.b_function = sp.lambdify((z1, z2), self.b)

        self.alpha_r = self.expression.subs([(a, self.a), (b, self.b)])

    #for a pure fluid
    @classmethod
    def for_purefluid(cls, fluid: PureFluid, alpha_r_expr: sp.Expr, a_c_expr: sp.Expr, alpha_T_expr: sp.Expr, b_expr: sp.Expr):
        """Create the non-dimensional Helmholtz energy of a cubic equation of state for a specific pure fluid
        
        Parameters
        ----------
        fluid: PureFluid
            A fluid object from the PureFluid class.
        alpha_r_expr: sp.Expr
            A sympy expression for the non-dimensional helmholtz energy as a function 
            of the symbols T (temperature in [K]), rho (density in [mol/m3]),
            a (a parameter in []), b (b parameter in []).
        a_c_expr: sp.Expr
            A sympy expression for the a_c coefficient of the specific cubic EOS.
            This expression is most often dependent on the T_c (the critical temperature
            of the fluid in [K]) and P_c (the critical pressure of the fluid in [Pa]) 
            sympy symbols.
        alpha_T_expr: sp.Expr
            A sympy expression for the alpha_T coefficient of the specific cubic EOS.
            This expression can be a function of T (temperature in [K]), T_c (the critical
            temperature of the fluid in [K]) and omega (the accentric factor of the fluid).
            In the specific case of the van-der-Waals EOS, alpha_T is a simpified 1.
        b_expr: sp.Expr
            A sympy expression for the b coefficient of the specific cubic EOS.
            This expression is most often dependent on the T_c (the critical temperature
            of the fluid in [K]) and P_c (the critical pressure of the fluid in [Pa]) 
            sympy symbols.
        """
        a_c_value = a_c_expr.subs([ (T_c, fluid.T_c), (P_c, fluid.P_c)])
        alpha_T_value = alpha_T_expr.subs([(omega, fluid.omega), (T_c, fluid.T_c)])

        a_value = a_c_value*alpha_T_value

        b_value = b_expr.subs([(T_c, fluid.T_c), (P_c, fluid.P_c)])

        return cls(alpha_r_expr, a_value, b_value)




    


