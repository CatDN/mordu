#06/05/2025

#sympy symbols
from symbols import *

# cubic non dimensional residual Helmholts energy class
# necessary to create any cubic EOS for a pure fluid
class AlphaRCubic():
    """
    Class for calculating the residual Helmholtz energy of cubic equations of state.

    This class provides methods to compute the residual Helmholtz energy as a function of temperature and density.
    The expression for the Helmholtz energy is created by substituting the symbols for values included in the input dictionary.

    Attributes:
        expression (sp.core.add.Add): The expression for the residual Helmholtz energy.
        a (float): The parameter 'a' related to the equation of state.
        b (float): The parameter 'b' related to the equation of state.
        alpha_r (sp.core.add.Add): The evaluated residual Helmholtz energy after substituting 'a' and 'b' into the expression.

    Methods:
        for_purefluid(cls, fluid, alpha_r_expr, a_c_expr, alpha_T_expr, b_expr):
            Class method to create an instance for a pure fluid using the provided expressions and fluid properties.

        for_mixture(cls, mix, mixture_rule, cubic_alpha_r_list, alpha_r_expr, **kwargs):
            Class method to create an instance for a mixture using the provided mixture rule and expressions.
    """

    def __init__(self, alpha_r_expr: sp.core.add.Add, a_value,  b_value):
        self.expression = alpha_r_expr

        self.a = a_value
        self.b = b_value

        self.a_function = sp.lambdify((T, z1, z2), self.a)
        self.b_function = sp.lambdify((z1, z2), self.b)

        self.alpha_r = self.expression.subs([(a, self.a), (b, self.b)])

    #for a pure fluid
    @classmethod
    def for_purefluid(cls, fluid, alpha_r_expr: sp.core.add.Add, a_c_expr: sp.core.mul.Mul, alpha_T_expr: sp.core.power.Pow, b_expr: sp.core.mul.Mul):
        a_c_value = a_c_expr.subs([ (T_c, fluid.T_c), (P_c, fluid.P_c)])
        alpha_T_value = alpha_T_expr.subs([(omega, fluid.omega), (T_c, fluid.T_c)])

        a_value = a_c_value*alpha_T_value

        b_value = b_expr.subs([(T_c, fluid.T_c), (P_c, fluid.P_c)])

        return cls(alpha_r_expr, a_value, b_value)

    #for a mixture
    @classmethod
    def for_mixture(cls, mix, cubic_alpha_r_list: list, mixture_rule, **kwargs):
        alpha_r_expr = cubic_alpha_r_list[0].expression
        a_mix, b_mix = mixture_rule(mix, cubic_alpha_r_list, **kwargs)

        return cls(alpha_r_expr, a_mix, b_mix)


    


