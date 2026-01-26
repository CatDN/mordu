#19/08/2025

from .symbols import *
from .purefluid import PureFluid

# non dimensional residual Helmholtz energy for Helmholtz EOS
class AlphaRHelmholtz():
    """The non-dimensional residual Helmholtz energy for a multiparameter equation of state.

    Multiparameter equations of state are also sometimes called Helmholtz equations of state.

    Attributes
    ----------
    expression: sp.Expr
        Expression of the residual non-dimensional Helmholtz energy , `alpha^r`, as a 
        function of only temperature in [K], `T`, and density in [mol/m3], `rho`.
    alpha_r: sp.Expr
        Expression of the residual non-dimensional Helmholtz energy , `alpha^r`, as a 
        function of only temperature in [K], `T`, and density in [mol/m3], `rho`.
    
    """
    
    def __init__(self, alpha_r_expr: sp.Expr ):
        """Create the non-dimensional residual Helmholtz energy for a multiparameter equation of state.

        Parameters
        ----------
        alpha_r_expr: sp.Expr
            Sympy expression of the residual non-dimensional Helmholyz energy, `alpha^r`,
            as a function of only temperature in [K], `T`, and density in [mol/m3], `rho`.
        """
        self.expression = alpha_r_expr  
        
        self.alpha_r = alpha_r_expr
        

    # for a pure fluid (fluid is just a placeholder argument)
    @classmethod
    def for_purefluid(cls, fluid: PureFluid, alpha_r_expr: sp.Expr):
        """Create the non-dimensional residual Helmholtz energy for a pure fluid.

        Parameters
        ----------
        fluid: PureFluid
            The `PureFluid` object for which the non-dimensional residual Helmholtz energy is being created.
            Although this parameter is not used, its presence allows for the equal usage of the `AlphaRHelmholtz` class
            compared to the other non-dimensional residual Helmholtz energy classes.
        alpha_r_expr: sp.Expr
            Sympy expression of the residual non-dimensional Helmholyz energy, `alpha^r`,
            as a function of only temperature in [K], `T`, and density in [mol/m3], `rho`.
        """
        # sub delta and tau for their pure fluid values using density and temperature
        return cls(alpha_r_expr)
