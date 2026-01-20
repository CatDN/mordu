#19/08/2025

from .symbols import *

# non dimensional residual Helmholtz energy for Helmholtz EOS
class AlphaRHelmholtz():
    
    def __init__(self, alpha_r_expr: sp.Expr ):
        self.expression = alpha_r_expr  
        
        self.alpha_r = alpha_r_expr
        

    # for a pure fluid (fluid is just a placeholder argument)
    @classmethod
    def for_purefluid(cls, fluid, alpha_r_expr: sp.Expr):
        # sub delta and tau for their pure fluid values using density and temperature
        return cls(alpha_r_expr)
