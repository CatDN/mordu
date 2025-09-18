#27/02/2025
#EOS parent class
import sympy as sp
from alpha_0 import Alpha0

from symbols import *

class EOS:

    def __init__(self, name:str, fluid, cp0_object, AlphaRClass, **kwargs):
        self.fluid = fluid
        
        self.alpha_0 = Alpha0(self.fluid, cp0_object)
        self.alpha_r = AlphaRClass.for_purefluid(self.fluid, **kwargs)

        self.alpha = (self.alpha_0.alpha_0 + self.alpha_r.alpha_r)

        self.name = name

        self.pressure_equation = sp.lambdify((rho, P, T), P- self.pressure)

    @property
    def pressure(self):
        pressure_expression = rho*self.fluid.R*T*(1+ rho*sp.diff(self.alpha_r.alpha_r, rho))
        return pressure_expression

    @property
    def compressibility(self):
        compressibility_expression = (1+ rho*sp.diff(self.alpha_r.alpha_r, rho))
        return compressibility_expression

    @property
    def specific_heat(self):
        alpha_0 = self.alpha_0.alpha_0.subs([(T, 1/T_inv)])
        alpha_r = self.alpha_r.alpha_r.subs([(T, 1/T_inv)])

        cp_expression = self.fluid.R*(-T_inv**2*(sp.diff(alpha_0, T_inv, T_inv) + sp.diff(alpha_r, T_inv, T_inv)) + \
                        (1+rho*sp.diff(alpha_r, rho) - rho*T_inv*sp.diff(alpha_r, rho, T_inv))**2/ \
                        (1+2*rho*sp.diff(alpha_r, rho) +rho**2*sp.diff(alpha_r,rho, rho)))
        cp_expression = cp_expression.subs([(T_inv, 1/T)])

        return cp_expression

    @property
    def joule_thomson_coefficient(self):
        alpha_0 = self.alpha_0.alpha_0.subs([(T, 1/T_inv)])
        alpha_r = self.alpha_r.alpha_r.subs([(T, 1/T_inv)])

        muJT_expression = 1/(self.fluid.R*rho)*-(rho*sp.diff(alpha_r, rho) +rho**2*sp.diff(alpha_r, rho, rho)+rho*T_inv*sp.diff(alpha_r, rho, T_inv))/ \
                            ((1+rho*sp.diff(alpha_r,rho) -rho*T_inv*sp.diff(alpha_r, rho, T_inv))-T_inv**2*(sp.diff(alpha_0, T_inv, T_inv)+sp.diff(alpha_r,T_inv,T_inv))*\
                                (1+2*rho*sp.diff(alpha_r, rho)+rho**2*sp.diff(alpha_r, rho, rho)))

        muJT_expression = muJT_expression.subs([(T_inv, 1/T)])

        return muJT_expression

    @property
    def gibbs_free_energy(self):
        alpha_0 = self.alpha_0.alpha_0.subs([(T, 1/T_inv)])
        alpha_r = self.alpha_r.alpha_r.subs([(T, 1/T_inv)])

        g_expression = self.fluid.R*T*(1+ alpha_0 + alpha_r + rho*sp.diff(alpha_r, rho))

        g_expression = g_expression.subs([(T_inv, 1/T)])

        return g_expression

    #see log page 163
    @property
    def fugacity_coefficient(self):

        phi_expression = sp.exp(self.alpha_r.alpha_r + rho*sp.diff(self.alpha_r.alpha_r, rho) - sp.log(1+rho*sp.diff(self.alpha_r.alpha_r, rho)))

        return phi_expression  #fugacity coefficient phi as a function of density and temperature
    
