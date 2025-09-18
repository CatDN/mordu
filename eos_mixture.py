#04/03/2025
#05/03/2025

#class for mixture EOS
from symbols import *

class EOSMixture():

    def __init__(self, name, mixture, alpha_r_class, EOS_1, EOS_2, **kwargs):
        """
        **kwargs:
            mixture_rule
        """
        #check whether EOS1 uses fluid1 and EOS2 uses fluid2
        if mixture.fluid_1 != EOS_1.fluid or mixture.fluid_2 != EOS_2.fluid:
            raise AttributeError("The fluid order used in the mixture object does not match the fluid order of the EOS")

        self.mixture = mixture
        self.EOS_1 = EOS_1
        self.EOS_2 = EOS_2

        # calculate alpha0
        self.alpha_0 = self.mixture.z1*(self.EOS_1.alpha_0.alpha_0 + sp.log(self.mixture.z1)) + self.mixture.z2*(self.EOS_2.alpha_0.alpha_0 + sp.log(self.mixture.z2))

        # calculate alpha_r
        self.alpha_r = alpha_r_class.for_mixture(self.mixture, [self.EOS_1.alpha_r, self.EOS_2.alpha_r], **kwargs)



        # #residual temperature and density parameters, harcoded for now
        # beta_T_12 = 1/0.98824 # because [0301] take ammonia to be the first component while I use hydrogen
        # gamma_T_12 = 1.1266
        # beta_V_12 = 1/1.0103 # because [0301] take ammonia to be the first component while I use hydrogen
        # gamma_V_12 = 0.7298

        beta_T_12 = 0.98824 
        gamma_T_12 = 1.1266
        beta_V_12 = 1.0103 
        gamma_V_12 = 0.7298

        # beta_T_12 = 1 
        # gamma_T_12 = 0.5 * (self.mixture.fluid_1.T_c + self.mixture.fluid_2.T_c)/(self.mixture.fluid_1.T_c * self.mixture.fluid_2.T_c)**0.5
        # beta_V_12 = 1 
        # gamma_V_12 = 4 * (1/self.mixture.fluid_1.rho_c + 1/self.mixture.fluid_2.rho_c) * (1/self.mixture.fluid_1.rho_c**(1/3) + 1/self.mixture.fluid_2.rho_c**(1/3))**(-3)
        
        # residual density
        self.rho_r = 1/(self.mixture.z1**2/self.mixture.fluid_1.rho_c + self.mixture.z2**2/self.mixture.fluid_2.rho_c + 2*self.mixture.z1*self.mixture.z2*beta_V_12*gamma_V_12*(self.mixture.z1+self.mixture.z2)/(beta_V_12**2*self.mixture.z1+self.mixture.z2)*1/8*(1/self.mixture.fluid_1.rho_c**(1/3) + 1/self.mixture.fluid_2.rho_c**(1/3))**3)

        # residual temperature
        self.T_r = 1*self.mixture.z1**2*self.mixture.fluid_1.T_c + self.mixture.z2**2*self.mixture.fluid_2.T_c + 2*self.mixture.z1*self.mixture.z2*beta_T_12*gamma_T_12*(self.mixture.z1+self.mixture.z2)/(beta_T_12**2*self.mixture.z1+self.mixture.z2)*(self.mixture.fluid_1.T_c*self.mixture.fluid_2.T_c)**0.5

        #ensure alpha_r is in terms of T and rho and not tau and delta
        self.alpha_r.alpha_r = self.alpha_r.alpha_r.subs([(delta, rho/rho_r), (tau, T_r/T), (rho_r, self.rho_r), (T_r, self.T_r)])

        #total alpha
        self.alpha = self.alpha_0 + self.alpha_r.alpha_r
        
        # optional name of EOS
        self.name = name
        
        # pressure equation, included as attribute because it is frequently necessary
        self.pressure_equation = sp.lambdify((rho, P, T, z1, z2), P- self.pressure)

    

    #pressure
    @property
    def pressure(self):
        pressure = rho*R*T*(1 + rho*sp.diff(self.alpha_r.alpha_r, rho))
        
        return pressure

    #compressibility factor
    @property
    def compressibility(self):
        compressibility = 1 + rho*sp.diff(self.alpha_r.alpha_r, rho)
        
        return compressibility

    # fugacity coefficients
    @property
    def fugacity_coefficients(self):
        alpha_r = self.alpha_r.alpha_r
        Z = self.compressibility


        ln_phi_i = [
            alpha_r + sp.diff(alpha_r, z1) - z2*sp.diff(alpha_r, z2)  - z1*sp.diff(alpha_r, z1) - sp.log(Z) + Z -1,
            alpha_r + sp.diff(alpha_r, z2) - z1*sp.diff(alpha_r, z1)  - z2*sp.diff(alpha_r, z2) - sp.log(Z) + Z -1
        ]

        phi_i = [sp.exp(ln_phi_i[i]) for i in range(len(ln_phi_i))]

        return phi_i
    





