#04/03/2025
#05/03/2025

#class for mixture EOS
from symbols import *

class EOSMixture():

    def __init__(self, name, mixture, alpha_r_class, EOS_1, EOS_2, beta_T:list = None, gamma_T:list = None, beta_V: list = None, gamma_V: list = None , **kwargs):
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

        # for making the calculation of other properties easier
        fluids = [self.mixture.fluid_1, self.mixture.fluid_2]
        z = [self.mixture.z1, self.mixture.z2]

        # calculate alpha0
        self.alpha_0 = self.mixture.z1*(self.EOS_1.alpha_0.alpha_0 + sp.log(self.mixture.z1)) + self.mixture.z2*(self.EOS_2.alpha_0.alpha_0 + sp.log(self.mixture.z2))

        # calculate alpha_r
        self.alpha_r = alpha_r_class.for_mixture(self.mixture, [self.EOS_1.alpha_r, self.EOS_2.alpha_r], **kwargs)

        # generalised binary coefficients according to [0315] page 69 on pdf
        if beta_T == None or gamma_T == None or beta_V == None or gamma_V == None:
            gamma_V = [[4*(1/fluids[i].rho_c + 1/fluids[j].rho_c) / (1/fluids[i].rho_c**(1/3) + 1/fluids[j].rho_c**(1/3))**3 for i in [0,1]] for j in [0,1]]
            gamma_T = [[0.5 * (fluids[i].T_c + fluids[j].T_c) / (fluids[i].T_c*fluids[j].T_c)**0.5 for i in [0,1]] for j in [0,1]]
            beta_T = [[1, 1],[1, 1]]
            beta_V = [[1, 1],[1, 1]]
            

        # residual density
        self.rho_r = 1/sum([sum([z[i]*z[j]*beta_V[i][j]*gamma_V[i][j] * (z[i] + z[j])/ (beta_V[i][j]**2 * z[i] + z[j]) * 1/8 * (fluids[i].rho_c**(-1/3) + fluids[j].rho_c**(-1/3))**3 for i in [0,1]]) for j in [0,1]])

        # residual temperature
        self.T_r = sum([sum([z[i]*z[j]*beta_T[i][j]*gamma_T[i][j] * (z[i] + z[j])/ (beta_T[i][j]**2 * z[i] + z[j]) * (fluids[i].T_c * fluids[j].T_c)**0.5 for i in [0,1]]) for j in [0,1]])

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
    





