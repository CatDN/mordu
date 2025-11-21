#19/08/2025

from symbols import *

class AlphaRHelmholtz():
    
    def __init__(self, alpha_r_expr: sp.core.add.Add ):
        self.expression = alpha_r_expr  
        
        self.alpha_r = alpha_r_expr
        

    # for a pure fluid (fluid is just a placeholder argument)
    @classmethod
    def for_purefluid(cls, fluid, alpha_r_expr: sp.core.add.Add):
        # sub delta and tau for their pure fluid values using density and temperature
        return cls(alpha_r_expr)

    # for a binary mixture
    @classmethod
    def for_mixture(cls, mix, alpha_r_list: list, 
                    F:float =0, n_ij=[], t_ij=[], d_ij=[], eta_ij=[], beta_ij=[], gamma_ij=[], epsilon_ij=[], 
                    k_pol=0, k_exp=0, k_gbs=0, 
                    beta_T:list = None, gamma_T:list = None, beta_V: list = None, gamma_V: list = None ,
                    **kwargs):
        
        z = [mix.z1, mix.z2]
        fluids = [mix.fluid_1, mix.fluid_2]

        # for each pure fluid alpha_r substitute rho and T for delta and tau with the critical values of each respective fluid
        alpha_rs = [alpha_r_list[i].alpha_r.subs([(rho, delta*fluids[i].rho_c), (T, fluids[i].T_c/tau)]) for i in range(len(fluids))]

        # now each alpha_r is a function of delta and tau instead of rho and T

        if F != 0:
            # general formula from [0451]
            alpha_r_ij =    sum([n_ij[k] * delta**d_ij[k] * tau**t_ij[k] for k in range(0,k_pol)]) + \
                            sum([n_ij[k] * delta**d_ij[k] * tau**t_ij[k] * sp.exp(- eta_ij[k] * (delta - epsilon_ij[k])**2 - beta_ij[k] * (delta - gamma_ij[k]) ) for k in range(k_pol, k_pol+k_exp)]) + \
                            sum([n_ij[k] * delta**d_ij[k] * tau**t_ij[k] * sp.exp(- eta_ij[k] * (delta - epsilon_ij[k])**2 - beta_ij[k] * (tau - gamma_ij[k])**2) for k in range(k_pol+k_exp, k_exp+k_gbs)])           
            
            Delta_alpha_r = z[0]*z[1] * F * alpha_r_ij
        else:
            Delta_alpha_r = 0

        # alpha_r should be a function of only delta and tau and z1 and z2
        alpha_r  = sum([z[i]*alpha_rs[i] for i in range(len(z))]) + Delta_alpha_r


        # generalised binary coefficients according to [0315] page 69 on pdf
        if beta_T == None or gamma_T == None or beta_V == None or gamma_V == None:
            gamma_V = [[4*(1/fluids[i].rho_c + 1/fluids[j].rho_c) / (1/fluids[i].rho_c**(1/3) + 1/fluids[j].rho_c**(1/3))**3 for i in [0,1]] for j in [0,1]]
            gamma_T = [[0.5 * (fluids[i].T_c + fluids[j].T_c) / (fluids[i].T_c*fluids[j].T_c)**0.5 for i in [0,1]] for j in [0,1]]
            beta_T = [[1, 1],[1, 1]]
            beta_V = [[1, 1],[1, 1]]
            
        # residual density
        rho_r = 1/sum([sum([z[i]*z[j]*beta_V[i][j]*gamma_V[i][j] * (z[i] + z[j])/ (beta_V[i][j]**2 * z[i] + z[j]) * 1/8 * (fluids[i].rho_c**(-1/3) + fluids[j].rho_c**(-1/3))**3 for i in [0,1]]) for j in [0,1]])

        # residual temperature
        T_r = sum([sum([z[i]*z[j]*beta_T[i][j]*gamma_T[i][j] * (z[i] + z[j])/ (beta_T[i][j]**2 * z[i] + z[j]) * (fluids[i].T_c * fluids[j].T_c)**0.5 for i in [0,1]]) for j in [0,1]])

        #ensure alpha_r is in terms of T and rho and not tau and delta
        alpha_r = alpha_r.subs([(delta, rho/rho_r), (tau, T_r/T)])

        return cls(alpha_r)