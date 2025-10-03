# 30/09/2025

from symbols import *

class AlphaRSAFT():
    
    def __init__(self):

        pass
        

    # for a pure fluid (fluid is just a placeholder argument)
    @classmethod
    def for_purefluid(cls, fluid, 
                      epsilon:float = 0, sigma:float = 0, m:float = 0, epsilon_AB:float = 0, k_AB:float = 0, M:int = 0, x_p:float =0, 
                      a:list=[], b:list=[], 
                      association_scheme:str = "", **kwargs):
        
        # hard sphere
        alpha_hs = cls.alpha_hs(sigma, epsilon)

        # chain
        alpha_chain = cls.alpha_chain(sigma, epsilon, m)

        # hard chain
        alpha_hc = m*alpha_hs + alpha_chain

        # dispersion
        alpha_disp = cls.alpha_disp(sigma, epsilon, m, a, b)


        return 

    # for a binary mixture
    @classmethod
    def for_mixture(cls):


        return
    


    ################################################# Static methods
    # hard sphere
    @staticmethod
    def alpha_hs(sigma, epsilon, m):
        d = sigma * (1 - 0.12 * sp.exp(-3*epsilon/T))

        zeta_n = [sp.pi/6*rho*m*d**n for n in range(0, 4)]

        alpha_hs = m * 1/zeta_n[0] * ( 3*zeta_n[1]*zeta_n[2]/(1-zeta_n[3]) + zeta_n[2]**3/(zeta_n[3]*(1-zeta_n[3])**2) + (zeta_n[2]**3/zeta_n[3]**2 -zeta_n[0])*sp.log(1-zeta_n[3]))
        
        return alpha_hs
    
    # chain
    @staticmethod
    def alpha_chain(sigma, epsilon, m):
        d = sigma * (1 - 0.12 * sp.exp(-3*epsilon/T))

        zeta_n = [sp.pi/6*rho*m*d**n for n in range(0, 4)]

        g_hs = 1/(1-zeta_n[3]) + d/2* 3*zeta_n[2]/(1-zeta_n[3])**2 + d**2/4 * 2*zeta_n[2]**2/(1-zeta_n[3])**3

        alpha_chain = -(m-1)*sp.log(g_hs)
        
        return alpha_chain
    
    # dispersion
    @staticmethod
    def alpha_disp(sigma, epsilon, m, a:list, b:list):
        d = sigma * (1 - 0.12 * sp.exp(-3*epsilon/T))

        eta = sp.pi/6*rho*m*d**3

        C1 = (1 + m*(8*eta - 2*eta**2)/(1-eta)**4 + (1-m)* (20*eta-27*eta**2 + 12*eta**3 -2*eta**4)/((1-eta)*(2-eta))**2)**(-1)

        a_i = [a[i][0] + (m-1)/m*a[i][1] + (m-1)/m * (m-2)/m * a[i][2] for i in range(0,6+1)]

        I1 = sum([a_i[i] *eta**i for i in range (0, 6)])

        b_i = [b[i][0] + (m-1)/m*b[i][1] + (m-1)/m * (m-2)/m * b[i][2] for i in range(0,6+1)]

        I2 = sum([b_i[i] *eta**i for i in range (0, 6+1)])

        alpha_disp = -2*sp.pi*rho*I1* m**2*epsilon*1/T*sigma**3 - sp.pi*rho*m*C1*I2* m**2*(epsilon*1/T)**2*sigma**3
        
        return alpha_disp
    
    # association
    @staticmethod
    def alpha_assoc(association_scheme:str ="", epsilon:float = 0, sigma:float = 0, m:float = 0, epsilon_AB:float = 0, k_AB:float = 0, M:int = 0):
        ########################################### define association scheme methods
        def assoc_2B(Delta):
            # source = [0328]
            X_A = (-1 + (1+4*rho*Delta)**0.5)/(4*rho*Delta)
            X_B = X_A
            return [X_A, X_B]

        def assoc_4B(Delta):
            # source = [0328]
            X_A = (-(1 - 2 * rho * Delta) + ((1 + 2 * rho * Delta)**2 + 4 * rho *Delta)**0.5)/(6 * rho * Delta)
            X_B = X_A
            X_C = X_A
            X_D = 3*X_A -2
            return [X_A, X_B, X_C, X_D]            

        def assoc_3B(Delta):
            # source = [0328]
            X_A = (-(1 - rho * Delta) + ((1 + rho * Delta)**2 + 4 * rho *Delta)**0.5)/(4 * rho * Delta)
            X_B = X_A
            X_C = 2 * X_A -1
            return [X_A, X_B, X_C]            

        def assoc_4C(Delta):
            # source = [0328]
            X_A = (-1 + (1 + 8 * rho * Delta)**0.5)/(4 * rho * Delta)
            X_B = X_A
            X_C = X_A
            X_D = X_A
            return [X_A, X_B, X_C, X_D]            


        ########################################## select association scheme
        if association_scheme=="2B":
            scheme = assoc_2B

        elif association_scheme=="4B":
            scheme = assoc_4B

        elif association_scheme=="3B":
            scheme = assoc_3B 

        elif association_scheme=="4C":
            scheme = assoc_4C

        else:
            raise(ValueError("Please select a valid association scheme..."))

        ########################################## calculate the association alpha

        g_hs = 1/(1-zeta_n[3])

        Delta = g_hs*(sp.exp(epsilon_AB/T)-1) * (sigma**3 * k_AB)
        
        X = scheme(Delta)
        
        alpha_assoc = sum([sp.log(X_A) - X_A/2 for X_A in X]) + 0.5*M

        return alpha_assoc
    
    # multipolar
    def alpha_multipolar():
        return
    
    