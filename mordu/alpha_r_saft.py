# 30/09/2025

from .symbols import *

# non dimensional Helmholtz energy for SAFT EOS
class AlphaRSAFT():
    
    def __init__(self, epsilon:float, sigma:float, m:float, epsilon_AB:float, k_AB:float, M:int, x_p:float, assocation_scheme:str,
                 alpha_hs, alpha_chain, alpha_disp, alpha_assoc, alpha_multipolar):

        self.epsilon = epsilon

        self.sigma = sigma
        
        self.m = m

        self.epsilon_AB = epsilon_AB

        self.k_AB = k_AB

        self.M = M

        self.x_p = x_p

        self.assocation_scheme = assocation_scheme

        self.alpha_hs = alpha_hs
        self.alpha_chain = alpha_chain
        self.alpha_disp = alpha_disp
        self.alpha_assoc = alpha_assoc
        self.alpha_multipolar = alpha_multipolar

        self.alpha_r = alpha_hs + alpha_chain + alpha_disp + alpha_assoc + alpha_multipolar

        

    # for a pure fluid (fluid is just a placeholder argument)
    @classmethod
    def for_purefluid(cls, fluid, 
                      epsilon:float = 0, sigma:float = 0, m:float = 0, 
                      epsilon_AB:float = 0, k_AB:float = 0, M:int = 0, 
                      x_p:float = 0, 
                      a:list=[], b:list=[], 
                      association_scheme:str = "", 
                      alpha_hs = None, alpha_chain=None, alpha_disp = None, alpha_assoc=None, alpha_multipolar=None,
                      **kwargs):
        
        # values used in most terms
        d = sigma * (1 - 0.12 * sp.exp(-3*epsilon/T))

        zeta_n = [pi/6*rho*m*d**n for n in range(0, 4)]

        eta = zeta_n[3]

        # print(epsilon, sigma, m, epsilon_AB, k_AB, M, a, b)
        # hard sphere
        if alpha_hs == None:
            alpha_hs = cls.alpha_hs(m, zeta_n).subs([(rho, rho*N_av*1e-30)])

        # chain
        if alpha_chain ==None:
            alpha_chain = cls.alpha_chain(m, d, zeta_n).subs([(rho, rho*N_av*1e-30)])

        # dispersion
        if alpha_disp == None:
            m2e1sigma3 = m**2*epsilon/T*sigma**3
            m2e2sigma3 = m**2*(epsilon/T)**2*sigma**3

            alpha_disp = cls.alpha_disp(m, eta, a, b, m2e1sigma3, m2e2sigma3).subs([(rho, rho*N_av*1e-30)])

        # association
        if alpha_assoc == None:
            g_hs = 1/(1-zeta_n[3]) + d/2* 3*zeta_n[2]/(1-zeta_n[3])**2 + d**2/4 * 2*zeta_n[2]**2/(1-zeta_n[3])**3

            Delta = g_hs*(sp.exp(epsilon_AB/T)-1) * (sigma**3 * k_AB)
            
            alpha_assoc = cls.alpha_assoc(association_scheme, Delta).subs([(rho, rho*N_av*1e-30)])

        # multipolar
        if alpha_multipolar ==None:
            alpha_multipolar = cls.alpha_multipolar(fluid, 1, sigma, m, x_p).subs([(rho, rho*1e-6)])

        return cls(epsilon, sigma, m, epsilon_AB, k_AB, M, x_p, association_scheme, alpha_hs, alpha_chain, alpha_disp, alpha_assoc, alpha_multipolar)

   
    ################################################# Static methods
    # hard sphere
    @staticmethod
    def alpha_hs(m, zeta_n):
        # from [0329]
        alpha_hs = m / zeta_n[0] * ( 3*zeta_n[1]*zeta_n[2]/(1-zeta_n[3]) + zeta_n[2]**3/(zeta_n[3]*(1-zeta_n[3])**2) + (zeta_n[2]**3/zeta_n[3]**2 -zeta_n[0])*sp.log(1-zeta_n[3]))
        # print(f"alpha_hs = {alpha_hs}")
        return alpha_hs
    
    # chain
    @staticmethod
    def alpha_chain(m, d, zeta_n):
        g_hs = 1/(1-zeta_n[3]) + d/2* 3*zeta_n[2]/(1-zeta_n[3])**2 + d**2/4 * 2*zeta_n[2]**2/(1-zeta_n[3])**3

        alpha_chain = -(m-1)*sp.log(g_hs)
        # print(f"alpha_chain = {alpha_chain}")
        return alpha_chain
    
    # dispersion
    @staticmethod
    def alpha_disp(m, eta, a:list, b:list, m2e1sigma3, m2e2sigma3):

        C1 = (1 + m*(8*eta - 2*eta**2)/(1-eta)**4 + (1-m)* (20*eta-27*eta**2 + 12*eta**3 -2*eta**4)/((1-eta)*(2-eta))**2)**(-1)
        # print(f"C1 = {C1}")

        a_i = [a[i][0] + (m-1)/m*a[i][1] + (m-1)/m * (m-2)/m * a[i][2] for i in range(0,7)]

        I1 = sum([a_i[i] *eta**i for i in range (0, 7)])
        # print(f"I1 = {I1}")

        b_i = [b[i][0] + (m-1)/m*b[i][1] + (m-1)/m * (m-2)/m * b[i][2] for i in range(0,7)]

        I2 = sum([b_i[i] *eta**i for i in range (0, 7)])
        # print(f"I2 = {I2}")

        alpha_disp = -2*pi*rho*I1* m2e1sigma3 - pi*rho*m*C1*I2* m2e2sigma3
        # print(f"alpha_disp = {alpha_disp}")
        return alpha_disp
    
    # association
    @staticmethod
    def alpha_assoc(association_scheme: str ="", Delta: sp.core.add.Add = 0):
        ########################################### define association scheme methods
        def assoc_2B(Delta):
            # source = [0328]
            X_A = (-1 + (1+4*rho*Delta)**0.5)/(2*rho*Delta)
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
            print("association scheme 2B was selected")

        elif association_scheme=="4B":
            scheme = assoc_4B
            print("association scheme 4B was selected")

        elif association_scheme=="3B":
            scheme = assoc_3B 
            print("association scheme 3B was selected")

        elif association_scheme=="4C":
            scheme = assoc_4C
            print("association scheme 4C was selected")

        elif association_scheme=="":   # for if the molecule does not associate, like hydrogen
            print("no association scheme was selected")
            return sp.simplify(0)

        else:
            raise(ValueError("Please select a valid association scheme..."))

        ########################################## calculate the association alpha
           
        X = scheme(Delta)
        
        alpha_assoc = sum([sp.log(X_A) - X_A/2 + 0.5 for X_A in X])
        # print(f"alpha_assoc = {alpha_assoc}")
        return alpha_assoc
    
    # multipolar
    @staticmethod
    def alpha_multipolar(fluid: object , z: float , sigma: float, m: float, x_p: float):
        if x_p == 0:
            return sp.sympify(0)

        # # new multipolar variation based on [0324], [0326] and [0333] only
        # # see logging 2025-12-20
        mu = fluid.dipole   # in Debyes

        # non dimensional variables
        rho_star = rho*N_av*(sigma*1e-8)**3                                         # non dimensional number density, [0334]

        mu_star = mu*1e-18/(k_b_Gaussian*T*(sigma*1e-8)**3)**0.5        #non dimensional dipole moment []
        T_star = 1/mu_star**2  

        J_6 = sp.exp(-0.488498 * rho_star**2 * sp.log(T_star) + 
                     0.863195*rho_star**2 + 
                     0.761344*rho_star*sp.log(T_star) +
                     -0.750086*rho_star +
                     -0.218562*sp.log(T_star) +
                     -0.538463)

        alpha_2 = -2/3 * pi * rho * N_av/(k_b_Gaussian**2 * T**2) * z**2 * x_p**2 * m**2 * (mu*1e-18)**4 / (sigma*1e-8)**3 *J_6

        K = sp.exp( -1.050534 * rho_star**2 * sp.log(T_star) + 
                    1.747476 * rho_star**2 + 
                    1.749366 * rho_star*sp.log(T_star) +
                    -1.999227 * rho_star +
                    -0.661046 *sp.log(T_star) +
                    -3.028720)

        alpha_3 = 32/135 * pi**3 * (14*pi/5)**0.5 * N_av**2 *rho**2 / (k_b_Gaussian * T)**3 * z**3 * x_p**3 * m**3 * (mu*1e-18)**6 * 1/(sigma*1e-8)**3 * K

        # A_multipolar = A_2/(1- A_3/A_2)
        alpha_multipolar = alpha_2/(1-alpha_3/alpha_2) 

        return alpha_multipolar
    
    
    
    