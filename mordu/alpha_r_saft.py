# 30/09/2025

from .symbols import *
from .purefluid import PureFluid

# non dimensional Helmholtz energy for SAFT EOS
class AlphaRSAFT():
    """Non-dimensional residual Helmholtz energy for SAFT (statistically associating fluid theory)
    equations of state.

    Attributes
    ----------
    epsilon: float
        Segment dispersion energy divided by the Boltzmann constant in [K].
    sigma: float
        Segment diameter in [A] (angstroms).
    m: float
        Number of segments (dimensionless).
    epsilon_AB: float
        Energy of association divided by the Boltzmann constant in [K].
    k_AB: float
        Volume of association between sites A and B (dimensionless).
    M: int
        Number of associations per segment. This number will correspond to the first character
        of the association_scheme string.
    x_p: float
        Fraction of dipolar segments, dependent on chain length `m`. Some sources will define `x_p^{\mu}`
        as the fraction of dipolar segments while other will define `x_p^{\mu}` as the number of polar
        spherical segments in the molecule which is independent of chain length.
    association_scheme: str
        String representing the association scheme ofthe SAFT EOS, where the first character is an
        integer representing the number of associations per molecule and the second character is
        a letter characterising the association strengths of each association [0328].
    alpha_hs: sp.Expr
        SymPy expression for the hard-sphere term of the SAFT equation as a function of temperature
        in [K], `T`, and density in [mol/m3], `rho`. 
    alpha_chain: sp.Expr
        SymPy expression for the chain term of the SAFT equation as a function of temperature
        in [K], `T`, and density in [mol/m3], `rho`. 
    alpha_disp: sp.Expr
        SymPy expression for the dispersion term of the SAFT equation as a function of temperature
        in [K], `T`, and density in [mol/m3], `rho`. 
    alpha_assoc: sp.Expr
        SymPy expression for the association term of the SAFT equation as a function of temperature
        in [K], `T`, and density in [mol/m3], `rho`. If the EOS does not have a multipolar term then
        alpha_assoc will be a sympyfied 0.
    alpha_multipolar: sp.Expr
        SymPy expression for the multipolar term of the SAFT equation as a function of temperature
        in [K], `T`, and density in [mol/m3], `rho`. If the EOS does not have a multipolar term then
        alpha_multipolar will be a sympyfied 0.
    alpha_r: sp.Expr
        SymPy expression for the non-dimensional Helmholtz energy as a function of temperature in [K],
        `T`, and density in [mol/m3], `rho`.

    References
    ----------
    [0328] Huang, Stanley H., and Maciej Radosz. "Equation of state for small, large, polydisperse, and associating molecules." Industrial & Engineering Chemistry Research 29, no. 11 (1990): 2284-2294.

    """
    
    def __init__(self, epsilon:float, sigma:float, m:float, epsilon_AB:float, k_AB:float, M:int, x_p:float, assocation_scheme:str,
                 alpha_hs: sp.Expr, alpha_chain: sp.Expr, alpha_disp: sp.Expr, alpha_assoc: sp.Expr, alpha_multipolar: sp.Expr):
        """Creating the non-dimensional residual Helmholtz energy for a SAFT EOS.

        Parameters
        ----------
        epsilon: float
            Segment dispersion energy divided by the Boltzmann constant in [K].
        sigma: float
            Segment diameter in [A] (angstroms).
        m: float
            Number of segments (dimensionless).
        epsilon_AB: float
            Energy of association divided by the Boltzmann constant in [K].
        k_AB: float
            Volume of association between sites A and B (dimensionless).
        M: int
            Number of associations per segment. This number will correspond to the first character
            of the association_scheme string.
        x_p: float
            Fraction of dipolar segments, dependent on chain length `m`. Some sources will define `x_p^{\mu}`
            as the fraction of dipolar segments while other will define `x_p^{\mu}` as the number of polar
            spherical segments in the molecule which is independent of chain length.
        association_scheme: str
            String representing the association scheme ofthe SAFT EOS, where the first character is an
            integer representing the number of associations per molecule and the second character is
            a letter characterising the association strengths of each association [0328].
        alpha_hs: sp.Expr
            SymPy expression for the hard-sphere term of the SAFT equation as a function of temperature
            in [K], `T`, and density in [mol/m3], `rho`. 
        alpha_chain: sp.Expr
            SymPy expression for the chain term of the SAFT equation as a function of temperature
            in [K], `T`, and density in [mol/m3], `rho`. 
        alpha_disp: sp.Expr
            SymPy expression for the dispersion term of the SAFT equation as a function of temperature
            in [K], `T`, and density in [mol/m3], `rho`. 
        alpha_assoc: sp.Expr
            SymPy expression for the association term of the SAFT equation as a function of temperature
            in [K], `T`, and density in [mol/m3], `rho`. If the EOS does not have a multipolar term then
            alpha_assoc will be a sympyfied 0.
        alpha_multipolar: sp.Expr
            SymPy expression for the multipolar term of the SAFT equation as a function of temperature
            in [K], `T`, and density in [mol/m3], `rho`. If the EOS does not have a multipolar term then
            alpha_multipolar will be a sympyfied 0.
        """
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
    def for_purefluid(cls, fluid: PureFluid, 
                      epsilon:float = 0, sigma:float = 0, m:float = 0, 
                      epsilon_AB:float = 0, k_AB:float = 0, M:int = 0, 
                      x_p:float = 0, 
                      a:list=[], b:list=[], 
                      association_scheme:str = "", 
                      alpha_hs: sp.Expr = None, alpha_chain: sp.Expr = None, alpha_disp: sp.Expr = None, alpha_assoc: sp.Expr = None, alpha_multipolar: sp.Expr = None,
                      **kwargs):
        """Create the non-dimensional residual Helmholtz energy of a SAFT EOS for a pure fluid.

        Parameters
        ----------
        fluid: PureFluid
            PureFluid object corresponding to the relevant fluid.
        epsilon: float
            Segment dispersion energy divided by the Boltzmann constant in [K].
        sigma: float
            Segment diameter in [A] (angstroms).
        m: float
            Number of segments (dimensionless).
        epsilon_AB: float
            Energy of association divided by the Boltzmann constant in [K].
        k_AB: float
            Volume of association between sites A and B (dimensionless).
        M: int
            Number of associations per segment. This number will correspond to the first character
            of the association_scheme string.
        x_p: float
            Fraction of dipolar segments, dependent on chain length `m`. Some sources will define `x_p^{\mu}`
            as the fraction of dipolar segments while other will define `x_p^{\mu}` as the number of polar
            spherical segments in the molecule which is independent of chain length.
        a: list
            List of a coefficients used in the definition of the dispersion term, `alpha_disp`.
            The values for theses coefficients can be taken either from [0325] or [0329] depending
            on the EOS.
        b: list
            List of b coefficients used in the definition of the dispersion term, `alpha_disp`.
            The values for theses coefficients can be taken either from [0325] or [0329] depending
            on the EOS.
        association_scheme: str
            String representing the association scheme ofthe SAFT EOS, where the first character is an
            integer representing the number of associations per molecule and the second character is
            a letter characterising the association strengths of each association [0328].
        alpha_hs: sp.Expr
            If not given, then it will be calculated.
            SymPy expression for the hard-sphere term of the SAFT equation as a function of temperature
            in [K], `T`, and density in [mol/m3], `rho`. 
        alpha_chain: sp.Expr
            If not given, then it will be calculated.
            SymPy expression for the chain term of the SAFT equation as a function of temperature
            in [K], `T`, and density in [mol/m3], `rho`. 
        alpha_disp: sp.Expr
            If not given, then it will be calculated.
            SymPy expression for the dispersion term of the SAFT equation as a function of temperature
            in [K], `T`, and density in [mol/m3], `rho`. 
        alpha_assoc: sp.Expr
            If not given, then it will be calculated.
            SymPy expression for the association term of the SAFT equation as a function of temperature
            in [K], `T`, and density in [mol/m3], `rho`. If the EOS does not have a multipolar term then
            alpha_assoc will be a sympyfied 0.
        alpha_multipolar: sp.Expr
            If not given, then it will be calculated.
            SymPy expression for the multipolar term of the SAFT equation as a function of temperature
            in [K], `T`, and density in [mol/m3], `rho`. If the EOS does not have a multipolar term then
            alpha_multipolar will be a sympyfied 0.
        kwargs: dict, optional
            Extra keyword arguments.

        References
        ----------
        [0325] Gross, J., & Sadowski, G. (2000). Application of perturbation theory to a hard-chain reference fluid: an equation of state for square-well chains. Fluid Phase Equilibria, 168(2), 183-199.

        [0328] Huang, S. H., & Radosz, M. (1990). Equation of state for small, large, polydisperse, and associating molecules. Industrial & Engineering Chemistry Research, 29(11), 2284-2294.

        [0329] Gross, J., & Sadowski, G. (2001). Perturbed-chain SAFT: An equation of state based on a perturbation theory for chain molecules. Industrial & engineering chemistry research, 40(4), 1244-1260.
        """        
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
    def alpha_hs(m: float, zeta_n: list[sp.Expr, sp.Expr, sp.Expr, sp.Expr]) -> sp.Expr:
        r"""Hard-sphere term as a function of temperature in [K], `T`, and number density in 
        [molecules/A], `rho_N`.

        Defined in [0329] as:

        .. math:
        \alpha^{hs}= m \frac{1}{\zeta_0} \left(\frac{3 \zeta_1 \zeta_2}{1-\zeta_3}  + \frac{\zeta_2^3}{\zeta_3(1-\zeta_3)^2} + \left(\frac{\zeta_2^3}{\zeta_3^2} - \zeta_0 \right) \ln (1-\zeta_3)   \right)

        where:

        .. math:
        \zeta_n = \frac{\pi}{6} \rho_N m d^n
    
        where:

        .. math:
        d = \sigma (1 -0.12 \exp(-3 \frac{\epsilon}{k} \frac{1}{T}))


        Parameters
        ----------
        m: float
            Number of segments (dimensionless).
        zeta_n: list
            List of `zeta` terms between zeta_0 and zeta_3.

        References
        ----------
        [0329] Gross, J., & Sadowski, G. (2001). Perturbed-chain SAFT: An equation of state based on a perturbation theory for chain molecules. Industrial & engineering chemistry research, 40(4), 1244-1260.
        """
        # from [0329]
        alpha_hs = m / zeta_n[0] * ( 3*zeta_n[1]*zeta_n[2]/(1-zeta_n[3]) + zeta_n[2]**3/(zeta_n[3]*(1-zeta_n[3])**2) + (zeta_n[2]**3/zeta_n[3]**2 -zeta_n[0])*sp.log(1-zeta_n[3]))
        # print(f"alpha_hs = {alpha_hs}")
        return alpha_hs
    
    # chain
    @staticmethod
    def alpha_chain(m: float, d: sp.Expr, zeta_n: list[sp.Expr, sp.Expr, sp.Expr, sp.Expr]) -> sp.Expr:
        r"""Chain term as a function of temperature in [K], `T`, and number density in 
        [molecules/A], `rho_N`.

        Defined in [0329] as:

        .. math:
        \alpha^{chain}=-(m-1)\ln(g^{hs})

        where:

        .. math:
        g^{hs} = \frac{1}{1-\zeta_3}+\frac{d}{2}\frac{3 \zeta_2}{(1-\zeta_3)^2} + \frac{d}{2}\frac{2\zeta_2^2}{(1-\zeta_3)^3}

        where:

        .. math:
        \zeta_n = \frac{\pi}{6} \rho_N m d^n
    
        where:

        .. math:
        d = \sigma (1 -0.12 \exp(-3 \frac{\epsilon}{k} \frac{1}{T}))

        Parameters
        ----------
        m: float
            Number of segments (dimensionless).
        d: sp.Expr
            The temperature dependent segment diameter in [A] (angstroms).
        zeta_n: list[sp.Expr, sp.Expr, sp.Expr, sp.Expr]
            List of `zeta` terms between zeta_0 and zeta_3.

        References
        ----------
        [0329] Gross, J., & Sadowski, G. (2001). Perturbed-chain SAFT: An equation of state based on a perturbation theory for chain molecules. Industrial & engineering chemistry research, 40(4), 1244-1260.     
        """

        g_hs = 1/(1-zeta_n[3]) + d/2* 3*zeta_n[2]/(1-zeta_n[3])**2 + d**2/4 * 2*zeta_n[2]**2/(1-zeta_n[3])**3

        alpha_chain = -(m-1)*sp.log(g_hs)
        # print(f"alpha_chain = {alpha_chain}")
        return alpha_chain
    
    # dispersion
    @staticmethod
    def alpha_disp(m: float, eta: sp.Expr, a:list, b:list, m2e1sigma3: sp.Expr, m2e2sigma3:sp.Expr) -> sp.Expr:
        r"""Dispersion term as a function of temperature in [K], `T`, and number density in 
        [molecules/A], `rho_N`.

        Defined by [0329] as:

        .. math:
        \alpha^{disp}=-2 \pi \rho_N I_1 m^2 \frac{\epsilon}{kT}\sigma^3 - \pi \rho_N m C_1 I_2 m^2 \left(\frac{\epsilon}{kT}\right)^2 \sigma^3
        
        where

        .. math:
        I_1 = \sum_{i=0}^{6} a_i \eta^i

        .. math:
        I_2 = \sum_{i=0}^{6} b_i \eta^i

        .. math:
        C_1 = \left(1 + m \frac{8\eta-2\eta^2}{(1-\eta)^4} +(1-m)\frac{20\eta-27\eta^2+12\eta^3-2\eta^4}{((1-\eta)(2-\eta))^2} \right)^{-1}

        where

        .. math:
        \eta = \zeta_3 = \frac{\pi}{6}\rho_N d^3

        where

        .. math:
        d = \sigma (1 -0.12 \exp(-3 \frac{\epsilon}{k} \frac{1}{T}))


        .. math:
        a_i = a_{0i} + \frac{m-1}{m}a_{1i} + \frac{m-1}{m}\frac{m-2}{m}a_{2i}

        .. math:
        b_i = b_{0i} + \frac{m-1}{m}b_{1i} + \frac{m-1}{m}\frac{m-2}{m}b_{2i}

        Parameters
        ----------
        m: float
            Number of segments.
        eta: sp.Expr
            SymPy expression of the packing factor as a function of temperature in [K], `T`,
            and number density in [molecules/m3], `rho_N`. The packing fraction is defined as
            `pi/6*rho*m*d**n`, which is equivalent to `zeta_3`.
        a: list
            List of coefficients.
            The values for theses coefficients can be taken either from [0325] or [0329] depending
            on the EOS.
        b: list
            List of coefficients.
            The values for theses coefficients can be taken either from [0325] or [0329] depending
            on the EOS.        
        m2e1sigma3: sp.Expr
            SymPy expression defined as `m**2*epsilon/T*sigma**3`
        m2e2sigma3: sp.Expr
            SymPy expression defined as `m**2*(epsilon/T)**2*sigma**3`

        References
        ----------
        [0325] Gross, J., & Sadowski, G. (2000). Application of perturbation theory to a hard-chain reference fluid: an equation of state for square-well chains. Fluid Phase Equilibria, 168(2), 183-199.

        [0329] Gross, J., & Sadowski, G. (2001). Perturbed-chain SAFT: An equation of state based on a perturbation theory for chain molecules. Industrial & engineering chemistry research, 40(4), 1244-1260.
        """

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
    def alpha_assoc(association_scheme: str ="", Delta: sp.Expr = 0) -> sp.Expr:
        r"""Association term as a function of temperature in [K], `T`, and number density in 
        [molecules/A], `rho_N`.

        The association term is defined in [0328] as:

        .. math:
          \alpha^{assoc} = \sum_A \left(\ln X^A - \frac{X^A}{2}\right) + \frac{1}{2} M

        Where the definition of the fraction of segments not connected at site A, `X_A`,
        is dependent on both the assocition scheme chosen as well as the association strength,
        `Delta`. The association strength is defined as:

        .. math:
        \Delta=g^{hs}\left(\exp\left(\frac{\epsilon^{AB}}{kT}\right) -1\right)\sigma^3 \kappa^{AB}

        Parameters
        ----------
        association_scheme: str
            String defining the association scheme, where the first character is the number of
            association sites per molecule and the second character defines the associations
            themselves.
        Delta: sp.Expr
            Association strength as defined in [0328].

        References
        ----------
        [0328] Huang, S. H., & Radosz, M. (1990). Equation of state for small, large, polydisperse, and associating molecules. Industrial & Engineering Chemistry Research, 29(11), 2284-2294.
        """

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
    def alpha_multipolar(fluid: PureFluid , z: float , sigma: float, m: float, x_p: float) -> sp.Expr:
        r"""Multipolar term as a function of temperature in [K], `T`, and number density in 
        [molecules/A], `rho_N`.

        Currently can only be applied to dipolar and apolar fluid.
        This term is calculated using Gaussian units (also known as centimeter-gram-second).
        Defined using the equations from [0324], [0326] and [0333] as:

        .. math:
        \alpha^{multipolar}=\frac{\alpha_2}{1-\alpha_3/\alpha_2}

        where

        .. math:
        \alpha_2 = \alpha_2^{mult}(112)=-\frac{2}{3}\frac{\pi \rho_N}{(k_BT)^2} (x^{\mu}_{p}m)^2  m^2 \frac{\mu^4}{\sigma^3}J^{(6)}

        .. math:
        \alpha_3 = \alpha_{3B}^{mult}(112;112;112)=\frac{32 \pi^3}{135} \left(\frac{14 \pi}{5} \right)^{1/2} \frac{\rho_N^2}{(k_B T)^3} (x_p^{\mu}m)^3 m^3 \frac{\mu^6}{ \sigma^3} K(222;333) 

        where

        .. math:
        \ln|J^{(n)}| = A_n \rho^{*2}\ln T^* + B_n \rho^{*2} + C_n \rho^* \ln T^* + D_n \rho^* +E_n \ln T^* + F_n


        .. math:
        \ln|K^{(n)}| = A_n \rho^{*2}\ln T^* + B_n \rho^{*2} + C_n \rho^* \ln T^* + D_n \rho^* +E_n \ln T^* + F_n


        The non-dimentional density, `rho*`, and non-dimensional temperature, `T*`, are defined as:

        .. math:
        \rho^*=\rho_N\sigma^3

        .. math:
        T^*=\frac{kT\sigma^3}{\mu^2}

        
        Parameters
        ----------
        fluid: PureFluid
            PureFluid object corresponding to the relevant fluid. This is used to obtain the dipole
            moment of the fluid.
        z: float
            Mole fraction of relevant fluid. For a pure fluid EOS `z = 1`.
        sigma: float
            Segment length in [A] (angstroms).
        m: float
            Number of segments (dimensionless).
        x_p: float
            Fraction of dipolar segments.

        References
        ----------
        [0324] Grandjean, L., de Hemptinne, J. C., & Lugo, R. (2014). Application of GC-PPC-SAFT EoS to ammonia and its mixtures. Fluid Phase Equilibria, 367, 159-172.

        [0326] NguyenHuynh, D., Passarello, J. P., Tobaly, P., & de Hemptinne, J. C. (2008). Application of GC-SAFT EOS to polar systems using a segment approach. Fluid Phase Equilibria, 264(1-2), 62-75.

        [0333] Gubbins, K. E., & Twu, C. H. (1978). Thermodynamics of polyatomic fluid mixtures—I theory. Chemical Engineering Science, 33(7), 863-878.
        """
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
    
    
    
    