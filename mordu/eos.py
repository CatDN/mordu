#27/02/2025

from .alpha_0 import Alpha0
from .symbols import *
from .purefluid import PureFluid
from .cp0 import Cp0

#EOS parent class
class EOS:
    """Class for the creation of all pure fluid equations of state (EOS).

    Attributes
    ----------
    name: str
        Name of the EOS.
    fluid: PureFluid
        PureFluid object corresponding to the fluid for which the EOS is beeing created.
    cp0_object: Cp0
        Cp0 object corresponding to the ideal gas isobaric heat capacity for the specific EOS.
    AlphaRClass: AlphaRCubic or AlphaRHelmholtz or AlphaRSAFT
        Non-dimensional residual Helmholtz energy class which depends on the type of EOS
        beeing created. Currently the only options available are: AlphaRCubic or 
        AlphaRHelmholtz or AlphaRSAFT.
    ref: str, optional
        Reference of where the equation of state was first defined (such as an academic 
        paper or textbook).
    kwargs: 
        Optional keyword arguments for the AlphaRClass chosen and then applied using the following
        class method:

        * mordu.alpha_r_cubi.AlphaRCubic.for_purefluid()

        * mordu.alpha_r_helmholtz.AlphaRHelmholtz.for_purefluid()

        * mordu.alpha_r_saft.AlphaRSAFT.for_purefluid()
    

    """

    def __init__(self, name:str, fluid: PureFluid, cp0_object: Cp0, AlphaRClass, ref:str = "", **kwargs):
        """Creating an equation of state for a single pure fluid.

        Parameters
        ----------
        name: str
            Name of the equation of state
        fluid: PureFluid
            PureFluid object corresponding to the fluid for which the EOS is beeing created.
        cp0_object: Cp0
            Cp0 object corresponding to the ideal gas isobaric heat capacity for the specific EOS.
        AlphaRClass: AlphaRCubic or AlphaRHelmholtz or AlphaRSAFT
            Non-dimensional residual Helmholtz energy class which depends on the type of EOS
            beeing created. Currently the only options available are: AlphaRCubic or 
            AlphaRHelmholtz or AlphaRSAFT.
        ref: str, optional
            Reference of where the equation of state was first defined (such as an academic 
            paper or textbook).
        kwargs: 
            Optional keyword arguments for the AlphaRClass chosen and then applied using the following
            class method:

            * mordu.alpha_r_cubi.AlphaRCubic.for_purefluid()

            * mordu.alpha_r_helmholtz.AlphaRHelmholtz.for_purefluid()

            * mordu.alpha_r_saft.AlphaRSAFT.for_purefluid()            

        """
        self.fluid = fluid
        
        self.alpha_0 = Alpha0(self.fluid, cp0_object)
        self.alpha_r = AlphaRClass.for_purefluid(self.fluid, **kwargs)

        self.alpha = (self.alpha_0.alpha_0 + self.alpha_r.alpha_r)

        self.name = name

        self.pressure_equation = sp.lambdify((rho, P, T), P- self.pressure)

        self.ref = ref

    @property
    def pressure(self) -> sp.Expr:
        r"""SymPy expression for pressure in [Pa], `P`, as a function of temperature in [K], `T`, and
        density in [mol/m3], `rho`.

        .. math:

        P=\rho R T (1+\rho \frac{\partial \alpha^r}{\partial \rho})

        Notes
        -----
        See [0315], page 107, for this formula and others.

        References
        ----------
        [0315] Kunz, Oliver, Reinhard Klimeck, Wolfgang Wagner, and Manfred Jaeschke. "The GERG-2004 wide-range equation of state for natural gases and other mixtures." (2007).

        """
        pressure_expression = rho*self.fluid.R*T*(1+ rho*sp.diff(self.alpha_r.alpha_r, rho))
        return pressure_expression

    @property
    def compressibility(self):
        r"""SymPy expression for compressibility, `Z`, as a function of temperature in [K], `T`, and
        density in [mol/m3], `rho`.

        .. math:

        Z= 1+ \rho \frac{\partial \alpha^r}{\partial \rho})

        Notes
        -----
        See [0315], page 107, for this formula and others.

        References
        ----------
        [0315] Kunz, Oliver, Reinhard Klimeck, Wolfgang Wagner, and Manfred Jaeschke. "The GERG-2004 wide-range equation of state for natural gases and other mixtures." (2007).

        """
        compressibility_expression = (1+ rho*sp.diff(self.alpha_r.alpha_r, rho))
        return compressibility_expression

    @property
    def specific_heat(self):
        r"""SymPy expression for isobaric heat capacity (also known as specific heat), `c_p`,
        as a function of temperature in [K], `T`, and density in [mol/m3], `rho`.

        .. math:
        \tau = \frac{1}{T}

        c_p = R(-tau^2 (\alpha^0_{\tau \tau} + \alpha^r_{\tau \tau}) + \frac{(1+ \rho \alpha^r_{\rho} - \rho \tau \alpha^r_{\rho \tau})^2}{1 + 2 \rho \alpha^r_{\rho} + \rho^2 \alpha^r_{\rho \rho}})
        
        Notes
        -----
        See [0315], page 107, for this formula and others.

        References
        ----------
        [0315] Kunz, Oliver, Reinhard Klimeck, Wolfgang Wagner, and Manfred Jaeschke. "The GERG-2004 wide-range equation of state for natural gases and other mixtures." (2007).

        """
        alpha_0 = self.alpha_0.alpha_0.subs([(T, 1/T_inv)])
        alpha_r = self.alpha_r.alpha_r.subs([(T, 1/T_inv)])

        cp_expression = self.fluid.R*(-T_inv**2*(sp.diff(alpha_0, T_inv, T_inv) + sp.diff(alpha_r, T_inv, T_inv)) + \
                        (1+rho*sp.diff(alpha_r, rho) - rho*T_inv*sp.diff(alpha_r, rho, T_inv))**2/ \
                        (1+2*rho*sp.diff(alpha_r, rho) +rho**2*sp.diff(alpha_r,rho, rho)))
        cp_expression = cp_expression.subs([(T_inv, 1/T)])

        return cp_expression

    @property
    def joule_thomson_coefficient(self):
        r"""SymPy expression for Joule-Thomson coefficient, `mu_JT`,
        as a function of temperature in [K], `T`, and density in [mol/m3], `rho`.

        .. math:
        \tau = \frac{1}{T}

        \mu_{JT} = \frac{1}{R \rho} \frac{- (\rho \alpha^r_{\rho} + \rho^2 \alpha^r_{\rho \rho} + \rho \tau \alpha^r_{\rho \tau})}{(1 + \rho \alpha^r_{\rho} - \rho \tau \alpha^r_{\rho \tau})^2 - \tau^2(\alpha^0_{\tau \tau} + \alpha^r_{\tau \tau}) (1 + 2\rho \alpha^r_{\rho} + \rho^2 \alpha^r_{\rho \rho})}
        
        Notes
        -----
        See [0315], page 107, for this formula and others.

        References
        ----------
        [0315] Kunz, Oliver, Reinhard Klimeck, Wolfgang Wagner, and Manfred Jaeschke. "The GERG-2004 wide-range equation of state for natural gases and other mixtures." (2007).

        """
        alpha_0 = self.alpha_0.alpha_0.subs([(T, 1/T_inv)])
        alpha_r = self.alpha_r.alpha_r.subs([(T, 1/T_inv)])

        muJT_expression = 1/(self.fluid.R*rho)*-(rho*sp.diff(alpha_r, rho) +rho**2*sp.diff(alpha_r, rho, rho)+rho*T_inv*sp.diff(alpha_r, rho, T_inv))/ \
                            ((1+rho*sp.diff(alpha_r,rho) -rho*T_inv*sp.diff(alpha_r, rho, T_inv))-T_inv**2*(sp.diff(alpha_0, T_inv, T_inv)+sp.diff(alpha_r,T_inv,T_inv))*\
                                (1+2*rho*sp.diff(alpha_r, rho)+rho**2*sp.diff(alpha_r, rho, rho)))

        muJT_expression = muJT_expression.subs([(T_inv, 1/T)])

        return muJT_expression

    @property
    def gibbs_free_energy(self):
        r"""SymPy expression for Gibbs free energy, `g`,
        as a function of temperature in [K], `T`, and density in [mol/m3], `rho`.

        .. math:
        \tau = \frac{1}{T}

        g = 1 + \alpha^0 + \alpha^r + \rho \alpha^r_{\rho}
                        
        Notes
        -----
        See [0315], page 107, for this formula and others.

        References
        ----------
        [0315] Kunz, Oliver, Reinhard Klimeck, Wolfgang Wagner, and Manfred Jaeschke. "The GERG-2004 wide-range equation of state for natural gases and other mixtures." (2007).

        """
        alpha_0 = self.alpha_0.alpha_0.subs([(T, 1/T_inv)])
        alpha_r = self.alpha_r.alpha_r.subs([(T, 1/T_inv)])

        g_expression = self.fluid.R*T*(1+ alpha_0 + alpha_r + rho*sp.diff(alpha_r, rho))

        g_expression = g_expression.subs([(T_inv, 1/T)])

        return g_expression

    #see log page 163
    @property
    def fugacity_coefficient(self):
        r"""SymPy expression for the fugacity coefficient, `phi`,
        as a function of temperature in [K], `T`, and density in [mol/m3], `rho`.

        .. math:

        \phi = \exp(\alpha^r + \rho \alpha^r_{\rho} - \log( 1 + \rho \alpha^r_{\rho}))
                        
        Notes
        -----
        See [0315], page 107, for this formula and others.

        References
        ----------
        [0315] Kunz, Oliver, Reinhard Klimeck, Wolfgang Wagner, and Manfred Jaeschke. "The GERG-2004 wide-range equation of state for natural gases and other mixtures." (2007).

        """

        phi_expression = sp.exp(self.alpha_r.alpha_r + rho*sp.diff(self.alpha_r.alpha_r, rho) - sp.log(1+rho*sp.diff(self.alpha_r.alpha_r, rho)))

        return phi_expression  #fugacity coefficient phi as a function of density and temperature
    
