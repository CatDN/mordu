#class for mixtures of fluids
from dataclasses import dataclass
import sympy as sp
from symbols import *

@dataclass
class MixtureFluid():
    fluid_1: object
    fluid_2: object
    z1: float = z1
    z2: float = z2

    def __post_init__(self):
        if self.z1+self.z2 !=1 and (self.z1!=z1 and self.z2!=z2):
            raise ValueError("The sum of the molar fractions must be 1, or the mole fractions must be symbolic")
        
        self.M = self.z1*self.fluid_1.M + self.z2*self.fluid_2.M

    def __str__(self):
        return f"MixtureFluid({self.fluid_1.formula}, {self.z1}, {self.fluid_2.formula}, {self.z2})"


# #example
# from PureFluid import H2, NH3
# eta = 0.6525
# mix = MixtureFluid(H2, NH3)
# print(mix.M)