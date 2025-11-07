# 07/11/2025

# create presets for all the pure fluid EOS relevant to the paper

from purefluid import NH3, H2

# classes
from eos import EOS
from alpha_r_cubic import AlphaRCubic
from alpha_r_helmholtz import AlphaRHelmholtz
from alpha_r_saft import AlphaRSaft

# objects
from purefluid import NH3, H2
from cp0 import NH3_cp0_NIST, H2_cp0_NIST 
from presets_alpha_r import alpha_r_0290, alpha_r_0300, alpha_r_0313


# algorithms and other functions and variables
from symbols import *
from mixture_rules import one_fluid_theory
from other_functions import multi_root
