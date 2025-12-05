__all__ = ["PureFluid",
           "symbols",
           "storeroom"]

import importlib

from .purefluid import PureFluid
from .symbols import *

# Lazy imports
def __getattr__(name):
    # 2. Check if the user is asking for 'storeroom'
    if name == "storeroom":
        # 3. Only NOW do we import the subpackage
        # print("Lazy loading 'storeroom' subpackage...") # Debug print
        return importlib.import_module(f".{name}", __name__)
    
    # 4. Standard error if they ask for something that doesn't exist
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


