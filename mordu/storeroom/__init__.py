"""Lazy-load storeroom submodules/objects on attribute access (PEP 562)."""

import importlib
from typing import Dict

# map attribute names to module paths under mordu.storeroom
_submodule_names: Dict[str, str] = {
    "fluids": "mordu.storeroom.fluids",
    "eos_purefluid": "mordu.storeroom.eos_purefluid",
    "cp0s": "mordu.storeroom.cp0s",
    "alpha_r_func": "mordu.storeroom.alpha_r_func",
}

# exported names visible from `from mordu.storeroom import ...`
__all__ = list(_submodule_names.keys())


def __getattr__(name: str):
    """Lazy import a submodule when accessed (e.g. storeroom.fluids)."""
    if name in _submodule_names:
        module = importlib.import_module(_submodule_names[name])
        # cache on the package module to avoid repeated imports
        globals()[name] = module
        return module
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def __dir__():
    return sorted(list(globals().keys()) + __all__)


