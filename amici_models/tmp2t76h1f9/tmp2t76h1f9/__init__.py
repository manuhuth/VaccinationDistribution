"""AMICI-generated module for model tmp2t76h1f9"""

import amici

# Ensure we are binary-compatible, see #556
if "0.11.16" != amici.__version__:
    raise RuntimeError(
        "Cannot use model tmp2t76h1f9, generated with AMICI "
        "version 0.11.16, together with AMICI version"
        f" {amici.__version__} which is present in your "
        "PYTHONPATH. Install the AMICI package matching the "
        "model version or regenerate the model with the AMICI "
        "currently in your path."
    )

from tmp2t76h1f9._tmp2t76h1f9 import *

__version__ = "0.1.0"
