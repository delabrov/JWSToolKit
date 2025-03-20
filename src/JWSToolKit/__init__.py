try:
    from .Cube import Cube
except Exception as e:
    print("Error importing Cube module:", e)
    raise

try:
    from .Image import Image
except Exception as e:
    print("Error importing Image module:", e)
    raise

try:
    from .Spec import Spec
except Exception as e:
    print("Error importing Spec module:", e)
    raise

import os
import sys

sys.path.insert(0, os.path.dirname(__file__))

__version__ = "1.1.2"
__author__ = "Valentin Delabrosse"
__email__ = "valentin.delabrosse@univ-grenoble-alpes.fr"

__all__ = ["Cube", "Image", "Spec"]

