import os
print(f"current running path is {os.getcwd()}")

import sys
sys.path.append(os.getcwd())

from build import pyxtract