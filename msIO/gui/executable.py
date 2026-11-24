import os
import sys

# Windowed PyInstaller builds have no console: sys.stdout / sys.stderr are None, so print(), tqdm and logging
# would crash with "'NoneType' object has no attribute 'write'".
if sys.stdout is None:
    sys.stdout = open(os.devnull, "w")
if sys.stderr is None:
    sys.stderr = open(os.devnull, "w")