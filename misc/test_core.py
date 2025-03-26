import sys
import time
import numpy as np
from ardal import Ardal
import _ardal

data = ["./data/usher_matrix.npy", "./data/usher_headers.json"]
ard = Ardal(data)

ard.stats()

core = ard.core(["B.1.1.406", "B.1.1.406"])

print(core)