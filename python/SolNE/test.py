from fd import *
from ud import *

print("\n\nSteffensen's: ")
print(sne_fd_1("x^2 - 3", 2, 0.000001))
print("\n\nHalley's: ")
print(sne_ud_1("x^2 - 3", 2, 0.000001))
print("\n\n")
