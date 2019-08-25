from fd import *
from ud import *

expression = "x^2 - 3"
x0 = 2
tol = 0.000001

print("\n\nSteffensen's: ")
print(sne_fd_1(expression, x0, tol))
print("\n\nHalley's: ")
print(sne_ud_1(expression, x0, tol))
print("\n\nFrontini's y Sormani's: ")
print(sne_ud_2(expression, x0, tol))
print("\n\n")
