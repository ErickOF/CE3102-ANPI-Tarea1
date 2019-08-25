from fd import *
from ud import *

# expression = "x^3-x-2"                  # Division by zero
expression = "exp(x) - 3*x"
x0 = 1
tol = 0.000001

print("\n\nSteffensen's: ")
print(sne_fd_1(expression, x0, tol))
print("\n\nHalley's: ")
print(sne_ud_1(expression, x0, tol))
print("\n\nFrontini's y Sormani's: ")
print(sne_ud_2(expression, x0, tol))
print("\n\n")
