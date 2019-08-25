from fd import *
from ud import *

expression = "9*x+3"                  # Division by zero
# expression = "exp(x) - 3*x"
x0 = 0.5
tol = 0.000001

print("\n\nSteffensen's: ")
print(sne_fd_1(expression, x0, tol))
print("\n\nHalley's: ")
print(sne_ud_1(expression, x0, tol))
print("\n\nFrontini's y Sormani's: ")
print(sne_ud_2(expression, x0, tol))
print("\n\n")
