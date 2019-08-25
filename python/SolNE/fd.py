# =========== Important: ¡must install Sympy! (pip install sympy) ============
from sympy import *


# ========================== Global variables ================================
global x                                    # x (is a global symbol)
global ITER_LIMIT                           # Limit of iterations
global DECIMAL_PRECISION                    # Amount of decimal values

x = symbols('x')
ITER_LIMIT = 10000
DECIMAL_PRECISION = 100


# ============================== Method 1 ====================================
def sne_fd_1(expr, x0, tol):
    """Steffensen's Method

    Arguments:
        expr {string} -- polynomial whose solution must be found
        x0 {float, int} -- initial value to start iterations
        tol {float, int} -- tolerance that indicates the stop condition

    Returns:
        xn {float} -- root approximation
        itera {int} -- amount of iterations required
        graph {int} -- flag that indicates if a graph must be done

    """

    # -------------------------- Validations ---------------------------------
    if (validator(expr, x0, tol) != True):
        return x0, 0, 0
    # ------------------------ Local variables -------------------------------
    f = sympify(expr)                       # Transforms string to function
    graph = 0                               # Wether the graph will be shown
    itera = 0                               # Amount of iterations
    xn = sympify(x0)                        # xn is a Sympy variable
    xNext = sympify(0)                      # xNext (x_(n+1))
    error = abs(f.subs(x, xn))              # calculates the error of x0

    try:
        # -------------------- Steffensen's Method ---------------------------
        while (error > tol):
            if(itera >= ITER_LIMIT):
                print("WARNING: Iteration limit reached")
                return Float(xn, DECIMAL_PRECISION), itera, graph
            xNext = \
                xn - (f.subs(x, xn) * f.subs(x, xn)) / \
                (f.subs(x, xn + f.subs(x, xn)) - f.subs(x, xn))
            xn = Float(xNext, DECIMAL_PRECISION)
            error = abs(f.subs(x, xn))
            itera += 1                      # New iteration
        graph = 1                           # Graph can be displayed

    except Exception as exception:
        print("WARNING: [Math error]", type(exception).__name__)

    return Float(xn, DECIMAL_PRECISION), itera, graph


# ============================== Method 2 ====================================
def sne_fd_2(*args, **kwargs):
    pass


# ============================== Method 3 ====================================
def sne_fd_3(*args, **kwargs):
    pass


# ============================== Method 4 ====================================
def sne_fd_4(*args, **kwargs):
    pass


# ============================== Method 5 ====================================
def sne_fd_5(*args, **kwargs):
    pass


# ============================== Method 6 ====================================
def sne_fd_6(*args, **kwargs):
    pass


# =========================== Auxiliary function =============================
def validator(expr, x0, tol):
    """Auxiliary function that validates inputs

    Arguments:
        expr {string} -- polynomial whose solution must be found
        x0 {float, int} -- initial value to start iterations
        tol {float, int} -- tolerance that indicates the stop condition

    Returns:
        flag {boolean} -- true if all inputs are valid

    """

    # -------------------------- Validations ---------------------------------
    if (type(expr) != str):
        print("WARNING: invalid expression")
        return False
    if (type(x0) != float and type(x0) != int):
        print("WARNING: x_0 must be a number")
        return False
    if (type(tol) != float and type(tol) != int):
        print("WARNING: tolerance must be a number")
        return False
    try:
        f = sympify(expr)
    except Exception as exception:
        print("WARNING: [invalid expression]: ", type(exception).__name__)
        return False
    return True
