# =========== Important: ¡must install Sympy! (pip install sympy) ============
import matplotlib.pyplot as plt
import numpy as np
from scipy import misc
from sympy import *


# ========================== Global variables ================================
global x                                    # x (is a global symbol)
global ITER_LIMIT                           # Limit of iterations
global DECIMAL_PRECISION                    # Amount of decimal values

x = symbols('x')
ITER_LIMIT = 10000
DECIMAL_PRECISION = 100


# ============================== Method 1 ====================================
def sne_ud_1(expr, x0, tol):
    """Halley's Method

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
                return N(xn, DECIMAL_PRECISION), itera, graph
            fDiff = diff(f, x)              # First derivative
            fDiff2 = diff(fDiff, x)         # Second derivative
            div = (2*(fDiff.subs(x, xn))*(fDiff.subs(x, xn)) -
                   f.subs(x, xn)*fDiff2.subs(x, xn))
            if (div != 0):
                xNext = xn - (2*(f.subs(x, xn))*(fDiff.subs(x, xn))) / div
            else:
                print("WARNING: [Math error] Division by zero")
                return N(xn, DECIMAL_PRECISION), itera, graph
            xn = N(xNext, DECIMAL_PRECISION)
            error = abs(f.subs(x, xn))
            itera += 1                      # New iteration
        graph = 1                           # Graph can be displayed

    except Exception as exception:
        print("WARNING: [Math error]", type(
            exception).__name__, "-", str(exception))

    return N(xn, DECIMAL_PRECISION), itera, graph

# ============================== Method 2 ====================================
def sne_ud_2(expr, x0, tol):
    """Frontini's y Sormani's Method

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
                return N(xn, DECIMAL_PRECISION), itera, graph
            fDiff = diff(f, x)
            div = fDiff.subs(x, xn)
            if (div != 0):
                div2 = fDiff.subs(x, (xn - 0.5 * f.subs(x, xn)/div))
                if (div2 != 0):
                    xNext = xn - (f.subs(x, xn) / div2)
                else:
                    print("WARNING: [Math error] Division by zero")
                    return N(xn, DECIMAL_PRECISION), itera, graph
            else:
                print("WARNING: [Math error] Division by zero")
                return N(xn, DECIMAL_PRECISION), itera, graph

            xn = N(xNext, DECIMAL_PRECISION)
            error = abs(f.subs(x, xn))
            itera += 1                      # New iteration
        graph = 1                           # Graph can be displayed

    except Exception as exception:
        print("WARNING: [Math error]", type(
            exception).__name__, "-", str(exception))

    return N(xn, DECIMAL_PRECISION), itera, graph

# ============================== Method 3 ====================================
def sne_ud_3(f, x0, tol, graf=1):
    if (not isinstance(f, str)):
        raise ValueError('f must be a string')
    
    if (not isinstance(x0, (int, float))):
        raise ValueError('x0 must be a int or float')
    
    if (not isinstance(tol, (int, float))):
        raise ValueError('tol must be a int or float')
    
    if (graf != 0 and graf != 1):
        raise ValueError('graf must be 0 or 1')

    xAprox = np.array([x0])
    _iter = 0
    
    try:
        fx = lambda x: eval(f)
        error = np.array([abs(fx(xAprox[-1]))])

        while (abs(fx(xAprox[-1])) > tol):
            xk = xAprox[-1]

            y = fx(xk)
            df = misc.derivative(fx, xk, dx=1e-6)
            df2 = misc.derivative(fx, xk, n=2, dx=1e-6)

            Lf = y * df2 / df**2

            xk_next = xk - (1 + 0.5 * Lf) * y / df

            xAprox = np.append(xAprox, xk_next)
            error = np.append(error, abs(fx(xk_next)))

            _iter += 1
        
        if graf == 1:
            k = np.linspace(0, _iter, _iter + 1)
            plotFunction(k, error, 'Chebyshev Method')

        return xAprox[-1], _iter
    except AttributeError as e:
        raise ValueError('f has an unknown function. ' + str(e).capitalize())
    except TypeError as e:
        raise ValueError('f has an unknown symbol. ' + str(e).capitalize())

# ============================== Method 4 ====================================
def sne_ud_4(f, x0, tol, graf=1):
    if (not isinstance(f, str)):
        raise ValueError('f must be a string')
    
    if (not isinstance(x0, (int, float))):
        raise ValueError('x0 must be a int or float')
    
    if (not isinstance(tol, (int, float))):
        raise ValueError('tol must be a int or float')
    
    if (graf != 0 and graf != 1):
        raise ValueError('graf must be 0 or 1')

    xAprox = np.array([x0])
    _iter = 0
    
    try:
        fx = lambda x: eval(f)
        error = np.array([abs(fx(xAprox[-1]))])

        while (abs(fx(xAprox[-1])) > tol):
            xk = xAprox[-1]

            y = fx(xk)
            df = misc.derivative(fx, xk, dx=1e-6)
            yk = xk - y / df

            xk_next = xk - (y / (y - fx(yk))) * y / df

            xAprox = np.append(xAprox, xk_next)
            error = np.append(error, abs(fx(xk_next)))

            _iter += 1
        
        if graf == 1:
            k = np.linspace(0, _iter, _iter + 1)
            plotFunction(k, error, 'Newton-Secant Method')

        return xAprox[-1], _iter
    except AttributeError as e:
        raise ValueError('f has an unknown function. ' + str(e).capitalize())
    except TypeError as e:
        raise ValueError('f has an unknown symbol. ' + str(e).capitalize())

# ============================== Method 5 ====================================
def sne_ud_5(f, x0, tol, graf=1):
    if (not isinstance(f, str)):
        raise ValueError('f must be a string')
    
    if (not isinstance(x0, (int, float))):
        raise ValueError('x0 must be a int or float')
    
    if (not isinstance(tol, (int, float))):
        raise ValueError('tol must be a int or float')
    
    if (graf != 0 and graf != 1):
        raise ValueError('graf must be 0 or 1')

    xAprox = np.array([x0])
    _iter = 0
    
    try:
        fx = lambda x: eval(f)
        error = np.array([abs(fx(xAprox[-1]))])

        while (abs(fx(xAprox[-1])) > tol):
            xk = xAprox[-1]

            yk = fx(xk)
            df = misc.derivative(fx, xk, dx=1e-6)
            df2 = misc.derivative(fx, xk, n=2, dx=1e-6)
            df3 = misc.derivative(fx, xk, n=3, order=5, dx=1e-6)
            
            Lf = yk * df2 / df**2
            Ldf = df * df3 / df2**2

            xk_next = xk - 1.5 * ((2 - Lf)**2 / \
                      (6 - 9 * Lf + 3 * Lf**2 + Lf**2 * Ldf)) * (yk / df)

            xAprox = np.append(xAprox, xk_next)
            error = np.append(error, abs(fx(xk_next)))

            _iter += 1
        
        if graf == 1:
            k = np.linspace(0, _iter, _iter + 1)
            plotFunction(k, error, 'Danby Burkardt Method')

        return xAprox[-1], _iter
    except AttributeError as e:
        raise ValueError('f has an unknown function. ' + str(e).capitalize())
    except TypeError as e:
        raise ValueError('f has an unknown symbol. ' + str(e).capitalize())

# ============================== Method 6 ====================================
def sne_ud_6(f, x0, tol, graf=1):
    if (not isinstance(f, str)):
        raise ValueError('f must be a string')
    
    if (not isinstance(x0, (int, float))):
        raise ValueError('x0 must be a int or float')
    
    if (not isinstance(tol, (int, float))):
        raise ValueError('tol must be a int or float')
    
    if (graf != 0 and graf != 1):
        raise ValueError('graf must be 0 or 1')

    xAprox = np.array([x0])
    _iter = 0
    
    try:
        fx = lambda x: eval(f)
        error = np.array([abs(fx(xAprox[-1]))])

        while (abs(fx(xAprox[-1])) > tol):
            xk = xAprox[-1]

            yk = fx(xk)
            df = misc.derivative(fx, xk, dx=1e-6)
            df2 = misc.derivative(fx, xk, n=2, dx=1e-6)
            zk = -yk / df

            xk_next = xk - yk * (df + 0.5 * zk * df2)**-1

            xAprox = np.append(xAprox, xk_next)
            error = np.append(error, abs(fx(xk_next)))

            _iter += 1
        
        if graf == 1:
            k = np.linspace(0, _iter, _iter + 1)
            plotFunction(k, error, 'Newton-Secant Method')

        return xAprox[-1], _iter
    except AttributeError as e:
        raise ValueError('f has an unknown function. ' + str(e).capitalize())
    except TypeError as e:
        raise ValueError('f has an unknown symbol. ' + str(e).capitalize())

# =========================== Auxiliary functions =============================
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
        print("WARNING: [invalid expression]: ", type(exception).__name__,
              "-", str(exception))
        return False
    return True

def plotFunction(k, error, title):
    plt.title(title)
    plt.xlabel('Iterations k')
    plt.ylabel('Error |f(xk)|')
    plt.plot(k, error)
    plt.show()
