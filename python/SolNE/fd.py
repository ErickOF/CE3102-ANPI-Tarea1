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
DECIMAL_PRECISION = 50


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
                return N(xn, DECIMAL_PRECISION), itera, graph
            div = f.subs(x, xn + f.subs(x, xn)) - f.subs(x, xn)
            if (div != 0):
                xNext = xn - (f.subs(x, xn) * f.subs(x, xn)) / div
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
def sne_fd_2(f, x0, a0, b0, tol, graf=1):
    """
    Yun-Petkovic Method

    Métodos iterativos aplicados a la ecuación de Kepler. Page 112.

    Arguments:

        f  {string} - polynomial whose solution must be found
        x0 {float, int} - initial value to start iterations
        a0 {float, int} - interval low value
        b0 {float, int} - interval high value
        tol {float, int} - tolerance that indicates the stop condition
        graph {int} - flag that indicates if a graph must be done

    Returns:

        xn {float} - root approximation
        _iter {int} - amount of iterations required
    """
    if (not isinstance(f, str)):
        raise ValueError('f must be a string')
    
    if (not isinstance(x0, (int, float))):
        raise ValueError('x0 must be a int or float')
    
    if (not isinstance(tol, (int, float))):
        raise ValueError('tol must be a int or float')
    
    if (graf != 0 and graf != 1):
        raise ValueError('graf must be 0 or 1')

    xAprox = np.array([0, x0])
    _iter = 0
    
    try:
        fx = lambda x: eval(f)
        error = np.array([abs(fx(xAprox[-1]))])

        hk = 1
        ak = a0
        bk = b0

        while (abs(fx(xAprox[-1])) > tol):
            xk = xAprox[-1]

            xk_next = xk - fx(xk)*(2 * hk / (fx(bk) - fx(ak)))

            xAprox = np.append(xAprox, xk_next)
            error = np.append(error, abs(fx(xk_next)))

            hk = xk_next - xk
            ak = xk_next - hk
            bk = xk_next + hk

            _iter += 1
        
        if graf == 1:
            k = np.linspace(0, _iter, _iter + 1)
            plotFunction(k, error, 'Yun-Petkovic Method')

        return xAprox[-1], _iter
    except AttributeError as e:
        raise ValueError('f has an unknown function. ' + str(e).capitalize())
    except TypeError as e:
        raise ValueError('f has an unknown symbol. ' + str(e).capitalize())

# ============================== Method 3 ====================================
def sne_fd_3(f, x0, tol, graf=1):
    """
    Jain Method
    
    Metodos iterativos optimos para la resolucion de ecuaciones no lineales. Page 1. Equation 1.

    Arguments:

        f  {string} - polynomial whose solution must be found
        x0 {float, int} - initial value to start iterations
        tol {float, int} - tolerance that indicates the stop condition
        graph {int} - flag that indicates if a graph must be done

    Returns:

        xn {float} - root approximation
        _iter {int} - amount of iterations required
    """
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
            yk = steffensen_method(f, x0, _iter)

            xk_next = xk - y**3 / ((fx(xk + y) - y) * (y - fx(yk)))

            xAprox = np.append(xAprox, xk_next)
            error = np.append(error, abs(fx(xk_next)))

            _iter += 1
        
        if graf == 1:
            k = np.linspace(0, _iter, _iter + 1)
            plotFunction(k, error, 'Jain Method')

        return xAprox[-1], _iter
    except AttributeError as e:
        raise ValueError('f has an unknown function. ' + str(e).capitalize())
    except TypeError as e:
        raise ValueError('f has an unknown symbol. ' + str(e).capitalize())

# ============================== Method 4 ====================================
def sne_fd_4(f, x0, tol, graf=1):
    """
    Liu Method
    
    Metodos iterativos optimos para la resolucion de ecuaciones no lineales. Page 1. Equation 2.

    Arguments:

        f  {string} - polynomial whose solution must be found
        x0 {float, int} - initial value to start iterations
        tol {float, int} - tolerance that indicates the stop condition
        graph {int} - flag that indicates if a graph must be done

    Returns:

        xn {float} - root approximation
        _iter {int} - amount of iterations required
    """
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

            yk = steffensen_method(f, x0, _iter)
            zk = xk + fx(xk)

            f_xk_yk = (fx(yk) - fx(xk)) / (yk - xk)
            f_yk_zk = (fx(zk) - fx(yk)) / (zk - yk)
            f_xk_zk = (fx(zk) - fx(xk)) / (zk - xk)

            xk_next = yk - fx(yk)*((f_xk_yk - f_yk_zk + f_xk_zk) / f_xk_yk**2)

            xAprox = np.append(xAprox, xk_next)
            error = np.append(error, abs(fx(xk_next)))

            _iter += 1
        
        if graf == 1:
            k = np.linspace(0, _iter, _iter + 1)
            plotFunction(k, error, 'Liu Method')

        return xAprox[-1], _iter
    except AttributeError as e:
        raise ValueError('f has an unknown function. ' + str(e).capitalize())
    except TypeError as e:
        raise ValueError('f has an unknown symbol. ' + str(e).capitalize())

# ============================== Method 5 ====================================
def sne_fd_5(f, x0, tol, graf=1):
    """
    Ren Method
    
    Metodos iterativos optimos para la resolucion de ecuaciones no lineales. Page 2. Equation 2.

    Arguments:

        f  {string} - polynomial whose solution must be found
        x0 {float, int} - initial value to start iterations
        tol {float, int} - tolerance that indicates the stop condition
        graph {int} - flag that indicates if a graph must be done

    Returns:

        xn {float} - root approximation
        _iter {int} - amount of iterations required
    """
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

            yk = steffensen_method(f, x0, _iter)
            zk = xk + fx(xk)

            f_xk_yk = (fx(yk) - fx(xk)) / (yk - xk)
            f_yk_zk = (fx(zk) - fx(yk)) / (zk - yk)
            f_xk_zk = (fx(zk) - fx(xk)) / (zk - xk)

            div = f_xk_yk + f_yk_zk - f_xk_zk + (yk - xk) * (yk - zk)
            xk_next = yk - fx(yk) / div

            xAprox = np.append(xAprox, xk_next)
            error = np.append(error, abs(fx(xk_next)))

            _iter += 1
        
        if graf == 1:
            k = np.linspace(0, _iter, _iter + 1)
            plotFunction(k, error, 'Ren Method')

        return xAprox[-1], _iter
    except AttributeError as e:
        raise ValueError('f has an unknown function. ' + str(e).capitalize())
    except TypeError as e:
        raise ValueError('f has an unknown symbol. ' + str(e).capitalize())

# ============================== Method 6 ====================================
def sne_fd_6(f, x0, tol, graf=1):
    """
    Free Derivative Ostrowski Method
    
    Journal of Computational and Applied mathematics. Equation 4

    Arguments:

        f  {string} - polynomial whose solution must be found
        x0 {float, int} - initial value to start iterations
        tol {float, int} - tolerance that indicates the stop condition
        graph {int} - flag that indicates if a graph must be done

    Returns:

        xn {float} - root approximation
        _iter {int} - amount of iterations required
    """
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
            yk = xk - (2 * y**2) / (fx(xk + y) - fx(xk - y))

            xk_next = yk * (fx(yk) - y) / (2 * fx(yk) - y)

            xAprox = np.append(xAprox, xk_next)
            error = np.append(error, abs(fx(xk_next)))

            _iter += 1
        
        if graf == 1:
            k = np.linspace(0, _iter, _iter + 1)
            plotFunction(k, error, 'Free Derivative Ostrowski Method')

        return xAprox[-1], _iter
    except AttributeError as e:
        raise ValueError('f has an unknown function. ' + str(e).capitalize())
    except TypeError as e:
        raise ValueError('f has an unknown symbol. ' + str(e).capitalize())

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
        print("WARNING: [invalid expression]: ", type(
            exception).__name__, "-", str(exception))
        return False
    return True

def steffensen_method(expr, x0, n):
    """
    Steffensen Method
    
    This function is used to calculate some necessary values for other functions

    Arguments:

        exp {string} - polynomial whose solution must be found
        x0 {float, int} - initial value to start iterations
        n {int} - number of iterations

    Returns:

        xn {float} - root approximation
    """
    f = lambda x: eval(expr)
    itera = 0
    xn = x0

    while (itera <= n):
        df = misc.derivative(f, xn, dx=1e-6)
        df2 = misc.derivative(f, xn, n=2, dx=1e-6)
        div = (2 * df**2 - f(xn) * df2)

        xn = xn - (2 * f(xn) * df) / div
        itera += 1
    return xn

def plotFunction(k, error, title):
    """
    This function is used to plot iterations vs error

    Arguments:

        k {iterable} - an iterable with x axis values
        error {iterable} - an iterable with y axis values
        title {string} - plot title

    Returns:

        This function doesn't return
    """
    plt.title(title)
    plt.xlabel('Iterations k')
    plt.ylabel('Error |f(xk)|')
    plt.plot(k, error)
    plt.show()
