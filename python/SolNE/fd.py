# Importante: ¡debe instalarse!
# https://pypi.org/project/Equation/#description
from Equation import Expression


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

    # ========================== Validations =================================
    if (validator(expr, x0, tol) != True):
        return
    # ======================= Local variables ================================
    ITER_LIMIT = 1000
    graph = 0
    itera = 0
    xn = x0
    error = abs(f(xn))
    xNext = 0

    try:
        # ================== Steffensen's Method =============================
        while (error > tol):
            if(itera >= ITER_LIMIT):
                print("WARNING: Iteration limit reached")
                break
            xNext = xn - (f(xn) * f(xn))/(f(xn + f(xn)) - f(xn))
            xn = xNext
            error = abs(f(xn))
            itera += 1

        graph = 1                           # Graph can be displayed

    except:
        print("WARNING: Math error")

    return xn, itera, graph


def sne_fd_2(*args, **kwargs):
    pass


def sne_fd_3(*args, **kwargs):
    pass


def sne_fd_4(*args, **kwargs):
    pass


def sne_fd_5(*args, **kwargs):
    pass


def sne_fd_6(*args, **kwargs):
    pass


def validator(expr, x0, tol):
    """Auxiliary function that validates inputs

    Arguments:
        expr {string} -- polynomial whose solution must be found
        x0 {float, int} -- initial value to start iterations
        tol {float, int} -- tolerance that indicates the stop condition

    Returns:
        flag {boolean} -- true if all inputs are valid

    """

    # ========================== Validations =================================
    if (type(expr) != str):
        print("WARNING: invalid expression")
        return False
    if (type(x0) != float and type(x0) != int):
        print("WARNING: x0 must be a number")
        return False
    if (type(tol) != float and type(tol) != int):
        print("WARNING: tolerance must be a number")
        return False
    try:
        f = Expression(expr, ["x"])
    except:
        print("WARNING: invalid expression")
        return False
    return True
