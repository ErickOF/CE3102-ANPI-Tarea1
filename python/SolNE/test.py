from fd import *
from ud import *


def test_sne_ud_3():
    x0 = 3 / 4
    tol = 0.0000000000000001
    func = 'np.cos(2 * x)**2 - x**2'
    graf = 1

    xAprox, _iter = sne_ud_3(func, x0, tol, graf)
    print('xAprox = {}\nIteraciones = {}'.format(xAprox, _iter))

def test_sne_ud_4():
    x0 = 3 / 4
    tol = 0.0000000000000001
    func = 'np.cos(2 * x)**2 - x**2'
    graf = 1

    xAprox, _iter = sne_ud_4(func, x0, tol, graf)
    print('xAprox = {}\nIteraciones = {}'.format(xAprox, _iter))

def test_sne_ud_5():
    x0 = 3 / 4
    tol = 0.0000000000000001
    func = 'np.cos(2 * x)**2 - x**2'
    graf = 1

    xAprox, _iter = sne_ud_5(func, x0, tol, graf)
    print('xAprox = {}\nIteraciones = {}'.format(xAprox, _iter))

def test_sne_ud_6():
    x0 = 3 / 4
    tol = 0.0000000000000001
    func = 'np.cos(2 * x)**2 - x**2'
    graf = 1

    xAprox, _iter = sne_ud_6(func, x0, tol, graf)
    print('xAprox = {}\nIteraciones = {}'.format(xAprox, _iter))


if __name__ == '__main__':
    expression = '9*x+3'                  # Division by zero
    # expression = 'exp(x) - 3*x'
    x0 = 0.5
    tol = 0.000001

    print('\n\nSteffensen\'s: ')
    print(sne_fd_1(expression, x0, tol))
    print('\n\nHalley\'s: ')
    print(sne_ud_1(expression, x0, tol))
    print('\n\nFrontini\'s y Sormani\'s: ')
    print(sne_ud_2(expression, x0, tol))
    print('\n\n')

    test_sne_ud_3()
    test_sne_ud_4()
    test_sne_ud_5()
    test_sne_ud_6()
