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
    x0 = 3 / 4
    tol = 0.000001
    f = "cos(2*x)^2-x^2"
    graf = 1

    print('\n\nSteffensen\'s: ')
    print(sne_fd_1(f, x0, tol, graf))
    print('\n\nHalley\'s: ')
    print(sne_ud_1(f, x0, tol, graf))
    print('\n\nFrontini\'s y Sormani\'s: ')
    print(sne_ud_2(f, x0, tol, graf))
    print('\n\n')

    # test_sne_ud_3()
    # test_sne_ud_4()
    # test_sne_ud_5()
    # test_sne_ud_6()
