import numpy as np
import scipy
import scipy.special
import matplotlib.pyplot as plt

import mpmath

mpmath.mp.dps = 64
print(mpmath.mp)

import plot_settings


def rewrap(x):
    y = []
    for xi in x:
        y.append(mpmath.mpf(xi))
    return y


# it is noted in the documentation that
# return scipy.special.legendre_p(n, x, diff_n=0)
# would be equivalent, but it may be numerically unstable for values n >= 20
def legendre_eval(n, x):
    return mpmath.mpf(scipy.special.eval_legendre(n, float(x)))


# first derivative of Pn(x)
def deriv_legendre_eval(n, x):
    return mpmath.mpf(scipy.special.legendre_p(n, float(x), diff_n=1)[-1])


def legendre_zeros(n):
    xi, wi = scipy.special.roots_legendre(n)
    return rewrap(xi)


# I found this formula on the Wikipedia page for Gauss-Legendre quadrature rules.
# There is a detailed citation pointing to p. 887 of Abramowitz & Stegun.
#
# Abramowitz, Milton; Stegun, Irene Ann, eds. (1983) [June 1964].
# "Chapter 25.4, Integration". Handbook of Mathematical Functions with Formulas,
# Graphs, and Mathematical Tables. Applied Mathematics Series.
# Vol. 55 (Ninth reprint with additional corrections of tenth original printing with corrections (December 1972);
# first ed.). Washington D.C.; New York:
# United States Department of Commerce, National Bureau of Standards;
# Dover Publications.
# ISBN 978-0-486-61272-0. LCCN 64-60036. MR 0167642. LCCN 65-12253.
#
# I found a similar citation on my favorite Gauss-Legendre data page
# https://pomax.github.io/bezierinfo/legendre-gauss.html
# where the author points to Wolfram MathWorld
# https://mathworld.wolfram.com/Legendre-GaussQuadrature.html
#
# It seems like I have some learning to do here...
#
def legendre_weights(n, xi):
    wi = []
    for xii in xi:
        wi.append(
            mpmath.mpf(2)
            / ((mpmath.mpf(1) - xii**2) * deriv_legendre_eval(n, xii) ** 2)
        )
    return wi


if __name__ == "__main__":

    NMAX = 8
    N = range(2, NMAX + 1, 2)
    for n in N:
        xi = legendre_zeros(n)
        wi = legendre_weights(n, xi)
        print("n=", n)
        print(np.array(xi))
        print(np.array(wi))

    x = np.linspace(-1.0, 1.0, 1024)

    plt.figure()
    for ell in range(7):
        y = np.zeros_like(x)
        for i in range(len(x)):
            y[i] = float(legendre_eval(ell, x[i]))
        plt.plot(x, y, label="$P_{:d}$".format(ell))
    plt.legend()
    plt.xlabel("$x$")
    plt.ylabel("$P_n(x)$")
    plt.title("Legendre Polynomials")
    plt.tight_layout()

    plt.show()
