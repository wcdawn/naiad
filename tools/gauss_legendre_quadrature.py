import numpy as np
import scipy
import scipy.special
import scipy.optimize
import matplotlib.pyplot as plt
import sys

import mpmath

# TODO probably too fine
mpmath.mp.dps = 64
print(mpmath.mp)

import plot_settings


def rewrap(x):
    y = []
    for xi in x:
        y.append(mpmath.mpf(xi))
    return y


def bisect(f, a, b, tol):
    xlo = a
    xhi = b
    flo = f(xlo)
    fhi = f(xhi)
    if mpmath.sign(flo) == mpmath.sign(fhi):
        print("Bisection failed, same sign.")
        sys.exit(1)
    while True:
        flo = f(xlo)
        fhi = f(xhi)
        xmid = mpmath.mpf(0.5) * (xlo + xhi)
        fmid = f(xmid)
        if mpmath.fabs(fmid) < tol:
            return xmid
        if mpmath.sign(flo) == mpmath.sign(fmid):
            xlo = xmid
        else:
            xhi = xmid


def legendre_eval(n, x):
    if n == 0:
        return mpmath.mpf(1.0)
    elif n == 1:
        return mpmath.mpf(x)
    else:
        return (
            (2 * n - 1) * x * legendre_eval(n - 1, x)
            - (n - 1) * legendre_eval(n - 2, x)
        ) / n


# first derivative of Pn(x)
def deriv_legendre_eval(n, x):
    if n == 0:
        return mpmath.mpf(0)
    elif n == 1:
        return mpmath.mpf(1)
    else:
        return n * legendre_eval(n - 1, x) + x * deriv_legendre_eval(n - 1, x)


def legendre_zeros_old(n):
    xi, wi = scipy.special.roots_legendre(n)
    return rewrap(xi)


def legendre_zeros(n):
    xi = []

    if (n + 1) % 2 == 0:
        xi.append(mpmath.mpf(0.0))
    if n == 1:
        return xi

    f = lambda x: legendre_eval(n, x)
    atol = mpmath.mpf(1e-32)

    # NOTE: this is actually too fine!
    # there are only n/2 zeros in this range since we're searching the half-range
    equal_width = mpmath.mpf(1.0) / n
    regions = [mpmath.mpf(0.0)]
    for i in range(n):
        regions.append((i + 1) * equal_width)
    regions[-1] = mpmath.mpf(1.0)
    regions[0] = atol

    for i in range(n):
        a = regions[i]
        b = regions[i + 1]
        if np.sign(f(a)) == np.sign(f(b)):
            continue
        x = bisect(f, a, b, atol)
        if x == mpmath.mpf(0):
            continue
        xi.append(x)

    for x in xi:
        if x > mpmath.mpf(0):
            xi.append(-x)

    if len(xi) != n:
        print("n=", n, "found=", len(xi))
        print("xi=", xi)
        print("Failed to find all roots")
        sys.exit(1)

    xi = sorted(xi)
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
    # note: only need to evaluate half
    for i in range(int((len(xi) + 1) / 2)):
        xii = xi[i]
        wi.append(
            mpmath.mpf(2)
            / ((mpmath.mpf(1) - xii**2) * deriv_legendre_eval(n, xii) ** 2)
        )
    flip = list(reversed(wi))
    if n % 2 == 1:
        flip = flip[1:]
    wi.extend(flip)
    return wi


if __name__ == "__main__":

    NMAX = 8
    fname = "gauss_legendre.txt"

    open(fname, "w")
    for n in range(NMAX):
        xi = legendre_zeros(n + 1)
        wi = legendre_weights(n + 1, xi)
        with open(fname, "a") as f:
            f.write("  std::vector<Quadrature_point>{{ // n = {:d}\n".format(n + 1))
            for i in range(len(xi)):
                f.write(
                    "    {{.x = {:23.16e}, .w = {:.16e}}},\n".format(
                        float(xi[i]), float(wi[i])
                    )
                )
            f.write("  },\n")

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
