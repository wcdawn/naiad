import numpy as np
import matplotlib.pyplot as plt
import scipy
import scipy.optimize
import matplotlib
import mpmath

mpmath.mp.dps = 64
print(mpmath.mp)

matplotlib.rcParams["lines.linewidth"] = 2
# matplotlib.rcParams["mathtext.fontset"] = "stix"
# matplotlib.rcParams["font.family"] = "STIXGeneral"
matplotlib.rcParams["font.size"] = 16


def bisection_search(f, xlo, xhi, tol):
    xlo = mpmath.mpf(xlo)
    xhi = mpmath.mpf(xhi)
    tol = mpmath.mpf(tol)
    while (xhi - xlo) >= tol:
        flo = f(xlo)
        fhi = f(xhi)

        xmid = mpmath.mpf(0.5) * (xlo + xhi)
        fmid = f(xmid)

        if mpmath.sign(flo) == mpmath.sign(fmid):
            xlo = xmid
        else:
            xhi = xmid

    return mpmath.mpf(0.5) * (xhi + xlo)


if __name__ == "__main__":

    extension = "png"
    dpi = 600

    nusf = mpmath.mpf(0.02)
    sigma_r_fuel = mpmath.mpf(0.02)
    D_fuel = mpmath.mpf(1.2)

    sigma_r_refl = mpmath.mpf(0.015)
    D_refl = mpmath.mpf(0.7)

    kappa_refl = mpmath.sqrt(sigma_r_refl / D_refl)

    LF = mpmath.mpf(80.0)
    LR = mpmath.mpf(100.0)

    f = lambda bf: (D_refl * kappa_refl) / mpmath.tanh(
        kappa_refl * (LR - LF)
    ) - D_fuel * bf * mpmath.tan(bf * LF)
    BF = bisection_search(f, 0.017, 0.018, tol=1e-16)
    print("BF={:.20f}".format(float(BF)))

    keff = nusf / (D_fuel * BF**2 + sigma_r_fuel)
    print("keff={:.20f}".format(float(keff)))

    phi0 = 1.0

    phi_fuel = lambda x: phi0 * mpmath.cos(BF * x)
    phi_refl = (
        lambda x: phi0
        * mpmath.cos(BF * LF)
        * (
            mpmath.cosh(kappa_refl * (x - LF))
            - mpmath.sinh(kappa_refl * (x - LF)) / mpmath.tanh(kappa_refl * (LR - LF))
        )
    )

    def phi(x):
        p = np.zeros_like(x)
        for i in range(len(x)):
            if x[i] < LF:
                p[i] = float(phi_fuel(x[i]))
            else:
                p[i] = float(phi_refl(x[i]))
        return p

    x = np.linspace(0.0, LR, 1024)

    plt.figure()
    plt.plot(x, phi(x))
    plt.title("Flux")
    plt.xlabel("x [cm]")
    plt.ylabel("$\\phi(x)$")
    plt.tight_layout()
    plt.savefig("analytic_phi." + extension, dpi=dpi)

    plt.show()
