import numpy as np
import matplotlib.pyplot as plt

if __name__ == "__main__":

    fname = "analytic_g0.csv"

    dat = np.loadtxt(fname, delimiter=",", skiprows=1)
    x = dat[:,0]
    exact = dat[:,1]
    calculated = dat[:,2]

    plt.figure()
    plt.plot(x, exact, label="exact")
    plt.plot(x, calculated, label="calculated")
    plt.legend()
    plt.xlabel("x [cm]")
    plt.ylabel("phi")
    plt.title("Result")
    plt.tight_layout()

    plt.figure()
    plt.plot(x, exact-calculated)
    plt.xlabel("x [cm]")
    plt.ylabel("phi diff")
    plt.title("Absolute Difference")
    plt.tight_layout()

    plt.figure()
    plt.plot(x, (exact-calculated) / exact)
    plt.xlabel("x [cm]")
    plt.ylabel("phi diff")
    plt.title("Relative Difference")
    plt.tight_layout()

    plt.show()
