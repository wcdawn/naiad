import numpy as np
import matplotlib.pyplot as plt
import sys

import plot_settings

if __name__ == "__main__":

    fname = sys.argv[1]
    extension = "png"
    resolution = 600

    dat = np.loadtxt(fname, delimiter=",", skiprows=1)

    x = dat[:, 0]
    flux = dat[:, 1:]

    print(x)

    plt.figure()
    for g in range(flux.shape[1]):
        plt.plot(x, flux[:, g], label="g={:d}".format(g))
    if flux.shape[1] <= 10:
        plt.legend()
    plt.xlabel("x [cm]")
    plt.ylabel("Flux")
    plt.title("NAIAD Flux")
    plt.tight_layout()
    plt.savefig("flux." + extension, dpi=resolution)

    plt.show()
