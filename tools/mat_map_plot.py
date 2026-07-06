import numpy as np
import matplotlib.pyplot as plt
import sys

import plot_settings


def load(fname):
    dat = np.loadtxt(fname, delimiter=",", skiprows=1, dtype=str)
    xcenter = dat[:, 0].astype(float)
    mat_map = dat[:, 1]

    mat_names = set(list(mat_map))
    mat_dict = {}
    i = 0
    for mname in mat_names:
        mat_dict[mname] = i
        i += 1

    mat_idx = np.zeros_like(mat_map, dtype=int)
    for i in range(len(mat_idx)):
        mat_idx[i] = mat_dict[mat_map[i]]

    return xcenter, mat_idx, mat_dict


if __name__ == "__main__":

    extension = "pdf"
    resolution = 600

    fname = sys.argv[1]

    xcenter, mat_idx, mat_dict = load(fname)

    plt.figure()
    for mname in mat_dict:
        m = mat_idx[mat_idx == mat_dict[mname]]
        x = xcenter[mat_idx == mat_dict[mname]]
        plt.plot(x, m, "o", label=mname)
    plt.legend()
    plt.xlabel("x [cm]")
    plt.ylabel("Material Index")
    plt.title("Naiad Material Map")
    plt.tight_layout()
    plt.savefig(fname.replace(".csv", "." + extension), dpi=resolution)

    plt.show()
