import numpy as np

from lgvalues_weights import legendre_weights
from lgvalues_abscissa import legendre_roots


if __name__ == "__main__":

    NMAX = 64
    fname = "gauss_legendre.txt"

    open(fname, "w")
    for n in range(NMAX):
        xi = legendre_roots[n + 1]
        wi = legendre_weights[n + 1]
        idx = np.argsort(xi)
        xi = xi[idx]
        wi = wi[idx]
        with open(fname, "a") as f:
            f.write("  std::vector<Quadrature_point>{{ // n = {:d}\n".format(n + 1))
            for i in range(len(xi)):
                f.write(
                    "    {{.x = {:23.16e}, .w = {:.16e}}},\n".format(
                        float(xi[i]), float(wi[i])
                    )
                )
            f.write("  },\n")
