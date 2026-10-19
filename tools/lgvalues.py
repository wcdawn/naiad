import numpy as np

from lgvalues_weights import legendre_weights
from lgvalues_abscissa import legendre_roots


def print_cpp(NMAX, fname):
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


def print_f90(NMAX, fname):
    open(fname, "w")  # clear file
    for n in range(NMAX):
        xi = legendre_roots[n + 1]
        wi = legendre_weights[n + 1]
        idx = np.argsort(xi)
        xi = xi[idx]
        wi = wi[idx]
        with open(fname, "a") as f:
            f.write("      case ({:d})\n".format(n + 1))
            f.write("        xi = (/ &\n")
            for i in range(len(xi)):
                x = xi[i]
                if i < len(xi) - 1:
                    f.write("          {:23.16e}_rk, &\n".format(x))
                else:
                    f.write("          {:23.16e}_rk  &\n".format(x))
            f.write("        /)\n")
            f.write("        wi = (/ &\n")
            for i in range(len(wi)):
                x = wi[i]
                if i < len(wi) - 1:
                    f.write("          {:23.16e}_rk, &\n".format(x))
                else:
                    f.write("          {:23.16e}_rk  &\n".format(x))
            f.write("        /)\n")


if __name__ == "__main__":

    NMAX = 64
    fname = "gauss_legendre.txt"

    print_cpp(NMAX, fname)
    print_f90(NMAX, fname)
