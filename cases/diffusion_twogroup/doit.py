import numpy as np
import subprocess
import matplotlib
import matplotlib.pyplot as plt
import sys

matplotlib.rcParams["lines.linewidth"] = 2
# matplotlib.rcParams["mathtext.fontset"] = "stix"
# matplotlib.rcParams["font.family"] = "STIXGeneral"
matplotlib.rcParams["font.size"] = 16


def set_input(txt, refine):
    out = []
    for line in txt:
        if "refine" in line:
            out.append("refine {:d}\n".format(refine))
        else:
            out.append(line)
    return out


def run(exe, inp):
    result = subprocess.run([exe, inp], stdout=subprocess.PIPE)
    return result.stdout.decode("utf-8")


def get_keyword(lines, keyword):
    for line in lines:
        if (keyword + " = " in line) or (keyword + "  = " in line):
            line = line.split()
            return float(line[2])


def get_keff(lines):
    return get_keyword(lines, "keff")


def get_nx(lines):
    ready = False
    for line in lines:
        if "after refinement" in line:
            ready = True
        elif "nx=" in line:
            nx = int(line.split()[1])
            if ready:
                return nx
    # in case refinement not performed
    return nx


if __name__ == "__main__":

    executable = "/Users/williamdawn/work/naiad/src/naiad.x"
    fname_base = "twogroup.inp"
    max_refine = 9

    fname_run = fname_base.replace(".inp", "_run.inp")

    runtxt = open(fname_base, "r").readlines()

    linferr = np.zeros((max_refine, 2))
    keff_diff = np.zeros(max_refine)
    nx = np.zeros(max_refine, dtype=int)

    for r in range(0, max_refine):
        inp = set_input(runtxt, r)
        open(fname_run, "w").writelines(inp)

        out = run(executable, fname_run)
        if "CONVERGENCE!" not in out:
            print("failed to converge r=", r)
            sys.exit(1)
        out = out.split("\n")

        keff = get_keff(out)
        keff_diff[r] = get_keyword(out, "keff_diff")
        for g in range(linferr.shape[1]):
            linferr[r, g] = get_keyword(out, "linf_g{:d}".format(g))

        nx[r] = get_nx(out)

        print("r=", r, "keff= {:.16f}".format(keff))

    with open("result.csv", "w") as f:
        f.write(" , keff_diff [pcm]")
        for g in range(linferr.shape[1]):
            f.write(" , linferr_g{:d}".format(g))
        f.write("\n")
        for ridx in range(max_refine):
            f.write("r{:d}".format(ridx))
            f.write(" , {:.16f}".format(keff_diff[ridx]))
            for g in range(linferr.shape[1]):
                f.write(" , {:.16e}".format(linferr[ridx, g]))
            f.write("\n")

    plt.figure()
    plt.loglog(nx, np.abs(keff_diff), "-o")
    plt.loglog(nx, np.abs(keff_diff[0]) / nx[0] / nx**2, "-k", lw=1, label="_hide")
    plt.xlabel("$N_X$")
    plt.ylabel("keff error [pcm]")
    plt.title("Spatial Refinement")
    plt.tight_layout()

    plt.figure()
    for g in range(linferr.shape[1]):
        plt.loglog(nx, linferr[:, g], "-o", label="g={:d}".format(g))
    plt.loglog(nx, np.abs(linferr[0, 0]) / nx[0] / nx**2, "-k", lw=1, label="_hide")
    plt.legend()
    plt.xlabel("$N_X$")
    plt.ylabel("$\\| \\phi - \\phi_{\\text{exact}} \\|$")
    plt.title("Spatial Refinement")
    plt.tight_layout()

    plt.show()
