import numpy as np
import subprocess
import matplotlib
import matplotlib.pyplot as plt
import sys

matplotlib.rcParams["lines.linewidth"] = 2
# matplotlib.rcParams["mathtext.fontset"] = "stix"
# matplotlib.rcParams["font.family"] = "STIXGeneral"
matplotlib.rcParams["font.size"] = 14


def set_input(txt, snorder, refine):
    out = []
    for line in txt:
        if "refine" in line:
            out.append("refine {:d}\n".format(refine))
        elif "snorder" in line:
            out.append("snorder {:d}\n".format(snorder))
        else:
            out.append(line)
    return out


def run(exe, inp):
    result = subprocess.run([exe, inp], stdout=subprocess.PIPE)
    return result.stdout.decode("utf-8")


def get_keff(lines):
    for line in lines:
        if "keff = " in line:
            line = line.split()
            return float(line[2])


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

    extension = "png"
    resolution = 600

    executable = "/Users/williamdawn/work/naiad/src/naiad.x"
    fname_base = "pin_slab_uo2.inp"
    arr_snorder = [2, 4, 8, 16, 32, 64, 128, 256]
    max_refine = 10

    fname_run = fname_base.replace(".inp", "_run.inp")

    runtxt = open(fname_base, "r").readlines()

    table = np.zeros((len(arr_snorder), max_refine))

    nx = np.zeros(max_refine, dtype=int)
    for snidx in range(len(arr_snorder)):
        sn = arr_snorder[snidx]
        for r in range(max_refine):
            inp = set_input(runtxt, sn, r)
            open(fname_run, "w").writelines(inp)

            out = run(executable, fname_run)
            if "WARNING" in out:
                print("Encountered warning sn=", sn, "r=", r)
                # sys.exit(1)
            out = out.split("\n")
            keff = get_keff(out)
            table[snidx, r] = keff

            # note: overwritten several times
            nx[r] = get_nx(out)

            print("sn=", sn, "r=", r, "keff= {:.16f}".format(keff))

        with open("result.csv", "w") as f:
            for snidx in range(table.shape[0]):
                f.write(" , S{:d}".format(arr_snorder[snidx]))
            f.write("\n")
            for ridx in range(table.shape[1]):
                f.write("r{:d}".format(ridx))
                for snidx in range(table.shape[0]):
                    f.write(" , {:.16f}".format(table[snidx, ridx]))
                f.write("\n")

    # perform Richardson extrapolation
    extrap_space = np.zeros(table.shape[0])
    # second order in space
    for i in range(table.shape[0]):
        extrap_space[i] = table[i, -1] + (table[i, -1] - table[i, -2]) / 3.0
    # diff is the difference
    diff = np.zeros_like(table)
    for i in range(table.shape[0]):
        diff[i, :] = extrap_space[i] - table[i, :]

    # dump table in LaTeX format
    with open("result.tex", "w") as f:
        for i in range(table.shape[0]):
            f.write(" & {{$S_{{{:d}}}$}}".format(arr_snorder[i]))
        f.write("\\\\\n")
        f.write("\\midrule\n")
        for j in range(table.shape[1]):
            f.write("r{:d}".format(j))
            for i in range(table.shape[0]):
                f.write(" & {:.6f}".format(table[i, j]))
            f.write("\\\\\n")
        f.write("\\midrule\n")
        f.write("Extrap.")
        for e in extrap_space:
            f.write(" & {:.6f}".format(e))
        f.write("\\\\\n")

    plt.figure()
    for snidx in range(table.shape[0]):
        plt.semilogx(
            nx, table[snidx, :], "-o", label="$S_{{{:d}}}$".format(arr_snorder[snidx])
        )
    plt.legend()
    plt.xlabel("NX")
    plt.ylabel("keff")
    plt.title("Spatial Refinement")
    plt.tight_layout()
    plt.savefig("refinement_space." + extension, dpi=resolution)

    plt.figure()
    for snidx in range(table.shape[0]):
        plt.loglog(
            nx,
            np.abs(diff[snidx, :]) * 1e5,
            "-o",
            label="$S_{{{:d}}}$".format(arr_snorder[snidx]),
        )
    plt.loglog(nx, 1e5 * nx.astype(float) ** (-2), "-k", lw=1, label="_hide")
    plt.legend()
    plt.xlabel("NX")
    plt.ylabel("|extrapolation - keff| [pcm]")
    plt.title("Spatial Refinement to Extrapolation")
    plt.tight_layout()
    plt.savefig("refinement_conv." + extension, dpi=resolution)

    snorder = np.array(arr_snorder)

    plt.figure()
    for ridx in range(table.shape[1]):
        plt.plot(snorder, table[:, ridx], "-o", label="NX={:d}".format(nx[ridx]))
    plt.legend()
    plt.xlabel("$S_N$")
    plt.gca().set_xticks(arr_snorder)
    plt.ylabel("keff")
    plt.title("Angular Refinement")
    plt.tight_layout()
    plt.savefig("refinement_moment." + extension, dpi=resolution)

    plt.show()
