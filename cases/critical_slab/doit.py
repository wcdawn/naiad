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
        elif "nx= " in line:
            nx = int(line.split()[1])
            if ready:
                return nx
    # in case refinement not performed
    return nx


if __name__ == "__main__":

    executable = "/Users/williamdawn/work/naiad/src/naiad.x"
    fname_base = "slab.inp"
    max_snorder = 64
    max_refine = 16

    fname_run = fname_base.replace(".inp", "_run.inp")

    runtxt = open(fname_base, "r").readlines()

    table = np.zeros((int(max_snorder / 2), max_refine))
    nx = np.zeros(max_refine, dtype=int)

    sidx = 0
    for sn in range(2, max_snorder + 1, 2):
        ridx = 0
        for r in range(0, max_refine):
            inp = set_input(runtxt, sn, r)
            open(fname_run, "w").writelines(inp)

            out = run(executable, fname_run)
            if "WARNING" in out:
                print("Encountered warning sn=", sn, "r=", r)
                sys.exit(1)
            out = out.split("\n")
            keff = get_keff(out)
            table[sidx, ridx] = keff

            # note: overwritten several times
            nx[ridx] = get_nx(out)

            print("sn=", sn, "r=", r, "keff= {:.16f}".format(keff))
            ridx += 1
        sidx += 1

    with open("result.csv", "w") as f:
        for sidx in range(table.shape[0]):
            f.write(" , S{:d}".format((sidx + 1) * 2))
        f.write("\n")
        for ridx in range(table.shape[1]):
            f.write("r{:d}".format(ridx))
            for sidx in range(table.shape[0]):
                f.write(" , {:.16f}".format(table[sidx, ridx]))
            f.write("\n")

    diff = 1.0 - table

    plt.figure()
    for sidx in range(table.shape[0]):
        plt.semilogx(
            nx, diff[sidx, :] * 1e5, "-o", label="$S_{{{:d}}}$".format(2 * (sidx + 1))
        )
    plt.legend()
    plt.xlabel("NX")
    plt.ylabel("Error [pcm]")
    plt.title("Spatial Refinement")
    plt.tight_layout()

    plt.figure()
    for sidx in range(table.shape[0]):
        plt.semilogx(
            nx, table[sidx, :], "-o", label="$S_{{{:d}}}$".format(2 * (sidx + 1))
        )
    plt.legend()
    plt.xlabel("NX")
    plt.ylabel("keff")
    plt.title("Spatial Refinement")
    plt.tight_layout()

    snorder = np.zeros(table.shape[0])
    for sidx in range(table.shape[0]):
        snorder[sidx] = 2 * (sidx + 1)

    plt.figure()
    for ridx in range(table.shape[1]):
        plt.plot(snorder, diff[:, ridx] * 1e5, "-o", label="NX={:d}".format(nx[ridx]))
    plt.legend()
    plt.xlabel("$S_N$")
    plt.ylabel("Error [pcm]")
    plt.title("Angle Refinement")
    plt.tight_layout()

    plt.figure()
    for ridx in range(table.shape[1]):
        plt.plot(snorder, table[:, ridx], "-o", label="NX={:d}".format(nx[ridx]))
    plt.legend()
    plt.xlabel("$S_N$")
    plt.ylabel("keff")
    plt.title("Angle Refinement")
    plt.tight_layout()

    plt.show()
