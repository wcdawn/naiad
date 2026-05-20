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
    fname_base = "assembly_mix_uniform.inp"
    start_snorder = 8
    start_refine = 0
    max_refine = 7

    fname_run = fname_base.replace(".inp", "_run.inp")

    runtxt = open(fname_base, "r").readlines()

    keff = np.zeros(max_refine)
    nx = np.zeros(max_refine, dtype=int)
    for i in range(max_refine):
        sn = start_snorder * 2**i
        r = start_refine + i

        inp = set_input(runtxt, sn, r)
        open(fname_run, "w").writelines(inp)

        out = run(executable, fname_run)
        if "WARNING" in out:
            print("Encountered warning sn=", sn, "r=", r)
        out = out.split("\n")
        keff[i] = get_keff(out)
        nx[i] = get_nx(out)
        print("sn=", sn, "r=", r, "keff= {:.16f}".format(keff[i]))

    for k in keff:
        print("{:.16f}".format(k))
