import numpy as np
import scipy
import scipy.optimize
import subprocess


def get_keff(fname_out):
    with open(fname_out, "r") as f:
        for line in f:
            if "keff = " in line:
                line = line.split()
                return float(line[2])


def run(exe, inp):
    result = subprocess.run([exe, inp], stdout=subprocess.PIPE)
    return result.stdout.decode("utf-8")


def make_input(fname_inp, snorder, refinement, width):
    lines = open(fname_inp, "r").readlines()
    out = []
    for line in lines:
        if "snorder" in line:
            out.append("snorder {:d}\n".format(snorder))
        elif "refine" in line:
            out.append("refine {:d}\n".format(refinement))
        elif "dx" in line:
            out.append("dx {:.8e}\n".format(width))
        else:
            out.append(line)
    fname_run = fname_inp.replace(".inp", "_run.inp")
    open(fname_run, "w").writelines(out)
    return fname_run


def do_it(exe, fname_inp, snorder, refinement, width):
    if hasattr(width, "__len__"):
        width = width[0]
    fname_run = make_input(fname_inp, snorder, refinement, width)
    data = run(exe, fname_run)
    fname_out = fname_run.replace(".inp", ".out")
    keff = get_keff(fname_out)
    return keff


def bisection(f, xlo, xhi, ftarget, ftol):
    lo = xlo
    hi = xhi
    fcalc = ftarget + 2.0 * ftol
    while abs(fcalc - ftarget) > ftol:
        flo = f(lo)
        fhi = f(hi)
        fcalc = 0.5 * (fhi + flo)
        mid = 0.5 * (lo + hi)
        if ftarget >= fcalc:
            hi = hi
            lo = mid
        else:
            hi = mid
            lo = lo
    return 0.5 * (hi + lo)


if __name__ == "__main__":

    myexec = "/Users/williamdawn/work/naiad/src/naiad.x"
    fname_base = "slab.inp"

    k_search = 1.0
    k_tol = 1e-8

    # bounding guesses
    xlo = 1.0
    xhi = 5.0

    search_method = "fsolve"  # "fsolve" or "bisection"

    snmax = 32
    refine_max = 10
    result = np.zeros((refine_max, int(snmax / 2)))

    snidx = 0
    for sn in range(2, snmax + 1, 2):
        for refine in range(refine_max):

            if search_method == "bisection":
                f = lambda x: do_it(myexec, fname_base, sn, refine, x)
                width = bisection(f, xlo, xhi, k_search, k_tol)
            elif search_method == "fsolve":
                f = lambda x: do_it(myexec, fname_base, sn, refine, x) - k_search
                calc = scipy.optimize.fsolve(f, 0.5 * (xlo + xhi), xtol=k_tol)
                width = calc[0]

            print("refine", refine, "sn", sn, "width", width)
            result[refine, snidx] = width

        snidx += 1

    with open("result.csv", "w") as f:
        s = "refinement"
        for sn in range(2, snmax + 1, 2):
            s += "," + "S{:d}".format(sn)
        f.write(s + "\n")
        for refine in range(refine_max):
            s = "r{:0d}".format(refine)
            for snidx in range(result.shape[1]):
                s += ",{:.8e}".format(result[refine,snidx])
            f.write(s + "\n")
