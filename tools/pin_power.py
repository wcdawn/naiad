import numpy as np
import sys


def compute_dx(x):
    dx = np.zeros_like(x)
    xleft = 0.0
    for i in range(len(dx)):
        dx[i] = 2.0 * (x[i] - xleft)
        xleft += dx[i]
    return dx


if __name__ == "__main__":

    fname = sys.argv[1]
    pin_pitch = 1.26

    dat = np.loadtxt(fname, delimiter=",", skiprows=1)

    x = dat[:, 0]
    power = dat[:, 1]

    dx = compute_dx(x)

    pin_edges = [0.0]
    xmax = x[-1] + 0.5 * dx[-1]
    npin = 0
    while pin_edges[-1] < xmax:
        pin_edges.append(pin_pitch * (npin + 1))
        npin += 1

    pin_power = np.zeros(npin)
    for i in range(npin):
        f1 = x > pin_edges[i]
        f2 = x < pin_edges[i + 1]
        filt = f1 & f2
        pin_power[i] = np.sum(power[filt] * dx[filt])

    xsum = np.sum(pin_power)
    npin_hot = np.sum(pin_power > 0)
    xmean = xsum / npin_hot
    pin_power /= xmean

    count = 0
    for i in range(npin):
        if pin_power[i] > 0.0:
            count += 1
            # print("{:d} {:.3f}".format(count, pin_power[i]))
            print("{:.3f}".format(pin_power[i]))
