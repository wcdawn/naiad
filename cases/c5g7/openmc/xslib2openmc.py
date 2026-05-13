import sys
import xslib
import openmc
import openmc.mgxs
import numpy as np


def dump(xs, fname):

    ngroup = len(xs[list(xs.keys())[0]]["sigma_t"])
    print("ngroup=", ngroup)

    groups = openmc.mgxs.EnergyGroups(np.logspace(-5, 7, ngroup + 1))
    nubar = 2.0  # arbitrarily set nu=2 everywhere since I only store nusf

    xsdata = []
    for mat in xs:
        xsthis = openmc.XSdata(mat, groups)
        xsthis.order = xs[mat]["scatter"].shape[0] - 1

        xsthis.set_total(xs[mat]["sigma_t"])

        # absorption is required by OpenMC
        if "sigma_a" in xs[mat]:
            xsthis.set_absorption(xs[mat]["sigma_a"])
        else:
            sigma_a = xs[mat]["sigma_t"].copy()
            for g in range(ngroup):
                sigma_a[g] -= np.sum(xs[mat]["scatter"][0, :, g])
            xsthis.set_absorption(sigma_a)

        scatter_matrix = xs[mat]["scatter"].copy()
        for ell in range(xs[mat]["scatter"].shape[0]):
            scatter_matrix[ell, :, :] = scatter_matrix[ell, :, :].transpose()
        scatter_matrix = np.rollaxis(scatter_matrix, 0, 3)
        xsthis.set_scatter_matrix(scatter_matrix)

        if "nusf" in xs[mat]:
            xsthis.set_nu_fission(xs[mat]["nusf"])
            xsthis.set_chi(xs[mat]["chi"])

        xsdata.append(xsthis)

    mgxs_file = openmc.MGXSLibrary(groups)
    mgxs_file.add_xsdatas(xsdata)
    mgxs_file.export_to_hdf5(fname)


if __name__ == "__main__":

    fname = sys.argv[1]
    xs = xslib.load(fname)

    dump(xs, fname.replace(".xs", ".h5"))
