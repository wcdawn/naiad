import numpy as np


class Problem:
    def __init__(self, number, name, kref=1.0):
        self.number = number
        self.name = name
        self.kref = kref
        self.keff = 0.0
        self.diff = 0.0

        s = self.name.split("-")

        self.ng = int(s[-3])
        self.pn = int(s[-2])

    def ngroup(self):
        return self.ng


cases = [
    Problem(1, "PUa-1-0-IN", 2.612903),
    Problem(2, "PUa-1-0-SL"),
    Problem(3, "PUa-H2O(1)-1-0-SL"),
    Problem(4, "PUa-H2O(0.5)-1-0-SL"),
    Problem(5, "PUb-1-0-IN", 2.290323),
    Problem(6, "PUb-1-0-SL"),
    Problem(7, "PUb-1-0-CY"),
    Problem(8, "PUb-1-0-SP"),
    Problem(9, "PUb-H2O(1)-1-0-CY"),
    Problem(10, "PUb-H2O(10)-1-0-CY"),
    Problem(11, "Ua-1-0-IN", 2.25),
    Problem(12, "Ua-1-0-SL"),
    Problem(13, "Ua-1-0-CY"),
    Problem(14, "Ua-1-0-SP"),
    Problem(15, "Ub-1-0-IN", 2.330917),
    Problem(16, "Ub-H2O(1)-1-0-SP"),
    Problem(17, "Uc-1-0-IN", 2.256083),
    Problem(18, "Uc-H2O(2)-1-0-SP"),
    Problem(19, "Ud-1-0-IN", 2.232667),
    Problem(20, "Ud-H2O(3)-1-0-SP"),
    Problem(21, "UD2O-1-0-IN", 1.1333333),
    Problem(22, "UD2O-1-0-SL"),
    Problem(23, "UD2O-1-0-CY"),
    Problem(24, "UD2O-1-0-SP"),
    Problem(25, "UD2O-H2O(1)-1-0-SL"),
    Problem(26, "UD2O-H2O(10)-1-0-SL"),
    Problem(27, "UD2O-H2O(1)-1-0-CY"),
    Problem(28, "UD2O-H2O(10)-1-0-CY"),
    Problem(29, "Ue-1-0-IN", 2.1806667),
    Problem(30, "Ue-Fe-Na-1-0-SL"),
    Problem(31, "PU-1-1-IN", 2.5),
    Problem(32, "PUa-1-1-SL"),
    Problem(33, "PUa-1-2-SL"),
    Problem(34, "PUb-1-1-SL"),
    Problem(35, "PUb-1-2-SL"),
    Problem(36, "Ua-1-1-CY"),
    Problem(37, "Ub-1-1-CY"),
    Problem(38, "UD2Oa-1-1-IN", 1.205587),
    Problem(39, "UD2Oa-1-1-SP"),
    Problem(40, "UD2Ob-1-1-IN", 1.227391),
    Problem(41, "UD2Ob-1-1-SP"),
    Problem(42, "UD2Oc-1-1-IN", 1.130933),
    Problem(43, "UD2Oc-1-1-SP"),
    Problem(44, "PU-2-0-IN", 2.683767),
    Problem(45, "PU-2-0-SL"),
    Problem(46, "PU-2-0-SP"),
    Problem(47, "U-2-0-IN", 2.216349),
    Problem(48, "U-2-0-SL"),
    Problem(49, "U-2-0-SP"),
    Problem(50, "UAL-2-0-IN", 2.661745),
    Problem(51, "UAL-2-0-SL"),
    Problem(52, "UAL-2-0-SP"),
    Problem(53, "URRa-2-0-IN", 1.631452),
    Problem(54, "URRa-2-0-SL"),
    Problem(55, "URRa-2-0-SP"),
    Problem(56, "URRb-2-0-IN", 1.365821),
    Problem(57, "URRc-2-0-IN", 1.633380),
    Problem(58, "URRb-H2Oa(1)-2-0-SL"),
    Problem(59, "URRb-H2Oa(5)-2-0-SL"),
    Problem(60, "URRb-H2Oa(IN)-2-0-SL"),
    Problem(61, "URRc-H2Oa(IN)-2-0-SL"),
    Problem(62, "URRd-2-0-IN", 1.034970),
    Problem(63, "URRd-H2Ob(1)-2-0-ISLC"),
    Problem(64, "URRd-H2Ob(10)-2-0-ISLC"),
    Problem(65, "URRd-H2Oc(1)-2-0-ISLC"),
    Problem(66, "URRd-H2Oc(10)-2-0-ISLC"),
    Problem(67, "UD2O-2-0-IN", 1.000196),
    Problem(68, "UD2O-2-0-SL"),
    Problem(69, "UD2O-2-0-SP"),
    Problem(70, "URRa-2-1-IN", 1.631452),
    Problem(71, "URRa-2-1-SL"),
    Problem(72, "UD2O-2-1-IN", 1.000196),
    Problem(73, "UD2O-2-1-SL"),
    Problem(74, "URR-3-0-IN", 1.60),
    Problem(75, "URR-6-0-IN", 1.60),
]


def get_fname(prob):
    fname = prob.name
    fname = fname.replace("(", "_")
    fname = fname.replace(")", "")
    fname = fname.replace(".", "")
    return fname


def print_csv(cases, fname):
    with open(fname, "w") as f:
        f.write("Number , Case , kref , keff , Diff. [pcm]\n")
        for case in cases:
            if case.keff != 0.0:
                f.write(
                    "{:d} , {:s} , {:.6f} , {:.6f} , {:.2f}\n".format(
                        case.number,
                        case.name,
                        case.kref,
                        case.keff,
                        case.diff,
                    )
                )
    return


def print_tex_onegroup(cases, fname):
    # NOTE: I have omitted the "header" for now
    # I expect that whoever is using this will write that themselves
    with open(fname, "w") as f:
        for case in cases:
            if (case.keff != 0.0) and (case.ngroup() == 1):
                f.write(
                    "{:d} & {:s} & {:.6f} & {:.6f} & {:.2f} \\\\\n".format(
                        case.number,
                        case.name,
                        case.kref,
                        case.keff,
                        case.diff,
                    )
                )
    return


def print_tex_multigroup(cases, fname):
    # NOTE: I have omitted the "header" for now
    # I expect that whoever is using this will write that themselves
    with open(fname, "w") as f:
        for case in cases:
            if (case.keff != 0.0) and (case.ngroup() != 1):
                f.write(
                    "{:d} & {:s} & {:.6f} & {:.6f} & {:.2f} \\\\\n".format(
                        case.number,
                        case.name,
                        case.kref,
                        case.keff,
                        case.diff,
                    )
                )
    return


def get_keff(fname):
    with open(fname, "r") as f:
        for line in f:
            if "keff = " in line:
                line = line.split()
                return float(line[2])


if __name__ == "__main__":

    evaluated = 0
    skipped = 0

    for case in cases:
        fname_stub = get_fname(case)
        fname_full = "./" + fname_stub + "/" + fname_stub + ".out"
        try:
            keff = get_keff(fname_full)
            print(case.name, fname_stub, "{:.6f}".format(keff))
            case.keff = keff
            evaluated += 1
        except FileNotFoundError:
            print("skipping case ", case.name)
            skipped += 1

    print("Evaluated: ", evaluated)
    print("Skipped: ", skipped)

    maxdiff = 0.0
    rmsdiff = 0.0
    count = 0

    for case in cases:
        if case.keff != 0.0:
            case.diff = (case.kref - case.keff) * 1e5
            maxdiff = np.max((maxdiff, np.abs(case.diff)))
            rmsdiff += case.diff**2
            count += 1

    print("Maximum difference {:.2f} [pcm]".format(maxdiff))
    rmsdiff = np.sqrt(rmsdiff / count)
    print("RMS difference {:.2f} [pcm]".format(rmsdiff))

    print_csv(cases, "summary.csv")
    print_tex_onegroup(cases, "summary_onegroup.tex")
    print_tex_multigroup(cases, "summary_multigroup.tex")
