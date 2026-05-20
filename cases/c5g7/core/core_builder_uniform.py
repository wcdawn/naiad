import sys

if __name__ == "__main__":

    mat_uo2 = [
        "moderator",
        "uo2",
        "uo2",
        "uo2",
        "uo2",
        "uo2",
        "uo2",
        "uo2",
        "uo2",
        "uo2",
        "uo2",
        "uo2",
        "uo2",
        "moderator",
    ]
    dx_uo2 = [
        0.09,
        0.09,
        0.09,
        0.09,
        0.09,
        0.09,
        0.09,
        0.09,
        0.09,
        0.09,
        0.09,
        0.09,
        0.09,
        0.09,
    ]

    mat_mox43 = [
        "moderator",
        "mox43",
        "mox43",
        "mox43",
        "mox43",
        "mox43",
        "mox43",
        "mox43",
        "mox43",
        "mox43",
        "mox43",
        "mox43",
        "mox43",
        "moderator",
    ]
    dx_mox43 = [
        0.09,
        0.09,
        0.09,
        0.09,
        0.09,
        0.09,
        0.09,
        0.09,
        0.09,
        0.09,
        0.09,
        0.09,
        0.09,
        0.09,
    ]

    mat_mox70 = [
        "moderator",
        "mox70",
        "mox70",
        "mox70",
        "mox70",
        "mox70",
        "mox70",
        "mox70",
        "mox70",
        "mox70",
        "mox70",
        "mox70",
        "mox70",
        "moderator",
    ]
    dx_mox70 = [
        0.09,
        0.09,
        0.09,
        0.09,
        0.09,
        0.09,
        0.09,
        0.09,
        0.09,
        0.09,
        0.09,
        0.09,
        0.09,
        0.09,
    ]

    mat_mox87 = [
        "moderator",
        "mox87",
        "mox87",
        "mox87",
        "mox87",
        "mox87",
        "mox87",
        "mox87",
        "mox87",
        "mox87",
        "mox87",
        "mox87",
        "mox87",
        "moderator",
    ]
    dx_mox87 = [
        0.09,
        0.09,
        0.09,
        0.09,
        0.09,
        0.09,
        0.09,
        0.09,
        0.09,
        0.09,
        0.09,
        0.09,
        0.09,
        0.09,
    ]

    mat_gdt = [
        "moderator",
        "guide_tube",
        "guide_tube",
        "guide_tube",
        "guide_tube",
        "guide_tube",
        "guide_tube",
        "guide_tube",
        "guide_tube",
        "guide_tube",
        "guide_tube",
        "guide_tube",
        "guide_tube",
        "moderator",
    ]
    dx_gdt = [
        0.09,
        0.09,
        0.09,
        0.09,
        0.09,
        0.09,
        0.09,
        0.09,
        0.09,
        0.09,
        0.09,
        0.09,
        0.09,
        0.09,
    ]

    mat_fiss = [
        "moderator",
        "fission_chamber",
        "fission_chamber",
        "fission_chamber",
        "fission_chamber",
        "fission_chamber",
        "fission_chamber",
        "fission_chamber",
        "fission_chamber",
        "fission_chamber",
        "fission_chamber",
        "fission_chamber",
        "fission_chamber",
        "moderator",
    ]
    dx_fiss = [
        0.09,
        0.09,
        0.09,
        0.09,
        0.09,
        0.09,
        0.09,
        0.09,
        0.09,
        0.09,
        0.09,
        0.09,
        0.09,
        0.09,
    ]

    mat_mod = [
        "moderator",
        "moderator",
        "moderator",
        "moderator",
        "moderator",
        "moderator",
        "moderator",
        "moderator",
        "moderator",
        "moderator",
        "moderator",
        "moderator",
        "moderator",
        "moderator",
    ]
    dx_mod = [
        0.09,
        0.09,
        0.09,
        0.09,
        0.09,
        0.09,
        0.09,
        0.09,
        0.09,
        0.09,
        0.09,
        0.09,
        0.09,
        0.09,
    ]

    assembly_uo2 = [
        "uo2",
        "uo2",
        "gdt",
        "uo2",
        "uo2",
        "gdt",
        "uo2",
        "uo2",
        "gdt",
        "uo2",
        "uo2",
        "gdt",
        "uo2",
        "uo2",
        "gdt",
        "uo2",
        "uo2",
    ]

    assembly_mox = [
        "mox43",
        "mox70",
        "gdt",
        "mox87",
        "mox87",
        "gdt",
        "mox87",
        "mox87",
        "gdt",
        "mox87",
        "mox87",
        "gdt",
        "mox87",
        "mox87",
        "gdt",
        "mox70",
        "mox43",
    ]

    assembly_mod = [
        "mod",
        "mod",
        "mod",
        "mod",
        "mod",
        "mod",
        "mod",
        "mod",
        "mod",
        "mod",
        "mod",
        "mod",
        "mod",
        "mod",
        "mod",
        "mod",
        "mod",
    ]

    core = assembly_mox + assembly_uo2 + assembly_mod

    dx = []
    mat = []
    for pin in core:
        if pin == "uo2":
            dx += dx_uo2
            mat += mat_uo2
        elif pin == "mox43":
            dx += dx_mox43
            mat += mat_mox43
        elif pin == "mox70":
            dx += dx_mox70
            mat += mat_mox70
        elif pin == "mox87":
            dx += dx_mox87
            mat += mat_mox87
        elif pin == "gdt":
            dx += dx_gdt
            mat += mat_gdt
        elif pin == "mod":
            dx += dx_mod
            mat += mat_mod
        elif pin == "fiss":
            dx += dx_fiss
            mat += mat_fiss
        else:
            print("unknown pin: ", pin)
            sys.exit(1)

    print("nx ", len(dx))
    print("dx")
    for x in dx:
        print("{:.2f}".format(x))
    print("mat_map")
    for m in mat:
        print(m)
