import sys
import openmc
import h5py

if __name__ == "__main__":

    mat_uo2 = ["moderator", "uo2", "moderator"]
    dx_uo2 = [0.09, 1.08, 0.09]

    mat_mox43 = ["moderator", "mox43", "moderator"]
    dx_mox43 = [0.09, 1.08, 0.09]

    mat_mox70 = ["moderator", "mox70", "moderator"]
    dx_mox70 = [0.09, 1.08, 0.09]

    mat_mox87 = ["moderator", "mox87", "moderator"]
    dx_mox87 = [0.09, 1.08, 0.09]

    mat_gdt = ["moderator", "guide_tube", "moderator"]
    dx_gdt = [0.09, 1.08, 0.09]

    mat_mod = ["moderator"]
    dx_mod = [1.26]

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

    model = openmc.Model()

    fname_xs = "../c5g7.h5"
    # For every cross section data set in the library, assign an openmc.Macroscopic object to a material
    material_dict = {}
    with h5py.File(fname_xs, "r") as h5:
        for xs in h5:
            xsname = xs
            if xsname == "moderator":
                xsname = "water"
            elif xsname == "mox70":
                xsname = "mox7"
            material_dict[xsname] = openmc.Material(name=xs)
            material_dict[xsname].set_density("macro", 1.0)
            material_dict[xsname].add_macroscopic(xs)

    # Instantiate a Materials collection, register all Materials, and export to XML
    materials = openmc.Materials(material_dict.values())

    # Set the location of the cross sections file to our pre-written set
    materials.cross_sections = fname_xs
    materials.export_to_xml()

    surf = []
    cell = []

    surf.append(openmc.XPlane(0.0))
    surf[0].boundary_type = "reflective"

    xx = 0.0

    for i in range(len(dx)):
        xx += dx[i]
        surf.append(openmc.XPlane(xx))

        c = openmc.Cell()
        c.region = +surf[-2] & -surf[-1]

        name = mat[i]
        if name == "mox70":
            name = "mox7"
        elif name == "moderator":
            name = "water"
        c.fill = material_dict[name]

        cell.append(c)

    surf[-1].boundary_type = "vacuum"

    geom = openmc.Geometry(cell)
    geom.export_to_xml()

    settings = openmc.Settings()
    settings.energy_mode = "multi-group"
    settings.particles = 10_000_000
    settings.batches = 50 + 1000
    settings.inactive = 50

    # Create an initial uniform spatial source distribution over fissionable zones
    bounds = [0.0, 0.0, 0.0, 21.42, 0.0, 0.0]
    uniform_dist = openmc.stats.Box(bounds[:3], bounds[3:])
    settings.source = openmc.IndependentSource(space=uniform_dist)

    settings.run_mode = "eigenvalue"

    settings.export_to_xml()

    openmc.run(threads=12)
