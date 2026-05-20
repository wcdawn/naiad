import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np

import openmc
import h5py

# create a model to tie together geometry, materials, settings, and tallies
model = openmc.Model()

fname_xs = "../c5g7.h5"

# For every cross section data set in the library, assign an openmc.Macroscopic object to a material
material_dict = {}
with h5py.File(fname_xs, "r") as h5:
    for xs in h5:
        xsname = xs
        if (xsname == "moderator"):
            xsname = "water"
        elif (xsname == "mox70"):
            xsname = "mox7"
        material_dict[xsname] = openmc.Material(name=xs)
        material_dict[xsname].set_density("macro", 1.0)
        material_dict[xsname].add_macroscopic(xs)

# Instantiate a Materials collection, register all Materials, and export to XML
materials = openmc.Materials(material_dict.values())

# Set the location of the cross sections file to our pre-written set
materials.cross_sections = fname_xs
materials.export_to_xml()

left = openmc.XPlane(0.0)
left.boundary_type = "reflective"

interface = openmc.XPlane(0.54)

right = openmc.XPlane(0.63)
right.boundary_type = "reflective"

fuel_cell = openmc.Cell()
fuel_cell.region = +left & -interface
fuel_cell.fill = material_dict["uo2"]

mod_cell = openmc.Cell()
mod_cell.region = +interface & -right
mod_cell.fill = material_dict["water"]

geom = openmc.Geometry([fuel_cell, mod_cell])
geom.export_to_xml()

settings = openmc.Settings()
settings.energy_mode = "multi-group"
settings.particles = 4_000_000
settings.batches = 1_000 + 50
settings.inactive = 50

# Create an initial uniform spatial source distribution over fissionable zones
bounds = [0.0, 0.0, 0.0, 0.54, 0.0, 0.0]
uniform_dist = openmc.stats.Box(bounds[:3], bounds[3:])
settings.source = openmc.IndependentSource(space=uniform_dist)

settings.run_mode = "eigenvalue"

settings.export_to_xml()

openmc.run(threads=12)
