import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np

import openmc

# create a model to tie together geometry, materials, settings, and tallies
model = openmc.Model()

# For every cross section data set in the library, assign an openmc.Macroscopic object to a material
material_dict = {}
for xs in ["uo2", "mox43", "mox7", "mox87", "fiss_chamber", "guide_tube", "water"]:
    material_dict[xs] = openmc.Material(name=xs)
    material_dict[xs].set_density("macro", 1.0)
    material_dict[xs].add_macroscopic(xs)

# Instantiate a Materials collection, register all Materials, and export to XML
materials = openmc.Materials(material_dict.values())

# Set the location of the cross sections file to our pre-written set
# NOTE: these cross sections were produced by the OpenMC developers (not me).
# So, hopefully they serve as an independent check.
materials.cross_sections = "../c5g7.h5"
materials.export_to_xml()

left = openmc.XPlane(0.0)
left.boundary_type = "reflective"

interface = openmc.XPlane(0.54)

right = openmc.XPlane(0.63)
right.boundary_type = "reflective"

fuel_cell = openmc.Cell()
fuel_cell.region = +left & -interface
fuel_cell.fill = material_dict["mox7"]

mod_cell = openmc.Cell()
mod_cell.region = +interface & -right
mod_cell.fill = material_dict["water"]

geom = openmc.Geometry([fuel_cell, mod_cell])
geom.export_to_xml()

settings = openmc.Settings()
settings.energy_mode = "multi-group"
settings.particles = 1_000_000
settings.batches = 150
settings.inactive = 50

# Create an initial uniform spatial source distribution over fissionable zones
bounds = [0.0, 0.0, 0.0, 0.54, 0.0, 0.0]
uniform_dist = openmc.stats.Box(bounds[:3], bounds[3:])
settings.source = openmc.IndependentSource(space=uniform_dist)

settings.run_mode = "eigenvalue"

settings.export_to_xml()

openmc.run(threads=8)
