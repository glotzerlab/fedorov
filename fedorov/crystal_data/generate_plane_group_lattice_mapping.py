# Copyright (c) 2019-2020 The Regents of the University of Michigan
# This file is part of the fedorov project, released under the BSD 3-Clause
# License.

# Maintainer: Pengji Zhou

# NOTE: this is the code for record that generates the plane group number -
# lattice type mapping.  The use of this code is not required to use this
# package

# generate plane group and lattice type mapping
import json

plane_group_lattice_mapping = {}
for i in range(1, 18):
    if i <= 2:
        plane_group_lattice_mapping[i] = "oblique"
    elif i <= 9:
        plane_group_lattice_mapping[i] = "rectangular"
    elif i <= 12:
        plane_group_lattice_mapping[i] = "square"
    else:
        plane_group_lattice_mapping[i] = "hexagonal"

with open("plane_group_lattice_mapping.json", "w") as f:
    json.dump(plane_group_lattice_mapping, f)
