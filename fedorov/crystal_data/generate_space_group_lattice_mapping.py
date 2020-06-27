# Copyright (c) 2019-2020 The Regents of the University of Michigan
# This file is part of the fedorov project, released under the BSD 3-Clause License.

# Maintainer: Pengji Zhou

# NOTE: this is the code for record that generates the space group number - lattice type mapping.
# The use of this code is not required to use this package

# generate spacegroup and lattice type mapping
import json

space_group_lattice_mapping = {}
for i in range(1, 231):
    if i in [146, 148, 155, 160, 161, 166, 167]:
        space_group_lattice_mapping[i] = 'rhombohedral'
    elif i <= 2:
        space_group_lattice_mapping[i] = 'triclinic'
    elif i <= 15:
        space_group_lattice_mapping[i] = 'monoclinic'
    elif i <= 74:
        space_group_lattice_mapping[i] = 'orthorhombic'
    elif i <= 142:
        space_group_lattice_mapping[i] = 'tetragonal'
    elif i <= 194:
        space_group_lattice_mapping[i] = 'hexagonal'
    elif i <= 230:
        space_group_lattice_mapping[i] = 'cubic'

with open('space_group_lattice_mapping.json', 'w') as f:
    json.dump(space_group_lattice_mapping, f)
