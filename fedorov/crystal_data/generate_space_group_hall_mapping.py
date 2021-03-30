# Copyright (c) 2019-2020 The Regents of the University of Michigan
# This file is part of the fedorov project, released under the BSD 3-Clause
# License.

# Maintainer: Pengji Zhou

# NOTE: this is the code for record that generates the hall number - space group
# number mapping.  The use of this code is not required to use this package

import json

# hall_spacegroup_mapping
import spglib as spg

two_origin_choice_list = [
    48,
    50,
    59,
    68,
    70,
    85,
    86,
    88,
    125,
    126,
    129,
    130,
    133,
    134,
    137,
    138,
    141,
    142,
    201,
    203,
    222,
    224,
    227,
    228,
]
H_R_choice_list = [146, 148, 155, 160, 161, 166, 167]

spacegroup_hall_mapping = {}
for hall_number in range(1, 531):
    spacegroup_type = spg.get_spacegroup_type(hall_number)
    if spacegroup_type["number"] not in spacegroup_hall_mapping:
        if (
            spacegroup_type["number"]
            in two_origin_choice_list + H_R_choice_list
        ):
            shift = 1
        else:
            shift = 0
        spacegroup_hall_mapping[spacegroup_type["number"]] = hall_number + shift

with open("space_group_hall_mapping.json", "w") as f:
    json.dump(spacegroup_hall_mapping, f)
