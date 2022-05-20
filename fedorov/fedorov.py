# Copyright (c) 2019-2020 The Regents of the University of Michigan
# This file is part of the fedorov project, released under the BSD 3-Clause
# License.

# Maintainer: Pengji Zhou

import copy
import json
import os
import re

import numpy as np
import pandas as pd

from . import data, space_group

_WYCKOFF_FILE = "space_group_{}_Wyckoff_site_data.json"


class Prototype:
    """Crystal prototype class.

    This class uses the minimal necessay information needed to fully define a
    crystal structures with space group number, wyckoff postions(in letter name
    convention) and free parameters for each relavent wyckoff postions.

    :param space_group_number:
        space group number between 1 and 230
    :type space_group_numbers:
        int
    :param wyckoff_site:
        wyckoff site letters included in the prototype
    :type wyckoff_site:
        str
    :param type_by_site:
        type name letter for each site set in wyckoff_sites
    :type type_by_site:
        str
    """

    def __init__(
        self,
        space_group_number=1,
        wyckoff_site="",
        type_by_site="",
    ):
        if space_group_number > 230 or space_group_number < 1:
            raise ValueError(
                "space_group_number must be an integer between 0 and 230, "
                "default = 1"
            )

        if not isinstance(wyckoff_site, str) or not wyckoff_site.isalpha():
            raise ValueError(
                "wyckoff_postions must be string consists of all the Wyckoff "
                "postions letters, e.g. 'abcc' denotes one set of Wyckoff "
                "postions for both a and b, and two sets at Wyckoff postion c"
            )

        if type_by_site == "":
            type_by_site = "A" * len(wyckoff_site)
        elif (
            not isinstance(type_by_site, str)
            or len(type_by_site) != len(wyckoff_site)
            or not type_by_site.isalpha()
        ):
            raise ValueError(
                "type_by_site must be string consists of type name (A/B/C, etc)"
                "for each Wyckoff site, default all particles with same type A "
                "if not provided"
            )

        wyckoff_site_list = list(wyckoff_site.lower())
        type_by_site = list(type_by_site.upper())

        wyckoff_data_dir = os.path.join(
            data._DATA_PATH, _WYCKOFF_FILE.format(space_group_number)
        )
        with open(wyckoff_data_dir, "r") as f:
            full_wyckoff_positions = json.load(f)

        basis_params_list = []
        order = 1
        for site in wyckoff_site_list:
            pos = copy.deepcopy(full_wyckoff_positions[site])
            pos = "".join(pos)
            for letter in ("x", "y", "z"):
                if letter in pos:
                    basis_params_list.append(letter + str(order))
            order += 1
        basis_params_value_list = [None] * len(basis_params_list)

        self.space_group_number = space_group_number
        self.space_group = space_group.SpaceGroup(space_group_number)
        self.wyckoff_site_list = wyckoff_site_list
        self.full_wyckoff_positions = full_wyckoff_positions
        self.type_by_site = type_by_site
        self.lattice_params = self.space_group.lattice.lattice_params
        self.basis_params = dict(
            zip(basis_params_list, basis_params_value_list)
        )

    def print_info(self):
        print(
            f"Wyckoff sites:{self.wyckoff_site_list}\n",
            f"Particle type for each Wyckoff sites:{self.type_by_site}\n"
            f"lattice parameters list:{list(self.lattice_params)}\n"
            f"basis parameters list:{self.basis_params}",
        )

    def update_basis_params(self, user_basis_params):
        params = copy.deepcopy(self.basis_params)
        for param, value in user_basis_params.items():
            if param in params:
                if value is not None:
                    params[param] = value
            else:
                print(
                    "warning: {} is not required and not used to define this "
                    "structure".format(param)
                )
        for value in params.values():
            if value is None:
                raise ValueError(
                    "not all necessary parameters were provided! Call "
                    "print_info() to see the full list of necessary parameters"
                )
        return params

    def get_basis_vectors(self, **user_basis_params):
        """Initialize fractional coordinates of the particles in the unitcell.

        :param user_basis_params:
            user defined parameters for different Wyckoff site degree of
            freedom, when applicable
        :type user_basis_params:
            float
        :return:
            basis_vectors
        :rtype:
            np.ndarray
        """
        basis_params = self.update_basis_params(user_basis_params)
        base_positions = np.zeros((0, 3))

        order = 1
        for site in self.wyckoff_site_list:
            pos = copy.deepcopy(self.full_wyckoff_positions[site])
            for letter in ("x", "y", "z"):
                if letter + str(order) in basis_params.keys():
                    exec(
                        "{} = {}".format(
                            letter, basis_params[letter + str(order)]
                        )
                    )
            for i in range(0, 3):
                # add * back for eval
                target = re.findall(r"(\d[xyz])", pos[i])
                for item in target:
                    pos[i] = pos[i].replace(
                        item, "*".join(re.findall(r"(\d)([xyz])", item)[0])
                    )
                pos[i] = eval(pos[i])
            base_positions = np.append(
                base_positions, np.array(pos).reshape(1, -1), axis=0
            )
            order += 1

        return self.space_group.get_basis_vectors(
            data.wrap(base_positions), base_type=self.type_by_site
        )

    def update_lattice_params(self, user_lattice_params):
        params = copy.deepcopy(self.lattice_params)
        for param, value in user_lattice_params.items():
            if param in params:
                if value is not None:
                    params[param] = value
            else:
                print(
                    "warning: {} is not required and not used to define this "
                    "structure".format(param)
                )
        return params

    def get_lattice_vectors(self, **user_lattice_params):
        """Initialize the unitcell and return lattice vectors [a1, a2, a3]

        :param user_lattice_params:
            unit cell parameters, provide a, b, c, alpha, beta, gamma where
            applicable
        :type user_lattice_params:
            float
        :return:
            lattice_vectors
        :rtype:
            np.ndarray
        """
        lattice_params = self.update_lattice_params(user_lattice_params)
        return self.space_group.lattice.get_lattice_vectors(**lattice_params)


class AflowPrototype(Prototype):
    """Aflow prototype class.

    This class uses the crystal prototypes in Aflow database to initialize
    crystal structures.

    :param prototype_index:
        prototype index [0, 589] for all 590 prototypes in AFLOW.
    :type prototype_index:
        int
    :param set_type:
        allow setting different type name(in A, B, C order) for different atoms
        in AFLOW prototype
    :type set_type:
        bool
    """

    _Aflow_database = pd.read_csv(
        os.path.join(data._DATA_PATH, "Aflow_processed_data.csv"), index_col=0
    )
    _name_regex = re.compile(r"'(.*?)'")

    def __init__(self, prototype_index=0, set_type=False):
        if prototype_index < 0 or prototype_index >= 590:
            raise ValueError(
                "prototype_index must be an integer between 0 and 590."
            )
        entry = self._Aflow_database.iloc[prototype_index]

        def get_values(value_str: str):
            try:
                return [float(i) for i in value_str.strip("[]").split(",")]
            except Exception:
                return []

        lattice_params = self._name_regex.findall(entry["lattice_params"])
        lattice_params_value = get_values(entry["lattice_params_value"])
        lattice_params = dict(zip(lattice_params, lattice_params_value))

        basis_params = self._name_regex.findall(entry["basis_params"])
        basis_params_value = get_values(entry["basis_params_value"])
        basis_params = dict(zip(basis_params, basis_params_value))

        # convert Aflow angle unit from degree to rad
        for key in {"alpha", "beta", "gamma"} & lattice_params.keys():
            lattice_params[key] = lattice_params[key] / 180 * np.pi

        space_group_number = entry["space_group"]
        # process proper unitcell params
        if space_group_number in {146, 148, 155, 160, 161, 166, 167}:
            a = lattice_params.pop("a")
            c = lattice_params.pop("c/a") * a
            lattice_params["a"] = np.sqrt(a**2 / 3 + c**2 / 9)
            lattice_params["alpha"] = np.arccos(
                (2 * c**2 - 3 * a**2) / (2 * (c**2 + 3 * a**2))
            )
        else:
            a = lattice_params["a"]
            if "b/a" in lattice_params:
                lattice_params["b"] = lattice_params.pop("b/a") * a
            if "c/a" in lattice_params:
                lattice_params["c"] = lattice_params.pop("c/a") * a

        wyckoff_sites, wyckoff_positions, types = self._get_wyckoff_sites(
            entry, space_group_number, set_type
        )

        self.id = entry["id"]
        self.pearson_symbol = entry["pearson_symbol"]
        self.prototype = entry["prototype"]
        self.space_group_number = space_group_number
        self.space_group = space_group.SpaceGroup(space_group_number)
        self.wyckoff_site_list = wyckoff_sites
        self.full_wyckoff_positions = wyckoff_positions
        self.type_by_site = types
        self.lattice_params = lattice_params
        self.basis_params = basis_params

    def print_info(self):
        print(
            f"Info for the chosen crystal structure prototype:\n"
            f"id: {self.id}, (Pearson-Chemistry-SpaceGroup)\n"
            f"Wyckoff sites: {self.wyckoff_site_list}\n"
            f"available lattice parameters: {self.lattice_params}\n"
            f"available basis parameters: {self.basis_params}"
        )

    def _get_wyckoff_sites(self, entry, space_group_number, set_type):
        wyckoff_sites_by_type = self._name_regex.findall(entry["wyckoff_sites"])
        wyckoff_sites = sorted("".join(wyckoff_sites_by_type))

        wyckoff_data_path = os.path.join(
            data._DATA_PATH, _WYCKOFF_FILE.format(space_group_number)
        )

        with open(wyckoff_data_path, "r") as f:
            wyckoff_positions = json.load(f)

        # get type label
        if not set_type:
            return wyckoff_sites, wyckoff_positions, ["A"] * len(wyckoff_sites)

        type_by_site = list("A" * len(wyckoff_sites))
        sorted_site_string = "".join(wyckoff_sites)
        base = ord("A")
        for wyckoffs in wyckoff_sites_by_type:
            for site in wyckoffs:
                order = sorted_site_string.find(site)
                sorted_site_string = sorted_site_string.replace(site, "0", 1)
                type_by_site[order] = chr(base)
            base += 1
        return wyckoff_sites, wyckoff_positions, type_by_site

    @classmethod
    def from_query(
        cls,
        pearson_symbol: "str | None" = None,
        space_group: "int | None" = None,
        prototype: "str | None" = None,
        set_type: bool = False,
    ):
        """Create all `AflowPrototype` matching the given query.

        Args:
            pearson_symbol (`str`, optional): The Pearson symbol to search for,
                defaults to ``None`` which accepts any Pearson symbol.
            space_group (`int`, optional): The space group to search for,
                defaults to ``None`` which accepts any space group.
            prototype (`str`, optional): The chemical prototype to search for,
                defaults to ``None`` which accepts any prototype.
            set_type (`bool`, optional): Set different type name (in alphabetic
                order starting with "A") for different atoms in AFLOW prototype.

        Returns:
            lattices (list[`AflowPrototype`]): The list of all
                `AflowPrototype`'s with a given Pearson symbol.
        """

        def query(row):
            return (
                (prototype is None or row.prototype == prototype)
                and (
                    pearson_symbol is None
                    or row.pearson_symbol == pearson_symbol
                )
                and (space_group is None or row.space_group == space_group)
            )

        return [
            cls(i, set_type)
            for i, row in cls._Aflow_database.iterrows()
            if query(row)
        ]
