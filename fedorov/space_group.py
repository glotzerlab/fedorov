# Copyright (c) 2019-2020 The Regents of the University of Michigan
# This file is part of the fedorov project, released under the BSD 3-Clause
# License.

import copy
import json
import os
import pickle
import warnings

import numpy as np
import rowan
import spglib as spg

from . import data, lattice, util


class PlaneGroup:
    """A class for plane group symmetry operation.

    This class provides method to initialize a crystal unit cell with plane
    group and Wyckoff possition information.

    :param plane_group_number:
        Plane group number between 1 and 17.
    :type plane_group_number:
        int
    """

    plane_group_info_dir = os.path.join(
        data._DATA_PATH, "plane_group_info.pickle"
    )
    with open(plane_group_info_dir, "rb") as f:
        plane_group_info_dict = pickle.load(f)
    plane_group_lattice_mapping_dir = os.path.join(
        data._DATA_PATH, "plane_group_lattice_mapping.json"
    )
    with open(plane_group_lattice_mapping_dir, "r") as f:
        plane_group_lattice_mapping = json.load(
            f, object_hook=util.json_key_to_int
        )

    def __init__(self, plane_group_number=1):
        if plane_group_number <= 0 or plane_group_number > 17:
            raise ValueError(
                "plane_group_number must be an integer between 1 and 17"
            )

        self.plane_group_number = plane_group_number
        self.lattice_type = self.plane_group_lattice_mapping[
            self.plane_group_number
        ]
        self.lattice = lattice.lattice_system_dict_2D[self.lattice_type]
        info = self.plane_group_info_dict[plane_group_number]
        self.translations = info["translations"]
        self.rotations = info["rotations"]

    def print_info(self):
        print(
            f"Plane group number: {self.plane_group_number}\n"
            f"lattice type: {self.lattice_type}\n"
            f"Default parameters for lattice: {self.lattice.lattice_params}"
        )

    def get_basis_vectors(
        self,
        base_positions,
        base_type=[],
        base_quaternions=None,
        is_complete=False,
        apply_orientation=False,
    ):
        """Get the basis vectors for the defined crystall structure.

        :param base_positions:
            N by 2 np array of the Wyckoff postions
        :type base_positions:
            np.ndarray
        :param base_type:
            a list of string for particle type name
        :type base_type:
            list
        :param base_quaternions:
            N by 4 np array of quaternions, default None
        :type base_quaternions:
            np.ndarray
        :param is_complete:
            bool value to indicate if the positions are complete postions in a
            unitcell
        :type is_complete:
            bool
        :param apply_orientations:
            bool value to indicate if the space group symmetry should be applied
            to orientation
        :type apply_orientations:
            bool
        :return:
            basis_vectors
        :rtype:
            np.ndarray
        """

        # check input accuracy
        if (
            not isinstance(base_positions, np.ndarray)
            or len(base_positions.shape) == 1
            or base_positions.shape[1] != 2
        ):
            raise ValueError(
                "base_positions must be an numpy array of shape Nx3"
            )
        if apply_orientation:
            if (
                not isinstance(base_quaternions, np.ndarray)
                or len(base_quaternions.shape) == 1
                or base_quaternions.shape[1] != 4
            ):
                raise ValueError(
                    "base_quaternions must be an numpy array of shape Nx4"
                )
        if len(base_type):
            if (
                not isinstance(base_type, list)
                or len(base_type) != base_positions.shape[0]
                or not all(isinstance(i, str) for i in base_type)
            ):
                raise ValueError(
                    "base_type must contain a list of type name the same length"
                    "as the number of basis positions"
                )
        else:
            base_type = ["A"] * base_positions.shape[0]

        def _expand2d(matrix2d):
            matrix3d = np.identity(3)
            matrix3d[0:2, 0:2] = matrix2d
            return matrix3d

        threshold = 1e-6
        reflection_exist = False
        for i in range(0, len(self.rotations)):
            # Generate the new set of positions from the base
            pos = data.wrap(
                self.rotations[i].dot(base_positions.T).T + self.translations[i]
            )
            if apply_orientation:
                if np.linalg.det(self.rotations[i]) == 1:
                    quat_rotate = rowan.from_matrix(
                        _expand2d(self.rotations[i]), require_orthogonal=False
                    )
                    quat = rowan.multiply(quat_rotate, base_quaternions)
                else:
                    quat = base_quaternions
                    reflection_exist = True

            if i == 0:
                positions = pos
                type_list = copy.deepcopy(base_type)
                if apply_orientation:
                    quaternions = quat
            else:
                pos_comparisons = pos - positions[:, np.newaxis, :]
                norms = np.linalg.norm(pos_comparisons, axis=2)
                comps = np.all(norms > threshold, axis=0)
                positions = np.append(positions, pos[comps], axis=0)
                type_list += [x for x, y in zip(base_type, comps) if y]
                if apply_orientation:
                    quaternions = np.append(quaternions, quat[comps], axis=0)
                    if norms.min() < threshold:
                        warnings.warn(
                            "Orientation quaterions may have multiple values "
                            "for the same particle postion under the symmetry "
                            "operation for this space group and is not well "
                            "defined, only the first occurance is used."
                        )
        if reflection_exist:
            warnings.warn(
                "Reflection operation is included in this space "
                "group, and is ignored for quaternion calculation."
            )

        if is_complete and len(positions) != len(base_positions):
            raise ValueError(
                "the complete basis postions vector does not match the space "
                "group chosen. More positions are generated based on the "
                "symmetry operation within the provided space group"
            )

        if apply_orientation:
            return data.wrap(positions), type_list, quaternions
        else:
            return data.wrap(positions), type_list

    def get_lattice_vectors(self, **user_lattice_params):
        """Initialize the unitcell and return lattice vectors [a1, a2].

        :param user_lattice_params:
            unit cell parameters, provide a, b, theta where applicable
        :type user_lattice_params:
            float
        :return: lattice_vectors
        :rtype: np.ndarray
        """
        return self.lattice.get_lattice_vectors(**user_lattice_params)


class SpaceGroup:
    """A class for space group symmetry operation.

    This class provides method to initialize a crystal unit cell with space
    group and Wyckoff possition information.

    :param space_group_number:
        Space group number between 1 and 230.
    :type space_group_number:
        int
    """

    space_group_hall_mapping_dir = os.path.join(
        data._DATA_PATH, "space_group_hall_mapping.json"
    )
    with open(space_group_hall_mapping_dir, "r") as f:
        space_group_hall_mapping = json.load(
            f, object_hook=util.json_key_to_int
        )
    space_group_lattice_mapping_dir = os.path.join(
        data._DATA_PATH, "space_group_lattice_mapping.json"
    )
    with open(space_group_lattice_mapping_dir, "r") as f:
        space_group_lattice_mapping = json.load(
            f, object_hook=util.json_key_to_int
        )

    def __init__(self, space_group_number=1):
        if space_group_number <= 0 or space_group_number > 230:
            raise ValueError(
                "space_group_number must be an integer between 1 and 230"
            )

        self.space_group_number = space_group_number
        self.lattice_type = self.space_group_lattice_mapping[
            self.space_group_number
        ]
        self.lattice = lattice.lattice_system_dict_3D[self.lattice_type]
        info = spg.get_symmetry_from_database(
            self.space_group_hall_mapping[space_group_number]
        )
        self.translations = info["translations"]
        self.rotations = info["rotations"]

    def print_info(self):
        print(
            f"Space group number: {self.space_group_number}\n"
            f"lattice type: {self.lattice_type}\n"
            f"Default parameters for lattice: {self.lattice.lattice_params}"
        )

    def get_basis_vectors(
        self,
        base_positions,
        base_type=[],
        base_quaternions=None,
        is_complete=False,
        apply_orientation=False,
    ):
        """Get the basis vectors for the defined crystall structure.

        :param base_positions:
            N by 3 np array of the Wyckoff postions
        :type base_positions:
            np.ndarray
        :param base_type:
            a list of string for particle type name
        :type base_type:
            list
        :param base_quaternions:
            N by 4 np array of quaternions, default None
        :type base_quaternions:
            np.ndarray
        :param is_complete:
            bool value to indicate if the positions are complete postions in a
            unitcell
        :type is_complete:
            bool
        :param apply_orientations:
            bool value to indicate if the space group symmetry should be applied
            to orientatioin
        :type apply_orientations:
            bool
        :return:
            basis_vectors
        :rtype:
            np.ndarray
        """

        # check input accuracy
        if (
            not isinstance(base_positions, np.ndarray)
            or len(base_positions.shape) == 1
            or base_positions.shape[1] != 3
        ):
            raise ValueError(
                "base_positions must be an numpy array of shape Nx3"
            )
        if apply_orientation:
            if (
                not isinstance(base_quaternions, np.ndarray)
                or len(base_quaternions.shape) == 1
                or base_quaternions.shape[1] != 4
            ):
                raise ValueError(
                    "base_quaternions must be an numpy array of shape Nx4"
                )
        if len(base_type):
            if (
                not isinstance(base_type, list)
                or len(base_type) != base_positions.shape[0]
                or not all(isinstance(i, str) for i in base_type)
            ):
                raise ValueError(
                    "base_type must contain a list of type name the same length"
                    "as the number of basis positions"
                )
        else:
            base_type = ["A"] * base_positions.shape[0]

        threshold = 1e-6
        reflection_exist = False
        for i in range(0, len(self.rotations)):
            # Generate the new set of positions from the base
            pos = data.wrap(
                self.rotations[i].dot(base_positions.T).T + self.translations[i]
            )
            if apply_orientation:
                if np.linalg.det(self.rotations[i]) == 1:
                    quat_rotate = rowan.from_matrix(
                        self.rotations[i], require_orthogonal=False
                    )
                    quat = rowan.multiply(quat_rotate, base_quaternions)
                else:
                    quat = base_quaternions
                    reflection_exist = True

            if i == 0:
                positions = pos
                type_list = copy.deepcopy(base_type)
                if apply_orientation:
                    quaternions = quat
            else:
                pos_comparisons = pos - positions[:, np.newaxis, :]
                norms = np.linalg.norm(pos_comparisons, axis=2)
                comps = np.all(norms > threshold, axis=0)
                positions = np.append(positions, pos[comps], axis=0)
                type_list += [x for x, y in zip(base_type, comps) if y]
                if apply_orientation:
                    quaternions = np.append(quaternions, quat[comps], axis=0)
                    if norms.min() < threshold:
                        print(
                            "Orientation quaterions may have multiple values "
                            "for the same particle postion under the symmetry "
                            "operation for this space group and is not well "
                            "defined, only the first occurance is used."
                        )
        if reflection_exist:
            print(
                "Warning: reflection operation is included in this "
                "space group, and is ignored for quaternion calculation."
            )

        if is_complete and len(positions) != len(base_positions):
            raise ValueError(
                "the complete basis postions vector does not match the space "
                "group chosen. More positions are generated based on the "
                "symmetry operation within the provided space group"
            )

        if apply_orientation:
            return data.wrap(positions), type_list, quaternions
        else:
            return data.wrap(positions), type_list

    def get_lattice_vectors(self, **user_lattice_params):
        """Initialize the unitcell and return lattice vectors [a1, a2, a3].

        :param user_lattice_params:
            unit cell parameters, provide a, b, c, alpha, beta, gamma where
            applicable
        :type user_lattice_params:
            float
        :return: lattice_vectors
        :rtype: np.ndarray
        """
        return self.lattice.get_lattice_vectors(**user_lattice_params)


class PointGroup:
    """A class to access all point group symmetry operations.

    This class provides method to access all point group symmetry operation in
    both rotational matrix form or quaternion form.

    :param point_group_number:
        Point group number between 1 and 32.
    :type point_group_number:
        int
    """

    point_group_rotation_matrix_dir = os.path.join(
        data._DATA_PATH, "point_group_rotation_matrix_dict.pickle"
    )
    with open(point_group_rotation_matrix_dir, "rb") as f:
        point_group_rotation_matrix_dict = pickle.load(f)

    point_group_quat_dir = os.path.join(
        data._DATA_PATH, "point_group_quat_dict.json"
    )
    with open(point_group_quat_dir, "r") as f:
        point_group_quat_dict = json.load(f, object_hook=util.json_key_to_int)

    point_group_name_mapping_dir = os.path.join(
        data._DATA_PATH, "point_group_name_mapping.json"
    )
    with open(point_group_name_mapping_dir, "r") as f:
        point_group_name_mapping = json.load(
            f, object_hook=util.json_key_to_int
        )

    def __init__(self, point_group_number=1):
        if point_group_number <= 0 or point_group_number > 32:
            raise ValueError(
                "point_group_number must be an integer between 1 and 32"
            )

        self.point_group_number = point_group_number
        self.point_group_name = self.point_group_name_mapping[
            point_group_number
        ]
        self.rotation_matrix = self.point_group_rotation_matrix_dict[
            point_group_number
        ]["rotations"]
        self.quaternion = self.point_group_quat_dict[point_group_number]

    def print_info(self):
        print(
            f"Point group number: {self.point_group_number}\n"
            f"Name {self.point_group_name}"
        )

    def get_quaternion(self):
        """Get the quaternions for the point group symmetry.

        :return:
            list of quaternions
        :rtype:
            list
        """

        return self.quaternion

    def get_rotation_matrix(self):
        """Get the rotation matrixes for the point group symmetry.

        :return:
            n by 3 by 3 numpy array containing n rotational matrixes
        :rtype:
            numpy.ndarray
        """

        return self.rotation_matrix
