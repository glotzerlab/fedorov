# Copyright (c) 2019-2020 The Regents of the University of Michigan
# This file is part of the fedorov project, released under the BSD 3-Clause License.

# Maintainer: Pengji Zhou

import numpy as np
import pandas as pd
import json
import pickle
import spglib as spg
import copy
import re
import rowan
import os


def json_key_to_int(x):
    if isinstance(x, dict):
        return {int(k): v for k, v in x.items()}
    return x


def wrap(basis_vectors):
    """Wrap fractional coordinates within a unitcell based on periodic boundary.

    :param basis_vectors:
        fractional coordinates for particle positions in N by 3 numpy array
    :type basis_vectors:
        np.ndarray
    :return:
        wrapped basis_vectors
    :rtype:
        np.ndarray
    """
    basis_vectors[basis_vectors >= 1.0] -= 1
    basis_vectors[basis_vectors < 0.0] += 1
    return basis_vectors


def convert_to_box(lattice_vectors):
    """Convert lattice vectors to box parameter: Lx, Ly, Lz, xy, xz, yz.

    :param lattice_vectors:
        3 by 3 numpy array of lattice vectors [a1, a2, a3]
    :type lattice_vectors:
        np.ndarray
    :return:
        Lx, Ly, Lz, xy, xz, yz
    :rtype:
        float
    """
    v0 = lattice_vectors[0, :]
    v1 = lattice_vectors[1, :]
    v2 = lattice_vectors[2, :]
    Lx = np.sqrt(np.dot(v0, v0))
    a2x = np.dot(v0, v1) / Lx
    Ly = np.sqrt(np.dot(v1, v1) - a2x * a2x)
    xy = a2x / Ly
    v0xv1 = np.cross(v0, v1)
    v0xv1mag = np.sqrt(np.dot(v0xv1, v0xv1))
    Lz = np.dot(v2, v0xv1) / v0xv1mag
    a3x = np.dot(v0, v2) / Lx
    xz = a3x / Lz
    yz = (np.dot(v1, v2) - a2x * a3x) / (Ly * Lz)
    return Lx, Ly, Lz, xy, xz, yz


def convert_to_vectors(Lx, Ly, Lz, xy, xz, yz):
    """Convert box parameter: Lx, Ly, Lz, xy, xz, yz to lattice vectors [a1, a2, a3].

    :param Lx:
    :type Lx:
        float
    :param Ly:
    :type Ly:
        float
    :param Lz:
    :type Lz:
        float
    :param xy:
    :type xy:
        float
    :param xz:
    :type xz:
        float
    :param yz:
    :type yz:
        float
    :return:
        lattice_vectors in form [a1, a2, a3]
    :rtype:
        np.ndarray
    """
    lattice_vectors = np.array([[Lx, 0, 0],
                                [xy * Ly, Ly, 0],
                                [xz * Lz, yz * Lz, Lz]])
    return lattice_vectors


def get_volumn(lattice_vectors):
    """Calculate volume of the unitcell

    :param lattice_vectors:
        3 by 3 numpy array of lattice vectors [a1, a2, a3]
    :type lattice_vectors:
        np.ndarray
    :return:
        volume
    :rtype:
        float
    """
    a1 = lattice_vectors[0, :]
    a2 = lattice_vectors[1, :]
    a3 = lattice_vectors[2, :]
    return abs(np.cross(a1, a2).dot(a3))


def fractional_to_cartesian(basis_vectors, lattice_vectors):
    """Convert fractional coordinates to cartesian coordinates.

    :param basis_vectors:
        N by 3 numpy array of basis vectors
    :type basis_vectors:
        np.ndarray
    :param lattice_vectors:
        3 by 3 numpy array of lattice vectors [a1, a2, a3]
    :type lattice_vectors:
        np.ndarray
    :return:
        N by 3 numpy array of cartesiann coordinates
    :rtype:
        np.ndarray
    """
    return basis_vectors.dot(lattice_vectors)


def translate_to_vector(a=1, b=1, c=1, alpha=np.pi / 2, beta=np.pi / 2, gamma=np.pi / 2):
    """Convert box parameters a, b, c, alpha, beta, gamma to lattice vectors [a1, a2, a3].

    :param a:
    :type a:
        float
    :param b:
    :type b:
        float
    :param c:
    :type c:
        float
    :param alpha:
    :type alpha:
        float
    :param beta:
    :type beta:
        float
    :param gamma:
    :type gamma:
        float
    :return:
        lattice_vectors
    :rtype:
        np.ndarray
    """
    # https://en.wikipedia.org/wiki/Fractional_coordinates
    cg = np.cos(gamma)
    sg = np.sin(gamma)
    ca = np.cos(alpha)
    cb = np.cos(beta)
    cy = (ca - cb * cg) / sg
    if (1 - ca * ca - cb * cb - cg * cg + 2 * ca * cb * cg) < 0:
        raise ValueError('Error: the box length and angle parameters provided are not feasible. '
                         'Please not the unit used for angle paramters should be in unit of rad')
    cz = np.sqrt(1 - ca * ca - cb * cb - cg * cg + 2 * ca * cb * cg) / sg
    lattice_vectors = np.array([
                               [a, 0, 0],
                               [b * cg, b * sg, 0],
                               [c * cb, c * cy, c * cz]
                               ])
    return lattice_vectors


def translate_to_vector_2D(a=1, b=1, theta=np.pi / 2):
    """Convert box parameters a, b, theta to lattice vectors [a1, a2].

    :param a:
    :type a:
        float
    :param b:
    :type b:
        float
    :param theta:
    :type theta:
        float
    :return:
        lattice_vectors
    :rtype:
        np.ndarray
    """
    lattice_vectors = np.array([
                               [a, 0],
                               [b * np.cos(theta), b * np.sin(theta)]
                               ])
    return lattice_vectors


# 2D systems
class Oblique2D(object):
    """A class for constructing a 2D oblique unitcell

    This class provides method to initialize a 2D oblique unitcell
    """
    lattice_params = {'a': 1, 'b': 1, 'theta': np.pi / 2}

    @classmethod
    def update_lattice_params(cls, user_lattice_params):
        params = copy.deepcopy(cls.lattice_params)
        for param, value in user_lattice_params.items():
            if param in params:
                if value is not None:
                    params[param] = value
            else:
                print('warning: {} is not required and not used to define this '
                      'structure'.format(param))
        return params

    @classmethod
    def get_lattice_vectors(cls, **user_lattice_params):
        """Initialize a 2D oblique unitcell and return lattice vectors [a1, a2].

        :param user_lattice_params:
            unit cell parameters, provide a, b, theta where applicable
        :type user_lattice_params:
            float
        :return:
            lattice_vectors
        :rtype:
            np.ndarray
        """
        params = cls.update_lattice_params(user_lattice_params)
        lattice_vectors = translate_to_vector_2D(**params)
        return lattice_vectors


class Rectangular2D(Oblique2D):
    """A class for constructing a 2D rectangular unitcell

    This class provides method to initialize a 2D rectangular unitcell
    """
    lattice_params = {'a': 1, 'b': 1}

    @classmethod
    def get_lattice_vectors(cls, **user_lattice_params):
        """Initialize a 2D rectangular unitcell and return lattice vectors [a1, a2].

        :param user_lattice_params:
            unit cell parameters, provide a, b, theta where applicable
        :type user_lattice_params:
            float
        :return:
            lattice_vectors
        :rtype:
            np.ndarray
        """
        params = cls.update_lattice_params(user_lattice_params)
        lattice_vectors = np.array([
                                   [params['a'], 0.0],
                                   [0.0, params['b']]
                                   ])
        return lattice_vectors


class Hexagonal2D(Oblique2D):
    """A class for constructing a 2D hexagonal unitcell

    This class provides method to initialize a 2D hexagonal unitcell
    """
    lattice_params = {'a': 1}

    @classmethod
    def get_lattice_vectors(cls, **user_lattice_params):
        """Initialize a 2D hexagonal unitcell and return lattice vectors [a1, a2].

        :param user_lattice_params:
            unit cell parameters, provide a, b, theta where applicable
        :type user_lattice_params:
            float
        :return:
            lattice_vectors
        :rtype:
            np.ndarray
        """
        params = cls.update_lattice_params(user_lattice_params)
        lattice_vectors = np.array([
                                   [params['a'], 0.0],
                                   [-0.5 * params['a'], params['a'] * np.sqrt(3) / 2]
                                   ])
        return lattice_vectors


class Square2D(Rectangular2D):
    """A class for constructing a 2D square unitcell

    This class provides method to initialize a 2D square unitcell
    """
    lattice_params = {'a': 1}

    @classmethod
    def get_lattice_vectors(cls, **user_lattice_params):
        """Initialize a 2D square unitcell and return lattice vectors [a1, a2].

        :param user_lattice_params:
            unit cell parameters, provide a, b, theta where applicable
        :type user_lattice_params:
            float
        :return:
            lattice_vectors
        :rtype:
            np.ndarray
        """
        params = cls.update_lattice_params(user_lattice_params)
        lattice_vectors = np.array([
                                   [params['a'], 0.0],
                                   [0.0, params['a']]
                                   ])
        return lattice_vectors


lattice_system_dict_2D = {'oblique': Oblique2D, 'rectangular': Rectangular2D,
                          'hexagonal': Hexagonal2D, 'square': Square2D}


class PlaneGroup(object):
    """A class for plane group symmetry operation.

    This class provides method to initialize a crystal unit cell with plane group and Wyckoff
    possition information.

    :param plane_group_number:
        Plane group number between 1 and 17.
    :type plane_group_number:
        int
    :param print_info:
        Print plane group information upon initialization.
    :type print_info:
        bool
    """

    dir_path = os.path.dirname(__file__)
    plane_group_info_dir = os.path.join(dir_path,
                                        'crystal_data/plane_group_info.pickle')
    with open(plane_group_info_dir, 'rb') as f:
        plane_group_info_dict = pickle.load(f)
    plane_group_lattice_mapping_dir = os.path.join(dir_path,
                                                   'crystal_data/plane_group_lattice_mapping.json')
    with open(plane_group_lattice_mapping_dir, 'r') as f:
        plane_group_lattice_mapping = json.load(f, object_hook=json_key_to_int)

    def __init__(self, plane_group_number=1, print_info=False):
        if plane_group_number <= 0 or plane_group_number > 17:
            raise ValueError('plane_group_number must be an integer between 1 and 17')

        self.plane_group_number = plane_group_number
        self.lattice_type = self.plane_group_lattice_mapping[self.plane_group_number]
        self.lattice = lattice_system_dict_2D[self.lattice_type]
        info = self.plane_group_info_dict[plane_group_number]
        self.translations = info['translations']
        self.rotations = info['rotations']

        if print_info:
            print('Plane group number: {}\n'.format(plane_group_number),
                  'lattice type: {}\n'.format(self.lattice_type),
                  'Default parameters for lattice: {}'.format(self.lattice.lattice_params))

    def get_basis_vectors(self, base_positions, base_type=[], base_quaternions=None,
                          is_complete=False, apply_orientation=False):
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
            bool value to indicate if the positions are complete postions in a unitcell
        :type is_complete:
            bool
        :param apply_orientations:
            bool value to indicate if the space group symmetry should be applied to orientatioin
        :type apply_orientations:
            bool
        :return:
            basis_vectors
        :rtype:
            np.ndarray
        """

        # check input accuracy
        if not isinstance(base_positions, np.ndarray) or len(base_positions.shape) == 1 \
           or base_positions.shape[1] != 2:
            raise ValueError('base_positions must be an numpy array of shape Nx3')
        if apply_orientation:
            if not isinstance(base_quaternions, np.ndarray) or len(base_quaternions.shape) == 1 \
               or base_quaternions.shape[1] != 4:
                raise ValueError('base_quaternions must be an numpy array of shape Nx4')
        if len(base_type):
            if not isinstance(base_type, list) or len(base_type) != base_positions.shape[0] \
               or not all(isinstance(i, str) for i in base_type):
                raise ValueError('base_type must contain a list of type name the same length'
                                 'as the number of basis positions')
        else:
            base_type = ['A'] * base_positions.shape[0]

        def _expand2d(matrix2d):
            matrix3d = np.identity(3)
            matrix3d[0:2, 0:2] = matrix2d
            return matrix3d

        threshold = 1e-6
        reflection_exist = False
        for i in range(0, len(self.rotations)):
            # Generate the new set of positions from the base
            pos = wrap(self.rotations[i].dot(base_positions.T).T + self.translations[i])
            if apply_orientation:
                if np.linalg.det(self.rotations[i]) == 1:
                    quat_rotate = rowan.from_matrix(_expand2d(self.rotations[i]),
                                                    require_orthogonal=False)
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
                        print('Orientation quaterions may have multiple values for the same '
                              'particle postion under the symmetry operation for this space group '
                              'and is not well defined, only the first occurance is used.')
        if reflection_exist:
            print('Warning: reflection operation is included in this space group, '
                  'and is ignored for quaternion calculation.')

        if is_complete and len(positions) != len(base_positions):
            raise ValueError('the complete basis postions vector does not match the space group '
                             'chosen. More positions are generated based on the symmetry operation '
                             'within the provided space group')

        if apply_orientation:
            return wrap(positions), type_list, quaternions
        else:
            return wrap(positions), type_list

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


# 3D systems
class Triclinic(object):
    """A class for constructing a triclinic unitcell."""
    Pearson = 'a'
    dimensions = 3
    lattice_params = {'a': 1, 'b': 1, 'c': 1, 'alpha': np.pi / 2, 'beta': np.pi / 2,
                      'gamma': np.pi / 2}

    @classmethod
    def update_lattice_params(cls, user_lattice_params):
        params = copy.deepcopy(cls.lattice_params)
        for param, value in user_lattice_params.items():
            if param in params:
                if value is not None:
                    params[param] = value
            else:
                print('warning: {} is not required and not used to define this '
                      'structure'.format(param))
        return params

    @classmethod
    def get_lattice_vectors(cls, **user_lattice_params):
        """Initialize a triclinic unitcell and return lattice vectors [a1, a2, a3].

        :param user_lattice_params:
            unit cell parameters, provide a, b, c, alpha, beta, gamma where applicable
        :type user_lattice_params:
            float
        :return:
            lattice_vectors
        :rtype:
            np.ndarray
        """
        params = cls.update_lattice_params(user_lattice_params)
        lattice_vectors = translate_to_vector(**params)
        return lattice_vectors


class Monoclinic(Triclinic):
    """A class for constructing a monoclinic unitcell

    This class provides method to initialize a monoclinic unitcell
    """
    Pearson = 'm'
    lattice_params = {'a': 1, 'b': 1, 'c': 1, 'beta': np.pi / 2}

    @classmethod
    def get_lattice_vectors(cls, **user_lattice_params):
        """Initialize a monoclinic unitcell and return lattice vectors [a1, a2, a3].

        :param user_lattice_params:
            unit cell parameters, provide a, b, c, alpha, beta, gamma where applicable
        :type user_lattice_params:
            float
        :return:
            lattice_vectors
        :rtype:
            np.ndarray
        """
        params = cls.update_lattice_params(user_lattice_params)
        lattice_vectors = translate_to_vector(a=params['a'], b=params['b'],
                                              c=params['c'], beta=params['beta'])
        return lattice_vectors


class Orthorhombic(Monoclinic):
    """A class for constructing a orthorhombic unitcell."""
    Pearson = 'o'
    lattice_params = {'a': 1, 'b': 1, 'c': 1}

    @classmethod
    def get_lattice_vectors(cls, **user_lattice_params):
        """Initialize a orthorhombi unitcell and return lattice vectors [a1, a2, a3].

        :param user_lattice_params:
            unit cell parameters, provide a, b, c, alpha, beta, gamma where applicable
        :type user_lattice_params:
            float
        :return:
            lattice_vectors
        :rtype:
            np.ndarray
        """
        params = cls.update_lattice_params(user_lattice_params)
        lattice_vectors = np.array([
                                   [params['a'], 0.0, 0.0],
                                   [0.0, params['b'], 0.0],
                                   [0.0, 0.0, params['c']]
                                   ])
        return lattice_vectors


class Tetragonal(Orthorhombic):
    """A class for constructing a tetragonal unitcell."""
    Pearson = 't'
    lattice_params = {'a': 1, 'c': 1}

    @classmethod
    def get_lattice_vectors(cls, **user_lattice_params):
        """Initialize a tetragona unitcell and return lattice vectors [a1, a2, a3].

        :param user_lattice_params:
            unit cell parameters, provide a, b, c, alpha, beta, gamma where applicable
        :type user_lattice_params:
            float
        :return:
            lattice_vectors
        :rtype:
            np.ndarray
        """
        params = cls.update_lattice_params(user_lattice_params)
        lattice_vectors = np.array([
                                   [params['a'], 0.0, 0.0],
                                   [0.0, params['a'], 0.0],
                                   [0.0, 0.0, params['c']]
                                   ])
        return lattice_vectors


class Hexagonal(Triclinic):
    """A class for constructing a hexagonal unitcell."""
    Pearson = 'hP'
    lattice_params = {'a': 1, 'c': 1}

    @classmethod
    def get_lattice_vectors(cls, **user_lattice_params):
        """Initialize a hexagonal unitcell and return lattice vectors [a1, a2, a3].

        :param user_lattice_params:
            unit cell parameters, provide a, b, c, alpha, beta, gamma where applicable
        :type user_lattice_params:
            float
        :return:
            lattice_vectors
        :rtype:
            np.ndarray
        """
        params = cls.update_lattice_params(user_lattice_params)
        lattice_vectors = np.array([
                                   [params['a'], 0.0, 0.0],
                                   [-0.5 * params['a'], np.sqrt(3.0) / 2.0 * params['a'], 0.0],
                                   [0.0, 0.0, params['c']]
                                   ])
        return lattice_vectors


class Rhombohedral(Triclinic):
    """A class for constructing a rhombohedral unitcell."""
    Pearson = 'hR'
    lattice_params = {'a': 1, 'alpha': np.pi / 2}

    @classmethod
    def get_lattice_vectors(cls, **user_lattice_params):
        """Initialize a rhombohedral unitcell and return lattice vectors [a1, a2, a3].

        :param user_lattice_params:
            unit cell parameters, provide a, b, c, alpha, beta, gamma where applicable
        :type user_lattice_params:
            float
        :return:
            lattice_vectors
        :rtype:
            np.ndarray
        """
        params = cls.update_lattice_params(user_lattice_params)
        lattice_vectors = translate_to_vector(a=params['a'], b=params['a'], c=params['a'],
                                              alpha=params['alpha'], beta=params['alpha'],
                                              gamma=params['alpha'])
        return lattice_vectors


class Cubic(Tetragonal):
    """A class for constructing a cubic unitcell."""
    Pearson = 'c'
    lattice_params = {'a': 1}

    @classmethod
    def get_lattice_vectors(cls, **user_lattice_params):
        """Initialize a cubicc unitcell and return lattice vectors [a1, a2, a3].

        :param user_lattice_params:
            unit cell parameters, provide a, b, c, alpha, beta, gamma where applicable
        :type user_lattice_params:
            float
        :return:
            lattice_vectors
        :rtype:
            np.ndarray
        """
        params = cls.update_lattice_params(user_lattice_params)
        lattice_vectors = np.array([
                                   [params['a'], 0.0, 0.0],
                                   [0.0, params['a'], 0.0],
                                   [0.0, 0.0, params['a']]
                                   ])
        return lattice_vectors


lattice_system_dict_3D = {'triclinic': Triclinic, 'monoclinic': Monoclinic,
                          'orthorhombic': Orthorhombic, 'tetragonal': Tetragonal,
                          'hexagonal': Hexagonal, 'rhombohedral': Rhombohedral, 'cubic': Cubic}


class SpaceGroup(object):
    """A class for space group symmetry operation.

    This class provides method to initialize a crystal unit cell with space group and Wyckoff
    possition information.

    :param space_group_number:
        Space group number between 1 and 230.
    :type space_group_number:
        int
    :param print_info:
        Print space group information upon initialization.
    :type print_info:
        bool
    """
    dir_path = os.path.dirname(__file__)
    space_group_hall_mapping_dir = os.path.join(dir_path,
                                                'crystal_data/space_group_hall_mapping.json')
    with open(space_group_hall_mapping_dir, 'r') as f:
        space_group_hall_mapping = json.load(f, object_hook=json_key_to_int)
    space_group_lattice_mapping_dir = os.path.join(dir_path,
                                                   'crystal_data/space_group_lattice_mapping.json')
    with open(space_group_lattice_mapping_dir, 'r') as f:
        space_group_lattice_mapping = json.load(f, object_hook=json_key_to_int)

    def __init__(self, space_group_number=1, print_info=False):
        if space_group_number <= 0 or space_group_number > 230:
            raise ValueError('space_group_number must be an integer between 1 and 230')

        self.space_group_number = space_group_number
        self.lattice_type = self.space_group_lattice_mapping[self.space_group_number]
        self.lattice = lattice_system_dict_3D[self.lattice_type]
        info = spg.get_symmetry_from_database(self.space_group_hall_mapping[space_group_number])
        self.translations = info['translations']
        self.rotations = info['rotations']

        if print_info:
            print('Space group number: {}\n'.format(space_group_number),
                  'lattice type: {}\n'.format(self.lattice_type),
                  'Default parameters for lattice: {}'.format(self.lattice.lattice_params))

    def get_basis_vectors(self, base_positions, base_type=[], base_quaternions=None,
                          is_complete=False, apply_orientation=False):
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
            bool value to indicate if the positions are complete postions in a unitcell
        :type is_complete:
            bool
        :param apply_orientations:
            bool value to indicate if the space group symmetry should be applied to orientatioin
        :type apply_orientations:
            bool
        :return:
            basis_vectors
        :rtype:
            np.ndarray
        """

        # check input accuracy
        if not isinstance(base_positions, np.ndarray) or len(base_positions.shape) == 1 \
           or base_positions.shape[1] != 3:
            raise ValueError('base_positions must be an numpy array of shape Nx3')
        if apply_orientation:
            if not isinstance(base_quaternions, np.ndarray) or len(base_quaternions.shape) == 1 \
               or base_quaternions.shape[1] != 4:
                raise ValueError('base_quaternions must be an numpy array of shape Nx4')
        if len(base_type):
            if not isinstance(base_type, list) or len(base_type) != base_positions.shape[0] \
               or not all(isinstance(i, str) for i in base_type):
                raise ValueError('base_type must contain a list of type name the same length'
                                 'as the number of basis positions')
        else:
            base_type = ['A'] * base_positions.shape[0]

        threshold = 1e-6
        reflection_exist = False
        for i in range(0, len(self.rotations)):
            # Generate the new set of positions from the base
            pos = wrap(self.rotations[i].dot(base_positions.T).T + self.translations[i])
            if apply_orientation:
                if np.linalg.det(self.rotations[i]) == 1:
                    quat_rotate = rowan.from_matrix(self.rotations[i], require_orthogonal=False)
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
                        print('Orientation quaterions may have multiple values for the same '
                              'particle postion under the symmetry operation for this space group '
                              'and is not well defined, only the first occurance is used.')
        if reflection_exist:
            print('Warning: reflection operation is included in this space group, '
                  'and is ignored for quaternion calculation.')

        if is_complete and len(positions) != len(base_positions):
            raise ValueError('the complete basis postions vector does not match the space group '
                             'chosen. More positions are generated based on the symmetry operation '
                             'within the provided space group')

        if apply_orientation:
            return wrap(positions), type_list, quaternions
        else:
            return wrap(positions), type_list

    def get_lattice_vectors(self, **user_lattice_params):
        """Initialize the unitcell and return lattice vectors [a1, a2, a3].

        :param user_lattice_params:
            unit cell parameters, provide a, b, c, alpha, beta, gamma where applicable
        :type user_lattice_params:
            float
        :return: lattice_vectors
        :rtype: np.ndarray
        """
        return self.lattice.get_lattice_vectors(**user_lattice_params)


class Prototype(object):
    """Crystal prototype class.

    This class uses the minimal necessay information needed to fully define a crystal structures
    with space group number, wyckoff postions(in letter name convention) and free parameters for
    each relavent wyckoff postions.

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
    :param print_info:
        print space group information upon initialization
    :type print_info:
        bool
    """

    dir_path = os.path.dirname(__file__)

    def __init__(self, space_group_number=1, wyckoff_site='', type_by_site='', print_info=False):
        if space_group_number > 230 or space_group_number < 1:
            raise ValueError('space_group_number must be an integer between 0 and 230, '
                             'default = 1')

        if not isinstance(wyckoff_site, str) or not wyckoff_site.isalpha():
            raise ValueError('wyckoff_postions must be string consists of all the Wyckoff '
                             'postions letters, e.g. \'abcc\' denotes one set of Wyckoff postions '
                             'for both a and b, and two sets at Wyckoff postion c')

        if type_by_site == '':
            type_by_site = 'A' * len(wyckoff_site)
        elif not isinstance(type_by_site, str) or len(type_by_site) != len(wyckoff_site) or \
             not type_by_site.isalpha():
            raise ValueError('type_by_site must be string consists of type name (A/B/C, etc)'
                             'for each Wyckoff site, default all particles with same type A '
                             'if not provided')

        wyckoff_site_list = list(wyckoff_site.lower())
        type_by_site = list(type_by_site.upper())

        wyckoff_data_dir = os.path.join(self.dir_path, 'crystal_data/space_group_{}_Wyckoff'
                                        '_site_data.json'.format(space_group_number))
        with open(wyckoff_data_dir, 'r') as f:
            full_wyckoff_positions = json.load(f)

        basis_params_list = []
        order = 1
        for site in wyckoff_site_list:
            pos = copy.deepcopy(full_wyckoff_positions[site])
            pos = ''.join(pos)
            for letter in ('x', 'y', 'z'):
                if letter in pos:
                    basis_params_list.append(letter + str(order))
            order += 1
        basis_params_value_list = [None] * len(basis_params_list)

        self.space_group_number = space_group_number
        self.space_group = SpaceGroup(space_group_number, print_info)
        self.wyckoff_site_list = wyckoff_site_list
        self.full_wyckoff_positions = full_wyckoff_positions
        self.type_by_site = type_by_site
        self.lattice_params = self.space_group.lattice.lattice_params
        self.basis_params = dict(zip(basis_params_list, basis_params_value_list))

        if print_info:
            print('Wyckoff sites:{}\n'.format(wyckoff_site_list),
                  'Particle type for each Wyckoff sites:{}\n'.format(type_by_site),
                  'lattice parameters list:{}\n'.format(list(self.lattice_params)),
                  'basis parameters list:{}'.format(basis_params_list))

    def update_basis_params(self, user_basis_params):
        params = copy.deepcopy(self.basis_params)
        for param, value in user_basis_params.items():
            if param in params:
                if value is not None:
                    params[param] = value
            else:
                print('warning: {} is not required and not used to define this '
                      'structure'.format(param))
        for value in params.values():
            if value is None:
                raise ValueError('not all necessary parameters were provided! Set print_info=True '
                                 'at prototype initialization to see the full list of necessary '
                                 'parameters')
        return params

    def get_basis_vectors(self, **user_basis_params):
        """Initialize fractional coordinates of the particles in the unitcell.

        :param user_basis_params:
            user defined parameters for different Wyckoff site degree of freedom, when applicable
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
            for letter in ('x', 'y', 'z'):
                if letter + str(order) in basis_params.keys():
                    exec('{} = {}'.format(letter, basis_params[letter + str(order)]))
            for i in range(0, 3):
                # add * back for eval
                target = re.findall(r"(\d[xyz])", pos[i])
                for item in target:
                    pos[i] = pos[i].replace(item,
                                            '*'.join(re.findall(r"(\d)([xyz])", item)[0]))
                pos[i] = eval(pos[i])
            base_positions = np.append(base_positions, np.array(pos).reshape(1, -1), axis=0)
            order += 1

        return self.space_group.get_basis_vectors(wrap(base_positions), base_type=self.type_by_site)

    def update_lattice_params(self, user_lattice_params):
        params = copy.deepcopy(self.lattice_params)
        for param, value in user_lattice_params.items():
            if param in params:
                if value is not None:
                    params[param] = value
            else:
                print('warning: {} is not required and not used to define this '
                      'structure'.format(param))
        return params

    def get_lattice_vectors(self, **user_lattice_params):
        """Initialize the unitcell and return lattice vectors [a1, a2, a3]

        :param user_lattice_params:
            unit cell parameters, provide a, b, c, alpha, beta, gamma where applicable
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

    This class uses the crystal prototypes in Aflow database to initialize crystal structures.

    :param prototype_index:
        prototype index [0, 589] for all 590 prototypes in AFLOW.
    :type prototype_index:
        int
    :param set_type:
        allow setting different type name(in A, B, C order) for different atoms in AFLOW prototype
    :type set_type:
        bool
    :param print_info:
        Print prototype information upon initialization
    :type print_info:
        bool
    """

    dir_path = os.path.dirname(__file__)
    Aflow_data_dir = os.path.join(dir_path,
                                  'crystal_data/Aflow_processed_data.csv')
    Aflow_database = pd.read_csv(Aflow_data_dir, index_col=0)

    def __init__(self, prototype_index=0, set_type=False, print_info=True):
        # should do search, return best match and all options, define a Structure
        # name search should support one or any combination of: Pearson symbol, space_group number
        # and chemistry
        # must define unitcell type and space_group now
        if prototype_index < 0 or prototype_index >= 590:
            raise ValueError('prototype_index must be an integer between 0 and 590, default '
                             'value of 0 will skip search by index and use Pearson symbol or '
                             'chemistry for search')
        if prototype_index + 1:
            entry = self.Aflow_database.iloc[prototype_index]
        # TODO: a search and use best match feature with pearson and chemistry input

        space_group_number = entry['space_group_number']
        lattice_params_list = re.findall(r"'(.*?)'", entry['lattice_params_list'])
        try:
            lattice_params_value_list = [float(i) for i in
                                         entry['lattice_params_value_list'].strip('[]').split(',')]
        except BaseException:
            lattice_params_value_list = []

        basis_params_list = re.findall(r"'(.*?)'", entry['basis_params_list'])
        try:
            basis_params_value_list = [float(i) for i in
                                       entry['basis_params_value_list'].strip('[]').split(',')]
        except BaseException:
            basis_params_value_list = []

        lattice_params = dict(zip(lattice_params_list, lattice_params_value_list))
        basis_params = dict(zip(basis_params_list, basis_params_value_list))

        # convert Aflow angle unit from degree to rad
        for key in ['alpha', 'beta', 'gamma']:
            if key in lattice_params.keys():
                lattice_params[key] = lattice_params[key] / 180 * np.pi

        # process proper unitcell params
        if space_group_number in [146, 148, 155, 160, 161, 166, 167]:
            a = lattice_params.pop('a')
            c = lattice_params.pop('c/a') * a
            lattice_params['a'] = np.sqrt(a ** 2 / 3 + c ** 2 / 9)
            lattice_params['alpha'] = np.arccos((2 * c ** 2 - 3 * a ** 2) / (
                                                2 * (c ** 2 + 3 * a ** 2)))
        else:
            # for others
            a = lattice_params['a']
            try:
                lattice_params['b'] = lattice_params.pop('b/a') * a
            except BaseException:
                pass
            try:
                lattice_params['c'] = lattice_params.pop('c/a') * a
            except BaseException:
                pass

        wyckoff_site_list_by_type = re.findall(r"'(.*?)'", entry['Wyckoff_site'])
        wyckoff_site_list = list(''.join(wyckoff_site_list_by_type))
        wyckoff_site_list.sort()

        wyckoff_data_dir = os.path.join(self.dir_path, 'crystal_data/space_group_{}_Wyckoff'
                                        '_site_data.json'.format(space_group_number))

        with open(wyckoff_data_dir, 'r') as f:
            full_wyckoff_positions = json.load(f)

        # get type label
        type_by_site = list('A' * len(wyckoff_site_list))
        if set_type:
            sorted_site_string = ''.join(wyckoff_site_list)
            base = ord('A')
            for wyckoffs in wyckoff_site_list_by_type:
                for site in wyckoffs:
                    order = sorted_site_string.find(site)
                    sorted_site_string = sorted_site_string.replace(site, '0', 1)
                    type_by_site[order] = chr(base)
                base += 1

        self.space_group_number = space_group_number
        self.space_group = SpaceGroup(space_group_number)
        self.wyckoff_site_list = wyckoff_site_list
        self.full_wyckoff_positions = full_wyckoff_positions
        self.type_by_site = type_by_site
        self.lattice_params = lattice_params
        self.basis_params = basis_params

        if print_info:
            print('Info for the chosen crystal structure prototype:\n',
                  'id: {}, (Pearson-Chemistry-SpaceGroup)\n'.format(entry['id']),
                  'Wyckoff sites: {}\n'.format(entry['Wyckoff_site']),
                  'available lattice parameters: {}\n'.format(lattice_params),
                  'available basis parameters: {}'.format(basis_params))


# Point Group operations
class PointGroup(object):
    """A class to access all point group symmetry operations.

    This class provides method to access all point group symmetry operation in both rotational
    matrix form or quaternion form.

    :param point_group_number:
        Point group number between 1 and 32.
    :type point_group_number:
        int
    :param print_info:
        Print point group information upon initialization.
    :type print_info:
        bool
    """

    dir_path = os.path.dirname(__file__)
    point_group_rotation_matrix_dir = os.path.join(dir_path,
                                                   'crystal_data/'
                                                   'point_group_rotation_matrix_dict.pickle')
    with open(point_group_rotation_matrix_dir, 'rb') as f:
        point_group_rotation_matrix_dict = pickle.load(f)

    point_group_quat_dir = os.path.join(dir_path,
                                            'crystal_data/point_group_quat_dict.json')
    with open(point_group_quat_dir, 'r') as f:
        point_group_quat_dict = json.load(f, object_hook=json_key_to_int)

    point_group_name_mapping_dir = os.path.join(dir_path,
                                                   'crystal_data/point_group_name_mapping.json')
    with open(point_group_name_mapping_dir, 'r') as f:
        point_group_name_mapping = json.load(f, object_hook=json_key_to_int)

    def __init__(self, point_group_number=1, print_info=False):
        if point_group_number <= 0 or point_group_number > 32:
            raise ValueError('point_group_number must be an integer between 1 and 32')

        self.point_group_number = point_group_number
        self.point_group_name = self.point_group_name_mapping[point_group_number]
        self.rotation_matrix = \
            self.point_group_rotation_matrix_dict[point_group_number]['rotations']
        self.quaternion = self.point_group_quat_dict[point_group_number]

        if print_info:
            print('Point group number: {}\n'.format(point_group_number),
                  'Name {}'.format(self.point_group_name))

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
