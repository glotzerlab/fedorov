# Copyright (c) 2019-2020 The Regents of the University of Michigan
# This file is part of the fedorov project, released under the BSD 3-Clause
# License.

import copy
import warnings

import numpy as np

from . import data


class Lattice:
    """Base class for 2D and 3D lattices."""

    pass


class Oblique2D:
    """A class for constructing a 2D oblique unitcell

    This class provides method to initialize a 2D oblique unitcell
    """

    lattice_params = {"a": 1, "b": 1, "theta": np.pi / 2}

    @classmethod
    def update_lattice_params(cls, user_lattice_params):
        params = copy.deepcopy(cls.lattice_params)
        for param, value in user_lattice_params.items():
            if param in params:
                if value is not None:
                    params[param] = value
            else:
                warnings.warn(
                    f"warning: {param} is not used to define this structure",
                    RuntimeWarning,
                )
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
        lattice_vectors = data.translate_to_vector_2D(**params)
        return lattice_vectors


class Rectangular2D(Oblique2D):
    """A class for constructing a 2D rectangular unitcell

    This class provides method to initialize a 2D rectangular unitcell
    """

    lattice_params = {"a": 1, "b": 1}

    @classmethod
    def get_lattice_vectors(cls, **user_lattice_params):
        """Initialize a 2D rectangular unitcell and return lattice vectors.

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
        lattice_vectors = np.array([[params["a"], 0.0], [0.0, params["b"]]])
        return lattice_vectors


class Hexagonal2D(Oblique2D):
    """A class for constructing a 2D hexagonal unitcell

    This class provides method to initialize a 2D hexagonal unitcell
    """

    lattice_params = {"a": 1}

    @classmethod
    def get_lattice_vectors(cls, **user_lattice_params):
        """Initialize a 2D hexagonal unitcell and return lattice vectors.

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
        lattice_vectors = np.array(
            [
                [params["a"], 0.0],
                [-0.5 * params["a"], params["a"] * np.sqrt(3) / 2],
            ]
        )
        return lattice_vectors


class Square2D(Rectangular2D):
    """A class for constructing a 2D square unitcell

    This class provides method to initialize a 2D square unitcell
    """

    lattice_params = {"a": 1}

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
        lattice_vectors = np.array([[params["a"], 0.0], [0.0, params["a"]]])
        return lattice_vectors


lattice_system_dict_2D = {
    "oblique": Oblique2D,
    "rectangular": Rectangular2D,
    "hexagonal": Hexagonal2D,
    "square": Square2D,
}


# 3D systems
class Triclinic:
    """A class for constructing a triclinic unitcell."""

    Pearson = "a"
    dimensions = 3
    lattice_params = {
        "a": 1,
        "b": 1,
        "c": 1,
        "alpha": np.pi / 2,
        "beta": np.pi / 2,
        "gamma": np.pi / 2,
    }

    @classmethod
    def update_lattice_params(cls, user_lattice_params):
        params = copy.deepcopy(cls.lattice_params)
        for param, value in user_lattice_params.items():
            if param in params:
                if value is not None:
                    params[param] = value
            else:
                warnings.warn(
                    f"warning: {param} is not used to define this structure",
                    RuntimeWarning,
                )
        return params

    @classmethod
    def get_lattice_vectors(cls, **user_lattice_params):
        """Initialize a triclinic unitcell and return lattice vectors.

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
        params = cls.update_lattice_params(user_lattice_params)
        lattice_vectors = data.translate_to_vector(**params)
        return lattice_vectors


class Monoclinic(Triclinic):
    """A class for constructing a monoclinic unitcell

    This class provides method to initialize a monoclinic unitcell
    """

    Pearson = "m"
    lattice_params = {"a": 1, "b": 1, "c": 1, "beta": np.pi / 2}

    @classmethod
    def get_lattice_vectors(cls, **user_lattice_params):
        """Initialize a monoclinic unitcell and return lattice vectors.

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
        params = cls.update_lattice_params(user_lattice_params)
        lattice_vectors = data.translate_to_vector(
            a=params["a"], b=params["b"], c=params["c"], beta=params["beta"]
        )
        return lattice_vectors


class Orthorhombic(Monoclinic):
    """A class for constructing a orthorhombic unitcell."""

    Pearson = "o"
    lattice_params = {"a": 1, "b": 1, "c": 1}

    @classmethod
    def get_lattice_vectors(cls, **user_lattice_params):
        """Initialize a orthorhombi unitcell and return lattice vectors.

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
        params = cls.update_lattice_params(user_lattice_params)
        lattice_vectors = np.array(
            [
                [params["a"], 0.0, 0.0],
                [0.0, params["b"], 0.0],
                [0.0, 0.0, params["c"]],
            ]
        )
        return lattice_vectors


class Tetragonal(Orthorhombic):
    """A class for constructing a tetragonal unitcell."""

    Pearson = "t"
    lattice_params = {"a": 1, "c": 1}

    @classmethod
    def get_lattice_vectors(cls, **user_lattice_params):
        """Initialize a tetragona unitcell and return lattice vectors.

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
        params = cls.update_lattice_params(user_lattice_params)
        lattice_vectors = np.array(
            [
                [params["a"], 0.0, 0.0],
                [0.0, params["a"], 0.0],
                [0.0, 0.0, params["c"]],
            ]
        )
        return lattice_vectors


class Hexagonal(Triclinic):
    """A class for constructing a hexagonal unitcell."""

    Pearson = "hP"
    lattice_params = {"a": 1, "c": 1}

    @classmethod
    def get_lattice_vectors(cls, **user_lattice_params):
        """Initialize a hexagonal unitcell and return lattice vectors.

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
        params = cls.update_lattice_params(user_lattice_params)
        lattice_vectors = np.array(
            [
                [params["a"], 0.0, 0.0],
                [-0.5 * params["a"], np.sqrt(3.0) / 2.0 * params["a"], 0.0],
                [0.0, 0.0, params["c"]],
            ]
        )
        return lattice_vectors


class Rhombohedral(Triclinic):
    """A class for constructing a rhombohedral unitcell."""

    Pearson = "hR"
    lattice_params = {"a": 1, "alpha": np.pi / 2}

    @classmethod
    def get_lattice_vectors(cls, **user_lattice_params):
        """Initialize a rhombohedral unitcell and return lattice vectors.

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
        params = cls.update_lattice_params(user_lattice_params)
        lattice_vectors = data.translate_to_vector(
            a=params["a"],
            b=params["a"],
            c=params["a"],
            alpha=params["alpha"],
            beta=params["alpha"],
            gamma=params["alpha"],
        )
        return lattice_vectors


class Cubic(Tetragonal):
    """A class for constructing a cubic unitcell."""

    Pearson = "c"
    lattice_params = {"a": 1}

    @classmethod
    def get_lattice_vectors(cls, **user_lattice_params):
        """Initialize a cubicc unitcell and return lattice vectors.

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
        params = cls.update_lattice_params(user_lattice_params)
        lattice_vectors = np.array(
            [
                [params["a"], 0.0, 0.0],
                [0.0, params["a"], 0.0],
                [0.0, 0.0, params["a"]],
            ]
        )
        return lattice_vectors


lattice_system_dict_3D = {
    "triclinic": Triclinic,
    "monoclinic": Monoclinic,
    "orthorhombic": Orthorhombic,
    "tetragonal": Tetragonal,
    "hexagonal": Hexagonal,
    "rhombohedral": Rhombohedral,
    "cubic": Cubic,
}
