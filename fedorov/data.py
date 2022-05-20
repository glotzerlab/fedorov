# Copyright (c) 2019-2020 The Regents of the University of Michigan
# This file is part of the fedorov project, released under the BSD 3-Clause
# License.

import os

import numpy as np

_DATA_PATH = os.path.join(os.path.dirname(__file__), "crystal_data")


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
    """Convert box parameter: Lx, Ly, Lz, xy, xz, yz to lattice vectors.

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
    lattice_vectors = np.array(
        [[Lx, 0, 0], [xy * Ly, Ly, 0], [xz * Lz, yz * Lz, Lz]]
    )
    return lattice_vectors


def get_volume(lattice_vectors):
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


def translate_to_vector(
    a=1, b=1, c=1, alpha=np.pi / 2, beta=np.pi / 2, gamma=np.pi / 2
):
    """Convert box parameters a, b, c, alpha, beta, gamma to lattice vectors.

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
        raise ValueError(
            "Error: the box length and angle parameters provided are not "
            "feasible. Please not the unit used for angle paramters should be "
            "in unit of rad"
        )
    cz = np.sqrt(1 - ca * ca - cb * cb - cg * cg + 2 * ca * cb * cg) / sg
    lattice_vectors = np.array(
        [[a, 0, 0], [b * cg, b * sg, 0], [c * cb, c * cy, c * cz]]
    )
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
    lattice_vectors = np.array([[a, 0], [b * np.cos(theta), b * np.sin(theta)]])
    return lattice_vectors


__all__ = [
    "wrap",
    "convert_to_box",
    "convert_to_vectors",
    "translate_to_vector",
    "fractional_to_cartesian",
    "get_volume",
    "translate_to_vector_2D",
]
