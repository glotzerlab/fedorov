# pytest
import json
import os
import re

import numpy as np
import pandas as pd

import fedorov
from fedorov import (
    AflowPrototype,
    PlaneGroup,
    PointGroup,
    Prototype,
    SpaceGroup,
    convert_to_box,
    convert_to_vectors,
    translate_to_vector,
    wrap,
)


# test all available structures from Aflow database
def test_aflow_prototype():
    dir_path = os.path.dirname(fedorov.__file__)
    Aflow_data_dir = os.path.join(dir_path, "crystal_data/Aflow_raw_data.csv")
    database = pd.read_csv(Aflow_data_dir, index_col=0)

    # test all run
    for i in range(0, 590):
        structure_test = AflowPrototype(prototype_index=i, print_info=False)
        basis_vectors, type_list = structure_test.get_basis_vectors()
        lattice_vectors = structure_test.get_lattice_vectors()
        N = int(re.findall(r"\d+", database.loc[i, "Pearson Symbol"])[0])
        assert not np.isnan(lattice_vectors).any()
        assert len(basis_vectors) == N
        assert len(type_list) == N
        assert "".join(type_list).isalpha()

    # detailed test result for different case
    # hexagonal
    structure_test = AflowPrototype(prototype_index=0, print_info=False)
    basis_vectors, type_list = structure_test.get_basis_vectors()
    lattice_vectors = structure_test.get_lattice_vectors()
    lattice_vectors_reference = np.array(
        [[2.508, 0, 0], [-1.254, 2.1719917, 0], [0, 0, 4.183]]
    )
    basis_vectors_reference = np.array(
        [
            [0.33333333333333, 0.66666666666667, 0.05995000000000],
            [0.66666666666667, 0.33333333333333, -0.05995000000000],
            [0.66666666666667, 0.33333333333333, 0.55995000000000],
            [0.33333333333333, 0.66666666666667, 0.44005000000000],
        ]
    )
    assert np.allclose(lattice_vectors, lattice_vectors_reference)
    assert np.allclose(basis_vectors, wrap(basis_vectors_reference))

    # Rhombohedral
    structure_test = AflowPrototype(prototype_index=10, print_info=False)
    basis_vectors, type_list = structure_test.get_basis_vectors()
    lattice_vectors = structure_test.get_lattice_vectors()
    basis_vectors_reference = np.array(
        [
            [0, 0, 0],
            [0.5, 0.5, 0.5],
            [0.2667, 0.2667, 0.2667],
            [-0.2667, -0.2667, -0.2667],
        ]
    )
    a = 6.773648
    alpha = 0.531214
    lattice_vectors_reference = translate_to_vector(
        a=a, b=a, c=a, alpha=alpha, beta=alpha, gamma=alpha
    )
    assert np.allclose(lattice_vectors, lattice_vectors_reference)
    assert np.allclose(basis_vectors, wrap(basis_vectors_reference))

    # tetragonal and type check
    structure_test = AflowPrototype(
        prototype_index=5, print_info=False, set_type=True
    )
    basis_vectors, type_list = structure_test.get_basis_vectors()
    lattice_vectors = structure_test.get_lattice_vectors()
    Lx, Ly, Lz, xy, xz, yz = convert_to_box(lattice_vectors)
    basis_vectors_reference = np.array(
        [
            [0, 0, 0],
            [0.5, 0.5, 0.8973],
            [0.5, 0.5, 0.4517],
            [0.5, 0, 0.3785],
            [0, 0.5, 0.3785],
        ]
    )
    lattice_vectors_reference = np.array(
        [[4.046, 0, 0], [0, 4.046, 0], [0, 0, 4.1394]]
    )
    assert np.allclose(lattice_vectors, lattice_vectors_reference)
    assert np.allclose(basis_vectors, wrap(basis_vectors_reference))
    assert type_list == ["B", "A", "C", "A", "A"]


def test_prototype():
    structure_test = Prototype(
        space_group_number=230, wyckoff_site="hh", type_by_site="ac"
    )
    basis_params = {
        "x1": 0.12,
        "y1": 0.13,
        "z1": 0.14,
        "x2": -0.125,
        "y2": -0.135,
        "z2": -0.145,
    }

    basis_vectors, type_list = structure_test.get_basis_vectors(**basis_params)
    lattice_vectors = structure_test.get_lattice_vectors(a=4)
    assert type_list == ["A", "C"] * 96
    assert np.allclose(basis_vectors[0, :], np.array([0.12, 0.13, 0.14]))
    assert np.allclose(
        basis_vectors[1, :], wrap(np.array([-0.125, -0.135, -0.145]))
    )
    assert np.allclose(
        lattice_vectors, np.array([[4, 0, 0], [0, 4, 0], [0, 0, 4]])
    )


# test math for basic crystal calculation methods
def test_fractional_coordinates():
    # test_fractional_coordinates:
    lattice_vectors = np.array(
        [
            [1, 0, 0],
            [1.08060461, 1.68294197, 0],
            [2.19506661, 1.7193084, 1.10709584],
        ]
    )
    assert np.allclose(
        translate_to_vector(1, 2, 3, 0.5, 0.75, 1.0), lattice_vectors
    )


def test_vector_to_box():
    lattice_vectors = np.array(
        [
            [1, 0, 0],
            [1.08060461, 1.68294197, 0],
            [2.19506661, 1.7193084, 1.10709584],
        ]
    )
    assert np.allclose(
        lattice_vectors, convert_to_vectors(*convert_to_box(lattice_vectors))
    )


# test generation from SpaceGroup
def test_space_group():
    spg_test = SpaceGroup(220)
    basis_positions = np.array([[0.1, 0.12, 0.13], [0.14, 0.15, 0.17]])
    base_quaternions = np.array([[0, 0, 0, 1], [0, 0, 1, 0]])
    basis_vectors, type_list, quaternions = spg_test.get_basis_vectors(
        basis_positions,
        base_type=["B", "A"],
        base_quaternions=base_quaternions,
        apply_orientation=True,
    )
    lattice_vectors = spg_test.get_lattice_vectors()
    basis_vector_last = np.array([[0.61, 0.92, 0.1]])
    assert np.allclose(basis_vectors[-1, :], basis_vector_last)
    assert basis_vectors.shape[0] == 96
    assert np.allclose(
        lattice_vectors, np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
    )
    assert type_list == ["B", "A"] * 48
    # TODO test quaternion accuracy


# test all values against structures_for_test.txt
def test_prototype_accuracy():
    reference_file = os.path.join(
        os.path.dirname(__file__), "structures_for_test.txt"
    )
    print(reference_file)
    with open(reference_file, "r") as f:
        for line in f:
            struct = json.loads(line)
            spg_num = struct["space_group_number"]
            wyckoff_site = struct["Wyckoff_site"]
            test_prototype = Prototype(
                space_group_number=spg_num, wyckoff_site=wyckoff_site
            )
            basis_vectors, type_list = test_prototype.get_basis_vectors(
                **struct["basis_params"]
            )
            lattice_vectors = test_prototype.get_lattice_vectors(
                **struct["lattice_params"]
            )

            assert np.allclose(
                lattice_vectors, np.asarray(struct["lattice_vectors"])
            )
            pos_comparisons = (
                basis_vectors
                - wrap(np.asarray(struct["basis_vectors"]))[:, np.newaxis, :]
            )
            norms = np.linalg.norm(pos_comparisons, axis=2)
            norms = np.amin(norms, axis=1)
            assert np.linalg.norm(norms) < 1e-6


# test aflow_database parameter is complete
def test_aflow_database_accuracy():
    for i, name in AflowPrototype.Aflow_database["id"].iteritems():
        cdbs = AflowPrototype(i, print_info=False)
        keys = set(cdbs.lattice_params.keys())
        if name.startswith("a"):  # triclinic
            assert keys == set(("a", "b", "c", "alpha", "beta", "gamma"))
        elif name.startswith("m"):  # monoclinic
            assert keys == set(("a", "b", "c", "beta"))
        elif name.startswith("o"):  # orthorhombic
            assert keys == set(("a", "b", "c"))
        elif name.startswith("t"):  # tetragonal
            assert keys == set(("a", "c"))
        elif name.startswith("hR"):  # rhombohedral
            assert keys == set(("a", "alpha"))
        elif name.startswith("h"):  # hexagonal
            assert keys == set(("a", "c"))
        elif name.startswith("c"):  # cubic
            assert keys == set(("a"))
        else:
            raise (name + "not start with a Pearson symbol")


# test generation from plane group
def test_plane_group():
    pg_test = PlaneGroup(9)
    basis_positions = np.array([[0.1, 0.12], [0.14, 0.15]])
    base_quaternions = np.array([[0, 0, 0, 1], [0, 0, 1, 0]])
    basis_vectors, type_list, quaternions = pg_test.get_basis_vectors(
        basis_positions,
        base_type=["B", "A"],
        base_quaternions=base_quaternions,
        apply_orientation=True,
    )
    lattice_vectors = pg_test.get_lattice_vectors(a=1, b=2)
    basis_vector_last = np.array([0.64, 0.35])
    assert np.allclose(basis_vectors[-1, :], basis_vector_last)
    assert basis_vectors.shape[0] == 16
    assert np.allclose(lattice_vectors, np.array([[1, 0], [0, 2]]))
    assert type_list == ["B", "A"] * 8
    # TODO test quaternion accuracy


# test generation from point group symmetry access
def test_point_group():
    point_group_test = PointGroup(29)
    assert point_group_test.point_group_name == "m-3"
    assert point_group_test.get_quaternion()[2] == [0.0, 0.0, -1.0, 0.0]
    assert np.allclose(
        point_group_test.get_rotation_matrix()[2],
        np.array([[-1, 0, 0], [0, 1, 0], [0, 0, -1]]),
    )
