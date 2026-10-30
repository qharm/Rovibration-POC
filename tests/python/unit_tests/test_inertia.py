# Copyright 2026 QHARM Team
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import unittest

import numpy as np

from rovib import (
    Molecule,
    inertia_tensor,
    rotational_constants,
    water_molecule,
)
from rovib.constants import HARTREE_TO_CM


def _planar_molecule():
    """Three atoms in the yz-plane whose principal axes are x, y, and z.

    The moments of inertia are I_xx = 14, I_yy = 12, and I_zz = 2.
    """
    coords = [[0.0, 0.0, -1.0], [0.0, -1.0, 2.0], [0.0, 1.0, 2.0]]
    masses = [4.0, 1.0, 1.0]
    return Molecule(("X", "Y", "Y"), coords, masses)


class TestInertiaTensor(unittest.TestCase):
    def test_planar(self):
        np.testing.assert_allclose(
            inertia_tensor(_planar_molecule()), np.diag([14.0, 12.0, 2.0])
        )

    def test_off_diagonal(self):
        mol = Molecule(("X",), [[1.0, 2.0, 3.0]], [2.0])
        expected = 2.0 * np.array(
            [[13.0, -2.0, -3.0], [-2.0, 10.0, -6.0], [-3.0, -6.0, 5.0]]
        )
        np.testing.assert_allclose(inertia_tensor(mol), expected)


class TestRotationalConstants(unittest.TestCase):
    def setUp(self):
        # The planar molecule's rotational constants in ascending order.
        self.ascending = [1.0 / 28.0, 1.0 / 24.0, 1.0 / 4.0]

    def test_default_ordering(self):
        # x is the smallest, y the largest, and z the middle one.
        np.testing.assert_allclose(
            rotational_constants(_planar_molecule()),
            [1.0 / 28.0, 1.0 / 4.0, 1.0 / 24.0],
        )

    def test_ordering(self):
        mol = _planar_molecule()
        np.testing.assert_allclose(
            rotational_constants(mol, (1, 2, 3)), self.ascending
        )
        np.testing.assert_allclose(
            rotational_constants(mol, (3, 2, 1)), self.ascending[::-1]
        )

    def test_string_ordering(self):
        mol = _planar_molecule()
        np.testing.assert_array_equal(
            rotational_constants(mol, "231"),
            rotational_constants(mol, (2, 3, 1)),
        )

    def test_bad_ordering(self):
        mol = _planar_molecule()
        for ordering in ((1, 2), (1, 1, 2), (0, 1, 2), (1, 2, 4), "1234"):
            with self.subTest(ordering=ordering):
                with self.assertRaises(ValueError):
                    rotational_constants(mol, ordering)

    def test_water(self):
        # The paper's rotational constants (x, y, z = a, b, c) in cm^-1.
        constants = rotational_constants(water_molecule()) * HARTREE_TO_CM
        np.testing.assert_allclose(constants, [9.49, 27.24, 14.57], atol=0.01)

    def test_rotated(self):
        # The rotational constants do not depend on the orientation.
        water = water_molecule()
        angle = np.radians(10.0)
        rotation = np.array(
            [
                [1.0, 0.0, 0.0],
                [0.0, np.cos(angle), -np.sin(angle)],
                [0.0, np.sin(angle), np.cos(angle)],
            ]
        )
        rotated = Molecule(
            water.symbols, water.coords @ rotation.T, water.masses
        )
        np.testing.assert_allclose(
            rotational_constants(rotated), rotational_constants(water)
        )

    def test_linear(self):
        masses = [16.0, 12.0, 16.0]
        on_axis = [[0.0, 0.0, -1.0], [0.0, 0.0, 0.0], [0.0, 0.0, 1.0]]
        off_axis = np.outer([-1.0, 0.0, 1.0], [1.0, 1.0, 1.0]) / np.sqrt(3.0)
        for coords in (on_axis, off_axis):
            mol = Molecule(("O", "C", "O"), coords, masses)
            with self.assertRaises(ValueError):
                rotational_constants(mol)
