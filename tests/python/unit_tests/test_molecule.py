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

from rovib import Molecule


class TestMolecule(unittest.TestCase):
    def setUp(self):
        self.symbols = ("O", "H", "H")
        self.coords = np.arange(9.0).reshape(3, 3)
        self.masses = np.array([16.0, 1.0, 1.0])
        self.mol = Molecule(self.symbols, self.coords, self.masses)

    def test_state(self):
        self.assertEqual(self.mol.symbols, self.symbols)
        np.testing.assert_array_equal(self.mol.coords, self.coords)
        np.testing.assert_array_equal(self.mol.masses, self.masses)

    def test_natoms(self):
        self.assertEqual(self.mol.natoms, 3)

    def test_nvib(self):
        self.assertEqual(self.mol.nvib, 3)
        four_atoms = Molecule("ABCD", np.zeros((4, 3)), np.ones(4))
        self.assertEqual(four_atoms.nvib, 6)

    def test_arrays_are_read_only(self):
        with self.assertRaises(ValueError):
            self.mol.coords[0, 0] = 1.0
        with self.assertRaises(ValueError):
            self.mol.masses[0] = 1.0

    def test_input_is_copied(self):
        self.coords[0, 0] = 100.0
        self.assertEqual(self.mol.coords[0, 0], 0.0)

    def test_bad_coords_shape(self):
        with self.assertRaises(ValueError):
            Molecule(self.symbols, np.zeros((2, 3)), self.masses)
        with self.assertRaises(ValueError):
            Molecule(self.symbols, np.zeros(9), self.masses)

    def test_bad_masses_shape(self):
        with self.assertRaises(ValueError):
            Molecule(self.symbols, self.coords, np.ones(2))

    def test_equality(self):
        same = Molecule(self.symbols, self.coords, self.masses)
        self.assertEqual(self.mol, same)

        moved = Molecule(self.symbols, self.coords + 1.0, self.masses)
        self.assertNotEqual(self.mol, moved)

        heavier = Molecule(self.symbols, self.coords, self.masses * 2.0)
        self.assertNotEqual(self.mol, heavier)

        renamed = Molecule(("O", "H", "D"), self.coords, self.masses)
        self.assertNotEqual(self.mol, renamed)

        self.assertNotEqual(self.mol, "not a molecule")
