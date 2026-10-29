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
from rovib.constants import HARTREE_TO_CM
from rovib.water_params import (
    harmonic_analysis,
    is_water,
    water_harmonic_frequencies,
    water_molecule,
    water_normal_modes,
)

#: The normal modes are only given to four digits.
L_TOL = 2.0e-3


class TestWaterMolecule(unittest.TestCase):
    def setUp(self):
        self.water = water_molecule()

    def test_atoms(self):
        self.assertEqual(self.water.symbols, ("O", "H", "H"))
        self.assertEqual(self.water.nvib, 3)

    def test_isotopic_masses(self):
        masses_in_amu = self.water.masses / 1822.888486209
        np.testing.assert_allclose(
            masses_in_amu, [15.99491461957, 1.00782503223, 1.00782503223]
        )

    def test_center_of_mass(self):
        com = self.water.masses @ self.water.coords / self.water.masses.sum()
        np.testing.assert_allclose(com, 0.0, atol=1.0e-12)

    def test_bond_lengths_and_angle(self):
        o, h1, h2 = self.water.coords
        self.assertAlmostEqual(np.linalg.norm(h1 - o), 1.8121, places=12)
        self.assertAlmostEqual(np.linalg.norm(h2 - o), 1.8121, places=12)
        cos_angle = np.dot(h1 - o, h2 - o) / 1.8121**2
        self.assertAlmostEqual(
            np.degrees(np.arccos(cos_angle)), 104.368, places=10
        )

    def test_orientation(self):
        o, h1, h2 = self.water.coords
        np.testing.assert_allclose(self.water.coords[:, 0], 0.0, atol=1e-15)
        self.assertLess(o[2], 0.0)
        self.assertLess(h1[1], 0.0)
        self.assertGreater(h2[1], 0.0)

    def test_principal_axis_frame(self):
        r = self.water.coords
        m = self.water.masses
        r2 = np.einsum("n,ni,ni->", m, r, r)
        inertia = r2 * np.eye(3) - np.einsum("n,ni,nj->ij", m, r, r)
        off_diagonal = inertia - np.diag(np.diag(inertia))
        np.testing.assert_allclose(off_diagonal, 0.0, atol=1.0e-10)

    def test_rotational_constants(self):
        # Checks the masses and geometry against the paper's rotational
        # constants (x, y, z = a, b, c), which are given in cm^-1.
        r = self.water.coords
        m = self.water.masses
        r2 = np.einsum("n,ni,ni->", m, r, r)
        inertia = r2 * np.eye(3) - np.einsum("n,ni,nj->ij", m, r, r)
        constants = HARTREE_TO_CM / (2.0 * np.diag(inertia))
        np.testing.assert_allclose(constants, [9.49, 27.24, 14.57], atol=0.01)


class TestWaterNormalModes(unittest.TestCase):
    def setUp(self):
        self.L = water_normal_modes()
        self.water = water_molecule()

    def test_shape(self):
        self.assertEqual(self.L.shape, (3, 3, 3))

    def test_orthonormal(self):
        L = self.L.reshape(9, 3)
        np.testing.assert_allclose(L.T @ L, np.eye(3), atol=L_TOL)

    def test_eckart_translation(self):
        sqrt_m = np.sqrt(self.water.masses)
        translation = np.einsum("n,nak->ak", sqrt_m, self.L)
        self.assertLess(np.abs(translation).max(), L_TOL * sqrt_m.max())

    def test_eckart_rotation(self):
        sqrt_m = np.sqrt(self.water.masses)
        weighted = sqrt_m[:, None, None] * self.L
        rotation = np.cross(self.water.coords[:, :, None], weighted, axis=1)
        rotation = rotation.sum(axis=0)
        scale = np.abs(sqrt_m[:, None] * self.water.coords).max()
        self.assertLess(np.abs(rotation).max(), L_TOL * scale)


class TestWaterHarmonicFrequencies(unittest.TestCase):
    def test_frequencies(self):
        frequencies = water_harmonic_frequencies() * HARTREE_TO_CM
        np.testing.assert_allclose(frequencies, [3830.87, 1649.74, 3940.46])


class TestHarmonicAnalysis(unittest.TestCase):
    def test_water(self):
        normal_modes, frequencies = harmonic_analysis(water_molecule())
        np.testing.assert_array_equal(normal_modes, water_normal_modes())
        np.testing.assert_array_equal(
            frequencies, water_harmonic_frequencies()
        )

    def test_not_water(self):
        water = water_molecule()
        mol = Molecule(("N", "H", "H"), water.coords, water.masses)
        with self.assertRaises(ValueError):
            harmonic_analysis(mol)


class TestIsWater(unittest.TestCase):
    def setUp(self):
        self.water = water_molecule()

    def test_water(self):
        self.assertTrue(is_water(self.water))

    def test_wrong_symbols(self):
        mol = Molecule(("O", "H", "D"), self.water.coords, self.water.masses)
        self.assertFalse(is_water(mol))

    def test_wrong_masses(self):
        masses = self.water.masses.copy()
        masses[2] *= 2.0
        mol = Molecule(self.water.symbols, self.water.coords, masses)
        self.assertFalse(is_water(mol))

    def test_perturbed_geometry(self):
        coords = self.water.coords.copy()
        coords[1, 1] += 1.0e-3
        mol = Molecule(self.water.symbols, coords, self.water.masses)
        self.assertFalse(is_water(mol))

    def test_within_tolerance(self):
        coords = self.water.coords + 1.0e-8
        mol = Molecule(self.water.symbols, coords, self.water.masses)
        self.assertTrue(is_water(mol))

    def test_wrong_number_of_atoms(self):
        mol = Molecule(("H", "H"), np.zeros((2, 3)), np.ones(2))
        self.assertFalse(is_water(mol))
