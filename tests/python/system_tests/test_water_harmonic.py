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
"""Reproduces the "HO" row of Table I of doi:10.1063/5.0333947."""

import unittest

import numpy as np

from rovib import (
    RovibrationalAnsatz,
    dominant_basis_function,
    hartree_to_wavenumber,
    rovibrational_energies,
    water_molecule,
)

#: The paper's ground-state energy, E_000, in cm^-1. (RMR verified)
GROUND_STATE = 4710.5

#: The paper's excitation energies, E_v - E_000, in cm^-1. (RMR verified)
EXCITATIONS = {
    (0, 1, 0): 1649.7,
    (0, 2, 0): 3299.5,
    (1, 0, 0): 3830.9,
    (0, 0, 1): 3940.5,
}

#: The paper reports energies to 0.1 cm^-1.
TOL = 0.05


class TestWaterHarmonic(unittest.TestCase):

    def test_table1_harmonic(self):
        water = water_molecule()
        energies, eigenvectors = rovibrational_energies(water, vmax=3, J=0)
        energies = hartree_to_wavenumber(energies)
        ansatz = RovibrationalAnsatz(water.nvib, vmax=3, J=0)
        levels = {
            dominant_basis_function(eigenvectors[:, i], ansatz).v: energy
            for i, energy in enumerate(energies)
        }

        ground_state = levels[(0, 0, 0)]
        self.assertAlmostEqual(ground_state, energies[0])
        self.assertAlmostEqual(ground_state, GROUND_STATE, delta=TOL)
        for v, excitation in EXCITATIONS.items():
            with self.subTest(v=v):
                self.assertAlmostEqual(levels[v] - ground_state,
                                       excitation,
                                       delta=TOL)

    def test_rotational_degeneracy(self):
        J = 2
        energies, _ = rovibrational_energies(water_molecule(), vmax=2, J=J)
        self.assertEqual(energies.size, 3**3 * (2 * J + 1))
        multiplets = energies.reshape(-1, 2 * J + 1)
        np.testing.assert_allclose(multiplets, multiplets[:, :1] * np.ones(
            (1, 2 * J + 1)))
