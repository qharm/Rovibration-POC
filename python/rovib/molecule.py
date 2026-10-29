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

"""Defines the Molecule class."""

import numpy as np


class Molecule:
    """The user's description of a (non-linear) molecule.

    This class is a stand-in for SimDE's ``ChemicalSystem`` class. The
    coordinates and masses are stored as read-only arrays.

    :param symbols: The atomic symbol of each atom.
    :type symbols: tuple of str
    :param coords: Cartesian coordinates of each atom in bohr. The shape is
        (number of atoms, 3).
    :type coords: numpy.ndarray
    :param masses: Isotopic mass of each atom in electron masses. The shape is
        (number of atoms,).
    :type masses: numpy.ndarray
    :raises ValueError: If ``coords`` or ``masses`` does not have the shape
        implied by ``symbols``.
    """

    def __init__(self, symbols, coords, masses):
        self._symbols = tuple(symbols)
        self._coords = np.array(coords, dtype=float)
        self._masses = np.array(masses, dtype=float)

        natoms = len(self._symbols)
        if self._coords.shape != (natoms, 3):
            raise ValueError(
                f"Expected coordinates with shape ({natoms}, 3), got "
                f"{self._coords.shape}."
            )
        if self._masses.shape != (natoms,):
            raise ValueError(
                f"Expected masses with shape ({natoms},), got "
                f"{self._masses.shape}."
            )

        self._coords.setflags(write=False)
        self._masses.setflags(write=False)

    @property
    def symbols(self):
        """The atomic symbol of each atom.

        :rtype: tuple of str
        """
        return self._symbols

    @property
    def coords(self):
        """Cartesian coordinates of each atom in bohr.

        :rtype: numpy.ndarray
        """
        return self._coords

    @property
    def masses(self):
        """Isotopic mass of each atom in electron masses.

        :rtype: numpy.ndarray
        """
        return self._masses

    @property
    def natoms(self):
        """The number of atoms in the molecule.

        :rtype: int
        """
        return len(self._symbols)

    @property
    def nvib(self):
        """The number of vibrational degrees of freedom.

        The molecule is assumed to be non-linear, so this is
        :math:`3\\Natom - 6`.

        :rtype: int
        """
        return 3 * self.natoms - 6

    def __eq__(self, other):
        if not isinstance(other, Molecule):
            return NotImplemented
        return (
            self.symbols == other.symbols
            and np.array_equal(self.coords, other.coords)
            and np.array_equal(self.masses, other.masses)
        )
