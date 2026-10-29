.. Copyright 2026 QHARM Team
..
.. Licensed under the Apache License, Version 2.0 (the "License");
.. you may not use this file except in compliance with the License.
.. You may obtain a copy of the License at
..
.. http://www.apache.org/licenses/LICENSE-2.0
..
.. Unless required by applicable law or agreed to in writing, software
.. distributed under the License is distributed on an "AS IS" BASIS,
.. WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
.. See the License for the specific language governing permissions and
.. limitations under the License.

##############
Library Design
##############

**************
Considerations
**************

#. Should be extensible to other molecules. But will be hard-coded to water
   for now.

   - Main feature gap is getting the anharmonic force constants.

#. First pass is a standalone Python module. Eventually will want to go back 
   and modularize it for use with SimDE.

   - Envisioning a top-level driver where each subfunction will eventually be a
     PluginPlay module.

#. First pass will build entire Hamiltonian matrix.

   - This is not scalable. Eventually want to do Davidson diagonalization so 
     will need to be able to do sigma builds.

#. Will use atomic units for all calculations.

#. Need complex arithmetic because matrix elements of P and J are complex.

#. Object-oriented.

   - Classes discussed in :ref:`classes` section.

#. One class per Hamiltonian term.

   - Like SimDE, terms should carry state needed to evaluate matrix elements,
     but serve primarily as a DSL for determining what terms are in the
     Hamiltonian.

#. Need a class for describing the space in which the Hamiltonian will be
   diagonalized.

   - Carries parameters to denote how the basis is truncated.
   - Also DSL-like. Storing full basis is not scalable.

#. Driver responsible for going from a molecule to a Hamiltonian operator to
   the eigenvalues of the Hamiltonian operator.


*******************
Driver Control Flow
*******************

#. Input is a minimum energy ``Molecule`` object.

   - Hard-code to water in the system test.
   - Assuming minimum energy allows focusing on Taylor expanding the PES.

#. Do harmonic analysis to get |Lvec| and |omegavec|.

#. Fill ``HarmonicHamiltonian`` object.

   - Hard-code for water. Throw if input ``Molecule`` is not water.
   - A general approach needs second derivatives of the electronic energy with 
     respect to displacements from |rzvec|.

#. Fill ``AnharmonicPotential`` object.

   - Hard-code for water. Throw if input ``Molecule`` is not water.
   - A general approach needs third and fourth derivatives of the electronic
     energy with respect to displacements from |rzvec|.

#. Fill ``RigidRotorHamiltonian`` object.

   - Build the inertia tensor and diagonalize it.
   - Test against hard-coded values for water.

#. Build the eta tensors.

   #. Compute |avec|.

      - |Nvib| by 3 by 3
      - Needs |rzvec|, |Mvec|, and |Lvec|.

   #. Compute |Bvec| by transforming |avec| using |Izvec|.

      - |Nvib| by 3 by 3

   #. Build the five eta tensors from |Bvec| and |Izvec|.

      - n = 0 is 3 by 3
      - n = 1 is 3 by 3 by |Nvib|
      - n = 2 is 3 by 3 by |Nvib| by |Nvib|
      - n = 3 is 3 by 3 by |Nvib| by |Nvib| by |Nvib|
      - n = 4 is 3 by 3 by |Nvib| by |Nvib| by |Nvib| by |Nvib| 

#. Build ``WatsonHamiltonian`` object.

#. Compute |zetavec|.

   - |Nvib| by 3 by |Nvib|
   - Involves cross-products of |Lvec|.

#. Build ``RovibrationalHamiltonian`` and ``VibrationalCoriolisHamiltonian``

   - Input is |zetavec| and |etavec|.

#. Build ``RovibrationalAnsatz`` object.

   - Input is the maximum vibrational quantum number and total angular momentum.

#. Build the Hamiltonian in the ``RovibrationalBasisFunction`` basis.

   - Details are in :ref:`matrix_elements` section.

#. Diagonalize the Hamiltonian matrix to get eigenvalues/eigenvectors.

   - Driver returns the eigenvalues and eigenvectors
