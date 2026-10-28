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

.. _classes:

##############
Classes Design
##############

#. ``Molecule`` class.

   - Motivation: Describe the user's input.
   - State is the atomic symbols, XYZ coordinates, and isotopic masses.
   - Eventually replaced by SimDE's ``ChemicalSystem`` class.

#. ``Hamiltonian`` class.

   - Container for all the Hamiltonian terms.
   - Will eventually get replaced by SimDE's ``Hamiltonian`` class.

#. ``HamiltonianTerm`` class.

   - Base class for all Hamiltonian terms.
   - Will eventually get replaced by a SimDE class.

#. ``HarmonicHamiltonian`` class.

   - |omegavec| which is |Nvib| element vector 

#. ``AnharmonicPotential`` class.

   - cubic force constants (|Nvib| by |Nvib| by |Nvib|) and 
   - quartic force constants (|Nvib| by |Nvib| by |Nvib| by |Nvib|)

#. ``RigidRotorHamiltonian`` class.

   - Rotational constants

#. ``VibrationalCoriolisHamiltonian`` class.

   - |zetavec| which is |Nvib| by 3 by |Nvib|
   - |etavec| list

#. ``RovibrationalHamiltonian`` class.

    - |zetavec| which is |Nvib| by 3 by |Nvib|
    - |etavec| list

#. ``WatsonHamiltonian`` class.

    - |etavec| list.
    - Technically only need diagonal elements of |etavec|.

#. ``RovibrationalBasisFunction`` class.

   - Product of |Nvib| harmonic oscillator wavefunctions and an angular 
     momentum eigenstate.
   - Class state is a vector of vibrational quantum numbers, J, and K.

#. ``RovibrationalAnsatz`` class.

   - State is conceptually a vector of ``RovibrationalBasisFunction`` objects.
   - In practice, basis functions are generated on-the-fly.