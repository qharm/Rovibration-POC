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

.. _angular_wavefunction:

####################
Angular Wavefunction
####################

The point of this page is to derive relationships involving the angular
component of the rovibrational wavefunction, :math:`\ket{J,K,M}`.

:math:`\ket{J,K,M}` is an eigenstate of :math:`\jop^2`, and the molecule-fixed
z-component of the angular momentum operator, :math:`\jop_z` with eigenvalues:

.. math::

   \jop^2\ket{J,K,M} =& \hbar^2 J(J+1)\ket{J,K,M}\\
   \jop_z\ket{J,K,M} =& \hbar K\ket{J,K,M}.

We additionally note that :math:`\jop^2 = \jop_x^2 + \jop_y^2 + \jop_z^2`.

The molecule-fixed components of |jvec| obey the anomalous commutation
relation :math:`\left[\jop_\alpha, \jop_\beta\right] = -i\hbar \jop_\gamma`. As
a result :math:`\jop_\pm = \jop_x \pm i\jop_y` lower and raise :math:`K`
respectively. Note that the sign of the commutator and the function of the
ladder operators are the **reverse** of the space-fixed case.

The action of the operators on the angular momentum eigenstates are given by:
.. math::

   \jop^2\ket{J,K} =& \hbar^2 J(J+1)\ket{J,K}\\
   \jop_z\ket{J,K} =& \hbar K\ket{J,K}\\
   \jop_\pm\ket{J,K} =& \hbar C^\mp_{J,K}\ket{J,K\mp1}

where:

.. math::

   C^\pm_{J,K} =& \sqrt{J(J+1) - K(K\pm1)}
