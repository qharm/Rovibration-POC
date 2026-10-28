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

.. _rigid_rotor:

#######################
Rigid Rotor Hamiltonian
#######################

The rigid rotor Hamiltonian is then defined as:

.. math::

   \HRR \equiv& \frac{1}{2}\jvec^T\muvec^{(0)}\jvec\\
               =& \frac{1}{2}\sum_{\alpha,\beta=1}^{3}
                  \jop_\alpha\muop^{(0)}_{\alpha,\beta}\jop_\beta\\
               =& \frac{1}{2}\sum_{\alpha=1}^{3}
                  \frac{\jop_\alpha^2}{I^0_{\alpha,\alpha}}

fixed frame. Using :math:`\gamma_x = \frac{\hbar^2}{2I^0_{xx}}`,
:math:`\gamma_y = \frac{\hbar^2}{2I^0_{yy}}`, and
:math:`\gamma_z = \frac{\hbar^2}{2I^0_{zz}}`, we can rewrite:

.. math::
   \HRR =&\frac{1}{\hbar^2}\left[\gamma_x \jop_x^2 + \gamma_y \jop_y^2 +
                                  \gamma_z \jop_z^2\right]\\
          =&\frac{1}{\hbar^2}\left[
            \frac{\gamma_x+\gamma_y}{2}\left(\jop_x^2 + \jop_y^2\right) +
            \frac{\gamma_x-\gamma_y}{2}\left(\jop_x^2 - \jop_y^2\right) +
            \gamma_z \jop_z^2\right]\\
          =&\frac{1}{\hbar^2}\left[
            \frac{\gamma_x+\gamma_y}{2}\left(\jop^2 - \jop_z^2\right) +
            \frac{\gamma_x-\gamma_y}{4}\left(\jop_+^2 + \jop_-^2\right) +
            \gamma_z \jop_z^2\right]

where we used the fact that :math:`\muvec^{(0)}` is diagonal in the molecule-

|HRR| is not diagonal in this basis, because the rigid rotor
eigenstates are linear combinations of the :math:`\ket{J,K}` states. The matrix
elements of |HRR| in this basis are given by:

.. math::

    \braket{\vvec',J',K'\mid \HRR\mid\vvec,J,K} =&
       \delta_{\vvec',\vvec}\delta_{J',J}\left[
       \left(\frac{\gamma_x+\gamma_y}{2}\left(J(J+1) - K^2\right) +
             \gamma_z K^2\right)\delta_{K',K}\right.\\
       &\left.+ \frac{\gamma_x-\gamma_y}{4}\left(
          C^-_{J,K}C^-_{J,K-1}\,\delta_{K',K-2} +
          C^+_{J,K}C^+_{J,K+1}\,\delta_{K',K+2}\right)\right]

where the :math:`\hbar^2` from the two ladder operators cancels the
:math:`1/\hbar^2` in |HRR|.