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

.. _vibration_coriolis:

################################
Vibrational-Coriolis Hamiltonian
################################

The vibrational Coriolis Hamiltonian is purely vibrational. Writing out the
matrix products and inserting the expansions above:

.. math::

   \HvibCor =& \frac{1}{2}\sum_{n=0}^{2}\sum_{\alpha,\beta=1}^{3}
      \piop_\alpha\muop^{(n)}_{\alpha\beta}\piop_\beta\\
     =& \frac{1}{2}\sum_{n=0}^{2}\sum_{\alpha,\beta=1}^{3}
      \sum_{i,j,k,\ell=1}^{\Nvib}\sum_{\midx}
      \zeta^\alpha_{i,j}\zeta^\beta_{k,\ell}
      \eta^{(n)}_{\alpha\beta,\midx}
      \qop_i \pop_j \qop^{\midx} \qop_k \pop_\ell 

The operator ordering must be kept, because :math:`\pop_j` does not commute
with :math:`\qop^{\midx}\qop_k`. The :math:`n=0` terms reduce to

.. math::

   \frac{1}{2}\sum_{\alpha=1}^{3}
      \frac{\piop_\alpha^2}{I^0_{\alpha\alpha}} -
   \frac{\hbar^2}{8}\sum_{\alpha=1}^{3}
      \frac{1}{I^0_{\alpha\alpha}}

the second of which is a constant energy shift.

|HvibCor| is purely vibrational, so it is diagonal in :math:`J` and
:math:`K`:

.. math::

   \braket{\vvec',J',K'\mid \HvibCor\mid\vvec,J,K} =&
      \delta_{J',J}\delta_{K',K}\left[
      \frac{1}{2}\sum_{n=0}^{2}\sum_{\alpha,\beta=1}^{3}
      \sum_{i,j,k,\ell=1}^{\Nvib}\sum_{\midx}
      \zeta^\alpha_{i,j}\zeta^\beta_{k,\ell}
      \eta^{(n)}_{\alpha\beta,\midx}
      \braket{\vvec'\mid \qop_i \pop_j \qop^{\midx} \qop_k \pop_\ell
              \mid\vvec}\right.\\
      &\left.- \frac{\hbar^2}{8}\sum_{n=0}^{4}\sum_{\alpha=1}^{3}
      \sum_{\midx}\eta^{(n)}_{\alpha\alpha,\midx}
      \braket{\vvec'\mid \qop^{\midx}\mid\vvec}\right]

The :math:`n`\ th Coriolis term is a product of :math:`n+4` vibrational
operators and the :math:`n`\ th Watson term is a product of :math:`n`, which
sets the range of :math:`\vvec'` each term couples to.