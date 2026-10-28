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

.. _rovibrational_coupling:

##################################
Rovibrational Coupling Hamiltonian
##################################

For the rovibrational coupling Hamiltonian, note that the molecule-fixed
components of |jvec| act only on the rotational coordinates, so they
commute with |qvec|, |pvec|, and therefore
:math:`\muvec^{(n)}` and |pivec|. Combined with the symmetry of
:math:`\muvec^{(n)}`, this lets the two Coriolis terms be combined:

.. math::

   \pivec^T\muvec^{(n)}\jvec =&
      \sum_{\alpha,\beta=1}^{3}
      \piop_\alpha\muop^{(n)}_{\alpha\beta}\jop_\beta
    = \sum_{\alpha,\beta=1}^{3}
      \jop_\alpha\piop_\beta\muop^{(n)}_{\alpha\beta}\\
   \jvec^T\muvec^{(n)}\pivec + \pivec^T\muvec^{(n)}\jvec =&
      \sum_{\alpha,\beta=1}^{3}\jop_\alpha\left(
        \muop^{(n)}_{\alpha\beta}\piop_\beta +
        \piop_\beta\muop^{(n)}_{\alpha\beta}\right)

where the first line relabeled :math:`\alpha\leftrightarrow\beta`. Inserting
this and the expansions above:

.. math::

   \Hrovib =& \frac{1}{2}\sum_{n=1}^{4}\sum_{\alpha,\beta=1}^{3}
      \muop^{(n)}_{\alpha\beta}\jop_\alpha \jop_\beta -
      \frac{1}{2}\sum_{n=0}^{3}\sum_{\alpha,\beta=1}^{3}\jop_\alpha\left(
        \muop^{(n)}_{\alpha\beta}\piop_\beta +
        \piop_\beta\muop^{(n)}_{\alpha\beta}\right)\\
     =& \frac{1}{2}\sum_{n=1}^{4}\sum_{\alpha,\beta=1}^{3}\sum_{\midx}
      \eta^{(n)}_{\alpha\beta,\midx}\qop^{\midx}\jop_\alpha \jop_\beta\\
      &- \frac{1}{2}\sum_{n=0}^{3}\sum_{\alpha,\beta=1}^{3}
      \sum_{k,\ell=1}^{\Nvib}\sum_{\midx}
      \zeta^\beta_{k,\ell}\eta^{(n)}_{\alpha\beta,\midx}
      \left(\qop^{\midx}\qop_k \pop_\ell + \qop_k \pop_\ell \qop^{\midx}\right)
      \jop_\alpha

Because :math:`\muvec^{(0)}` is a diagonal constant, the :math:`n=0`
Coriolis term reduces to

.. math::

   -\sum_{\alpha=1}^{3}
      \frac{\jop_\alpha\piop_\alpha}{I^0_{\alpha\alpha}}

|Hrovib| couples vibration and rotation, and its matrix elements
factor into a vibrational and a rotational part:

.. math::

   \braket{\vvec',J',K'\mid \Hrovib\mid\vvec,J,K} =&
      \delta_{J',J}\left[
      \frac{1}{2}\sum_{n=1}^{4}\sum_{\alpha,\beta=1}^{3}\sum_{\midx}
      \eta^{(n)}_{\alpha\beta,\midx}
      \braket{\vvec'\mid \qop^{\midx}\mid\vvec}
      \braket{J,K'\mid \jop_\alpha \jop_\beta\mid J,K}\right.\\
      &- \frac{1}{2}\sum_{n=0}^{3}\sum_{\alpha,\beta=1}^{3}
      \sum_{k,\ell=1}^{\Nvib}\sum_{\midx}
      \zeta^\beta_{k,\ell}\eta^{(n)}_{\alpha\beta,\midx}\\
      &\left.\times\braket{\vvec'\mid
         \qop^{\midx}\qop_k \pop_\ell + \qop_k \pop_\ell \qop^{\midx}
         \mid\vvec}
      \braket{J,K'\mid \jop_\alpha\mid J,K}\right]

The first (centrifugal distortion) term couples :math:`K` to
:math:`K, K\pm1, K\pm2`. The second (Coriolis) term couples :math:`K` to
:math:`K, K\pm1`.
