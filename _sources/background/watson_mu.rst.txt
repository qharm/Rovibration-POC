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

.. _watson_mu_tensor:

################
Watson mu Tensor
################

The problem description notes that the Watson mu tensor, |muvec|, is expanded as
a fourth-order Taylor series about the minimum energy geometry,
:math:`\muvec \approx \sum_{n=0}^{4}\muvec^{(n)}`. This page works out the
details.

*******************
General Formulation
*******************

The :math:`n`-th term in the Taylor expansion of |muvec| is given by:

.. math::

   \muvec^{(n)} =& c_n \Izmhalf \bvec^n \Izmhalf

with:

- :math:`c_n = (n+1)\left(-\frac{1}{2}\right)^n`
- :math:`\bvec =& \sum_{k=1}^{\Nvib} \Izmhalf \avec_k \Izmhalf \qop_k`

and :math:`\avec_k` is a 3 by 3 matrix with elements:

.. math::

    a_{k,\alpha\beta} =&
       2 \sum_{\ell=1}^{\Natom} \sum_{m,n=1}^{3}
       \left(\delta_{\alpha\beta}\delta_{mn} -
             \delta_{\alpha n}\delta_{\beta m}\right)
       \sqrt{M_\ell}\,
       r_{\ell,m}^0
       L_{\ell,k,n}

***************
0-th Order Term
***************

The 0-th order term is straightforward:

.. math::

   \muvec^{(0)} =& c_0 \Izmhalf \bvec^0 \Izmhalf\\
                 =& 1 \cdot \Izmhalf \Izmhalf\\
                 =& \Ivec
