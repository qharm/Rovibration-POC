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

.. _harmonic_oscillator:

###############################
Harmonic Oscillator Hamiltonian
###############################


|HHO| is diagonal in this basis with matrix elements:

.. math::

    \braket{\vvec',J',K'\mid\HHO\mid\vvec,J,K} =
       \delta_{\vvec',\vvec}\delta_{J',J}\delta_{K',K}
       \sum_{i=1}^{\Nvib}\hbar\omega_i\left(v_i + \frac{1}{2}\right)