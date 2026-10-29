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

#####################
Problem Specification
#####################

The goal of rovib is to model the rovibrational energy levels of a non-linear
molecule containing |Natom| atoms. The location of the :math:`i`\ th atom in the
molecule is given by the vector :math:`\rvec_i` and the location of the
molecule as a whole is given by the 3|Natom|-dimensional vector
:math:`\rvec = \left(\rvec_1, \rvec_2, \ldots, \rvec_\Natom\right)`. Similarly,
|Mvec| is the diagonal matrix of atomic masses such that :math:`M_{ii}` is the
mass of the :math:`i`\ th atom. We assume that the Born-Oppenheimer
approximation is valid. In this case, we can assume that the nuclei move on a
single potential energy surface (PES), :math:`\Vop\left(\rvec\right)`, defined
by the electronic Hamiltonian. We will also assume that we have located a
minimum energy geometry on the PES, :math:`\rzvec`.

Our original goal can be restated as solving the time-independent Schrödinger
equation for the motion of the nuclei on the PES.

*****************
Coordinate System
*****************

:math:`\rvec` will in general be a somewhat arbitrary space-fixed reference
frame. In this coordinate system the molecule can translate, rotate, and
vibrate. The rovibrational Hamiltonian is simplified by defining a molecular
reference frame that maximally decouples the :math:`\Nvib = (3\Natom - 6)`
internal degrees of freedom of the molecule (vibrations) from the 6 external
degrees of freedom (translations and rotations). In practice, rotational-
vibrational coupling is only completely eliminated if the PES is harmonic. The
conditions under which the internal and external degrees of freedom are
maximally decoupled are known as the Eckart conditions.

To eliminate the translational degrees of freedom, we require the molecular
reference frame to move with the molecule. This is accomplished by requiring
the center of mass of the molecule to be at the origin of the molecular
reference frame.

We assume we have a non-linear molecule with |Natom| atoms. This means that
the molecule has  vibrational degrees of freedom.
We will define a molecular reference frame shortly. In the molecular reference
frame, |rvec| is the vector containing the Cartesian coordinates of all |Natom|
atoms such that :math:`\rvec_i` is the position of the :math:`i`\ th atom in
the molecule.  |Ivec| is the inertia
tensor of the molecule with coordinates defined by |rvec|.

We assume that  is the potential energy surface
(PES) of the molecule and that |rzvec| is a minimum energy geometry on the PES.
The molecular reference frame is defined so that its origin is the center of
mass of |rzvec| and its axes are aligned with the principal axes of inertia of
|rzvec|. To satisfy the Eckart conditions, we will use displacements from
|rzvec| to define the geometry of the molecule instead of |rvec|.
The displacements are defined as :math:`\rhovec = \rvec - \rzvec`.

The orientation and location of the molecular reference frame with respect to
the space-fixed reference frame can be fixed by requiring the center of mass to
be at the origin of the space-fixed reference frame and the Hessian of the PES
at |rzvec| to be diagonal in mass-weighted Cartesian coordinates. |Lvec| is the
transformation
matrix that diagonalizes the Hessian and |Izvec| is the diagonal inertia tensor
of the molecule at |rzvec|. :math:`\qvec = \Lvec^T\Mvec^{1/2}\rhovec` are the
resulting mass-weighted normal coordinates of the molecule describing the
vibrational degrees of freedom. The corresponding vibrational momenta are
:math:`\pvec = -i\hbar\frac{\partial}{\partial \qvec}`. As the molecule vibrates
and rotations it will be assumed to have total angular momentum |jvec| and
vibrational angular momentum |pivec|.

Other quantities used in the rovibrational Hamiltonian are:

- |Ieffvec| is the effective inertia tensor of the molecule with
  coordinates defined by |rvec|; it is the difference between |Ivec| and the
  vibrational contribution to the inertia tensor.
- |muvec| is the inverse of the effective inertia tensor of the
  molecule with coordinates defined by |rzvec|.
- |zetavec| is the Coriolis coupling tensor.
- |gammavec| are the rotation constants of the molecule.
- |omegavec| are the harmonic frequencies of the vibrational
  degrees of freedom of the molecule.
- |Hvib| is the vibrational Hamiltonian of the molecule.

***********
Definitions
***********

|zetavec| is given by:

.. math::

    \zeta_{i,j}^\alpha = \sum_{k=1}^{\Natom}
        \left(\Lvec_{k,i}\times\Lvec_{k,j}\right)_\alpha

|Ivec| is given by:

.. math::

    \Iop_{\alpha,\beta} = \sum_{i=1}^{\Natom} M_i
        \left(\left(\rvec_i\cdot\rvec_i\right)\delta_{\alpha,\beta} -
              \rop_{i,\alpha}\rop_{i,\beta}\right)

|Izvec| is given by:

.. math::

    I^0_{\alpha,\alpha} = \frac{\hbar^2}{2\gamma_\alpha}

|Ieffvec| is given by:

.. math::
    \Ieffop_{\alpha,\beta} = \Iop_{\alpha,\beta} -
       \sum_{i,j,k}^{\Nvib} \zeta^{\alpha}_{i,k}\zeta^\beta_{j,k}
       \qop_i \qop_j

|pivec| is given by:

.. math::

   \piop_\alpha = \sum_{i,j}^{\Nvib} \zeta_{i,j}^\alpha \qop_i \pop_j

**************************************
Deriving the Rovibrational Hamiltonian
**************************************

|Hvib| is given by:

.. math::

    \Hvib =
       \frac{1}{2}\left(\jvec-\pivec\right)^T\muvec\left(\jvec-\pivec\right) +
       \frac{1}{2}\sum_{i=1}^{\Nvib}\pop_i^2 +
       \Vop(\qvec) -
       \frac{\hbar^2}{8}\tr(\muvec)

where :math:`\frac{\hbar^2}{8}\tr(\muvec)` is the Watson term for a
non-linear molecule and is a type of centrifugal distortion term.


:math:`\Vop(\qvec)` is then expanded as a Taylor series about the minimum
energy geometry (which means the zeroth and 1st order terms are zero).
Truncating at the quartic term gives:

.. math::

    \Vop(\qvec) \approx&  \Vharm + \Vanharm\\
    \Vharm \equiv& \frac{1}{2}\sum_{i}^{\Nvib}\omega_i^2 \qop_i^2 \\
    \Vanharm \equiv& \frac{1}{6}\sum_{i,j,k}^{\Nvib}g_{i,j,k}
                       \qop_i \qop_j \qop_k +
               \frac{1}{24}\sum_{i,j,k,\ell}^{\Nvib}F_{i,j,k,\ell}
                 \qop_i \qop_j \qop_k \qop_\ell

where we used the fact that the harmonic force constants are diagonal in
mass-weighted coordinates and defined the third- and fourth-order force
constants as :math:`g_{i,j,k}` and :math:`F_{i,j,k,\ell}` respectively.

For convenience and consistency, |muvec| is also expanded as a
fourth-order Taylor series about the minimum energy geometry:

.. math::

    \muvec \approx& \sum_{n=0}^{4}\muvec^{(n)}\\
    \muvec^{(n)} =& c_n \Izmhalf \bvec^n \Izmhalf\\
    c_n =& (n+1)\left(-\frac{1}{2}\right)^n\\
    \bvec =& \sum_{k=1}^{\Nvib} \Izmhalf \avec_k \Izmhalf \qop_k\\
    a_{k,\alpha\beta} =&
       2 \sum_{\ell=1}^{\Natom} \sum_{m,n=1}^{3}
       \left(\delta_{\alpha\beta}\delta_{mn} -
             \delta_{\alpha n}\delta_{\beta m}\right)
       \sqrt{M_\ell}\,
       r_{\ell,m}^0
       L_{\ell,k,n}

Inserting the Taylor expansions of :math:`\Vop(\qvec)` and |muvec| into the
vibrational Hamiltonian, and retaining terms up to fourth order in |qvec|
gives:

.. math::

       \Hvib =&
         \HRR
         + \HHO + \Vanharm  + \HvibCor + \Hrovib + \Hwatson +
         \mathcal{O}(\qvec^5)\\
       \HHO \equiv&  \frac{1}{2}\sum_{i=1}^{\Nvib}\pop_i^2 + \Vharm\\
       \HvibCor\equiv& \frac{1}{2}\sum_{n=0}^{2}
         \pivec^T\muvec^{(n)}\pivec\\
       \Hwatson \equiv& -\sum_{n=0}^{4}\frac{\hbar^2}{8}\tr(\muvec^{(n)})\\
       \Hrovib \equiv& \frac{1}{2}\sum_{n=1}^{4}
         \jvec^T\muvec^{(n)}\jvec -
         \frac{1}{2}\sum_{n=0}^{3}\jvec^T\muvec^{(n)}\pivec-
         \frac{1}{2}\sum_{n=0}^{3}\pivec^T\muvec^{(n)}\jvec\\
         \mathcal{O}(\qvec^5) =&
                \frac{1}{2}\sum_{n=3}^{4}\pivec^T\muvec^{(n)}\pivec -
                \frac{1}{2}\jvec^T\muvec^{(4)}\pivec -
                \frac{1}{2}\pivec^T\muvec^{(4)}\jvec\\


- See :ref:`rigid_rotor` for details on the rigid rotor Hamiltonian,
  :math:`\HRR`.
- See :ref:`harmonic_oscillator` for details on the harmonic oscillator
  Hamiltonian, :math:`\HHO`.
- See :ref:`vibration_coriolis` for details on the vibrational Coriolis
  Hamiltonian, :math:`\HvibCor`.
- See :ref:`watson_hamiltonian` for details on the Watson Hamiltonian,
  :math:`\Hwatson`.
- See :ref:`rovibrational_coupling` for details on the rovibrational coupling
  Hamiltonian, :math:`\Hrovib`.

TODO: Stuff below here needs a new home.

The remaining Hamiltonians are easier to work with once
:math:`\muvec^{(n)}` and |pivec| are written explicitly in terms
of |qvec| and |pvec|.

.. math::

   \Bvec_k \equiv \Izmhalf \avec_k \Izmhalf

so that :math:`\bvec = \sum_{k}\Bvec_k \qop_k`, and using the fact that
the :math:`\qop_k` commute with one another, :math:`\muvec^{(n)}` is a
polynomial of degree :math:`n` in |qvec|:

.. math::

   \muop^{(n)}_{\alpha\beta} =&
      \sum_{\midx}\eta^{(n)}_{\alpha\beta,\midx}\qop^{\midx}\\
   \eta^{(n)}_{\alpha\beta,\midx} \equiv& c_n\left[
      \Izmhalf
      \Bvec_{m_1}\Bvec_{m_2}\cdots\Bvec_{m_n}
      \Izmhalf\right]_{\alpha\beta}\\
   \qop^{\midx} \equiv& \qop_{m_1}\qop_{m_2}\cdots \qop_{m_n}

where :math:`\midx = \left(m_1, m_2, \ldots, m_n\right)` is a multi-index
and :math:`\sum_{\midx}` runs each :math:`m_t` from 1 to
|Nvib|. For :math:`n=0`, :math:`\qop^{\midx} = 1` and
:math:`\eta^{(0)}_{\alpha\beta} = \delta_{\alpha\beta}/I^0_{\alpha\alpha}`.
Two properties are used below:

- :math:`\avec_k` is symmetric (a consequence of the Eckart conditions),
  so :math:`\muvec^{(n)}` is symmetric.
- :math:`\zeta^\alpha_{i,j} = -\zeta^\alpha_{j,i}`, so
  :math:`\zeta^\alpha_{i,i} = 0` and each :math:`\piop_\alpha` is Hermitian.

***************
Matrix Elements
***************

Let :math:`\ket{v}` be an eigenfunction of |HHO| with eigenvalue
:math:`E_v=\hbar\omega\left(v + \frac{1}{2}\right)` and :math:`\ket{J,K}` be an eigenfunction of
:math:`\jop^2`
and :math:`\jop_z` with eigenvalues :math:`\hbar^2 J(J+1)` and
:math:`K\hbar` respectively (here :math:`K` is the projection of
|jvec| on to the molecule's z-axis). A basis for the rovibrational
Hamiltonian is then given by the direct
product of one :math:`\ket{v}` per vibrational degree of freedom and
:math:`\ket{J,K}`, i.e., :math:`\ket{\vvec,J,K}`.

The molecule-fixed components of |jvec| obey the anomalous commutation
relation :math:`\left[\jop_x, \jop_y\right] = -i\hbar \jop_z`. As a result
:math:`\jop_\pm = \jop_x \pm i\jop_y` lower and raise :math:`K` respectively
(the reverse of the space-fixed case):

.. math::

   \jop_\pm\ket{J,K} = \hbar C^\mp_{J,K}\ket{J,K\mp1},\qquad
   C^\pm_{J,K} \equiv \sqrt{J(J+1) - K(K\pm1)}




|HvibCor| and |Hrovib| are sums of products of
:math:`\qop_k`, :math:`\pop_k`, and :math:`\jop_\alpha`, so their matrix
elements follow from those of the individual operators.

In mass-weighted normal coordinates
:math:`\qop_k = \sqrt{\frac{\hbar}{2\omega_k}}
\left(\aop_k + \aop_k^\dagger\right)` and
:math:`\pop_k = i\sqrt{\frac{\hbar\omega_k}{2}}
\left(\aop_k^\dagger - \aop_k\right)`,
where :math:`\aop_k` and :math:`\aop_k^\dagger` are the ladder operators of
mode :math:`k`. The single-mode matrix elements are therefore:

.. math::

   \braket{v'_k\mid \qop_k\mid v_k} =& \sqrt{\frac{\hbar}{2\omega_k}}
      \left(\sqrt{v_k+1}\,\delta_{v'_k,v_k+1} +
            \sqrt{v_k}\,\delta_{v'_k,v_k-1}\right)\\
   \braket{v'_k\mid \pop_k\mid v_k} =& i\sqrt{\frac{\hbar\omega_k}{2}}
      \left(\sqrt{v_k+1}\,\delta_{v'_k,v_k+1} -
            \sqrt{v_k}\,\delta_{v'_k,v_k-1}\right)

Operators on different modes commute, so the matrix element of a product of
:math:`\qop` and :math:`\pop` operators factors into one matrix element per
mode:

.. math::

   \braket{\vvec'\mid \op{O}_1 \op{O}_2\cdots \op{O}_p\mid\vvec} =
      \prod_{k=1}^{\Nvib}\braket{v'_k\mid \op{O}^{(k)}\mid v_k}

where :math:`\op{O}^{(k)}` is the product, in the original order, of the
factors that act on mode :math:`k` (the identity if there are none).
Single-mode products are evaluated by inserting
:math:`\sum_{v''_k}\ket{v''_k}\bra{v''_k}` between factors. Each factor changes
one quantum number by :math:`\pm 1`, so a product of :math:`p` factors only
connects states with
:math:`\sum_k\left|v'_k - v_k\right| \le p` and
:math:`\sum_k\left(v'_k - v_k\right) \equiv p \pmod 2`.

The reversed action of :math:`\jop_\pm` does not change the magnitude of the
|HRR| matrix elements, which only involve :math:`\jop_+^2 + \jop_-^2`, but it
fixes the signs below. Using
:math:`\jop_x = \frac{1}{2}\left(\jop_+ + \jop_-\right)` and
:math:`\jop_y = \frac{1}{2i}\left(\jop_+ - \jop_-\right)`:

.. math::

   \braket{J',K'\mid \jop_x\mid J,K} =& \frac{\hbar}{2}\delta_{J',J}
      \left(C^+_{J,K}\,\delta_{K',K+1} + C^-_{J,K}\,\delta_{K',K-1}\right)\\
   \braket{J',K'\mid \jop_y\mid J,K} =& \frac{i\hbar}{2}\delta_{J',J}
      \left(C^+_{J,K}\,\delta_{K',K+1} - C^-_{J,K}\,\delta_{K',K-1}\right)\\
   \braket{J',K'\mid \jop_z\mid J,K} =& \hbar K\,\delta_{J',J}\delta_{K',K}

Every :math:`\jop_\alpha` conserves :math:`J`, so products are evaluated within
a single :math:`J` block:

.. math::

   \braket{J,K'\mid \jop_\alpha \jop_\beta\mid J,K} =
      \sum_{K''=K-1}^{K+1}\braket{J,K'\mid \jop_\alpha\mid J,K''}
                          \braket{J,K''\mid \jop_\beta\mid J,K}

which is nonzero only for :math:`\left|K'-K\right| \le 2`.
