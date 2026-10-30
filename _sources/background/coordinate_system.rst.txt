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

.. _coordinate_system:

#################
Coordinate System
#################

The instantaneous geometry of the molecule, :math:`\tens{r}`, and the minimum
energy geometry, :math:`\tens{r}^0`, will in general be expressed with respect
to a somewhat arbitrary space-fixed reference frame. Describing nuclear motion
in the space-fixed coordinate system will in general couple the molecule's
internal (vibrational) and external (translational and rotational) degrees of
freedom. In practice, it is not possible to completely decouple the internal and
external degrees of freedom (e.g., centrifugal distortion caused by rotation
will affect the vibrational motion). Our goal is to instead find a coordinate
system in which the internal and external degrees of freedom are as decoupled as
possible. This coordinate system is known as the Eckart frame and the Eckart
conditions provide the mechanism for decoupling the internal and external
degrees of freedom. The goal of this section is to re-express quantities in
terms of the Eckart frame.

The minimum energy geometry of the molecule is a fixed reference geometry that
does not change as the molecule moves. It will follow the Eckart frame as the
molecule moves. At each instant we can position the minimum energy geometry in
the Eckart frame however we want. We will choose to position it so that its
center of mass lies at the origin and the principal axes are aligned with the
Eckart frame's axes. Since the center of mass of the minimum energy geometry is
at the origin, we know:

.. math::

   \sum_{i=1}^{\Natom} M_{ii}\vec{f}^0_i = \vec{0}.

where the use of :math:`f` instead of :math:`r` is to indicate that the
components of the vector are expressed in the Eckart frame. Similarly, we know
that in the Eckart frame the inertia tensor of the minimum energy geometry,
:math:`\tens{I}^0`, is diagonal and the elements are given by:

.. math::

   I^0_{\alpha\beta} = \delta_{\alpha,\beta}\sum_{i=1}^{\Natom} M_{ii}
       \left[\left(\vec{f}^0_i\cdot\vec{f}^0_i\right) -
             f^0_{i,\alpha}f^0_{i,\beta}\right]

Now consider the instantaneous geometry of the molecule. In the Eckart reference
frame, the instantaneous center of mass, :math:`\vec{f}_{\text{COM}}`, of the
molecule is given by:

.. math::

   \vec{f}_{\text{COM}} = \frac{\sum_{i=1}^{\Natom} M_{ii}\vec{f}_i}
       {\sum_{i=1}^{\Natom} M_{ii}}

and the instantaneous inertia tensor of the molecule, :math:`\tens{I}`, is
given by:

.. math::

   I_{\alpha\beta} = \sum_{i=1}^{\Natom} M_{ii}
       \left[\left(\vec{f}_i\cdot\vec{f}_i\right)\delta_{\alpha,\beta} -
             f_{i,\alpha}f_{i,\beta}\right].

Since the PES will be Taylor expanded it is convenient to describe the geometry
of the molecule in terms of displacements from the minimum energy geometry.
Since both the instantaneous and minimum energy geometries move with the Eckart
frame, so too will the displacements. To that end, we define
:math:`\vec{\rho}_i` as the displacement of the :math:`i`\ th atom from
:math:`\vec{f}^0_i`.

The definition of the displacements are then:

.. math::

   \vec{\rho}_i = \vec{f}_i - \vec{f}^0_i

For convenience we also define mass weighted displacements,
:math:`\tens{P} = \Mvec^{1/2}\tens{\rho}`.

*****************
Eckart Conditions
*****************

The Eckart frame will be defined in terms of three vectors,
:math:`\hat{e}_x`, :math:`\hat{e}_y`, and :math:`\hat{e}_z`, termed the Eckart
vectors. Again, for convenience, we define mass-weighted counterparts
:math:`\hat{E}_x`, :math:`\hat{E}_y`, and :math:`\hat{E}_z`. Since the Eckart
frame moves with the molecule, the Eckart vectors will coincide with the
translational vectors of the molecule. Projecting the translations onto the
displacements gives:

.. math::

   \hat{E}_\alpha\cdot\tens{P} =& \sum_{i=1}^{\Natom} M_{ii}
       \left(\hat{e}_\alpha\cdot\vec{\rho}_i\right)\\
       =& \sum_{i=1}^{\Natom} M_{ii}\rho_{i\alpha}

To define internal degrees of freedom, we want the displacements to be
orthogonal to the Eckart vectors. This leads to the *translational Eckart
conditions*:

.. math::

   \sum_{i=1}^{\Natom} M_{ii}\vec{\rho}_i = \vec{0}

Inserting the definition of the displacements gives:

.. math::

   \sum_{i=1}^{\Natom} M_{ii}\left(\vec{f}_i - \vec{f}^0_i\right) =&
       \sum_{i=1}^{\Natom} M_{ii}\vec{f}_i -
       \sum_{i=1}^{\Natom} M_{ii}\vec{f}^0_i \\
       =& \sum_{i=1}^{\Natom} M_{ii}\vec{f}_i = \vec{0}

where we have used the fact that the minimum energy geometry is defined so that
its center of mass is at the origin. Therefore, the translational Eckart
conditions are satisfied if the instantaneous center of mass of the molecule is
also at the origin of the Eckart frame.

We can define a basis for the rotational motion of the molecule as
:math:`\vec{s}_\alpha = \hat{e}_\alpha \times \vec{f}^0`. Defining
:math:`\vec{S}_\alpha = \tens{M}^{1/2}\vec{s}_\alpha`, we can project the
rotations onto the displacements:

.. math::

   \vec{S}_\alpha\cdot\tens{P} =& \sum_{i=1}^{\Natom} M_{ii}
       \left(\hat{e}_\alpha\times\vec{f}^0_i\right)\cdot\vec{\rho}_i\\
       =& \sum_{i=1}^{\Natom} M_{ii}
          \hat{e}_\alpha\cdot\left(\vec{f}^0_i\times\vec{\rho}_i\right)\\
       =& \sum_{i=1}^{\Natom}
           M_{ii}\left(\vec{f}^0_i\times\vec{\rho}_i\right)_\alpha

Yielding the *rotational Eckart conditions*:

.. math::

   \sum_{i=1}^{\Natom} M_{ii}\vec{f}^0_i\times\vec{\rho}_i = \vec{0}

Inserting the definition of the displacements gives:

.. math::

   \vec{0} =& \sum_{i=1}^{\Natom}
      M_{ii}\vec{f}^0_i\times\left(\vec{f}_i - \vec{f}^0_i\right) \\
       =& \sum_{i=1}^{\Natom} M_{ii}\vec{f}^0_i\times\vec{f}_i -
          \sum_{i=1}^{\Natom} M_{ii}\vec{f}^0_i\times\vec{f}^0_i \\
       =& \sum_{i=1}^{\Natom} M_{ii}\vec{f}^0_i\times\vec{f}_i

In other words, the rotational Eckart conditions are satisfied when the
instantaneous geometry of the molecule has no net rotation with respect to the
minimum energy geometry.

*******************************************
Conversion from Space-Fixed to Eckart Frame
*******************************************

The Eckart conditions lead to a projection operator:

.. math::

   \tens{P}_{\text{Ext}} =
     \frac{1}{\sum_{i=1}^{\Natom} M_{ii}}
     \sum_{\alpha=1}^{3}\hat{E}_\alpha\hat{E}_\alpha^T +
     \sum_{\alpha=1}^{3}\frac{\vec{S}_\alpha\vec{S}_\alpha^T}
                             {I^0_{\alpha\alpha}}

If :math:`\bar{H}` is the mass-weighted Hessian of the PES with respect to
displacements from the minimum energy geometry, then the mass-weighted Hessian
in the Eckart frame is given by:

.. math::

   \bar{\tens{H}}_{\text{Eckart}} =
      \left(\tens{1} - \tens{P}_{\text{Ext}}\right)^T
      \bar{\tens{H}}\left(\tens{1} - \tens{P}_{\text{Ext}}\right)

Solving:

.. math::

   \bar{\tens{H}}_{\text{Eckart}}\bar{\tens{L}} =
   \bar{\tens{L}}\tens{\Omega}^2

gives a diagonal matrix :math:`\tens{\Omega}` containing the vibrational
frequencies and the mass-weighted vibrational modes, :math:`\bar{\tens{L}}`.
The transformation from the space-fixed frame to the Eckart frame is then
given by:

.. math::

   \tens{L} = \tens{M}\bar{\tens{L}}
