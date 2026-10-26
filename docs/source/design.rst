#############
Module Design
#############

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

#. First pass will build entire Hamiltonian matrix. This is not scalable for
   large molecules so the design should have a mechanism for selecting subsets
   of the Hamiltonian matrix to build.

*****************
Classes Necessary
*****************

#. ``Molecule`` class.
   - State is the XYZ coordinates and masses.
   - Eventually replaced by SimDE's ``ChemicalSystem`` class.

#.  ``Qk`` and ``Pk`` classes.
   - Operators for normal modes and their conjugate momenta.
   - State is the mode index ``k``.

#. ``RovibrationalBasisFunction`` class.
   - Product of |Nvib| harmonic oscillator wavefunctions and an angular momentum
     eigenstate.
   - Class state is a vector of vibrational quantum numbers, J, and K.

#. ``RovibrationalAnsatz`` class.
   - State is a vector of ``RovibrationalBasisFunction`` objects.

#. ``TaylorExpansion`` class.
   - Represents the Taylor expansion of the PES about the minimum energy 
     geometry.
   - Minimum energy geometry is set to zero.  
   - Gradient is zero for minimum energy geometry.
   - Only harmonic and higher order terms are stored.
   - State:
     - |Lvec| which is |Natom| by 3 by |Nvib|
     - |omegavec| which is |Nvib| element vector 
     - cubic force constants (|Nvib| by |Nvib| by |Nvib|) and 
     - quartic force constants (|Nvib| by |Nvib| by |Nvib| by |Nvib|)

*******************
Driver Control Flow
*******************

#. Input is a minimum energy ``Molecule`` object.
   - Hard-code to water in the system test.
   - Assuming minimum energy allows focusing on Taylor expanding the PES.

#. Fill ``TaylorExpansion`` object.
   - Hard-code for water. Throw if input ``Molecule`` is not water.
   - A general approach needs second, third, and fourth derivatives of the
     electronic energy with respect to displacements from |rzvec|.

#. Compute |Izvec| and rotational parameters.
   - Result is a three-element vector of rotational constants.
   - Hard-code for water. Throw if input ``Molecule`` is not water.
   - A general approach would form the inertia tensor and diagonalize it.

#. Compute |zvec|.
    - |Nvib| by 3 by |Nvib|
    - Involves cross-products of |Lvec|.

#. Compute |bvec|.
   - |Nvib| by 3 by 3
   - Compute |avec|.
     - |Nvib| by 3 by 3
     - Needs |rzvec|, |Mvec|, and |Lvec|.
   - Transform |avec| to |bvec| using |Izvec|.

#. Build ``RovibrationaAnsatz`` object.
   - Input is the maximum vibrational quantum number and total angular momentum.

#. Build the Hamiltonian in the ``RovibrationalBasisFunction`` basis.
   - Loop over all pairs of basis functions and compute the matrix element.
     - One function per term (RR, HO, vibCor, rovib) in the Hamiltonian.
     - Functions needed to build the matrix elements are described in the
       next section.

#. Diagonalize the Hamiltonian matrix to get eigenvalues and eigenvectors.

************************
Building Matrix Elements
************************

We'll wrap the logic for building the matrix elements appearing in the
Hamiltonian terms. Each function will take the parameters needed to compute
a single element of the matrix element and will return a scalar. We need
functions to compute:

- Ladder coefficients
- Elements of eta
- One angular momentum operator.
- Two angular momenta operators.
- A single vibrational mode operator.
- The vibrational Coriolis term (QPQQP term)
- The rovibrational Coriolis term (QQP + QPQ term)
- An element of the rigid rotor Hamiltonian
- An element of the harmonic oscillator Hamiltonian
- An element of the vibrational Coriolis Hamiltonian
- An element of the rovibrational Coriolis Hamiltonian