
.. _matrix_elements:

##########################
Evaluating Matrix Elements
##########################

The bulk of the work of this package is evaluating the matrix elements of the
Hamiltonian. This page focuses on the details of how the matrix elements will
be evaluated.

**************
Considerations
**************

#. Hamiltonian is a sum of terms.

   - The design calls for one class per Hamiltonian term.
   - The matrix elements for the individual terms will be evaluated term by 
     term, and then summed to get the total Hamiltonian matrix element.

#. Basis functions are assumed to be products of harmonic oscillator
   wavefunctions and angular momentum eigenstates.

#. The Hamiltonian will be sparse in the basis set.

   - Whether it's sparse enough to use sparse matrix libraries is TBD.
   - Initial implementation will focus on filling in the non-zero elements of
     a dense matrix.
   - Need a mechanism to get non-zero indices.

****************
Helper Functions
****************

#. ``non_zero_indices(Bra, t, Ket)``. Takes the basis set for the bra, a
   Hamiltonian term, and the basis set for the ket. Returns a generator that
   runs over bra-ket pairs for which the corresponding matrix element of the 
   Hamiltonian term is non-zero.

   - The idea to take a bra and a ket is to allow for the using different basis
     sets or batches of basis functions.
   - One common API. Internally can branch based on the input types.

***********
Algorithm 1
***********

This algorithm assumes we are building the full Hamiltonian matrix in a dense
format.

#. Input is the Hamiltonian ``H`` and the basis set ``B``.

#. Loop over Hamiltonian terms, ``t``.

#. Call  ``non_zero_indices(B,t,B)`` to get a generator ``G``.

#. Loop over ``(b,k)`` pairs in ``G``.

#. Compute the matrix element ``H[b,k] = <b|t|k>``.
