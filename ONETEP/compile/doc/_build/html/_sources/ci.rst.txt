==============================
Configuration Interaction (CI)
==============================

:Author: David Turban, University of Cambridge

Introduction
============

The aim of the CI functionality is to evaluate the electronic
Hamiltonian with respect to a set of reference configurations obtained
from constrained DFT (cDFT) and LR-TDDFT. This makes it possible to
obtain transition rates between different electronic states using
Fermi’s Golden Rule (e.g. rate of charge transfer at
donor/acceptor-interface in an organic solar cell). Also, CI can help to
construct eigenstates of a system in situations where the basic ground-
and excited state methods fail due to deficiencies of the approximate
exchange-correlation functional.

A classic example for such a case is the binding curve of
:math:`\text{H}_{2}^{+}`. The LDA functional gives an incorrect
dissociation limit with a binding energy that is significantly too
small. The reason is that there is a spurious self-interaction of the
single electron with itself, such that a delocalisation of the electron
over both atoms will always yield a lower energy, even at infinite
separation. Even more advanced hybrid functionals (like B3LYP) do not
solve this problem since they only include part of the exact
Hartree-Fock exchange which would cancel the self-interaction. In this
scenario CI offers a way forward. By choosing the two states with the
electron fully localised on either atom as references one ensures the
correct long-range limit. Finally, CI is used to evaluate the
Hamiltonian in the reference basis and obtain approximate eigenstates.
This gives a very good match with LDA at short range and also retains
the physical dissociation limit.

Theory
======

As a first step we assume that one intends to find the Hamiltonian
matrix element :math:`\langle\Psi_{B}|\hat{H}|\Psi_{A}\rangle` between
two cDFT states. In cDFT the constrained solutions are obtained as
ground states of the electronic Hamiltonian augmented with a
constraining potential that pushes charge (and/or spin) around the
system:

.. math:: (\hat{H}+\hat{V}_{c})|\Psi_{c}\rangle=(\hat{H}+V_{c}\hat{w}_{c})|\Psi_{c}\rangle=F|\Psi_{c}\rangle.

Here the potential is written as the product of its magnitude
:math:`V_{c}` and a weighting operator :math:`\hat{w}_{c}` which
specifically acts on the donor and acceptor regions of the system with
appropriate signs. In ONETEP the weighting operator is built from local
orbitals (like PAOs or NGWFs) that define the donor and acceptor
regions. The eigenvalue :math:`F` is the energy :math:`E` of the
constrained solution plus a correction due to the (unphysical)
constraining potential:

.. math:: F=\langle\Psi_{c}|\hat{H}+V_{c}\hat{w}_{c}|\Psi_{c}\rangle=E[\rho_{c}]+V_{c}\int d\mathbf{r}w_{c}(\mathbf{r})\rho_{c}(\mathbf{r})=E+V_{c}N_{c}.

The magnitude of the potential :math:`V_{c}` takes the role of a
Lagrange multiplier that is chosen such that the amount of displaced
charge/spin matches the population target :math:`N_{c}`. It should be
noted that cDFT states are generally not eigenstates of the electronic
Hamiltonian :math:`\hat{H}`.

Using the cDFT potentials we now obtain

.. math::

   \begin{aligned}
   \langle\Psi_{B}|\hat{H}|\Psi_{A}\rangle & = & \langle\Psi_{B}|\hat{H}+\hat{V}_{A}-\hat{V}_{A}|\Psi_{A}\rangle\nonumber \\
    & = & F_{A}\langle\Psi_{B}|\Psi_{A}\rangle-\langle\Psi_{B}|\hat{V}_{A}|\Psi_{A}\rangle,\end{aligned}

which reduces the problem to calculating overlaps of states and matrix
elements of the cDFT potentials. At this point the expression still
contain the full many-body wave functions :math:`\Psi`. To obtain a
practical computational scheme we replace them with Kohn-Sham (KS)
determinants :math:`\Phi` and approximate

.. math:: \langle\Psi_{B}|\Psi_{A}\rangle\approx\langle\Phi_{B}|\Phi_{A}\rangle,\;\;\langle\Psi_{B}|\hat{V}|\Psi_{A}\rangle\approx\langle\Phi_{B}|\hat{V}|\Phi_{A}\rangle.

We will assume that all orbitals are chosen to be real functions, and
not worry about complex conjugation in inner products. The standard
result for the overlap of two Slater determinants is given by

.. math::

   \begin{aligned}
   \langle\Phi_{B}|\Phi_{A}\rangle & = & \int d\mathbf{r}_{1}\ldots d\mathbf{r}_{N}\frac{1}{\sqrt{N!}}\left|\begin{array}{ccc}
   \psi_{1}^{B}(\mathbf{r}_{1}) & \cdots & \psi_{1}^{B}(\mathbf{r}_{N})\\
   \vdots &  & \vdots\\
   \psi_{N}^{B}(\mathbf{r}_{1}) & \cdots & \psi_{N}^{B}(\mathbf{r}_{N})
   \end{array}\right|\times\frac{1}{\sqrt{N!}}\left|\begin{array}{ccc}
   \psi_{1}^{A}(\mathbf{r}_{1}) & \cdots & \psi_{1}^{A}(\mathbf{r}_{N})\\
   \vdots &  & \vdots\\
   \psi_{N}^{A}(\mathbf{r}_{1}) & \cdots & \psi_{N}^{A}(\mathbf{r}_{N})
   \end{array}\right|\nonumber \\
    & = & \int d\mathbf{r}_{1}\ldots d\mathbf{r}_{N}\psi_{1}^{B}(\mathbf{r}_{1})\ldots\psi_{N}^{B}(\mathbf{r}_{N})\left|\begin{array}{ccc}
   \psi_{1}^{A}(\mathbf{r}_{1}) & \cdots & \psi_{1}^{A}(\mathbf{r}_{N})\\
   \vdots &  & \vdots\\
   \psi_{N}^{A}(\mathbf{r}_{1}) & \cdots & \psi_{N}^{A}(\mathbf{r}_{N})
   \end{array}\right|\nonumber \\
    & = & \left|\begin{array}{ccc}
   \langle\psi_{1}^{A}|\psi_{1}^{B}\rangle & \cdots & \langle\psi_{1}^{A}|\psi_{N}^{B}\rangle\\
   \vdots &  & \vdots\\
   \langle\psi_{N}^{A}|\psi_{1}^{B}\rangle & \cdots & \langle\psi_{N}^{A}|\psi_{N}^{B}\rangle
   \end{array}\right|\;\;\;=\;\;\;\det\left(S_{AB}\right),\end{aligned}

where the functions :math:`\psi` denote KS orbitals. The result is
simply the determinant of the matrix :math:`S_{AB}` of overlaps between
KS orbitals of states A and B. In a similar fashion we can evaluate
matrix elements of a potential operator:

.. math::

   \begin{aligned}
   \langle\Phi_{B}|\hat{V}|\Phi_{A}\rangle & = & \int d\mathbf{r}_{1}\ldots d\mathbf{r}_{N}\psi_{1}^{B}(\mathbf{r}_{1})\ldots\psi_{N}^{B}(\mathbf{r}_{N})\left[\sum_{i}V(\mathbf{r}_{i})\right]\times\left|\begin{array}{ccc}
   \psi_{1}^{A}(\mathbf{r}_{1}) & \cdots & \psi_{1}^{A}(\mathbf{r}_{N})\\
   \vdots &  & \vdots\\
   \psi_{N}^{A}(\mathbf{r}_{1}) & \cdots & \psi_{N}^{A}(\mathbf{r}_{N})
   \end{array}\right|\nonumber \\
    & = & \sum_{i}\left|\begin{array}{ccccc}
   \langle\psi_{1}^{A}|\psi_{1}^{B}\rangle & \cdots & \langle\psi_{1}^{A}|\hat{V}|\psi_{i}^{B}\rangle & \cdots & \langle\psi_{1}^{A}|\psi_{N}^{B}\rangle\\
   \vdots &  & \vdots &  & \vdots\\
   \langle\psi_{N}^{A}|\psi_{1}^{B}\rangle & \cdots & \langle\psi_{N}^{A}|\hat{V}|\psi_{i}^{B}\rangle & \cdots & \langle\psi_{N}^{A}|\psi_{N}^{B}\rangle
   \end{array}\right|\nonumber \\
    & = & \sum_{ij}\langle\psi_{j}^{A}|\hat{V}|\psi_{i}^{B}\rangle\cdot C_{j,i}\nonumber \\
    & = & \det\left(S_{AB}\right)\times\text{tr}\left(V_{AB}\cdot S_{AB}^{-1}\right).\end{aligned}

In the third line the determinant is expanded along the :math:`i`-th
column. :math:`C_{j,i}` denotes cofactors of :math:`S_{AB}` which are
sign-adapted determinants of the submatrices formed by deleting the
:math:`j`-th row and :math:`i`-th column of :math:`S_{AB}`. A well-known
theorem in linear algebra states that the matrix of cofactors of an
invertible matrix is equal to the transpose of the inverse of the matrix
times its determinant.

Next, we discuss how a constrained reference state can be coupled to an
excited state from LR-TDDFT. In LR-TDDFT the excited states are
represented as superpositions of single-particle excitations from an
occupied to an unoccupied orbital. This information is contained in the
response density matrix :math:`R_{jb}`. A particular non-zero entry
indicates that a transition from valence orbital :math:`\psi_{j}` to
conduction orbital :math:`\psi_{b}` contributes to the excited state. In
the following indices :math:`i,j,k,\ldots` will denote valence orbitals
and indices :math:`b,c,d,\ldots` conduction orbitals. A natural choice
for a DFT wave function of such an excitation that retains the response
density by construction is

.. math:: |\Phi\rangle=\sum_{jb}R_{jb}|\Phi_{j}^{b}\rangle.

:math:`|\Phi_{j}^{b}\rangle` denotes a Slater determinant constructed
from the valence orbitals, except for the single valence orbital
:math:`j` replaced with conduction orbital :math:`b`. For the following
we assume that state B was obtained as a LR-TDDFT excitation, and A is a
(constrained) ground state as before. For the overlap we calculate

.. math::

   \begin{aligned}
   \langle\Phi_{B}|\Phi_{A}\rangle & = & \int d\mathbf{r}_{1}\ldots d\mathbf{r}_{N}\sum_{jb}R_{jb}\psi_{1}^{B}(\mathbf{r}_{1})\ldots\psi_{b}^{B}(\mathbf{r}_{j})\ldots\psi_{N}^{B}(\mathbf{r}_{N})\times\left|\begin{array}{ccc}
   \psi_{1}^{A}(\mathbf{r}_{1}) & \cdots & \psi_{1}^{A}(\mathbf{r}_{N})\\
   \vdots &  & \vdots\\
   \psi_{N}^{A}(\mathbf{r}_{1}) & \cdots & \psi_{N}^{A}(\mathbf{r}_{N})
   \end{array}\right|\nonumber \\
    & = & \sum_{jb}R_{jb}\underset{\uparrow j}{\left|\begin{array}{ccccc}
   \langle\psi_{1}^{A}|\psi_{1}^{B}\rangle & \cdots & \langle\psi_{1}^{A}|\psi_{b}^{B}\rangle & \cdots & \langle\psi_{1}^{A}|\psi_{N}^{B}\rangle\\
   \vdots &  & \vdots &  & \vdots\\
   \langle\psi_{N}^{A}|\psi_{1}^{B}\rangle & \cdots & \langle\psi_{N}^{A}|\psi_{b}^{B}\rangle & \cdots & \langle\psi_{N}^{A}|\psi_{N}^{B}\rangle
   \end{array}\right|}\nonumber \\
    & = & \sum_{ijb}R_{jb}\langle\psi_{i}^{A}|\psi_{b}^{B}\rangle\cdot C_{i,j}\nonumber \\
    & = & \det\left(S_{AB}\right)\times\text{tr}\left(T_{AB}^{vc}\cdot R^{\top}\cdot S_{AB}^{-1}\right).\end{aligned}

The matrix :math:`T_{AB}^{vc}` represents the overlap of the valence
orbitals of state A with the conduction orbitals of state B. The
derivation of the overlap with a potential operator is a bit more
involved but proceeds along similar lines:

.. math::

   \begin{aligned}
   \langle\Phi_{B}|\hat{V}|\Phi_{A}\rangle & = & \int d\mathbf{r}_{1}\ldots d\mathbf{r}_{N}\sum_{jb}R_{jb}\psi_{1}^{B}(\mathbf{r}_{1})\ldots\psi_{b}^{B}(\mathbf{r}_{j})\ldots\psi_{N}^{B}(\mathbf{r}_{N})\left[\sum_{i}V(\mathbf{r}_{i})\right]\left|\begin{array}{ccc}
   \psi_{1}^{A}(\mathbf{r}_{1}) & \cdots & \psi_{1}^{A}(\mathbf{r}_{N})\\
   \vdots &  & \vdots\\
   \psi_{N}^{A}(\mathbf{r}_{1}) & \cdots & \psi_{N}^{A}(\mathbf{r}_{N})
   \end{array}\right|\nonumber \\
   \nonumber \\
    & = & \sum_{i\ne j}\sum_{b}R_{jb}\left|\begin{array}{ccccccc}
   \langle\psi_{1}^{A}|\psi_{1}^{B}\rangle & \cdots & \langle\psi_{1}^{A}|\hat{V}|\psi_{i}^{B}\rangle & \cdots & \langle\psi_{1}^{A}|\psi_{b}^{B}\rangle & \cdots & \langle\psi_{1}^{A}|\psi_{N}^{B}\rangle\\
   \vdots &  & \vdots &  & \vdots &  & \vdots\\
   \langle\psi_{N}^{A}|\psi_{1}^{B}\rangle & \cdots & \langle\psi_{N}^{A}|\hat{V}|\psi_{i}^{B}\rangle & \cdots & \langle\psi_{N}^{A}|\psi_{b}^{B}\rangle & \cdots & \langle\psi_{N}^{A}|\psi_{N}^{B}\rangle
   \end{array}\right|\nonumber \\
    &  & +\sum_{jb}R_{jb}\left|\begin{array}{ccccc}
   \langle\psi_{1}^{A}|\psi_{1}^{B}\rangle & \cdots & \langle\psi_{1}^{A}|\hat{V}|\psi_{b}^{B}\rangle & \cdots & \langle\psi_{1}^{A}|\psi_{N}^{B}\rangle\\
   \vdots &  & \vdots &  & \vdots\\
   \langle\psi_{N}^{A}|\psi_{1}^{B}\rangle & \cdots & \langle\psi_{N}^{A}|\hat{V}|\psi_{b}^{B}\rangle & \cdots & \langle\psi_{N}^{A}|\psi_{N}^{B}\rangle
   \end{array}\right|\nonumber \\
   \nonumber \\
    & = & \sum_{ijb}R_{jb}\sum_{kl}\langle\psi_{k}^{A}|\hat{V}|\psi_{i}^{B}\rangle\langle\psi_{l}^{A}|\psi_{b}^{B}\rangle\cdot\epsilon_{kl}\epsilon_{ij}C_{kl,ij}\;\;+\;\;\sum_{ijb}R_{jb}\langle\psi_{i}^{A}|\hat{V}|\psi_{b}^{B}\rangle\cdot C_{i,j}.\end{aligned}

In the first determinant two columns are distinct from the overlap
:math:`S_{AB}`, we therefore expand along both. This leads to an
expression including the second cofactors :math:`C_{kl,ij}`. It follows
from Jacobi’s theorem that

.. math::

   \begin{aligned}
   \epsilon_{kl}\epsilon_{ij}C_{kl,ij} & = & \det\left(S_{AB}\right)\times\left[\left(S_{AB}^{-1}\right)_{ik}\left(S_{AB}^{-1}\right)_{jl}-\left(S_{AB}^{-1}\right)_{il}\left(S_{AB}^{-1}\right)_{jk}\right].\end{aligned}

Putting everything together we finally obtain

.. math::

   \begin{aligned}
   \langle\Phi_{B}|\hat{V}|\Phi_{A}\rangle & = & \det\left(S_{AB}\right)\times\left[\text{tr}\left(V_{AB}\cdot S_{AB}^{-1}\right)\text{tr}\left(T_{AB}^{vc}\cdot R^{\top}\cdot S_{AB}^{-1}\right)-\text{tr}\left(V_{AB}\cdot S_{AB}^{-1}\cdot T_{AB}^{vc}\cdot R^{\top}\cdot S_{AB}^{-1}\right)\right]\nonumber \\
    &  & +\det\left(S_{AB}\right)\times\text{tr}\left(W_{AB}^{vc}\cdot R^{\top}\cdot S_{AB}^{-1}\right),\end{aligned}

where :math:`W_{AB}^{vc}` refers to matrix elements of :math:`\hat{V}`
between valence orbitals of state A and conduction orbitals of state B.

We note that in general the Hamiltonian matrix obtained in the way shown
is not symmetric due to the approximations inherent in the DFT
formalism. Hence, the Hamiltonian must be symmetrised before eigenstates
can be obtained.

Implementation
==============

The CI functionality is implemented in ``couplings_mod``. For each
reference state the density kernel and NGWFs are read from the
corresponding files. Additionally, the cDFT-potentials are read from
file for a cDFT reference state. For an excited state from LR-TDDFT,
conduction kernel, conduction NGWFs and the response kernel are read. A
set of orthonormal orbitals representing the valence space is obtained
from the NGWF representation by solving the eigenvalue problem

.. math:: \sum_{\beta\gamma}K^{\alpha\beta}S_{\beta\gamma}x^{\gamma}=n\cdot x^{\alpha},

and restricting to the occupied subspace. Here :math:`K^{\alpha\beta}`
is the valence density kernel and :math:`S_{\beta\gamma}` the overlap
matrix of valence NGWFs. Orthonormal conduction orbitals are obtained in
an equivalent manner. The actual CI calculations then proceeds in this
basis as outlined in the theory section. It should be noted that the
orbitals obtained in this way generally do not correspond to the KS
orbitals (they do not result from a diagonalisation of the Hamiltonian).
However, both are related through an orthogonal transformation. Hence,
the determinants are identical, therefore all results are unaffected by
this choice of basis.

The transformation to a orthonormal basis comes with an inherent
:math:`N^{3}` scaling of the method. The computational effort is
expected to be comparable with a properties calculation (which involves
a diagonalisation of the Hamiltonian).

Performing a calculation
========================

This section explains how to set up a CI calculation, and points out a
couple of important things to look out for.

-  First perform calculations for desired reference states. For each
   state the density kernel and NGWFs have to be written to files
   (``.dkn`` and ``.tightbox_ngwfs``). For cDFT reference states the
   potentials and projectors are required (``.cdft`` and
   ``.tightbox_hub_projs``). For LR-TDDFT states conduction kernel and
   NGWFs are required (``.dkn_cond`` and ``.tightbox_ngwfs_cond``), as
   well as the response density matrix.

-  It is currently required that all reference states use the same unit
   cell, grid, geometry and identical atomic species with the same
   number of NGWFs and the same pseudopotentials. **NOTE:** The current
   implementation is not compatible with PAW!

-  Now set up a new input file for the CI calculation. It is recommended
   to copy the input file of one of the cDFT reference calculations.
   This ensures that the setup of the CI run is consistent with the
   reference calculations (in particular with the correct projectors).
   If a LR-TDDFT reference state is used, also copy the conduction
   species block into the file. Set ``TASK`` to ``COUPLINGS``.

-  Add the block ``couplings_states``. This tells the CI calculation
   which reference states to use. Here is an example:Each line
   corresponds to one reference state. The first column is short
   identifier for the state (currently unused). The second column
   indicates whether it is a cDFT or LR-TDDFT state, the third column
   contains the root name (i.e. name of original input file without
   extension). For a LR-TDDFT state, the fourth column determines the
   index of the excitation to be used (set to 0 for cDFT states).
   Finally, the fifth column is the energy in Hartree. It should be made
   sure that all energies are referenced to the same zero point.

-  The output is written to matrix files (using ``dense_write``). The
   names of these files consist of the root name of the CI run with
   extensions ``_ci_ham`` and ``_ci_ham_sym`` for the CI Hamiltonian and
   its symmetrised version, respectively. Eigenvalues and -states of the
   (symmetrised) CI Hamiltonian are written to files with extensions
   ``_ci_eigvals`` and ``_ci_eigstates`` (column-wise). If
   ``output_detail : VERBOSE`` is chosen, the results are also written
   to standard output.

References
==========

-  Extracting electron transfer coupling elements from constrained
   density functional theory, Q. Wu and T. Van Voorhis, J. Chem. Phys.
   **125**, 164105 (2006)

-  Exciton/Charge-Transfer Electronic Couplings in Organic
   Semiconductors, S. Difley and T. Van Voorhis, JCTC **7**, 594 (2011)

-  Determinants and matrices, A.C. Aitken, University mathematical texts
   vol. 20, Oliver and Boyd (1958)
