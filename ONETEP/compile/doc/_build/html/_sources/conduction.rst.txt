=====================================================================
Conduction NGWF optimisation and optical absorption spectra
=====================================================================

:Author: Laura E. Ratcliff, Imperial College London
	 
:Date:   July 2011

Conduction calculations
=======================

As a consequence of the NGWF optimisation process in ONETEP the occupied
(valence) Kohn-Sham states are well represented by the NGWFs, but the
unoccupied (conduction) NGWFs are not, so that upon diagonalisation of
the Hamiltonian at the end of a calculation, if one were to compare the
resulting eigenvalues with a conventional cubic-scaling DFT code such as
CASTEP [Clark2005]_, the ONETEP conduction states would be
higher in energy than the CASTEP states, and some conduction states
might be missing [Skylaris2005]_. In order to correct this
problem, a method has been implemented whereby a second set of NGWFs
(referred to as the conduction NGWFs) are optimised to accurately
represent the Kohn-Sham conduction states. This is done with asymptotic
linear-scaling computational effort by constructing an idempotent
density matrix representing the manifold of conduction states (rather
than by solving for them explicitly). Optimisation of the NGWFs
describing the conduction states proceeds as with the valence states,
using a dual-loop system by which the conduction density kernel and the
conduction NGWF coefficients are simultaneously optimised.

It should be noted that the Kohn-Sham eigenvalues will of course not be
expected to exactly correspond to the true quasi-particle energies,
however in practice reasonable agreement with experiment has been seen
to occur in a number of systems, particularly when using the scissor
operator [Godby1986]_, [Gygi1989]_.

The conduction NGWF optimisation takes the form of a non-self-consistent
calculation following a ground-state calculation, where the density and
potential calculated in the ground-state calculation are re-used. A
projected Hamiltonian is then constructed in the conduction NGWF basis,
using the density operator as a projection operator. This projected
Hamiltonian is modified to avoid problems which might occur if the
Hamiltonian and density operators do not commute perfectly.
Additionally, the valence states are shifted up in energy by some amount
:math:`w`, such that they become higher in energy than the conduction
states. The projected conduction Hamiltonian is thus written:

.. math::

   \begin{aligned}
   \left(H_\chi^{\textrm{proj}}\right)_{\alpha\beta}&=&\langle \chi_\alpha|\hat{H}-\hat{\rho}\left(\hat{H}-w\right)\hat{\rho}|\chi_\beta\rangle\\ \nonumber
   &=&\left(H_\chi\right)_{\alpha\beta} -\left(T^\dag K H_\phi KT\right)_{\alpha\beta}\\ \nonumber
   &&+w\left(T^\dag K S_\phi KT\right)_{\alpha\beta}, \end{aligned}

where :math:`\{|\phi_{\alpha}\rangle\}` is the set of valence NGWFs and
:math:`\{|\chi_{\alpha}\rangle\}` the set of conduction NGWFs.
:math:`{\bm{\rho}}` is the valence density matrix, :math:`\mathbf{K}` is
the valence density kernel, :math:`\mathbf{S_{\phi}}` is the valence
overlap matrix and :math:`\mathbf{H_{\phi}}` is the valence Hamiltonian.
:math:`\mathbf{S_{\chi}}` is the conduction overlap matrix,
:math:`\mathbf{T}` is the valence-conduction cross overlap matrix
defined as
:math:`T_{\alpha\beta}=\langle \phi_{\alpha} | \chi_{\beta} \rangle`,
:math:`\mathbf{H_{\chi}}` is the (unprojected) conduction Hamiltonian,
:math:`\mathbf{H_\chi^{\textrm{proj}}}` is the projected conduction
Hamiltonian, :math:`\mathbf{Q}` is the conduction density matrix and
:math:`\mathbf{M}` is the conduction density kernel. The conduction
NGWFs and kernel are then minimised with respect to the energy
expression
:math:`E=\text{tr}\left[\mathbf{Q}\mathbf{H_\chi^{\textrm{proj}}}\right]`,
following the same procedure as in a standard ONETEP calculation. The
shift can either be set to a constant value, or updated during a
calculation, by setting it to be higher than the highest eigenvalue as
calculated in the conduction NGWF basis.

At the end of the conduction NGWF optimisation process, the valence and
conduction NGWF basis sets are combined into a new ‘joint’ basis, which
will be capable of accurately representing both the occupied and
unoccupied Kohn-Sham states. Other properties such as optical absorption
spectra can then be calculated in this joint basis.

For further information see Ratcliff *et
al*. [Ratcliff2011]_.

Performing conduction calculations in ONETEP
============================================

In order to optimise a set of NGWFs capable of accurately representing
the Kohn-Sham conduction states in ONETEP, it is first necessary to have
performed a standard ONETEP ground-state calculation and have retained
the density kernel and NGWF output files. No special parameter values
are required for this stage, although it may be worth setting
ODD\_PSINC\_GRID to true, as conduction NGWF radii generally need to be
larger than valence NGWF radii in order to achieve large convergence,
and so it is more likely that the FFT box will be required to be equal
to the psinc grid, and as both stages of the calculation must have the
same cut-off energy and therefore grid size, it is desirable to have an
odd grid for both the cell size and FFT box.

Once a ground-state calculation has been performed, a conduction
calculation can be performed by setting TASK=COND. The number of
conduction NGWFs per species and their radii must then be specified in
the SPECIES\_COND block, which follows the same pattern as the species
block. The initial NGWFs, if not specified, will be equal to the initial
NGWFs used for the valence density matrix. This choice can be overridden
by specifying different choices in a SPECIES\_ATOMIC\_SET\_COND block,
which can be set to use the pseudoatomic solver by setting “SOLVE” for
each species. One useful option is to specify that certain valence
states, particularly those known to be fully filled and thus not
expected to contribute to the manifold of conduction states, should be
included in the calculation of the ground state of the pseudoatom but
left out of the conduction NGWF set. For example, when generating
conduction NGWFs for Cadmium, one might want to include the 10 filled 4d
states in the atom calculation, but since they are not expected to
contribute to the unoccupied states, they can be omitted from the
conduction NGWFs by setting a splitnorm of “-1” for them, through the
following solver string: “SOLVE conf=4d10:-1”.

The conduction density kernel must contain a specific number of occupied
states, and only these states will contribute to the NGWF gradient:
COND\_NUM\_STATES is used to specify the number of conduction states to
be optimised. In principle this could be any number, but in practice the
higher energy conduction states converge rather slowly with respect to
conduction NGWF radii, and in particular completely delocalised
conduction states are very hard to represent using localised basis
functions. Therefore results should be treated with caution when
optimising high energy conduction states. A good rule-of-thumb is, if at
all possible, to try to choose a manifold of conduction states above
which there is a reasonable sized gap in the density-of-states (as
calculated with the valence NGWFs).

For truly asymptotically linear-scaling calculations, one must truncate
the conduction density kernel. The cutoff for this truncation is
specified using COND\_KERNEL\_CUTOFF, although it is expected that high
levels of kernel truncation will significantly limit the accuracy of the
calculated conduction states. If unsure, do not use kernel truncation
unless you are confident you have verified that the properties you are
interested in are unaffected by the truncation.

At the end of a conduction calculation, diagonalisations are
automatically performed of the valence Hamiltonian, both the projected
and unprojected conduction Hamiltonians and the joint valence-conduction
Hamiltonian. The eigenvalues are written to the corresponding .bands
files. However, no joint basis density kernel is generated and so the
occupancies are not calculated within this basis. The unprojected
conduction eigenvalues are of limited use to most users, as it is
difficult to determine which are conduction states and which are poorly
represented valence states. For the projected conduction eigenvalues,
the gap referred to in the output is not the usual gap, rather it is the
gap between the highest optimised conduction state and the lowest
unoptimised conduction state. If required, it is also possible to plot
the orbitals in either the valence and conduction NGWF basis sets,
and/or in the joint basis set, using the keywords
COND\_PLOT\_VC\_ORBITALS and COND\_PLOT\_JOINT\_ORBITALS.

A standalone properties calculation can also be performed on the basis
of sets of valence and conduction NGWFs and kernels which have already
been calculated. To enable this, set TASK=PROPERTIES\_COND: the options
COND\_READ\_TIGHTBOX\_NGWFS and COND\_READ\_DENSKERN will automatically
be enabled, the NGWF optimisation will skipped and the calculation will
proceed straight to the stage of diagonalisation of the Hamiltonian and
plotting of the orbitals.

Automatic setup of COND\_NUM\_STATES
------------------------------------

Since it is not always straightforward to guess a sensible number of
conduction states to converge, the code will by default attempt to
choose an appropriate number of states for the user. By default, at the
start of a conduction optimisation, the code will count the number of
unoccupied states of the valence Hamiltonian with negative eigenvalues
to arrive at a guess of the number of bound states in a finite system.
The code will also check for any degeneracy of the highest unoccupied
state included in the calculation and automatically include more states
until an energy gap of at least 0.001 Ha between states that get
optimised and states that are unoptimised is achieved.

Since counting the number of eigenstates with negative eigenvalues in
order to obtain the number of bound states is only strictly valid in
finite systems, it is possible for the user to define an energy range,
as measured from the HOMO energy level, and the code will attempt to
optimise all states within that energy range. The two keywords
controlling the automatic conduction state setup are given by
COND\_ENERGY\_RANGE and COND\_ENERGY\_GAP. The first keyword defines the
desired energy range in Hartree while COND\_ENERGY\_GAP defines the
minimum required energy gap between the highest conduction state that is
optimised and the lowest conduction state that stays unoptimised.

Setting the shift
-----------------

There are a number of parameters relating to the shift, :math:`w`, used
in the projected conduction Hamiltonian. It is possible to keep the
shift at some fixed value (defined using COND\_INIT\_SHIFT) during the
calculation, by setting COND\_FIXED\_SHIFT to true. Alternatively, it
can be automatically updated during the calculation, which is usually
the safest way to proceed. This is achieved by calculating the highest
eigenvalue within the conduction NGWF basis at the start of each NGWF
iteration (providing COND\_CALC\_MAX\_EIGEN is set to true), and
comparing the current shift to this eigenvalue. Providing the shift is
higher than the highest eigenvalue, it remains unchanged, but if the
maximum eigenvalue has become greater than the current shift, it is
updated to equal the maximum eigenvalue plus some extra buffer value
(defined by COND\_SHIFT\_BUFFER).

Local minima
------------

In practice, it is sometimes possible to become trapped in local minima,
where the ordering of states within the initial unoptimised basis
doesn’t correspond to the correct order, and so sometimes states are
missed. The problem can be identified by decreasing
NGWF\_THRESHOLD\_ORIG and seeing if the gradient stagnates while the
energy continues to decrease, or by plotting convergence graphs with
conduction NGWF radii where sharp changes in energy are sometimes
observed with small changes in conduction NGWF radii. In practice it is
therefore very important to systematically converge calculations with
respect to the conduction NGWF radii, which might require larger values
than ground-state ONETEP calculations. This problem can typically be
avoided by optimising some extra states (COND\_NUM\_EXTRA\_STATES) above
the required number of conduction states for a few iterations
(COND\_NUM\_EXTRA\_ITS) (typically 5-10 iterations). Selecting the
required number of extra states to include is mostly a trial and error
process whereby the number of extra states should be increased until no
changes are seen in the calculated conduction energy.

Additional notes on input parameters
------------------------------------

As the ground-state NGWFs and density kernel are required for the
conduction calculation, READ\_TIGHTBOX\_NGWFS and READ\_DENSKERN are
automatically set to true. There are separate variables for the
corresponding conduction quantities (COND\_READ\_TIGHTBOX\_NGWFS and
COND\_READ\_DENSKERN) which can be set to true for restarting conduction
calculations. The parameters WRITE\_TIGHTBOX\_NGWFS and WRITE\_DENSKERN
are not independently specified for the conduction and valence NGWF
basis sets.

Conduction calculations in implicit solvent
-------------------------------------------

| Some care has to be taken when performing a conduction optimisation
  for a system embedded in an implicit solvent (see the separate
  Implicit Solvation documentation on how to perform a ground state
  calculation in implicit solvent). The reason for this is that the
  ground state in the implicit solvation model is often computed in a
  two step process. In the first step the solvation cavity is computed
  as an isosurface of the ground state density in vacuum, while in the
  second step the ground state of the system is evaluated for that fixed
  cavity. In order for the conduction calculation to be consistent with
  the ground state calculation, the same solvation cavity has to be used
  in both cases. To ensure that the code uses the correct restart files
  when setting up the ground state Hamiltonian at the beginning of a
  conduction optimisation, the following sets of keywords should be
  used:

| ``Task : Singlepoint Cond``
| ``is_implicit_solvent : T``
| ``is_auto_solvation : T``
| ``is_smeared_ion_rep : T``

| The code will then automatically perform a SinglePoint and a
  conduction calculation in the implicit solvent, using the same
  solvation cavity for both the ground state and the conduction state
  calculation. This is achieved by writing .vacuum\_dkn and
  .vacuum\_tightbox\_ngwfs files that are used to set up the solvation
  cavity.

If further conduction calculations are required using the same ground
state, for example in order to change the number of conduction states
converged, it is possible to change the Task to COND and include the
keyword is\_separate\_restart\_files: T. This triggers the use of the
.vacuum files to set up the correct solvation cavity at the beginning of
the COND calculation.

Optical absorption spectra
==========================

The calculation of matrix elements for the generation of optical
absorption spectra using Fermi’s golden rule has been implemented in
ONETEP following the method used in CASTEP, as outlined by
Pickard [Pickard1997]_. Using the dipole approximation, the
imaginary component of the dielectric function is defined as

.. math::
   :label: imag_diel

   \varepsilon_2\left(\omega\right)=\frac{2e^2\pi}{\Omega\epsilon_0}\sum_{\mathbf{k},v,c}\left|\langle\psi_{\mathbf{k}}^{c}|\mathbf{\hat{q}}\cdot\mathbf{r}|\psi_{\mathbf{k}}^{v}\rangle\right|^2\delta\left(E_{\mathbf{k}}^{c}-E_{\mathbf{k}}^{v}-\hbar\omega\right) ,

where :math:`v` and :math:`c` denote valence and conduction bands
respectively, :math:`|\psi_{\mathbf{k}}^{n}\rangle` is the :math:`n`\ th
eigenstate at a given :math:`\mathbf{k}`-point with a corresponding
energy :math:`E_{\mathbf{k}}^n`, :math:`\Omega` is the cell volume,
:math:`\mathbf{\hat{q}}` is the direction of polarization of the photon
and :math:`\hbar\omega` its energy. Currently only the :math:`\Gamma`
point is included in the sum over :math:`\mathbf{k}`-points.

As the position operator is ill-defined in periodic boundary conditions,
this should instead be calculated using a momentum operator formalism,
where the two are related via [Read1991]_:

.. math:: \langle\phi_f|\mathbf{r}|\phi_i\rangle = \frac{1}{i\omega m}\langle\phi_f|\mathbf{p}|\phi_i\rangle + \frac{1}{\hbar\omega}\langle\phi_f|\left[\hat{V}_{nl},\mathbf{r}\right]|\phi_i\rangle .

The commutator term can then be found using the
identity [Motta2010]_:

.. math::

   \begin{aligned}
   &&\left(\nabla_\mathbf{k}+\nabla_\mathbf{k'}\right)\left[\int e^{-i\mathbf{k}\cdot\mathbf{r}} V_{nl}\left(\mathbf{r},\mathbf{r'}\right) e^{i\mathbf{k'}\cdot\mathbf{r'}} d\mathbf{r}\ d\mathbf{r'}\right] \\
   &=&i\int e^{-i\mathbf{k}\cdot\mathbf{r}}\left[V_{nl}\left(\mathbf{r},\mathbf{r'}\right)\mathbf{r'}-\mathbf{r}V_{nl}\left(\mathbf{r},\mathbf{r'}\right)\right] e^{i\mathbf{k'}\cdot\mathbf{r'}} d\mathbf{r}\ d\mathbf{r'} \nonumber,\end{aligned}

where the derivative can either be calculated directly or using finite
differences in reciprocal space. Once the matrix elements have been
calculated in this manner, they can then be used to form a weighted
density of states according to equation :eq:`imag_diel`.

Calculating optical absorption spectra in ONETEP
================================================

The calculation of matrix elements for optical absorption spectra is
activated by setting COND\_CALC\_OPTICAL\_SPECTRA to true. The matrix
elements are then calculated at the end of a conduction calculation in
both the valence and joint valence-conduction NGWF basis sets. Various
options can be modified, including the choice of calculating the matrix
elements in either the position or momentum representation, using the
parameter COND\_SPEC\_CALC\_MOM\_MAT\_ELS. For accurate results, the
position operator should only be used for molecules, where the
conduction NGWFs are sufficiently small compared to the size of the unit
cell that they do not overlap with any periodic copies. If using the
momentum formulation, the default behaviour is to also calculate the
commutator between the nonlocal potential and the position operator,
although setting COND\_SPEC\_CALC\_NONLOC\_COMM will switch this off.
Additionally, the method of calculation of the commutator can be
specified using COND\_SPEC\_CONT\_DERIV, so that either a continuous
derivative or finite difference method is employed. If using the finite
difference option, the finite difference shifting parameter can also be
specified using COND\_SPEC\_NONLOC\_COMM\_SHIFT.

Outputs
-------

If the input filename is seed.dat then the matrix elements will be
written to seed\_val\_OPT\_MAT\_ELS.txt and
seed\_joint\_OPT\_MAT\_ELS.txt. These contain the matrix elements
between all states in the :math:`x`, :math:`y` and :math:`z` directions,
and the energies of each state, as well as the transition energy, are
also printed. For calculations in the momentum representation, the real
and imaginary components of the matrix element are printed in the
additional two columns at the end.

[Clark2005] S. J. Clark, M. D. Segall, C. J. Pickard, P. J. Hasnip, M. J. Probert, K. Refson and M. C. Payne, Z. Kristallogr. **220**, 567 (2005).

[Skylaris2005] C.-K. Skylaris, P. D. Haynes, A. A. Mostofi and M. C. Payne, J. Phys. Condens. Matter **17**, 5757 (2005).

[Godby1986] R. W. Godby, M. Schlüter and L. J. Sham, Phys. Rev. Lett **56**, 2415 (1986).

[Gygi1989] F. Gygi and A. Baldereschi, Phys. Rev. Lett **62**, 2160 (1989).

[Ratcliff2011] L. E. Ratcliff, N. D. M. Hine and P. D. Haynes *In Preparation* (2011).

[Pickard1997] C. J. Pickard, Ph.D. thesis, University of Cambridge (1997).

[Read1991] A. J. Read and R. J. Needs, Phys. Rev. B **44**, 13071 (1991).

[Motta2010] C. Motta, M. Giantomassi, M. Cazzaniga, K. Gal-Nagy and X. Gonze, Comput. Mater. Sci. **50**, 698 (2010).
