===================================================================
Linear response time-dependent density-functional theory (LR-TDDFT)
===================================================================

:Author: Tim Zuehlsdorff, Imperial College London

Linear Response TDDFT
=====================

The linear response TDDFT (LR-TDDFT) functionality in ONETEP allows the
calculation of the low energy excited states of a system in linear
scaling effort. In contrast to time-evolution TDDFT, where the density
matrix of the system is propagated explicitly in time, LR-TDDFT recasts
the problem of finding TDDFT excitation energies into an effective
non-hermitian eigenvalue equation of the form:

.. math::
   :label: full_tddft

   \begin{pmatrix} \textbf{A} & \textbf{B} \\ \textbf{B} & \textbf{A}\end{pmatrix}\begin{pmatrix} \vec{\textbf{X}} \\ \vec{\textbf{Y}} \end{pmatrix} = \omega \begin{pmatrix} 1 & 0 \\ 0 & -1 \end{pmatrix}\begin{pmatrix} \vec{\textbf{X}} \\ \vec{\textbf{Y}} \end{pmatrix}

where the elements of the block matrices :math:`\textbf{A}` and
:math:`\textbf{B}` can be expressed in canonical Kohn-Sham
representation as

.. math::

   \begin{aligned}
   A_{cv,c'v'}&=&\delta_{c,c'}\delta_{v,v'}(\epsilon^{\textrm{\scriptsize{KS}}}_{c}-\epsilon^{\textrm{\scriptsize{KS}}}_{v})+K_{cv,c'v'} \\
   B_{cv,c'v'}&=&K_{cv, v'c'}\end{aligned}

Here, :math:`c` and :math:`v` denote Kohn-Sham conduction and valence
states and **K** is the coupling matrix with elements given by

.. math::

   \begin{aligned}
    \nonumber
   K_{cv,c'v'}=2\int \mathrm{d} ^3 r \mathrm{d} ^3 r'\left[\frac{1}{|\textbf{r}-\textbf{r}'|}+\left. \frac{\delta^2 E_{\textrm{\scriptsize{xc}}}}{\delta\rho(\textbf{r})\delta\rho(\textbf{r}')}\right|_{\rho^{\{0\}}}\right]  \\
   \times \psi^{\textrm{\scriptsize{KS}}*}_{c}(\textbf{r})\psi^{\textrm{\scriptsize{KS}}}_{v}(\textbf{r})\psi^{\textrm{\scriptsize{KS}}*}_{v'}(\textbf{r}')\psi^{\textrm{\scriptsize{KS}}}_{c'}(\textbf{r}').\end{aligned}

with :math:`E_{\textrm{\scriptsize{xc}}}` being the
exchange-correlation energy. Its second derivative, evaluated at the
ground-state density :math:`\rho^{\{0\}}` of the system, is normally
referred to as the exchange-correlation kernel.

The above equation can be understood as an effective 2-particle
Hamiltonian consisting of a diagonal part of conduction-valence
eigenvalue differences and a coupling term :math:`K_{cv,c'v'}`
connecting individual Kohn-Sham excitations.

In ONETEP, LR-TDDFT is implemented both in terms of the full TDDFT
eigenvalue equation (Eqn. [full\_tddft]) and in the Tamm-Dancoff
approximation, a commonly used simplification to the full non-hermitian
eigenvalue equation, where the off diagonal elements :math:`\textbf{B}`
are set to zero. The problem of calculating the TDDFT excitation
energies thus becomes equivalent to solving the hermitian eigenvalue
equation

.. math:: \textbf{A}\vec{\textbf{X}}=\omega \vec{\textbf{X}}

The Tamm-Dancoff approximation violates time-reversal symmetry and
oscillator strength sum rules and can blue-shift strong peaks in the
spectrum by up to 0.3 eV, however, dark states are typically left almost
unaltered from their corresponding states in the Tamm-Dancoff
approximation.

In the ONETEP code, the Tamm-Dancoff eigenvalue equation is re-expressed
in terms of two sets of NGWFs, one optimised for the valence space
(denoted as :math:`\{ \phi_\alpha\}`) and one optimised for a low energy
subspace of the conduction manifold (denoted as :math:`\{\chi_\beta \}`,
see the documentation of the conduction NGWF optimiation functionality).
Furthermore, the eigenvalue equation is solved iteratively for the
lowest few eigenvalues using a conjugate gradient method. In order to do
so we define the action :math:`\textbf{q}` of operator
:math:`\textbf{A}` acting :math:`\vec{\textbf{X}}` in conduction-valence
NGWF space as

.. math::
   :label: operator

   (q^{\chi\phi})^{\alpha\beta}=(P^{\{\mathrm{c}\}} H^{\chi}P^{\{1\}}-P^{\{1\}} H^{\phi}P^{\{\mathrm{v}\}})^{\alpha\beta}
   +(P^{\{\mathrm{c}\}} V^{\{1\}\chi\phi}_{\textrm{\scriptsize{SCF}}}P^{\{\mathrm{v}\}})^{\alpha\beta}.

where :math:`\textbf{H}^{\chi}` and :math:`\textbf{H}^\phi` are the
Hamiltonians in conduction and valence NGWF representation respectively,
:math:`\textbf{P}^{\{c\}}` and :math:`\textbf{P}^{\{v\}}` denote the
conduction and valence density matrices and :math:`\textbf{P}^{\{1\}}`
is the response density matrix, a representation of the trial vector
:math:`\vec{\textbf{X}}` in conduction-valence NGWF space.
:math:`V^{\{1\}}_{\textrm{\scriptsize{SCF}}}` is the first order
response of the system due to the density
:math:`\rho^{\{1\}}(\textbf{r})` associated with
:math:`\textbf{P}^{\{1\}}`. Under this redefinition of the action
:math:`\textbf{A}` in conduction-valence NGWF space, finding the lowest
:math:`N_\omega` excitation energies is equivalent to minimising

.. math:: \Omega=\sum_i^{N_\omega}\omega_i=\sum_i^{N_{\omega}}\left[  \frac{\textrm{Tr}\left[\textbf{P}^{\{1\}\dagger}_i\textbf{S}^{\chi}\textbf{q}^{\chi\phi}_i\textbf{S}^\phi\right]}{\textrm{Tr}\left[\textbf{P}^{\{1\}^\dagger}_i\textbf{S}^{\chi}\textbf{P}^{\{1\}}_i\textbf{S}^\phi\right]}\right]

with respect to :math:`\left\{ \textbf{P}^{\{1\}}_i\right\}` under the
constraint

.. math::
   :label: ortho
	   
    \textrm{Tr}\left[\textbf{P}^{\{1\}\dagger}_i\textbf{S}^{\chi}\textbf{P}^{\{1\}}_j\textbf{S}^\phi\right]=\delta_{ij}.

If all density matrices involved in the above expressions, ie.
:math:`\textbf{P}^{\{1\}}`, :math:`\textbf{P}^{\{c\}}` and
:math:`\textbf{P}^{\{v\}}` are truncated and thus become sparse, the
algorithm scales as :math:`O(N)` with system size for a fixed number of
excitation energies :math:`N_\omega` and as :math:`O(N_\omega^2)` with
the number of excitation energies required.

A similar algorithm can be derived for the full TDDFT eigenvalue
equation, where we make use of the change of variables
:math:`\textbf{p}=\vec{\textbf{X}}+\vec{\textbf{Y}}` and
:math:`\textbf{q}=\vec{\textbf{X}}-\vec{\textbf{Y}}`. Each TDDFT
excitation then has two effective density matrices,
:math:`\textbf{P}^{\{p\}}` and :math:`\textbf{P}^{\{q\}}`, associated
with it that have the same structure as :math:`\textbf{P}^{\{1\}}` in
the Tamm-Dancoff approximation. The density matrices do obey an updated
orthonormality constraint of the form

.. math::
   :label: ortho2

    \frac{1}{2}\left(\textrm{Tr}\left[\textbf{P}^{\{p\}\dagger}_i\textbf{S}^{\chi}\textbf{P}^{\{q\}}_j\textbf{S}^\phi\right]+ \textrm{Tr}\left[\textbf{P}^{\{q\}\dagger}_i\textbf{S}^{\chi}\textbf{P}^{\{p\}}_j\textbf{S}^\phi\right]\right)=\delta_{ij}

and an analogous expression for the total energy :math:`\Omega` in full
TDDFT can be derived.

Performing a LR-TDDFT calculation
=================================

The LR-TDDFT calculation in ONETEP is enabled by setting the task flag
to TASK=LR\_TDDFT. The LR-TDDFT calculation mode reads in the density
kernels and NGWFs of a converged ground state and conduction state
calculation, so the .dkn, .dkn\_cond, .tightbox\_ngwfs and
.tightbox\_ngwfs\_cond files all need to be present. The most important
keywords in a TDDFT calculation are:

-  | :math:`\tt{lr\_tddft\_RPA}`: T/F.
   | Boolean, default :math:`\tt{lr\_tddft\_RPA}`\ =F. If set to T, the
     code performs a full TDDFT calculation without relying on the
     simplified Tamm-Dancoff approximation.

-  | :math:`\tt{lr\_tddft\_num\_states}`: n
   | Integer, default :math:`\tt{lr\_tddft\_num\_states}=1`.
   | The keyword specifies how many excitations we want to converge. If
     set to a positive integer n, the TDDFT algorithm will converge the
     n lowest excitations of the system.

-  | :math:`\tt{lr\_tddft\_cg\_threshold}`: x
   | Real, default :math:`\tt{lr\_tddft\_cg\_threshold}=10^{-6}`.
   | The keyword specifies the convergence tolerance on the sum of the n
     TDDFT excitation energies. If the sum of excitation energies
     changes by less than x in two consecutive iterations, the
     calculation is taken to be converged.

-  | :math:`\tt{lr\_tddft\_maxit\_cg}`: n
   | Integer, default :math:`\tt{lr\_tddft\_maxit\_cg}=60`.
   | The maximum number of conjugate gradient iterations the algorithm
     will perform.

-  | :math:`\tt{lr\_tddft\_triplet}`: T/F
   | Boolean, default :math:`\tt{lr\_tddft\_triplet}=F`.
   | Flag that decides whether the :math:`\tt{lr\_tddft\_num\_states}=n`
     states to be converged are singlet or triplet states.

-  | :math:`\tt{lr\_tddft\_write\_kernels}`: T/F
   | Boolean, default :math:`\tt{lr\_tddft\_write\_kernels}=T`.
   | If the flag is set to T, the TDDFT response density kernels are
     printed out at every conjugate gradient iteration. These files are
     necessary to restart a LR\_TDDFT calculation.

-  | :math:`\tt{lr\_tddft\_restart}`: T/F
   | Boolean, default :math:`\tt{lr\_tddft\_trestart}=F`.
   | If the flag is set to T, the algorithm reads in
     :math:`\tt{lr\_tddft\_num\_states}=n` response density kernels in
     .dkn format and uses them as initial trial vectors for a restarted
     LR\_TDDFT calculation.

-  | :math:`\tt{lr\_tddft\_restart\_from\_TDA}`: T/F
   | Boolean, default :math:`\tt{lr\_tddft\_trestart\_from\_TDA}=F`.
   | If the flag is set to T and :math:`\tt{lr\_tddft\_RPA}`: T, the
     code will read in already converged density kernels
     :math:`\left\{\textbf{P}^{\{1\}}_i\right\}` and use them as a
     starting guess for a full TDDFT calculation such that
     :math:`\textbf{P}^{\{p\}}_i=\textbf{P}^{\{q\}}_i=\textbf{P}^{\{1\}}`.
     In many cases, the full TDDFT results are similar to the
     Tamm-Dancoff results and this strategy of starting the full TDDFT
     calculation leads to a rapid convergence.

-  | :math:`\tt{lr\_tddft\_init\_random}` T/F
   | Boolean, default :math:`\tt{lr\_tddft\_init\_random}` T.
   | By default, the initial TDDFT eigenvector guesses are initialised
     to random matrices. This yields an unbiased convergence of the
     TDDFT algorithm but can mean that one starts the optimisation
     relatively far away from the minimum. If
     :math:`\tt{lr\_tddft\_init\_random}`\ =F, the code instead computes
     the :math:`n` minimum energy pure Kohn-Sham transitions in
     linear-scaling effort and initialises the :math:`n` TDDFT response
     density matrices to the pure Kohn-Sham transition density matrices.
     In many small to medium sized systems this leads to initial states
     much closer to the TDDFT minimum and rapid convergence. In large
     extended systems this can yield states that are spurious charge
     transfer states that are not ideal, especially if a more advanced
     density matrix truncation scheme is used. In this case it is
     possible to set the keyword
     :math:`\tt{lr\_tddft\_init\_max\_overlap}` T, in which, rather than
     choosing the lowest Kohn-Sham transitions, the code picks the
     lowest few transitions that also have a maximum overlap of electron
     and hole densities.

-  | :math:`\tt{lr\_tddft\_kernel\_cutoff}`: x
   | Real, default :math:`\tt{lr\_tddft\_kernel\_cutoff}=1000 a_0`.
   | Keyword sets a truncation radius on all response density kernels in
     order to achieve linear scaling computational effort with system
     size.

While the LR\_TDDFT calculation can be made to scale linearly for a
fixed number of excitations converged, it should be kept in mind that
the algorithm needs to perform orthonormalisation procedures and thus
scales as :math:`O(N^2)` with :math:`\tt{lr\_tddft\_num\_states}`.

Truncation of the Response density matrix
=========================================

To run a fully linear scaling TDDFT calculation the response density
matrix has to be truncated by setting
:math:`\tt{lr\_tddft\_kernel\_cutoff}`. This truncation introduces
numerical errors into the calculation, which mainly manifest themselves
in the form that the response density matrices do no longer exactly obey
a first order idempotency constraint that is placed on them. The
idempotency constraint can be written in form of an invariance equation:

.. math:: \textbf{P}^{\{1\}'}=\textbf{P}^{\{\mathrm{c}\}}\textbf{S}^{\chi}\textbf{P}^{\{1\}}\textbf{S}^{\phi}\textbf{P}^{\{\mathrm{v}\}}=\textbf{P}^{\{1\}}

To measure the degree to which the invariance relation is violated we
make use of a penalty functional
:math:`Q\left[ \textbf{P}^{\{1\}}\right]` given by:

.. math:: Q\left[\textbf{P}^{\{1\}}\right]=\textrm{Tr}\left[\left(\textbf{P}^{\{1\}\dagger}\textbf{S}^{\chi}\textbf{P}^{\{1\}}\textbf{S}^{\phi}-\textbf{P}^{\{1\}' \dagger}\textbf{S}^{\chi}\textbf{P}^{\{1\}'}\textbf{S}^{\phi} \right)^2\right].

For truncated :math:`\textbf{P}^{\{1\}}`,
:math:`Q\left[\textbf{P}^{\{1\}}\right]\neq 0` which can lead to
problems in the convergence of the conjugate gradient algorithm. In
order to avoid these issues, the TDDFT routines perform the minimisation
of the energy in an analogous form to the LNV method in ground-state
calculations: The auxiliary density kernel :math:`\textbf{P}^{\{1\}'}`
is used instead of :math:`\textbf{P}^{\{1\}}` for the minimisation of
:math:`\Omega`. While :math:`\textbf{P}^{\{1\}'}` is much less sparse
than :math:`\textbf{P}^{\{1\}}` it preserves idempotency to the same
degree as the conduction and valence density kernel, yielding a
stabilised convergence.

However, should :math:`Q\left[\textbf{P}^{\{1\}}\right]` diverge
significantly from 0 during the calculation, there are routines in place
similar to the kernel purification schemes in ground state DFT that
force the kernel towards obeying its idempotency constraint. The keyword
controlling these routines are given below:

-  | :math:`\tt{lr\_tddft\_penalty\_tol}`: x
   | Real, default :math:`\tt{lr\_tddft\_penalty\_tol}=10^{-8}`.
   | Keyword sets a tolerance for the penalty functional. If
     :math:`Q\left[\textbf{P}^{\{1\}}\right]` is larger than
     :math:`\tt{lr\_tddft\_penalty\_tol}` the algorithm will perform
     purification iterations in order to decrease the penalty value and
     force :math:`\textbf{P}^{\{1\}}` towards the correct idempotency
     behaviour.

-  | :math:`\tt{lr\_tddft\_maxit\_pen}`: n
   | Integer, default :math:`\tt{lr\_tddft\_maxit\_pen}=20`.
   | The maximum number purification iterations performed per conjugate
     gradient step.

More advanced TDDFT kernel truncation schemes
=============================================

| There are many situations where physical intuition allows one to
  specify a more sophisticated sparsity pattern than a uniform spherical
  kernel cutoff on :math:`\textbf{P}^{\{1\}}` (or
  :math:`\textbf{P}^{\{p\}}` and :math:`\textbf{P}^{\{q\}}` for full
  TDDFT). For example, in pigment-protein complexes the excitations of
  interest retain a relative localisation on the pigment and one would
  ideally converge these states directly, without obtaining any spurious
  charge transfer states from the pigment to far away regions of the
  protein, that can arise due to failures in semi-local exchange
  correlation functionals. This can be achieved by introducing a new
  block into the input file of the form

::

   %block species_tddft_kernel
     label1 label 2 label3 ...
     label5 ...
     ...
   %endblock species_tddft_kernel

| where the labels refer to atom labels. As an example, consider a
  pigment protein complex, where the pigment atoms are labelled H1, C1
  etc. while the protein atoms are labelled H, C, etc. Then we can force
  the excitations of the system to be fully localised on the pigment by
 including

::

   %block species_tddft_kernel
     C1 H1 ...
   %endblock species_tddft_kernel

| This has the effect of setting all elements of
  :math:`\textbf{P}^{\{1\}}` to zero that correspond to conduction or
  valence NGWFs centered on atoms of the environment. In this way the
  electrostatic effects of the environment are treated fully quantum
  mechanically, while no delocalisation into the protein is allowed. If
  one would like to introduce a coupling to the environment but wants to
  suppress any charge transfer coupling between the pigment and its
  environment, it is possible to specify

::

   %block species_tddft_kernel
     C1 H1 ...
     C H ...
   %endblock species_tddft_kernel

| It is possible to specify an arbitrary number of subregions in the
  system in this way. It is also possible to list the same species in
  different lines, allowing for charge transfer interactions between
  some atom types of two regions but not others.

| Rather than having the off-diagonal charge-transfer blocks defined in
  ``%block species_tddft_kernel`` set exactly to zero,
  it is also possible to give these blocks a more realistic sparsity
  pattern, for example that of the overlap matrix. While this process
  still suppresses any significant amount of charge transfer between
  TDDFT regions, it can be used to allow overlapping NGWFs from
  different TDDFT regions to contribute to the TDDFT transition density.
  In order to do so, set the block

::

   %block species_tddft_ct
     C1 H1 ...
     C2 H2 ...
   %endblock species_tddft_ct

| and set ``lr_tddft_ct_length`` to a chosen cutoff length
  for the charge-transfer interaction between the specified blocks. For
  example, if the off-diagonal blocks of the response density matrix
  (corresponding to charge-transfer excitations between TDDFT regions)
  should have the same sparsity pattern as the overlap matrix, set
  ``lr_tddft_ct_length`` to twice the NGWF localisation
  radius.

Preconditioning
===============

The TDDFT eigenvalue problem is generally ill-conditioned, which can
lead to a relatively slow convergence. For this reason, it is possible
to precondition the eigenvalue problem, which is achieved by solving a
linear system iteratively to a certain tolerance at each conjugate
gradient step. Solving the linear system only requires matrix-matrix
multiplications and is very cheap for small and medium sized systems,
however, it can get more costly for very large systems, especially when
no kernel truncation is used. In these cases, it can be necessary to
reduce the number of default iterations of the preconditioner. The main
keywords controlling the preconditioner are

-  | :math:`\tt{lr\_tddft\_precond}`: T/F
   | Boolean, default :math:`\tt{lr\_tddft\_precond}=T`.
   | Flag that decides whether the preconditioner is switched on or off.

-  | :math:`\tt{lr\_tddft\_precond\_iter}`: n
   | Integer, default :math:`\tt{lr\_tddft\_precond\_iter}=20`.
   | Maximum number of iterations in the linear system solver applying
     the preconditioner.

-  | :math:`\tt{lr\_tddft\_precond\_tol}`: x
   | Real, default :math:`\tt{lr\_tddft\_precond\_tol}=10^{-8}`.
   | The tolerance to which the linear system is solved in the
     preconditioner. Choosing a large tolerance means that the
     preconditioner is only applied approximately during each iteration.

Representation of the unoccupied subspace
=========================================

In the LR\_TDDFT method as implemented in ONETEP, the user has two
options regarding the representation of the unoccupied subspace. The
first option is to define the active unoccupied subspace of the
calculation to only contain the Kohn-Sham states that were explicitly
optimised in the COND calculation. The other is to make use of a
projector onto the entire unoccupied subspace, where we redefine the
conduction density matrix as:

.. math:: \textbf{P}^{\{\textrm{c}\}}=\left(\left(\textbf{S}^{\chi}\right)^{-1} -\left(\textbf{S}^{\chi}\right)^{-1}\textbf{S}^{\chi\phi}\textbf{P}^{\{\textrm{v}\}}\left(\textbf{S}^{\chi\phi}\right)^\dagger \left(\textbf{S}^{\chi}\right)^{-1}\right) .

The first option has the advantage that we only include states for
which the NGWFs are well optimised, but has the drawback that some
excitations converge very slowly with the size of the unoccupied
subspace and thus a good convergence with the number of conduction
states optimised is hard to reach. The second method implicitly includes
the entire unoccupied subspace (to the extent that it is representable
by a small, localised NGWF representation), but has the disadvantage
that now states are included in the calculation for which the NGWFs are
not optimised. Furthermore, the density matrix defined above is no
longer strictly idempotent, leading to violations of the idempotency
condition and thus a non-vanishing penalty functional
:math:`Q\left[\textbf{P}^{\{1\}}\right]`, requiring kernel purification
iterations as described in the previous section.

The problem of loss of idempotency can be avoided by using the joint
NGWF set to represent the conduction space when using the projector.
While this increases the computational cost of the LR\_TDDFT calculation
by a factor of 2, it preserves the idempotency of
:math:`\textbf{P}^{\{\textrm{c}\}}` and is the recommended option when
using the projector onto the unoccupied subspace.

The keywords controlling the use of the projector are

-  | :math:`\tt{lr\_tddft\_projector}`: T/F
   | Boolean, default :math:`\tt{lr\_tddft\_projector}=T`.
   | If the flag is set to T, the conduction density matrix
     :math:`\textbf{P}^{\{\textrm{c}\}}` is redefined to be a projector
     onto the entire unoccupied subspace.

-  | :math:`\tt{lr\_tddft\_joint\_set}`: T/F
   | Boolean, default :math:`\tt{lr\_tddft\_joint\_set}=T`.
   | If the flag is set to T, the joint NGWF set is used to represent
     the conduction space in the LR\_TDDFT calculation.

Calculations in implicit solvent
================================

A TDDFT calculation in implicit solvent is performed in an analogous way
to the implicit solvent calculation in combination with a conduction
optimisation (see the documentation of the conduction optimisation for
further details). By default, the implicit solvent only acts on the
ground state of the system and thus influences the conduction and
valence Kohn-Sham states mixed into the TDDFT calculation. However, a
screening of the response density due to a dynamic dielectric constant
is not included in the calculation. In order to activate dynamic
screening effects in TDDFT, the user can set the keyword
:math:`\tt{lr\_optical\_permittivity}` to the effective dynamic
dielectric constant :math:`\epsilon_\infty` of the system in question.

Outputs
=======

The LR\_TDDFT calculation will produce a number of outputs. At the end
of the calculation, the individual excitation energies and oscillator
strengths will be computed and printed in the main ONETEP output file.
Furthermore, the energies and oscillator strengths are used to generate
a excitation spectrum written to the textfile rootname.tddft\_spectrum.
The peaks in the spectrum are Gaussians of a width controlled by
:math:`\tt{lr\_tddft\_spectrum\_smear}`. Furthermore, by default,
density cube files of the response density, the electron and the hole
density for each excitation are printed out. The LR\_TDDFT code can also
perform an analysis of individual excitations, where the response
density matrix is decomposed into dominant Kohn-Sham transitions. Since
this analysis requires the Kohn-Sham eigenstates and thus a
diagonalisation of the Hamiltonian, it scales as :math:`O(N^3)` and
should not be performed for very large system sizes.

The keywords controlling these outputs are:

-  | :math:`\tt{lr\_tddft\_write\_densities}`: T/F
   | Boolean, default :math:`\tt{lr\_tddft\_write\_densities}=T`.
   | If the flag is set to T, the response density, electron density and
     hole density for each excitation is computed and written into a
     .cube file.

-  | :math:`\tt{lr\_tddft\_analysis}`: T/F
   | Boolean, default :math:`\tt{lr\_tddft\_analysis}=F`.
   | If the flag is set to T, a full :math:`O(N^3)` analysis of each
     TDDFT excitation is performed in which the response density is
     decomposed into dominant Kohn-Sham transitions.

Good practices and common problems
==================================

-  The quality of the TDDFT excitation energies critically depends on
   the representation of the conduction space manifold. Any excitation
   that has a large contribution from an unoccupied state that is not
   explicitly optimised in the COND calculation is not expected to be
   represented correctly in the LR\_TDDFT calculation. In general it is
   advisable to optimise as many conduction states as possible. However,
   high energy conduction states are often very delocalised and only
   representable if the conduction NGWF radius is increased
   significantly, thus leading to poor computational efficiency. In
   practice, there is a tradeoff between computational efficiency and
   the representation of the conduction state manifold (see also the
   documentation on conduction state optimisation on this issue).
   Generally, TDDFT excitations should be converged with respect to both
   the conduction NGWF radius and the number of conduction states
   explicitly optimised.

-  Since the ground state and conduction density kernels are used as
   projectors onto the occupied and unoccupied subspace in LR\_TDDFT,
   one often finds that the inner loop of the SINGLEPOINT and COND
   optimisation has to be converged to a higher degree of accuracy to
   achieve well behaved TDDFT results. It is therefore recommended to
   increase MAXIT\_LNV and MINIT\_LNV from their default value in the
   SINGLEPOINT and COND calculation. If no density kernel cutoff is
   used, the penalty functional value in the LR\_TDDFT calculation
   should be vanishingly small. If the number increases significantly
   during a calculation or if the code begins to perform penalty
   optimisation steps, that is a clear sign that the initial conduction
   and valence density kernels are not converged well enough.

-  In order to perform a LR\_TDDFT calculation that scales fully
   linearly with system size, all density matrices involved have to be
   sparse and thus a KERNEL\_CUTOFF has to be set for both the
   SINGLEPOINT and COND calculation. Using a density matrix truncation
   on the conduction states can sometimes be difficult depending on how
   the subspace of optimised conduction states is chosen and care has to
   be taken to prevent unphysical results.

-  When running calculations in full linear scaling mode, the ground
   state and conduction density kernels are no longer strictly
   idempotent, which means that the penalty functional in LR\_TDDFT will
   no longer be strictly zero. The code might perform penalty functional
   optimisation steps to keep the idempotency error small. However,
   these idempotency corrections can cause the conjugate gradient
   algorithm to stagnate and can even cause the energy to increase. If
   this happens, it is an indication that the minimum energy and maximum
   level of convergence for this truncation of the density kernel has
   been reached.

-  When placing a truncation onto the the response density kernels it
   should be kept in mind that this may cause the optimisation to miss
   certain low energy excitations completely. Very long range
   charge-transfer type excitations cannot be represented by a truncated
   response density kernel and will thus be missing from the spectrum of
   excitations converged. However, well localised excitations should be
   unaffected. In a similar way, if the TDDFT kernel is limited to a
   certain region, it should be checked whether increasing the region
   leads to a smooth convergence of the energy of the localised state
   within the region.

Reference
=========

For further background regarding the theory behind the LR\_TDDFT method
in ONETEP, as well as a number of benchmark tests, see

-  Linear-scaling time-dependent density-functional theory in the linear
   response formalism, T. J. Zuehlsdorff, N. D. M. Hine, J. S. Spencer,
   N. M. Harrison, D. J. Riley, and P. D. Haynes, J. Chem. Phys.
   **139**, 064104 (2013)
