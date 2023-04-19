==========================================================
DFT+\ :math:`U`\ (+\ :math:`J`)
==========================================================

:Author: David D. O'Regan, Trinity College Dublin
	 
:Date:   July 2015

DFT+\ :math:`U` is fully and self-consistently implemented in ONETEP,
together with a number of advanced ancillary functionalities. The method
is linear-scaling with respect to system size, exhibiting no systematic
tendency to slow convergence to the ground-state. DFT+\ :math:`U` in its
conventional fixed-projector form introduces only a small increase in
computational pre-factor with respect to the underlying
exchange-correlation functional [O-Regan2012]_.

**PLEASE NOTE: Seven columns are now required in the Hubbard block in
order to allow for** :math:`+J` **calculations. Older input files with six
columns will not yield incorrect results, but the code will exit.**

A very short introduction to DFT+\ :math:`U`
============================================

DFT+\ :math:`U` [Anisimov1991]_, [Anisimov1997]_, [Dudarev1998]_, also
known as LDA+\ :math:`U` or LSDA+\ :math:`U`, is a method used to
improve the description of so-called strongly correlated materials
offered by DFT within conventional approximations for
exchange-correlation (XC) such as the LSDA and :math:`\sigma`-GGA,
quantitatively and even qualitatively. These functionals, based on the
locally-evaluated density and its gradients, can sometimes fail to
reproduce the physics associated with localised orbitals of :math:`3d`
and :math:`4f` character characteristic of conventionally-classed
strongly correlated materials, a category consisting of not only
first-row transition metals and their oxides, but also lanthanoid oxide
materials, and other materials such as certain magnetic semiconductors
and organometallic molecules.

Typically, the LDA and its extensions underestimate local magnetic
moments and the tendency to favour high-spin ground-states in such
materials, and the insulating gap in cases where it is related to
electron localisation. Underestimation of the gap due to the absence, in
the LDA, of the derivative discontinuity with respect to orbital
occupancy in the exact XC-funtional may be confounded by an
underestimation of the exchange splitting induced by local magnetic
moments.

The DFT+U correction term is usually thought of as an explicit
mean-field treatment of the exchange-correlation energy contributed by
the correlated sites (subspaces projected out with functions of
:math:`3d` and or :math:`4f` character) within the Hubbard model,
including a double-counting correction for that contribution already
included in the LDA term. The flavour implemented in ONETEP is the
basis-set independent, rotationally invariant quadratic penalty
functional of Ref [Cococcioni2005]_, defined by the
additive energy correction

.. math::

   E_{DFT+U} \left[ n^{(I) (\sigma)} \right] =  \sum_{I \sigma} \frac{U^{(I)}}{2} \rm{Tr} 
   \left[  n^{(I) (\sigma)} \left( 1 -  n^{(I) (\sigma)} \right)\right].

Here, :math:`U` is an estimate of the scalar screened density-density
Coulomb repulsion between localised orbitals. The occupancy matrix of
the correlated site :math:`I`, for spin channel :math:`\sigma`, is
defined, in the case of orthonormal projector functions :math:`\lbrace \lvert \varphi^{(I)}_m \rangle \rbrace`, and density-matrix
:math:`\hat{\rho}^{(\sigma)}`, by

.. math::

   n^{(I)(\sigma)}_{m m'} = \langle \varphi_m^{(I)} \rvert \hat{\rho}^{(\sigma)} 
   \lvert \varphi_{m'}^{(\sigma)} \rangle.

Put simply, if the system under study comprises open :math:`3d` or
:math:`4f` sub-shells, then there is a good chance that the LDA will
find a minimum energy by partly occupying and leaving degenerate the
Kohn-Sham orbitals strongly overlapping with these states, rather than
splitting them into occupied and virtual Hubbard bands. This leads to an
underestimation of the insulating gap and any associated magnetic order.
In this case, the DFT+\ :math:`U` method can be used to penalise the
non-integer occupancy of these orbitals, tending to fill states with
occupancy greater than :math:`0.5` and empty states with occupancy less
than :math:`0.5`, as can be seen from the expression for the
DFT+\ :math:`U` potential

.. math::

   \hat{V}^{(\sigma)}_{DFT+U} = \sum_{I}  U^{(I)} 
    \lvert \varphi_m^{(I)} \rangle 
   \left( \frac{1}{2} \delta_{m m'} - n^{(I) (\sigma)}_{m m'} \right)  \langle 
   \varphi_{m'}^{(I)} \rvert .

The DFT+\ :math:`U` term may be considered as a correction which cancels
the contribution to the energy arising due to the spurious
self-interaction of a partially occupied
orbital [Cococcioni2005]_. In this case, the :math:`U`
parameter is the curvature of the total energy with respect to the
occupancy of the correlated manifold - which should be a piece-wise
linear curve were Janak’s theorem satisfied [Janak1978]_ –
which can be computed using linear-response theory (among other methods
such as constrained DFT) according to the prescription given in
Refs. [Cococcioni2005]_, [Kulik2006]_.

How to activate DFT+\ :math:`U` in ONETEP
=========================================

In order to activate the DFT+\ :math:`U` functionality, the **hubbard**
block is added to the input file. For example, in the case of a system
containing iron and cerium atoms incorrectly described by the
exchange-correlation functional, which we suspect could benefit from the
DFT+\ :math:`U` correction to improve the description of localisation,
we might use the ``hubbard`` block:

::

   % block hubbard
     Fe1   2   4.0   0.0  -10.0   0.00   1.0
     Fe2   2   4.0   0.0  -10.0   0.00  -1.0
     Ce1   3   6.0   0.0  -10.0   0.50   0.0
   % endblock hubbard

The columns of the ``hubbard`` block are described as follows:

#. The species label, e.g. :math:`Fe1` for iron atoms of a first type in
   the cell, :math:`Fe2` for iron atoms of a second kind, etc. Only the
   species to which orbitals to be corrected by DFT+\ :math:`U` are
   assigned should be listed in this block.

#. The angular momentum channel of the projectors used to delineate the
   strongly correlated sites on Hubbard atoms this type, e.g.
   :math:`l=2` for :math:`Fe1`. Conventionally, the radial quantum
   number :math:`r=l+1` is used to generate atom-centred atomic
   projectors, so that :math:`l=2` gives :math:`3d` orbitals,
   :math:`l=3` gives :math:`4f` orbitals etc. (please get in contact if
   you need to use a :math:`r \ne l+1` combination, or multiple
   sub-shells per atom).

#. The value of the Hubbard :math:`U` for this sub-shell, in
   electron-volts. Most users will simply work with the value for
   :math:`U` that they find corrects the band-gap or bond-lengths in the
   system they wish to study. Methods do, however, exist to estimate its
   value, for example the linear-response technique
   [Cococcioni2005]_, [Kulik2006]_, which is implemented in
   ONETEP.

#. The value of the Hund’s exchange :math:`J` for this sub-shell, in
   electron-volts. The rotationally invariant exchange corrective term
   described in detail in Ref. [Himmetoglu2011]_ is fully
   implemented in ONETEP (including forces etc), and activated for any
   :math:`J \ne 0`.

#. This number :math:`Z` selects how the radial part of the projector
   functions used to describe the :math:`1s`, :math:`2p`, :math:`3d` or
   :math:`4f` atomic orbitals entering the DFT+\ :math:`U` functional
   are defined. In the case that :math:`\mathbf{ Z < 0}`, a subset of
   the orbitals generated by solving the atomic problem subject to the
   pseudopotential for the species in question are chosen (in which case
   the projectors form a subset of the initial guesses for the ONETEP
   NGWFs); here the magnitude of the negative Z makes no difference. In
   the case that :math:`\mathbf{ Z > 0}`, for more advanced users, this
   number is the effective charge divided by the ratio of effective
   masses used to generate projectors in the form of solutions to the
   hydrogenic Schrödinger equation. A good guess for this number might
   be the Clementi-Raimondi effective charge, tabulated in
   Refs. [Clementi1963]_, [Clementi1967]_, and the choice of
   radial profile does matter [O-Regan2010]_. In both
   cases, the projectors are effectively renormalised within an
   atom-centred sphere with the same radius as the NGWFs on that atom.

#. An additional potential acting on the subspace in question, the
   prefactor :math:`\alpha` is here entered in electron-volts. This is
   needed, for example, in order to locally vary the potential in order
   to determine the value of :math:`U` which is consistent with the
   screened response in the system with linear-response
   theory [Cococcioni2005]_, [Kulik2006]_, or to break a
   spatial symmetry, such as in a mixed-valence system. In the example
   given, we are additionally penalising the occupancy on cerium
   :math:`4f` atomic orbitals.

#. The spin-splitting factor, in electron-volts, which is deducted from
   the :math:`\alpha` factor for the spin-up channel and added to
   :math:`\alpha` for the spin-down channel. In the example shown here
   we’re promoting spin-up magnetisation for iron atoms :math:`Fe1`, and
   spin-down for :math:`Fe2`. This can be very useful for appropriately
   breaking magnetic symmetries in antiferromagnetic solids or
   open-shell singlet molecules, or for estimating the magnetic
   susceptibility or exchange coupling.

   **N.B.** Users may find the DFT+\ :math:`U` functionality useful in
   cases of systems even when the DFT+\ :math:`U` correction is not
   needed (setting the all :math:`U` parameters to zero does not disable
   the functionality). The implementation offers a very inexpensive
   method for carrying out carefully-defined atom-centred atomic
   population analysis, or breaking symmetries in spin or charge ordered
   systems.

Compatibility
=============

The DFT+\ :math:`U` functionality is fully compatible with almost all
other parts of the ONETEP code, such as listed below, since it simply
involves an additional term in the Hamiltonian and ionic forces. Please
get in touch first if you would like to use a more exotic combination of
these functionalities:

#. Total-energy minimisation and ionic forces

#. Geometry optimisation, molecular dynamics and phonon calculations

#. All other functionals including hybrids and Van der Waals functionals

#. Implicit solvation

#. The PAW formalism and ultrasoft pseudopotentials

#. Constrained DFT

#. Local density of states (including a correlated subspace
   decomposition)

#. Natural bond orbital calculations

#. Conduction-band optimisation and Fermi’s Golden Rule spectra

#. Calculations of changes in electric polarisation

#. Time-dependent DFT

#. Electronic transmission calculations

The extension of the DFT+\ :math:`U` implementation to cluster Dynamical
mean-field theory has also been implemented in ONETEP; for an example of
its capabilities see Ref. [Weber2012]_.

Using NGWFs and projector self-consistency
==========================================

Any reasonable set of localised atomic-like functions may, in principle,
be used for the projectors defining the correlated subspaces in
DFT+\ :math:`U`; the choice is somewhat arbitrary and the description
“atomic orbitals" does not uniquely define them. One possible approach
is to use Wannier functions for the Kohn-Sham orbitals, so that the
correlated subspaces are proper subspaces of the Kohn-Sham Hilbert
space. Indeed, there is numerical evidence to suggest that Maximally
Localised Wannier Functions (MLWFs) [Marzari1997]_, [Souza2001]_,
in particular, provide a basis that maximises a particular measure of
the on-site Coulomb repulsion [Miyake2008]_, and MLWFs are
in common use as a minimal basis with which to construct tight-binding
models from first-principles.

In ONETEP, a set of variationally-optimised nonorthogonal generalised
Wannier functions (NGWFs) are generated as a by-product of total-energy
minimisation. NGWFs exhibit some similar properties to MLWFs and other
flavours of localised Wannier functions, and, for example, can be used
to calculate finite-difference response properties in a similar
way [O-Regan2012-2]_. As they are conveniently available in
ONETEP, we have made it possible to re-use the NGWFs from the end of a
ground-state calculation as a set of Hubbard projectors with which to
define the DFT+\ :math:`U` correction. For this, it was necessary to
develop a tensorially-consistent formulation of DFT+\ :math:`U` in order
to accommodate nonorthogonal projector
functions [O-Regan2011]_; projector nonorthogonality
for a given subspace is automatically compensated for.

In order to ensure that NGWFs with appropriate symmetry are chosen as
Hubbard projectors for a given atom, those :math:`n` NGWFs
:math:`\lvert \phi_\alpha \rangle` that maximise :math:`\sum^n_{m,\alpha }\langle \varphi_m  \rvert  \phi^\alpha \rangle \langle \phi_\alpha \rvert \varphi_m \rangle`, for a given set of
:math:`n` hydrogenic orbitals :math:`\lvert \varphi_m \rangle`, defined
in the ``hubbard`` block, are selected for the task. The keyword
``hubbard_max_iter``, (defaulting to :math:`0`), sets the task to
``HUBBARDSCF``, which performs a self-consistency cycle over the Hubbard
projectors, demonstrated in
Refs. [O-Regan2010]_, [O-Regan2011]_. The density from one
minimisation is re-used at the beginning of the next, and setting
``hubbard_max_iter`` to :math:`2` one can carry out a DFT+\ :math:`U`
calculation using the LDA NGWFs as projectors.

The keywords ``hubbard_energy_tol``, ``hubbard_conv_win``, and
``hubbard_proj_mixing`` are used to manage the Hubbard projector
self-consistency cycle. For convergence, the ground state energy must
deviate less than ``hubbard_energy_tol`` (defaulting to
:math:`10^{-8}Ha`) from one ``HUBBARDSCF`` iteration to the next, over
``hubbard_conv_win`` (defaulting to :math:`2`) iterations. A fraction
``hubbard_proj_mixing`` (defaulting to :math:`0.0`) of the previous
Hubbard projectors may be mixed with the new ones in order to accelerate
the procedure, although this has never been found to be necessary.
Setting ``hubbard_proj_mixing`` to a negative value causes the
projectors to be read in from a ``.tightbox_hub_projs`` file, for
restarting a ``HUBBARDSCF`` calculation or for a variety of
post-processing tasks.

[O-Regan2012] D. D. O’Regan, N. D. M. Hine, M. C. Payne and A. A. Mostofi, Phys. Rev. B **85**, 085107 (2012).

[Anisimov1991] J. Z. V. I. Anisimov and O. K. Andersen, Phys. Rev. B **44**, 943 (1991).

[Anisimov1997] V. I. Anisimov, F. Aryasetiawan, and A. I. Liechtenstein, J. Phys.: Condens. Matter **9**, 767 (1997).

[Dudarev1998] S. L. Dudarev, Phys. Rev. B **57**, 3 (1998).

[Cococcioni2005] M. Cococcioni and S. de Gironcoli, Phys. Rev. B **71**, 035105 (2005).

[Janak1978] J. F. Janak, Phys. Rev. B **18**, 12 (1978).

[Kulik2006] H. J. Kulik, M. Cococcioni, D. A. Scherlis and N. Marzari, Phys. Rev. Lett. **97**, 103001 (2006).

[Himmetoglu2011] B. Himmetoglu, R. M. Wentzcovitch, and M. Cococcioni, Phys. Rev. B,\ **84**, 115108 (2011).

[Clementi1963] E. Clementi and D.L. Raimondi, J. Chem. Phys. **38**, 2686 (1963).

[Clementi1967] E. Clementi, D.L. Raimondi, and W.P. Reinhardt, J. Chem. Phys. **47**, 1300 (1967).

[O-Regan2010] D. D. O’Regan, N. D. M. Hine, M. C. Payne and A. A. Mostofi, Phys. Rev. B **82**, 081102 (2010).

[Weber2012] C. Weber, D. D. O’Regan, N. D. M. Hine, M. C. Payne, G. Kotliar and P. B. Littlewood, Phys. Rev. Lett. **108**, 256402 (2012).

[Marzari1997] N. Marzari and D. Vanderbilt, Phys. Rev. B **56**, 12847 (1997).

[Souza2001] I. Souza, N. Marzari and D. Vanderbilt, Phys. Rev. B **65**, 035109 (2001).

[Miyake2008] T. Miyake and F. Aryasetiawan, Phys. Rev. B **77**, 085122 (2008).

[O-Regan2012-2] D. D. O’Regan, M. C. Payne, and A. A. Mostofi, Phys. Rev. B **85**, 193101 (2012).

[O-Regan2011] D. D. O’Regan, M. C. Payne and A. A. Mostofi, Phys. Rev. B **83**, 245124 (2011).
