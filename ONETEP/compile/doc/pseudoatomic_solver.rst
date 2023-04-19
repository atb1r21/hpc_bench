=========================================================
Using the Pseudoatomic Solver to Generate NGWFs
=========================================================

:Author: Nicholas D.M. Hine, University of Warwick
:Date:   September 2011
:Date:   Updated February 2016

What is being calculated?
=========================

When the atomic solver is used [Ruiz-Serrano2012]_, a
Kohn-Sham DFT calculation is performed for a “pseudoatom”. This means
that the pseudopotential of a single isolated ion is used as the
external potential, and the single-electron Kohn-Sham states are solved
self-consistently, for a given XC functional. The resulting states can
be expected to form an ‘ideal’ atomic orbital basis for a given
calculation [Soler2002]_, [Artacho1999]_, [Blum2009]_, [Tarralba2008]_, [Chen2010]_.
In practice they are a good starting point for initial NGWFs, or a
passable fixed basis for calculations without NGWF optimisation, in
which context they should be at least comparable to the basis sets
generated in SIESTA, for example (though note that a high cutoff energy
is required for accurate representation on the grid in such cases).

The pseudoatomic orbitals we are looking for solve the Kohn-Sham
equation:\ 

.. math:: \hat{H}_{\mathrm{KS}}\psi_{n}(\mathbf{r})=\epsilon_{n}\psi_{n}(\mathbf{r})

where\ 

.. math:: \hat{H}_{\mathrm{KS}}=-\frac{1}{2}\nabla^{2}+V_{\mathrm{loc}}(r)+\sum_{i}|p_{i}\rangle D_{ij}\langle p_{j}|

is the Hamiltonian in terms of kinetic, local potential and nonlocal
potential contributions. For a norm conserving pseudopotential the
matrix :math:`D_{ij}` is diagonal and there is one term per angular
momentum channel, so there is generally only one contributing nonlocal
projector for each wavefunction. In PAW and USPs, there may be more, and
nonzero cross-terms :math:`D_{ij}` for :math:`i\neq j`. Because it is a
spherical potential, solutions of this form of the Kohn-Sham equation
are easily separable into an angular part, solved by the spherical
harmonics (:math:`Y_{lm}(\hat{\mathbf{r}})`, or equivalently the real
spherical harmonics :math:`S_{lm}(\hat{\mathbf{r}})`) multiplied by
radial part, which we will be solving explicitly in terms of a basis.

The atomic orbitals are solved in a sphere of radius :math:`R_{c}`. The
most appropriate basis to use therefore consists of normalised spherical
Bessel functions of given angular momentum :math:`l`. The radial parts
of the basis functions, :math:`B_{l,\nu}(r)` are\ 

.. math:: B_{l,\nu}(r)=j_{l}(q_{l,\nu}r)\:/\,\Big[\int_{0}^{R_{c}}|j_{l}(q_{l,\nu}r)|^{2}r^{2}dr\Big]^{\frac{1}{2}}\;,

with :math:`q_{l,\nu}` such that :math:`q_{l,\nu}R_{c}` are the zeros
of the spherical Bessel functions. Thence :math:`B_{l,\nu}(R_{c})=0` and
:math:`\int_{0}^{R_{c}}|B_{l,\nu}(r)|^{2}r^{2}dr=1` for all :math:`\nu`.
The basis would be complete if :math:`\nu` were infinite: in practice it
must be truncated, and the number of functions included is determined by
a kinetic energy cutoff :math:`E_{\mathrm{cut}}`. The criterion
:math:`\frac{1}{2}q_{l,\nu}^{2}<E_{\mathrm{cut}}` determines the largest
:math:`\nu` for each :math:`l`.

In the pseudoatom calculation we are therefore calculating Kohn-Sham
states of the form

.. math:: \psi_{n}(\mathbf{r})=\sum_{\nu}c_{n,\nu}\, B_{l_{n},\nu}(r)\, S_{l_{n}m_{n}}(\hat{\mathbf{r}})\;,

which have eigenvalues :math:`\epsilon_{n}` and occupancies
:math:`f_{n}` which include spin-degeneracy. The occupancies are fixed,
and determined before the main calculation, such that they sum to the
number of valence electrons. Spherical symmetry of the density is
assumed, so the occupancies of all members of a given set of
:math:`m`-degenerate orbitals are always equal — and in fact in practice
the :math:`m` states for a given :math:`l`, :math:`n` pair are
amalgamated into one state with :math:`f_{n}` summed over all the
degenerate :math:`m`\ ’s. For example, for a nitrogen ion with valence
configuration :math:`2s^{2}\,2p^{3}`, we would have :math:`f_{2s}=2`,
:math:`f_{2p}=3`. For this, we therefore need to find the lowest-energy
self-consistent eigenstate of each of :math:`l=0` and :math:`l=1`.
Henceforth we will only consider the radial dependence
:math:`\psi_{n}(r)`. All radial quantities will be considered to have
been integrated over solid angle already, so factors of :math:`4\pi` are
omitted and :math:`\int|\psi_{n}(r)|^{2}r^{2}dr=1` for a normalised
orbital.

We define the local potential through\ 

.. math:: V_{\mathrm{loc}}(r)=V_{\mathrm{psloc}}(r)+V_{H}[n](r)+V_{XC}[n](r)+V_{\mathrm{conf}}(r)

where for a spherical charge distribution
:math:`n(r)=\sum_{n}f_{n}|\psi_{n}(r)|^{2}`, the Hartree potential is
given by\ 

.. math:: V_{H}(r)=\frac{1}{r}\int_{0}^{r}n(r')r'^{2}dr'+\int_{r}^{\infty}n(r')r'\, dr'\;.

and the XC potential is
:math:`V_{XC}[n](r)=\frac{\partial E_{XC}[n]}{\partial n(r)}`.
:math:`V_{\mathrm{conf}}(r)` is an optional confining potential whose
specific form will be discussed later [Blum2009]_.

For each :math:`l` we can define the Hamiltonian matrix\ 

.. math:: H_{\nu,\nu'}^{l}=\int_{0}^{R_{c}}B_{l,\nu}(r)\left[\hat{H}B_{l,\nu'}(r)\right]r^{2}dr

and the overlap matrix\ 

.. math:: S_{\nu,\nu'}^{l}=\int_{0}^{R_{c}}B_{l,\nu}(r)B_{l,\nu'}(r)r^{2}dr

We then solve the secular equation\ 

.. math:: \mathbf{H}^{l}.\mathbf{c}_{n}=\epsilon_{n}\mathbf{S}^{l}.\mathbf{c}_{n}

to give the coefficients :math:`c_{n,\nu}` which describe the orbitals.
The orbitals are generated on the real-space grid and density mixing
with a variable mixing parameter :math:`\alpha` is then used until
self-consistency is obtained. The result is deemed to be converged once
a) the Harris-Foulkes estimate of the total energy (the bandstructure
energy) matches the total energy as determined from the density to
within a given tolerance (:math:`10^{-5}` Ha) and the energy has stopped
changing each iteration, to within a given tolerance (:math:`10^{-7}`
Ha).

Performing a Calculation with the Pseudoatomic Solver
=====================================================

The atomic solver is the default approach to NGWF initialisation, so if
you do not need to change any settings for any species, simply omit the
``%block species_atomic_set`` block.

If there are any tweaks to be made to the default, this block is
required, and for each element symbol the string “SOLVE” should appear
in the entry. If you want to use automatic initialisation of the number
of NGWFs, then specify that to be ``-1`` in the species block. The code
will attempt to determine how many orbitals to use, which orbitals
constitute the valence, and what their default occupancies should be. To
illustrate what will happen, we present some simple examples.

Let us imagine setting up a calculation with only nitrogen, for which
:math:`Z_{\mathrm{atom}}=7` and :math:`Z_{\mathrm{ion}}=5`. The valence
manifold consists of :math:`4` NGWFs of radius :math:`R_{c}=8.0` per
atom, so we would have the following blocks in our input file:

````

``%block species``

``N N 7 4 8.0``

``%endblock species``

``%block species_atomic_set``

``N “SOLVE”``

``%endblock species_atomic_set``

Note that because we may well want to add extra options to this string
later, it’s best to always use the “” quotes around SOLVE. These
settings will activate the pseudoatomic solver and it will attempt to
guess a default configuration for the atom. Since
:math:`Z_{\mathrm{ion}}=5`, the code will count back five electrons from
the end of the default neutral atom occupancy, which is
:math:`1s^{2}\,2s^{2}\,2p^{3}`, and will discover that the valence
states are :math:`2s^{2}\,2p^{3}`. Since we have asked for :math:`N=4`
NGWFs, the solver will then count forward from the start of the valence
states and determine that by including the whole first set of :math:`s`
and :math:`p` states it has enough to span the valence space and create
four orbitals (and thus four NGWFs). The solver will therefore solve for
one state with :math:`l=0`, :math:`f=2` and one state with :math:`l=1`,
:math:`f=3`, all with radius :math:`R_{c}=8.0`, and from these states
will produce one :math:`s`-like NGWF and the three degenerate
:math:`p_{x}`, :math:`p_{y}` and :math:`p_{z}` NGWFs.

A slightly more complex example would be if we were generating orbitals
for iron (:math:`Z_{\mathrm{atom}}=26`, :math:`Z_{\mathrm{ion}}=8)`:

````

``%block species``

``Fe Fe 26 9 10.0``

``%endblock species``

``%block species_atomic_set``

``Fe “SOLVE”``

``%endblock species_atomic_set``

This time, to find the default configuration, the solver initialisation
routines will count back 8 electrons from the neutral atom configuration
of :math:`1s^{2}\,2s^{2}\,2p^{6}\,3s^{2}\,3p^{6}\,3d^{6}\,4s^{2}` and
thus will determine that the valence states are :math:`3d^{6}\,4s^{2}`.
However, this time we have asked for 9 NGWFs, so it will then count
forward from :math:`3d`, include the fivefold-degenerate lowermost
:math:`d`-like state and the lowest :math:`s`-like state. This only
makes six, so it then will also have to include the threefold-degenerate
:math:`4p`-like state. The solver will have to solve for one unoccupied
:math:`p`-like orbital, which will have :math:`f=0` throughout the
calculation.

Controlling the configuration
-----------------------------

The default neutral-atom configurations for all the elements up to
:math:`Z=92` are included in the code, and will be used by default to
generate the configuration. However, it is also possible to override
these default configurations. For example, to generate NGWFs for iron in
the 3+ state, we might want to set the occupancies to
:math:`3d^{5}\,4s^{0}`. To do this we use the “conf=” directive after
the SOLVE string:

````

``%block species_atomic_set``

``Fe “SOLVE conf=3d5 4s0”``

``%endblock species_atomic_set``

Any terms in the configuration which are not overridden are left at
their default values. Another example might be if we wanted to force the
partial occupation of more higher-lying states than would otherwise be
occupied for the neutral atom:

````

``%block species_atomic_set``

``C “SOLVE conf=2s1.5 2p2.5”``

``%endblock species_atomic_set``

Note that the solver counts through the configuration terms strictly in
the order :math:`n,l`, i.e. \ :math:`n` is looped over outermost, then
:math:`l=0` to :math:`l=n-1` for each :math:`n` innermost. This means
that sometimes a little thought is required to get the terms one
actually wants, and not spurious extra ones. For example, if we wanted
to run a calculation of oxygen with 9 NGWFs per atom, what we probably
wanted would be to run with 1 :math:`s`-like NGWF, 3 :math:`p`-like
NGWFs and 5 :math:`d`-like NGWFs. However, this is not by default what
one will get if one asks for

````

``%block species``

``O O 8 9 9.0``

``%endblock species``

``%block species_atomic_set``

``O “SOLVE”``

``%endblock species_atomic_set``

This will identify :math:`2s^{2}\,2p^{4}` as the valence orbitals, and
counting forward will identify :math:`2s`, :math:`2p`, :math:`3s`,
:math:`3p` and just 1 of the 5 degenerate :math:`3d` states as the NGWFs
required. Therefore, we must instruct the atomsolver to ignore the
unwanted excited :math:`3s` and :math:`3p` terms. We do this with an
“X”, which instructs the solver to knock out this term:

````

``%block species_atomic_set``

``O “SOLVE conf=2s2 2p4 3sX 3pX 3d0”``

``%endblock species_atomic_set``

Strictly speaking, the :math:`2s`, :math:`2p` and :math:`3d` strings are
not needed, as they are the default values anyway, but they are left in
for clarity. I find it advisable, so that I can keep track of the terms
which will generate the NGWFs, to add explicitly the terms with zero
occupancy to the conf string.

Generating larger, non-minimal bases
------------------------------------

ONETEP is generally used to create an *in-situ-optimised*, minimal basis
(eg 4 NGWFs/atom for C, N, O etc). However, it is also possible to fix
the NGWFs and run with a much larger, unoptimised basis, in a manner
akin to other DFT codes designed for large-scale simulations (eg
SIESTA). One would then normally want to use multiple NGWFs for each
angular momentum channel corresponding to the valence orbitals. This is
known as using a “multiple-zeta” basis set, where zeta refers to the
radial part of the valence atomic orbitals. For example, a “triple-zeta”
basis for carbon would have :math:`3` :math:`s`-like functions and
:math:`3` of each of :math:`p_{x}`, :math:`p_{y}`, and
:math:`p_{z}`-like functions. There are two approaches to generating
these extra radial functions. This simplest is just generate the
higher-lying orbitals of a given angular momentum. For carbon, for
example, a double-zeta basis in this scheme would include :math:`3s` and
:math:`3p`-like states. This approach, however, is not very quick to
converge with basis size.

It is often better to apply the commonly-used “split-valence” approach.
This allows the orbitals that have been generated to be “split” into
multiple functions, so as to generate so-called “split-valence
multiple-zeta” basis sets. In this formalism, one function :math:`f(r)`
can be split into two functions :math:`g_{1}(r)` and :math:`g_{2}(r)`
according to the following:

#. A matching radius :math:`r_{m}` is chosen: for :math:`r>r_{m}`, we
   set :math:`g_{2}(r)=f(r)`. For :math:`r\leq r_{m}`, we set
   :math:`g_{2}(r)=r^{l}(a_{l}-b_{l}r^{2})`, where :math:`a_{l}` and
   :math:`b_{l}` are chosen such that :math:`g_{2}(r_{m})=f(r_{m})` and
   :math:`g_{2}'(r_{m})=f'(r_{m})`.

#. The other function, :math:`g_{1}(r)`, is set to
   :math:`f(r)-g_{2}(r)`, so :math:`g_{1}(r)=0` for :math:`r\geq r_{m}`.

#. Both functions are renormalised, so
   :math:`\int_{0}^{R_{c}}|g_{1}(r)|^{2}r^{2}dr=1` and
   :math:`\int_{0}^{R_{c}}|g_{2}(r)|^{2}r^{2}dr=1`.

Splitting of a term is activated by adding a colon after the term and
specifying the “split norm” value. This is the fraction :math:`p` of the
total norm of the orbital which is beyond the matching radius
:math:`r_{m}`, such that
:math:`\int_{r_{m}}^{R_{c}}|f(r)|^{2}r^{2}dr=p`. If this colon is
present, the solver will take into account the total number of orbitals
which will result from this term *after splitting,* when counting
forward in the configuration terms to determine which orbitals to solve.
For example, if we wished to generate a Double-Zeta Polarisation (DZP)
basis for oxygen (:math:`2\times1\times s`,
:math:`2\times3\times p`,\ :math:`1\times5\times d`), where the last 15%
of the norm was matched for the :math:`s` and :math:`p`-orbitals, we
would use the following:

````

``%block species``

``O O 8 13 9.0``

``%endblock species``

``%block species_atomic_set``

``O “SOLVE conf=2s2:0.15 2p4:0.15 3sX 3pX 3d0”``

``%endblock species_atomic_set``

So that you can tell that it is happening, the code will output a
message along the following lines when splitting a given orbital.

````

``Splitting orbital 1, splitnorm= 0.150000000``

``Splitting orbital 1, splitnorm= 0.060000000``

The result of the splitting can be viewed in the
“initial\_rad\_ngwf\_xx” files.

Obtaining Polarisation Orbitals through Perturbation
----------------------------------------------------

As well as including more flexibility for the valence orbitals, in the
form of multiple-zeta basis sets, one frequently also wants to expand
the basis by extending it to higher angular momentum channels. This can
be done by simply increasing the number of NGWFs requested and ensuring
through the ’conf=’ string that the extra functions obtained are of the
right angular momentum. However, the resulting high-\ :math:`l` states
tend to be unbound in the free atom, and therefore do not necessarily
add anything particularly useful to the basis.

There is an alternative means to generate higher-\ :math:`l` states,
using perturbation theory. In this, one effectively applies an electric
field to the valence states of angular momentum :math:`l` and polarises
them, resulting in a set of states of angular momentum :math:`l+1`.
Imagine we have an orbital :math:`\psi_{0}(\mathbf{r})` of angular
momentum :math:`l`, :math:`m` which is an eigenstates of the original
Hamiltonian :math:`\hat{H}_{0}` with eigenvalue :math:`\epsilon_{0}`:

.. math:: \psi_{0}(\mathbf{r})=\sum_{\nu}c_{0,\nu}\, B_{l,\nu}(r)\, S_{lm}(\hat{\mathbf{r}})\;.

We wish to polarise this orbital by applying an electric field
:math:`\mathcal{E}` in the :math:`z`-direction:

.. math:: \hat{H}_{1}=\mathcal{E}r\, S_{10}(\hat{\mathbf{r}})\;,

(since :math:`S_{10}(\hat{\mathbf{r}})=z/r`). Perturbing
:math:`\psi_{0}` with :math:`\hat{H}_{1}` gives us no shift in energy to
first-order, since the perturbation is an odd multiplicative function of
:math:`z`, meaning :math:`\epsilon_{1}=0`. What about the change in the
wavefunction? This obeys\ 

.. math:: (\hat{H}_{0}-\epsilon_{0})\psi_{1}(\mathbf{r})=-(\hat{H}_{1}-\epsilon_{1})\psi_{0}(\mathbf{r})\label{eq:pert}

In principle, :math:`\psi_{1}(\mathbf{r})` could have any angular
momentum components, but we can see that in practice it only contains
:math:`L=l\pm1`, since the dipole selection rule excludes all other
channels. We already have :math:`l-1` states in our basis, so we
conclude that :math:`\psi_{1}(\mathbf{r})` need only include
:math:`l+1`, and we can expand :math:`\psi_{1}` in terms of the
:math:`l+1` basis functions:\ 

.. math:: \psi_{1}(\mathbf{r})=\sum_{\nu}c_{1,\nu}\, B_{l+1,\nu}\, S_{l+1,m}(\hat{\mathbf{r}})

Therefore we can generate a shifted Hamiltonian\ 

.. math::
   :label: pert

   H_{\nu,\nu'}^{l+1}=\int_{0}^{R_{c}}B_{l+1,\nu}(r)\left[(\hat{H}^{l+1}-\epsilon_{0})B_{l+1,\nu'}(r)\right]r^{2}dr\;,

and the components of the RHS of Eq. :eq:`pert`\ 

.. math:: D_{\nu}=-\int_{0}^{R_{c}}B_{l+1,\nu}(r)r\psi_{0}(r)r^{2}dr\;.

To solve for :math:`c_{1,\nu}` we just need to invert
:math:`H_{\nu,\nu'}^{l+1}`\ and apply it to :math:`D_{\nu}`, and then
renormalise the result to have a norm of 1.

In practice, polarisation of a given configuration term of angular
momentum :math:`l`, to form a perturbative polarisation orbital for
:math:`l+1`, is achieved by adding “\|P” to the term, for example:

````

``%block species``

``O O 8 13 9.0``

``%endblock species``

``%block species_atomic_set``

``O “SOLVE conf=2s2:0.15 2p4:0.15|P”``

``%endblock species_atomic_set``

So that you can tell that it is happening, the code will output a
message along the following lines when polarising a given orbital:

````

``Polarising orbital 1 to generate l= 2 function (orbital 3)``

Again, the result can be viewed by plotting the relevant
“initial\_rad\_ngwf\_xx” file.

Overriding radii
----------------

By default, the cutoff radius used for all the orbitals of an atom is
the same :math:`R_{c}` as defined in the ``%block species`` entry for
that element. However, we can override this, either for all orbitals, or
for certain angular momentum channels.

To override the radius for all channels, for example to
7.0\ :math:`a_{0}`, would we add the flag “R=7.0” to the SOLVE string:

````

``%block species_atomic_set``

``O “SOLVE conf=2s2:0.15 2p4:0.15 3sX 3pX 3d0 R=7.0”``

``%endblock species_atomic_set``

Or leave the default values for all other channels, but override the
:math:`d`-channel only to 5.0\ :math:`a_{0}`, we would use

````

``%block species_atomic_set``

``O “SOLVE conf=2s2:0.15 2p4:0.15 3sX 3pX 3d0 Rd=5.0”``

``%endblock species_atomic_set``

Adjusting confining potentials
------------------------------

By default, a confining potential is applied, of the form:\ 

.. math:: V_{\mathrm{conf}}(r)=S\,\exp[-w_{l}/(r-R_{c}+w_{l})]/(r-R_{c})^2

where :math:`S` is the maximum height of the confining potential (at
:math:`r=R_{c}`), and :math:`w_{l}` is the width of the region over
which it is applied. By default, :math:`S=100` Ha,
and\ :math:`w_{l}=3.0a_{0}` for all :math:`l`-channels used. These can
also be overridden, either all at once or for specific :math:`l`-values
in the case of :math:`w`.

For example, to set no confining potential on the confined
:math:`d`-orbitals in Zinc, but to keep the default one on all the other
orbitals, we could set :math:`w_{d}=0`:

````

``%block species_atomic_set``

``Zn “SOLVE conf=3d10 4s2 Rd=5.0 wd=0.0”``

``%endblock species_atomic_set``

Or to turn off confinement potentials entirely, and generate
:math:`R_{c}=15.0a_{0}` orbitals to match those generated by CASTEP’s
atomsolver (this should allow direct comparison of energies, given
suitable tweaking of the energy cutoffs so that they exactly match:

````

``%block species_atomic_set``

``O “SOLVE R=15.0 S=0.0”``

``%endblock species_atomic_set``

Note that there can be problems with convergence for certain choices of
confining potential. In particular, if you apply different confining
potentials to different *occupied* orbitals, there will be problems
obtaining agreement between the Harris-Foulkes estimator and the actual
total energy - because the band energy will incorporate the different
confining potentials, but the total energy cannot. The confining
potential on angular momentum channels with no occupied orbitals can
therefore be whatever you like, but those of occupied orbitals must all
match. The exception to this is if the cutoff radius ment of one channel
is less than the onset radius for the others. In this case, there is no
need to apply a confinement to the lower-cutoff channel at all (eg in
the example above for Zinc).

Initial Guess Density: Setting initial charges and spins
--------------------------------------------------------

The atomsolver solutions are by default also used to provide an initial
guess for the valence electron density. This is used to generate the
initial Hamiltonian, which determines the occupation of the orbitals via
Palser-Manolopoulos or other kernel optimisation schemes. Therefore it
is important that this initial density is a reasonably good guess to the
real density.

A superposition of atomic densities is usually fine for neutral systems
without large magnetic moments. Sometimes, however, one needs to adjust
the charges and spins of this initial density. It appears to be a rather
bad idea to actually adjust the orbital occupations of the pseudoatoms
self-consistently: it becomes impossible to disentangle the effect of
changing the orbital from that of changing the density.

A better approach, therefore, is to retain the same pseudoatomic
solutions for the neutral atom, but adjust the orbital occupations only
at the point where they are used to generate the density.

This can be done by specifying “INIT CHARGE=X SPIN=Y” in the SOLVE
string. Either CHARGE or SPIN can be omitted if they are not required.
For example, for a manganese ion with charge +3 and spin 4, we might set

````

``%block species_atomic_set``

``Mn “SOLVE conf=3d5 4s2 wd=7.0 INIT SPIN=4 CHARGE=+3”``

``%endblock species_atomic_set``

The default occupation for the neutral atom is :math:`3d^5`
:math:`4s^2`. However, we ask it to apply a charge of +3, and this will
remove occupation number from the orbitals with the highest energy until
the right charge is obtained. In this case the resulting occupation will
be :math:`3d^4` :math:`4s^0`. Next, the spin is applied, resulting in
:math:`3d_{\uparrow}^4` :math:`3d_{\downarrow}^0`. Note that the charge
is applied first, followed by the spin.

[Ruiz-Serrano2012] A. Ruiz-Serrano, N.D.M. Hine and C.-K. Skylaris, *in press*, (2012).

[Soler2002] J.M. Soler, E. Artacho, J.D. Gale, A. Garcia, J. Junquera, P. Ordejon,
and D. Sanchez-Portal,\ *The SIESTA method for ab initio order-N
materials simulation*, J. Phys. Condens. Matter 14, (2002).

[Artacho1999] E. Artacho, D. Sanchez-Portal, P. Ordejon, A. Garca, and J. M. Soler, *Linear-scaling ab-initio calculations for large and complex systems*, Phys. Status Solidi B 215, 809 (1999).

[Blum2009] V. Blum, R. Gehrke, F. Hanke, P. Havu, V. Havu, X. Ren, K. Reuter, and M. Scheffler: *Ab initio molecular simulations with numeric atom-centered orbitals*, Comput. Phys. Commun. 180, 2175 (2009).

[Tarralba2008] A. S. Torralba, M. Todorovic, V. Brazdova, R. Choudhury, T. Miyazaki, M. J. Gillan, and D. R. Bowler. *Pseudo-atomic orbitals as basis sets for the O(N) DFT code CONQUEST*, J. Phys. Condens. Matt. 20(29), (2008).

[Chen2010] M. Chen, G.-C. Guo, and L. He, *Systematically improvable optimized atomic basis sets for ab initio calculations*, J. Phys. Condens. Matt. 22, 445501, (2010).
