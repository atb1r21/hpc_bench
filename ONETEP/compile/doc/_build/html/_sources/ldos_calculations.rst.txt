==========================================================================================================
Calculating the Local/Partial Density of States and Angular Momentum Projected Density of States
==========================================================================================================

:Author: Nicholas D.M. Hine, University of Warwick (originally Imperial College London)
:Author: Jolyon Aarons, University of Warwick

:Date: June 2019 (Updated by Jolyon Aarons to add angular momentum PDOS information).
:Date: Originally written by Nicholas D.M. Hine April 2012. 

What is being calculated?
=========================

The local density of states (LDOS) provides a description encompassing
both the spatial and energetic distribution of the single-particle
eigenstates simultaneously. The angular momentum projected density of
states (PDOS) decomposes the density of states energetic distribution
into angular momentum components. Both decompositions may be combined
into an LPDOS. The LDOS and PDOS can thus be two valuable sources of
information for understanding and interpreting electronic structure
calculations.

In the local-orbital framework of ONETEP
[Skylaris2005]_, both decompositions are
achieved by first performing a diagonalisation of the Hamiltonian matrix
in the NGWF basis. This is a post-processing step performed at the end
of the calculation, once NGWF and density kernel convergence have been
achieved. While this comes with a :math:`O(N^{3})` computational cost,
the prefactor is low because the NGWF basis is generally quite small.
Therefore, such a diagonalisation remains fast up to quite large system
sizes, particularly if a parallel eigensolver such as ScaLAPACK is used.

The generalised eigenproblem that needs to be solved to provide the
eigenvalues and eigenvectors is:\ 

.. math::
   :label: gen_eig_prob

   \sum_{\beta}H_{\alpha\beta}M_{\phantom{\beta}n}^{\beta}=\epsilon_{n}\sum_{\beta}S_{\alpha\beta}M_{\phantom{\beta}n}^{\beta}

The matrix :math:`M_{\phantom{\beta}n}^{\beta}` describes the
eigenvectors, which take the form
:math:`|\psi_{n}\rangle=\sum_{\beta}|\phi_{\beta}\rangle
M_{\phantom{\beta}n}^{\beta}`. In ONETEP, Eq. :eq:`gen_eig_prob` can be
solved using either the LAPACK routine DSYGVX or (preferably) the
ScaLAPACK routine PDSYGVX, depending on whether -DSCALAPACK has been
specified at compile time, and ScaLAPACK libraries have been provided.

The result is the eigenvalues :math:`\{\epsilon_{n}\}` and eigenvectors
:math:`M_{\phantom{\beta}n}^{\beta}`. From these, the total density of
states can be obtained as

.. math::
   :label: DOS

   D(\epsilon)=\sum_{n}\delta(\epsilon-\epsilon_{n})\;.

In practice the delta function is replaced with a Gaussian with
broadening :math:`\sigma`, typically of the order of
:math:`0.1\,\mathrm{eV}`:

.. math:: \delta(\epsilon-\epsilon_{n})\approx\sqrt{\frac{\log(2)}{\pi\sigma^2}}*\exp{\left(\frac{-\log(2)(\epsilon-\epsilon_n)^2}{\pi\sigma^2}\right)}.

Local Density of States
-----------------------

The local density of states in a given region :math:`I` is calculated by
projecting each eigenstate onto the local orbitals of region :math:`I`,
as

.. math::
   :label: LDOS

   D_{I}(\epsilon)=\sum_{n}\delta(\epsilon-\epsilon_{n})\;\langle\psi_{n}|\sum_{\alpha\in I}\left(|\phi^{\alpha}\rangle\langle\phi_{\alpha}|\right)|\psi_{n}\rangle.

Here, the non-orthogonality of the NGWFs when used as projectors means
that the ket must be contravariant. We are therefore implicitly using
the contravariant dual of the NGWF (as in DFT+U
[O-Regan2011]_, [O-Regan2012]_).
Fortunately, the functions :math:`|\phi^{\alpha}\rangle` need not be
explicitly constructed in real space: the relationship
:math:`\langle\phi_{\alpha}|\phi^{\beta}\rangle=\delta_{\alpha\beta}`
implies that we can re-write Eq. :eq:`LDOS` as:

.. math::
   :label: LDOS2

   \begin{aligned}
   D_{I}(\epsilon) & = & \sum_{n}\!\delta(\epsilon-\epsilon_{n})\!\!\!\!\!\!\sum_{\beta,\gamma,\alpha\in I}\!\!\!\!(M^{\dagger})_{n}^{\phantom{n}\gamma}\langle\phi_{\gamma}|\left(|\phi^{\alpha}\rangle\langle\phi_{\alpha}|\right)|\phi_{\beta}\rangle M_{\phantom{\beta}n}^{\beta}\nonumber \\
    & = & \sum_{n}\delta(\epsilon-\epsilon_{n})\sum_{\alpha\in I}(M^{\dagger})_{n}^{\phantom{n}\alpha}(\sum_{\beta}S_{\alpha\beta}M_{\phantom{\beta}n}^{\beta})\end{aligned}

What is therefore obtained is a series of functions
:math:`D_{I}(\epsilon)` for each of the chosen regions :math:`I`, which
may be the NGWFs of a single atom, or those of a group of atom types.

Angular Momentum Projected Density of States
--------------------------------------------

To calculate the PDOS, we need to insert an additional resolution of the
identity into equation :eq:`LDOS2` using a basis of angular momentum
resolved functions on which to project our NGWFS,
:math:`| \chi_{l,m}\rangle`:

.. math::
   :label: DOS_identity_operator
	   
   D_{l,I}(\epsilon) \approx \sum_n  \delta(\epsilon-\epsilon_n) \sum_{\alpha,l\in I}(M^{\dagger})_n^{\,\,\,\,\alpha} \sum_{m \in l}\langle{\phi_\alpha
   | \chi'_{\alpha l m}}\rangle \sum_{l'm'} \Lambda^{ l m, l'm'} \sum_\beta \left(\langle{ \chi'_{
   l' m'} |\phi_\beta}\rangle M^\beta_{\ \, n}   \, \right),

where we need to include the overlap matrix of angular momentum
resolved functions, :math:`\Lambda`, since this basis is also
non-orthogonal.

We have considerable scope in which basis we choose for the angular
momentum resolved functions. Effectively, this is a set of spherical
harmonics multiplied by some radial term. In ONETEP, we currently have
two options implemented for the radial term: either spherical waves, or
pseudo-atomic functions, as used to initialise the NGWFs, before
optimisation in the NGWF SCF loop. More details about the theory behind
these options as well as tests and comparisons can be found in our paper
[Aarons2019]_.

Performing an LDOS Calculation
==============================

An LDOS calculation is performed as part of the optional post-processing
activated using ``do_properties: T`` or using ``task: PROPERTIES``.
To activate LDOS we then need to specify the Gaussian broadening, such
as ``dos_smear : 0.1 eV``. The default value of
``dos_smear : -0.1 eV`` disables LDOS.

Then, we need to specify the groups of atom types. This is done via a
block, with each line listing a group of atoms. For example, in a
benzene ring, we might use the following to find the contributions of
the carbon and hydrogen atoms respectively:

::

   %block species_ldos_groups
     C
     H
   %endblock species_ldos_groups

A more complex example would be for a GaAs nanorod with hydrogen
termination on the faces. If we wished to see the LDOS varying over 5
layers, labelled 1-5, we could use:

::

   %block species_ldos_groups
     Ga1 As1 H1
     Ga2 As2 H2
     Ga3 As3 H3
     Ga4 As4 H4
     Ga5 As5 H5
   %endblock species_ldos_groups

Examples of the use of LDOS analysis, including example plots, can be
found in several recent papers employing ONETEP
[Avraam2011]_, [Avraam2012]_, [Hine2012]_.

Performing a PDOS Calculation
=============================

The default settings in ONETEP for PDOS calculations are to use the
pseudo-atomic states as the angular momentum resolved projection basis
with a Löwdin orthogonalisation. For most applications, the spilling
parameter associated with this basis will be sufficiently small.
PDOS calculations are enabled as part of the optional post-processing by
writing ``do_properties : T`` into the ONETEP input file, along with a
Gaussian smearing width, such as ``dos_smear : 0.1 eV`` and a maximum
angular momentum in the angular momentum resolved projection basis, such
as ``pdos_max_l : 2`` to include up to d-states (when using the
default pseudo-atomic basis, this will also be limited by the maximum
angular momentum state in each species, as calculated by the
pseudo-atomic solver).

Local PDOS (LPDOS) Calculations
-------------------------------

If you intend to calculate an angular momentum projected DOS on a subset
of atoms, this can be achieved by specifying a block in the input file,
in the same way as for LDOS. The block can be set up by using the
``species_pdos_groups`` keyword. For example, in a benzene ring
calculation, if you want to find the contributions to the PDOS coming
from solely carbon atoms and solely hydrogen atoms, you could write:

::

   %block species_pdos_groups
     C
     H
   %endblock species_pdos_groups

If you also want the combined contribution from carbon and hydrogen to
each angular momentum channel, you should add a line for this:

::

   %block species_pdos_groups
     C
     H
     C H
   %endblock species_pdos_groups

This will calculate PDOS histogram data up to ``pdos_max_l`` for each
line. As many or as few combinations of species as you require can be
calculated by adding extra lines.

If you instead want a specific subset of atoms of a particular species,
this can be achieved easily by labelling this subset differently to the
others in its species. For example, if you have the following species
specification:

::

   %block species
     Pt Pt 78 9 9.0
   %endblock species

   %block species_atomic_set
     Pt SOLVE conf=5d9 6s1 6p0
   %endblock species_atomic_set

   %block species_pot
     Pt platinum.paw
   %endblock species_pot

   %block species_pdos_groups
     Pt
   %endblock species_pdos_groups

   %block positions_abs
     Pt 8.8292 12.2847 8.7330
     Pt 9.2819 11.1839 11.1325
   %endblock positions_abs

| Then you may wish to duplicate the platinum species definitions to
  label a subset of the atoms in the ``positions_abs`` block, as shown
  here for example:

::
  
   %block species
     Pt Pt 78 9 9.0
     Pt1 Pt 78 9 9.0
   %endblock species

   %block species_atomic_set
     Pt SOLVE conf=5d9 6s1 6p0
     Pt1 SOLVE conf=5d9 6s1 6p0
   %endblock species_atomic_set

   %block species_pot
     Pt platinum.paw
     Pt1 platinum.paw
   %endblock species_pot

   %block species_pdos_groups
     Pt
     Pt1
     Pt Pt1
   %endblock species_pdos_groups

   %block positions_abs
     Pt  8.8292 12.2847 8.7330
     Pt1 9.2819 11.1839 11.1325
   %endblock positions_abs

| and hence calculate PDOS contributions for subsets of atoms of a
  single species.

Expert PDOS Options
-------------------

Further options are available in the ONETEP PDOS functionality to
control the quality of the projection. These will be unneeded in most
cases, but if, for instance, you are observing larger spilling
parameters than your requirements permit, you may wish to enable some of
these options.

The most reliable way we have found to reduce the spilling parameter is
to use a spherical-wave basis rather than the pseudo-atomic basis as the
angular momentum resolved projection basis. To do this in ONETEP, add
``pdos_pseudoatomic : F`` to your input file. By default, this will
create a set of contracted spherical waves by fitting spherical waves to
your converged NGWFs, via the contraction coefficients.

The spherical wave basis is contracted by default to reduce the memory
requirements of the code. You may, however, not see an improvement in
the spilling parameter by using this set. To be certain of reducing the
spilling parameter, you should also opt to use the full, non-contracted
spherical wave basis, by setting ``pdos_reduce_sws`` : T in your input
file, along with an adequately large ``pdos_max_l``. For ``pdos_max_l``
you can start by running with 2 and increase to 3 if required. If you
choose to take this approach, beware of the memory requirements, which
can be *up to* 10 times greater.

If you choose to use a contracted set, then you almost certainly want to
use the default fitting coefficients (fitted to NGWFs). These can be
changed to unity by setting ``pdos_construct_basis : F``, however this
is not likely to improve your results, and is likely to be removed in a
future version of ONETEP due to it being mainly of use for debugging
purposes. To reduce the spilling parameter with the contracted set, we
recommend increasing the ``pdos_max_l`` parameter.

In specialised cases, you may also wish to *not* sum over the magnetic
quantum number. This can be achieved by setting ``pdos_sum_mag : F``.
This will give histogram data for every magnetic quantum number of every
angular momentum channel of each atom group.

Interpreting Outputs
--------------------

ONETEPs PDOS outputs are written to several files as well as to stdout,
which is itself usually redirected to the main log/output file. The
PDOS output to stdout will look something like the following (for an
input file with 3 ``pdos_groups`` and ``pdos_max_l=2``):

::

      
    ================ Projected Density of States (pDOS) calculation ================

    Constructing AM resolved functions  ...... done

    Performing overlap integrals ...  done

    Computing pDOS weights ...  done

    All bands spilling parameter =   2.16 %
    Occupancy-weighted spilling parameter =   0.30 %

     => Outputting data for OptaDOS <=

    Writing pDOS weights to file "Pt3O.val_pdos_bin" ... done

    Writing band gradients to file "Pt3O.val_dome_bin" ... done

    Writing Castep output cell file to "Pt3O-out.cell" ... done

     => Computing Gaussian smeared pDOS <=
    Writing "Pt3O_PDOS.txt" ...  done

     => Computing Occupancy-weighted Gaussian smeared pDOS <=
    Writing "Pt3O_occ_PDOS.txt" ...  done
      => Band centres:
     S band centre of group 1:  -10.784858 eV
     P band centre of group 1:   -6.380333 eV
     D band centre of group 1:   -1.992269 eV
     S band centre of group 2:   -3.492084 eV
     P band centre of group 2:   -5.494629 eV
     D band centre of group 2:   -1.992269 eV
     S band centre of group 3:  -20.033217 eV
     P band centre of group 3:   -6.607254 eV
      Band centres done. <=
      => Integrated number of electrons in each AM band:
     S num electrons of group 1:    3.769061
     P num electrons of group 1:    5.624284
     D num electrons of group 1:   26.497454
     S num electrons of group 2:    2.107330
     P num electrons of group 2:    1.147080
     D num electrons of group 2:   26.497454
     S num electrons of group 3:    1.661731
     P num electrons of group 3:    4.477204
      Integrated number of electrons done. <=
    ================================================================================

First, we can see the spilling parameters – effectively how well the
angular momentum resolved basis is able to represent the NGWFs. A lower
value is better; if you are running production calculations, you should
want a value lower than a few percent. If not, consider making some of
the changes suggested in the expert options section above. The all-bands
value includes un-occupied bands as well as valence states, which will
be the same as the occupancy-weighted version unless you are using
EDFT with a finite electronic temperature.

Following this are the files for use with **OptaDOS**. OptaDOS is a
freely available piece of software for plotting various DOS projections.
If you intend to use it, then please follow the CASTEP section of the
OptaDOS manual with these files, as they are compatible.

The histogram files are then written, “\*\_PDOS.txt”. These also come in
all-bands and occupancy weighted flavours (the occupancy weighted
variant is more reliable for a usual ground state calculation with
ONETEP as the conduction states are not well described. Only if you are
doing LPDOS on the output of a ONETEP conduction calculation the
occupancy un-weighted LPDOS outputs will be meaningful). The order of
columns is firstly the energy column, followed by the angular momentum
component columns (i.e. s,p,d...) for each pdos group. This can be
plotted trivially with xmgrace, or any other plotting tool.

Finally ONETEP reports the energy and occupancy weighted averages of the
PDOS, so called-band centres, useful in catalysis (e.g. the value
“d-band centre” is a very useful decsriptor about the ability of a metal
surface to bind atomic oxygen and other types of adsorbates) and the
integrated number of electrons in each component.

[Skylaris2005] C.-K. Skylaris, P. D. Haynes, A. A. Mostofi, and M. C. Payne, J. Chem. Phys. **122**, 084119 (2005).

[O-Regan2011] D. D. O’Regan, M. C. Payne and A. A. Mostofi, Phys. Rev. B **83**, 245124 (2011).

[O-Regan2012] D. D. O’Regan, N. D. M. Hine, M. C. Payne and A. A. Mostofi, Phys. Rev. B **85**, 085107 (2012).

[Avraam2011] P. W. Avraam, N. D. M. Hine, P. Tangney, and P. D. Haynes, Phys. Rev. B **83**, 241402(R) (2011).

[Avraam2012] P. W. Avraam, N. D. M. Hine, P. Tangney, and P. D. Haynes, Phys. Rev. B **85**, 115404 (2012).

[Hine2012] N. D. M. Hine, P. W. Avraam, P. Tangney, and P. D. Haynes, J. Phys. Conf. Ser. (2012).

[Aarons2019] J. A. Aarons, L. G. Verga, N. D. M. Hine, and C.-K. Skylaris, Submitted (2019).
