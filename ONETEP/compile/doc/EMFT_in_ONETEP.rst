====================================
Embedded Mean-Field Theory
====================================

:Author: Robert J. Charlton, Imperial College London
:Author: Joseph C.A. Prentice, Imperial College London

:Date:   January 2019
:Date:	 Updated by J.C.A. Prentice October 2021

Embedded Mean-Field Theory (EMFT)
=================================

Often in simulations of materials we wish to consider the impact of a
host environment on a system of interest, such as chromophores in
solvent or doped molecular crystals. While the interesting physics or
chemistry may be associated with the subsystem, the effects of the
environment can be significant and warrant description at the quantum
level of theory. However, the cost of applying accurate quantum methods
such as hybrid functionals to potentially large environments can be
restricted by the cost of such methods. Quantum embedding
[Huang2008]_, [Gomes2012]_ methods are intended to combine
an accurate, high-level description of the subsystem of interest
(active) with a cheaper, low-level method for the host environment.

Embedded mean field theory (EMFT) [Fornace2015]_ is an
approach to quantum embedding based on the one-electron density matrix.
To begin, we partition the density matrix into subsystem components,

.. math::
   :label: emft_dm

   \rho = \begin{pmatrix}
           \rho_\text{AA} & \rho_\text{AB}\\
           \rho_\text{BA} & \rho_\text{BB}
       \end{pmatrix},

where A is the active region and B is the inactive environment. The
total energy can be written as

.. math:: E{\left[\rho\right]}=\text{tr}{\left[\rho H_0\right]}+G{\left[\rho\right]},

where :math:`H_0` contains the one-electron terms of the Hamiltonian
and :math:`G{\left[\rho\right]}` contains all two-electron terms (local,
Hartree and exchange-correlation effects). In embedded mean-field theory
(EMFT), the two-electron interaction for the active subsystem A is
constructed at a higher level of theory to the rest of the system,

.. math::
   :label: emft_en
   
   E^\text{EMFT}{\left[\rho\right]}=
           \text{tr}{\left[\rho H_0\right]}+G^\text{low}{\left[\rho\right]}+\left(G^\text{high}{\left[\rho_\text{AA}\right]}-G^\text{low}{\left[\rho_\text{AA}\right]}\right),

where :math:`G^\text{low}` and :math:`G^\text{high}` are the
two-electron interaction energies at the lower and higher levels of
theory, respectively. For example, the low level theory could be LDA
while the higher level uses a hybrid functional such as B3LYP. We assume
here that the core Hamiltonian :math:`H_0` is the same at both levels of
theory, though this need not necessarily be the case. The ground state
of the embedded system can thus be obtained by minimising
:eq:`emft_en` with respect to the elements of the density
matrix.

Block orthogonalisation
-----------------------

Normalisation is maintained provided the trace of the density matrix
with the overlap matrix satisfies

.. math:: \text{Tr}\left[\rho\textbf{S}\right]=N.

EMFT partitioning can result in unrealistic charge spillover from the
low-level to the high-level region, producing large negative results for
the off-diagonal terms
:math:`\text{Tr}\left[\rho_\text{AB}\textbf{S}_\text{BA}\right]` and
:math:`\text{Tr}\left[\rho_\text{BA}\textbf{S}_\text{AB}\right]`. One
possible remedy is to impose a block-orthogonalisation (BO) between the
subsystem orbitals [Ding2017]_,

.. math::

   \begin{aligned}
       {\lvert\tilde{\phi_i^\text{B}}\rangle}&
           =\left(1-\hat{P}^\text{A}\right){\lvert\phi_i^\text{B}\rangle}, \\
       \hat{P}^\text{A}&
           =\sum_{j,k\in\text{A}}{\lvert\phi_j^\text{A}\rangle}
           \left(\textbf{S}^\text{AA}\right)_{jk}^\text{-1}{\langle\phi_k^\text{A}\rvert}.\end{aligned}

By construction
:math:`\text{Tr}\left[\rho_\text{AB}\textbf{S}_\text{BA}\right]` and
:math:`\text{Tr}\left[\rho_\text{BA}\textbf{S}_\text{AB}\right]` are
strictly zero and all electrons are associated with the diagonal blocks.

Implementation in ONETEP
========================

Quantum embedding as implemented in ONETEP is based around
EMFT [Prentice2020]_. Here we denote the active system
NGWFs as :math:`{\lvert\chi_i^\text{A}\rangle}` and the environment NGWFs as
:math:`{\lvert\phi_j^\text{B}\rangle}`. The fundamental quantity of interest is
the Hamiltonian,

.. math::

   \textbf{H}^\text{EMFT}=\begin{pmatrix}
           \textbf{H}^\text{high}_\text{AA} & \textbf{H}^\text{low}_\text{AB} \\
           \textbf{H}^\text{low}_\text{BA} & \textbf{H}^\text{low}_\text{BB}
       \end{pmatrix},
       \label{eq:emft_ham}

where the high- and low-level Hamiltonian operators are given as

.. math::

   \begin{aligned}
       \hat{H}^\text{high}=&\hat{T}+\hat{V}_\text{local}+\hat{V}_\text{Hartree}+\hat{V}_\text{XC}^\text{high},\\
       \hat{H}^\text{low}=&\hat{T}+\hat{V}_\text{local}+\hat{V}_\text{Hartree}+\hat{V}_\text{XC}^\text{low}.\end{aligned}

The total energy can thus be found by minimising the quantity

.. math::

   E^\text{EMFT}=
       \min_{\left\{K^{\alpha\beta}\right\},\left\{\chi_\alpha\right\}}
       \text{Tr}\left[\textbf{K}\textbf{H}^\text{EMFT}\right],
       \label{eq:emft_energy}

with respect to the NGWFs and elements of the density kernel
:math:`\textbf{K}`, using the conventional methods available in ONETEP.
The Hamiltonian is constructed as follows,

#. The total electron density :math:`n{\left(\mathbf{r}\right)}` is
   constructed from the full system NGWFs and kernel, from which
   :math:`V_\text{XC}^\text{low}{\left(\mathbf{r}\right)}` is
   calculated.

#. The active subsystem density
   :math:`n^\text{AA}{\left(\mathbf{r}\right)}` is constructed using the
   subsystem terms and the subsystem XC potentials
   :math:`V_\text{XC}^\text{low,A}{\left(\mathbf{r}\right)}` and
   :math:`V_\text{XC}^\text{high,A}{\left(\mathbf{r}\right)}`
   calculated.

#. Final EMFT potential can be written as

   .. math::

      V_\text{XC}^\text{high}{\left(\mathbf{r}\right)}
                  =V_\text{XC}^\text{low}{\left(\mathbf{r}\right)}+\left(V_\text{XC}^\text{high,A}{\left(\mathbf{r}\right)}-V_\text{XC}^\text{low,A}{\left(\mathbf{r}\right)}\right).

   with which we can construct the high-level Hamiltonian.

Although block orthogonalisation is found to work when just the density
kernel is being optimised, it does not do so generally for the
optimisation of the NGWFs. Because of this, there is an option to
optimise the NGWFs at the lower level of theory first (in all regions),
and then fix them for an optimisation of the kernel under EMFT.

If you would like to use hybrid functionals with embedding, there are
two things to bear in mind. Firstly, only hybrid-in-semi local DFT
calculations are currently supported – hybrid-in-hybrid calculations are
not possible. Secondly, the species in the ``species_swri-[swri name]``
block must match the species in the active region exactly. Anything else
will give incorrect results. Otherwise, the set-up of the hybrid
functional calculation is identical to a normal ONETEP calculation.

LR-TDDFT calculations can be performed with embedding
(TD-EMFT) [Ding2017-2]_, and this is also implemented
within ONETEP [Prentice2022]_. If you would like to
perform a TD-EMFT calculation, it may be advisable to restrict the
excitations to the active region, using the ``species_tddft_kernel``
block.

It is also possible to place the quantum embedding system within
implicit solvent, giving multi-level embedding
capability [Prentice2022]_. This should work very
similarly to standard implicit solvent calculations, although there are
a couple of additional keywords (see below). The main difference is
whether the cavity is constructed using the density kernel optimised
solely at the low level of theory, or optimised using EMFT. This only
makes a difference if the active region is close to the edge of the
cavity.

Keywords
========

-  ``species_ngwf_regions`` (block): This block defines which species
   are in which region. Each line of the block corresponds to a distinct
   region. The species within each region do not necessarily have to be
   physically next to one another. If this block is not defined, it is
   assumed that there is only one region, containing all the species in
   the system.

-  ``do_fandt`` (logical): Controls whether a freeze-and-thaw (F+T)
   optimisation of the NGWFs is performed or not. This is a cruder form
   of embedding, where all regions are treated at the same level of
   theory, but each region’s NGWFs are optimised in turn, with the
   others frozen. Default ``F``.

-  ``freeze_switch_steps`` (integer): How many NGWF CG optimisation
   steps should be spent on each region before moving onto the next in a
   F+T calculation. ``maxit_ngwf_cg`` represents the total number of
   NGWF optimisation steps across all regions. A value less than 0 means
   that all NGWFs are optimised together i.e. no F+T takes place.
   Default ``-1``.

-  ``use_emft`` (logical): Controls whether an EMFT calculation is
   performed, as described above. Default ``F``.

-  ``active_region`` (integer): Defines which region is the active
   region – 1 means the species on the first line in the
   ``species_ngwf_regions`` block constitute the active region, 2 means
   the second line, and so on. Default ``1``.

-  ``active_xc_functional`` (string): Defines what functional is used as
   the higher level of theory within EMFT. Default is the value of
   ``xc_functional`` i.e. no difference between the regions.

-  ``freeze_envir_ngwfs`` (logical): Controls whether the environment
   NGWFs should ever be optimised or not. Default ``F``.

-  ``use_emft_follow`` (logical): Controls whether the EMFT calculation
   is only performed after a regular calculation, so the NGWFs are
   optimised at the lower level of theory first, before applying EMFT.
   Default ``F``.

-  ``use_emft_lnv_only`` (logical): Controls whether only the kernel is
   optimised within EMFT, with the NGWFs optimised at the lower level of
   theory and then fixed. Usually used in conjunction with
   ``use_emft_follow``. Default ``F``.

-  ``emft_lnv_steps`` (integer): Controls the number of LNV kernel
   optimisation steps to be used in conjunction with
   ``use_emft_lnv_only``. Default ``10``.

-  ``block_orthogonalise`` (logical): Controls whether the environment
   NGWFs are orthogonalised with respect to the active region NGWFs, as
   described above. Default ``F``.

-  ``parallel_scheme`` (string): Defines the parallel scheme used for
   the calculation. See Appendix for more information. Default ``NONE``.

-  ``read_sub_denskern`` (logical): Controls whether only diagonal
   blocks of the density kernel are read in when restarting. This is
   useful for starting an embedding calculation from two separate
   calculations on the individual regions, so you only have the diagonal
   blocks of the density kernel. Default ``F``.

-  ``embed_debug`` (logical): Turns on verbose printing for debugging of
   embedding functionalities. Default ``F``.

-  ``is_restart_vac_from_vac`` (logical): Decides whether the vacuum
   calculation in an autosolvation implicit solvent calculation should
   be restarted from the vacuum\_ files or not. Useful for restarting
   autosolvation calculations if they time-out or similar. Default
   ``F``.

-  ``is_emft_cavity`` (logical): Decides whether the cavity used in
   implicit solvent calculations is determined using the low-level
   density kernel (``F``), or the EMFT-optimised density kernel (``T``).
   Default ``F``.

The most reliable way to run EMFT calculations is to have ``use_emft``,
``use_emft_follow``, ``use_emft_lnv_only`` and ``block_orthogonalise``
all set to ``T``. These can be set to ``F`` (most sensibly in reverse
order i.e. ``block_orthogonalise`` first), but the calculation may
become more unstable, depending on the system, the regions chosen and
the functionals chosen.

Example input file
==================

::

    !====================================================!
    ! Input for calculation with the ONETEP program      !
    !                                                    !
    ! O2 and H2 form the embedded system to be treated   !
    ! at the higher level of theory, O1 and H1 are the   !
    ! environment treated at the low-level.              !
    !====================================================!

    %block species_ngwf_regions
    O2 H2
    O1 H1
    %endblock species_ngwf_regions

    task: SINGLEPOINT
    cutoff_energy 1000 eV
    write_forces: T
    xc_functional: LDA
    active_xc_functional: PBE

    use_emft: T
    use_emft_follow: T
    use_emft_lnv_only: T
    block_orthogonalise : T
    parallel_scheme: HOUSE

    %block species_atomic_set
    H1 "SOLVE"
    O1 "SOLVE"
    H2 "SOLVE"
    O2 "SOLVE"
    %endblock species_atomic_set

    %block species
    H1 H 1 1 7.0
    O1 O 8 4 7.0
    H2 H 1 1 7.0
    O2 O 8 4 7.0
    %endblock species

    %block species_pot
    H1 "pseudo/hydrogen.recpot"
    O1 "pseudo/oxygen.recpot"
    H2 "pseudo/hydrogen.recpot"
    O2 "pseudo/oxygen.recpot"
    %endblock species_pot

    %block lattice_cart
         30.000000000       0.000000000       0.000000000
          0.000000000      30.000000000       0.000000000
          0.000000000       0.000000000      30.000000000
    %endblock lattice_cart

    %block positions_abs
    O1       16.203224001     15.100000000     11.536063353
    H1       15.100000000     15.100000000     10.100000000
    H1       15.100000000     15.100000000     12.991451046
    O2       12.600158789     15.100000000     17.306583960
    H2       13.051873252     13.656398529     18.308114239
    H2       13.051873252     16.543601471     18.308114239
    %endblock positions_abs

Interaction with other functionalities
======================================

Fully tested
------------

-  Energy and forces calculations

-  Hybrid-in-semi local DFT

-  Restarting calculations

-  LR-TDDFT

-  Implicit solvent

Should work, not thoroughly tested
----------------------------------

-  Geometry optimisation

-  Finite displacement phonons

-  Molecular dynamics

-  Conduction NGWF optimisation

-  Ensemble DFT

-  Kernel DIIS

-  QNTO

-  NAO

-  Cutoff Coulomb

-  Spin polarised calculations

-  Some properties calculations (eigenstates, Mulliken charges,
   plotting, DoS)

Not compatible with embedding
-----------------------------

-  Hubbard calculations

-  DMFT

-  PAW

-  cDFT

-  Bandstructure calculations

-  DMA

-  EDA

-  Electronic transport

-  Hybrid-in-hybrid DFT

-  NEB

-  EELS

-  Polarisable embedding

-  Transition state searching

-  DDEC

Any functionalities missed above are likely to not work with embedding.

Appendix: Parallel strategies with embedding
============================================

In a normal ONETEP calculation, atoms are distributed across the
available MPI processes according to a ‘parallel strategy’. This
determines how resources such as matrix elements will be spread across
the MPI environment in order to reduce the communication between nodes
and maximise the efficiency of the calculation. Details on maximising
parallel efficiency are available via the ONETEP documentation and
website.

As part of the embedding infrastructure, each subsystem is given its own
parallel strategy. This contains all information relating to the
distribution of resources across the MPI nodes available to the
calculation, which are determined by the parameter ``PARALLEL_SCHEME``.
There are three settings for the distribution of resources during an
embedding calculation:

-  ``NONE``: All subsystems are treated completely independently, with
   atoms distributed across all available processors as though the other
   subsystems do not exist. The number of MPI processes cannot be
   greater than the number of atoms in the smallest subsystem. For
   example, if there are 8 processors available then each will hold
   atoms and data from all subsystems, though the calculation will fail
   if any subsystem has less than 8 atoms (or possibly slightly more if
   the space-filling curve is in use). This is the default setting for
   testing but is not recommended for practical calculations due to the
   constraint on the number of processors.

-  ``SENATE``: Nodes are partitioned evenly between all subsystems. For
   example, if there are 8 processors and 2 subsystems, then each will
   be allocated 4 processors, regardless of the number of atoms in each
   subsystem. Unlike the ``NONE`` setting, there is no upper bound on
   the number of processors which may be used, so user discretion is
   advised.

-  ``HOUSE``: Divides the processors proportionally between all
   subsystems, with a minimum of 1 processor per subsystem. For example,
   if we have two subsystems consisting of 15 and 5 atoms each, then
   with 8 processors each subsystem will be allocated 6 and 2 nodes
   respectively. At a minimum all subsystems are granted 1 processor —
   if we had two subsystems with 1 and 100 atoms in our 8 processor
   example, then they will receive 1 and 7 processor respectively. Like
   ``SENATE``, there is no upper bound on the number of processors that
   can be allocated and finding a sensible setting is left to the user.

``HOUSE`` is the recommended setting for running calculations, the
others are mainly of use for testing. Since they should all produce the
same results, any significant differences may be a sign of an underlying
problem, so comparing them is a useful consistency check.

[Ding2017-2] F. Ding, T. Tsuchiya, F. R. Manby and T. F. Miller, *J. Chem. Theory Comput.*, **13**, 4216–4227, (2017).

[Prentice2020] J. C. A. Prentice, R. J. Charlton, A. A. Mostofi and P. D. Haynes, *J. Chem. Theory Comput.*, **16**, 354–365, (2020).

[Prentice2022] J. C. A. Prentice, *J. Chem. Theory Comput.*, **18**, 1542-1554 (2022).

[Huang2008] P. Huang and E. M. Carter, *Annu. Rev. Phys. Chem.*, **59**, 261–290, (2008).

[Gomes2012] A. S. P. Gomes and C. R. Jacob, *Annu. Rep. Prog. Chem., Sect. C: Phys. Chem.*, **108**, 222–277, (2012).

[Fornace2015] M. E. Fornace, J. Lee, M. Kaito, F. R. Manby, T. F. Miller, *J. Chem. Theory Comput.*, **11**, 568–580, (2015).

[Ding2017] F. Ding, F. R. Manby and T. F. Miller, *J. Chem. Theory Comput.*, **13**, 1605–1615, (2017).
