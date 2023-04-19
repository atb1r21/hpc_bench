============================================
Constrained Density Functional Theory (cDFT)
============================================

:Author: Gilberto Teobaldi, University of Liverpool (g.teobaldi@liv.ac.uk)

:Date:   November 2013

cDFT input check-list
=====================

This is a short check-list meant to help the setting up of a
constrained-DFT (cDFT) simulation. ONETEP implements some rather
extensive check of the input (as specified in .dat file). In spite of
this, users may be still capable of creating erroneous inputs which I
could not think of. If you experience disaster (commutators larger than
1.E-3, NWGFs- and cDFT-optimisation stuck at the same value of the
gradient for many iterations, or NaN (Not A Number) and ‘\*\*\*’ outputs
from ONETEP), please do get in touch.

Can I start by running a cDFT geometry-optimisation and/or Molecular Dynamics straightaway?
-------------------------------------------------------------------------------------------

Currently, this is **not** possible. ONETEP will not allow you to run a
geometry-optimisation/molecular dynamics simulation unless a .cdft file
(containing the cDFT-potentials, Uq/s) is present. This is meant to
force you to first make sure that your cDFT-settings let you obtain a
good convergence in a single-point calculation for the system and
constrained solution you are interested in.

What projectors am I telling ONETEP to use?
-------------------------------------------

Is cdft\_read\_projector=F (this is the default) and cdft\_multi\_proj=F
(this is the default) in your input?

If so, you will be using as cDFT-projectors the orbitals, of angular
momentum *L-projectors*, of a hydrogenic atom of atomic number
*Z-projector*, provided *Z-projector* is positive in the %block
constrained\_dft (Z>0). Conversely, if you have entered a negative
*Z-projector* value in the %block constrained\_dft (Z<0), you will be
using the numerical orbital (of angular momentum *L-projectors*) for the
given cDFT-atom as cDFT-projectors.

Is cdft\_multi\_proj=T (this is NOT the default) in your input?

If so, you will be using as cDFT-projectors the numerical orbitals,
obtained from the pseudo-atom solver. You will be using as many
cDFT-projectors as NGWFs for the given cDFT-site. Thus, to know the
number and angular momentum of your cDFT-projectors, you need to refer
to the specific settings of the cDFT-site in the %block
species\_atom\_set and the relevant output in the standard .out file.

Is cdft\_read\_projector=T in your input?

If so, you will use the projectors stored in the .tightbox\_hub\_projs
file (or experience a crash if this file is not present). Mind that the
file will contain different classes of projectors (hydrogenic orbitals,
numerical orbitals, self-consistent projectors) depending on the
settings you used for the cDFT-run which generated the
.tightbox\_hub\_projs file. This means that, unless one wants to
experience a crash, the projectors in the .tightbox\_hub\_projs file
should have been generated with the same entries for *Z-projector* *and
L-projector* in %block constrained\_dft **OR** cdft\_multi\_proj and
%block species\_atom\_set as in the current calculation.

Is this the first single-point cDFT submission or do I want to restart a single-point calculation?
--------------------------------------------------------------------------------------------------

For submission from scratch the user can chose between starting the
cDFT-run from previously obtained DFT-NGWFs and Density Kernel (DKN), by
setting (read\_denskern=T and read\_tightbox\_ngwfs=T), **or not**
(read\_denskern=F and read\_tightbox\_ngwfs=F).

**To restart a cDFT**-single point run (from the latest Uq/s-potentials)
or to run cDFT-geometry optimisations/Molecular Dynamics it is necessary
to activate cdft\_continuation=T. **Mind** that if, instead of
hydrogenic or numerical orbitals as cDFT-projectors (see above) one
wants to use/keep using pre-optimised projectors (stored in
.tightbox\_hub\_projs file) it is necessary to activate also
cdft\_read\_projector=T in a restart.

Are the columns of the %block constrained\_dft requesting the right targets for the given cDFT-flavour?
-------------------------------------------------------------------------------------------------------

**MIND** that, as explained in the description of the keywords,
different columns of the %block constrained\_dft are used depending on
the selected cDFT-modes (cdft\_atom\_charge, cdft\_atom\_spin,
cdft\_group\_charge/spin\_acceptor/donor,
cdft\_group\_charge/spin\_diff). Are you requesting the right targets in
terms of atomic-population, group-populations, atomic magnetic-moments
and group magnetic-moments?

Mind also that the cdft\_group\_charge/spin\_acceptor/donor,
cdft\_group\_charge/spin\_diff cDFT-modes require that targeted
population, magnetic-moment and differences are entered explicitly with
the corresponding \_TARGET keyword. If you have forgotten entering the
relevant \_TARGET keyword in your input, the simulation will stop and
ONETEP will tell you about it.

Have I set CDFT\_GROUP\_CHARGE\_UP/DOWN\_ONLY=T?
------------------------------------------------

Remember that, for CDFT\_GROUP\_CHARGE\_UP/DOWN\_ONLY=T, only the
corresponding spin-channel will be constrained in
cdft\_group\_charge\_acceptor/donor and cdft\_group\_charge\_diff runs.
Accordingly, you should target populations for one-spin channel only
(**not** for the UP+DOWN channels).

Is the sign of U\ :sub:`q/s` correct for the group\_charge/spin\_acceptor/donor/diff run?
-----------------------------------------------------------------------------------------

Remember that atoms in acceptor- and donor-group are identified by mean
of the sign of the Uq/s in the %block constrained\_dft. Thus,
**acceptor-group atoms** need to be assigned negative (e-attractive)
U\ :sub:`q/s` (**U\ :sub:`q/s`\ <0**), whereas **donor-group atoms**
need to be assigned positive (e-repulsive) U\ :sub:`q/s`
(**U\ :sub:`q/s`>0**).

Are my constraints compatible with the spin (=N\ :sub:`UP`\ -N\ :sub:`DOWN`\ ) keyword in the input (.dat) file?
----------------------------------------------------------------------------------------------------------------

Remember that ONETEP optimises the Density Kernel keeping the number of
UP (N\ :sub:`UP`\ ) and DOWN (N\ :sub:`DOWN`\ ) electrons fixed. As a result,
the net magnetization of the system (N\ :sub:`UP`\ -N\ :sub:`DOWN`\ ) is also
fixed. Have you introduced any incongruence between the spin=0,1,2,etc
keyword and the cDFT-ones in your input? If your targeted cDFT-solution
results in a high-value (i.e. > 0 μ\ :sub:`B`) magnetization of the
**total** system, the spin keyword should reflect it.

What values of min/maxit\_lnv, maxit\_ngwf\_cg and maxit\_cdft\_u\_cg should I use?
-----------------------------------------------------------------------------------

Based on the tests run so far, and considering that the cDFT-potentials
(Uq/s) are iteratively optimised at each (1:sup:`st`) step of the
NGWFs-optimisation line-search, the following choices seem to be
reasonable:

::

   maxit_palser_mano : -1
   kerfix : 2
   maxit_pen : 0
   minit_lnv : 5
   maxit_lnv : 10
   maxit_ngwf_cg : 60
   lnv_check_trial_steps : T
   lnv_threshold_orig : 1.0e-10
   maxit_cdft_u_cg : 5

Increasing maxit\_cdft\_u\_cg above 5 is hardly going to accelerate the
convergence of the cDFT-run (the NGWFs will change at the following step
of NGWFs-optimisation), decreasing minit\_lnv below 5 might be risky.

To what value should I initialise the cDFT-potentials (Uq/s)?
-------------------------------------------------------------

In the current implementation, unless cdft\_guru=T, in which case the
Uq/s entered in the %block constrained\_dft will be used, the absolute
value of the Uq/s cDFT-potentials are internally initialised to 1 eV
(their original sign is, of course, maintained). For spin-excitation
\|Us\|=1eV may be too large, and initialising the \|Us\| with 0.1-0.3 eV
may accelerate convergence (this has been tested only on triplet
excitations in benzene dimers).

What is the difference between a cdft\_atom\_charge run and a cdft\_group\_charge\_acceptor/donor one with one-atom group?
--------------------------------------------------------------------------------------------------------------------------

Whereas the cdft\_atom\_charge=T mode allows independent (i.e.
potentially different) constraining potentials to be applied to the UP
and DOWN spin-channels, for **one-atom**
cdft\_group\_charge\_acceptor/donor=T the same Uq will be applied to
both the UP and DOWN spin-channel. Activation of
cdft\_group\_charge\_up/down\_only=T in cdft\_group\_acceptor/donor
modes allows to optimise Uq for only a spin-channel, leaving the other
spin-channel unconstrained.

Have I chosen a meaningful cdft\_cg\_max for the cDFT-mode I wish to use?
-------------------------------------------------------------------------

For cDFT-runs with only one cDFT-group in the system
(cdft\_group\_charge/spin\_acceptor/donor modes) and
group\_charge/spin\_diff runs, **only one** cDFT-potential (Uq/s) will
be optimised in the cDFT-loop. Accordingly, for these cases it is
recommended to perform the cDFT-optimisation via a steepest descendent
algorithm (cdft\_cg\_max=1).

How do I obtain the population of the cDFT-sites for a standard DFT-run?
------------------------------------------------------------------------

Performing a single-point (task=singlepoint) fixed-Uq/s
(maxit\_cdft\_u\_cg=0) cDFT-run using very small (e.g. 1.E-60)
Uq/s-potentials in the %block constrained\_dft and setting
output\_detail=VERBOSE will result in the cDFT-population of all the
cDFT-sites being printed in the standard output (.out) file. To obtain
atom-specific (instead of atomic\_species-specific) information on the
population of the cDFT-sites, it is necessary to set
CDFT\_PRINT\_ALL\_OCC=T.

**Mind.** The population of a given cDFT-site depends critically on the
projector used. Make sure you decide your cDFT-targets from the
DFT-populations obtained with the same set of projectors!

How do I optimise self-consistently the projectors for a given geometry?
------------------------------------------------------------------------

In analogy with DFT+U simulations, self-consistent optimisation of the
cDFT-projectors (for a fixed geometry) is activated by setting
task=HUBBARDSCF in the .dat file. It is recommended to start the
task=HUBBARDSCF run from pre-optimised PAO-cDFT projectors (Z<0 in
%block constrained\_dft), cDFT-potentials (Uq/s), NGWFs and DKN. This is
accomplished, regardless of the keywords specific to the chosen
cDFT-flavour, by making sure the input file (.dat) contains:

::

   task : HUBBARDSCF
   hubbard_max_iter : 40 # perform 40 Hubbard-SCF iterations
   read_denskern : T
   read_tightbox_ngwfs : T
   cdft_read_projectors : T
   cdft_continuation : T

The percentage of the latest-optimised cDFT-NGWFs to be used as new
cDFT-projectors is controlled by the hubbard\_proj\_mixing keyword [0 <
\|hubbard\_proj\_mixing\| < 1], with hubbard\_proj\_mixing=1 meaning
that the latest optimised-NGWFs are used entirely as new
cDFT-projectors.

**Mind.** Provided cdft\_read\_projector=T, the cDFT-projectors in the
.tightbox\_hub\_projs file are used (entirely) as new projector during
the 1\ :sup:`st` HUBBARDSCF iteration.

Is self-consistent optimisation of the cDFT-projectors at the DFT-geometry a good idea?
---------------------------------------------------------------------------------------

From the tests so far, this is **not a good idea**. For tightly
constrained systems (for instance an hypothetical
N\ :sup:`(+)`\ =N\ :sup:`(-)` excitation) a cDFT task=HUBBARDSCF run at
the DFT-optimised geometry may result in very slow convergence. The
recommended procedure is to **first** optimise the cDFT-geometry using
numerical orbitals (PAO) as projectors (Z-projector <0 in %block
constrained\_dft) and **then** optimise the cDFT-projectors at the
PAO-cDFT optimised geometry.

cDFT Keywords
=============

Intermediate Keywords
---------------------

**CDFT\_ATOM\_CHARGE**

Syntax: CDFT\_ATOM\_CHARGE [Logical]

Description: Activate atom charge-constrained-DFT mode. This mode is
incompatible with any other cDFT-mode.

Default: False

Example: CDFT\_ATOM\_CHARGE T

**CDFT\_ATOM\_SPIN**

Syntax: CDFT\_ATOM\_SPIN [Logical]

Description: Activate atom magnetic-moment-constrained-DFT mode. This
mode is incompatible with any other cDFT-mode.

Default: False

Example: CDFT\_ATOM\_SPIN T

**CDFT\_CG\_MAX**

Syntax: CDFT\_CG\_MAX [Real]

Description: Specifies the maximum number of constraining potential
(Uq/s) conjugate gradient iterations between resets.

Default: Number of independent Uq/s for cdft\_guru=F

Example: CDFT\_CG\_MAX 1 #Perform steepest descents optimisation

**CDFT\_CG\_MAX\_STEP**

Syntax: CDFT\_CG\_MAX\_STEP [Real]

Description: Maximum length of trial step for the constraining potential
(Uq/s) optimisation line search.

Default: 50.0

Example: CDFT\_CG\_MAX\_STEP 10.0

**CDFT\_CG\_THRESHOLD**

Syntax: CDFT\_CG\_THRESHOLD [Real]

Description: Specifies the convergence threshold for the RMS gradient of
the constraining potentials (Uq/s).

Default: 1.0E-3

Example: CDFT\_CG\_THRESHOLD 0.01

**CDFT\_CG\_TYPE**

Syntax: CDFT\_CG\_TYPE [Text]

Description: Specifies the variant of the conjugate gradients algorithm
used for the optimization of the constraining potentials (Uq/s),
currently either NGWF\_FLETCHER for Fletcher-Reeves or NGWF\_POLAK for
Polak-Ribiere.

Default: NGWF\_FLETCHER

Example: CDFT\_CG\_TYPE NGWF\_POLAK

**CDFT\_CHARGE\_ACCEPTOR\_TARGET**

Syntax: CDFT\_CHARGE\_ACCEPTOR\_TARGET [Real]

Description: Targeted acceptor-group electron population for
acceptor-group charge-constrained-DFT mode
[CDFT\_GROUP\_CHARGE\_ACCEPTOR=T].

Default: 0.

Example: CDFT\_CHARGE\_ACCEPTOR\_TARGET 17 #Constrain Nup+Ndown=17 e in
subspace

**CDFT\_CHARGE\_DONOR\_TARGET**

Syntax: CDFT\_CHARGE\_DONOR\_TARGET [Real]

Description: Targeted donor-group electron population for donor-group
charge-constrained-DFT mode [CDFT\_GROUP\_CHARGE\_DONOR=T].

Default: 0.

Example: CDFT\_CHARGE\_DONOR\_TARGET 17 #Constrain Nup+Ndown=17 e in
subspace

**CDFT\_CONTINUATION**

Syntax: CDFT\_CONTINUATION [Logical]

Description: Continue a constraining potential (Uq/s) optimisation from
a previous run using the .cdft file with the latest cDFT-potentials.
CDFT\_CONTINUATION=T allows also to perform single-point cDFT runs
(MAXIT\_CDFT\_U\_CG=0) reading atom-specific constraining potentials
from .cdft file (instead of species-specific ones from the
CONSTRAINED\_DFT block). For cdft\_continuation=T, the constraining
potentials (Uq/s) are read from the .cdft file no matter the setting of
cdft\_guru.

Default: False

Example: CDFT\_CONTINUATION T

**CDFT\_GROUP\_CHARGE\_ACCEPTOR**

Syntax: CDFT\_GROUP\_CHARGE\_ACCEPTOR [Logical]

Description: Activate acceptor-group charge-constrained-DFT mode. This
mode is compatible with CDFT\_GROUP\_CHARGE\_DONOR and
CDFT\_GROUP\_SPIN\_ACCEPTOR/DONOR cDFT-modes, and incompatible with
CDFT\_ATOM\_CHARGE/SPIN and CDFT\_GROUP\_CHARGE/SPIN\_DIFF cDFT modes.

Default: False

Example: CDFT\_GROUP\_CHARGE\_ACCEPTOR T

**CDFT\_GROUP\_CHARGE\_DIFF**

Syntax: CDFT\_GROUP\_CHARGE\_DIFF [Logical]

Description: Activate group charge-difference constrained-DFT mode. This
mode is compatible with CDFT\_GROUP\_SPIN\_DIFF cDFT mode only. Thus, it
is incompatible with any other CDFT\_ATOM\_CHARGE/SPIN and
CDFT\_GROUP\_CHARGE/SPIN\_ACCEPTOR/DONOR cDFT modes.

Default: False

Example: CDFT\_GROUP\_CHARGE\_DIFF T

**CDFT\_GROUP\_CHARGE\_DIFF\_TARGET**

Syntax: CDFT\_CHARGE\_DIFF\_TARGET [Real]

Description: Targeted electron population difference between acceptor
and donor group for -group charge-difference constrained-DFT mode
[CDFT\_GROUP\_CHARGE\_DIFF=T].

Default: 0.

Example: CDFT\_CHARGE\_ACCEPTOR\_TARGET 2

    #Constrain [Nup+Ndown]\_ACC - [Nup+Ndown]\_DON to 2 e.

**CDFT\_GROUP\_CHARGE\_DONOR**

Syntax: CDFT\_GROUP\_CHARGE\_DONOR [Logical]

Description: Activate donor-group charge-constrained-DFT mode. This mode
is compatible with CDFT\_GROUP\_CHARGE\_ACCEPTOR and
CDFT\_GROUP\_SPIN\_ACCEPTOR/DONOR cDFT-modes, and incompatible with
CDFT\_ATOM\_CHARGE/SPIN and CDFT\_GROUP\_CHARGE/SPIN\_DIFF cDFT modes.

Default: False

Example: CDFT\_GROUP\_CHARGE\_DONOR T

**CDFT\_GROUP\_CHARGE\_DOWN\_ONLY**

Syntax: CDFT\_GROUP\_CHARGE\_DOWN\_ONLY [Logical]

Description: Constrain only SPIN-DOWN channel in
CDFT\_GROUP\_CHARGE\_ACCEPTOR, CDFT\_GROUP\_CHARGE\_DONOR and
CDFT\_GROUP\_CHARGE\_DIFF modes. To avoid disaster, make sure the
specified CDFT\_CHARGE\_ACCEPTOR/DONOR\_TARGET or
CDFT\_CHARGE\_DIFF\_TARGET keywords are consistent with the fact only
one spin channel is being constrained. This functionality is NOT
compatible with CDFT\_GROUP\_CHARGE\_UP\_ONLY, CDFT\_ATOM\_CHARGE/SPIN,
and CDFT\_GROUP\_SPIN\_ACCEPTOR/DONOR and CDFT\_GROUP\_SPIN\_DIFF cDFT
modes.

Default: False

Example: CDFT\_GROUP\_CHARGE\_DOWN\_ONLY T

**CDFT\_GROUP\_CHARGE\_UP\_ONLY**

Syntax: CDFT\_GROUP\_CHARGE\_UP\_ONLY [Logical]

Description: Constrain only SPIN-UP channel in
CDFT\_GROUP\_CHARGE\_ACCEPTOR, CDFT\_GROUP\_CHARGE\_DONOR and
CDFT\_GROUP\_CHARGE\_DIFF modes. To avoid disaster, make sure the
specified CDFT\_CHARGE\_ACCEPTOR/DONOR\_TARGET or
CDFT\_CHARGE\_DIFF\_TARGET keywords are consistent with the fact only
one spin channel is being constrained. This functionality is NOT
compatible with CDFT\_GROUP\_CHARGE\_DOWN\_ONLY,
CDFT\_ATOM\_CHARGE/SPIN, and CDFT\_GROUP\_SPIN\_ACCEPTOR/DONOR and
CDFT\_GROUP\_SPIN\_DIFF cDFT modes.

Default: False

Example: CDFT\_GROUP\_CHARGE\_UP\_ONLY T

**CDFT\_GROUP\_SPIN\_ACCEPTOR**\ Syntax: CDFT\_GROUP\_SPIN\_ACCEPTOR
[Logical]

Description: Activate acceptor-group magnetic-moment constrained-DFT
mode. This mode is compatible with CDFT\_GROUP\_SPIN\_DONOR and
CDFT\_GROUP\_CHARGE\_ACCEPTOR/DONOR cDFT-modes, and incompatible with
CDFT\_ATOM\_CHARGE/SPIN and CDFT\_GROUP\_CHARGE/SPIN\_DIFF cDFT modes.

Default: False

Example: CDFT\_GROUP\_SPIN\_ACCEPTOR T

**CDFT\_GROUP\_SPIN\_DIFF**

Syntax: CDFT\_GROUP\_SPIN\_DIFF [Logical]

Description: Activate group magnetic-moment-difference constrained-DFT
mode. This mode is compatible with CDFT\_GROUP\_CHARGE\_DIFF cDFT mode
only. Thus, it is incompatible with any other CDFT\_ATOM\_CHARGE/SPIN
and CDFT\_GROUP\_CHARGE/SPIN\_ACCEPTOR/DONOR cDFT modes.

Default: False

Example: CDFT\_GROUP\_CHARGE\_DIFF T

**CDFT\_GROUP\_SPIN\_DIFF\_TARGET**

Syntax: CDFT\_SPIN\_DIFF\_TARGET [Real]

Description: Targeted magnetic-moment difference between acceptor and
donor group for group magnetic-moment-difference constrained-DFT mode
[CDFT\_GROUP\_SPIN\_DIFF=T].

Default: 0.

Example: CDFT\_CHARGE\_ACCEPTOR\_TARGET 2

    #Constrain [Nup-Ndown]\_ACC - [Nup-Ndown]\_DON to 2 e.

**CDFT\_GROUP\_SPIN\_DONOR**

Syntax: CDFT\_GROUP\_SPIN\_DONOR [Logical]

Description: Activate donor-group magnetic-moment constrained-DFT mode.
This mode is compatible with CDFT\_GROUP\_SPIN\_ACCEPTOR and
CDFT\_GROUP\_CHARGE\_ACCEPTOR/DONOR cDFT-modes, and incompatible with
CDFT\_ATOM\_CHARGE/SPIN and CDFT\_GROUP\_CHARGE/SPIN\_DIFF cDFT modes.

Default: False

Example: CDFT\_GROUP\_SPIN\_DONOR T

**CDFT\_GURU**

Syntax: CDFT\_GURU [Logical]

Description: Tell ONETEP you are a cDFT-expert and prevent it from
initialising the active \|Uq/s\| to failsafe value of 1 eV overwriting
the values entered in the %block constrained\_dft (Uq/s).

Default: False

Example: CDFT\_GURU T

**CDFT\_HUBBARD**

Syntax: CDFT\_HUBBARD [Logical]

Description: Activate the constrained-DFT+U functionality. It requires
specifications of a positive value for the Hubbard correction (Uh) in
the CONSTRAINED\_DFT Block.

Default: False

Example: CDFT\_HUBBARD T

**CDFT\_MAX\_GRAD**

Syntax: CDFT\_MAX\_GRAD [Real]

Description: Specifies the convergence threshold for the maximum value
of the constraining-potential (Uq/s) gradient at any cDFT-site

Default: 1.0E-3

Example: CDFT\_MAX\_GRAD 0.01

**CDFT\_MULTI\_PROJ**

Syntax: CDFT\_MULTI\_PROJ [Logical]

Description: Activate the “as many cDFT-projectors as NGWFs” cDFT-mode.
In this mode, the number of cDFT-projectors for a given cDFT-atom equals
the number of NWGFs for that atom as specified in the %block species.
Both the cDFT-projectors and the NGWFs are localised within spheres of
the same radius. When activated, this mode overwrites the L-projectors
and Z-projectors settings in %block constrained\_dft, and the
cDFT-projectors are built according to the settings in %block
species\_atomic\_set for that atom=cDFT-site.

Default: False

Example: CDFT\_MULTI\_PROJ T

**CDFT\_PRINT\_ALL\_OCC**

Syntax: CDFT\_PRINT\_ALL\_OCC [Logical]

Description: Print detailed information of occupancies for al the
cDFT-sites, for OUTPUT\_DETAIL = VERBOSE.

Default: False

Example: CDFT\_PRINT\_ALL\_OCC T

**CDFT\_READ\_PROJ**

Syntax: CDFT\_READ\_PROJ [Logical]

Description: Read cDFT-projectors from .tightbox\_hub\_proj file.
Activation of this keyword overwrites any Z-projector setting in %block
constrained\_dft. It also makes not necessary to set
hubbard\_proj\_mixing<0 to have task=HUBBARDSCF runs with projectors
read in from file.

Default: False

Example: CDFT\_READ\_PROJ T

**CDFT\_SPIN\_ACCEPTOR\_TARGET**

Syntax: CDFT\_SPIN\_ACCEPTOR\_TARGET [Real]

Description: Targeted group magnetic-moment for acceptor-group
magnetic-moment constrained-DFT mode [CDFT\_GROUP\_SPIN\_ACCEPTOR=T].

Default: 0.

Example: CDFT\_SPIN\_ACCEPTOR\_TARGET -2 #Constrain Nup-Ndown=-2 in
subspace

**CDFT\_SPIN\_DONOR\_TARGET**

Syntax: CDFT\_SPIN\_DONOR\_TARGET [Real]

Description: Targeted group magnetic-moment for donor-group
magnetic-moment constrained-DFT mode [CDFT\_GROUP\_SPIN\_DONOR=T].

Default: 0.

Example: CDFT\_SPIN\_DONOR\_TARGET -2 #Constrain Nup-Ndown=-2 in
subspace

**CDFT\_TRIAL\_LENGTH**

Syntax: CDFT\_TRIAl\_LENGTH [Real]

Description: Specifies initial trial length for first step of
constraining-potential (Uq/s) conjugate gradients optimisation.

Default: 0.1

Example: CDFT\_TRIAL\_LENGTH 1.0

**CI\_CDFT**

Syntax: CI\_CDFT [Logical]

Description: Perform a Configuration Interaction calculation based on
constrained-DFT configurations.

Default: False

Example: CI\_CDFT T

**CI\_CDFT\_NUM\_CONF**

Syntax: CDFT\_MAX\_GRAD [Integer]

Description: Specifies the number of constrained-DFT configuration
available for a CI\_CDFT=T simulation

Default: 0

Example: CI\_CDFT\_NUM\_CONF 4

**CONSTRAINED\_DFT**

Syntax: CONSTRAINED\_DFT [Block]

Syntax: %BLOCK CONSTRAINED\_DFT

S1 L1 Z1 Uh1 Uq1(UP) Uq1(DOWN) Us1 N1(UP) N1(DOWN) [N1(UP)-N1(DOWN)] S2
L2 Z2 Uh2 Uq2(UP) Uq2(DOWN) Us2 N2(UP) N2(DOWN) [N2(UP)-N2(DOWN)] . . .
. .

. . . . .

SM LM ZM UhM UqM(UP) UqM(DOWN) UsM NM(UP) NM(DOWN) [NM(UP)-NM(DOWN)]

%ENDBLOCK CONSTRAINED\_DFT

Description: Manages constrained-DFT simulations. Provided
cdft\_multi\_proj=F, for species S and subspace of angular momentum
channel L (with principal quantum number n=L+1) we apply charge
spin-specific [Uq(UP), Uq(DOWN)] or magnetic-moment-specific (Us)
constraining potentials (eV). For cdft\_atom\_charge=T, N(UP) and
N(DOWN) indicate the targeted e-population for spin-channel UP and DOWN,
respectively. For cdft\_atom\_spin=T, [N1(UP)-N1(DOWN)] indicates the
targeted e-population difference (i.e. local magnetic moment). Uh
indicates the optional Hubbard parameter (U, eV) to be applied for
cdft\_hubbard=T. An effective nuclear charge Z defines the hydrogenic
orbitals spanning the subspace unless a negative value is given, e.g.,
Z=-10, in which case the NGWFs initial guess orbitals (numerical atomic
orbitals) are used. Depending on the activated cDFT-mode, different
columns of the block are used. These are:

    S, L, Z, (Uh), Uq(UP), Uq(DOWN), N(UP), N(DOWN) for
    cdft\_atom\_charge=T

    S, L, Z, (Uh), Us, [N(UP)-N(DOWN)] for cdft\_atom\_spin=T

    S, L, Z, (Uh), Uq(UP), Uq(DOWN) for cdft\_group\_charge\_acceptor=T,
    cdft\_group\_charge\_donor=T, or cdft\_group\_charge\_diff=T. In
    this case, Uq(UP) must be equal to Uq(DOWN). Acceptor and donor
    atoms are differentiated by mean of negative [Uq(UP/DOWN)<0] and
    positive [Uq(UP/DOWN)>0] constraining-potentials, respectively.
    Setting Uq=0 in the %block constrained\_dft will result in the given
    cDFT-atom being excluded from the list of the atoms in a given
    cdft\_group\_charge\_donor/acceptor/diff group.

    S, L, Z, (Uh), and Us for cdft\_group\_spin\_acceptor=T,
    cdft\_group\_spin\_donor=T, or cdft\_group\_spin\_diff=T. In this
    case, Acceptor and donor atoms are differentiated by mean of
    negative (Us<0) and positive (Us>0) constraining-potentials,
    respectively. Setting Us=0 in the %block constrained\_dft will
    result in the given cDFT-atom being excluded from the list of the
    atoms in a given cdft\_group\_spin\_donor/acceptor/diff group.

    cdft\_group\_spin\_acceptor=T, cdft\_group\_spin\_donor=T,
    cdft\_group\_charge\_acceptor=T and cdft\_group\_charge\_donor=T are
    all compatible one with another. Charge- and magnetic-moment
    acceptor- and donor-groups may or may not be the same group. Thus,
    besides simultaneously constraining the charge and magnetic-moment
    on a given group, it is also possible (by setting the appropriate
    sign of Uq and Us in the %block constrained\_dft) to create, within
    the same input and system, a charge\_acceptor group-A, a
    charge\_donor group-B, a spin\_acceptor group-C and a spin\_donor
    group-C. Similar considerations apply also for simultaneous
    activation of group\_charge\_diff and group\_spin\_diff cDFT-modes.
    In sum,

    Activation of cdft\_group\_charge\_up(down)\_only=T for
    cdft\_group\_charge\_acceptor/donor or cdft\_group\_charge\_diff
    modes leads to optimisation of the Uq potentials only for the
    selected spin-channel i.e. Uq(UP) only for
    cdft\_group\_charge\_up\_only=T, and Uq(DOWN) only for
    cdft\_group\_charge\_down\_only=T, leaving the other spin channel
    unconstrained.

    For cdft\_multi\_proj=T the L-projector and Z-projector columns in
    the %block constrained\_dft are read but NOT used. The
    cDFT-projectors are set on the basis of the %block
    species\_atomic\_set and taken as the NGWFs initial guess (numerical
    atomic orbital). This leads to as many cDFT-projectors as NGWFs for
    the cDFT-atom being used. In the current implementation, the same
    Uq/s is applied to all the projectors of a given cDFT-atom
    regardless of their principal quantum number and angular momentum.

    For all the cDFT-modes, unless maxit\_cdft\_u\_cg=0, and depending
    of the specific cDFT-mode, the constraining potentials (Uq,Us) will
    be automatically optimised. Note that, unless cdft\_guru=T, the
    constraining potentials (Uq/s) will be initialised to 1 eV. Thus, to
    perform fixed-Uq/s cDFT-runs or to initialise Uq/s with values
    different from 1 eV (useful for low-energy spin-excitation), it is
    necessary to set cdft\_guru=T.

    The CONSTRAINED\_DFT Block is incompatible with the HUBBARD Block.
    To perform a constrained-DFT+U simulation with Hubbard (Uh)
    correction applied to the subspace in addition to the constraining
    potentials (Uq/s) it is necessary to set cdft\_hubbard=T. For
    cdft\_hubbard=F (which is the default), the Hubbard correction will
    NOT be applied to the subspace.

Example: %BLOCK CONSTRAINED\_DFT

# L Z Uh Uq(UP) Uq(DOWN) Us N(UP) N1(DOWN) [N1(UP)-N1(DOWN)]

N1 1 -5. 0.0 11.0 11.0 0.0 2.3 1.3 0.

N2 1 -5. 0.0 -26.0 -26.0 0.0 2.7 2.7 0.

%ENDBLOCK CONSTRAINED\_DFT

**MAXIT\_CDFT\_U\_CG**

Syntax: MAXIT\_CDFT\_U\_CG [Integer]

Description: Specifies the maximum number of iterations for the
constraining potentials (Uq/s) conjugate gradients optimisation.

Default: 60

Example: MAXIT\_CDFT\_U\_CG 5

**HUBBARD\_TENSOR\_CORR**

Syntax: HUBBARD\_TENSOR\_CORR [Integer]

Description:

    1: Correct tensorially for the slight nonorthogonality between DFT+U or constrained DFT (cDFT) projectors of numerical pseudoatomic orbital form, individually on each atom, which arise due to finite psinc sampling. See details see Phys. Rev. B 83, 245124 (2011).

    2: Use the full simulation-cell overlap matrix of the initial numerical pseudoatomic orbitals (whether or not they are selected as projectors) to form the nonorthogonality correction, in the vein of Mulliken analysis. Experimental, non-Hermitian at present, and not recommended.

    3: Do not correct for the slight nonorthogonality between DFT+U or constrained DFT (cDFT) projectors on a given atom. This is standard in many codes, and currently necessary to choose when using USP/PAW.

    4: This is a reasonable, but not necessary choice, when using cDFT with constraints based on atom group populations. The non-negligible nonorthogonality between projectors on different atoms in the group is accounted for tensorially. This also activates multi-site Pulay force terms in constrained DFT (cDFT) that account for varying inter-atom nonorthogonality. For details see Phys. Rev. B 97, 205120 (2018).

    5: This is also an arguably reasonable, but not necessary choice, when using cDFT with constraints based on, e.g. the difference of atom or atom group populations, that is source-drain cDFT. The non-negligible nonorthogonality between projectors on different atoms in a group is accounted for tensorially, and also the possible nonorthogonality between the source and drain atoms or atom-groups. For an application see Phys. Rev. B 93, 165102 (2016). This also activates (in principle) multi-site Pulay force terms in constrained DFT (cDFT) that account for varying inter-subspace nonorthogonality. This also activates multi-site Pulay force terms in constrained DFT (cDFT) that account for varying inter-atom nonorthogonality. This also activates multi-site Pulay force terms in constrained DFT (cDFT) that account for varying inter-atom nonorthogonality. For details see Phys. Rev. B 97, 205120 (2018).

   The HUBBARD_TENSOR_CORR functionality is not activated in the rare case that analytical hydrogenic projectors are instead of the default numerical pseudoatomic ones.

Default: 1

Example: HUBBARD\_TENSOR\_CORR 4
