=================================
Electron Localisation Descriptors
=================================

:Author: Rebecca J. Clements, University of Southampton
:Author: James C. Womack, University of Southampton
:Author: Chris-Kriton Skylaris, University of Southampton
	 
:Date:   July 2019

ONETEP now has the electron localisation descriptors:

-  Electron Localisation Function (ELF) [Becke1990]_

-  Localised Orbital Locator (LOL) [Becke2000]_

These descriptors provide a visualisation tool for indicating the
location of electron pairs, including bonding and lone pairs, and
distinguishing between :math:`\sigma` and :math:`\pi` bonds. It is a
dimensionless quantity in the range of :math:`0` and :math:`1`. The
plots can be output in :math:`.cube` or :math:`.dx` format for
visualisation as isosurfaces or volume slices.

The ELF was first introduced by Edgecombe and Becke in 1990 for
Hartree-Fock Theory, and later updated for Density Functional Theory by
Savin *et al.*.
[Becke1990]_, [Becke2000]_, [Savin1992]_ It is based
on the Hartree Fock probability of finding two particles of the same
spin :math:`\sigma` at two different positions of a multielectron
system. From this, Edgecombe and Becke obtained the conditional
probability of an electron existing in the proximity of a reference
electron, which relates to electron localisation. The smaller the
probability, the increased likelihood that the reference electron is
localised, provided both the electrons are of the same spin. This
probability vanishes to zero when the two electrons have the same
position, in agreement with the Pauli principle.

This probability is not upper-bounded and so for more convenient
graphical interpretation, it is inverted. Localisation is represented at
unity. Details are below.

The LOL is similar to the ELF, but with a simpler representation, and produces cleaner plots in some cases. [Schmider2004]_

The ELF provides quantum Valence Shell Electron Pair Repulsion (VSEPR)
representation of coordination compounds, and the identification of
covalent bonding across crystalline solids and surfaces. This provides a
useful tool for the main applications of ONETEP; biomolecular
simulations, catalysis, and the design of nanostructured materials.

Electron Localisation Function
==========================================

The starting point to the ELF formula is the exact kinetic energy
density for spin :math:`\sigma`,

.. math::
   :label: kedensity

   \tau_{\sigma}^{exact}(\textbf{r}) = \sum_{\alpha} \tau(\textbf{r};\alpha),

where
:math:`\tau(\textbf{r};\alpha)=(\nabla \psi_{\alpha}(\textbf{r})) \cdot \bigg( \nabla \sum_{\beta}K^{\alpha\beta}\psi_{\beta}(\textbf{r}) \bigg)`

defined in ONETEP [Womack2016]_ in terms of NGWFs
where :math:`\beta` are all which overlap with :math:`\psi_{\alpha}`.
Note that the standard :math:`\frac{1}{2}` coefficient in ONETEP is
omitted for the purposes of the ELF, to follow the definition by Becke.
The term :math:`\tau_{\sigma}^{exact}` is extended into a formula to
describe electron localisation, a *non-negative* probability density.
This will be called the Pauli kinetic energy, :math:`D_{\sigma}`. For
the individual spin, the Pauli kinetic energy takes the form:

.. math:: D_{\sigma} = \tau_{\sigma}^{exact} - \frac{1}{4} \frac{\left( \nabla n_{\sigma} \right) ^{2}} {n_{\sigma}}.

where :math:`n_{\sigma}(\textbf{r})` is the charge density for each
spin :math:`\sigma`, and :math:`\nabla n_{\sigma}(\textbf{r})` is its
gradient. :math:`D_{\sigma}` is compared with the uniform electron gas
as a reference, by taking the ratio:

.. math:: \chi_{\sigma} = \frac {D_{\sigma}} {D_{\sigma}^{0}}.

where

.. math:: D_{\sigma}^{0} = \frac{3}{5} \left( 6\pi^{2} \right) ^{\frac{2}{3}} n_{\sigma}^{\frac{5}{3}}.

The charge density is the local value of
:math:`n_{\sigma} \left( \textbf{r} \right)` here and the coefficient is
the Fermi constant.

Hence, the measure of electron localisation has now become
dimensionless. :math:`\chi_{\sigma}` is then reformulated to avoid the
open bounds of the above formula, limiting the ELF to a more desirable
finite range of values of :math:`0` to :math:`1` for visual
representation:

.. math:: ELF = \frac{1}{1+\chi_{\sigma}^{2}}

Localised Orbital Locator
=====================================

The LOL is similar to the ELF. Again, the starting point is the exact
kinetic energy for spin :math:`\sigma`, as in Equation :eq:`kedensity`,
with the :math:`\frac{1}{2}` coefficient omitted.

The uniform electron gas reference, :math:`D_{\sigma}^{0}`, is used
again here, also known as the local spin density approximation (LSDA):

.. math:: \tau_{\sigma}^{LSDA} = \frac{3}{5} \left( 6\pi^{2} \right) ^{\frac{2}{3}} n_{\sigma}^{\frac{5}{3}}.

The ratio follows a similar structure to the ELF, except for it is
inverted, which again, makes the quantity dimensionless:

.. math:: t_{\sigma} = \frac {\tau_{\sigma}^{LSDA}} {\tau_{\sigma}^{exact}}.

The quantity is reformulated to change the infinite range into a
:math:`0` to :math:`1` range like before:

.. math:: v_{\sigma} = \frac{t_{\sigma}} {1 + t_{\sigma}}

How to use
==========

Keywords
--------

The keywords related to the implementation of the electron localisation
descriptors are as follows:

::

    Keyword:                Options (default):
    eld_calculate           T/F (F)
    eld_function            ELF/LOL (ELF)
    ke_density_calculate    T/F (F)
    do_properties           T/F (F)
    cube_format             T/F (T)
    dx_format               T/F (F)

-  Setting :math:`eld\_calculate` to true turns on the calculation. The
   calculation will not proceed if this keyword is missing or if it set
   to false.

-  The keyword :math:`eld\_function` determines which of the ELF or LOL
   ONETEP is to calculate, by specifying either string. The default here
   is the ELF, provided the keyword :math:`eld\_calculate` has been
   specified.

-  As part of this implementation, the kinetic energy density can now
   also be output, using the logical keyword
   :math:`ke\_density\_calculate`. This does not automatically output
   with :math:`eld\_calculate`.

-  Electron localisation descriptors and kinetic energy density are
   available in the formats of :math:`.cube` or :math:`.dx` files.

-  For spin polarised systems, there will be an ELF output for each spin
   individually, showing the electron localisation for one of the spins.

In order to use any of the above keywords, ONETEP’s properties
calculation must be enabled, using :math:`do\_properties` or setting the
task to :math:`properties`, if reading in density results of an energy
minimisation calculation. To produce the density plot during the
original energy calculation, the input should include:

::

    task                singlepoint
    write_density_plot  T

Example input file
------------------

Below is an example input for using the ELF, for the water molecule:

::

    task            singlepoint
    cutoff_energy   900.0 eV
    maxit_ngwf_cg   100
    output_detail   verbose
    do_properties   T
    cube_format     T
    dx_format       F
    grd_format      F
    eld_calculate   T
    eld_function    ELF

    %block lattice_cart
    40.000000000000    0.000000000000    0.000000000000
     0.000000000000   40.000000000000    0.000000000000
     0.000000000000    0.000000000000   40.000000000000
    %endblock lattice_cart

    %block positions_abs
    O 20.000000000  20.000000000  20.000000000
    H 18.565580829  18.889354011  20.000000000
    H 21.434419171  18.889354011  20.000000000
    %endblock positions_abs

    %block species
    O O 8 4 8.0
    H H 1 1 8.0
    %endblock species

    %block species_pot
    O <path to oxygen.recpot>
    H <path to hydrogen.recpot>
    %endblock species_pot

[Becke1990] A. D. Becke and K. E. Edgecombe. A simple measure of electron localization in atomic and molecular systems. *J. Chem. Phys.*, 92(9):5397-5403, 1990. 

[Becke2000] A. D. Becke and H. L. Schmider. Chemical content of the kinetic energy density. *Journal of Molecular Structure (Theochem)*, 527:51–61, 2000.

[Savin1992] A. Savin, O. Jepsen, J. Flad, O. K. Andersen, H. Preuss, and H. G. von Schnering. Electron localization in solid-state structures of the elements: the diamond structure. *Angewandte Chemie International Edition in English*, 31(2):187–188, 1992.

[Schmider2000] H. Schmider and A. Becke. Chemical content of the kinetic energy density. *Journal of Molecular Structure (Theochem)*, 527(1):51 – 61, 2000.

[Womack2016] J. C. Womack, N. Mardirossian, M. Head-Gordon, and C.-K. Skylaris. Self-consistent implementation of meta-gga functionals for the ONETEP linear-scaling electronic structure package. *The Journal of Chemical Physics*, 145(20):204114, 2016.
