===============================================================================
Density derived electrostatic and chemical (DDEC) electron density partitioning
===============================================================================

:Author: Louis P. Lee, University of Cambridge
:Author: Daniel J. Cole, University of Cambridge

Atoms-in-molecule electron density partitioning is a useful
post-processing analysis tool for computing atomic charges (as well as
higher order atomic multipoles) from the total electron density. ONETEP
uses the DDEC3 method [1,3] for this purpose, as the computed charges
are both chemically meaningful and reproduce the electrostatic potential
of the underlying QM calculation. Options are also available for
computing *Hirshfeld* and *iterated stockholder atoms* (ISA) charges
[3,4].

A DDEC3 calculation to partition the electron density and output atomic
charges, multipoles and volumes is performed by specifying:

::

   ddec_calculate : T
   ddec_multipole : T
   ddec_moment : 3

along with the ddec\_rcomp block for your system (see below). Iterated
stockholder atoms (ISA) partitioning may be performed instead by
additionally specifying:

::

   ddec_IH_fraction : 0.00

Classical Hirshfeld partitioning may be performed instead by
additionally specifying:

::

   ddec_classical_hirshfeld : T
   ddec_IH_fraction : 1.00
   ddec_maxit : 1

| The reference ion densities for use with DDEC3 are read in from an
  external library kindly provided by Thomas A. Manz and Nidia Gabaldon
  Limas (please cite Refs. [1,2]), and are available for download from
  the ONETEP website:
| http://www.onetep.org/pmwiki/uploads/Main/Utilities/ddec\_atomic\_densities.tar.gz
| The paths to the reference densities are specified in the block
  ``ddec_rcomp``. Specify one core and one total density file for each
  species in your system (except for hydrogen and helium which do not
  require a core density file). The example below is for methanol:

::
  
   %block ddec_rcomp
     H ALL “H_c2.refconf”
     O ALL “O_c2.refconf”
     O CORE “O_c2.coreconf”
     C ALL “C_c2.refconf”
     C CORE “C_c2.coreconf”
   %endblock ddec_rcomp

References
==========

| For the development of the DDEC method:
| :math:`[1]` T.A. Manz and D.S. Sholl, “Improved Atoms-in-Molecule
  Charge Partitioning Functional for Simultaneously Reproducing the
  Electrostatic Potential and Chemical States in Periodic and
  Non-Periodic Materials,” J. Chem. Theory Comput. 8 (2012) 2844-2867.
| :math:`[2]` T. A. Manz and D. S. Sholl, “Chemically Meaningful Atomic
  Charges that Reproduce the Electrostatic Potential in Periodic and
  Nonperiodic Materials”, J. Chem. Theory Comput. 6 (2010) 2455-2468.
| And its implementation in ONETEP:
| :math:`[3]` L. P. Lee, N. Gabaldon Limas, D. J. Cole, M. C. Payne,
  C.-K. Skylaris, T. A. Manz, “Expanding the Scope of Density Derived
  Electrostatic and Chemical Charge Partitioning to Thousands of Atoms”,
  J. Chem. Theory Comput., 10 (2014) 5377.
| :math:`[4]` L. P. Lee, D. J. Cole, C.-K. Skylaris, W. L. Jorgensen, M.
  C. Payne, “Polarized Protein-Specific Charges from Atoms-in-Molecule
  Electron Density Partitioning”, J. Chem. Theory Comput., 9 (2013),
  2981.
