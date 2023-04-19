=============================================
Empirical Dispersion Correction
=============================================

:Author: Max Phipps, University of Southampton
:Author: Arihant Bhandari, University of Southampton
	 
:Date:   November 2020

Theory
======

Energy
------

The modification to the total DFT energy when correcting for
dispersion is given by,

  .. math::
     :label: disp_corr
	     
       E_{DFT-D} = E_{KS-DFT} + E_{disp}

where the dispersion energy correction is given by
raw-latex:[Vasp]_, [Hill2009]_:

  .. math::
     :label: E_disp
     
       E_{disp} = -\frac{1}{2}\cdot s_6\cdot \sum^{N_{at}}_{i=1} \sum^{N_{at}}_{j=1}
     \sum_{L}^* \dfrac{C_{6,ij}}{R^6_{ij}} f_{damp}( R_{ij} ) g_{smooth}(R_{ij})

where:

-  :math:`s_6` is a DF-dependent global scaling factor,

-  :math:`C_{6,ij}` is a dispersion coefficient for the atom pair
   :math:`ij`,

-  :math:`L` denotes all translations of the unit cell within the van
   der Waals radial cutoff :math:`R_{cut}`,

-  :math:`*` denotes :math:`i\ne j` in :math:`L=0`

-  :math:`R_{ij}` is the distance :math:`| \overrightarrow{R^{0}_{i}} -
   \overrightarrow{R^L_{j}} |` between atom i in the parent cell
   :math:`L=0` and the atom j in all possible translations :math:`L`,

-  :math:`f_{damp}(R_{ij})` is a damping function that is unity at large
   distances and zero at small
   distances [Hill2009]_, [Grimme2006]_,

-  :math:`g_{smooth}(R_{ij})` is smoothening function for truncating the
   interactions beyond van der Waals radial cutoff.

-  The dispersion coefficients and the form of the damping function are
   dependent upon the empirical vdW correction model adopted. The
   damping functions for the different models are given below:

.. list-table::
   :widths: 8 40 32 20
   :header-rows: 1

   * - Index
     - Description
     - Damping function :math:`f_{damp}( R_{ij} )`
     - vdW radii (:math:`R_{0,ij}`)
   * - 1
     - Damping function of Elstner [Elstner2001]_.
     - :math:`(1-\exp(-c_{damp}(R_{ij}/R_{0,ij})^7))^4`
     - :math:`\frac{{R_{0,i}}^3+{R_{0,j}}^3}{{R_{0,i}}^2+{R_{0,j}}^2}`
   * - 2
     - First damping function of Wu and Yang [Wu2002]_.
     - :math:`(1-\exp(-c_{damp}(R_{ij}/R_{0,ij})^3))^2`
     - :math:`\frac{{R_{0,i}}^3+{R_{0,j}}^3}{{R_{0,i}}^2+{R_{0,j}}^2}`
   * - 3
     - Second damping function of Wu and Yang [Wu2002]_.
     - :math:`\frac{1}{1+\exp(-c_{damp}(R_{ij}/R_{0,ij}-1))}`
     - :math:`\frac{{R_{0,i}}^3+{R_{0,j}}^3}{{R_{0,i}}^2+{R_{0,j}}^2}`
   * - 4
     - Damping function of D2 correction of Grimme [Grimme2006]_.
     - :math:`\frac{1}{1+\exp(-c_{damp}(R_{ij}/R_{0,ij}-1))}`
     - :math:`{R_{0,i}}+{R_{0,j}}`

where :math:`c_{damp}` is a damping constant (referred to as :math:`d`
within the literature of Grimme [Grimme2006]_ and
:math:`c_{damp}` by Hill [Hill2009]_) and
:math:`R_{0,ij}` is determined by the vdW radii of the atomic pair
:math:`i` and :math:`j`.

-  The smoothening function has the following form for all dispersion
   models:

   .. math:: g_{smooth}(R_{ij})=1-e^{-\left(R_{ij}-R_{cut}\right)^2}

-  | where, :math:`R_{cut}` is the radial cutoff for van der Waals
     interactions, and can be set by the following keyword in the input
     file:
   | ``vdw_radial_cutoff: x units`` [Real physical, default
     ``vdw_radial_cutoff: 100 bohr``\ ].

Forces
------

During geometry optimization, the van der Waals forces are calculated
from the derivative of the dispersion energy with respect to the ionic
coordinates :math:`\{s_i\}`:

.. math::

   \frac{\partial E_{disp}}{\partial s_i}  = -\frac{1}{2}\cdot s_6\cdot
   \sum^{N_{at}}_{j=1} \sum_{L}^* \dfrac{C_{6,ij}}{R^6_{ij}} \left[f( R_{ij} )
   g'(R_{ij})+g(R_{ij})\left(f'(R_{ij})-\frac{6f(R_{ij})}{R_{ij}}\right)\right]\frac{\partial
   R_{ij}}{\partial s_i}

Activating the dispersion corrections
=====================================

Four vdW correction options have been implemented within ONETEP.
Activation of the vdW corrections within ONETEP is achieved using the

``DISPERSION``

keyword followed by the dispersion index. eg. For the D2 correction of
Grimme [Grimme2006]_,

``DISPERSION 4``

The exchange-correlation functionals available with optimized parameters
for the dispersion models are given below:

.. list-table::
   :widths: 10 20 70
   :header-rows: 1

   * - Index
     - Keyword
     - Available XC functionals
   * - 1
     - ELSTNER
     - BLYP, PBE, PW91, REVPBE, RPBE, XLYP.
   * - 2
     - WUYANG1
     - BLYP, PBE, PW91, REVPBE, RPBE, XLYP.
   * - 3
     - WUYANG2
     - BLYP, PBE, PW91, REVPBE, RPBE, XLYP.
   * - 4
     - GRIMMED2
     - BLYP, PBE, B3LYP.
    
The :math:`C_{6,ij}` coefficients and the coefficients of the damping
function of dispersion corrections 1–3 are optimized for each xc
functional to minimize the root mean square difference of the
interaction energy from the literature values for a selection of
dispersion-dominant biological complexes from the JSCH-2005 and S22
sets [Jurecka2006]_` and the database of
Morgado [Morgado2007]_. For the D2 correction of
Grimme the parameters are optimized as described in the
literature [Grimme2006]_.

The :math:`s_6` parameters are unity for dispersion corrections 1–3 and
are fitted by least squares optimization of interaction energy error for
40 noncovalently bound complexes for the D2 correction of Grimme.

In the case of using an unoptimized functional not given within the
list, default unoptimized :math:`c_{damp}`, :math:`C_{6,ij}`) and
:math:`R_{0,i}` values are adopted for the damping functions of Elstner
or Wu and Yang as described by Hill [Hill2009]_ and a
default :math:`s_6` value of 1.00 is adopted for Grimme’s D2 correction
model.

Overriding dispersion correction parameters
===========================================

It is possible to override the default parameters of the dispersion
damping functions. This option allows the user to specify parameters for
elements and functionals for which values are not given. The
atom-dependent variables :math:`C_{6,i}` (used to calculate
:math:`C_{6,ij}`), :math:`R_{0,i}` (related to the atomic vdW radius of
an atom :math:`i`), and :math:`n_{eff}` (used in the calculation of
:math:`C_{6,ij}` for all damping functions excluding the D2 correction
of Grimme) are modified using the

``vdw_params``

block. This override block applies the parameter changes to atoms by
their atomic number (nzatom). eg. To override the dispersion parameters
associated with nitrogen,

::

    %block vdw_params
    ! nzatom, c6coeff, radzero, neff
      7       21.1200  2.6200   2.5100
    %endblock vdw_params

To override the damping constant :math:`c_{damp}` associated with a
damping function, the keyword

``vdw_dcoeff``

followed by the modified damping constant parameter is used. eg.

``vdw_dcoeff 11.0``

Boundary Conditions
===================

The boundary conditions are set by the following keyword:

``vdw_bc P/O P/O P/O``

which accepts a string which should contain three characters (which may
be separated by spaces), specifying the BCs along the :math:`x`,
:math:`y` and :math:`z` directions of the simulation cell. ‘P’ for
periodic and ‘O’ for open. In case, the keyword is not specified, BCs
set same as

``ion_ion_bc``.


[Vasp] https://www.vasp.at/wiki/index.php/DFT-D2

[Hill2009] Q. Hill and C-K Skylaris, *Proc. R. Soc. A* 465(2103):669-683, 2009

[Grimme2006] S. Grimme, *J. Comput. Chem.* 27(15):1787-1799, 2006

[Elstner2001] M. Elstner, P. Hobza, T. Frauenheim, S. Suhai and E. Kaxiras, *J. Chem. Phys.* 114(12):5149-5155, 2001

[Wu2002] Q. Wu and W. Yang, *J. Chem. Phys.*, 116(2):515-524, 2002

[Jurecka2006] P. Jurečka, J. Šponer, J. Černýa and P. Hobza, *Phys. Chem. Chem. Phys.* 8:1985-1993, 2006

[Morgado2007] C. A. Morgado, J. P. McNamara, I. H. Hillier, N. A. Burton and M. A. Vincent *J. Chem. Theory Comput.* 3(5):1656-1664, 2007
