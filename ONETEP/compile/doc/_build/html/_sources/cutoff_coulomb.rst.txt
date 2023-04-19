=========================
Cut-off Coulomb
=========================

:Author: Nick Hine, University of Warwick (implementation)
:Author: Gabriel Bramley, University of Southampton (documentation)
	 
:Date:   July 2019

Theory
======

The plane wave approach inherently involves the use of periodic boundary
conditions (PBC). In order to model isolated molecules, the supercell
technique was developed, which involves adding large quantities of
vacuum to the simualtion cell in order to separate their periodic
images. Although this method is adequate for neutral molecules, charged
systems and systems with significant multipoles require additional
considerations. The potential of the monopole and multipole moments
decay in accordance to the power law, where point charges decay with
:math:`\frac{1}{r^{1}}`, dipoles with :math:`\frac{1}{r^{2}}` etc..
Using the supercell method for systems with a net charge or significant
dipoles require large volumes of vacuum to eliminate the electrostatic
interactions between the system of the unit cell and its periodic
images. However, as traditional plane wave codes extend across the
entire cell, adding large quantities of vacuum is computationally
costly. Various dipole corrections such as the cut-off Coulomb (CC)
[Jarvis1997]_, Continuum Screening Method
[Otani2006]_ and Gaussian
Counter Charge model
[Dabo2008]_ have been
developed in order to isolate the electrostatic interactions of the unit
cell from its periodic images.

The cut-off Coulomb (or Coulomb cut-off) approach, as implemented in
ONETEP by Nick Hine [Hine2011]_, achieves this by only
allowing Coulombic interactions within a specified region, thereby
emulating the electrostatics of an isolated system with a periodic
representation of the charge density [Jarvis1997]_. By
selecting an appropriate region, one can eliminate spurious
electrostatic interactions between charges in the periodic replicas and
the home cell, while also retaining the correct description of the
potential between charges in the isolated system. One can choose to
maintain periodic Coulombic interactions along specified axes, thereby
representing electrostatics in systems with either 1D (wire), 2D (slab)
or 0D (sphere) periodicity.

3D Periodic Coulomb Interaction
-------------------------------

In a standard periodic calculation, the electrostatic potential is
defined as through the Hartree potential:

.. math::
   :label: hartree

   V_{H}(\mathbf{r}) = \int_{space} {\frac{n(\mathbf{r'})}{|\mathbf{r} -  \mathbf{r'}|}}  d\mathbf{r'}

Where :math:`n(\mathbf{r'})` describes the charge density of the system.
This can be equivalently represented as:

.. math:: V_{H}(\mathbf{r}) = \int_{space} v(\mathbf{r},\mathbf{r'})n(\mathbf{r'}) d\mathbf{r'}

Where :math:`v(\mathbf{r})` represents the Coulomb interaction,
:math:`\frac{1}{|\mathbf{r} -  \mathbf{r'}|}`. This is more conveniently
calculated in reciprocal space, which is achieved through a Fourier
transformation of the :math:`V_{H}` in accordance with convolution
theory:

.. math:: V_{H}(\mathbf{G}) = v(\mathbf{G})n(\mathbf{G})

In the fully periodic case, the above expressions are taken over all
reciprocal space (:math:`\infty` to :math:`-\infty`), which yields a
Coulomb interaction of:

.. math:: v(\mathbf{G}) = \frac{4 \pi}{|\mathbf{G}|^2}

Coulomb Cut-off in 3D
---------------------

The standard periodic approach takes the integration of the Coulomb
interaction over all space, which allows interaction between the
periodic images and the original unit cell. In contrast, Coulomb cut-off
sets the Coulomb interaction between charges to zero beyond a specified
radius. In the simplest case, one can define this cut-off as a sphere,
which assumes a system without periodicity. By selecting an appropriate
cut-off radius, :math:`R_C`, one can retain the correct electrostatic
interactions of between charges in the unit cell, while eliminating
spurious interactions the periodic images:

.. math::

   v^{3D}(\mathbf{r},\mathbf{r'}) =
        \begin{cases}
         \frac{1}{|\mathbf{r} -  \mathbf{r'}|} & \text{for $R_C < |\mathbf{r} -  \mathbf{r'}|$}\\
         0 & \text{for $R_C > |\mathbf{r} -  \mathbf{r'}|$}\\
        \end{cases}

Performing an analytic Fourier Transformation of
:math:`v(\mathbf{r},\mathbf{r'})` with modified boundaries yields the
following reciprocal space representation of the Coulomb interaction:

.. math:: v^{3D}(\mathbf{G}) =  \frac{4 \pi}{\mathbf{G}^2}(1 - \cos(\mathbf{G}R_C))

Coulomb Cut-off in 1D
---------------------

Although this approach is satisfactory for systems with no periodicity,
additional considerations must be made if periodicity is required in
either 1D (ie. a wire) or in 2D (ie. an infinitely extended plane). In
the 1D case, the Coulomb cut-off is defined as a cylinder, where
periodicity is retained in the z-axis and the Coulomb interaction is
applied in the xy-direction. In theory, this redefines the Coulomb
interaction as:

.. math::
   :label: 1DCC

   \begin{aligned}
       v^{1D}(\mathbf{G_{\bot}, G}_x) =  \frac{4 \pi}{\mathbf{G}^2}  [1 + \mathbf{G}_{\bot}R_C J_1 (\mathbf{G}_{\bot}R_C) K_{0}(\mathbf{G}_x R_C)
       - \mathbf{G}_x R_C J_{0}(\mathbf{G}_{\bot}R_C) K_{1}(\mathbf{G}_x R_C) ]
       \end{aligned}

Where :math:`J` and :math:`K` are modified Bessel functions,
:math:`\mathbf{G}_x`, :math:`\mathbf{G}_y` and :math:`\mathbf{G}_z` the
reciprocal lattice vectors in :math:`x`, :math:`y` and :math:`z` and
:math:`\mathbf{G}_{\bot} = \sqrt{\mathbf{G}_x^2 + \mathbf{G}_z^2}`.

However, a divergence occurs at :math:`G_x = 0`, which is handled by
re-casting Equation :eq:`1DCC` into an analytical expression solved through
:

.. math:: v^{1D}(\mathbf{G_{\bot}, G}_x = 0) =  - 4 \pi  \int_{0}^{R} r J_{0}(\mathbf{G}_{\bot})\ln{(\mathbf{r})} d \mathbf{r}

Coulomb Cut-off in 2D
---------------------

For systems where periodicity is maintained in 2D, the Coulomb cut-off
must only be applied in the out-of-plane direction, while retaining PBC
in the xy-plane. Originally, this was implemented in ONETEP from the
formulation of *Rozzi et al.*
[Rozzi2006]_, where the Coulomb
interaction :math:`v^{3D}(\mathbf{G})` is re-cast to the following
expression:

.. math:: v^{2D}(\mathbf{G_{\|}},\mathbf{G}_{z}) = \frac{4 \pi}{\mathbf{G}^2} \bigg \lbrack 1 + e^{-\mathbf{G}_{\|}R_C}\frac{\mathbf{G}_z}{\mathbf{G}_{\|}}\sin(\mathbf{G}_z R_C) - e^{-\mathbf{G}_{\|}R_C}\cos{|\mathbf{G}_z|R_C}) \bigg \rbrack

Where :math:`\mathbf{G}_{\|} = \sqrt{\mathbf{G}_{x}^2+\mathbf{G}_{y}^2}`
and :math:`\mathbf{G}_{z}` represent the in-plane and out-of-plane
reciprocal space vectors respectively. However, as described by *Sohier
et al.* [Sohier2017]_, if
:math:`R_C = \frac{L}{2}`, where L represents the length of the
simulation cell, this expression simplifies to:

.. math:: v^{2D}(\mathbf{\mathbf{G}_{\|}},\mathbf{G}_{z}) = \frac{4 \pi}{\mathbf{G}^2} \bigg \lbrack 1 - e^{-\mathbf{G}_{\|}R_C}\cos({|\mathbf{G}_z|R_C}) \bigg \rbrack

Where :math:`G_z` is a multiple of :math:`\frac{2 \pi}{L}`. As with the
Coulomb interaction under periodic boundary conditions, this term
diverges at :math:`\mathbf{G} = 0`, and is therefore treated separately
and :math:`v^{2D}(\mathbf{\mathbf{G}}=0) = 0` as argued by *Sohier et
al.* [Sohier2017]_.

Performing a Calculation with Coulomb Cut-off
=============================================

To use Coulomb cut-off, the keyword ``COULOMB_CUTOFF_TYPE`` must be
inserted, with the input specifying the periodicity of the system:

-  1D - ``COULOMB_CUTOFF_TYPE: WIRE``

-  2D - ``COULOMB_CUTOFF_TYPE: SLAB``

-  3D - ``COULOMB_CUTOFF_TYPE: SPHERE``

In addition, the length/radius of the cut-off must be specified with
either ``COULOMB_CUTOFF_RADIUS`` or ``COULOMB_CUTOFF_LENGTH``:

-  1D & 2D - ``COULOMB_CUTOFF_LENGTH``

-  3D - ``COULOMB_CUTOFF_RADIUS``

As part of the Coulomb cut-off in ONETEP, the electron density
:math:`n(\mathbf{r})` in the original cell is placed into a larger,
padded cell in which :math:`n(\mathbf{r}) = 0`. This is determined in a
similar way as the original lattice block through a new block
``%BLOCK PADDED_LATTICE_CART``, which determines the size and dimensions
of the larger cell: ``%BLOCK PADDED_LATTICE_CART`` ``a11 a21 a31``
``a21 a22 a23`` ``a31 a32 a33`` ``%ENDBLOCK PADDED_LATTICE_CART``

This is automatically specified, so adjusting this block is not
recommended. The recommended set-ups for calculations of each
dimensionality are summarized in the table below:

+------------------------+---------------------------------------+------------------------------------------+
| Coulomb Cut-off Type   | Cut-off Length/Radius                 | Cell Dimensions                          |
+========================+=======================================+==========================================+
|                        |                                       | :math:`a_{11}^{pad} = 2a_{11}^{cell}`    |
+------------------------+---------------------------------------+------------------------------------------+
| Sphere\*               | :math:`R_C = \sqrt{3}a_{33}^{cell}`   | :math:`a_{22}^{pad} = 2a_{22}^{cell}`    |
+------------------------+---------------------------------------+------------------------------------------+
|                        |                                       | :math:`a_{33}^{pad} = 2a_{33}^{cell}`    |
+------------------------+---------------------------------------+------------------------------------------+
|                        |                                       | :math:`a_{11}^{pad} = 2a_{11}^{cell}`    |
+------------------------+---------------------------------------+------------------------------------------+
| Wire\*\*               | :math:`R_C = \sqrt{2}a_{33}^{cell}`   | :math:`a_{22}^{pad} = 2a_{22}^{cell}`    |
+------------------------+---------------------------------------+------------------------------------------+
|                        |                                       | :math:`a_{33}^{pad} = a_{33}^{cell}`     |
+------------------------+---------------------------------------+------------------------------------------+
|                        |                                       | :math:`a_{11}^{pad} = a_{11}^{cell}`     |
+------------------------+---------------------------------------+------------------------------------------+
| Slab\*\*\*             | :math:`R_C = a_{33}`                  | :math:`a_{22}^{pad} = a_{22}^{cell}`     |
+------------------------+---------------------------------------+------------------------------------------+
|                        |                                       | :math:`a_{33}^{pad} = 2 a_{33}^{cell}`   |
+------------------------+---------------------------------------+------------------------------------------+

Table: The recommended calculation parameters for each periodicity of
the Coulomb cut-off, where :math:`a_{ii}^{cell}` and
:math:`a_{ii}^{pad}` represent the diagonal components of the original
simulation cell specified in ``% BLOCK_LATTICE_CART`` and the padded
cell respectively. Assumed orthogonal cell in all cases.

| \* Assuming :math:`a_{11}^{cell} = a_{22}^{cell} = a_{33}^{cell}`.
| \*\* Assuming :math:`a_{11}^{cell} = a_{22}^{cell}`. :math:`a_{33}` being the periodic direction.
| \*\*\* :math:`a_{33}` defined as the non-periodic direction.

These choices of both the padded cell dimension and the cut-off
length/radius ensure two conditions are satisfied:

#. Charges within the original unit cell correctly interact with one
   another.

#. The interaction between charges of the periodic image and the
   original simulation cell are set to zero.

The first condition is satisfied by setting the cut-off distance equal
to or greater than the distance between any two non-zero charges within
the original unit cell. For the 3D case, this is typically satisfied by
:math:`R_C > \sqrt{3}L_{cell}`, while in 2D and 1D, this is satisfied by
:math:`R_C = \frac{L_{cell}}{2}`, where :math:`L_{cell}` is the cell
dimension in the non-periodic axis. The second condition requires that
the distance between non-zero charges of the simulation cell and the
periodic image must be greater than or equal to the cut-off length. In
ONETEP, this is achieved by placing the unit cell inside a larger padded
cell, in which the charge density :math:`\rho(\mathbf{r})=0`. The second
condition is satisfied when the total cell length,
:math:`L_{total} = L_{cell} + L_{pad} \geq R_C + L_{cell}`.

[Jarvis1997] M. R. Jarvis, I. D. White, R. W. Godby, and M. C. Payne, *Supercell technique for total-energy calculations of finite charged and polar systems*, Phys. Rev. B, **56**, (1997).

[Otani2006] M. Otani, and O. Sugino, *First-principles calculations of charged surfaces and interfaces: A plane-wave nonrepeated slab approach*, Phys. Rev. B, **73**, (2006).

[Dabo2008] I. Dabo, B. Kozinsky, N. E. Singh-Miller, N. Marzari, *Electrostatics in periodic boundary conditions and real-space corrections*, Phys. Rev. B, **77**, (2008).

[Hine2011] N. D. M. Hine, J. Dziedzic, P. D. Haynes, and C.-K. Skylaris, *Electrostatic interactions in finite systems treated with periodic boundary conditions: Application to linear-scaling density functional theory*, J. Chem. Phys. **135** (2011).

[Rozzi2006] C. A. Rozzi, D. Varsano, A. Marini, E. K. U. Gross, and A. Rubio, *Exact Coulomb cutoff technique for supercell calculations*, Phys. Rev. B **73** (2006).

[Sohier2017] T. Sohier, M. Calandra, and F. Mauri, *Density functional perturbation theory for gated two-dimensional heterostructures: Theoretical developments and application to flexural phonons in graphene*, Phys. Rev. B **96** (2017).
