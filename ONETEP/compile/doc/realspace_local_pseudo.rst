=========================================
Realspace local pseudopotential in ONETEP
=========================================

:Author: Jacek Dziedzic, University of Southampton
:Author: Chris-Kriton Skylaris, University of Southampton

Motivation
==========

In standard ONETEP the local pseudopotential is obtained in reciprocal
space by a discrete Fourier transform, by assuming the cell is
periodically repeated in space. However, there are certain use-cases,
where one is interested in the properties of an isolated (not
periodically repeated) system. This is especially true if other energy
terms, such the Hartree energy or the ion-ion energy are already
calculated with open boundary conditions, which is the case, e.g., for
implicit solvent calculations in ONETEP.

Theory
======

Assume that :math:`{v_{loc}\left(\vec{r}\right)}` is located on an atom
:math:`A` at a position :math:`\vec{R}_A` and we want to determine the
contribution to the local pseudopotential coming from this atom. Owing
to the spherical symmetry of the potential, we have

.. math:: {v_{loc,A}\left(\vec{r}\right)}= v_{loc}\left(\vec{r}-\vec{R}_A\right) = v_{loc}\left(\vert\vec{r}-\vec{R}_A\vert\right).

The local pseudopotential is given to us in terms of its continuous
Fourier coefficients,
:math:`{\tilde{v}_{loc}\left(\vert\vec{g}\vert\right)}`, read from a
recpot file. To generate the pseudopotential at a point :math:`\vec{r}`
in real space, we use the continuous Fourier transform:

.. math:: v_{loc}\left(\vec{r}-\vec{R}_A\right) = \frac{1}{{\left(2\pi\right)}^{3}}\int {\tilde{v}_{loc}\left(\vec{g}\right)}e^{i\vec{g}\cdot \left(\vec{r}-\vec{R}_A\right)}d\vec{g}=\int{\tilde{v}_{loc}\left(\vec{g}\right)}e^{i\vec{g}\cdot \vec{x}}d\vec{g},

where we have set :math:`\vec{x}=\vec{r}-\vec{R}_A`. Expanding the
plane wave :math:`e^{i\vec{g}\cdot\vec{x}}` in terms of localised
functions, we get

.. math::

   {v_{loc,A}\left(\vec{r}\right)}= 
   \frac{1}{{\left(2\pi\right)}^{3}}
   \int {\tilde{v}_{loc}\left(\vec{g}\right)}\cdot \left[ 4\pi 
   \sum_{l=0}^{\infty}
   \sum_{m=-l}^{l} i^l 
   j_l\left( gx\right)
   Z_{lm}\left(\Omega_{\vec{g}}\right)
   Z_{lm}\left(\Omega_{\vec{x}}\right)
   d\vec{g}
   \right],

.. math::
   :label: eq4

   {v_{loc,A}\left(\vec{r}\right)}= 
   \frac{1}{{\left(2\pi\right)}^{3}}
   4\pi
   \sum_{l=0}^{\infty}
   \sum_{m=-l}^{l} i^l 
   Z_{lm}\left(\Omega_{\vec{x}}\right)
   \underbrace{
   \int {\tilde{v}_{loc}\left(\vec{g}\right)}j_l\left( gx\right)
   Z_{lm}\left(\Omega_{\vec{g}}\right)
   d\vec{g}}_{I_1}.

The orthogonality of harmonics means that all of the terms, except for
that of :math:`l=m=0`, disappear and, after a change of coordinates
(:math:`g^2\sin{\theta}` being the Jacobian), we obtain a new expression
for the integral in (:eq:`eq4`):

.. math::

   I_1=
   \int\limits_{0}^{2\pi}
   \int\limits_{0}^{\pi}
   Z_{lm}\left(\Omega_{\vec{g}}\right)
   Z_{00}
   \sin{\theta}\,d\theta\,d\varphi
   \cdot\int\limits_{0}^{\infty}{\tilde{v}_{loc}\left(g\right)}{}j_l\left( gx\right)g^2 dg.

With :math:`Z_{00}=1/{\sqrt{4\pi}}`, the double integral simplifies to
1 and we obtain, after realizing that all terms except for :math:`l=0`
disappear,

.. math::

   {v_{loc,A}\left(\vec{r}\right)}= 
   \frac{1}{{\left(2\pi\right)}^{3}}
   4\pi
   \int {\tilde{v}_{loc}\left(g\right)}j_0\!\left( gx\right)g^2 dg
   =
   \frac{1}{{\left(2\pi\right)}^{3}}
   4\pi
   \int {\tilde{v}_{loc}\left(g\right)}\frac{\sin\left( g    x\right)}{gx}g^2\,dg.

ONETEP uses a convention where an additional factor of :math:`4\pi` is
needed when transforming between real and reciprocal space. Thus the
final formula for the local pseudopotential at a distance of :math:`x`
from an atom of species :math:`s` becomes

.. math::
   :label: eq7

   {v^s_{loc}\left(x\right)}= \frac{2}{\pi}\int\limits_0^{\infty}
   {\tilde{v}^s_{loc}\left(g\right)}\frac{\sin\left(gx\right)}{x}g\,dg.

Implementation
==============

In practice, however, it is not possible to evaluate the integral
:eq:`eq7` with :math:`\infty` as the upper limit, because
:math:`{\tilde{v}^s_{loc}\left(g\right)}` is defined in the recpot file
only up to a :math:`g_{max}` of 100 Å\ :math:`^{-1}`. Furthermore, to
ensure the results are consistent with standard ONETEP, we must lower
this limit even more, to prevent aliasing, as high :math:`g`\ ’s will
not be representable on our reciprocal space grid. Thus, in practice we
evaluate

.. math::
   :label: eq8

   {v^s_{loc}\left(x\right)}= \frac{2}{\pi}\int\limits_0^{g_{cut}}
   {\tilde{v}^s_{loc}\left(g\right)}\frac{\sin\left(gx\right)}{x}g\,dg,

where :math:`g_{cut}=2\pi\max{\left(d_1,d_2,d_3\right)}` (:math:`d_i`
being the grid spacings of ``pub_cell``) and will usually be in the
order of 20-30 :math:`a_0^{-1}`.

The integral is evaluated for :math:`x`\ ’s on a fine radial grid
running from :math:`0` to the maximum possible distance, which is the
magnitude of the cell diagonal. The calculation is distributed across
nodes (each node deals with a portion of the fine radial grid). The
total pseudopotential for any point on the real space fine grid is
evaluated by interpolation from the fine radial grid and by summing over
all atoms. This calculation is distributed across nodes as well (each
node deals with its own slabs of the real space fine grid). The default
number of points in the radial grid is 100000 and can be changed with
the directive ``openbc_pspot_finetune_nptsx``.

The integral :eq:`eq8` is difficult to evaluate numerically. One source of
difficulties is the oscillatory nature of :math:`\sin\left(gx\right)`.
For larger cells, where the maximum interesting :math:`x` is in the
order of :math:`100\,a_0`, this oscillates so rapidly that the
resolution of the recpot file (0.05 Å\ :math:`^{-1}`) is not enough and
it becomes necessary to interpolate
:math:`{\tilde{v}^s_{loc}\left(g\right)}`, and the whole integrand,
between the :math:`g`-points specified in the recpot file. The result of
the interpolation is stored on a fine radial :math:`g`-grid, which is
:math:`f` times as fine as the original radial :math:`g`-grid of the
recpot file. :math:`f` is determined automatically so that every full
period of :math:`\sin\left(gx\right)` is sampled by at least 50 points.
For typical cells, this yields :math:`f` in the order of 5-50, depending
on the cell size. Alternatively, :math:`f` may be specified manually by
the ``openbc_pspot_finetune_f`` directive.

Another difficulty is caused by the singularity in
:math:`{\tilde{v}^s_{loc}\left(g\right)}` as :math:`g\to0`, where the
behaviour of :math:`{\tilde{v}^s_{loc}\left(g\right)}` approaches that
of :math:`-Z_s/g^2`. Although the integral is convergent, this
singularity cannot be numerically integrated in an accurate fashion. The
singularity also presents problems when interpolating between the
:math:`g` points – the usual cubic interpolation of
``services_1d_interpolation`` becomes inaccurate at low :math:`g`\ ’s.
The second problem is solved by subtracting the Coulombic potential,
:math:`-Z_s/g^2`, before interpolation to the fine radial :math:`g`-grid
and then adding it back. The first problem is difficult to treat. An
approach where at low :math:`g`\ ’s
:math:`{\tilde{v}^s_{loc}\left(g\right)}` is assumed to be exactly equal
to :math:`-Z_s/g^2` (which allows the low-\ :math:`g` part of integral
(:eq:`eq8`) to be evaluated analytically) gives better results than
attempting to numerically integrate the singularity, but is not accurate
enough, leading to errors in the order of :math:`50-100\,\mu{}`\ Ha in
the energy for a hydrogen atom test-case (with a total energy of ca.
0.477Ha. Attempting to fit :math:`A/g^2+B/g+C` (which also allows
analytical integration at low :math:`g`\ ’s) gives similar results. The
numerical inaccuracy presents itself as a near-constant shift of the
obtained pseudopotential and clearly affects total energy.

To solve this problem, we observe that the local pseudopotential can be
split into a long-range part and a short-range part:

.. math:: {v^s_{loc}\left(x\right)}= {v^{s (long)}_{loc}\left(x\right)}+ {v^{s (short)}_{loc}\left(x\right)},

.. math:: {\tilde{v}^s_{loc}\left(g\right)}= {\tilde{v}^{s (long)}_{loc}\left(g\right)}+ {\tilde{v}^{s (short)}_{loc}\left(g\right)}.

Following [Martyna1999], we observe that
:math:`{\tilde{v}^{s (long)}_{loc}\left(g\right)}=\frac{4\pi}{g^2}\exp{\left(\frac{-g^2}{4\alpha^2}\right)}`
(where :math:`\alpha` is an adjustable parameter, controllable with
``openbc_pspot_finetune_alpha``) which easily transforms to real space
to give
:math:`{v^{s (long)}_{loc}\left(x\right)}=-\frac{\operatorname{erf}{\left(\alpha{}x\right)}}{x}`
and is conveniently calculated in real space. The short-range part
(corresponding to high :math:`g`\ ’s) is
:math:`{\tilde{v}^{s (long)}_{loc}\left(g\right)}={\tilde{v}^s_{loc}\left(g\right)}\cdot\left[1-\exp{\left(\frac{-g^2}{4\alpha^2}\right)}\right]`.
In this way, the integral (:eq:`eq8`) can be rewritten as

.. math::
   :label: eqsplit

   {v^s_{loc}\left(x\right)}= 
   -\frac{\operatorname{erf}{\left(\alpha{}x\right)}}{x}+
   \frac{2}{\pi}
   \underbrace{
   \int\limits_0^{g_{cut}}
   {\tilde{v}^s_{loc}\left(g\right)}\cdot \left[ 1 - exp\left(\frac{-g^2}{4\alpha^2}\right) \right]
   \cdot \frac{\sin\left(gx\right)}{x}g\,dg}_{I_s(x)}
   .

Owing to the
:math:`\left[1-\exp{\left(\frac{-g^2}{4\alpha^2}\right)}\right]` factor,
the integral :math:`I_s(x)` is no longer singular at :math:`g=0` and can
be accurately evaluated numerically, if :math:`\alpha` is large enough.
Small values of :math:`\alpha` make the numerical integration more
difficult (requiring larger values for :math:`f`), because the
oscillations at low :math:`g`\ ’s are large in magnitude. Larger values
of :math:`\alpha` allow for easy integration, but they cause the
long-range behaviour to “kick in” earlier. As this long-range behaviour
is calculated in real space, it lacks the oscillations that are present
in standard ONETEP because of a finite value for :math:`g_{cut}`. Even
though these oscillations are an artifact, obtaining a long-range
behaviour that is physically more correct, but without the oscillations,
leads to aliasing in reciprocal space and to a departure from the
results of standard ONETEP. For this reason we want :math:`\alpha` to be
as small as possible, without negatively impacting the numerical
integration. The accuracy of the obtained method can be judged by
comparing the real space tail of the obtained pseudopotential with the
Coulombic potential. Since we expect the obtained pseudopotential to
oscillate slightly around :math:`-Z_s/x`, a good measure of accuracy,
which we will call :math:`b`, is the average value of
:math:`\dfrac{{v^s_{loc}\left(x\right)}-(-Z_s/x)}{-Z_s/x}` over the tail
of the pseudopotential, from, say, 5\ :math:`a_0` to the maximum
:math:`x` for which :math:`{v^s_{loc}\left(x\right)}` is evaluated.
Ideally, :math:`b` should be zero. Numerical inaccuracies will cause a
shift in :math:`{v^s_{loc}\left(x\right)}` which will present itself as
a finite, non-zero value of :math:`b`. Naïve numerical integration by a
direct calculation of (:eq:`eq8`) led, for our test-case, to :math:`b` in
the order of 0.01, which can be reduced by an order of magnitude by
using a very fine radial :math:`g`-grid (high value of :math:`f`).
Subtracting out the Coulombic potential and integrating only the
difference between :math:`{\tilde{v}^s_{loc}\left(g\right)}` and the
Coulombic potential numerically, while integrating the remaining part
analytically reduced b to about 0.0005. Application of the proposed
formula (:eq:`eqsplit`) yielded :math:`b=5\cdot10^{-8}` for
:math:`\alpha=0.5/l` and :math:`b=3\cdot10^{-9}` for
:math:`\alpha=0.1/l` with a suitably large :math:`f` to ease the
numerical integration at low :math:`g` (:math:`l` is the box length).
With the default value for :math:`f`, the total energy is not sensitive
(to more than 0.0001%) to the choice of :math:`\alpha`, provided it is
in a resonable range of :math:`0.1/l - 2/l`. The value of :math:`0.3/l`
was chosen as a default.

The calculation of the realspace local pseudo is implemented in
``norm_conserv_pseudo.F90`` in the subroutine
``pseudo_local_on_grid_openbc`` and its internal subroutine
``internal_Is_of_x``, which evaluates :math:`I_s(x)`. A typical
calculation would use default values for all the parameters. The
realspace local pseudo is off by default and is turned on automatically
when smeared ions or implicit solvent is in use. It can also be forced
to be on (for development tests) by using ``openbc_pspot T``.

.. list-table:: Directives controlling the calculation of the realspace local pseudo
   :widths: 40 20 40
   :header-rows: 1

   * - Directive
     - Action
     - Rationale for use
   * - ``openbc_pspot T``
     - Forces the realspace pseudo to be used
     - Normally not needed, the realspace pseudo will be turned on when necessary. This directive allows turning it on even though the Hartree potential calculation and Ewald calculation proceed in reciprocal space, which might be useful for certain test calculations. A related directive, ``openbc_ion_ion T`` may be used in conjuction, to replace Ewald with a direct Coulombic sum.
   * - ``openbc_pspot_finetune_f`` :math:`value` [:math:`value` is an integer.]
     - Sets the fineness parameter, :math:`f`, to :math:`value`.
     - Default value of -1 causes :math:`f` to be determined automatically. Positive values can be used to increase :math:`f` to obtain extra accuracy. Decreasing :math:`f` will reduce accuracy and is not recommended.
   * - ``openbc_pspot_finetune_nptsx`` :math:`value` [:math:`value` is an integer.]
     - Sets the number of radial grid points (distinct values of :math:`x`) to :math:`value`.
     - The default of 100000 should be enough, unless huge boxes are used, where it might make sense to increase it. Decreasing this value is not recommended, as it will impact accuracy.
   * - ``openbc_pspot_finetune_alpha`` :math:`value` [:math:`value` is a real.]
     - Sets the short-range-long-range crossover parameter :math:`\alpha` to :math:`value/l`,  where :math:`l` is the maximum dimension of the cell.
     - A default value of 0.3 should be OK for most applications. Increasing :math:`\alpha` will reduce the numerical inaccuracy in :math:`I_s(x)`, but will cause the long-range behaviour to lack the oscillations of usual ONETEP and thus increase aliasing. Decreasing :math:`\alpha` will make :math:`I_s(x)` inaccurate, which can be helped, to a certain extent, by increasing :math:`f`.


[Martyna1999] G. J. Martyna and M. E. Tuckerman *J. Chem. Phys.* **110** (1999).
