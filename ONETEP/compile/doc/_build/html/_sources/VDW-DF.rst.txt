=================================================
Using van der Waals Density Functionals
=================================================

:Author: Lampros Andrinopoulos, Imperial College London
:Date:   November 2012

Activating vdW-DF
=================

The van der Waals energy is calculated in ONETEP using the van der Waals
Density Functional method, developed by Dion *et al*
[Dion2004]_.

The only input variable needed to activate the vdW-DF functional is to
set ``xc_functional VDWDF``. If a ``vdW_df_kernel`` file is not present
in the working directory, it will be automatically generated.

Theory
======

The form for the exchange-correlation functional proposed by Dion *et
al* is:

.. math::
   :label: E_xc
	   
   E_{xc} = E_{x}^{\rm{revPBE}} + E_c^{\rm{PW92}} + E_c^{\rm{nl}}

where the non-local exchange-correlation energy is given by:

.. math::
   :label: E_nl_dion
	   
   E_{c}^{\mathrm{nl}} = \frac{1}{2}\int\int d{\mathbf{r}}d{\mathbf{r}}'\rho({\mathbf{r}})\phi({\mathbf{r}},{\mathbf{r}}')\rho({\mathbf{r}}')

where :math:`\rho({\mathbf{r}})` is the electron density at
:math:`{\mathbf{r}}` and :math:`\phi({\mathbf{r}},{\mathbf{r}}')` is the
nonlocal exchange correlation kernel whose form is explained in
[Dion2004]_.

Non-local correlation energy
----------------------------

The direct calculation of the integral in the form of Eq. :eq:`E_nl_dion`
is very computationally expensive, as it involves a six-dimensional
spatial integral.

The algorithm proposed later by Roman-Perez and Soler
[Roman-Perez2009]_ improves the efficiency of the
calculation. They observed that with the form used by Dion *et al* for
:math:`\phi`, the above expression can be re-written as:

.. math:: E_{c}^{\mathrm{nl}} = \frac{1}{2}\int\int d{\mathbf{r}}d{\mathbf{r}}'\rho({\mathbf{r}})\phi(q,q',r)\rho({\mathbf{r}}')

where :math:`r=|{\mathbf{r}}-{\mathbf{r}}'|`, and :math:`q` and
:math:`q'` are the values of a universal function
:math:`q_0[\rho({\mathbf{r}}),|\nabla \rho({\mathbf{r}})|]` at
:math:`{\mathbf{r}}` and :math:`{\mathbf{r}}'`. They thus proposed a way
to expand the kernel :math:`\phi` using interpolating polynomials
:math:`p_\alpha(q)` for chosen values :math:`q_\alpha` of :math:`q`, and
tabulated functions :math:`\phi_{\alpha\beta}(r)` for the kernel
corresponding to each pair of interpolating polynomials. The
interpolating polynomials :math:`p_{\alpha}` are cubic splines that
evaluate to a Kronecker delta on each respective interpolating point. A
mesh of 20 interpolation points is used in Soler’s implementation. The
Soler form of the nonlocal energy can be written as:

.. math::
   :label: kernel
	   
   \phi(q_1,q_2,r) = \sum_{\alpha\beta}\phi_{\alpha\beta}(r) p_{\alpha}(q_1) p_{\beta}(q_2)

The universal function :math:`q_0({\mathbf{r}})` is in practice given
by:

.. math::
   :label: q0
	   
   q_0({\mathbf{r}}) = \Bigg(1 + \frac{\epsilon_c^{\rm{PW92}}}{\epsilon_x^{\rm{LDA}}} + 
   \frac{0.8491}{9}\Big(\frac{\nabla\rho}{2\rho k_F}\Big)^2\Bigg) k_F

with :math:`k_F=(3\pi^2\rho)^{1/3}`. The quantity :math:`q_0` is first
“saturated” to limit its maximum value, according to:

.. math:: q_0^{\text{sat}}(\rho,{|\nabla{\rho}|}) = q_c \Bigg(1-\exp\Big(-\sum_{m=1}^{m_c}\frac{(q/q_c)^m}{m}\Big)\Bigg)

where :math:`q_c` is the maximum value of the mesh of
:math:`q_{\alpha}`.

To evaluate this, we first define a quantity
:math:`\theta_{\alpha}({\mathbf{r}}) = \rho({\mathbf{r}}) p_{\alpha}(q(\rho({\mathbf{r}}),\nabla\rho({\mathbf{r}}))`
in real space. In terms of this, Eq. :eq:`E_nl_dion` can be written as:

.. math::
   :label: E_nl_real
	   
   E_c^{\mathrm{nl}} = \frac{1}{2} \sum_{\alpha\beta} \int \int d{\mathbf{r}}d{\mathbf{r}}'
   \theta_{\alpha}({\mathbf{r}}) \theta_{\beta}({\mathbf{r}}') \phi_{\alpha\beta}(r)

It can be shown that this can be written as a reciprocal space integral:

.. math::
   :label: E_nl
	   
       E_c^{\mathrm{nl}} = \frac{1}{2} \sum_{\alpha\beta}\int d\mathbf{k} 
       \theta^{*}_{\alpha}(\mathbf{k})\theta_{\beta}(\mathbf{k})\phi_{\alpha\beta}(k)

Since the kernel is radially dependent in real space, it is only
dependent on the magnitude of the G-vectors, hence the kernel need only
be evaluated as a 1-dimensional function :math:`\phi_{\alpha\beta}(k)`
for each :math:`\alpha`, :math:`\beta`.

The kernel :math:`\phi` and its second derivatives are tabulated for a
specific set of radial points and transformed to reciprocal space. These
values are then used to interpolate the kernel at every point
:math:`\mathbf{k}` in reciprocal space required to calculate Eq.
:eq:`E_nl`.

Kernel
------

This section details the evaluation of the NLXC kernel. The kernel
:math:`\phi({\mathbf{r}},{\mathbf{r}}')` as specified by Dion *et al*
[Dion2004]_ is given by (in atomic units):

.. math::

   \phi({\mathbf{r}},{\mathbf{r}}') = \frac{1}{\pi^2}\int_{0}^{\infty}a^2da 
       \int_0^{\infty}b^2db W(a,b) T(\nu(a),\nu(b),\nu'(a),\nu'(b))

where

.. math:: T(w,x,y,z) = \frac{1}{2}\Big[\frac{1}{w+x} + \frac{1}{y+z}\Big]\Big[\frac{1}{(w+y)(x+z)}+\frac{1}{(w+z)(y+x)}\Big],

and

.. math::

   \begin{aligned}
       W(a,b) = 2\Big[ & (3-a^2)b\cos b \sin a \\
                       + & (3-b^2)a\cos a \sin b   \\
                       + & (a^2+b^2-3)\sin a\sin b \\
                       - & 3ab\cos a \cos b \Big]/(a^3b^3),\end{aligned}

and

.. math:: \nu(y) = 1- e^{-\gamma y^2/d^2}; \quad \nu'(y) = 1- e^{-\gamma y^2/d'^2};

where :math:`d=|{\mathbf{r}}-{\mathbf{r}}'|q_0({\mathbf{r}})`,
:math:`d'=|{\mathbf{r}}-{\mathbf{r}}'|q_0(\mathbf{r'})`

Following this chain of logic, it is clear that this the kernel can in
fact be considered as a function only of
:math:`|{\mathbf{r}}-{\mathbf{r}}'|`, :math:`q_0({\mathbf{r}})` and
:math:`q_0({\mathbf{r}}')`, since all other variables are dummy
variables which are integrated over. The kernel can therefore be written
as

.. math::
   :label: phi_tab
	   
   \phi(r,q_0({\mathbf{r}}),q_0({\mathbf{r}}'))

This makes it possible to evaluate the integrals above so as to
tabulate the kernel values numerically for a pre-chosen set of radial
points and :math:`q_0` values.

Non-local potential
-------------------

Starting from :eq:`E_nl`, one can evaluate the potential
:math:`v^{\mathrm{nl}}({\mathbf{r}})` corresponding to this energy, by
evaluating all terms in :math:`\partial E_{\mathrm{nl}} /
\partial n({\mathbf{r}})`. The non-local potential
:math:`v_i^{\mathrm{nl}}` at point :math:`{\mathbf{r}}_i` on the grid is
thus written explicitly in terms of the derivatives of the
:math:`\theta_{\alpha}` quantities with respect to the values
:math:`\rho_j` at all other points on the grid:

.. math::
   :label: v_nl
	   
   v_i^{\mathrm{nl}} = \sum_{\alpha}(u_{\alpha i}{\frac{\partial{\theta_{\alpha i}}}{\partial{\rho_i}}}+\sum_j u_{\alpha j}
   {\frac{\partial{\theta_{\alpha j}}}{\partial{\nabla\rho_j}}}{\frac{\partial{\nabla\rho_j}}{\partial{\rho_i}}})

This makes use of the quantities
:math:`u_\alpha({\mathbf{r}})= \sum_{\beta}\mathcal{F}(\theta_{\beta}(\mathbf{k})\phi_{\alpha\beta}(k))`:
which are already calculated in the evaluation of the energy.

Using the White and Bird [White1994]_ approach, Eq.
:eq:`v_nl` can be written as:

.. math::
   :label: v_nl_wnb
	   
   v_{\mathrm{nl}}({\mathbf{r}}) = \sum_{\alpha} \Big( 
     u_{\alpha}({\mathbf{r}}){\frac{\partial{\theta_{\alpha}({\mathbf{r}})}}{\partial{\rho({\mathbf{r}})}}} 
     - \int\int i\mathbf{G}\cdot \frac{\nabla\rho({\mathbf{r}}')}{|\nabla\rho({\mathbf{r}}')|}
     {\frac{\partial{\theta_{\alpha}({\mathbf{r}}')}}{\partial{|\nabla\rho({\mathbf{r}}')|}}}e^{i\mathbf{G}\cdot ({\mathbf{r}}-{\mathbf{r}}')} d{\mathbf{r}}d\mathbf{G}
     \Big)

For this we need to calculate
:math:`{\frac{\partial{\theta}}{\partial{\rho}}}` and
:math:`{\frac{\partial{\theta}}{\partial{{|\nabla{\rho}|}}}}`:

.. math::
   :label: dtheta_drho

   \begin{aligned}
       {\frac{\partial{\theta_\alpha}}{\partial{\rho}}} &=p_\alpha + \rho {\frac{\partial{p_\alpha}}{\partial{\rho}}} \nonumber \\
                                                  &=p_\alpha + \rho {\frac{\partial{p_\alpha}}{\partial{q}}}{\frac{\partial{q}}{\partial{\rho}}} \nonumber \\
                                                  &=p_\alpha + \rho {\frac{\partial{p_\alpha}}{\partial{q}}} \frac{q}{k_F} {\frac{\partial{k_F}}{\partial{\rho}}} + \rho {\frac{\partial{p_\alpha}}{\partial{q}}} k_F ({\frac{\partial{{\varepsilon}_c}}{\partial{\rho}}}{\varepsilon}_x^{-1}-{\varepsilon}_c{\varepsilon}_x^{-2}{\frac{\partial{{\varepsilon}_x}}{\partial{\rho}}} - \frac{8}{3(3\pi^2)^{2/3}}\frac{Z}{4}(\nabla\rho)^2 \rho^{-11/3}) \nonumber \\
                                                  &=p_\alpha + q/3{\frac{\partial{p_\alpha}}{\partial{q}}} + k_F\rho {\frac{\partial{p_\alpha}}{\partial{q}}} ({\frac{\partial{{\varepsilon}_c}}{\partial{\rho}}}{\varepsilon}_x^{-1}-{\varepsilon}_c{\varepsilon}_x^{-2}{\frac{\partial{{\varepsilon}_x}}{\partial{\rho}}}- \frac{2Z}{3(3\pi^2)^{2/3}} {|\nabla{\rho}|}^2\rho^{-11/3})\end{aligned}

.. math::
   :label: dtheta_dgradrho
	   
       {\frac{\partial{\theta_\alpha}}{\partial{{|\nabla{\rho}|}}}} = \rho {\frac{\partial{p_\alpha}}{\partial{q}}} {\frac{\partial{q}}{\partial{{|\nabla{\rho}|}}}} = \frac{Z}{2\rho k_F} \rho {\frac{\partial{p_\alpha}}{\partial{q}}}{|\nabla{\rho}|}

Combining Eqs. :eq:`v_nl_wnb`, :eq:`dtheta_drho` and :eq:`dtheta_dgradrho` gives
us the final expression for the nonlocal potential.

Overview of computational algorithm
===================================

Module ``nlxc_mod``
-------------------

The main module for the calculation of the non-local energy and
potential is ``nlxc_mod``. The tabulation of the kernel :math:`\phi` is
performed only if a kernel file is not found, by ``vdwdf_kernel``.

The input required to calculate the non-local energy and potential is
essentially just the density and its gradient on the fine grid. The
calculation of :math:`q` and the Fourier transformed
:math:`\theta_\alpha` from Eq. :eq:`E_nl` is performed first, in the
routine ``nlxc_q0_theta``. The derivatives of the
:math:`\theta_\alpha`\ s with respect to the density and the module of
its gradient are calculated on-the-fly in the real-space loop for the
calculation of the non-local potential :math:`v_{nl}` in Eq. :eq:`v_nl`. This
is to avoid storing unnecessary arrays. From Eq. :eq:`v_nl_wnb` two
transforms are required per :math:`\alpha` value, a forward FFT,
followed by a backward FFT for calculating the non-local potential.

Subroutines to interpolate the polynomials as well as the kernel using
cubic splines are used (``spline_interp`` and ``interpolate``). The
interpolating polynomials :math:`p_\alpha` used are Kronecker deltas, so
they take the value 1 on the interpolating point and are zero at the
other points.

Module ``vdwdf_kernel``
-----------------------

The kernel :math:`\phi_{\alpha\beta}(k)` is tabulated for 1024 radial
reciprocal space points and 20 :math:`q_0` points. Gaussian quadrature
is used to calculate Eq. :eq:`phi_tab` and then the result is Fourier
transformed. The second derivatives of the kernel are calculated by
interpolation, and also tabulated. The default name of the file is
``vdw_df_kernel``. The program will first check if this file exists: if
it does, it will be loaded in and need not be calculated. If it does
not, it will be generated from scratch (which only takes a few minutes)
and then it is written out for future re-use in the current working
directory.

| The format of the ``vdw_df_kernel`` file is:

````

| ``N_alpha``  ``N_radial``
| ``max_radius``
| ``q_points(:)``
| ``kernel(0:N_radial,alpha,beta)``
| ``kernel2(0:N_radial,alpha,beta)`` 

where ``kernel2`` is the array of second derivatives of the kernel.

[Dion2004] M. Dion, H. Rydberg, E. Schröder, D. C. Langreth, and B. I. Lundqvist,
Phys. Rev. Lett. **92**, 246401 (2004).

[Roman-Perez2009] G. Román-Pérez and J. M. Soler, Phys. Rev. Lett. **103**, 096102 (2009).

[White1994] J. A. White and D. M. Bird, Phys. Rev. B **50**, 4954 (1994).
