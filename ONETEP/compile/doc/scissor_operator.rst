================================
Species Dependent Scissor Shifts
================================

:Author: Nelson Yeung, University of Warwick
:Author: Nicholas Hine, University of Warwick

Scissor Hamiltonian
===================

The scissor hamiltonian allows one to apply species-dependent and
subspace-dependent energy-level shifts to the hamiltonian, which has the
effect of shifting eigenvalues associated with specific layers of the
material. One can separately shift the valence and conduction subspaces
associated with each layer. This also affects the total energy, so must
be applied with great care if you are using the total energy for any
purpose. The idea is that a band-alignment correction can be applied to
the individual layers of non-covalently-bonded layered materials, though
there may well be many other applications as well. It would seem
“unwise” at best to apply this approach to different regions of the same
molecule or solid which are strongly bonded: results would be
unpredictable and likely unphysical. One ideal use would be to correct
the alignment of the band-edges of a layered material heterobilayer so
that the appropriate heterostructure type was realised, for example
straddled-gap rather than broken-gap, using shifts chosen by reference
to beyond-DFT accuracy calculations of the individual materials, or from
experimental techniques such as ARPES.

In order to apply this shift, we define a scissor Hamiltonian operator
as follows:

.. math::

   \hat{H}_\text{scissor} = \lvert\phi_\eta\rangle K_\text{shifted}^{\eta\delta}
       \langle\phi_\delta\rvert \, ,

where :math:`K_\text{shifted}` is the sum of the species-dependent
shifted valence and conduction density kernel, which is defined as

.. math:: K_\text{shifted} = \sum_L \left( \sigma_{v, L} K_L + \sigma_{c, L} (S^{-1}_L - K_L) \right)\,.

The :math:`\sigma_v` and :math:`\sigma_c` are the shifts for valence
and conduction states, respectively, and the sum is over layer
:math:`L`. The scissor shifted eigenvalues are simply

.. math:: H_\text{shifted} = H + S K_\text{shifted} S \, .

The gradient of the scissor energy with respect to the NGWFs can be
calculated using

.. math::

   \begin{aligned}
       \begin{split}
           \frac{\partial E_\text{scissor}}{\partial\phi_\gamma^*(\mathbf{r})}
           &=
           K_n^{\beta\alpha} \frac{\partial}{\partial\phi_\gamma^*(\mathbf{r})}
           \langle \phi_\alpha\rvert\hat{H}_\text{scissor}\lvert\phi_\beta\rangle \\
           &=
           K_n^{\beta\alpha} \frac{\partial}{\partial\phi_\gamma^*(\mathbf{r})}
           \langle\phi_\alpha\rvert\phi_\eta\rangle K_\text{shifted}^{\eta\delta}
           \langle\phi_\delta\rvert\phi_\beta\rangle \\
           &=
           K_\text{shifted}^{\beta\delta} S_{\delta\eta} K_n^{\eta\alpha} \left(
               \phi_\beta + \sum_{ij} \tilde{p}^i(\mathbf{r}) O_{ij} {R^j}_\alpha
           \right) \, .
       \end{split}\end{aligned}

The approach has been reasonably well tested in the context of LNV
calculations. As of April 2019 it has not been validated for EDFT,
conduction NGWF optimisation, TDDFT etc, but these would be reasonably
expected to work as well.

Performing a Species Dependent Scissor Calculation
==================================================

To activate the shift, the ``species_scissor`` block must be present.
The user must specify (on separate lines) groups of atom types, with
each line finishing with two numbers representing the valence and
conduction shifts to be applied to that group of species.

For example, for a :math:`\textrm{MoS}_2` / :math:`\textrm{MoSe}_2`
heterobilayer, we might use the following to correct the energies of the
individual layers:

::

        %block species_scissor
        Mo1 S  -0.5 0.5
        Mo2 Se -0.1 0.8
        %endblock species_scissor

This would have the effect of opening the gap of the
:math:`\textrm{MoS}_2` layer by 1eV and opening the gap of the
:math:`\textrm{MoSe}_2` layer by 0.9 eV, and shifting the alignment of
the valence bands by 0.4 eV.
