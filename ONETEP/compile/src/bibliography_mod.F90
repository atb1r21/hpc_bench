! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-
!=================================================================!
!                                                                 !
!                       Bibliography module                       !
!                                                                 !
! This module keeps track of what modules have been used by the   !
! code for the purposes of outputting a bibliography file at the  !
! end of the run, suggesting citations for work using the results !
! of that run.                                                    !
!-----------------------------------------------------------------!
! This module was created by Nicholas Hine in April 2012.         !
! Largely based on bib.F90 written for the CASTEP code by Matt    !
! Probert.                                                        !
!=================================================================!

module bibliography

  implicit none

  private

  ! Type to store bibliography entries in
  type :: bib_entry
    character (len=20) :: tag        ! lookup
    character (len=20) :: type       ! article or misc
    character (len=1000) :: author
    character (len=300) :: title
    character (len=300) :: journal   ! also doubles up as URL
    character (len=10) :: year
    character (len=10) :: volume     ! char not int as might be range
    character (len=20) :: pages
    character (len=1000) :: abstract ! comment to user as why to cite this paper
  end type bib_entry

  ! Module variables
  integer, parameter :: nbib_max = 100
  logical, save :: bib_init=.false.
  integer :: nbib
  type(bib_entry), save :: bib_library(1:nbib_max)
  logical, save :: bib_cited(1:nbib_max) = .false.

  ! Public subroutines
  public :: bibliography_cite
  public :: bibliography_output

contains

  subroutine bibliography_cite(tag)

    !=========================================================================!
    ! If tag matches an entry in the database then add it to the list of      !
    ! records to output.                                                      !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   tag intent(in): text label to check against database                  !
    !-------------------------------------------------------------------------!
    ! Adapted for ONETEP in April 2012 by Nicholas Hine from a routine        !
    ! written for CASTEP by Matt Probert, 26/01/2011                          !
    !=========================================================================!

    use utils, only: utils_abort

    implicit none

    ! Argument
    character(len=*), intent(in) :: tag

    ! Local variables
    integer :: ibib

    ! Catch first usage and set up database etc
    if (.not.bib_init) call bibliography_init

    ! Check if tag matches an entry
    match_loop: do ibib=1,nbib
       !if (io_uppercase(tag,' -_')==io_uppercase(bib_library(ibib)%tag,' -_')) then
       if (tag==bib_library(ibib)%tag) then
          bib_cited(ibib)=.true.
          exit match_loop
       end if
    end do match_loop

    ! Check we have found a match somewhere
    if (ibib==nbib+1) call utils_abort('Unrecognized value of tag '//tag// &
         ' in bibliography_cite')

    return

  end subroutine bibliography_cite

  subroutine bibliography_output(output_unit)

    !=========================================================================!
    ! Write out bibtex formatted bibliography based on all matched tags found.!
    !-------------------------------------------------------------------------!
    ! Adapted for ONETEP in April 2012 by Nicholas Hine from a routine        !
    ! written for CASTEP by Matt Probert, 26/01/2011                          !
    !=========================================================================!

    use comms, only: pub_on_root, comms_reduce
    use utils, only: utils_abort

    implicit none

    ! Argument
    integer, intent(in) :: output_unit

    ! Local Variables
    integer :: ibib

    ! Catch error
    if (.not.bib_init) &
         call utils_abort('Error in bibliography_output: module not initialised')

    call comms_reduce('OR',bib_cited,nbib_max)

    ! Only output from master root
    if (pub_on_root) then
       do ibib=1,nbib
          if (bib_cited(ibib)) then

             ! Write formatted bibtex entry to output_unit
             select case(trim(bib_library(ibib)%Type))
             case ('article')
                write (output_unit,1) trim(bib_library(ibib)%tag)
                write (output_unit,2) trim(bib_library(ibib)%Author)
                write (output_unit,3) trim(bib_library(ibib)%Title)
                write (output_unit,4) trim(bib_library(ibib)%Journal)
                write (output_unit,5) trim(bib_library(ibib)%Year)
                write (output_unit,6) trim(bib_library(ibib)%Volume)
                write (output_unit,7) trim(bib_library(ibib)%Pages)
                write (output_unit,8) trim(bib_library(ibib)%Abstract)
                write (output_unit,9)
                write (output_unit,*) ''
             case ('misc')
                write (output_unit,1) trim(bib_library(ibib)%tag)
                write (output_unit,2) trim(bib_library(ibib)%Author)
                write (output_unit,3) trim(bib_library(ibib)%Title)
                write (output_unit,10) trim(bib_library(ibib)%Journal)
                write (output_unit,8) trim(bib_library(ibib)%Abstract)
                write (output_unit,9)
                write (output_unit,*) ''
             case default
                call utils_abort('Unrecognised Type in bibliography_output &
                     &- only support article and misc so far')
             end select

          end if
       end do
    end if

1   format('@article{',A,',')
2   format('  Author    = {',A,'},')
3   format('  Title     = {',A,'},')
4   format('  Journal   = {',A,'},')
5   format('  Year      = {',A,'},')
6   format('  Volume    = {',A,'},')
7   format('  Pages     = {',A,'},')
8   format('  Abstract  = {',A,'}')
9   format('}')
10  format('  URL       = {',A,'}')

    return

  end subroutine bibliography_output

  subroutine bibliography_init()

    !=========================================================================!
    ! Create internal bibtex database                                         !
    !-------------------------------------------------------------------------!
    ! Adapted for ONETEP in April 2012 by Nicholas Hine from a routine        !
    ! written for CASTEP by Matt Probert, 26/01/2011                          !
    !=========================================================================!

    use utils, only: utils_abort

    implicit none

    !create bibtex data in standard layout:
    !  tag, Type, Author, &
    !& Title, &
    !& Journal, Year, Volume, Pages, &
    !& Abstract (=usage comment)

    bib_library(1) = bib_entry('DFT-HK', 'article', 'Hohenberg, P. and Kohn, W.', &
         & 'Inhomogeneous electron gas', &
         & 'Phys. Rev.', '1964', '136', 'B864-B871', &
         & 'The original DFT paper with Hohenberg-Kohn theorem')

    bib_library(2) = bib_entry('DFT-KS', 'article','Kohn, W. and Sham, L. J.', &
         & 'Self-consistent equations including exchange and correlation effects', &
         & 'Phys. Rev.', '1965', '140', 'A1133-A1138', &
         & 'The Kohn-Sham equations')

    bib_library(3) = bib_entry('non-local recpot', 'article', 'Kleinman, L. and Bylander, D. M.', &
         & 'Efficacious form for model pseudopotentials', &
         & 'Phys. Rev. Lett.', '1982', '48', '1425-1428', &
         & 'Non-local pseudopotentials')

    bib_library(4) = bib_entry('Ceperley-Alder', 'article', 'Ceperley, D. M. and Alder, B. J.', &
         & 'Ground state of the electron gas by a stochastic method', &
         & 'Phys. Rev. Lett.', '1980', '45', '566-569', &
         & 'QMC calculation of homogenous electron gas')

    bib_library(5) = bib_entry('PBE', 'article', 'Perdew, J. P. and Burke, K. and Ernzerhof, M.', &
         & 'Generalized Gradient Approximation Made Simple', &
         & 'Phys. Rev. Lett.', '1996', '77', '3865-3868', &
         & 'The PBE version of the Generalised Gradient Approximation XC functional')

    bib_library(6) = bib_entry('PZ-LDA', 'article', 'Perdew, J. P. and Zunger, A.', &
         & 'Self-interaction correction to density-functional approximations for many-electron systems', &
         & 'Phys. Rev. B', '1981', '23', '5048-5079', &
         & 'The Perdew-Zunger parameterization of the local density approximation')

    bib_library(7) = bib_entry('BFGS', 'article', 'Pfrommer, B. G. and Cote, M. and Louie, S. G. and Cohen, M. L.', &
         & 'Relaxation of crystals with the quasi-{N}ewton method', &
         & 'J. Comput. Phys.', '1997', '131', '233-240', &
         & 'The BFGS variant used in ONETEP for geometry optimisation')

    bib_library(8) = bib_entry('LBFGS', 'article', 'Byrd, R. H. and Nocedal, J. and Schnabel, R. B.', &
         & 'Representations of quasi-Newton matrices and their use in limited memory methods', &
         & 'Math. Prog.', '1994', '63', '129-156', &
         & 'The standard LBFGS reference')

    bib_library(9) = bib_entry('DI', 'article', 'Andzelm, J. and King-Smith, R. D. and Fitzgerald, G.', &
         & 'Geometry optimization of solids using delocalized internal coordinates', &
         & 'Chem. Phys. Lett.', '2001', '335', '321-326', &
         & 'The DI reference for ONETEP geometry optimisation')

    bib_library(10) = bib_entry('SUPERCELL_PHONON', 'article', ' W. Frank, C. Els\"{a}sser, and M. F\"{a}hnle', &
         & 'Ab initio Force-Constant Method for Phonon Dispersions in Alkali Metals,', &
         & 'Phys. Rev. Lett.', '1995', '74', '1791-1794', &
         & 'The reference for the supercell/finite displacement phonon method')

    bib_library(11) = bib_entry('ONETEP', 'article', 'C.-K. Skylaris, P. D. Haynes, A. A. Mostofi, M. C. Payne', &
         & 'Introducing {ONETEP}: Linear-scaling density functional simulations on parallel computers', &
         & 'J. Chem. Phys.', '2005', '122', '084119', &
         & 'The standard ONETEP reference')

    bib_library(12) = bib_entry('NGWFS', 'article', 'C.-K. Skylaris, A. A. Mostofi, P. D. Haynes, O. Dieguez, M. C. Payne.', &
         & 'Nonorthogonal generalized Wannier function pseudopotential  plane-wave method', &
         & 'Phys. Rev. B', '2002', '66', '035119', &
         & 'The reference for psinc function methods based on NGWFs')

    bib_library(13) = bib_entry('PRECOND', 'article', 'A. A. Mostofi, P. D. Haynes, C.-K. Skylaris, M. C. Payne', &
         & 'Preconditioned iterative minimization for linear-scaling electronic structure calculations', &
         & 'J. Chem. Phys.', '2003', '119', '8842', &
         & 'The reference for NGWF gradient KE preconditioning')

    bib_library(14) = bib_entry('TE_LOC', 'article', 'A. A. Mostofi, C.-K. Skylaris, P. D. Haynes, M. C. Payne', &
         & 'Total-energy calculations on a real space grid with localized functions and a plane-wave basis', &
         & 'Comput. Phys. Commun.', '2002', '147', '788', &
         & 'The reference for methods based on the FFT box')

    bib_library(15) = bib_entry('KINETIC', 'article', 'C.-K. Skylaris, A. A. Mostofi, P. D. Haynes,&
         &  C. J. Pickard, M. C. Payne', &
         & 'Accurate kinetic energy evaluation in electronic structure calculations with localized&
         &  functions on real space grids', &
         & 'Comput. Phys. Commun.', '2001', '140', '315', &
         & 'Evaluation of the Kinetic energy in the FFT box')

    bib_library(16) = bib_entry('FORCES', 'article', &
         & 'N. D. M. Hine, M. Robinson, P. D. Haynes, C.-K. Skylaris, M. C. Payne, and A. A. Mostofi', &
         & 'Accurate ionic forces and geometry optimization in linear-scaling density-functional theory with local orbitals', &
         & 'Phys. Rev. B', '2011', '83', '195102', &
         & 'The reference for calculation of forces in ONETEP')

    bib_library(17) = bib_entry('PBCs', 'article', 'N. D. M. Hine, J. Dziedzic, P. D. Haynes, C.-K. Skylaris', &
         & 'Electrostatic Interactions in Finite Systems treated with Periodic Boundary Conditions: Application to &
         &Linear-Scaling Density Functional Theory ', &
         & 'J. Chem. Phys.', '2011', '135', '204103', &
         & 'Electrostatics avoiding the effect of PBCs in ONETEP: Cutoff Coulomb, Martyna-Tuckerman, OBCs')

    bib_library(18) = bib_entry('IS', 'article', 'J. Dziedzic, H. H. Helal, C.-K. Skylaris, A. A. Mostofi, M. C. Payne', &
         & 'Minimal parameter implicit solvent model for ab initio electronic structure calculations ', &
         & 'Europhysics Letters', '2011', '95', '43001', &
         & 'Implicit Solvent calculations with ONETEP')

    bib_library(19) = bib_entry('COND', 'article', 'L. E. Ratcliff, N. D. M. Hine, and P. D. Haynes', &
         & 'Calculating optical absorption spectra for large systems using linear-scaling density functional theory', &
         & 'Phys. Rev. B', '2011', '84', '165131', &
         & 'Conduction NGWF optimisation with ONETEP')

    bib_library(20) = bib_entry('DISP', 'article', 'Q. Hill and C.-K. Skylaris', &
         & 'Including dispersion interactions in the ONETEP program for linear-scaling density functional theory calculations', &
         & 'Proc. R. Soc. A', '2009', '465', '669', &
         & 'Empirical Dispersion Corrections')

    bib_library(21) = bib_entry('DFTU', 'article', 'D. D. O`Regan, N. D. M. Hine, M. C. Payne and A. A. Mostofi', &
         & 'Linear-scaling DFT+U with full local orbital optimization', &
         & 'Phys. Rev. B', '2012', '85', '085107', &
         & 'Linear-Scaling DFT+U')

    bib_library(22) = bib_entry('PROJSC', 'article', 'D. D. O`Regan, N. D. M. Hine, M. C. Payne and A. A. Mostofi', &
         & 'Projector self-consistent DFT+U using non-orthogonal generalized Wannier functions', &
         & 'Phys. Rev. B', '2010', '82', '081102(R)', &
         & 'Projector-self consistency in Linear-Scaling DFT+U')

    bib_library(23) = bib_entry('NBO', 'article', 'L. P. Lee, D. J. Cole, M. C. Payne, and C.-K. Skylaris', &
         & 'Natural Bond Orbital Analysis in Linear-Scaling Density Functional Theory Calculations: application to a &
         &protein-protein interface', &
         & 'in press', '2012', 'xx', 'xx', &
         & 'Natural Bond Orbital Analysis in ONETEP')

    bib_library(24) = bib_entry('LDOS', 'article', 'N. D. M. Hine, P. W. Avraam, P. Tangney and P. D. Haynes', &
         & 'Linear-Scaling Density Functional Theory Simulations of Polar Semiconductor Nanorods', &
         & 'in press', '2012', 'xx', 'xx', &
         & 'Explanation of LDOS calculations in ONETEP')

    bib_library(25) = bib_entry('PAW', 'article', 'N. D. M. Hine', &
         & 'Linear-Scaling Density Functional Theory using the Projector Augmented Wave Method', &
         & 'J. Phys. Condens. Matter', '2017', '29', '024001', &
         & 'Explanation of PAW calculations in ONETEP')

    bib_library(26) = bib_entry('PULAYFORCES', 'article', 'A. Ruiz Serrano, N. D. M. Hine and C.-K. Skylaris', &
         & 'Pulay forces from localized orbitals optimized in situ using a psinc basis set', &
         & 'J. Chem. Phys.', '2012', '136', '234101', &
         & 'Description of Pulay forces in ONETEP')

    bib_library(27) = bib_entry('EMBEDDING', 'article', &
         & 'S. J. Fox., C. Pittock, T. Fox, C. Tautermann, N. Malcolm, and C.-K. Skylaris', &
         & 'Electrostatic embedding in large-scale first principles quantum mechanical calculations on biomolecules ', &
         & 'J. Chem. Phys.', '2011', '135', '224107', &
         & 'Description of Electrostatic Embedding techniques used in ONETEP')

    bib_library(28) = bib_entry('SUBSPACE', 'article', &
         & 'D. D. O`Regan, M. C. Payne and A. A. Mostofi,', &
         & 'Subspace representations in ab initio methods for strongly correlated systems', &
         & 'Phys. Rev. B', '2011', '83', '245124', &
         & 'Description of use of nonorthogonal functions in methods such as DFT+U, DMFT, cDFT etc.')

    bib_library(29) = bib_entry('CHRISTOFFEL', 'article', &
         & 'D. D. O`Regan', &
         & 'Optimised Projections for the Ab Initio Simulation of Large and Strongly Correlated Systems', &
         & 'Springer Theses, Springer, Berlin, Heidelberg', '2012', 'XVI', '225p.', &
         & 'Description of geometric "Christoffel" corrections preserving operators under representation update')

    bib_library(30) = bib_entry('HFX', 'article', &
         & 'J. Dziedzic, J. C. Womack, R. Ali and C.-K. Skylaris', &
         & 'Massively parallel linear-scaling Hartree-Fock exchange and hybrid &
         &exchange-correlation functionals with plane wave basis set accuracy', &
         & 'J. Chem. Phys.', '2021', '155', '224106', &
         & 'Hartree-Fock exchange and hybrid functional calculations with ONETEP')

    bib_library(31) = bib_entry('POLEMB', 'article', 'J. Dziedzic, Y. Mao, Y. &
         &Shao, J. Ponder, T. Head-Gordon, M. Head-Gordon and C.-K. Skylaris', &
         'TINKTEP: A fully self-consistent, mutually polarizable QM/MM approach&
         & based on the AMOEBA force field', 'J. Chem. Phys.','2016','145',&
         '124106','Polarisable embedding and QM/MM with ONETEP')

    bib_library(32) = bib_entry('ANHARM_AND_DMA', 'article', 'V. Vitale, J. &
         &Dziedzic, S. M.-M. Dubois, H. Fangohr and C.-K. Skylaris', &
         'Anharmonic Infrared Spectroscopy through the Fourier Transform of Time&
         & Correlation Function Formalism in ONETEP','JCTC','2015','11', '3321', &
         'Anharmonic spectra and distributed multipole analysis in ONETEP')

! lam81
    bib_library(33) = bib_entry('Precond-EXP', 'article', 'Packwood, D. and Kermode, J. and &
         & Mones, L. and Bernstein, N. and Woolley, J. and Gould, N. and Ortner, C. and Csanyi, G.', &
         & 'A universal preconditioner for simulating condensed phase materials', &
         & 'J. Chem. Phys.', '2016', '144', '164109', &
         & 'Exponential based preconditioner paper')

    bib_library(34) = bib_entry('Precond-FF', 'article', 'Mones, L. and Csanyi, G. and Ortner, C.', &
         & 'Molecular Mechanical Terms based Preconditioners for Geometry Optimisation and &
         & Saddle Point Search of Molecular Systems', &
         & 'J. Chem. Phys', '2017', '0', '0', &
         & 'Force field based preconditioner paper')

    bib_library(35) = bib_entry('Lindh-FF', 'article', 'Lindh, R. and Bernhardsson, A. and Karlstr\"om G. and Malmqvist, P.', &
         & 'On the use of a Hessian model function in molecular geometry optimizations', &
         & 'Chem. phys. lett.', '1995' ,'241', '423', &
         & 'Simple general force field paper')

    bib_library(36) = bib_entry('ELF', 'article', 'Becke, A. D. and Edgecombe, K. E.', &
         & 'A simple measure of electron localization in atomic and molecular systems', &
         & 'J. Chem. Phys.', '1990' ,'92(9)', '5397-5403', &
         & 'The method for electron localisation descriptor, ELF')

    bib_library(37) = bib_entry('ELF_DFT', 'article', 'Savin, A., Jepsen, O., Flad, J., Andersen, O. K., Preuss, H. &
         & and von Schnering, H. G.', &
         & 'Electron Localization in Solid-State Structures of the Elements: the Diamond Structure', &
         & 'Angew. Chem. Int. Ed. Engl.', '1992' ,'31(2)', '187-188', &
         & 'The method for electron localisation descriptor, ELF, in DFT')

    bib_library(38) = bib_entry('LOL', 'article', 'Schmider, H. L. and Becke, A. D.', &
         & 'Chemical content of the kinetic energy density', &
         & 'J. Mol. Struct (Theochem)', '2000' ,'527', '51-61', &
         & 'The method for electron localisation descriptor, LOL')

    bib_library(39) = bib_entry('Electrolyte', 'article', 'Dziedzic J., &
         &Bhandari A., Anton L., Peng C., Womack C. C., Famili M., Kramer D. &
         &and C.-K. Skylaris', 'Practical Approach to Large-Scale Electronic &
         &Structure Calculations in Electrolyte Solutions via Continuum-Embed&
         &ded Linear-Scaling Density Functional Theory', 'J. Phys. Chem. C', &
         '2020', '124', '7860-7872', 'Implicit solvation with Boltzmann ions')

    bib_library(40) = bib_entry('JCP_special_2020', 'article', &
         'Prentice, Joseph C. A. &
         &and Aarons, Jolyon and Womack, James C. and Allen, Alice E. A. and&
         & Andrinopoulos, Lampros and Anton, Lucian and Bell, Robert A. and &
         &Bhandari, Arihant and Bramley, Gabriel A. and Charlton, Robert J. &
         &and Clements, Rebecca J. and Cole, Daniel J. and Constantinescu, G&
         &abriel and Corsetti, Fabiano and Dubois, Simon M.-M. and Duff, Kev&
         &in K. B. and Escartin, Jose Maria and Greco, Andrea and Hill, Quin&
         &tin and Lee, Louis P. and Linscott, Edward and O`Regan, David D. a&
         &nd Phipps, Maximillian J. S. and Ratcliff, Laura E. and Ruiz Serra&
         &no, Alvaro and Tait, Edward W. and Teobaldi, Gilberto and Vitale, &
         &Valerio and Yeung, Nelson and Zuehlsdorff, Tim J. and Dziedzic, Ja&
         &cek and Haynes, Peter D. and Hine, Nicholas D. M. and Mostofi, Ara&
         &sh A. and Payne, Mike C. and Skylaris, Chris-Kriton', &
         'The {ONETEP} linear-scaling density functional theory program', &
         'J. Chem. Phys.', '2020', '152', '174111', &
         'ONETEP reference paper')

    bib_library(41) = bib_entry('Fisicaro_Soft_Sphere', 'article', &
         'Fisicaro G., Genovese L., Andreussi O., Mandal S., Nair N., &
         &Marzari N. and Goedecker S.', &
         'Soft-Sphere Continuum Solvation in Electronic-Structure Calculations', &
         'J. Chem. Theory Comput.', '2017', '13', '3829', &
         'Soft Sphere Cavity Reference')

    bib_library(42) = bib_entry('Alvarez_vdW_Radii', 'article', &
         'S. Alvaerez', &
         'A cartography of the van der Waals territories', &
         'Dalton Trans.', '2013', '42', '8617', &
         'Alvarez vdW radii for softer sphere')

    !bib_library(33) = bib_entry(tag, Type, Author, &
    !& Title, &
    !& Journal, Year, Volume, Pages, &
    !& Abstract)

    ! Keep track of number of entries
    nbib=42
    if (nbib>nbib_max) &
         call utils_abort('Error in bibliography_init - need to increase nbib_max')

    ! Set global flags
    bib_cited=.false.
    bib_init=.true.

    return

  end subroutine bibliography_init

end module bibliography
