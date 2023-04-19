!-*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-
module sparse

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!   The module in this file was originally created by Nicholas Hine in May    !
!   2009, based largely on the existing SPAM2 code, written by Peter Haynes   !
!   and modified by Nicholas Hine, between 2004 and 2009.                     !
!                                                                             !
!   TCM Group, Cavendish laboratory, University of Cambridge                  !
!   Madingley Road, Cambridge CB3 0HE, UK                                     !
!                                                                             !
!   Subsequent modifications and contributions by Nicholas Hine,              !
!   Jacek Dziedzic, Jose M Escartin, Andrea Greco and Robert Charlton.        !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!=============================================================================!
!                                                                             !
! Block sparse matrix module - documentation                                  !
!                                                                             !
!-----------------------------------------------------------------------------!
!                                                                             !
! The SPAM3 type and its associated module are an update to the SPAM2 type on !
! which ONETEP v2.0-2.3 relied, to improve both the performance and the       !
! capabilities of the matrix algebra. Most significantly, SPAM3 allows        !
! matrices of arbitrary (non-square) dimensions, for applications such as     !
! < NGWF | projector > type overlaps, < NGWF | Hubbard U projector >          !
! matrices, and many others. The new code also represents a redesign          !
! of the data structures for even more efficient parallelisation.             !
! Whereas in SPAM2, matrices were always either wholly-sparse or wholly-dense !
! SPAM3 divides the matrix into "segments", which are blocks corresponding to !
! the sections of the columns belonging to a given proc associated with rows  !
! of a second proc. Segments are assigned sparse or dense depending on the    !
! number of nonzero elements in that segment.                                 !
!                                                                             !
! Sparse matrices in ONETEP are divided column-wise over the different procs  !
! on which the simulation is running, attempting to share out the NGWFs       !
! equally over the procs. A given proc has a certain set of atoms (and thus   !
! NGWFs) associated with it, and only stores the data corresponding to those  !
! NGWF columns. When representing non-square matrices such as an              !
! NGWF-Projector overlap matrix, the same atom-blocking scheme is retained,   !
! with the sizes of the blocks modified accordingly.                          !
!                                                                             !
! The block sparse matrices in ONETEP arise because the structure of the      !
! sparse matrices containing matrix elements between NGWFs (the localised     !
! functions) depends upon whether the spherical atom-centred regions overlap. !
! The matrix elements are therefore grouped into blocks whose size is         !
! determined by the number of NGWFs on each atom.                             !
!                                                                             !
! Imagine a sparse matrix with the structure the overlap matrix of a          !
! simplified linear butane molecule, with 1 NGWF on each hydrogen atom and 4  !
! on each carbon atom                                                         !
!                                                                             !
!            H3   H6   H9   H12                                               !
!            |    |    |    |                                                 !
!       H1 - C2 - C5 - C8 - C11- H14                                          !
!            |    |    |    |                                                 !
!            H4   H7   H10  H13                                               !
!                                                                             !
! Assume the NGWF radii are small enough that only atoms associated with      !
! neighbouring carbon atoms overlap (for illustrative purposes only).         !
! Imagine simulating this on 4 procs of a parallel computer, numbered 0,1,2,3,!
! with 4,3,3, and 4 atoms respectively on them. The sparsity pattern of the   !
! resulting matrix would look something like this, with X's representing      !
! nonzero elements of the overlap matrix:                                     !
!                                                                             !
!  |--------------------||------------------------------------------||        !
!  |                Proc|| 0        || 1      || 2      || 3        ||        !
!  |--------------------||----------++--------++--------++----------||        !
!  |    |    |     |Atom||1| 2  |3|4|| 5  |6|7|| 8  |910|| 11 121314||        !
!  |    |    |     |Type||H|CCCC|H|H||CCCC|H|H||CCCC|H|H||CCCC|H|H|H||        !
!  |Proc|Atom| Type|NGWF||1|1234|1|1||1234|1|1||1234|1|1||1234|1|1|1||        !
!  |--------------------||----------++--------++--------++----------||        !
!  |    | 1  | H   | 1  ||X|XXXX|X|X||OOOO|O|O||OOOO|O|O||OOOO|O|O|O||        !
!  |    |    |     |    ||----------++--------++--------++----------||        !
!  |    |    |     | 1  ||X|XXXX|X|X||XXXX|X|X||OOOO|O|O||OOOO|O|O|O||        !
!  |    | 2  | C   | 2  ||X|XXXX|X|X||XXXX|X|X||OOOO|O|O||OOOO|O|O|O||        !
!  | 0  |    |     | 3  ||X|XXXX|X|X||XXXX|X|X||OOOO|O|O||OOOO|O|O|O||        !
!  |    |    |     | 4  ||X|XXXX|X|X||XXXX|X|X||OOOO|O|O||OOOO|O|O|O||        !
!  |    |    |     |    ||----------++--------++--------++----------||        !
!  |    | 3  | H   | 1  ||X|XXXX|X|X||XXXX|X|X||OOOO|O|O||OOOO|O|O|O||        !
!  |    |    |     |    ||----------++--------++--------++----------||        !
!  |    | 4  | H   | 1  ||X|XXXX|X|X||XXXX|X|X||OOOO|O|O||OOOO|O|O|O||        !
!  |--------------------||==========++========++========++==========||        !
!  |    |    |     | 1  ||O|XXXX|X|X||XXXX|X|X||XXXX|X|X||OOOO|O|O|O||        !
!  |    | 5  | C   | 2  ||O|XXXX|X|X||XXXX|X|X||XXXX|X|X||OOOO|O|O|O||        !
!  | 1  |    |     | 3  ||O|XXXX|X|X||XXXX|X|X||XXXX|X|X||OOOO|O|O|O||        !
!  |    |    |     | 4  ||O|XXXX|X|X||XXXX|X|X||XXXX|X|X||OOOO|O|O|O||        !
!  |    |    |     |    ||----------++--------++--------++----------||        !
!  |    | 6  | H   | 1  ||O|XXXX|X|X||XXXX|X|X||XXXX|X|X||OOOO|O|O|O||        !
!  |    |    |     |    ||----------++--------++--------++----------||        !
!  |    | 7  | H   | 1  ||O|XXXX|X|X||XXXX|X|X||XXXX|X|X||OOOO|O|O|O||        !
!  |--------------------||==========++========++========++==========||        !
!  |    |    |     | 1  ||O|OOOO|O|O||XXXX|X|X||XXXX|X|X||XXXX|X|X|O||        !
!  |    | 8  | C   | 2  ||O|OOOO|O|O||XXXX|X|X||XXXX|X|X||XXXX|X|X|O||        !
!  |    |    |     | 3  ||O|OOOO|O|O||XXXX|X|X||XXXX|X|X||XXXX|X|X|O||        !
!  | 2  |    |     | 4  ||O|OOOO|O|O||XXXX|X|X||XXXX|X|X||XXXX|X|X|O||        !
!  |    |    |     |    ||----------++--------++--------++----------||        !
!  |    | 9  | H   | 1  ||O|OOOO|O|O||XXXX|X|X||XXXX|X|X||XXXX|X|X|O||        !
!  |    |    |     |    ||----------++--------++--------++----------||        !
!  |    | 10 | H   | 1  ||O|OOOO|O|O||XXXX|X|X||XXXX|X|X||XXXX|X|X|O||        !
!  |--------------------||==========++========++========++==========||        !
!  |    |    |     | 1  ||O|OOOO|O|O||OOOO|O|O||XXXX|X|X||XXXX|X|X|X||        !
!  |    | 11 | C   | 2  ||O|OOOO|O|O||OOOO|O|O||XXXX|X|X||XXXX|X|X|X||        !
!  |    |    |     | 3  ||O|OOOO|O|O||OOOO|O|O||XXXX|X|X||XXXX|X|X|X||        !
!  | 3  |    |     | 4  ||O|OOOO|O|O||OOOO|O|O||XXXX|X|X||XXXX|X|X|X||        !
!  |    |    |     |    ||----------++--------++--------++----------||        !
!  |    | 12 | H   | 1  ||O|OOOO|O|O||OOOO|O|O||XXXX|X|X||XXXX|X|X|X||        !
!  |    |    |     |    ||----------++--------++--------++----------||        !
!  |    | 13 | H   | 1  ||O|OOOO|O|O||OOOO|O|O||XXXX|X|X||XXXX|X|X|X||        !
!  |    |    |     |    ||----------++--------++--------++----------||        !
!  |    | 14 | H   | 1  ||O|OOOO|O|O||OOOO|O|O||OOOO|O|O||XXXX|X|X|X||        !
!  |--------------------||----------++--------++--------++----------||        !
!  |--------------------||------------------------------------------||        !
!                                                                             !
! We can see that, for example, the columns stored on proc 0 have no nonzero  !
! elements in rows associated with procs 2 and 3, whereas it has some but not !
! all nonzero elements corresponding to rows associated with proc 1 and every !
! nonzero element corresponding to rows stored on row 0. Putting this into a  !
! matrix of size Nprocs x Nprocs, where 'D' represents a segment with dense   !
! storage, 'S' represents a segment with sparse storage and 'B' is a blank    !
! segment, we have                                                            !
!                                                                             !
! ++-+-------++                                                               !
! ||N|0|1|2|3||                                                               !
! ++-+-+-+-+-++                                                               !
! ||0|D|S|B|B||                                                               !
! ||1|S|D|D|B||                                                               !
! ||2|B|D|D|S||                                                               !
! ||3|B|B|S|D||                                                               !
! +++--------++                                                               !
!                                                                             !
! For the sparse segments, the segment is split into 'blocks', which are      !
! rectangular sections of the matrix corresponding to the overlap of all the  !
! elements associated with one atom with all the elements associated with     !
! another (or itself). The sizes of the blocks are determined by the number   !
! of elements on each atom, be they NGWFs, projectors, etc.                   !
!                                                                             !
! Matrices arising in ONETEP (in the limit of large system sizes) will contain!
! many blocks of zeroes, and therefore it is more efficient to index atomic   !
! blocks rather than individual matrix elements. However, if all the blocks   !
! within a given segment are nonzero, then the indexing is unnecessary and    !
! this part of the matrix may should be stored as a dense segment. Equally, if!
! there are no blocks within a given segment, there is no need to index them. !
!                                                                             !
!-----------------------------------------------------------------------------!
!                                                                             !
! Whether in sparse or dense format, the elements are stored in               !
! column-indexed format. This is because in Fortran, the first array index    !
! changes fastest as one steps through memory i.e. the order in which         !
! elements are stored for an array a(n,m) is:                                 !
!                                                                             !
! a(1,1) | a(2,1) | a(3,1) | ... | a(n,1) | a(1,2) | a(2,2) | ... | a(n,2) ...!
!                                                                             !
! The first index conventionally refers to the rows, and the second to the    !
! columns, of a matrix i.e. a matrix A may be represented by an array of      !
! dimension 2 where the elements are stored as:                               !
!                                                                             !
!    A    <->  a(i,j)                                                         !
!     ij                                                                      !
!                                                                             !
! Thus stepping through memory corresponds to moving down a column.           !
!                                                                             !
! For efficiency, one wishes to address memory sequentially if possible. So   !
! for a simple operation such as axpy one would code as follows:              !
!                                                                             !
!    A    := A   + alpha B                                                    !
!     ij      ij          ij                                                  !
!                                                                             !
!    do j=1,m                                ! loop over columns              !
!       do i=1,n                             ! loop over rows in column j     !
!          a(i,j) = a(i,j) + alpha * b(i,j)  ! element-wise axpy              !
!       end do                               ! loop over rows in column j     !
!    end do                                  ! loop over columns              !
!                                                                             !
! Matrix-matrix multiplication is more complicated. In Fortran maths libraries!
! matrices are represented by arrays in this column-wise fashion.             !
!                                                                             !
!-----------------------------------------------------------------------------!
!                                                                             !
! Sparse segments are represented using an atom-blocked column-indexed sparse !
! storage method. Two integer arrays and one double precision real array are  !
! used to store the matrices as follows:                                      !
!                                                                             !
! Let the number of block-rows and block-columns (i.e. atoms) be nblk, the    !
! number of segments with any nonzero blocks in be nseg, and the number of    !
! nonzero blocks to be stored be nzb.                                         !
!                                                                             !
!  integer :: seg_idx(0:nprocs)          ! Starting indices for each segment  !
!  integer :: blk_idx((nblk+1)*nseg+nzb) ! an index of nonzero blocks         !
!  integer :: seg_ptr(0:nprocs)          ! Pointers to each segment's data    !
!  integer :: blk_ptr((nblk+1)*nseg+nzb) ! pointers to the data in the blocks !
!  real(kind=DP) :: dmtx(nze)        ! the data (matrix elements) themselves  !
!                                                                             !
! The index array works as follows. The segment index points to the starting  !
! indices of each segment, and then the first nblk+1 entries for each         !
! index segment describe the index itself i.e. they point to a list of        !
! nonzero block-rows in each block-column for that segment.                   !
! This list begins at blk_idx(nblk+2) and there is one entry in this list for !
! each of the nzb nonzero blocks in the matrix. The list entries themselves   !
! are the indexes of the nonzero block-rows in a given block-column. For the  !
! carbon dioxide example above we have nprocs=4 (0-3), and on proc 0, we have !
! nblk=4 and nzb=25, with the '0' segment having 16 nonzero blocks and the    !
! '1' segment having 9. Hence, the segment index on proc 0 would be:          !
!                                                                             !
! seg_idx(0)  =                                 1                             !
! seg_idx(1)  = seg_idx(0) + nblk+1 + nzb(0) = 22                             !
! seg_idx(2)  = seg_idx(1) + nblk+1 + nzb(1) = 36                             !
! seg_idx(3)  = seg_idx(2) (as 2 is blank)   = 36                             !
! seg_idx(4)  = seg_idx(3) (as 3 is blank)   = 36                             !
!                                                                             !
!                                                                             !
! blk_idx(1)  = nblk+2 = 6      ! start of list of nonzero block-rows in col 1!
! blk_idx(2)  = 10              ! start of list of nonzero block-rows in col 2!
! blk_idx(3)  = 14              ! start of list of nonzero block-rows in col 3!
! blk_idx(4)  = 18              ! start of list of nonzero block-rows in col 4!
! blk_idx(5)  = nblk+nzb+2 = 22 ! end+1 of list of nonzero block-rows in col 4!
! blk_idx(6)  = 1               ! index of first nonzero block-row in col 1   !
! blk_idx(7)  = 2               ! index of second nonzero block-row in col 1  !
! blk_idx(8)  = 3               ! index of third nonzero block-row in col 1   !
! blk_idx(9)  = 4               ! index of fourth nonzero block-row in col 1  !
! blk_idx(10) = 1               ! index of first nonzero block-row in col 2   !
! blk_idx(11) = 2               ! index of second nonzero block-row in col 2  !
! ...                                                                         !
! blk_idx(21) = 4               ! index of fourth nonzero block-row in col 4  !
! blk_idx(22) = 27              ! start of list of nonzero block-rows in col 1!
! blk_idx(23) = 27              ! start of list of nonzero block-rows in col 2!
! blk_idx(24) = 30              ! start of list of nonzero block-rows in col 3!
! blk_idx(25) = 33              ! start of list of nonzero block-rows in col 4!
! blk_idx(26) = 36              ! end+1 of list of nonzero block-rows in col 4!
! blk_idx(27) = 5               ! index of first nonzero block-row in col 2   !
! blk_idx(28) = 6               ! index of second nonzero block-row in col 2  !
! blk_idx(29) = 7               ! index of third nonzero block-row in col 2   !
! blk_idx(30) = 5               ! index of first nonzero block-row in col 3   !
! blk_idx(31) = 6               ! index of second nonzero block-row in col 3  !
! ...                                                                         !
! blk_idx(35) = 7               ! index of third nonzero block-row in col 4   !
!                                                                             !
! Thus one can straightforwardly loop over all nonzero blocks in a given      !
! column of the whole matrix as follows:                                      !
!                                                                             !
! do seg=0,nprocs-1                                                           !
!   seg_start = seg_idx(seg)                                                  !
!   do blk_col=1,nblk                          ! loop over all columns on proc!
!      do idx=blk_idx(seg_start+blk_col-1), &  ! loop along list of nonzero   !
!             blk_idx(seg_start+blk_col)-1     ! block-rows in block-column   !
!                                              ! blk_col for this segment     !
!         blk_row = blk_idx(idx)               ! index of block-row           !
!         ...                                                                 !
!                                                                             !
! Thus the number of non-zero block-rows in block-column blk_col in segment   !
! seg is given by:                                                            !
!    blk_idx(seg_idx(seg)+blk_col) - blk_idx(seg_idx(seg)+blk_col-1)          !
! and this is why blk_idx(seg_idx(seg)+nblk) is set to nblk+nzb+2.            !
!                                                                             !
! The matrix elements in each nonzero block are stored sequentially as        !
! matrices in the data array in the same order in memory as they would be if  !
! stored as dense two-dimensional arrays **for each segment** as in section 1 !
! above. The blk_ptr points to the start of each block's matrix elements      !
! in the data array. The first nblk entries of blk_ptr contain additional     !
! pointers to the diagonal blocks. The entry blk_ptr(nblk+1) is unused. There !
! is also an extra pointer at the end to mark the end of the data. For the    !
! above example this works out as follows:                                    !
!                                                                             !
! blk_ptr(1)  = 1               ! pointer to diagonal block (1,1)             !
! blk_ptr(2)  = 9               ! pointer to diagonal block (2,2)             !
! blk_ptr(3)  = 41              ! pointer to diagonal block (3,3)             !
! blk_ptr(4)  = 49              ! pointer to diagonal block (4,4)             !
! blk_ptr(5)  = 0               ! unused                                      !
! blk_ptr(6)  = 1               ! pointer to 1x1 array for block (1,1)        !
! blk_ptr(7)  = 2               ! pointer to 1x4 array for block (2,1)        !
! blk_ptr(8)  = 6               ! pointer to 1x1 array for block (3,1)        !
! blk_ptr(9)  = 7               ! pointer to 1x1 array for block (4,1)        !
! blk_ptr(10) = 8               ! pointer to 4x1 array for block (1,2)        !
! blk_ptr(11) = 12              ! pointer to 4x4 array for block (2,2)        !
! ...                                                                         !
! blk_ptr(21) = 50              ! pointer to end+1 of data for block (4,4)    !
! blk_ptr(22) = 0               ! unused                                      !
! blk_ptr(23) = 0               ! unused                                      !
! blk_ptr(24) = 0               ! unused                                      !
! blk_ptr(25) = 0               ! unused                                      !
! blk_ptr(26) = 50              ! pointer to 4x4 array for block (5,2)        !
! blk_ptr(27) = 66              ! pointer to 4x1 array for block (6,2)        !
! blk_ptr(28) = 70              ! pointer to 4x1 array for block (7,2)        !
! blk_ptr(30) = 74              ! pointer to 1x4 array for block (5,3)        !
! blk_ptr(31) = 78              ! pointer to 1x1 array for block (6,3)        !
! ...                                                                         !
! blk_ptr(35) = 86              ! pointer to end+1 of data array for this proc!
!                                                                             !
! The columns of the  sparse matrix on this proc thus contain 85 nonzero      !
! elements and this is the length of the dmtx array nze on this proc.         !
!                                                                             !
! Note that on procs beyond proc zero, the list of diagonal blocks would      !
! appear later, in the section of the array corresponding to that segment of  !
! that proc.                                                                  !
!                                                                             !
!-----------------------------------------------------------------------------!


  use constants, only: DP, LONG, stdout
  use parallel_strategy, only: PARAL_INFO

  implicit none

  private

  ! Public subroutines and functions

  ! sparse_init routines
  public :: sparse_count_ss
  public :: sparse_index_ss
  public :: sparse_count_union
  public :: sparse_index_union

  ! sparse_base routines
  public :: sparse_create
  public :: sparse_destroy
  public :: sparse_index_length
  public :: sparse_generate_index
  public :: sparse_transpose_structure
  public :: sparse_get_element
  public :: sparse_put_element
  public :: sparse_get_block
  public :: sparse_put_block
  public :: sparse_get_col
  public :: sparse_put_col
  public :: sparse_clr_col

  ! sparse_inquiry routines
  public :: sparse_rms_element
  public :: sparse_max_abs_element
  public :: sparse_element_exists
  public :: sparse_fill_fac_denom
  public :: sparse_is_dense
  public :: sparse_num_element
  public :: sparse_proc_num_element
  public :: sparse_num_rows
  public :: sparse_num_cols
  public :: sparse_first_elem_on_proc
  public :: sparse_last_elem_on_proc
  public :: sparse_first_elem_on_atom
  public :: sparse_num_elems_on_proc
  public :: sparse_num_elems_on_atom
  public :: sparse_atom_of_elem
  public :: sparse_proc_of_elem
  public :: sparse_memory
  public :: sparse_show_memory_usage
  public :: sparse_max_imag_part
  public :: sparse_check_hermitian
  public :: sparse_get_par
  public :: sparse_rms_element_array

  ! sparse_ops routines
  public :: sparse_copy
  public :: sparse_scale
  public :: sparse_conjugate
  public :: sparse_axpy
  public :: sparse_product
  public :: sparse_trace
  public :: sparse_transpose
  public :: sparse_expand
  public :: sparse_extremal_eigenvalue
  public :: sparse_extremal_eigenvalue_array
  public :: sparse_outer_product
  public :: sparse_1norm
  public :: sparse_any_isnan
  public :: sparse_entrywise_norm
  public :: sparse_take_real_part
  public :: sparse_safe_axpy

  ! sparse_algor routines
  public :: sparse_solve2
  public :: sparse_hotelling_init
  public :: sparse_hotelling_invert

  ! sparse_utils routines
  public :: sparse_init_blocking_scheme
  public :: sparse_mod_init
  public :: sparse_exit
  public :: sparse_convert
#ifdef SCALAPACK
  public :: sparse_spam3toblacs
  public :: sparse_blacstospam3
#endif
  public :: sparse_show_matrix
  public :: sparse_show_segment_filling
  public :: sparse_show_network
  public :: sparse_write
  public :: sparse_read
  public :: sparse_convert_unsegment_real
  public :: sparse_convert_segment_real
  public :: sparse_convert_unsegment_complex
  public :: sparse_convert_segment_complex
  public :: sparse_set_to_garbage
  public :: sparse_diagnose
  public :: sparse_mat_in_library

  ! Enquiry routines
  interface sparse_num_element
     module procedure sparse_num_element_lib
     module procedure sparse_num_element_mat
     !module procedure sparse_num_element_mat_embed
  end interface
  interface sparse_num_elems_on_atom
     module procedure sparse_num_elems_on_atom_mat
     module procedure sparse_num_elems_on_atom_blks
  end interface

  ! Element operation routines
  interface sparse_get_element
     module procedure sparse_get_element_real
     module procedure sparse_get_element_complex
     module procedure sparse_get_element_coef
  end interface
  interface sparse_put_element
     module procedure sparse_put_element_real
     module procedure sparse_put_element_complex
     module procedure sparse_put_element_coef
  end interface

  ! Block operation routines
  interface sparse_get_block
     module procedure sparse_get_block_real
     module procedure sparse_get_block_complex
  end interface
  interface sparse_put_block
     module procedure sparse_put_block_real
     module procedure sparse_put_block_complex
  end interface

  ! Column operation routines
  interface sparse_get_col
     module procedure sparse_get_col_real
     module procedure sparse_get_col_complex
  end interface
  interface sparse_put_col
     module procedure sparse_put_col_real
     module procedure sparse_put_col_complex
  end interface
  interface sparse_clr_col
     module procedure sparse_clr_col_real
     module procedure sparse_clr_col_complex
  end interface


  ! AXPY routines
  interface sparse_axpy
     module procedure sparse_axpy_real
     module procedure sparse_axpy_complex
  end interface

  ! agrecocmplx
  ! safe AXPY routines
  interface sparse_safe_axpy
     module procedure sparse_safe_axpy_real
     module procedure sparse_safe_axpy_complex
  end interface

  ! Scale and shift routines
  interface sparse_scale
     module procedure sparse_scale_real
     module procedure sparse_scale_complex
     module procedure sparse_scale_complex2
  end interface

  ! Outer product routines
  interface sparse_outer_product
     module procedure sparse_outer_product_real
     module procedure sparse_outer_product_complex
  end interface

  ! Conversion routines
  interface sparse_convert
     module procedure sparse_spam3tofull_real
     module procedure sparse_spam3tofull_complex
     module procedure sparse_fulltospam3_real
     module procedure sparse_fulltospam3_complex
  end interface
#ifdef SCALAPACK
  interface sparse_spam3toblacs
     module procedure sparse_spam3toblacs_real
     module procedure sparse_spam3toblacs_complex
  end interface
  interface sparse_blacstospam3
     module procedure sparse_blacstospam3_real
     module procedure sparse_blacstospam3_complex
  end interface
#endif

  ! Read/write routines
  interface sparse_write
     module procedure sparse_write_scalar
     module procedure sparse_write_vector
     module procedure sparse_write_matrix
  end interface
  interface sparse_read
     module procedure sparse_read_scalar
     module procedure sparse_read_vector
     module procedure sparse_read_matrix
  end interface

  ! Type definition for the structure of a block sparse matrix
  type STRUC3

    character(len=30) :: structure ! The structure code for this matrix type
    character(len=30) :: transpose_structure ! The structure code for the
                                             ! transpose of this structure
    integer :: nrows          ! The total number of rows in the matrix
    integer :: mcols          ! The total number of cols in the matrix
    integer :: nblk           ! The total number of block-rows (atoms)
    integer :: mblk           ! The total number of block-columns (atoms)
    integer(kind=LONG) :: nze ! The total number of nonzero elements (<n*m)
    integer :: nzb            ! The total number of nonzero blocks stored
    integer :: my_mcols       ! The number of cols on this proc
    integer :: my_nblks       ! The number of blocks on this proc
    integer :: my_nze         ! The number of nonzero elements on this proc
    integer :: my_nzb         ! The number of nonzero blocks on this proc
    integer :: max_nze        ! Maximum number of nonzero elements on any proc
    integer :: max_nzb        ! Maximum number of nonzero blocks on any proc
    integer :: max_nblks      ! Maximum number of block-cols on any proc
    integer :: col_blks       ! Identifier for column blocking scheme
    integer :: row_blks       ! Identifier for row blocking scheme
    integer,allocatable :: idx_lens(:) ! Lengths of index arrays on all procs
    integer,allocatable :: idx_seg_lens(:) ! Lengths of the segments of idx
         ! array on other procs corresponding to the segment for this proc
    integer,allocatable :: mtx_seg_lens(:) ! Lengths of the segments of data
         ! array on other procs corresponding to the segment for this proc
    integer,allocatable :: idx_groupseg_lens(:) ! Total length of the segments
         ! of idx array on other procs corresponding to procs in this group
    integer,allocatable :: mtx_groupseg_lens(:) ! Total length of the segments
         ! of data on other procs corresponding to procs in this group
    ! Segment information: (s_type,:) elements store whether each segment is
    ! dense, sparse or blank. (s_idx,:) elements store the starting position of
    ! the index for each segment. (s_ptr,:) elements store the pointer to the
    ! start of the data for each segment in the dmtx/zmtx arrays.
    integer,allocatable :: seg_info(:,:)
    integer,allocatable :: blk_idx(:)   ! Column indexed list of nonzero blocks
    integer,allocatable :: blk_ptr(:)   ! Pointers to nonzero blocks
    type(PARAL_INFO), pointer :: row_par ! rc2013: pointer to row parallel strategy
    type(PARAL_INFO), pointer :: col_par ! rc2013: pointer to column parallel strategy
  end type STRUC3

  ! Type definition for a communicator to send matrix data between MPI processes
  type, public :: COM3
    integer :: lib
    integer :: buflen
    integer :: datlen
    integer :: num_handles
    integer :: idx_hdl       ! Index of index recv msg in handles array
    integer :: info_hdl      ! Index of seginfo recv msg in handles array
    integer :: ptr_hdl       ! Index of pointer recv msg in handles array
    integer :: data_hdl      ! Index of data recv msg in handles array
    integer :: req_send_hdl  ! Index of request send msg in handles array
    integer :: data_send_hdl ! Index of data send msg in handles array
    integer :: nreqs
    integer :: reqdatlen
    integer :: nbuf, cbuf, ubuf
    logical :: iscmplx
    logical :: cropped
    logical :: segonly
    logical :: groupsegonly
    logical :: got_index
    logical :: got_reqs
    logical :: got_data
    type(PARAL_INFO), pointer :: row_par ! rc2013: pointer to row parallel strategy
    type(PARAL_INFO), pointer :: col_par ! rc2013: pointer to column parallel strategy
    integer, allocatable :: seginfobuf(:,:,:)
    integer, allocatable :: ptrbuf(:,:)
    integer, allocatable :: idxbuf(:,:)
    real(kind=DP), allocatable :: dmtxrecvbuf(:,:)
    real(kind=DP), allocatable :: dmtxsendbuf(:,:)
    logical, allocatable :: send_buffer_free(:)
    complex(kind=DP), allocatable :: zmtxrecvbuf(:,:)
    complex(kind=DP), allocatable :: zmtxsendbuf(:,:)
    integer, allocatable :: ptrreqrecvbuf(:,:)
    integer, allocatable :: ptrreqsendbuf(:,:)
    integer, allocatable :: index_reqs(:) ! Receive buffer for index requests
    integer, allocatable :: data_reqs(:)  ! Receive buffer for data requests
    integer, allocatable :: handles(:)    ! MPI message handles for each recv
  end type COM3

  ! Type definition for data-sharing storage within a comms group
  type, public :: SHARE3

    integer :: lib
    logical :: iscmplx
    integer, allocatable :: seginfobuf(:,:,:) ! Recv buffer for segment info
    integer, allocatable :: idxbuf(:,:)   ! Recv buffer for index
    integer, allocatable :: ptrbuf(:,:)   ! Recv buffer for ptr list
    type(PARAL_INFO), pointer :: row_par ! rc2013: pointer to row parallel strategy
    type(PARAL_INFO), pointer :: col_par ! rc2013: pointer to column parallel strategy
    real(kind=DP), allocatable :: dmtxbuf(:,:)   ! Buffer for real data
    complex(kind=DP), allocatable :: zmtxbuf(:,:)! Buffer for cplx data

  end type SHARE3

  ! Type definition for the matrix elements themselves, with a pointer to
  ! the appropriate structure
  type, public :: SPAM3
    integer :: lib     ! The index for the structure in the library
    logical :: iscmplx ! TRUE if matrix is complex, otherwise FALSE
    real(kind=DP), allocatable :: dmtx(:)  ! If real, the non-zero elements of
                                           ! the matrix, stored in blocks
    complex(kind=DP), allocatable :: zmtx(:)  ! If complex, the non-zero elmnts
                                           ! of the matrix, stored in blocks
    character(len=30) :: structure         ! Character string identifying the
                                           ! sparse matrix structure
  end type SPAM3

  ! Segment info list components
  integer, parameter :: s_type = 1
  integer, parameter :: s_idx = 2
  integer, parameter :: s_ptr = 3

  ! Segment types
  integer, parameter :: SEG_DENSE      =  222
  integer, parameter :: SEG_SPARSE     =  333
  integer, parameter :: SEG_BLANK      = -111

  ! Block size list indices
  integer, public, parameter :: BLKS_NGWF      =  1
  integer, public, parameter :: BLKS_PROJ      =  2
  integer, public, parameter :: BLKS_HUB_PROJ  =  3
  integer, public, parameter :: BLKS_COND      =  4
  integer, public, parameter :: BLKS_JOINT     =  5
  integer, public, parameter :: BLKS_SW        =  6
  integer, public, parameter :: BLKS_CORE      =  7
  integer, public, parameter :: BLKS_AUX       =  8
  integer, public, parameter :: BLKS_PDOS      =  9
  integer, public, parameter :: BLKS_JOINTPH   =  10 ! gcc32

  ! Element patterns for sparse_expand
  integer, public, parameter :: PATTERN_LOWER     = 444
  integer, public, parameter :: PATTERN_ALTERNATE = 555
  integer, public, parameter :: PATTERN_FULL      = 666

  ! Allow creation of new matrix structures?
  logical, public :: pub_sparse_allow_new_matrices

  ! Tag identifiers and special send values for sparse_product
  integer, parameter :: SEGINFO_TAG  = 200000001
  integer, parameter :: BLKIDX_TAG   = 200000002
  integer, parameter :: BLKPTR_TAG   = 200000003
  integer, parameter :: DATA_TAG     = 200000004
  integer, parameter :: IDX_REQ_TAG  = 500000001
  integer, parameter :: DATA_REQ_TAG = 500000002
  integer, parameter :: PTR_REQ_TAG  = 500000003
  integer, parameter :: index_sent = -444
  integer, parameter :: index_not_sent = -555
  integer, parameter :: index_needed = -666
  integer, parameter :: data_sent = -777
  integer, parameter :: data_not_sent = -888

  ! Library of known sparse matrix structures
  ! rc2013: for testing let's make the library bigger
  integer, parameter :: max_library = 2500 ! Maximum number of indices in library
  integer :: num_library                 ! Number of indices in library
  type(STRUC3) :: library(max_library)   ! Library of sparse indices

  ! Record of matrices currently in use
  integer :: nalloc
  integer(kind=LONG) :: global_struc_mem
  integer(kind=LONG) :: local_struc_mem
  integer(kind=LONG) :: global_mat_mem
  integer(kind=LONG) :: local_mat_mem

  ! rc2013: parallelisation variables moved to PARAL_INFO

  ! File version number
  real(kind=DP), parameter :: file_version = 1.0_DP

contains

  subroutine sparse_init_blocking_scheme(scheme,num,num_funcs_on_proc, &
       num_funcs_on_atom, first_func_on_proc, first_func_on_atom, &
       atom_of_func, proc_of_func, par)

    use comms, only: pub_total_num_procs
    use rundat, only: pub_any_nl_proj, pub_PAW, pub_hubbard, &
         pub_cond_calculate, pub_eels_calculate, pub_use_aux_ngwfs, pub_use_swx, &
         & pub_dos_smear, pub_pdos_max_l, pub_pdos_reduce_sws, &
           pub_lr_phonons_calculate, pub_pdos_construct_basis, &
           pub_use_activeswx, pub_active_region
    use parallel_strategy, only: PARAL_INFO
    use utils, only: utils_alloc_check

    implicit none

    ! Arguments
    integer, intent(in) :: scheme
    integer, intent(in) :: num
    integer, intent(in) :: num_funcs_on_proc(0:pub_total_num_procs-1)
    type(PARAL_INFO), intent(inout) :: par ! rc2013: parallel strategy
                      ! associated with this blocking scheme (region).
    integer, intent(in) :: num_funcs_on_atom(1:par%nat)
    integer, intent(in) :: first_func_on_proc(0:pub_total_num_procs)
    integer, intent(in) :: first_func_on_atom(1:par%nat)
    integer, intent(in) :: atom_of_func(1:num)
    integer, intent(in) :: proc_of_func(1:num)

    ! Local Variables
    integer :: ierr
    integer :: num_elem_types
    integer :: max_elems

    ! Allocate arrays for local copies of all blocking information
    if (.not.allocated(par%first_elem_on_atom)) then

       ! Find number of blocking schemes
       num_elem_types = BLKS_NGWF
       if (pub_any_nl_proj.or.pub_paw) num_elem_types = BLKS_PROJ
       if (pub_hubbard) num_elem_types = BLKS_HUB_PROJ
       if (pub_cond_calculate) num_elem_types = BLKS_JOINT
       if (pub_use_swx.or.pub_use_activeswx) num_elem_types = BLKS_SW
       if (pub_eels_calculate) num_elem_types = BLKS_CORE
       if (pub_use_aux_ngwfs) num_elem_types = BLKS_AUX
       if (pub_dos_smear > 0.0_dp .and. pub_pdos_max_l >=0) num_elem_types = BLKS_PDOS
       if (pub_lr_phonons_calculate) num_elem_types = BLKS_JOINTPH

       ! Find max number of elements in any blocking scheme
       max_elems = par%num_ngwfs
       if (pub_any_nl_proj) max_elems = max(max_elems,par%num_projectors)
       if (pub_paw) max_elems = max(max_elems,par%num_pawpws)
       if (pub_hubbard) max_elems = max(max_elems,par%num_hub_proj)
       if (pub_cond_calculate) max_elems = max(max_elems, &
            par%num_ngwfs_cond+par%num_ngwfs)
       if (pub_use_swx .or. (pub_use_activeswx.and.(par%par_index==pub_active_region))) &
            max_elems = max(max_elems, par%num_sw)
       if (pub_eels_calculate) max_elems = max(max_elems,par%num_corewfs)
       if (pub_use_aux_ngwfs) max_elems = max(max_elems,par%num_ngwfs_aux)

       if (pub_dos_smear > 0.0_dp .and. pub_pdos_max_l >= 0) then
          if(pub_pdos_reduce_sws) then
             if(pub_pdos_construct_basis) then
                max_elems = max(max_elems,((pub_pdos_max_l+1)**2)*par%num_ngwfs)
             else
                max_elems = max(max_elems,((pub_pdos_max_l+1)**2)*par%nat)
             end if
          else
             max_elems = max(max_elems,((pub_pdos_max_l+1)**2)*par%nat*par%max_pdos_n)
          end if
       end if

       if (pub_lr_phonons_calculate) max_elems = max(max_elems, 2*par%num_ngwfs)

       ! rc2013: allocate to the parallel strategy
       allocate(par%first_elem_on_atom(par%nat,num_elem_types),stat=ierr)
       call utils_alloc_check('sparse_init','first_elem_on_atom',ierr)

       allocate(par%first_elem_on_proc(0:pub_total_num_procs,num_elem_types), &
            stat=ierr)
       call utils_alloc_check('sparse_init','first_elem_on_proc',ierr)

       allocate(par%num_elems_on_atom(par%nat,num_elem_types),stat=ierr)
       call utils_alloc_check('sparse_init','num_elems_on_atom',ierr)

       allocate(par%num_elems_on_proc(0:pub_total_num_procs-1,num_elem_types), &
            stat=ierr)
       call utils_alloc_check('sparse_init','num_elems_on_proc',ierr)

       allocate(par%atom_of_elem(max_elems,num_elem_types),stat=ierr)
       call utils_alloc_check('sparse_init','atom_of_elem',ierr)

       allocate(par%proc_of_elem(max_elems,num_elem_types),stat=ierr)
       call utils_alloc_check('sparse_init','proc_of_elem',ierr)
    end if

    ! rc2013: allocate to each par
    par%first_elem_on_atom(:,scheme) = first_func_on_atom(:)
    par%first_elem_on_proc(:,scheme) = first_func_on_proc(:)
    par%num_elems_on_atom(:,scheme) = num_funcs_on_atom(:)
    par%num_elems_on_proc(:,scheme) = num_funcs_on_proc(:)
    par%atom_of_elem(1:num,scheme) = atom_of_func(:)
    par%proc_of_elem(1:num,scheme) = proc_of_func(:)

  end subroutine sparse_init_blocking_scheme


  !==========================================================================!
  ! This routine performs standard initialisation tasks for the module, such !
  ! as allocating temporary arrays and setting up counters and shorthands.   !
  !--------------------------------------------------------------------------!
  ! Arguments:                                                               !
  ! par               :: parallel strategy                                   !
  ! first_sparse_call :: check if this is the 1st call to this module in     !
  !                      case of multiple parallel strategies                !
  !--------------------------------------------------------------------------!
  ! Written by Nicholas Hine in April 2011 from bits of sparse_init.         !
  ! Modified for multiple parallel strategies by Robert Charlton, 06/09/2016.!
  !==========================================================================!

  subroutine sparse_mod_init(par, first_sparse_call)

    use comms, only: pub_my_proc_id
    use parallel_strategy, only: PARAL_INFO
    use rundat, only: pub_timings_level
!$  use rundat, only: pub_threads_max
    use timer, only: timer_clock

    implicit none

    ! Arguments
    type(PARAL_INFO), pointer, intent(inout) :: par
    logical, intent(inout)                   :: first_sparse_call

    ! rc2013: initialise all these if this is the first call
    ! rc2013: otherwise simply iterate through the library
    if(first_sparse_call) then
       num_library = 0
       global_mat_mem = 0
       local_mat_mem = 0
       global_struc_mem = 0
       local_struc_mem = 0
       nalloc = 0
       first_sparse_call = .false.
    endif

    ! Set shorthand variables
    ! rc2013: allow for possibility that there are no atoms on some procs
    if(par%num_atoms_on_proc(pub_my_proc_id) == 0) then
       par%my_first_blk = 0
       par%my_last_blk = -1
    else
       par%my_first_blk = par%first_atom_on_proc(pub_my_proc_id)
       par%my_last_blk = par%first_atom_on_proc(pub_my_proc_id + 1) - 1
    endif

!$OMP PARALLEL NUM_THREADS(pub_threads_max) DEFAULT(NONE) &
!$OMP SHARED(pub_timings_level,pub_threads_max)

    ! Start and stop all the timers to avoid messing up the ordering later
    if (iand(pub_timings_level,8) /= 0) then
       call timer_clock('sparse_product',1)
       call timer_clock('sparse_product',2)
       call timer_clock('sparse_product_omp_barrier',1)
       call timer_clock('sparse_product_omp_barrier',2)
       call timer_clock('sparse_product_dense_dgemm',1)
       call timer_clock('sparse_product_dense_dgemm',2)
       call timer_clock('sparse_product_picks',1)
       call timer_clock('sparse_product_picks',2)
       call timer_clock('sparse_product_sparse',1)
       call timer_clock('sparse_product_sparse',2)
       call timer_clock('sparse_product_get_crop_step_data',1)
       call timer_clock('sparse_product_get_crop_step_data',2)
       call timer_clock('sparse_product_crop_indexB',1)
       call timer_clock('sparse_product_crop_indexB',2)
       call timer_clock('sparse_product_crop_indexL',1)
       call timer_clock('sparse_product_crop_indexL',2)
       call timer_clock('sparse_product_check_send_requests',1)
       call timer_clock('sparse_product_check_send_requests',2)
       call timer_clock('sparse_product_init_share',1)
       call timer_clock('sparse_product_init_share',2)
       call timer_clock('sparse_product_exit_comms',1)
       call timer_clock('sparse_product_exit_comms',2)
       call timer_clock('sparse_product_exit_share',1)
       call timer_clock('sparse_product_exit_share',2)
    end if

!$OMP END PARALLEL

  end subroutine sparse_mod_init


  !============================================================================!
  ! This subroutine performs all the finalisation routines involving the       !
  ! sparse matrices.                                                           !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   None                                                                     !
  !----------------------------------------------------------------------------!
  ! Written by Peter Haynes, March 2004.                                       !
  ! Updated for SPAM3 by Nicholas Hine, May 2009.                              !
  !============================================================================!

  subroutine sparse_exit(par)

    use utils, only: utils_dealloc_check

    implicit none

    ! Arguments:
    type(PARAL_INFO), intent(inout) :: par

    ! Local variables
    integer :: ierr       ! Status flag
    integer :: ilibrary   ! Counter over library entries

    ! Loop over library entries
    do ilibrary=num_library,1,-1

       ! Deallocate arrays within matrix structure
       call sparse_struc_dealloc(library(ilibrary))

    end do

    ! Library is empty
    num_library = 0

    ! Deallocate copies of parallel strategy arrays
    deallocate(par%proc_of_elem,stat=ierr)
    call utils_dealloc_check('sparse_exit','proc_of_elem',ierr)

    deallocate(par%atom_of_elem,stat=ierr)
    call utils_dealloc_check('sparse_exit','atom_of_elem',ierr)

    deallocate(par%num_elems_on_proc,stat=ierr)
    call utils_dealloc_check('sparse_exit','num_elems_on_proc',ierr)

    deallocate(par%num_elems_on_atom,stat=ierr)
    call utils_dealloc_check('sparse_exit','num_elems_on_atom',ierr)

    deallocate(par%first_elem_on_proc,stat=ierr)
    call utils_dealloc_check('sparse_exit','first_elem_on_proc',ierr)

    deallocate(par%first_elem_on_atom,stat=ierr)
    call utils_dealloc_check('sparse_exit','first_elem_on_atom',ierr)

  end subroutine sparse_exit

  !==========================================================================!
  ! This routine returns the numbers of matrix elements and blocks (total,   !
  ! per proc and per segment) which are nonzero due to overlap of pairs of   !
  ! spheres, as stored in the overlaps%overlap_list array.                   !
  !--------------------------------------------------------------------------!
  ! Arguments:                                                               !
  !   overlaps(input)  : Overlap list as provided by parallel_strategy_mod   !
  !   nze     (output) : Total number of nonzero elements                    !
  !   nzb     (output) : Total number of nonzero blocks                      !
  !   my_nze  (output) : Number of nonzero elements on this proc             !
  !   my_nzb  (output) : Number of nonzero blocks on this proc               !
  !   seg_nze (output) : Number of nonzero elements per segment on this proc !
  !   seg_nzb (output) : Number of nonzero blocks per segment on this proc   !
  !--------------------------------------------------------------------------!
  ! Written by Peter Haynes, March 2004.                                     !
  ! Revised for distributed data, June 2006.                                 !
  ! Updated for SPAM3 by Nicholas Hine, May 2009.                            !
  !==========================================================================!

  subroutine sparse_count_ss(overlaps,col_blks,row_blks,nze,nzb,my_nze,my_nzb, &
       seg_nze,seg_nzb, row_par, col_par)

    use comms, only: comms_reduce, comms_allgather, pub_total_num_procs
    use parallel_strategy, only: OVERLAP_LIST
    use utils, only: utils_abort, utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(OVERLAP_LIST), intent(inout) :: overlaps
    integer, intent(in) :: col_blks ! Number identifying col blk sizes
    integer, intent(in) :: row_blks ! Number identifying row blk sizes
    integer(kind=LONG), intent(out) :: nze ! Number of nonzero elements
    integer, intent(out) :: nzb     ! Number of nonzero blocks (atoms)
    integer, intent(out) :: my_nze  ! Number of nonzero elements on this proc
    integer, intent(out) :: my_nzb  ! Number of nonzero blocks on this proc
    integer, intent(out) :: seg_nze(0:pub_total_num_procs-1)
    integer, intent(out) :: seg_nzb(0:pub_total_num_procs-1)
    type(PARAL_INFO), intent(in), pointer :: row_par
    type(PARAL_INFO), intent(in), pointer :: col_par

    ! Local variables
    integer :: proc        ! Proc of each found block
    integer :: iblk, jblk  ! Atom counters
    integer :: iovlap      ! Atomic overlap counter
    integer :: ielems      ! Number of elements on atom iblk
    integer, allocatable :: ibuf(:) ! Communications buffer for comms_allgather
    integer :: ierr        ! Error flag for utils_{alloc,dealloc}_check
    real(kind=DP) :: nze_dp     ! Real nze for checking for integer overflow
    real(kind=DP) :: my_nze_dp  ! Real my_nze for checking for integer overflow
    character(len=:), allocatable :: error_string ! For output of error message
    logical :: integer_overflow ! Has an integer overflow been detected?

    ! No integer overflow... yet
    integer_overflow = .false.

    ! Zero counters
    my_nze = 0
    my_nze_dp = real(my_nze,kind=DP)
    my_nzb = 0
    seg_nzb(:) = 0
    seg_nze(:) = 0

    ! Loop over atoms on this proc (block-columns)
    do iblk=col_par%my_first_blk,col_par%my_last_blk

       ! Number of elements on atom iblk
       ielems = col_par%num_elems_on_atom(iblk,col_blks)

       ! Skip this atom if there are no elements on it - this is the case
       ! for atoms which are not Hubbard atoms when constructing V and W
       if (ielems==0) overlaps%num_overlaps(iblk) = 0

       ! Loop over all atoms jblk which overlap atom iblk (block-rows)
       do iovlap=1,overlaps%num_overlaps(iblk)
          jblk = overlaps%overlap_list(iovlap,iblk)

          ! Skip this atom if there are no elements on it
          if (row_par%num_elems_on_atom(jblk,row_blks)==0) cycle

          ! Count contribution to my_nze and my_nzb totals
          my_nze = my_nze + ielems * row_par%num_elems_on_atom(jblk,row_blks)
          my_nzb = my_nzb + 1

          ! Increment real version of my_nze, for integer overflow checking
          my_nze_dp = my_nze_dp + real(ielems,kind=DP) * &
               real(row_par%num_elems_on_atom(jblk,row_blks),kind=DP)

          ! Check for integer overflow
          if (my_nze_dp > HUGE(my_nze)) integer_overflow = .true.

          ! Find proc to which jblk belongs and record in segment nze and nzb
          proc = row_par%proc_of_elem(row_par%first_elem_on_atom(jblk,row_blks),row_blks)

          seg_nze(proc) = seg_nze(proc) + ielems * &
               row_par%num_elems_on_atom(jblk,row_blks)
          seg_nzb(proc) = seg_nzb(proc) + 1

       end do  ! Loop over atoms jblk which overlap atom iblk (block-rows)

    end do  ! Loop over atoms iblk on this proc (block-columns)

    ! Sum up over all procs
    ! JCW: Susceptible to integer overflow, since we have no comms_reduce
    !      procedures for integer(kind=LONG)
    !ibuf(1) = my_nze
    !ibuf(2) = my_nzb
    !call comms_reduce('SUM',ibuf)
    !nze = ibuf(1)
    !nzb = ibuf(2)
    ! JCW: In the absence of a comms_reduce for LONG integers, we avoid integer
    !      overflow by allgathering to an array of integers on each proc, then
    !      locally summing into nze, which is integer(kind=LONG)
    allocate(ibuf(pub_total_num_procs),stat=ierr)
    call utils_alloc_check('sparse_count_ss','ibuf',ierr)
    ! nze
    ibuf(:) = 0
    call comms_allgather(ibuf,my_nze)
    nze = sum(int(ibuf(:),kind=LONG))
    ! nzb
    ibuf(:) = 0
    call comms_allgather(ibuf,my_nzb)
    nzb = sum(ibuf(:))
    deallocate(ibuf,stat=ierr)
    call utils_dealloc_check('sparse_count_ss','ibuf',ierr)

    ! Sum over all procs for integer overflow checking
    nze_dp = my_nze_dp
    call comms_reduce('SUM',nze_dp)
    if (nze_dp > HUGE(nze)) integer_overflow = .true.

    if (integer_overflow) then
       allocate(character(len=80) :: error_string)
       ! len=80, because we want 4*20 width numeric values
       write(error_string,'(4es20.10)') &
            real(nze,kind=DP), nze_dp, real(my_nze,kind=DP), my_nze_dp
       call utils_abort('Error in sparse_count_ss: &
            &Integer overflow in sparse matrix index detected. Overflowed and &
            &expected values of nze and my_nze = '//error_string)
    end if

  end subroutine sparse_count_ss


  !==========================================================================!
  ! This routine generates a block index for a sparse matrix from a list of  !
  ! direct sphere-sphere overlaps and adds it to the sparse library.         !
  !--------------------------------------------------------------------------!
  ! Arguments:                                                               !
  !   nze     (input) : Number of nonzero elements                           !
  !   nzb     (input) : Number of nonzero blocks                             !
  !   my_nze  (inout) : Number of nonzero elements on this proc              !
  !   my_nzb  (input) : Number of nonzero blocks on this proc                !
  !   seg_nze (inout) : Number of nonzero elements per segment on this proc  !
  !   seg_nzb (input) : Number of nonzero blocks per segment on this proc    !
  !   name    (input) : Identifying name                                     !
  !--------------------------------------------------------------------------!
  ! Written by Peter Haynes, March 2004.                                     !
  ! Revised for distributed data, June 2006.                                 !
  ! Updated for SPAM3 by Nicholas Hine, May 2009.                            !
  !==========================================================================!

  subroutine sparse_index_ss(overlaps,col_blks,row_blks,name,nze,nzb,my_nze, &
       my_nzb,seg_nze,seg_nzb,rlib,transpose_name, row_par, col_par)

    use comms, only: pub_total_num_procs, pub_my_proc_id
    use parallel_strategy, only: OVERLAP_LIST
    use utils, only: utils_abort

    implicit none

    ! Arguments
    type(OVERLAP_LIST), intent(inout) :: overlaps
    integer, intent(in) :: col_blks       ! Number identifying col blk sizes
    integer, intent(in) :: row_blks       ! Number identifying row blk sizes
    character(len=*), intent(in) :: name  ! Identifying name
    integer(kind=LONG), intent(in) :: nze ! Number of nonzero elements
    integer, intent(in) :: nzb            ! Number of nonzero blocks
    integer, intent(inout) :: my_nze      ! Num of nonzero els on this proc
    integer, intent(in) :: my_nzb         ! Num of nonzero blocks on this proc
    integer, intent(inout) :: seg_nze(0:pub_total_num_procs-1)
    integer, intent(in) :: seg_nzb(0:pub_total_num_procs-1)
    integer, intent(out), optional :: rlib ! Library entry index return val
    character(len=*), intent(in), optional :: transpose_name ! Name of transpose
    type(PARAL_INFO), pointer, intent(in) :: row_par
    type(PARAL_INFO), pointer, intent(in) :: col_par

    ! Local variables
    integer :: my_nat     ! Number of atoms on this proc
    integer :: iblk,jblk  ! Atomic loop counters
    integer :: loc_iblk   ! Atomic loop counter in terms of local atoms
    integer :: ielems     ! Number of elements associated with atom iblk
    integer :: idx,ptr    ! Index and pointer
    integer :: dptr       ! Pointer to record (does not increment if dense)
    integer :: iovlap     ! Atomic overlap counter
    integer :: seg        ! Segment index
    integer :: seg_type   ! Segment type for this segment
    integer :: seg_start  ! Start position of this segment in the index

    ! Increment library counter
    num_library = num_library + 1
    if (num_library > max_library) &
         call sparse_library_full('sparse_index_ss')

    ! Set identification name
    library(num_library)%structure = name

    ! Set transpose name
    if (.not.present(transpose_name).and.(row_blks==col_blks)) then
       library(num_library)%transpose_structure = name
    else if (present(transpose_name)) then
       library(num_library)%transpose_structure = transpose_name
    else
       call utils_abort('Error in sparse_index_ss: &
            &Tranpose structure not specified')
    end if

    ! Set row and column block size list identifiers
    library(num_library)%col_blks = col_blks
    library(num_library)%row_blks = row_blks

    ! Set segment types and increase nze accordingly
    call sparse_segments_alloc(library(num_library),my_nze,seg_nze,seg_nzb,&
        row_par,col_par)

    ! Initialise dimensions of matrix
    call sparse_struc_alloc(library(num_library),nze,nzb,my_nze,my_nzb, &
        row_par, col_par)

    ! Ensure list of overlaps for each atom is in ascending order
    do iblk=col_par%my_first_blk,col_par%my_last_blk
       call internal_heapsort(overlaps%num_overlaps(iblk),overlaps%overlap_list(:,iblk))
    end do

    ! Initialise counters and pointers
    my_nat = col_par%my_last_blk - col_par%my_first_blk + 1
    idx = 1
    ptr = 1
    dptr = -1

    ! Loop over segments of the index
    do seg=0,pub_total_num_procs-1

       ! Set up entries in the seg_info(s_idx,:) and seg_info(s_ptr,:) arrays
       library(num_library)%seg_info(s_ptr,seg) = ptr
       library(num_library)%seg_info(s_idx,seg) = idx
       seg_type = library(num_library)%seg_info(s_type,seg)
       seg_start = idx

       ! Skip indexing of this segment if it contains no nonzero elements
       if (seg_type == SEG_BLANK) cycle

       if (seg/=pub_my_proc_id) then
          ! Set unused pointers to zero
          library(num_library)%blk_ptr(idx:idx+my_nat) = 0
       end if

       ! Go past list of column start positions for this segment
       idx = idx + my_nat + 1
       library(num_library)%blk_idx(seg_start) = idx

       ! Loop over atom block-columns iblk on this proc for this segment
       loc_iblk = 0
       do iblk=col_par%my_first_blk,col_par%my_last_blk
          loc_iblk = loc_iblk + 1

          ! Number of elements on atom iblk
          ielems = col_par%num_elems_on_atom(iblk,col_blks)

          ! Skip this atom if there are no elements on it - this is the case
          ! for atoms which are not Hubbard atoms when constructing V and W
          if (ielems==0) overlaps%num_overlaps(iblk) = 0

          ! Loop over overlapping block-rows (atoms)
          do iovlap=1,overlaps%num_overlaps(iblk)
             jblk = overlaps%overlap_list(iovlap,iblk)

             if (jblk < row_par%first_atom_on_proc(seg)) cycle
             if (jblk >= row_par%first_atom_on_proc(seg+1)) exit

             ! Skip this atom if there are no elements on it
             if (row_par%num_elems_on_atom(jblk,row_blks)==0) cycle

             ! Index contains block-row
             library(num_library)%blk_idx(idx) = jblk

             ! Pointer points to start of data
             if ( seg_type == SEG_SPARSE ) then
                dptr = ptr
                ! Increment block pointer by size of block
                ptr = ptr + ielems * row_par%num_elems_on_atom(jblk,row_blks)
             else if ( seg_type == SEG_DENSE ) then
                ! Find ptr to first element of this block in segment
                dptr = ptr + (col_par%first_elem_on_atom(iblk,col_blks) - &
                     col_par%first_elem_on_proc(pub_my_proc_id,col_blks)) * &
                     row_par%num_elems_on_proc(seg,row_blks) &
                     + row_par%first_elem_on_atom(jblk,row_blks) - &
                     row_par%first_elem_on_proc(seg,row_blks)
             end if

             library(num_library)%blk_ptr(idx) = dptr
             if (jblk == iblk) &
                  library(num_library)%blk_ptr(seg_start+loc_iblk-1) = dptr

             ! Increment block index by one
             idx = idx + 1

          end do  ! jblk

          ! Store idx in list of column start positions for this segment
          library(num_library)%blk_idx(seg_start+loc_iblk) = idx
          library(num_library)%blk_ptr(seg_start+loc_iblk) = 0

       end do  ! iblk

       ! For dense segments, move ptr on by total number of elements
       if ( seg_type == SEG_DENSE ) then
          ptr = ptr + col_par%num_elems_on_proc(pub_my_proc_id,col_blks) * &
               row_par%num_elems_on_proc(seg,row_blks)
       end if

       ! Set final pointer to end of data
       library(num_library)%blk_ptr(idx) = ptr

    end do  ! seg

    ! Record end of last segment
    library(num_library)%seg_info(s_idx,pub_total_num_procs) = idx
    library(num_library)%seg_info(s_ptr,pub_total_num_procs) = ptr

    ! Record library value of new matrix in rlib
    if (present(rlib)) rlib = num_library

  end subroutine sparse_index_ss


  !==========================================================================!
  ! This routine returns the numbers of matrix elements and blocks (total,   !
  ! per proc and per segment) in a sparse matrix whose sparsity pattern is   !
  ! the union of the nonzero elements of two already-defined sparse matrix   !
  ! indices.                                                                 !
  !--------------------------------------------------------------------------!
  ! Arguments:                                                               !
  !   alib    (input)  : Library entry for first predefined matrix           !
  !   blib    (input)  : Library entry for second predefined matrix          !
  !   nze     (output) : Total number of nonzero elements                    !
  !   nzb     (output) : Total number of nonzero blocks                      !
  !   my_nze  (output) : Number of nonzero elements on this proc             !
  !   my_nzb  (output) : Number of nonzero blocks on this proc               !
  !   seg_nze (output) : Number of nonzero elements per segment on this proc !
  !   seg_nzb (output) : Number of nonzero blocks per segment on this proc   !
  !--------------------------------------------------------------------------!
  ! Written by Nicholas Hine, June 2009.                                     !
  !==========================================================================!

  subroutine sparse_count_union(alib,blib,nze,nzb,my_nze,my_nzb, &
       seg_nze,seg_nzb)

    use comms, only: comms_allgather, comms_reduce, pub_total_num_procs
    use parallel_strategy, only: PARAL_INFO
    use utils, only: utils_abort, utils_assert, utils_alloc_check, &
         utils_dealloc_check

    implicit none

    ! Arguments
    integer, intent(in) :: alib
    integer, intent(in) :: blib
    integer(kind=LONG), intent(out) :: nze
    integer, intent(out) :: nzb
    integer, intent(out) :: my_nze
    integer, intent(out) :: my_nzb
    integer, intent(out) :: seg_nze(0:pub_total_num_procs-1)
    integer, intent(out) :: seg_nzb(0:pub_total_num_procs-1)

    ! Locals
    integer :: row_blks    ! Identifier for column blocking scheme
    integer :: col_blks    ! Identifier for row blocking scheme
    integer :: iblk,jblk   ! Atomic loop counters
    integer :: loc_iblk    ! Atomic loop counter in terms of local atoms
    integer :: iidx        ! Overlap index
    integer :: seg         ! Segment index
    integer :: aseg_start  ! Segment index start for library A
    integer :: bseg_start  ! Segment index start for library B
    integer :: aseg_type   ! Segment type for this segment of A
    integer :: bseg_type   ! Segment type for this segment of B
    integer :: num_ovlaps  ! Number of overlaps for this atom in this segment
    integer :: ielems      ! Number of elements associated with atom iblk
    integer, allocatable :: ibuf(:) ! Communications buffer for comms_allgather
    integer :: ierr        ! Error Flag
    logical,allocatable :: tfound_blk(:)
    integer,allocatable :: tfound_idx(:)
    type(PARAL_INFO), pointer :: row_par ! Parallel strategy for rows
    type(PARAL_INFO), pointer :: col_par ! Parallel strategy for columns
    real(kind=DP) :: nze_dp     ! Real nze for checking for integer overflow
    real(kind=DP) :: my_nze_dp  ! Real my_nze for checking for integer overflow
    character(len=:), allocatable :: error_string ! For output of error message
    logical :: integer_overflow ! Has an integer overflow been detected?

    ! No integer overflow... yet
    integer_overflow = .false.

    ! Zero counters
    my_nze = 0
    my_nzb = 0
    my_nze_dp = real(my_nze,kind=DP)

    ! Check inputs
    call utils_assert((alib>=1).and.(alib<=num_library), &
       'Error in sparse_count_union: invalid library entry provided for alib.')
    call utils_assert((blib>=1).and.(blib<=num_library), &
       'Error in sparse_count_union: invalid library entry provided for blib.')

    ! Find row and column blocking schemes
    row_blks = library(alib)%row_blks
    col_blks = library(alib)%col_blks

    ! rc2013: assign parallel strategy
    ! rc2013: count_union requires A and B have same parallel strategies
    row_par => library(blib)%row_par
    col_par => library(blib)%col_par

    ! rc2013: check parallel strategies
    call sparse_par_check(alib, blib, routine='sparse_count_union')

    ! Sanity check
    call utils_assert(library(blib)%row_blks == row_blks, 'Error in &
         &sparse_count_union: A and B row blocking schemes do not match')
    call utils_assert(library(blib)%col_blks == col_blks, 'Error in &
         &sparse_count_union: A and B col blocking schemes do not match')

    ! rc2013: this may need to be formatted for different system sizes
    allocate(tfound_idx(library(alib)%nblk),stat=ierr)
    call utils_alloc_check('sparse_count_union','tfound_idx',ierr)
    allocate(tfound_blk(library(alib)%nblk),stat=ierr)
    call utils_alloc_check('sparse_count_union','tfound_blk',ierr)
    tfound_blk = .false.

    ! Loop over segments of both structures
    do seg=0,pub_total_num_procs-1

       seg_nze(seg) = 0
       seg_nzb(seg) = 0
       aseg_start = library(alib)%seg_info(s_idx,seg)
       bseg_start = library(blib)%seg_info(s_idx,seg)
       aseg_type = library(alib)%seg_info(s_type,seg)
       bseg_type = library(blib)%seg_info(s_type,seg)

       ! Loop over atoms on this proc
       loc_iblk = 0
       do iblk=col_par%my_first_blk,col_par%my_last_blk
          loc_iblk = loc_iblk + 1

          ielems = col_par%num_elems_on_atom(iblk,col_blks)

          ! Set counter of number of overlaps found to zero
          num_ovlaps = 0
          if (aseg_type /= SEG_BLANK) then
             ! Loop over nonzero blocks in this column in structure A
             do iidx=library(alib)%blk_idx(aseg_start+loc_iblk-1), &
                  library(alib)%blk_idx(aseg_start+loc_iblk)-1
                jblk = library(alib)%blk_idx(iidx)
                tfound_blk(jblk) = .true.
                num_ovlaps = num_ovlaps + 1
                tfound_idx(num_ovlaps) = jblk
             end do  ! iidx
          end if

          if (bseg_type /= SEG_BLANK) then
             ! Loop over nonzero blocks in this column in structure B
             do iidx=library(blib)%blk_idx(bseg_start+loc_iblk-1), &
                  library(blib)%blk_idx(bseg_start+loc_iblk)-1
                jblk = library(blib)%blk_idx(iidx)

                ! Check if this block has already been found
                if (.not.tfound_blk(jblk)) then

                   tfound_blk(jblk) = .true.
                   num_ovlaps = num_ovlaps + 1
                   tfound_idx(num_ovlaps) = jblk

                end if

             end do  ! iidx
          end if

          ! Loop over all overlaps kat (block-rows)
          do iidx=1,num_ovlaps
             jblk = tfound_idx(iidx)

             ! Whole atomic block counts
             my_nze = my_nze + ielems * row_par%num_elems_on_atom(jblk,row_blks)

             ! Increment real version of my_nze, for integer overflow checking
             my_nze_dp = my_nze_dp + real(ielems,kind=DP) * &
                  real(row_par%num_elems_on_atom(jblk,row_blks),kind=DP)

             ! Check for integer overflow
             if (my_nze_dp > HUGE(my_nze)) integer_overflow = .true.

             my_nzb = my_nzb + 1
             seg_nze(seg) = seg_nze(seg) + &
                  ielems * row_par%num_elems_on_atom(jblk,row_blks)
             seg_nzb(seg) = seg_nzb(seg) + 1

             ! Reset flags
             tfound_blk(jblk) = .false.

          end do  ! iidx

       end do  ! iblk

    end do  ! seg

    ! Sum up over all procs
    ! JCW: Susceptible to integer overflow, since we have no comms_reduce
    !      procedures for integer(kind=LONG)
    !ibuf(1) = my_nze
    !ibuf(2) = my_nzb
    !call comms_reduce('SUM',ibuf)
    !nze = ibuf(1)
    !nzb = ibuf(2)
    ! JCW: In the absence of a comms_reduce for LONG integers, we avoid integer
    !      overflow by allgathering to an array of integers on each proc, then
    !      locally summing into nze, which is integer(kind=LONG)
    allocate(ibuf(pub_total_num_procs),stat=ierr)
    call utils_alloc_check('sparse_count_union','ibuf',ierr)
    ! nze
    ibuf(:) = 0
    call comms_allgather(ibuf,my_nze)
    nze = sum(int(ibuf(:),kind=LONG))
    ! nzb
    ibuf(:) = 0
    call comms_allgather(ibuf,my_nzb)
    nzb = sum(ibuf(:))
    deallocate(ibuf,stat=ierr)
    call utils_dealloc_check('sparse_count_union','ibuf',ierr)

    ! Sum over all procs for integer overflow checking
    nze_dp = my_nze_dp
    call comms_reduce('SUM',nze_dp)
    if (nze_dp > HUGE(nze)) integer_overflow = .true.

    if (integer_overflow) then
       allocate(character(len=80) :: error_string)
       ! len=80, because we want 4*20 width numeric values
       write(error_string,'(4es20.10)') &
            real(nze,kind=DP), nze_dp, real(my_nze,kind=DP), my_nze_dp
       call utils_abort('Error in sparse_count_union: &
            &Integer overflow in sparse matrix index detected. Overflowed and &
            &expected values of nze and my_nze = '//error_string)
    end if

    deallocate(tfound_blk,stat=ierr)
    call utils_dealloc_check('sparse_count_union','tfound_blk',ierr)
    deallocate(tfound_idx,stat=ierr)
    call utils_dealloc_check('sparse_count_union','tfound_idx',ierr)

  end subroutine sparse_count_union


  !==========================================================================!
  ! This routine generates a block index for a sparse matrix whose sparsity  !
  ! pattern is the union of the nonzero elements of two already-defined      !
  ! sparse matrix indices.                                                   !
  !--------------------------------------------------------------------------!
  ! Arguments:                                                               !
  !   alib    (input) : Library entry for first predefined matrix            !
  !   blib    (input) : Library entry for second predefined matrix           !
  !   nze     (input) : Total number of nonzero elements                     !
  !   nzb     (input) : Total number of nonzero blocks                       !
  !   my_nze  (inout) : Number of nonzero elements on this proc              !
  !   my_nzb  (input) : Number of nonzero blocks on this proc                !
  !   seg_nze (inout) : Number of nonzero elements per segment on this proc  !
  !   seg_nzb (input) : Number of nonzero blocks per segment on this proc    !
  !   name    (input) : Identifying name                                     !
  !--------------------------------------------------------------------------!
  ! Written by Nicholas Hine, June 2009.                                     !
  !==========================================================================!

  subroutine sparse_index_union(alib_in,blib_in,name,nze,nzb,my_nze,my_nzb, &
       seg_nze,seg_nzb,rlib,transpose_name)

    use comms, only: pub_my_proc_id, pub_total_num_procs
    use parallel_strategy, only: PARAL_INFO
    use utils, only: utils_assert, utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    integer, intent(in) :: alib_in
    integer, intent(in) :: blib_in
    character(len=*), intent(in) :: name
    integer(kind=LONG), intent(in) :: nze
    integer, intent(in) :: nzb
    integer, intent(inout) :: my_nze
    integer, intent(in) :: my_nzb
    integer, intent(inout) :: seg_nze(0:pub_total_num_procs-1)
    integer, intent(in) :: seg_nzb(0:pub_total_num_procs-1)
    integer, intent(out), optional :: rlib ! Library entry index return val
    character(len=*), intent(in), optional :: transpose_name  ! Name of transpose

    ! Locals
    integer :: row_blks      ! Identifier for column blocking scheme
    integer :: col_blks      ! Identifier for row blocking scheme
    integer :: iblk,jblk     ! Atomic loop counters
    integer :: loc_iblk      ! Atomic loop counter in terms of local atoms
    integer :: num_ovlaps    ! Number of overlaps for this atom in this segment
    integer :: iidx          ! Overlap index
    integer :: idx,ptr,dptr  ! Index, pointer and dense pointer
    integer :: my_nat        ! Number of atoms on this nod
    integer :: ielems        ! Number of elements associated with atom iblk
    integer :: seg           ! Segment index
    integer :: seg_type      ! Type for this segment of C
    integer :: seg_start     ! Index start position for this segment of C
    integer :: aseg_start    ! Segment index start for library A
    integer :: bseg_start    ! Segment index start for library B
    integer :: aseg_type     ! Segment type for this segment of A
    integer :: bseg_type     ! Segment type for this segment of B
    integer :: ulib          ! Library index of union
    integer :: temp_lib      ! Temporary library index
    integer :: alib          ! Local copy of alib variable
    integer :: blib          ! Local copy of alib variable
    integer :: ierr        ! Error Flag
    logical,allocatable :: tfound_blk(:)
    integer,allocatable :: tfound_idx(:)
    type(PARAL_INFO), pointer :: row_par ! parallel strategy for rows
    type(PARAL_INFO), pointer :: col_par ! parallel strategy for columns

    ! Check inputs
    call utils_assert((alib_in>=1).and.(alib_in<=num_library), &
       'Error in sparse_index_union: invalid library entry provided for alib.')
    call utils_assert((blib_in>=1).and.(blib_in<=num_library), &
       'Error in sparse_index_union: invalid library entry provided for blib.')

    alib = alib_in
    blib = blib_in

    ! Find row and column blocking schemes
    row_blks = library(alib)%row_blks
    col_blks = library(alib)%col_blks

    ! rc2013: find row and column parallel strategies
    ! rc2013: index_union should only be used where A and B have same parallel strategies
    row_par => library(alib)%row_par
    col_par => library(alib)%col_par

    ! Sanity check
    call utils_assert(library(blib)%row_blks == row_blks, 'Error in &
         &sparse_index_union: A and B row blocking schemes do not match')
    call utils_assert(library(blib)%col_blks == col_blks, 'Error in &
         &sparse_index_union: A and B row blocking schemes do not match')

    ! rc2013: check for consistency of parallel strategies
    call sparse_par_check(alib, blib, routine='sparse_index_union')

    allocate(tfound_idx(library(alib)%nblk),stat=ierr)
    call utils_alloc_check('sparse_index_union','tfound_idx',ierr)
    allocate(tfound_blk(library(alib)%nblk),stat=ierr)
    call utils_alloc_check('sparse_index_union','tfound_blk',ierr)
    tfound_blk = .false.

    ! Check if the new structure name already exists
    call sparse_search_library(name,ulib)
    temp_lib = 0

    ! If the structure code was not found, create a new one
    if (ulib == 0) then

       ! Increment library counter
       num_library = num_library + 1
       if (num_library > max_library) &
            call sparse_library_full('sparse_index_union')

       ! Set identification name
       library(num_library)%structure = name

       ! Set row and column block size list identifiers
       library(num_library)%col_blks = col_blks
       library(num_library)%row_blks = row_blks

       ulib = num_library

    else ! this structure already exists

       ! Check if ulib is also either alib or blib
       if (ulib==alib) then
          temp_lib = num_library + 1
          call sparse_struc_copy(library(temp_lib),library(alib))
          alib = temp_lib
          if (ulib==blib) blib = temp_lib
       else if (ulib==blib) then
          temp_lib = num_library + 1
          call sparse_struc_copy(library(temp_lib),library(blib))
          blib = temp_lib
       end if

       ! Deallocate previous storage associated with ulib
       call sparse_segments_dealloc(library(ulib))
       call sparse_struc_dealloc(library(ulib))

    end if

    ! Set transpose name
    if (.not.present(transpose_name).and.(row_blks==col_blks)) then
       library(ulib)%transpose_structure = name
    else if (present(transpose_name)) then
       library(ulib)%transpose_structure = transpose_name
    else
       !call utils_abort('Error in sparse_index_union: &
       !     &Tranpose structure not specified')
    end if

    ! Set segment types and increase nze accordingly
    call sparse_segments_alloc(library(ulib),my_nze,seg_nze,seg_nzb,&
        row_par,col_par)

    ! Initialise dimensions of matrix
    call sparse_struc_alloc(library(ulib),nze,nzb,my_nze,my_nzb, &
         row_par, col_par)

    ! Initialise counters and pointers
    my_nat = col_par%my_last_blk - col_par%my_first_blk + 1
    idx = 1
    ptr = 1
    dptr = -1

    ! Loop over segments of the index
    do seg=0,pub_total_num_procs-1

       ! Set up entries in the seg_info(s_idx,:) and seg_info(s_ptr,:) arrays
       library(ulib)%seg_info(s_ptr,seg) = ptr
       library(ulib)%seg_info(s_idx,seg) = idx
       seg_type = library(ulib)%seg_info(s_type,seg)
       seg_start = idx

       aseg_start = library(alib)%seg_info(s_idx,seg)
       bseg_start = library(blib)%seg_info(s_idx,seg)
       aseg_type = library(alib)%seg_info(s_type,seg)
       bseg_type = library(blib)%seg_info(s_type,seg)

       ! Skip indexing of this segment if it contains no nonzero elements
       if (seg_type == SEG_BLANK) cycle

       if (seg/=pub_my_proc_id) then
          ! Set unused pointers to zero
          library(ulib)%blk_ptr(idx:idx+my_nat) = 0
       end if

       ! Go past list of column start positions for this segment
       idx = idx + my_nat + 1
       library(ulib)%blk_idx(seg_start) = idx

       ! Loop over atom block-columns iblk on this proc for this segment
       loc_iblk = 0
       do iblk=col_par%my_first_blk,col_par%my_last_blk
          loc_iblk = loc_iblk + 1

          ! Number of elements on atom iblk
          ielems = col_par%num_elems_on_atom(iblk,col_blks)

          ! Set counter of number of overlaps found to zero
          num_ovlaps = 0

          if (aseg_type/=SEG_BLANK) then
             ! Loop over nonzero blocks in this column in structure A
             do iidx=library(alib)%blk_idx(aseg_start+loc_iblk-1), &
                  library(alib)%blk_idx(aseg_start+loc_iblk)-1
                jblk = library(alib)%blk_idx(iidx)
                tfound_blk(jblk) = .true.
                num_ovlaps = num_ovlaps + 1
                tfound_idx(num_ovlaps) = jblk
             end do  ! iidx
          end if

          if (bseg_type/=SEG_BLANK) then
             ! Loop over nonzero blocks in this column in structure B
             do iidx=library(blib)%blk_idx(bseg_start+loc_iblk-1), &
                  library(blib)%blk_idx(bseg_start+loc_iblk)-1
                jblk = library(blib)%blk_idx(iidx)

                ! Check if this block has already been found
                if (.not.tfound_blk(jblk)) then

                   tfound_blk(jblk) = .true.
                   num_ovlaps = num_ovlaps + 1
                   tfound_idx(num_ovlaps) = jblk

                end if

             end do  ! iidx
          end if

          ! Sort list of overlaps to ensure it is in ascending order
          call internal_heapsort(num_ovlaps,tfound_idx)

          do iidx=1,num_ovlaps
             jblk = tfound_idx(iidx)

             ! Index contains block-row
             library(ulib)%blk_idx(idx) = jblk

             ! Pointer points to start of data
             if ( seg_type == SEG_SPARSE ) then
                dptr = ptr
                ! Increment block pointer by size of block
                ptr = ptr + ielems * row_par%num_elems_on_atom(jblk,row_blks)
             else if ( seg_type == SEG_DENSE ) then
                ! Find ptr to first element of this block in segment
                dptr = ptr + (col_par%first_elem_on_atom(iblk,col_blks) - &
                     col_par%first_elem_on_proc(pub_my_proc_id,col_blks)) * &
                     row_par%num_elems_on_proc(seg,row_blks) &
                     + row_par%first_elem_on_atom(jblk,row_blks) - &
                     row_par%first_elem_on_proc(seg,row_blks)
             end if

             library(ulib)%blk_ptr(idx) = dptr
             if (jblk == iblk) &
                  library(ulib)%blk_ptr(seg_start+loc_iblk-1) = dptr

             ! Increment block index by one
             idx = idx + 1

            ! Reset flags
            tfound_blk(jblk) = .false.

          end do  ! jblk

          ! Store idx in list of column start positions for this segment
          library(ulib)%blk_idx(seg_start+loc_iblk) = idx
          library(ulib)%blk_ptr(seg_start+loc_iblk) = 0

       end do  ! iblk

       ! For dense segments, move ptr on by total number of elements
       if ( seg_type == SEG_DENSE ) then
          ptr = ptr + col_par%num_elems_on_proc(pub_my_proc_id,col_blks) * &
               row_par%num_elems_on_proc(seg,row_blks)
       end if

       ! Set final pointer to end of data
       library(ulib)%blk_ptr(idx) = ptr

    end do  ! seg

    ! Record end of last segment
    library(ulib)%seg_info(s_idx,pub_total_num_procs) = idx
    library(ulib)%seg_info(s_ptr,pub_total_num_procs) = ptr

    if (temp_lib > 0) then
       call sparse_struc_dealloc(library(temp_lib))
    end if

    ! Record library value of new matrix in rlib
    if (present(rlib)) rlib = ulib

    deallocate(tfound_blk,stat=ierr)
    call utils_dealloc_check('sparse_index_union','tfound_blk',ierr)
    deallocate(tfound_idx,stat=ierr)
    call utils_dealloc_check('sparse_index_union','tfound_idx',ierr)

  end subroutine sparse_index_union


  subroutine sparse_show_segment_filling

    use comms, only: pub_on_root

    implicit none

    ! Local variables
    type(SPAM3) :: q,h,s,k,ks,ksk,ksks

    q%structure = 'Q'
    call sparse_create(q)
    h%structure = 'H'
    call sparse_create(h)
    s%structure = 'S'
    call sparse_create(s)
    k%structure = 'K'
    call sparse_create(k)
    ks%structure = 'KS'
    call sparse_create(ks)
    ksk%structure = 'KSK'
    call sparse_create(ksk)
    ksks%structure = 'KSKS'
    call sparse_create(ksks)

    !if (pub_on_root) write(stdout,'(a)') 'Q'
    !call internal_show_seg_frac(q)
    !if (pub_on_root) write(stdout,'(a)') 'H'
    !call internal_show_seg_frac(h)
    if (pub_on_root) write(stdout,'(a)') 'S'
    call internal_show_seg_frac(s)
    if (pub_on_root) write(stdout,'(a)') 'K'
    call internal_show_seg_frac(k)
    if (pub_on_root) write(stdout,'(a)') 'KS'
    call internal_show_seg_frac(ks)
    if (pub_on_root) write(stdout,'(a)') 'KSK'
    call internal_show_seg_frac(ksk)
    !if (pub_on_root) write(stdout,'(a)') 'KSKS'
    !call internal_show_seg_frac(ksks)

    call sparse_destroy(ksks)
    call sparse_destroy(ksk)
    call sparse_destroy(ks)
    call sparse_destroy(k)
    call sparse_destroy(s)
    call sparse_destroy(h)
    call sparse_destroy(q)

  contains

    subroutine internal_show_segs(mat)

      use comms, only: pub_on_root,pub_my_proc_id,pub_total_num_procs, &
           comms_barrier,comms_bcast
      use utils, only: utils_flush, utils_alloc_check, utils_dealloc_check

      implicit none

      ! Arguments
      type(SPAM3),intent(in) :: mat

      ! Locals
      integer :: proc,seg
      integer :: ierr
      character(4), allocatable :: segchar(:,:),segline(:)
      character(20) :: fmt, tmp

      allocate(segchar(0:pub_total_num_procs-1,0:pub_total_num_procs-1), &
           stat=ierr)
      call utils_alloc_check('internal_show_segs','segchar',ierr)
      allocate(segline(0:pub_total_num_procs-1),stat=ierr)
      call utils_alloc_check('internal_show_segs','segline',ierr)

      do seg=0,pub_total_num_procs-1
         if (library(mat%lib)%seg_info(s_type,seg)==SEG_SPARSE) then
            segchar(seg,pub_my_proc_id) = '+'
         else if (library(mat%lib)%seg_info(s_type,seg)==SEG_DENSE) then
            segchar(seg,pub_my_proc_id) = '#'
         else
            segchar(seg,pub_my_proc_id) = '.'
         end if
      end do

      do proc=0,pub_total_num_procs-1
         call comms_bcast(proc,segchar(:,proc),pub_total_num_procs)
      end do

      if (pub_on_root) then
         write(stdout,'(a)',advance='no') '     '
         do proc=0,pub_total_num_procs-1
            write(stdout,'(i2)',advance='no') proc
         end do
         write(stdout,*)
         write(tmp,'(i10)') pub_total_num_procs
         write(fmt,'(a,a,a)')'(i5,',trim(adjustl(tmp)),'a2)'
         do seg=0,pub_total_num_procs-1
            segline(:) = segchar(seg,:)
            write(stdout,fmt) seg,segline(:)
         end do
      end if

      call comms_barrier
      call utils_flush
      call comms_barrier

      deallocate(segline,stat=ierr)
      call utils_dealloc_check('internal_show_segs','segline',ierr)
      deallocate(segchar,stat=ierr)
      call utils_dealloc_check('internal_show_segs','segchar',ierr)

    end subroutine internal_show_segs

    subroutine internal_show_seg_frac(mat)

      use comms, only: pub_on_root,pub_my_proc_id,pub_total_num_procs, &
           comms_barrier,comms_bcast
      use utils, only: utils_flush, utils_alloc_check, utils_dealloc_check

      implicit none

      ! Arguments
      type(SPAM3),intent(in) :: mat

      ! Locals
      integer :: proc,seg
      integer :: ierr
      real(kind=DP), allocatable :: segfrac(:,:),segline(:)
      character(20) :: fmt, tmp
      type(PARAL_INFO), pointer :: row_par
      type(PARAL_INFO), pointer :: col_par

      ! rc2013: assign parallel strategies
      row_par => library(mat%lib)%row_par
      col_par => library(mat%lib)%col_par

      allocate(segfrac(0:pub_total_num_procs-1,0:pub_total_num_procs-1), &
           stat=ierr)
      call utils_alloc_check('internal_show_seg_frac','segfrac',ierr)
      allocate(segline(0:pub_total_num_procs-1),stat=ierr)
      call utils_alloc_check('internal_show_seg_frac','segline',ierr)

      do seg=0,pub_total_num_procs-1
         segfrac(seg,pub_my_proc_id) = (library(mat%lib)%seg_info(s_ptr,seg+1)-&
              library(mat%lib)%seg_info(s_ptr,seg))/ &
              real(col_par%num_elems_on_proc(pub_my_proc_id, &
              library(mat%lib)%col_blks)*&
              row_par%num_elems_on_proc(seg,library(mat%lib)%row_blks),kind=DP)
      end do

      do proc=0,pub_total_num_procs-1
         call comms_bcast(proc,segfrac(:,proc),pub_total_num_procs)
      end do

      if (pub_on_root) then
         write(stdout,'(a)',advance='no') '     '
         do proc=0,pub_total_num_procs-1
            write(stdout,'(i5)',advance='no') proc
         end do
         write(stdout,*)
         write(tmp,'(i10)') pub_total_num_procs
         write(fmt,'(a,a,a)')'(i5,',trim(adjustl(tmp)),'f5.2)'
         do seg=0,pub_total_num_procs-1
            segline(:) = segfrac(seg,:)
            write(stdout,fmt) seg,segline(:)
         end do
      end if

      call comms_barrier
      call utils_flush
      call comms_barrier

      deallocate(segline,stat=ierr)
      call utils_dealloc_check('internal_show_seg_frac','segline',ierr)
      deallocate(segfrac,stat=ierr)
      call utils_dealloc_check('internal_show_seg_frac','segfrac',ierr)

    end subroutine internal_show_seg_frac

  end subroutine sparse_show_segment_filling

  !==========================================================================!
  ! This subroutine sorts a list of integers (to ensure that blocks are      !
  ! always indexed in ascending order.                                       !
  !--------------------------------------------------------------------------!
  ! Arguments:                                                               !
  !   n     (input) : The length of the list                                 !
  !   list  (inout) : The list to be sorted                                  !
  !--------------------------------------------------------------------------!
  ! Written by Peter Haynes, June 2006.                                      !
  !==========================================================================!

  subroutine internal_heapsort(n,list)

    implicit none

    ! Arguments
    integer, intent(in)    :: n
    integer, intent(inout) :: list(:)

    ! Local variables
    integer :: i,j                 ! Loop variables over list
    integer :: ic,is               ! Counters for heap creation/selection
    integer :: temp                ! Temporary copy of list item

    ! If there is only one item in the list, there's no sorting to be done!
    if (n <= 1) return

    ! Do the heap sort
    ic = n/2+1
    is = n
    do
       if (ic > 1) then ! in heap creation phase
          ic = ic - 1
          temp = list(ic)
       else ! in heap selection phase
          temp = list(is)
          list(is) = list(1)
          is = is - 1
          if (is == 1) then
             list(1) = temp
             exit
          end if
       end if
       ! Sift down temporary copy to correct level in heap
       i = ic
       j = 2 * ic
       do while (j <= is)
          if (j < is) then
             if (list(j) < list(j+1)) j = j + 1
          end if
          if (temp < list(j)) then ! demote temporary copy
             list(i) = list(j)
             i = j
             j = 2 * j
          else ! found level for temporary copy
             j = is + 1
          end if
       end do
       list(i) = temp
    end do

  end subroutine internal_heapsort


  !============================================================================!
  ! This subroutine creates a new sparse matrix by:                            !
  !   (i) looking its name up in the library of known structures, or           !
  !  (ii) copying the structure from another matrix, or                        !
  ! (iii) generating a new structure from the product of two sparse matrices.  !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   new_spam (output) : The sparse matrix to be created                      !
  !   spama    (input)  : An optional sparse matrix whose structure is used    !
  !   spamb    (input)  : A second (optional) sparse matrix whose structure    !
  !                       will be used                                         !
  !   iscmplx  (input)  : A flag to indicate the data stored in this matrix is !
  !                       complex                                              !
  ! embed_struc (input) : Structure code to be used for initialising new SPAM3 !
  !                       if this is an embedding matrix.                      !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, May 2009, based on SPAM2 sparse_create routine,  !
  ! Written by Peter Haynes, March 2004,                                       !
  ! Revised for distributed data, June 2006,                                   !
  ! Modified for dense matrices by Nicholas Hine, Dec 2007.                    !
  ! Minor edit for embedding structure codes by Robert Charlton, 16/08/2017.   !
  !============================================================================!

  subroutine sparse_create(new_spam,spama,spamb,iscmplx,rlib,embed_struc)

    use comms, only: comms_reduce, comms_send, comms_allgather, &
         comms_barrier, comms_free, comms_bcast, comms_probe, comms_wait, &
         pub_my_proc_id, pub_total_num_procs, &
         pub_first_proc_in_group, pub_last_proc_in_group, pub_group_comm, &
         pub_num_comms_groups, pub_comms_group_size
    use constants, only: stdout
!$  use comms, only: pub_my_rank_in_group
    use parallel_strategy, only: PARAL_INFO
!$  use rundat, only: pub_threads_max
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_abort, &
         utils_assert

    implicit none

    ! Arguments
    type(SPAM3), intent(inout) :: new_spam      ! Sparse matrix to be created
    type(SPAM3), optional, intent(in) :: spama  ! First (optional) matrix to use
    type(SPAM3), optional, intent(in) :: spamb  ! Other (optional) matrix to use
    logical, optional, intent(in) :: iscmplx    ! Flag for complex data (opt)
    integer, optional, intent(out) :: rlib      ! Optional library index
    character(len=*), optional, intent(in) :: embed_struc ! Optional embedding
    ! structure code (set in sparse_embed_create)

    ! Local variables
    logical :: loc_iscmplx ! Local copy of optional iscmplx argument
    integer :: ierr        ! Error flag
    integer :: ilib        ! Library entry
    integer(kind=LONG) :: nze ! Number of nonzero elements in new matrix
    integer :: nzb         ! Number of nonzero blocks in new matrix
    integer :: my_nze      ! Number of nonzero els in new matrix on this proc
    integer :: my_nzb      ! Number of nonzero blks in new matrix on this proc
    integer :: alib,blib   ! Library entries for spama and spamb
    integer,allocatable,save :: num_ovlaps(:) ! Num overlaps found for each atom
    integer,allocatable :: seg_nzb(:) ! Number of nonzero blocks on each seg
    integer,allocatable :: seg_nze(:) ! Number of nonzero elements on each seg
    logical,allocatable :: tfound_blk(:)
    integer,allocatable :: tfound_idx(:)
    type(COM3) :: acom
    type(SHARE3) :: bshare
    integer :: num_ovlaps_size
    type(PARAL_INFO), pointer :: row_par ! parallel strategy for rows
    type(PARAL_INFO), pointer :: col_par ! parallel strategy for columns
    integer :: first_atom_in_group       ! First atom in this proc's comms group
    integer :: last_atom_in_group        ! Last atom in this proc's comms group
    type(PARAL_INFO), pointer :: acol_par ! parallel strategy for A columns
    ! Switch for whether this group contains procs that are assigned
    ! to the region associated with the columns of this SPAM3
    logical :: relevant_group

    if (present(spamb) .and. (.not. present(spama))) then
       call utils_abort('Error in sparse_create: &
            &wrong combination of optional arguments (spamb requires spama).')
    end if

    if ((.not. present(spama)) .and. (.not. present(spamb))) then

       ! Check name of new matrix against library of known structures

       ! Search library for matching entry
       call sparse_search_library(new_spam%structure,ilib)

       ! Quit if not found
       call utils_assert(ilib /= 0, 'Error in sparse_create: &
               &no library structure found to match '//trim(new_spam%structure))

       ! Copy structure from this library entry
       new_spam%lib = ilib

       ! Set default for optional argument
       loc_iscmplx = .false.
       if (present(iscmplx)) loc_iscmplx = iscmplx

       ! Allocate space for data
       call sparse_data_alloc(new_spam,loc_iscmplx)

    else if (.not. present(spamb)) then

       ! Copy structure from given sparse matrix: spama
       new_spam%lib = spama%lib

       ! Set default for optional argument
       loc_iscmplx = spama%iscmplx
       if (present(iscmplx)) loc_iscmplx = iscmplx
       call sparse_par_check(new_spam%lib, spama%lib, routine='sparse_create')

       ! Allocate space for data
       call sparse_data_alloc(new_spam,loc_iscmplx)

    else
       ! Generate a new structure from the product of spama and spamb

       ! Set default for optional argument
       loc_iscmplx = (spama%iscmplx .or. spamb%iscmplx)
       if (present(iscmplx)) loc_iscmplx = iscmplx

       ! Generate the new structure identifier
       ! rc2013: use embedding structure instead to avoid spurious structures
       if(present(embed_struc)) then
          write(new_spam%structure,'(a28)') trim(embed_struc)
       else
          write(new_spam%structure,'(a28)') trim(spama%structure)// &
               trim(spamb%structure)
       end if
       new_spam%structure = adjustl(new_spam%structure)

       ! Check if this already exists in the library
       call sparse_search_library(new_spam%structure,ilib)

       ! If it does...
       if (ilib > 0) then

          ! Copy the known structure from the library
          new_spam%lib = ilib
          ! rc2013: check that the parallel strategies are consistent
          call sparse_par_check(new_spam%lib, spama%lib, spamb%lib, 'sparse_create')

          ! Allocate space for data
          call sparse_data_alloc(new_spam,loc_iscmplx)

       else    ! ... otherwise some work must be done!

          if (.not.pub_sparse_allow_new_matrices) then
             !if (pub_on_root) then
             !   write(stdout,'(a)') 'WARNING in sparse_create: &
             !        &Creation of new matrix structures is discouraged'
             !   write(stdout,'(a)') 'outside of the routine &
             !        &sparse_init_rep'
             !   write(stdout,'(6a)') 'Matrix Structures: ', &
             !        trim(spama%structure),',',trim(spamb%structure), ',', &
             !        trim(new_spam%structure)
             !end if
          end if

          ! Find library entries for spama and spamb
          alib = spama%lib ; blib = spamb%lib

          call utils_assert (library(alib)%col_blks == library(blib)%row_blks, &
               'Error in sparse_create: column block scheme of matrix structure&
               & '//trim(library(alib)%structure)//' does not match row block &
               &scheme of matrix structure '//trim(library(blib)%structure)//&
               'Matrix product cannot be indexed.')

          ! Increment library counter
          num_library = num_library + 1
          if (num_library > max_library) &
               call sparse_library_full('sparse_create')

          ! Enter structure into library
          library(num_library)%structure = new_spam%structure
          library(num_library)%col_blks = library(blib)%col_blks
          library(num_library)%row_blks = library(alib)%row_blks

          ! rc2013: extract the par strategies from spama and spamb
          row_par => library(alib)%row_par
          col_par => library(blib)%col_par
          ! rc2013: we also need the A column par for product
          acol_par => library(alib)%col_par

          ! jcap: check if this group contains any procs that are
          ! assigned to the appropriate region
          if ( ( (pub_first_proc_in_group.ge.col_par%first_proc).and.&
               (pub_first_proc_in_group.le.col_par%first_proc+col_par%num_procs-1) ) .or. &
               ( (pub_first_proc_in_group.lt.col_par%first_proc).and.&
               (pub_last_proc_in_group.ge.col_par%first_proc) ) ) then
             relevant_group=.true.
             ! rc2013: get indices of the first and last atoms in this proc's group
             first_atom_in_group = sparse_first_atom_in_group(col_par, &
                  pub_first_proc_in_group, pub_last_proc_in_group)
             last_atom_in_group = sparse_last_atom_in_group(col_par, &
                  pub_first_proc_in_group, pub_last_proc_in_group)
          else
             relevant_group=.false.
             first_atom_in_group = 0
             last_atom_in_group = 0
          end if

          ! Allocate workspace
          num_ovlaps_size = last_atom_in_group - first_atom_in_group + 1
          allocate(num_ovlaps(num_ovlaps_size),stat=ierr)
          call utils_alloc_check('sparse_create','num_ovlaps',ierr)
          allocate(seg_nzb(0:pub_total_num_procs-1),stat=ierr)
          call utils_alloc_check('sparse_create','seg_nzb',ierr)
          allocate(seg_nze(0:pub_total_num_procs-1),stat=ierr)
          call utils_alloc_check('sparse_create','seg_nze',ierr)

          ! Allocate arrays and initialise comms
          call sparse_com_allocate(acom,spama,2, &
               alloc_mtx=.false.,cropped=.false.,seg=.false.,groupseg=.false.)
          call sparse_init_share(bshare,spamb,share_data=.false., &
               idx_only=.true.)

          ! Count the number of nonzero blocks and elements using the library
          call internal_count_product(nze,nzb,my_nze,my_nzb)

          ! Allocate arrays for the structure segment arrays
          call sparse_segments_alloc(library(num_library),my_nze,seg_nze, &
               seg_nzb,row_par,col_par)

          ! Allocate arrays for matrix structure
          call sparse_struc_alloc(library(num_library),nze,nzb,my_nze,my_nzb, &
              row_par, col_par)

          ! Generate the new block sparse index
          call internal_index_product

          ! Deallocate arrays and exit comms
          call sparse_exit_share(bshare,.false.,.true.)
          call sparse_com_deallocate(acom,dealloc_mtx=.false.,cropped=.false.)

          ! Deallocate workspace
          deallocate(seg_nze,stat=ierr)
          call utils_dealloc_check('sparse_create','seg_nze',ierr)
          deallocate(seg_nzb,stat=ierr)
          call utils_dealloc_check('sparse_create','seg_nzb',ierr)
          deallocate(num_ovlaps,stat=ierr)
          call utils_dealloc_check('sparse_create','num_ovlaps',ierr)

          ! Point new_spam to new structure in the library
          new_spam%lib = num_library

          ! rc2013: final check before we leave
          call sparse_par_check(new_spam%lib, spama%lib, spamb%lib, 'sparse_create')

          ! Allocate space for data
          call sparse_data_alloc(new_spam,loc_iscmplx)

       end if

    end if

    ! Return library index if requested
    if (present(rlib)) rlib = new_spam%lib

  contains

    !==========================================================================!
    ! This subroutine counts the nonzero elements and blocks in a structure    !
    ! derived from the product of two existing sparse matrices.                !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   nze    (output) : Number of nonzero elements                           !
    !   nzb    (output) : Number of nonzero blocks                             !
    !   my_nze (output) : Number of nonzero elements on this proc              !
    !   my_nzb (output) : Number of nonzero blocks on this proc                !
    !--------------------------------------------------------------------------!
    ! Written by Nicholas Hine, May 2009, based on the SPAM3 routine:          !
    ! Original version written by Peter Haynes, March 2004,                    !
    ! Revised for distributed data, June 2006,                                 !
    ! Modified for dense matrices by Nicholas Hine, Dec 2007.                  !
    ! Revised to use comms groups by Nicholas Hine, May 2013.                  !
    !==========================================================================!

    subroutine internal_count_product(nze,nzb,my_nze,my_nzb)

      implicit none

      ! Arguments
      integer(kind=LONG), intent(out) :: nze ! Number of nonzero elements
      integer, intent(out) :: nzb    ! Number of nonzero blocks
      integer, intent(out) :: my_nze ! Number of nonzero elements on this proc
      integer, intent(out) :: my_nzb ! Number of nonzero blocks on this proc

      ! Local variables
      integer :: col_blks             ! Identifier for column blocking scheme
      integer :: row_blks             ! Identifier for row blocking scheme
      integer :: step                 ! Proc step loop counter
      integer :: reqproc              ! Proc to request index from for next step
      integer :: recvproc             ! Proc to receive index from this step
      integer :: group_nat            ! Number of atoms in this comms group
      integer :: iblk,jblk,kblk       ! Block counters
      integer :: loc_iblk,loc_jblk    ! Block counters local to proc
      integer :: grp_iblk             ! Block counter local to group
      integer :: iovlap,jovlap        ! Block overlap counters
      integer :: jovlap_lo,jovlap_hi  ! Block overlap binary search counters
      integer :: ielems               ! Number of elements on atom iblk
      integer :: iblk_ovlaps          ! Number of overlaps involving atom iblk
      integer :: max_ovlaps           ! Maximum num overlaps found for any atom
      integer :: a_row_seg            ! Segment index
      integer :: aseg_start           ! Start position of this segment in index
      integer :: b_seg_type           ! Segment type for this segment of B
      integer :: b_col_seg            ! Segment index within B shared cols
      integer :: b_col_seg_proc       ! Proc count of B shared col in full list
      integer :: bseg_start           ! Start position of this segment in index
      integer, allocatable :: ovlap_buf(:) ! Ovlap list buffer
      integer, allocatable :: ovlap_list(:,:) ! Ovlap list for atoms in group
      integer, allocatable :: ovlap_list_grp(:,:,:) ! Ovlap list (comms group)
      integer, allocatable :: proc_nze(:) ! Numbers of nonzero elements on procs
      integer :: first_atom_in_group ! rc2013: First atom in this proc's comms group
      integer :: last_atom_in_group  ! rc2013: Last atom in this proc's comms group
!$    integer, external :: omp_get_thread_num

      ! First obtain a quick over-estimate of the number of overlaps
      ! involving each local atom

      ! Zero counters
      num_ovlaps = 0

      call comms_barrier
      call comms_free

      ! Initialise amat comms
      acom%index_reqs(:) = index_not_sent
      call sparse_init_comms(acom,grouped=.true.)

!$OMP PARALLEL NUM_THREADS(pub_threads_max) DEFAULT(none) &
!$OMP PRIVATE(iblk,loc_iblk,iblk_ovlaps,iovlap,jblk,a_row_seg, &
!$OMP      loc_jblk,jovlap,kblk,grp_iblk,b_col_seg,b_col_seg_proc,bseg_start, &
!$OMP      b_seg_type,aseg_start,reqproc,step,tfound_idx,tfound_blk,ierr) &
!$OMP SHARED(acom,recvproc,pub_total_num_procs,row_par,col_par,bshare, &
!$OMP      num_ovlaps,pub_first_proc_in_group,pub_last_proc_in_group, &
!$OMP      pub_num_comms_groups,pub_my_proc_id,acol_par, &
!$OMP      alib,library,pub_comms_group_size,pub_threads_max, &
!$OMP      first_atom_in_group, last_atom_in_group, relevant_group)

      allocate(tfound_idx(library(alib)%nblk),stat=ierr)
      call utils_alloc_check('internal_count_product (sparse_create)', &
           'tfound_idx',ierr)
      allocate(tfound_blk(library(alib)%nblk),stat=ierr)
      call utils_alloc_check('internal_count_product (sparse_create)', &
           'tfound_blk',ierr)
      tfound_blk = .false.

      ! Loop over other processors for whole columns of A
      do step=0,pub_num_comms_groups-1

!$OMP BARRIER
!$       if (omp_get_thread_num()==0) then
         ! Receive data for this step
         recvproc = modulo(pub_my_proc_id-step*pub_comms_group_size + &
              pub_total_num_procs,pub_total_num_procs)
         call sparse_get_step_index(acom,recvproc)

         ! Request index and data for next step if required
         reqproc = modulo(pub_my_proc_id-(step+1)*pub_comms_group_size + &
              pub_total_num_procs,pub_total_num_procs)
         if (reqproc/=pub_my_proc_id) then
            ! Send request for index, pointers and data
            call comms_send(reqproc,index_needed,1,tag=IDX_REQ_TAG)
            ! Start asynchronous receives for index and pointers
            call sparse_recv_index(acom,reqproc,1,async=.true.)
         end if

!$       end if
!$OMP FLUSH(acom,recvproc)
!$OMP BARRIER

         ! jcap: check if this group contains any procs that are
         ! assigned to the appropriate region
         if (relevant_group) then
            ! rc2013: find the first and last atoms in this comms group
            first_atom_in_group = sparse_first_atom_in_group(col_par, &
                 pub_first_proc_in_group, pub_last_proc_in_group)
            last_atom_in_group = sparse_last_atom_in_group(col_par, &
                 pub_first_proc_in_group, pub_last_proc_in_group)
         else
            ! If not, we should skip the rest of this iteration
            first_atom_in_group = 0
            last_atom_in_group = 0
            cycle
         end if

         ! Loop over atom block-columns iblk of B in this group
!$OMP DO
         ! rc2013: construct the row and column pars of new_spam from
         ! A and B respectively
         do iblk=first_atom_in_group,last_atom_in_group

            ! Find the proc index in the full list for this b_col_seg
            b_col_seg_proc = col_par%proc_of_atom(col_par%orig_atom(iblk))
            b_col_seg = b_col_seg_proc - pub_first_proc_in_group
            bseg_start = bshare%seginfobuf(s_idx,recvproc,b_col_seg)

            ! Find segment type of this part of index of B and cycle if blank
            b_seg_type = bshare%seginfobuf(s_type,recvproc,b_col_seg)
            if (b_seg_type == SEG_BLANK) cycle

            loc_iblk = iblk - col_par%first_atom_on_proc(b_col_seg_proc) + 1
            grp_iblk = iblk - first_atom_in_group + 1

            ! Reset counter for number of overlaps for atom iblk
            iblk_ovlaps = 0

            ! Loop over atoms jblk overlapping atom iblk: block-rows in
            ! block-col iblk of B
            do iovlap=bshare%idxbuf(bseg_start+loc_iblk-1,b_col_seg), &
                 bshare%idxbuf(bseg_start+loc_iblk,b_col_seg)-1
               jblk = bshare%idxbuf(iovlap,b_col_seg)
               loc_jblk = jblk - acol_par%first_atom_on_proc(recvproc) + 1

               ! Loop over atoms kblk overlapping atom jblk: block-rows in
               ! block-column jblk of A
               do a_row_seg=0,pub_total_num_procs-1

                  ! If this segment of A is blank, move on to next
                  if (acom%seginfobuf(s_type,a_row_seg,2)==SEG_BLANK) cycle

                  aseg_start = acom%seginfobuf(s_idx,a_row_seg,2)

                  do jovlap=acom%idxbuf(aseg_start+loc_jblk-1,2), &
                       acom%idxbuf(aseg_start+loc_jblk,2)-1
                     kblk = acom%idxbuf(jovlap,2)

                     ! Check whether this sph-sph overlap has already been found
                     if (.not. tfound_blk(kblk)) then

                        ! Mark block as found
                        tfound_blk(kblk) = .true.
                        iblk_ovlaps = iblk_ovlaps + 1
                        tfound_idx(iblk_ovlaps) = kblk

                     end if

                  end do  ! kblk

               end do ! a_row_seg

            end do  ! jblk

            ! Keep count of total number of overlaps involving iblk found
            num_ovlaps(grp_iblk) = num_ovlaps(grp_iblk) + iblk_ovlaps

            ! Reset flags for overlaps found
            do iovlap=1,iblk_ovlaps
               kblk = tfound_idx(iovlap)
               tfound_blk(kblk) = .false.
            end do

!$          if (omp_get_thread_num()==0) then
            call sparse_check_send_requests(acom)
!$          end if

         end do  ! Loop over atoms iblk on this proc
!$OMP END DO

      end do  ! Loop over procs

      deallocate(tfound_blk,stat=ierr)
      call utils_dealloc_check('internal_count_product (sparse_create)', &
           'tfound_blk',ierr)
      deallocate(tfound_idx,stat=ierr)
      call utils_dealloc_check('internal_count_product (sparse_create)', &
           'tfound_idx',ierr)

!$OMP END PARALLEL

      ! Finalise amat comms
      call sparse_exit_comms(acom,spama)

      ! Reduce overlap count over all procs in comms group
      call comms_reduce('SUM',num_ovlaps,comm=pub_group_comm)

      ! Set up holding array to store overlap details
      max_ovlaps = min(maxval(num_ovlaps) + 1,library(alib)%nblk+1)
      group_nat = last_atom_in_group - first_atom_in_group + 1
      allocate(ovlap_list(max_ovlaps,group_nat),stat=ierr)
      call utils_alloc_check('internal_count_product (sparse_create)', &
           'ovlap_list',ierr)
      ovlap_list = 0

      ! Now count again, avoiding over-counting between processors

      ! Zero counters
      my_nze = 0
      my_nzb = 0
      num_ovlaps = 0
      seg_nzb(:) = 0
      seg_nze(:) = 0
      call comms_barrier

      ! Find block sizes from the product matrices
      col_blks = library(num_library)%col_blks
      row_blks = library(num_library)%row_blks

      ! Initialise amat comms
      acom%index_reqs(:) = index_not_sent
      call sparse_init_comms(acom,grouped=.true.)

!$OMP PARALLEL NUM_THREADS(pub_threads_max) DEFAULT(none) &
!$OMP PRIVATE(iblk,loc_iblk,iblk_ovlaps,iovlap,jblk,loc_jblk, &
!$OMP      a_row_seg,jovlap,kblk,ielems,grp_iblk,b_col_seg,b_col_seg_proc, &
!$OMP      bseg_start,b_seg_type,jovlap_hi,jovlap_lo,aseg_start,ovlap_buf, &
!$OMP      reqproc,tfound_blk,tfound_idx,ierr,step) &
!$OMP SHARED(acom,recvproc,pub_total_num_procs,alib,library, &
!$OMP      row_par,col_par,acol_par,num_ovlaps, &
!$OMP      pub_my_rank_in_group, &
!$OMP      row_blks,col_blks,ovlap_list,bshare,pub_my_proc_id, &
!$OMP      pub_first_proc_in_group,pub_last_proc_in_group, &
!$OMP      pub_num_comms_groups,pub_comms_group_size,max_ovlaps,pub_threads_max, &
!$OMP      first_atom_in_group, last_atom_in_group, relevant_group)

      allocate(tfound_idx(library(alib)%nblk),stat=ierr)
      call utils_alloc_check('internal_count_product (sparse_create)', &
           'tfound_idx',ierr)
      allocate(tfound_blk(library(alib)%nblk),stat=ierr)
      call utils_alloc_check('internal_count_product (sparse_create)', &
           'tfound_blk',ierr)
      allocate(ovlap_buf(max_ovlaps),stat=ierr)
      call utils_alloc_check('internal_count_product (sparse_create)', &
           'ovlap_buf',ierr)
      tfound_blk = .false.

      ! Loop over processors
      do step=0,pub_num_comms_groups-1

!$OMP BARRIER
!$       if (omp_get_thread_num()==0) then
         ! Receive data for this step
         recvproc = modulo(pub_my_proc_id-step*pub_comms_group_size + &
              pub_total_num_procs,pub_total_num_procs)
         call sparse_get_step_index(acom,recvproc)

         ! Request index and data for next step if required
         reqproc = modulo(pub_my_proc_id-(step+1)*pub_comms_group_size + &
              pub_total_num_procs,pub_total_num_procs)
         if (reqproc/=pub_my_proc_id) then
            ! Send request for index, pointers and data
            call comms_send(reqproc,index_needed,1,tag=IDX_REQ_TAG)
            ! Start asynchronous receives for index and pointers
            call sparse_recv_index(acom,reqproc,1,async=.true.)
         end if
!$       end if
!$OMP FLUSH(acom,recvproc)
!$OMP BARRIER

         ! jcap: Again, if this group does not contain any procs that
         ! are associated with the column region of A, cycle
         if (.not.relevant_group) cycle

         ! Loop over atom block-columns iblk on this proc
!$OMP DO
         do iblk=first_atom_in_group,last_atom_in_group

            ! Find the proc index in the full list for this b_col_seg
            b_col_seg_proc = col_par%proc_of_atom(col_par%orig_atom(iblk))
            b_col_seg = b_col_seg_proc - pub_first_proc_in_group
            bseg_start = bshare%seginfobuf(s_idx,recvproc,b_col_seg)

            ! Find segment type of this part of index of B and cycle if blank
            b_seg_type = bshare%seginfobuf(s_type,recvproc,b_col_seg)
            if (b_seg_type == SEG_BLANK) cycle

            loc_iblk = iblk - col_par%first_atom_on_proc(b_col_seg_proc) + 1
            grp_iblk = iblk - first_atom_in_group + 1

            ! Reset counter for number of overlaps for atom iblk
            iblk_ovlaps = 0

            ! Number of elements on atom iblk
            ielems = col_par%num_elems_on_atom(iblk,col_blks)

            ! Loop over atoms jblk overlapping atom iblk: block-rows in
            ! block-col iblk of B
            do iovlap=bshare%idxbuf(bseg_start+loc_iblk-1,b_col_seg), &
                 bshare%idxbuf(bseg_start+loc_iblk,b_col_seg)-1
               jblk = bshare%idxbuf(iovlap,b_col_seg)

               ! Loop over atoms kblk overlapping atom jblk: block-rows in
               ! block-column jblk of A
               do a_row_seg=0,pub_total_num_procs-1

                  ! If this segment of A is blank, move on to next
                  if (acom%seginfobuf(s_type,a_row_seg,2)==SEG_BLANK) cycle

                  aseg_start = acom%seginfobuf(s_idx,a_row_seg,2)

                  loc_jblk = jblk - acol_par%first_atom_on_proc(recvproc) + 1
                  do jovlap=acom%idxbuf(aseg_start+loc_jblk-1,2), &
                       acom%idxbuf(aseg_start+loc_jblk,2)-1
                     kblk = acom%idxbuf(jovlap,2)

                     ! Check whether this sph-sph overlap has already been found
                     if (.not. tfound_blk(kblk)) then

                        ! Mark block as found
                        tfound_blk(kblk) = .true.
                        iblk_ovlaps = iblk_ovlaps + 1
                        tfound_idx(iblk_ovlaps) = kblk

                     end if

                  end do  ! kblk

               end do ! a_row_seg

            end do  ! jblk

            ! Loop over overlaps found for block iblk (i.e. block-rows kblk in
            ! block-column iblk of the product AB)
            do iovlap=1,iblk_ovlaps
               kblk = tfound_idx(iovlap)

               ! Find position for kblk in the list using binary search
               jovlap = 1
               jovlap_lo = 1
               jblk = 0
               jovlap_hi = num_ovlaps(grp_iblk)
               do
                  if (jovlap_lo>jovlap_hi) exit
                  jovlap = (jovlap_hi + jovlap_lo)/2
                  jblk = ovlap_list(jovlap,grp_iblk)
                  if (jblk==kblk) exit
                  if (jblk>kblk) jovlap_hi = jovlap - 1
                  if (jblk<kblk) jovlap_lo = jovlap + 1
               end do
               if (jblk/=kblk) jovlap = jovlap_hi + 1

               ! Check that kblk does not already occur in the list
               if (.not. ovlap_list(jovlap,grp_iblk) == kblk) then

                  ! Increase number of overlaps found for iblk
                  num_ovlaps(grp_iblk) = num_ovlaps(grp_iblk) + 1

                  ! Shuffle up entries for atoms >kblk
                  ovlap_buf(jovlap+1:num_ovlaps(grp_iblk)) &
                       = ovlap_list(jovlap:num_ovlaps(grp_iblk)-1,grp_iblk)
                  ovlap_list(jovlap+1:num_ovlaps(grp_iblk),grp_iblk) &
                       = ovlap_buf(jovlap+1:num_ovlaps(grp_iblk))

                  ! Insert entry for kblk
                  ovlap_list(jovlap,grp_iblk) = kblk

               end if

               ! Reset flags
               tfound_blk(kblk) = .false.

            end do

!$          if (omp_get_thread_num()==0) then
            call sparse_check_send_requests(acom)
!$          end if

         end do  ! Loop over atoms iblk in this group
!$OMP END DO

      end do  ! Loop over procs

      deallocate(ovlap_buf,stat=ierr)
      call utils_dealloc_check('internal_count_product (sparse_create)', &
           'ovlap_buf',ierr)
      deallocate(tfound_blk,stat=ierr)
      call utils_dealloc_check('internal_count_product (sparse_create)', &
           'tfound_blk',ierr)
      deallocate(tfound_idx,stat=ierr)
      call utils_dealloc_check('internal_count_product (sparse_create)', &
           'tfound_idx',ierr)

!$OMP END PARALLEL

      ! Finalise amat comms
      call sparse_exit_comms(acom,spama)

      call comms_barrier

      allocate(ovlap_list_grp(max_ovlaps,group_nat,0:pub_comms_group_size-1), &
           stat=ierr)
      call utils_alloc_check('internal_count_product (sparse_create)', &
           'ovlap_list_grp',ierr)

      allocate(tfound_blk(library(alib)%nblk),stat=ierr)
      call utils_alloc_check('internal_count_product (sparse_create)', &
           'tfound_blk',ierr)
      allocate(tfound_idx(library(alib)%nblk),stat=ierr)
      call utils_alloc_check('internal_count_product (sparse_create)', &
           'tfound_idx',ierr)
      tfound_blk = .false.

      ! Gather overlap lists from all procs in group
      call comms_allgather(ovlap_list_grp(:,:,0),ovlap_list,comm=pub_group_comm)

      ! Loop over atoms on this proc, find overlaps from any rank in group
      num_ovlaps = 0
      ! jcap: Again, if this group does not contain any procs that
      ! are associated with the column region of A, skip this bit
      if (relevant_group) then

         do iblk=col_par%my_first_blk,col_par%my_last_blk
            grp_iblk = iblk - first_atom_in_group + 1

            ! Reset number of overlaps for this blk column
            iblk_ovlaps = 0

            ! Number of elements on atom iblk
            ielems = col_par%num_elems_on_atom(iblk,col_blks)

            ! Loop over procs in this comms group
            do recvproc=0,pub_comms_group_size-1

               ! Loop over overlaps found on this proc for this atom
               do iovlap=1,size(ovlap_list_grp,1)

                  if (ovlap_list_grp(iovlap,grp_iblk,recvproc)==0) exit

                  kblk = ovlap_list_grp(iovlap,grp_iblk,recvproc)

                  if (.not. tfound_blk(kblk)) then

                     ! Mark block as found
                     tfound_blk(kblk) = .true.
                     iblk_ovlaps = iblk_ovlaps + 1
                     tfound_idx(iblk_ovlaps) = kblk

                  end if
               end do
            end do

            ! Sort list of overlaps to ensure it is in ascending order
            call internal_heapsort(iblk_ovlaps,tfound_idx)

            ! Now loop over the overlaps we just found
            do iovlap=1,iblk_ovlaps

               kblk = tfound_idx(iovlap)

               ! Increase number of overlaps found for iblk
               num_ovlaps(grp_iblk) = num_ovlaps(grp_iblk) + 1

               ! Count elements and blocks
               my_nze = my_nze + ielems * row_par%num_elems_on_atom(kblk,row_blks)
               my_nzb = my_nzb + 1

               ! Count contributions to this segment of product AB
               a_row_seg = row_par%proc_of_elem(row_par%first_elem_on_atom(kblk,row_blks),row_blks)
               seg_nzb(a_row_seg) = seg_nzb(a_row_seg) + 1
               seg_nze(a_row_seg) = seg_nze(a_row_seg) + &
                    ielems * row_par%num_elems_on_atom(kblk,row_blks)

               ! Reset flags
               tfound_blk(kblk) = .false.

            end do
         end do

      end if

      ! Gather overlap counts from all procs in group and find max
      call comms_reduce('MAX',num_ovlaps,comm=pub_group_comm)

      deallocate(tfound_idx,stat=ierr)
      call utils_dealloc_check('internal_count_product (sparse_create)', &
           'tfound_idx',ierr)
      deallocate(tfound_blk,stat=ierr)
      call utils_dealloc_check('internal_count_product (sparse_create)', &
           'tfound_blk',ierr)

      ! Allocate array to calculate non-overflowing nze
      allocate(proc_nze(0:pub_total_num_procs-1),stat=ierr)
      call utils_alloc_check('internal_count_product (sparse_create)', &
           'proc_nze',ierr)

      ! ndmh: the following calculates nze allowing for it to exceed the
      ! length of a standard integer. Assumes only nze may require kind=LONG,
      ! but that the individual my_nze values can be stored as standard
      ! integers

      ! Get proc_nze from each proc
      proc_nze(:) = 0
      proc_nze(pub_my_proc_id) = my_nze
      call comms_reduce('SUM',proc_nze(0:pub_total_num_procs-1),&
           pub_total_num_procs)

      ! ndmh: non-overflowing way to add up nze (defined LONG)
      nze = 0
      do recvproc=0,pub_total_num_procs-1
         nze = nze + int(proc_nze(recvproc),kind=LONG)
      end do

      ! Dellocate array to calculate non-overflowing nze
      deallocate(proc_nze,stat=ierr)
      call utils_dealloc_check('internal_count_product (sparse_create)', &
           'proc_nze',ierr)

      ! Sum up nzb over all procs
      acom%idxbuf(1,1) = my_nzb
      call comms_reduce('SUM',acom%idxbuf,1)
      nzb = acom%idxbuf(1,1)

      ! Deallocate holding array to store overlap details
      deallocate(ovlap_list_grp,stat=ierr)
      call utils_dealloc_check('internal_count_product (sparse_create)', &
           'ovlap_list_grp',ierr)
      deallocate(ovlap_list,stat=ierr)
      call utils_dealloc_check('internal_count_product (sparse_create)', &
           'ovlap_list',ierr)

    end subroutine internal_count_product

    !==========================================================================!
    ! This subroutine indexes the nonzero blocks in a structure derived from   !
    ! the product of two existing sparse matrices.                             !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   None                                                                   !
    !--------------------------------------------------------------------------!
    ! Written by Nicholas Hine, May 2009, based on the SPAM3 routine:          !
    ! Written by Peter Haynes, March 2004,                                     !
    ! Revised for distributed data, June 2006,                                 !
    ! Modified for dense matrices by Nicholas Hine, Dec 2007.                  !
    !==========================================================================!

    subroutine internal_index_product

      implicit none

      ! Local variables
      integer :: col_blks             ! Identifier for column blocking scheme
      integer :: row_blks             ! Identifier for row blocking scheme
      integer :: step                 ! Proc step loop counter
      integer :: reqproc              ! Proc to request index from for next step
      integer :: recvproc             ! Proc to receive index from this step
      integer :: my_nat               ! Number of atoms on this proc
      integer :: group_nat            ! Number of atoms in this comms group
      integer :: iblk,jblk,kblk       ! Atom block counters
      integer :: loc_iblk,loc_jblk    ! Atom counters local to proc
      integer :: grp_iblk             ! Block counter local to group
      integer :: iovlap,jovlap        ! Block overlap counters
      integer :: jovlap_lo,jovlap_hi  ! Block overlap binary search counters
      integer :: ielems               ! Number of elements on atom iblk
      integer :: iblk_ovlaps          ! Number of overlaps found for this atom
      integer :: max_ovlaps           ! Maximum num overlaps found for any atom
      integer :: idx                  ! Index counter
      integer :: ptr                  ! Pointer counter
      integer :: dptr                 ! Dense segment pointer
      integer :: a_row_seg            ! Segment index
      integer :: aseg_start           ! Start position of this segment in index
      integer :: b_seg_type           ! Segment type for this segment of B
      integer :: b_col_seg            ! Segment index within B shared cols
      integer :: b_col_seg_proc       ! Proc count of B shared col in full list
      integer :: bseg_start           ! Start position of this segment in index
      integer :: c_seg_type           ! Segment type for this segment
      integer :: c_row_seg            ! Segment index
      integer :: cseg_start           ! Start position of this segment in index
      integer, allocatable :: ovlap_buf(:) ! Ovlap list buffer
      integer, allocatable :: ovlap_list(:,:) ! Ovlap list for atoms in group
      integer, allocatable :: ovlap_list_grp(:,:,:) ! Ovlap list (comms group)
!$    integer, external :: omp_get_thread_num

      ! Find block sizes from the product matrix
      col_blks = library(num_library)%col_blks
      row_blks = library(num_library)%row_blks

      ! Set up holding array to store overlap details
      max_ovlaps = min(maxval(num_ovlaps) + 1,library(alib)%nblk+1)
      group_nat = last_atom_in_group - first_atom_in_group +1
      allocate(ovlap_list(max_ovlaps,group_nat),stat=ierr)
      call utils_alloc_check('internal_index_product (sparse_create)', &
           'ovlap_list',ierr)
      ovlap_list = 0
      num_ovlaps = 0
      dptr = -1

      ! Initialise amat comms
      acom%index_reqs(:) = index_not_sent
      call sparse_init_comms(acom,grouped=.true.)

!$OMP PARALLEL NUM_THREADS(pub_threads_max) DEFAULT(none) &
!$OMP PRIVATE(iblk,loc_iblk,iblk_ovlaps,iovlap,jblk,a_row_seg,loc_jblk, &
!$OMP      jovlap,kblk,grp_iblk,b_col_seg,b_col_seg_proc,bseg_start, &
!$OMP      b_seg_type,jovlap_lo,jovlap_hi,aseg_start,reqproc,ovlap_buf, &
!$OMP      ierr,step,tfound_blk,tfound_idx) &
!$OMP SHARED(acom,pub_total_num_procs,alib,library, &
!$OMP      row_par,col_par,acol_par,num_ovlaps, &
!$OMP      row_blks,col_blks,ovlap_list,pub_my_rank_in_group,pub_my_proc_id, &
!$OMP      pub_first_proc_in_group,pub_last_proc_in_group, &
!$OMP      bshare,pub_num_comms_groups,pub_comms_group_size,max_ovlaps, &
!$OMP      recvproc,pub_threads_max, &
!$OMP      first_atom_in_group, last_atom_in_group, relevant_group)

      allocate(tfound_idx(library(alib)%nblk),stat=ierr)
      call utils_alloc_check('internal_index_product (sparse_create)', &
           'tfound_idx',ierr)
      allocate(tfound_blk(library(alib)%nblk),stat=ierr)
      call utils_alloc_check('internal_index_product (sparse_create)', &
           'tfound_blk',ierr)
      allocate(ovlap_buf(max_ovlaps),stat=ierr)
      call utils_alloc_check('internal_index_product (sparse_create)', &
           'ovlap_buf',ierr)
      tfound_blk = .false.

      ! Loop over processors
      do step=0,pub_num_comms_groups-1

!$OMP BARRIER
!$       if (omp_get_thread_num()==0) then
         ! Receive data for this step
         recvproc = modulo(pub_my_proc_id-step*pub_comms_group_size + &
              pub_total_num_procs,pub_total_num_procs)
         call sparse_get_step_index(acom,recvproc)

         ! Request index and data for next step if required
         reqproc = modulo(pub_my_proc_id-(step+1)*pub_comms_group_size + &
              pub_total_num_procs,pub_total_num_procs)
         if (reqproc/=pub_my_proc_id) then
            ! Send request for index, pointers and data
            call comms_send(reqproc,index_needed,1,tag=IDX_REQ_TAG)
            ! Start asynchronous receives for index and pointers
            call sparse_recv_index(acom,reqproc,1,async=.true.)
         end if
!$       end if
!$OMP FLUSH(acom,recvproc)
!$OMP BARRIER

         ! jcap: Skip this iteration if the group is not relevant
         if (.not.relevant_group) cycle

         ! Loop over atom block-columns iblk on this proc
!$OMP DO
         do iblk=first_atom_in_group,last_atom_in_group

            ! Find the proc index in the full list for this b_col_seg
            b_col_seg_proc = col_par%proc_of_atom(col_par%orig_atom(iblk))
            b_col_seg = b_col_seg_proc - pub_first_proc_in_group
            bseg_start = bshare%seginfobuf(s_idx,recvproc,b_col_seg)

            ! Find segment type of this part of index of B and cycle if blank
            b_seg_type = bshare%seginfobuf(s_type,recvproc,b_col_seg)
            if (b_seg_type == SEG_BLANK) cycle

            loc_iblk = iblk - col_par%first_atom_on_proc(b_col_seg_proc) + 1
            grp_iblk = iblk - first_atom_in_group + 1

            ! Reset counter for number of overlaps for atom iblk
            iblk_ovlaps = 0

            ! Loop over atoms jblk overlapping atom iblk according to spamb
            do iovlap=bshare%idxbuf(bseg_start+loc_iblk-1,b_col_seg), &
                 bshare%idxbuf(bseg_start+loc_iblk,b_col_seg)-1
               jblk = bshare%idxbuf(iovlap,b_col_seg)

               ! Loop over atoms kblk overlapping atom jblk: block-rows in
               ! block-column jblk of spama
               do a_row_seg=0,pub_total_num_procs-1
                  loc_jblk = jblk - acol_par%first_atom_on_proc(recvproc) + 1

                  ! If this segment of A is blank, move on
                  if (acom%seginfobuf(s_type,a_row_seg,2)==SEG_BLANK) cycle

                  aseg_start = acom%seginfobuf(s_idx,a_row_seg,2)

                  do jovlap=acom%idxbuf(aseg_start+loc_jblk-1,2), &
                       acom%idxbuf(aseg_start+loc_jblk,2)-1
                     kblk = acom%idxbuf(jovlap,2)

                     ! Check whether this sph-sph overlap has already been found
                     if (.not. tfound_blk(kblk)) then

                        ! Mark block as found
                        tfound_blk(kblk) = .true.
                        iblk_ovlaps = iblk_ovlaps + 1
                        tfound_idx(iblk_ovlaps) = kblk

                     end if

                  end do  ! kblk

               end do  ! a_row_seg

            end do  ! jblk

            ! Loop over overlaps found for atom iblk
            do iovlap=1,iblk_ovlaps
               kblk = tfound_idx(iovlap)

               ! Find position for kblk in the list using binary search
               jovlap = 1
               jovlap_lo = 1
               jblk = 0
               jovlap_hi = num_ovlaps(grp_iblk)
               do
                  if (jovlap_lo>jovlap_hi) exit
                  jovlap = (jovlap_hi + jovlap_lo)/2
                  jblk = ovlap_list(jovlap,grp_iblk)
                  if (jblk==kblk) exit
                  if (jblk>kblk) jovlap_hi = jovlap - 1
                  if (jblk<kblk) jovlap_lo = jovlap + 1
               end do
               if (jblk/=kblk) jovlap = jovlap_hi + 1

               ! Check that kblk does not already occur in the list
               if (.not. ovlap_list(jovlap,grp_iblk) == kblk) then

                  ! Increase number of overlaps found for iblk
                  num_ovlaps(grp_iblk) = num_ovlaps(grp_iblk) + 1

                  ! Shuffle up entries for atoms >kblk
                  ovlap_buf(jovlap+1:num_ovlaps(grp_iblk)) &
                       = ovlap_list(jovlap:num_ovlaps(grp_iblk)-1,grp_iblk)
                  ovlap_list(jovlap+1:num_ovlaps(grp_iblk),grp_iblk) &
                       = ovlap_buf(jovlap+1:num_ovlaps(grp_iblk))

                  ! Insert entry for kblk
                  ovlap_list(jovlap,grp_iblk) = kblk

               end if

               ! Reset flag
               tfound_blk(kblk) = .false.

            end do  ! Overlaps for atom iblk

!$          if (omp_get_thread_num()==0) then
            call sparse_check_send_requests(acom)
!$          end if

         end do  ! Loop over atoms iblk
!$OMP END DO
      end do  ! Loop over steps

      deallocate(ovlap_buf,stat=ierr)
      call utils_dealloc_check('internal_index_product (sparse_create)', &
           'ovlap_buf',ierr)
      deallocate(tfound_blk,stat=ierr)
      call utils_dealloc_check('internal_index_product (sparse_create)', &
           'tfound_blk',ierr)
      deallocate(tfound_idx,stat=ierr)
      call utils_dealloc_check('internal_index_product (sparse_create)', &
           'tfound_idx',ierr)
!$OMP END PARALLEL

      ! Finalise amat comms
      call sparse_exit_comms(acom,spama)

      allocate(ovlap_list_grp(max_ovlaps,group_nat,0:pub_comms_group_size-1), &
           stat=ierr)
      call utils_alloc_check('internal_index_product (sparse_create)', &
           'ovlap_list_grp',ierr)

      allocate(tfound_blk(library(alib)%nblk),stat=ierr)
      call utils_alloc_check('internal_index_product (sparse_create)', &
           'tfound_blk',ierr)
      allocate(tfound_idx(library(alib)%nblk),stat=ierr)
      call utils_alloc_check('internal_index_product (sparse_create)', &
           'tfound_idx',ierr)
      tfound_blk = .false.

      ! Gather overlap lists from all procs in group
      call comms_allgather(ovlap_list_grp(:,:,0),ovlap_list,comm=pub_group_comm)

      ! Initialise counters and pointers
      my_nat = col_par%my_last_blk - col_par%my_first_blk + 1
      idx = 1
      ptr = 1

      ! Loop over segments of the index of the product AB
      do c_row_seg=0,pub_total_num_procs-1

         ! Set up entries in the seg_info(s_idx,:) and seg_info(s_ptr,:) arrays
         library(num_library)%seg_info(s_ptr,c_row_seg) = ptr
         library(num_library)%seg_info(s_idx,c_row_seg) = idx
         c_seg_type = library(num_library)%seg_info(s_type,c_row_seg)
         cseg_start = idx

         ! Skip indexing of this segment if it contains no nonzero elements
         if (c_seg_type == SEG_BLANK) cycle

         if (c_row_seg/=pub_my_proc_id) then
            ! Set unused pointers to zero
            library(num_library)%blk_ptr(idx:idx+my_nat) = 0
         end if

         ! Go past list of column start positions for this segment
         idx = idx + my_nat + 1
         library(num_library)%blk_idx(cseg_start) = idx

         ! Loop over atom block-columns iblk on this proc for this segment
         do iblk=col_par%my_first_blk,col_par%my_last_blk
            loc_iblk = iblk - col_par%my_first_blk + 1
            grp_iblk = iblk - first_atom_in_group + 1

            ! Reset number of overlaps for this blk column
            iblk_ovlaps = 0

            ! Number of elements on atom iblk
            ielems = col_par%num_elems_on_atom(iblk,col_blks)

            ! Loop over procs in this comms group
            do recvproc=0,pub_comms_group_size-1

               ! Loop over overlaps found on this proc for this atom
               do iovlap=1,size(ovlap_list,1)

                  if (ovlap_list_grp(iovlap,grp_iblk,recvproc)==0) exit

                  kblk = ovlap_list_grp(iovlap,grp_iblk,recvproc)

                  !if (kblk < par%first_atom_on_proc(c_row_seg)) cycle
                  !if (kblk >= par%first_atom_on_proc(c_row_seg+1)) exit
                  if (kblk < row_par%first_atom_on_proc(c_row_seg)) cycle
                  if (kblk >= row_par%first_atom_on_proc(c_row_seg+1)) exit

                  if (.not. tfound_blk(kblk)) then

                     ! Mark block as found
                     tfound_blk(kblk) = .true.
                     iblk_ovlaps = iblk_ovlaps + 1
                     tfound_idx(iblk_ovlaps) = kblk

                  end if
               end do
            end do

            ! Sort list of overlaps to ensure it is in ascending order
            call internal_heapsort(iblk_ovlaps,tfound_idx)

            ! Loop over overlaps found with this atom
            do iovlap=1,iblk_ovlaps
               jblk = tfound_idx(iovlap)

               ! Index contains block-row
               library(num_library)%blk_idx(idx) = jblk

               ! Pointer points to start of data
               if ( c_seg_type == SEG_SPARSE ) then
                  dptr = ptr
                  ! Increment block pointer by size of block
                  ptr = ptr + ielems * row_par%num_elems_on_atom(jblk,row_blks)
               else if ( c_seg_type == SEG_DENSE ) then
                  ! Find ptr to first element of this block in segment
                  dptr = ptr + (col_par%first_elem_on_atom(iblk,col_blks) - &
                       col_par%first_elem_on_proc(pub_my_proc_id,col_blks)) * &
                       row_par%num_elems_on_proc(c_row_seg,row_blks) &
                       + row_par%first_elem_on_atom(jblk,row_blks) - &
                       row_par%first_elem_on_proc(c_row_seg,row_blks)
               end if

               library(num_library)%blk_ptr(idx) = dptr
               if (jblk == iblk) &
                    library(num_library)%blk_ptr(cseg_start+loc_iblk-1) = dptr

               ! Increment block index by one
               idx = idx + 1

               ! Reset flags
               tfound_blk(jblk) = .false.

            end do  ! jblk

            ! Store idx in list of column start positions for this segment
            library(num_library)%blk_idx(cseg_start+loc_iblk) = idx
            library(num_library)%blk_ptr(cseg_start+loc_iblk) = 0

         end do  ! iblk

         ! For dense segments, move ptr on by total number of elements
         if ( c_seg_type == SEG_DENSE ) then
            ptr = ptr + col_par%num_elems_on_proc(pub_my_proc_id,col_blks) * &
                 row_par%num_elems_on_proc(c_row_seg,row_blks)
         end if

         ! Set final pointer to end of data
         library(num_library)%blk_ptr(idx) = ptr

      end do  ! c_row_seg

      ! Record end of last segment
      library(num_library)%seg_info(s_idx,pub_total_num_procs) = idx
      library(num_library)%seg_info(s_ptr,pub_total_num_procs) = ptr

      deallocate(tfound_blk,stat=ierr)
      call utils_dealloc_check('internal_index_product (sparse_create)', &
           'tfound_blk',ierr)
      deallocate(tfound_idx,stat=ierr)
      call utils_dealloc_check('internal_index_product (sparse_create)', &
           'tfound_idx',ierr)

      ! Deallocate holding array to store overlap details
      deallocate(ovlap_list_grp,stat=ierr)
      call utils_dealloc_check('internal_index_product (sparse_create)', &
           'ovlap_list_grp',ierr)
      deallocate(ovlap_list,stat=ierr)
      call utils_dealloc_check('internal_index_product (sparse_create)', &
           'ovlap_list',ierr)

      ! Do not exit until all sends are completed
      call comms_free

      my_nat = pub_total_num_procs - &
           count(library(num_library)%seg_info(s_type,:)==SEG_BLANK)

      iovlap = (library(num_library)%my_nblks + 1)*my_nat + &
           library(num_library)%my_nzb
      call utils_assert(iovlap+1 == idx, 'Error in internal_index_product: &
           &Consistency check failed. The values of iovlap, idx-1 and proc &
           &number follow. ',iovlap, idx-1, pub_my_proc_id)

    end subroutine internal_index_product

  end subroutine sparse_create
  ! rc2013: end of sparse_create


  !============================================================================!
  ! This subroutine destroys a sparse matrix structure, freeing up the memory  !
  ! allocated to internal arrays.                                              !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   old_spam (inout) : The sparse matrix to be destroyed                     !
  !----------------------------------------------------------------------------!
  ! Written by Peter Haynes, March 2004.                                       !
  ! Revised for distributed data, June 2006.                                   !
  ! Updated for SPAM3 by Nicholas Hine, May 2009.                              !
  !============================================================================!

  subroutine sparse_destroy(old_spam)

    implicit none

    ! Arguments
    type(SPAM3), intent(inout) :: old_spam

    ! Deallocate memory assigned to internal structures
    call sparse_data_dealloc(old_spam)

  end subroutine sparse_destroy


  !============================================================================!
  ! This function returns the index length of a sparse matrix.                 !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   mat   (input)  : The matrix to be assessed                               !
  !----------------------------------------------------------------------------!
  ! Written by Peter Haynes, June 2006.                                        !
  !============================================================================!

  integer function sparse_index_length(mat)

    implicit none

    ! Arguments
    type(SPAM3), intent(in) :: mat   ! The sparse matrix

    ! Local variable
    integer :: ilib                  ! Library entry for mat

    ! Very simple...
    ilib = mat%lib
    sparse_index_length = library(ilib)%my_nblks + library(ilib)%my_nzb + 1

  end function sparse_index_length

  !============================================================================!
  ! This subroutine returns the index of a sparse matrix in a column-indexed,  !
  ! non-segmented form suitable for use with integrals_mod routines.           !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   idx   (output) : The array for the index                                 !
  !   mat   (input)  : The matrix to be assessed                               !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, May 2009.                                        !
  !============================================================================!

  subroutine sparse_generate_index(idx,mat)

    use comms, only: pub_total_num_procs
    use parallel_strategy, only: PARAL_INFO
    use utils, only: utils_assert

    implicit none

    ! Arguments
    integer, intent(out) :: idx(:)  ! The array for the index
    type(SPAM3), intent(in) :: mat  ! The sparse matrix

    ! Local variables
    integer :: ilib                 ! The library entry for mat
    integer :: ilen                 ! The length of the index
    integer :: seg                  ! Segment index
    integer :: seg_start            ! Start position of seg in index
    integer :: iovlap               ! Block overlap counter
    integer :: iblk                 ! Block-column counter
    integer :: loc_iblk             ! Block-column counter local to this proc
    integer :: jblk                 ! Block-row counter
    integer :: iidx                 ! Indices
    type(PARAL_INFO), pointer :: col_par ! rc2013: parallel strategy for columns

    ! Basic information about the matrix mat
    ilib = mat%lib
    ilen = library(ilib)%my_nblks + library(ilib)%my_nzb + 1
    col_par=>library(ilib)%col_par

    ! Check array is sufficiently long
    call utils_assert(size(idx) >= ilen, 'Error in sparse_generate_index: &
         &idx is not long enough for index.')

    ! Index is sorted into proc-proc segments but we want the whole index
    ! for a single proc's columns combined, so we need to reorganise it

    ! Loop over blocks on this proc
    iidx = col_par%my_last_blk - col_par%my_first_blk + 3
    loc_iblk = 0
    do iblk=col_par%my_first_blk,col_par%my_last_blk
       loc_iblk = loc_iblk + 1

       ! Store index counter in the initial list of column start positions
       idx(loc_iblk) = iidx

       ! Loop over proc-segments of the index
       do seg=0,pub_total_num_procs-1

          ! Nothing to add if this segment is blank
          if (library(ilib)%seg_info(s_type,seg)==SEG_BLANK) cycle

          ! Loop over overlaps for this atom in this segment
          seg_start = library(ilib)%seg_info(s_idx,seg)
          do iovlap=library(ilib)%blk_idx(seg_start+loc_iblk-1), &
               library(ilib)%blk_idx(seg_start+loc_iblk)-1

             ! Store this block in the new list
             jblk = library(ilib)%blk_idx(iovlap)
             idx(iidx) = jblk
             iidx = iidx + 1

          end do ! iovlap

       end do  ! seg

    end do  ! iblk
    idx(loc_iblk+1) = iidx

  end subroutine sparse_generate_index

  !============================================================================!
  ! This subroutine returns the structure code which corresponds to the        !
  ! transpose of the sparse matrix passed in. For symmetric structures this    !
  ! is the same as the structure of the matrix itself, but for rectangular     !
  ! matrices it may be different.                                              !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   mat   (input)  : The matrix to be transposed                             !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, March 2011                                       !
  !============================================================================!

  subroutine sparse_transpose_structure(struc_trans,mat)

    implicit none

    ! Arguments
    character(len=30), intent(out) :: struc_trans
    type(SPAM3), intent(in) :: mat   ! The sparse matrix

    ! Local variable
    integer :: ilib                  ! Library entry for mat

    ! Get transpose structure
    ilib = mat%lib
    struc_trans = library(ilib)%transpose_structure

  end subroutine sparse_transpose_structure

  !============================================================================!
  ! Returns the index of the first element on a given proc in the row or       !
  ! column blocking scheme of a given matrix.                                  !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   mat    (input) : The matrix                                              !
  !   proc   (input) : The proc for which the first element is required        !
  !   rowcol (input) : 'R' to use row blocks, 'C' to use column blocks         !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, June 2009.                                       !
  !============================================================================!

  integer function sparse_first_elem_on_proc(proc,mat,rowcol)

    implicit none

    ! Arguments
    type(SPAM3), intent(in) :: mat     ! The matrix
    integer, intent(in) :: proc        ! The proc
    character, intent(in) :: rowcol    ! Whether to use row or col blocks
    type(PARAL_INFO), pointer :: row_par
    type(PARAL_INFO), pointer :: col_par

    ! rc2013: assign parallel strategies
    row_par => library(mat%lib)%row_par
    col_par => library(mat%lib)%col_par

    if (rowcol == 'R') then
       sparse_first_elem_on_proc = &
            row_par%first_elem_on_proc(proc,library(mat%lib)%row_blks)
    else if (rowcol == 'C') then
       sparse_first_elem_on_proc = &
            col_par%first_elem_on_proc(proc,library(mat%lib)%col_blks)
    else
       sparse_first_elem_on_proc = 0
    end if

  end function sparse_first_elem_on_proc

  !============================================================================!
  ! Returns the index of the last element on a given proc in the row or column !
  ! blocking scheme of a given matrix.                                         !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   mat    (input) : The matrix                                              !
  !   proc   (input) : The proc for which the last element is required         !
  !   rowcol (input) : 'R' to use row blocks, 'C' to use column blocks         !
  !----------------------------------------------------------------------------!
  ! Written by Robert Charlton, 07/05/2018.                                    !
  !============================================================================!

  integer function sparse_last_elem_on_proc(proc,mat,rowcol)

    implicit none

    ! Arguments
    type(SPAM3), intent(in) :: mat     ! The matrix
    integer, intent(in) :: proc        ! The proc
    character, intent(in) :: rowcol    ! Whether to use row or col blocks

    ! rc2013: if there's nothing here then set last_elem to -1 (so we don't
    ! loop over anything)
    if(sparse_first_elem_on_proc(proc,mat,rowcol) == 0) then
       sparse_last_elem_on_proc = -1
    else
       sparse_last_elem_on_proc = sparse_first_elem_on_proc(proc,mat,rowcol) &
            + sparse_num_elems_on_proc(proc,mat,rowcol) - 1
    endif

  end function sparse_last_elem_on_proc

  !============================================================================!
  ! Returns the index of the first element on a given atom in the row or       !
  ! column blocking scheme of a given matrix.                                  !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   mat    (input) : The matrix                                              !
  !   atom   (input) : The atom for which the first element is required        !
  !   rowcol (input) : 'R' to use row blocks, 'C' to use column blocks         !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, June 2009.                                       !
  !============================================================================!

  integer function sparse_first_elem_on_atom(atom,mat,rowcol)

    implicit none

    ! Arguments
    type(SPAM3), intent(in) :: mat     ! The matrix
    integer, intent(in) :: atom        ! The atom
    character, intent(in) :: rowcol    ! Whether to use row or col blocks
    type(PARAL_INFO), pointer :: row_par
    type(PARAL_INFO), pointer :: col_par

    ! rc2013: assign parallel strategies
    row_par => library(mat%lib)%row_par
    col_par => library(mat%lib)%col_par

    if (rowcol == 'R') then
       sparse_first_elem_on_atom = &
            row_par%first_elem_on_atom(atom,library(mat%lib)%row_blks)
    else if (rowcol == 'C') then
       sparse_first_elem_on_atom = &
            col_par%first_elem_on_atom(atom,library(mat%lib)%col_blks)
    else
       sparse_first_elem_on_atom = 0
    end if

  end function sparse_first_elem_on_atom

  !============================================================================!
  ! Returns the number of elements on a given proc in the row or column        !
  ! blocking scheme of a given matrix.                                         !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   mat    (input) : The matrix                                              !
  !   proc   (input) : The proc for which the number of elements is required   !
  !   rowcol (input) : 'R' to use row blocks, 'C' to use column blocks         !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, June 2009.                                       !
  !============================================================================!

  integer function sparse_num_elems_on_proc(proc,mat,rowcol)

    implicit none

    ! Arguments
    type(SPAM3), intent(in) :: mat     ! The matrix
    integer, intent(in) :: proc        ! The proc
    character, intent(in) :: rowcol    ! Whether to use row or col blocks

    type(PARAL_INFO), pointer :: row_par
    type(PARAL_INFO), pointer :: col_par

    ! rc2013: assign parallel strategies
    row_par => library(mat%lib)%row_par
    col_par => library(mat%lib)%col_par

    if (rowcol == 'R') then
       sparse_num_elems_on_proc = &
            row_par%num_elems_on_proc(proc,library(mat%lib)%row_blks)
    else if (rowcol == 'C') then
       sparse_num_elems_on_proc = &
            col_par%num_elems_on_proc(proc,library(mat%lib)%col_blks)
    else
       sparse_num_elems_on_proc = 0
    end if

  end function sparse_num_elems_on_proc

  !============================================================================!
  ! Returns the number of elements on a given atom in the row or column        !
  ! blocking scheme of a given matrix.                                         !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   mat    (input) : The matrix                                              !
  !   atom   (input) : The atom for which the number of elements is required   !
  !   rowcol (input) : 'R' to use row blocks, 'C' to use column blocks         !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, June 2009.                                       !
  !============================================================================!

  integer function sparse_num_elems_on_atom_mat(atom,mat,rowcol)

    implicit none

    ! Arguments
    type(SPAM3), intent(in) :: mat     ! The matrix
    integer, intent(in) :: atom        ! The atom
    character, intent(in) :: rowcol    ! Whether to use row or col blocks

    type(PARAL_INFO), pointer :: row_par
    type(PARAL_INFO), pointer :: col_par

    ! rc2013: assign parallel strategies
    row_par => library(mat%lib)%row_par
    col_par => library(mat%lib)%col_par

    if (rowcol == 'R') then
       sparse_num_elems_on_atom_mat = &
            row_par%num_elems_on_atom(atom,library(mat%lib)%row_blks)
    else if (rowcol == 'C') then
       sparse_num_elems_on_atom_mat = &
            col_par%num_elems_on_atom(atom,library(mat%lib)%col_blks)
    else
       sparse_num_elems_on_atom_mat = 0
    end if

  end function sparse_num_elems_on_atom_mat

  !============================================================================!
  ! Returns the number of elements on a given atom in the row or column        !
  ! blocking scheme of a given matrix.                                         !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   mat    (input) : The matrix                                              !
  !   atom   (input) : The atom for which the number of elements is required   !
  !   rowcol (input) : 'R' to use row blocks, 'C' to use column blocks         !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, June 2009.                                       !
  !============================================================================!

  integer function sparse_num_elems_on_atom_blks(atom,blks,par)

    implicit none

    ! Arguments
    integer, intent(in) :: atom        ! The atom
    integer, intent(in) :: blks        ! The blocking scheme
    type(PARAL_INFO), intent(in), pointer :: par

    sparse_num_elems_on_atom_blks = par%num_elems_on_atom(atom,blks)

  end function sparse_num_elems_on_atom_blks

  !============================================================================!
  ! Returns the proc index of a given element in the row or column blocking    !
  ! scheme of a given matrix.                                                  !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   mat    (input) : The matrix                                              !
  !   elem   (input) : The element index required                              !
  !   rowcol (input) : 'R' to use row blocks, 'C' to use column blocks         !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, June 2009.                                       !
  !============================================================================!

  integer function sparse_proc_of_elem(elem,mat,rowcol)

    implicit none

    ! Arguments
    type(SPAM3), intent(in) :: mat     ! The matrix
    integer, intent(in) :: elem        ! The element
    character, intent(in) :: rowcol    ! Whether to use row or col blocks
    type(PARAL_INFO), pointer :: row_par
    type(PARAL_INFO), pointer :: col_par

    ! rc2013: assign parallel strategies
    row_par => library(mat%lib)%row_par
    col_par => library(mat%lib)%col_par

    if (rowcol == 'R') then
       sparse_proc_of_elem = &
            row_par%proc_of_elem(elem,library(mat%lib)%row_blks)
    else if (rowcol == 'C') then
       sparse_proc_of_elem = &
            col_par%proc_of_elem(elem,library(mat%lib)%col_blks)
    else
       sparse_proc_of_elem = 0
    end if

  end function sparse_proc_of_elem

  !============================================================================!
  ! Returns the atom index of a given element in the row or column blocking    !
  ! scheme of a given matrix.                                                  !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   mat    (input) : The matrix                                              !
  !   elem   (input) : The element index required                              !
  !   rowcol (input) : 'R' to use row blocks, 'C' to use column blocks         !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, June 2009.                                       !
  !============================================================================!

  integer function sparse_atom_of_elem(elem,mat,rowcol)

    implicit none

    ! Arguments
    integer, intent(in) :: elem        ! The element
    type(SPAM3), intent(in) :: mat     ! The matrix
    character, intent(in) :: rowcol      ! Whether to use row or col blocks
    type(PARAL_INFO), pointer :: row_par
    type(PARAL_INFO), pointer :: col_par

    ! rc2013: assign parallel strategies
    row_par => library(mat%lib)%row_par
    col_par => library(mat%lib)%col_par


    if (rowcol == 'R') then
       sparse_atom_of_elem = &
            row_par%atom_of_elem(elem,library(mat%lib)%row_blks)
    else if (rowcol == 'C') then
       sparse_atom_of_elem = &
            col_par%atom_of_elem(elem,library(mat%lib)%col_blks)
    else
       sparse_atom_of_elem = 0
    end if

  end function sparse_atom_of_elem

  !============================================================================!
  ! Tests whether a specified element is listed in the index of a sparse       !
  ! matrix.                                                                    !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   mat  (input)   : The matrix                                              !
  !   irow (input)   : The element index of the row required                   !
  !   jcol (input)   : The element index of the column required                !
  !----------------------------------------------------------------------------!
  ! Written by Peter Haynes, July 2006.                                        !
  ! Modified for SPAM3 by Nicholas Hine, May 2009.                             !
  !============================================================================!

  logical function sparse_element_exists(mat,irow,jcol)

    use comms, only: pub_my_proc_id
    use parallel_strategy, only: PARAL_INFO
    use rundat, only: pub_debug
    use utils, only: utils_assert

    implicit none

    ! Arguments
    type(SPAM3), intent(in) :: mat     ! The matrix
    integer, intent(in) :: irow        ! The element index of the row
    integer, intent(in) :: jcol        ! The element index of the column

    ! Local variables
    integer :: ilib      ! Library entry for mat
    integer :: iblk      ! The block-row to which the element belongs
    integer :: jblk      ! The block-col to which the element belongs
    integer :: kblk      ! Block counter
    integer :: iovlap    ! Block overlap counter
    integer :: loc_jblk  ! Block-column counter local to this proc
    integer :: seg       ! Segment index
    integer :: seg_start ! Start position of this segment in the index
    type(PARAL_INFO), pointer :: row_par ! parallel strategy for matrix rows
    type(PARAL_INFO), pointer :: col_par ! parallel strategy for matrix columns

    ! Get library entry
    ilib = mat%lib

    ! Check arguments
    if(pub_debug) then
       call utils_assert(irow >= 1 .and. irow <= library(ilib)%nrows, &
          'Error in sparse_element_exists: invalid row index.')
       call utils_assert(jcol >= 1 .and. jcol <= library(ilib)%mcols, &
          'Error in sparse_element_exists: invalid column index.')
    end if

    ! rc2013: assign the pointer for parallel strategies based on structure
    row_par => library(ilib)%row_par
    col_par => library(ilib)%col_par

    ! Obtain information about the relevant block for this element
    seg = row_par%proc_of_elem(irow,library(ilib)%row_blks)
    iblk = row_par%atom_of_elem(irow,library(ilib)%row_blks)
    jblk = col_par%atom_of_elem(jcol,library(ilib)%col_blks)

    ! rc2013: assign the pointer for col_par based on structure
    col_par => library(ilib)%col_par

    ! Check this block-column is local to this processor
    call utils_assert(jblk >= col_par%my_first_blk .and. jblk <= col_par%my_last_blk, &
         'Error in sparse_element_exists: requested column is not local to &
         &proc ', pub_my_proc_id)

    ! Assume not found
    sparse_element_exists = .false.

    ! Element definitely doesn't exist if segment is blank
    if (library(ilib)%seg_info(s_type,seg) == SEG_BLANK) return

    ! Find local index for the block-column
    loc_jblk = jblk - col_par%my_first_blk + 1

    ! Loop over the nonzero blocks in this segment
    seg_start = library(ilib)%seg_info(s_idx,seg)
    do iovlap=library(ilib)%blk_idx(seg_start+loc_jblk-1), &
         library(ilib)%blk_idx(seg_start+loc_jblk)-1
       kblk = library(ilib)%blk_idx(iovlap)   ! block-row
       if (kblk == iblk) then
          sparse_element_exists = .true.
          exit
       end if
    end do

  end function sparse_element_exists


  !============================================================================!
  ! This function returns true the denominator for the filling factor of a     !
  ! sparse matrix, based on the row and column blocking schemes.               !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   row_blks   (input)  : The row blocking scheme.                           !
  !   col_blks   (input)  : The col blocking scheme.                           !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine in April 2011.                                    !
  ! Edited for multiple parallel strategies by Robert Charlton, 10/12/2016.    !
  !============================================================================!

  real(kind=DP) function sparse_fill_fac_denom(row_blks,col_blks,row_par,col_par)

    use parallel_strategy, only: PARAL_INFO

    implicit none

    ! Arguments
    integer, intent(in) :: row_blks
    integer, intent(in) :: col_blks
    type(PARAL_INFO), intent(in), pointer :: row_par
    type(PARAL_INFO), intent(in), pointer :: col_par

    ! Local vaiables
    real(kind=DP) :: denom

    ! rc2013: allow for zeros in denominator
    denom = &
         real(sum(row_par%num_elems_on_proc(:,row_blks)),kind=DP) * &
         real(sum(col_par%num_elems_on_proc(:,col_blks)),kind=DP)
    if (denom>0.0_DP) then
       sparse_fill_fac_denom = 100.0_DP / (&
            real(sum(row_par%num_elems_on_proc(:,row_blks)),kind=DP) * &
            real(sum(col_par%num_elems_on_proc(:,col_blks)),kind=DP))
    else
       sparse_fill_fac_denom = 0.0_DP
    endif

  end function sparse_fill_fac_denom


  !============================================================================!
  ! This function returns true if every element in the matrix is nonzero, or   !
  ! false otherwise.                                                           !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   mat   (input)  : The matrix to be assessed                               !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine in December 2009.                                 !
  !============================================================================!

  logical function sparse_is_dense(mat)

    implicit none

    ! Arguments
    type(SPAM3), intent(in) :: mat   ! The sparse matrix

    ! Local Variables
    integer(kind=LONG) :: num_sq

    ! Check if the number of nonzero elements is equal to nrows*mcols
    num_sq = int(library(mat%lib)%nrows,kind=LONG) * &
         int(library(mat%lib)%mcols,kind=LONG)

    ! Return true if so, false otherwise
    if (library(mat%lib)%nze==num_sq) then
       sparse_is_dense = .true.
    else
       sparse_is_dense = .false.
    end if

  end function sparse_is_dense


  !============================================================================!
  ! This function calculates the RMS element value of a sparse matrix.         !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   mat   (input)  : The matrix to be assessed                               !
  !----------------------------------------------------------------------------!
  ! Written by Peter Haynes, June 2006.                                        !
  !============================================================================!

  real(kind=DP) function sparse_rms_element(mat)

    use comms, only: comms_reduce

    implicit none

    ! Arguments
    type(SPAM3), intent(in) :: mat   ! The sparse matrix

    ! Local variables
    integer :: ilib            ! Library entry for mat
    integer :: iel             ! Element loop counter
    real(kind=DP) :: loc_sum   ! Sum accumulated on this processor
    real(kind=DP) :: glob_sum  ! Sum accumulated across all processors

    ! Get library entry for mat
    ilib = mat%lib

    ! Sum modulus squared of each element on this processor
    loc_sum = 0.0_DP
    if (mat%iscmplx) then
       do iel=1,library(ilib)%my_nze
          loc_sum = loc_sum + real(mat%zmtx(iel)*conjg(mat%zmtx(iel)),kind=DP)
       end do
    else
       do iel=1,library(ilib)%my_nze
          loc_sum = loc_sum + mat%dmtx(iel)*mat%dmtx(iel)
       end do
    end if

    ! Sum up across processors
    glob_sum = loc_sum
    call comms_reduce('SUM',glob_sum)

    ! Calculate RMS value
    sparse_rms_element = sqrt(glob_sum / real(library(ilib)%nze,kind=DP))

  end function sparse_rms_element


  !============================================================================!
  ! This function calculates the maximum absolute element value of a sparse    !
  ! matrix.                                                                    !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   mat   (input)  : The matrix to be assessed                               !
  !----------------------------------------------------------------------------!
  ! Written by Peter Haynes, June 2006.                                        !
  !============================================================================!

  real(kind=DP) function sparse_max_abs_element(mat)

    use comms, only: comms_reduce

    implicit none

    ! Arguments
    type(SPAM3), intent(in) :: mat   ! The sparse matrix

    ! Local variables
    integer :: ilib            ! Library entry for mat
    integer :: iel             ! Element loop counter
    real(kind=DP) :: loc_max   ! Maximum value on this processor
    real(kind=DP) :: glob_max  ! Maximum value across all processors

    ! Get library entry for mat
    ilib = mat%lib

    ! Find maximum absolute value of every element on this processor
    loc_max = 0.0_DP
    if (mat%iscmplx) then
       do iel=1,library(ilib)%my_nze
          loc_max = max(loc_max, &
               real(mat%zmtx(iel) * conjg(mat%zmtx(iel)),kind=DP))
       end do
       loc_max = sqrt(loc_max)
    else
       do iel=1,library(ilib)%my_nze
          loc_max = max(loc_max, abs(mat%dmtx(iel)))
       end do
    end if

    !  across processors
    glob_max = loc_max
    call comms_reduce('MAX',glob_max)

    ! Calculate RMS value
    sparse_max_abs_element = glob_max

  end function sparse_max_abs_element


  !============================================================================!
  ! This function returns the number of non-zero elements in a sparse matrix   !
  ! in the columns stored on this proc.                                        !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   mat   (input)  : The matrix to be assessed                               !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, June 2009.                                       !
  !============================================================================!

  integer function sparse_proc_num_element(mat)

    implicit none

    ! Arguments
    type(SPAM3), intent(in) :: mat   ! The sparse matrix

    ! Very simple...
    sparse_proc_num_element = library(mat%lib)%my_nze

  end function sparse_proc_num_element


  !============================================================================!
  ! This function returns the number of non-zero elements in a sparse matrix   !
  ! specified by its library entry, in a floating-point real.                  !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   mat   (input)  : The matrix to be assessed                               !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, April 2011.                                      !
  !============================================================================!

  real(kind=DP) function sparse_num_element_lib(ilib)

    implicit none

    ! Arguments
    integer, intent(in) :: ilib   ! The library entry

    ! Very simple...
    sparse_num_element_lib = real(library(ilib)%nze,kind=DP)

  end function sparse_num_element_lib


  !============================================================================!
  ! This function returns the number of non-zero elements in a sparse matrix   !
  ! in a floating-point real (so as to avoid integer overflows).               !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   mat   (input)  : The matrix to be assessed                               !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, December 2009, based on sparse_num_element,      !
  ! originally written by Peter Haynes, June 2006.                             !
  !============================================================================!

  real(kind=DP) function sparse_num_element_mat(mat)

    implicit none

    ! Arguments
    type(SPAM3), intent(in) :: mat   ! The sparse matrix

    ! Very simple...
    sparse_num_element_mat = real(library(mat%lib)%nze,kind=DP)

  end function sparse_num_element_mat


  !============================================================================!
  ! This function returns the number of rows in a sparse matrix by consulting  !
  ! the relevant library entry.                                                !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   mat   (input)  : The matrix to be assessed                               !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, January 2010.                                    !
  ! Result type changed from real to integer by JM Escartin, Summer 2016.      !
  !============================================================================!

  integer function sparse_num_rows(mat)

    implicit none

    ! Arguments
    type(SPAM3), intent(in) :: mat   ! The sparse matrix

    ! Very simple...
    sparse_num_rows = library(mat%lib)%nrows

  end function sparse_num_rows


  !============================================================================!
  ! This function returns the number of cols in a sparse matrix by consulting  !
  ! the relevant library entry.                                                !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   mat   (input)  : The matrix to be assessed                               !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, January 2010.                                    !
  ! Result type changed from real to integer by JM Escartin, Summer 2016.      !
  !============================================================================!

  integer function sparse_num_cols(mat)

    implicit none

    ! Arguments
    type(SPAM3), intent(in) :: mat   ! The sparse matrix

    ! Very simple...
    sparse_num_cols = library(mat%lib)%mcols

  end function sparse_num_cols


  !============================================================================!
  ! This function returns the number of non-zero elements in a sparse matrix   !
  ! in the columns stored on this proc.                                        !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   mat   (input)  : The matrix to be assessed                               !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, February 2010.                                   !
  !============================================================================!

  subroutine sparse_memory(local_mem,global_mem,structure)

    use utils, only: utils_assert

    implicit none

    ! Arguments
    integer, intent(out) :: local_mem
    integer, intent(out) :: global_mem
    character(*), intent(in) :: structure   ! Sparse matrix structure code

    ! Local Variables
    integer :: ilib

    call sparse_search_library(structure,ilib)

    ! Quit if not found
    call utils_assert(ilib /= 0, 'Error in sparse_memory: &
            &no library structure found to match '//trim(structure))

    local_mem = library(ilib)%my_nze
    global_mem = int(library(ilib)%nze)

  end subroutine sparse_memory


  !============================================================================!
  ! This function displays the currently allocated memory of all the sparse    !
  ! matrices currently in existence.                                           !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, February 2010.                                   !
  !============================================================================!

  subroutine sparse_show_memory_usage

    use comms, only: pub_on_root

    implicit none

    if (pub_on_root) then
       write(stdout,'(a)') 'Currently allocated Sparse Matrix Memory'
       write(stdout,'(a,i20)') 'Global           :',global_mat_mem
       write(stdout,'(a,i20)') 'Local (root proc):',local_mat_mem
    end if

  end subroutine sparse_show_memory_usage

  !============================================================================!
  ! This subroutine displays the elements of the matrix smat for debugging     !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   smat  (inout) : The matrix to print                                      !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, May 2009.                                        !
  !============================================================================!

  subroutine sparse_show_matrix(smat,outunit,show_elems,matlab_format,labels)

    use comms, only: comms_barrier, comms_bcast, comms_send, &
         comms_recv, pub_on_root, pub_my_proc_id, pub_root_proc_id, &
         pub_total_num_procs, comms_free
    use parallel_strategy, only: PARAL_INFO
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_flush, &
         utils_assert

    implicit none

    ! Arguments
    type(SPAM3),intent(in) :: smat
    integer, intent(in), optional :: outunit
    logical, intent(in), optional :: show_elems
    logical, intent(in), optional :: matlab_format

    ! Locals
    integer :: loc_outunit        !
    logical :: loc_show_elems     !
    logical :: loc_matlab_format  !
    integer :: ierr               !
    integer :: seg                !
    integer :: proc               !
    integer :: nrows              ! Total number of rows in this block/segment
    integer :: mcols              !
    integer :: datlen             !
    integer :: ptr                !
    integer :: iorig              !
    integer :: row_start, row_end !
    integer :: row_len            !
    integer :: col_start, col_end !
    integer :: col_len            !
    integer :: col_blks           ! Identifier for col blocking scheme
    integer :: row_blks           ! Identifier for row blocking scheme
    integer :: ielem,global_ielem !
    integer :: iatom              !
    integer :: loc_iatom          !
    integer :: max_seg_rows       !
    integer :: predpwidth         ! Pre decimal point width
    integer :: width              ! Width of number for display
    real(kind=DP), allocatable :: dmtxbuf(:)    ! Buffer for real data
    complex(kind=DP), allocatable :: zmtxbuf(:) ! Buffer for complex data
    real(kind=DP), allocatable :: dmtxout(:)    !
    complex(kind=DP), allocatable :: zmtxout(:) !
    character(len=1000) :: fmt                  !
    character(len=12) :: tmp

    character(len=2),allocatable :: row_symbols(:) !
    character(len=2),allocatable :: col_symbols(:) !
    integer,allocatable :: row_elem_on_atom(:)     !
    integer,allocatable :: col_elem_on_atom(:)     !
    character(len=2),optional :: labels(library(smat%lib)%mcols)
    ! rc2013: local variables for row and column parallel strategies
    type(PARAL_INFO), pointer :: row_par
    type(PARAL_INFO), pointer :: col_par

    ! Process optional arguments
    if (present(outunit)) then
       loc_outunit = outunit
    else
       loc_outunit = 6
    end if
    if (present(show_elems)) then
       loc_show_elems = show_elems
    else
       loc_show_elems = .false.
    end if
    if (present(matlab_format)) then
       loc_matlab_format = matlab_format
    else
       loc_matlab_format = .false.
    end if

    ! Consistency Check
    call utils_assert(.not. loc_matlab_format.or. .not.loc_show_elems, &
         'Consistency check failed in sparse_show_matrix.')

    ! rc2013: assign the pointers for par based on library structure
    row_par => library(smat%lib)%row_par
    col_par => library(smat%lib)%col_par

    mcols = library(smat%lib)%mcols
    nrows = library(smat%lib)%nrows
    col_blks = library(smat%lib)%col_blks
    row_blks = library(smat%lib)%row_blks
    max_seg_rows = maxval(row_par%num_elems_on_proc(:,row_blks))

    ! rc2013: assign the pointers for par based on library structure
    row_par => library(smat%lib)%row_par
    col_par => library(smat%lib)%col_par

    ! Allocate communication buffer
    if (smat%iscmplx) then
       allocate(zmtxbuf(mcols*max_seg_rows),stat=ierr)
       call utils_alloc_check('sparse_show_matrix','zmtxbuf',ierr)
       allocate(zmtxout(mcols),stat=ierr)
       call utils_alloc_check('sparse_show_matrix','zmtxout',ierr)
    else
       allocate(dmtxbuf(mcols*max_seg_rows),stat=ierr)
       call utils_alloc_check('sparse_show_matrix','dmtxbuf',ierr)
       allocate(dmtxout(mcols),stat=ierr)
       call utils_alloc_check('sparse_show_matrix','dmtxout',ierr)
    end if

    ! Allocate storage for symbols and elem_on_atom lists
    if (loc_show_elems) then

       allocate(col_symbols(mcols),stat=ierr)
       call utils_alloc_check('sparse_show_matrix','col_symbols',ierr)
       allocate(row_symbols(nrows),stat=ierr)
       call utils_alloc_check('sparse_show_matrix','row_symbols',ierr)
       allocate(col_elem_on_atom(mcols),stat=ierr)
       call utils_alloc_check('sparse_show_matrix','col_elem_on_atom',ierr)
       allocate(row_elem_on_atom(nrows),stat=ierr)
       call utils_alloc_check('sparse_show_matrix','row_elem_on_atom',ierr)

       ! Collect this proc's symbol strings and element indices along cols
       ! jcap: if this proc is relevant
       if ( (pub_my_proc_id.ge.col_par%first_proc) .and. &
            (pub_my_proc_id.lt.col_par%first_proc+col_par%num_procs) ) then
          do ielem=col_par%first_elem_on_proc(pub_my_proc_id,col_blks), &
               col_par%first_elem_on_proc(pub_my_proc_id+1,col_blks)-1
             iatom = col_par%atom_of_elem(ielem,col_blks)
             loc_iatom = iatom - col_par%first_atom_on_proc(pub_my_proc_id) + 1
             col_symbols(ielem) = col_par%elements_on_proc(loc_iatom)%symbol
             col_elem_on_atom(ielem) = ielem - col_par%first_elem_on_atom(iatom,col_blks)+1
          end do
       end if

       ! Collect this proc's symbol strings and element indices along rows
       ! jcap: if this proc is relevant
       if ( (pub_my_proc_id.ge.row_par%first_proc) .and. &
            (pub_my_proc_id.lt.row_par%first_proc+row_par%num_procs) ) then
          do ielem=row_par%first_elem_on_proc(pub_my_proc_id,row_blks), &
               row_par%first_elem_on_proc(pub_my_proc_id+1,row_blks)-1
             iatom = row_par%atom_of_elem(ielem,row_blks)
             loc_iatom = iatom - row_par%first_atom_on_proc(pub_my_proc_id) + 1
             row_symbols(ielem) = row_par%elements_on_proc(loc_iatom)%symbol
             row_elem_on_atom(ielem) = ielem - row_par%first_elem_on_atom(iatom,row_blks)+1
          end do
       end if

       ! Share symbol strings and element indices
       do proc=col_par%first_proc,col_par%first_proc+col_par%num_procs-1
          col_start = col_par%first_elem_on_proc(proc,col_blks)
          col_end = col_par%first_elem_on_proc(proc+1,col_blks) - 1
          col_len = col_end - col_start + 1
          call comms_bcast(proc,col_symbols(col_start:col_end),col_len*2)
          call comms_bcast(proc,col_elem_on_atom(col_start:col_end),col_len)
       end do

       do proc=row_par%first_proc,row_par%first_proc+row_par%num_procs-1
          row_start = row_par%first_elem_on_proc(proc,row_blks)
          row_end = row_par%first_elem_on_proc(proc+1,row_blks) - 1
          row_len = row_end - row_start + 1
          call comms_bcast(proc,row_symbols(row_start:row_end),row_len*2)
          call comms_bcast(proc,row_elem_on_atom(row_start:row_end),row_len)
       end do
       ! CW
       if(present(labels)) then
          labels=col_symbols
       endif
    end if

    ! Only the root proc shows the matrix
    if (pub_on_root.and.loc_show_elems) then
       ! Write the line identifying the atoms and elements on each atom
       write(loc_outunit,"(12x)",advance='no')
       do ielem=1,mcols
          iorig = col_par%orig_atom(col_par%atom_of_elem(ielem,col_blks))
          if (smat%iscmplx) then
             write(loc_outunit,"(15x,a2,i4,2x,i2,1x)",advance='no') &
                  trim(col_symbols(ielem)),iorig,col_elem_on_atom(ielem)
          else
             write(loc_outunit,"(7x,a2,i4,2x,i2,1x)",advance='no') &
                  trim(col_symbols(ielem)),iorig,col_elem_on_atom(ielem)
          end if
       end do
       write(loc_outunit,*)
    end if

    ! qoh: Find width neccessary to show widest number in matrix
    ! dhpt: fix for case of elements exactly zero
    if (sparse_max_abs_element(smat)==0.0_DP) then
       predpwidth = 1
    else
       predpwidth = ceiling(log10(sparse_max_abs_element(smat)))
    endif
    predpwidth = max(predpwidth,1)

    if (pub_on_root) then
       ! Write the format string to use for the main data
       write(tmp,'(i10)') mcols
       if (smat%iscmplx) then
          if (loc_matlab_format) then
             write(fmt,"(a,a,a)")'(sp,',trim(adjustl(tmp)), &
                  "(e16.10,e16.10,'i',x),';')"
          else if (loc_show_elems) then
             write(fmt,"(a,a,a)")'(a2,i4,i4,2x,sp,',trim(adjustl(tmp)), &
                  "(f12.8,f12.8,'i',x))"
          else
             write(fmt,"(a,a,a)")'(sp,',trim(adjustl(tmp)), &
                  "(f20.14,f20.14,'i',x))"
          end if
       else
          if (loc_matlab_format) then
             write(fmt,"(a,a,a)")'(',trim(adjustl(tmp)),"f20.12,';')"
          else if (loc_show_elems) then
             write(fmt,"(a,a,a)")'(a2,i4,i4,2x,',trim(adjustl(tmp)),'f18.10)'
          else
             width = 17 + predpwidth
             write(fmt,"(a,a,a,i3,a)")'(',trim(adjustl(tmp)),'f',width,'.14)'
          end if
       end if

    end if

    ! Loop over segments
    if (pub_on_root.and.loc_matlab_format) &
         write(loc_outunit,'(a)',advance='no') '['
    do seg=0,pub_total_num_procs-1
       call comms_free

       ! Re-create the dense version of this segment if required
       if (smat%iscmplx) then
          call internal_dense_segment_complex(smat,seg,zmtxbuf)
       else
          call internal_dense_segment_real(smat,seg,dmtxbuf)
       end if

       nrows = row_par%num_elems_on_proc(seg,row_blks)
       datlen = nrows * col_par%num_elems_on_proc(pub_my_proc_id,col_blks)
       ptr = 1
       do proc=0,pub_total_num_procs-1
          datlen = nrows * col_par%num_elems_on_proc(proc,col_blks)
          if (pub_on_root.and.(proc/=pub_my_proc_id)) then
             if (smat%iscmplx) then
                call comms_recv(proc,zmtxbuf(ptr:ptr+datlen-1),datlen)
             else
                call comms_recv(proc,dmtxbuf(ptr:ptr+datlen-1),datlen)
             end if
          else if (.not.(pub_on_root).and.(proc==pub_my_proc_id)) then
             if (smat%iscmplx) then
                call comms_send(pub_root_proc_id,zmtxbuf(1:datlen),datlen)
             else
                call comms_send(pub_root_proc_id,dmtxbuf(1:datlen),datlen)
             end if
          end if
          ptr = ptr + datlen
       end do

       if (pub_on_root) then

          ! Loop over rows in this segment
          do ielem=1,row_par%num_elems_on_proc(seg,row_blks)
             global_ielem = ielem + row_par%first_elem_on_proc(seg,row_blks) - 1
             iorig = row_par%orig_atom(row_par%atom_of_elem(global_ielem,row_blks))
             if (smat%iscmplx) then
                ptr = ielem
                do proc=0,pub_total_num_procs-1
                   datlen = nrows * col_par%num_elems_on_proc(proc,col_blks)
                   zmtxout(col_par%first_elem_on_proc(proc,col_blks): &
                        col_par%first_elem_on_proc(proc+1,col_blks)-1) = &
                        zmtxbuf(ptr:ptr+datlen-ielem:nrows)
                   ptr = ptr + datlen
                end do
                if (loc_show_elems) then
                   write(loc_outunit,fmt) trim(row_symbols(global_ielem)), &
                        iorig,row_elem_on_atom(global_ielem),zmtxout(1:mcols)
                else
                   write(loc_outunit,fmt) zmtxout(1:mcols)
                end if
             else
                ptr = ielem
                do proc=0,pub_total_num_procs-1
                   datlen = nrows * col_par%num_elems_on_proc(proc,col_blks)
                   if (datlen == 0) cycle
                   dmtxout(col_par%first_elem_on_proc(proc,col_blks): &
                        col_par%first_elem_on_proc(proc+1,col_blks)-1) = &
                        dmtxbuf(ptr:ptr+datlen-ielem:nrows)
                   ptr = ptr + datlen
                end do
                if (loc_show_elems) then
                   write(loc_outunit,fmt) trim(row_symbols(global_ielem)), &
                        iorig,row_elem_on_atom(global_ielem),dmtxout(1:mcols)
                else
                   write(loc_outunit,fmt) dmtxout(1:mcols)
                end if
             endif
          end do

       end if

       call comms_barrier

    end do
    if (pub_on_root.and.loc_matlab_format) &
         write(loc_outunit,'(a)') ']'

    call utils_flush
    call comms_free
    call comms_barrier
    call utils_flush

    ! Dellocate storage for symbols and elem_on_atom list
    if (loc_show_elems) then
       deallocate(row_elem_on_atom,stat=ierr)
       call utils_dealloc_check('sparse_show_matrix','row_elem_on_atom',ierr)
       deallocate(col_elem_on_atom,stat=ierr)
       call utils_dealloc_check('sparse_show_matrix','col_elem_on_atom',ierr)
       deallocate(row_symbols,stat=ierr)
       call utils_dealloc_check('sparse_show_matrix','row_symbols',ierr)
       deallocate(col_symbols,stat=ierr)
       call utils_dealloc_check('sparse_show_matrix','col_symbols',ierr)
    end if

    ! Deallocate communication buffers
    if (smat%iscmplx) then
       deallocate(zmtxout,stat=ierr)
       call utils_dealloc_check('sparse_show_matrix','zmtxout',ierr)
       deallocate(zmtxbuf,stat=ierr)
       call utils_dealloc_check('sparse_show_matrix','zmtxbuf',ierr)
    else
       deallocate(dmtxout,stat=ierr)
       call utils_dealloc_check('sparse_show_matrix','dmtxout',ierr)
       deallocate(dmtxbuf,stat=ierr)
       call utils_dealloc_check('sparse_show_matrix','dmtxbuf',ierr)
    end if

contains

    !========================================================================!
    ! Collates the data of a real-valued segment of a matrix (sparse or      !
    ! dense) so that it can be printed.                                      !
    !------------------------------------------------------------------------!
    ! Written by Nicholas Hine, May 2009.                                    !
    !========================================================================!

    subroutine internal_dense_segment_real(mat,seg,dbuf)

      implicit none

      ! Arguments
      type(SPAM3),intent(in) :: mat
      integer,intent(in) :: seg
      real(kind=DP), intent(out) :: dbuf(:)

      ! Locals
      integer :: ilib          ! Library entry for mat
      integer :: iblk          ! Block-column counter
      integer :: loc_iblk      ! Block-column counter local to this proc
      integer :: jblk          ! Block-row counter
      integer :: iovlap        ! Block overlap counter
      integer :: ielem,jelem   ! Element row/col counters
      integer :: ptr           ! Pointer to data in block
      integer :: datlen        ! Length of data in block
      integer :: bufptr        ! Pointer to data in buffer
      integer :: row_blks      ! Identifier for column blocking scheme
      integer :: col_blks      ! Identifier for row blocking scheme
      integer :: seg_type      ! Segment type for this segment
      integer :: seg_start     ! Start position of this segment in the index

      ! Get library entry
      ilib = mat%lib
      row_blks = library(ilib)%row_blks
      col_blks = library(ilib)%col_blks

      seg_start = library(ilib)%seg_info(s_idx,seg)
      seg_type = library(ilib)%seg_info(s_type,seg)

      datlen = col_par%num_elems_on_proc(pub_my_proc_id,col_blks) * &
           row_par%num_elems_on_proc(seg,row_blks)
      dbuf(1:datlen) = 0.0_DP

      if (seg_type == SEG_DENSE) then

         ptr = library(ilib)%seg_info(s_ptr,seg)
         datlen = library(ilib)%seg_info(s_ptr,seg + 1) - ptr
         dbuf(1:datlen) = mat%dmtx(ptr:ptr+datlen-1)

      else if (seg_type == SEG_SPARSE) then

         ! Loop over atom block-columns iblk on this proc
         loc_iblk = 0
         seg_start = library(ilib)%seg_info(s_idx,seg)
         do iblk=col_par%my_first_blk,col_par%my_last_blk
            loc_iblk = loc_iblk + 1

            ! Loop over atoms jblk overlapping atom iblk: block-rows in
            ! block-col iblk of spam
            do iovlap=library(ilib)%blk_idx(seg_start+loc_iblk-1), &
                 library(ilib)%blk_idx(seg_start+loc_iblk)-1
               jblk = library(ilib)%blk_idx(iovlap)
               ptr = library(ilib)%blk_ptr(iovlap)

               bufptr = (col_par%first_elem_on_atom(iblk,col_blks) - &
                    col_par%first_elem_on_proc(pub_my_proc_id,col_blks))* &
                    row_par%num_elems_on_proc(seg,row_blks) + &
                    row_par%first_elem_on_atom(jblk,row_blks) - &
                    row_par%first_elem_on_proc(seg,row_blks) + 1
               do ielem=1,col_par%num_elems_on_atom(iblk,col_blks)
                  do jelem=1,row_par%num_elems_on_atom(jblk,row_blks)
                     dbuf(bufptr) = mat%dmtx(ptr)
                     bufptr = bufptr + 1
                     ptr = ptr + 1
                  end do ! jelem
                  bufptr = bufptr - row_par%num_elems_on_atom(jblk,row_blks)
                  bufptr =  bufptr + row_par%num_elems_on_proc(seg,row_blks)
               end do ! ielem

            end do ! iovlap

         end do ! iblk

      end if ! seg_type

    end subroutine internal_dense_segment_real

    !========================================================================!
    ! Collates the data of a complex-valued segment of a matrix (sparse or   !
    ! dense) so that it can be printed.                                      !
    !------------------------------------------------------------------------!
    ! Written by Nicholas Hine, May 2009.                                    !
    !========================================================================!

    subroutine internal_dense_segment_complex(mat,seg,zbuf)

      implicit none

      ! Arguments
      type(SPAM3),intent(in) :: mat
      integer,intent(in) :: seg
      complex(kind=DP), intent(out) :: zbuf(:)

      ! Locals
      integer :: ilib          ! Library entry for mat
      integer :: iblk          ! Block-column counter
      integer :: loc_iblk      ! Block-column counter local to this proc
      integer :: jblk          ! Block-row counter
      integer :: iovlap        ! Block overlap counter
      integer :: ielem,jelem   ! Element row/col counters
      integer :: ptr           ! Pointer to data in block
      integer :: datlen        ! Length of data in block
      integer :: bufptr        ! Pointer to data in buffer
      integer :: row_blks      ! Identifier for column blocking scheme
      integer :: col_blks      ! Identifier for row blocking scheme
      integer :: seg_type      ! Segment type for this segment
      integer :: seg_start     ! Start position of this segment in the index

      ! Get library entry
      ilib = mat%lib
      row_blks = library(ilib)%row_blks
      col_blks = library(ilib)%col_blks

      seg_start = library(ilib)%seg_info(s_idx,seg)
      seg_type = library(ilib)%seg_info(s_type,seg)

      datlen = col_par%num_elems_on_proc(pub_my_proc_id,col_blks) * &
           row_par%num_elems_on_proc(seg,row_blks)
      zbuf(1:datlen) = (0.0_DP,0.0_DP)

      if (seg_type == SEG_DENSE) then

         ptr = library(ilib)%seg_info(s_ptr,seg)
         datlen = library(ilib)%seg_info(s_ptr,seg + 1) - ptr
         zbuf(1:datlen) = mat%zmtx(ptr:ptr+datlen-1)

      else if (seg_type == SEG_SPARSE) then

         ! Loop over atom block-columns iblk on this proc
         loc_iblk = 0
         seg_start = library(ilib)%seg_info(s_idx,seg)
         do iblk=col_par%my_first_blk,col_par%my_last_blk
            loc_iblk = loc_iblk + 1

            ! Loop over atoms jblk overlapping atom iblk: block-rows in
            ! block-col iblk of spam
            do iovlap=library(ilib)%blk_idx(seg_start+loc_iblk-1), &
                 library(ilib)%blk_idx(seg_start+loc_iblk)-1
               jblk = library(ilib)%blk_idx(iovlap)
               ptr = library(ilib)%blk_ptr(iovlap)

               bufptr = (col_par%first_elem_on_atom(iblk,col_blks) - &
                    col_par%first_elem_on_proc(pub_my_proc_id,col_blks))* &
                    row_par%num_elems_on_proc(seg,row_blks) + &
                    row_par%first_elem_on_atom(jblk,row_blks) - &
                    row_par%first_elem_on_proc(seg,row_blks) + 1
               do ielem=1,col_par%num_elems_on_atom(iblk,col_blks)
                  do jelem=1,row_par%num_elems_on_atom(jblk,row_blks)
                     zbuf(bufptr) = mat%zmtx(ptr)
                     bufptr = bufptr + 1
                     ptr = ptr + 1
                  end do ! jelem
                  bufptr = bufptr - row_par%num_elems_on_atom(jblk,row_blks)
                  bufptr =  bufptr + row_par%num_elems_on_proc(seg,row_blks)
               end do ! ielem

            end do ! iovlap

         end do ! iblk

      end if ! seg_type

    end subroutine internal_dense_segment_complex

  end subroutine sparse_show_matrix


  subroutine sparse_show_network(mat,iunit)

    use comms, only: comms_free, pub_on_root, pub_my_proc_id, &
         pub_total_num_procs, pub_root_proc_id, comms_barrier
    use parallel_strategy, only: PARAL_INFO

    implicit none

    ! Arguments
    type(SPAM3), intent(inout) :: mat         ! The matrix
    integer, intent(in) :: iunit

    ! Local Variables
    type(COM3) :: com
    integer :: recvproc
    integer :: ilib
    integer :: iblk,loc_iblk
    integer :: jblk,loc_jblk
    integer :: kblk
    integer :: jidx
    integer :: seg, seg_start
    integer :: nrow,ncol
    logical :: found
    type(PARAL_INFO), pointer :: row_par ! parallel strategy for rows
    type(PARAL_INFO), pointer :: col_par ! parallel strategy for columns

    ! Get library entry
    ilib = mat%lib
    ! rc2013: assign the pointers for par based on library structure
    row_par => library(ilib)%row_par
    col_par => library(ilib)%col_par

    ! Allocate arrays and initialise comms
    call sparse_com_allocate(com,mat,1, &
         alloc_mtx=.false.,cropped=.false.,seg=.false.,groupseg=.false.)

    do recvproc=0,pub_total_num_procs-1

       call comms_barrier

       ! Check if this proc needs to send yet
       if ((recvproc==pub_my_proc_id).and.(.not.pub_on_root)) then
          call sparse_send_index(com,pub_root_proc_id)
       end if

       ! Root proc receives and writes
       if (pub_on_root) then
          call sparse_recv_index(com,recvproc,1,async=.false.)

          do iblk=row_par%first_atom_on_proc(recvproc), &
               row_par%first_atom_on_proc(recvproc+1)-1
             loc_iblk = iblk - row_par%first_atom_on_proc(recvproc) + 1
             ncol = row_par%num_elems_on_atom(iblk,library(mat%lib)%row_blks)

             do seg=0,pub_total_num_procs-1
                seg_start = com%seginfobuf(s_idx,seg,1)
                if (com%seginfobuf(s_type,seg,1)==SEG_BLANK) then
                   do loc_jblk=1,row_par%num_atoms_on_proc(seg)
                      write(iunit,'(i4)',advance='no') 0
                   end do
                   cycle
                end if
                do loc_jblk=1,row_par%num_atoms_on_proc(seg)
                   jblk = loc_jblk + row_par%first_atom_on_proc(seg) -1
                   found = .false.
                   ! Loop over block-rows kblk of A in column jblk
                   do jidx=com%idxbuf(seg_start+loc_iblk-1,1), &
                        com%idxbuf(seg_start+loc_iblk,1)-1
                      kblk = com%idxbuf(jidx,1)
                      if (kblk==jblk) then
                         found=.true.
                         exit
                      end if
                   end do
                   if (found) then
                      nrow = row_par%num_elems_on_atom(jblk,library(mat%lib)%row_blks)
                      write(iunit,'(i4)',advance='no') nrow*ncol
                   else
                      write(iunit,'(i4)',advance='no') 0
                   end if
                end do
             end do
             write(iunit,*)
          end do
       end if
    end do

    ! Deallocate arrays and exit comms
    call comms_free
    call comms_barrier
    call sparse_com_deallocate(com,dealloc_mtx=.false.,cropped=.false.)

  end subroutine sparse_show_network

  !============================================================================!
  ! This subroutine inserts an element into a given block sparse matrix local  !
  ! to this processor.                                                         !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   el (input)     : The element of the real matrix to insert                !
  !   mat (inout)    : The real matrix                                         !
  !   jrow (input)   : The element index of the row required                   !
  !   icol (input)   : The element index of the column required                !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, May 2009 based on SPAM2 version.                 !
  !----------------------------------------------------------------------------!
  ! SPAM2 version written by Peter Haynes, in June 2006.                       !
  !============================================================================!

  subroutine sparse_put_element_real(el,mat,jrow,icol)

    use comms, only: pub_my_proc_id
    use parallel_strategy, only: PARAL_INFO
    use rundat, only: pub_debug
    use utils, only: utils_assert

    implicit none

    ! Arguments
    real(kind=DP), intent(in) :: el           ! The element to insert
    type(SPAM3), intent(inout) :: mat         ! The matrix
    integer, intent(in) :: jrow               ! The index of the row
    integer, intent(in) :: icol               ! The index of the column

    ! Local variables
    integer :: ilib      ! Library entry for mat
    integer :: row_blks  ! Identifier for column blocking scheme
    integer :: col_blks  ! Identifier for row blocking scheme
    integer :: iblk      ! The block-row to which the element belongs
    integer :: jblk      ! The block-col to which the element belongs
    integer :: kblk      ! Block counter
    integer :: nrows     ! Total number of rows in this block/segment
    integer :: jelem0    ! The first element in the requested block-row
    integer :: ielem0    ! The first element in the requested block-column
    integer :: iovlap    ! Block overlap counter
    integer :: loc_iblk  ! Block-column counter local to this proc
    integer :: ptr       ! Pointer to data array
    integer :: seg       ! Segment index
    integer :: seg_type  ! Segment type for this segment
    integer :: seg_start ! Start position of this segment in the index
    type(PARAL_INFO), pointer :: row_par ! parallel strategy for rows
    type(PARAL_INFO), pointer :: col_par ! parallel strategy for columns

    ! Get library entry
    ilib = mat%lib
    row_blks = library(ilib)%row_blks
    col_blks = library(ilib)%col_blks
    nrows = 0

    ! rc2013: assign the pointers for par based on library structure
    row_par => library(ilib)%row_par
    col_par => library(ilib)%col_par

    ! Check arguments
    if(pub_debug) then
       call utils_assert(jrow >= 1 .and. jrow <= library(ilib)%nrows, &
            'Error in sparse_put_element_real: invalid row index')
       call utils_assert(icol >= 1 .and. icol <= library(ilib)%mcols, &
            'Error in sparse_put_element_real: invalid column index')
    end if

    ! Check this column is local to this processor
    call utils_assert (icol>=sparse_first_elem_on_proc(pub_my_proc_id,mat,'C').and.&
         icol<=sparse_last_elem_on_proc(pub_my_proc_id,mat,'C'), 'Error in sparse_&
         &put_element_real: requested column is not local to proc ', &
         pub_my_proc_id)

    ! Obtain information about the relevant block & segment for this element
    seg = row_par%proc_of_elem(jrow,row_blks)
    seg_type = library(ilib)%seg_info(s_type,seg)
    if (seg_type==SEG_BLANK) return
    seg_start = library(ilib)%seg_info(s_idx,seg)
    jblk = row_par%atom_of_elem(jrow,row_blks)
    iblk = col_par%atom_of_elem(icol,col_blks)
    jelem0 = row_par%first_elem_on_atom(jblk,row_blks)
    ielem0 = col_par%first_elem_on_atom(iblk,col_blks)
    if (seg_type==SEG_DENSE) then
       nrows = row_par%num_elems_on_proc(seg,row_blks)
    else if (seg_type==SEG_SPARSE) then
       nrows = row_par%num_elems_on_atom(jblk,row_blks)
    end if

    ! Find local index for the block-column
    loc_iblk = iblk - col_par%my_first_blk + 1

    do iovlap=library(ilib)%blk_idx(seg_start+loc_iblk-1), &
         library(ilib)%blk_idx(seg_start+loc_iblk)-1
       kblk = library(ilib)%blk_idx(iovlap)   ! block-row
       if (kblk < jblk) cycle
       if (kblk > jblk) exit
       ptr = library(ilib)%blk_ptr(iovlap)
       if (mat%iscmplx) then
          mat%zmtx(ptr + (icol - ielem0)*nrows + (jrow - jelem0)) = &
                     cmplx(el, 0.0_DP, kind=DP)
       else
          mat%dmtx(ptr + (icol - ielem0)*nrows + (jrow - jelem0)) = el
       end if
    end do

  end subroutine sparse_put_element_real

  !============================================================================!
  ! This subroutine retrieves an element into a given block sparse matrix      !
  ! local to this processor.                                                   !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   el (output)    : The element of the real matrix retrieved                !
  !   mat (inout)    : The real matrix                                         !
  !   jrow (input)   : The element index of the row required                   !
  !   icol (input)   : The element index of the column required                !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, May 2009 based on SPAM2 version.                 !
  !----------------------------------------------------------------------------!
  ! SPAM2 version written by Peter Haynes, in June 2006.                       !
  !============================================================================!

  subroutine sparse_get_element_real(el,mat,jrow,icol)

    use comms, only: pub_my_proc_id
    use parallel_strategy, only: PARAL_INFO
    use rundat, only: pub_debug
    use utils, only: utils_assert

    implicit none

    ! Arguments
    real(kind=DP), intent(out) :: el          ! The element to insert
    type(SPAM3), intent(in) :: mat            ! The matrix
    integer, intent(in) :: jrow               ! The index of the row
    integer, intent(in) :: icol               ! The index of the column

    ! Local variables
    integer :: ilib      ! Library entry for mat
    integer :: row_blks  ! Identifier for column blocking scheme
    integer :: col_blks  ! Identifier for row blocking scheme
    integer :: iblk      ! The block-row to which the element belongs
    integer :: jblk      ! The block-col to which the element belongs
    integer :: kblk      ! Block counter
    integer :: nrows     ! Total number of rows in this block/segment
    integer :: jelem0    ! The first element in the requested block-row
    integer :: ielem0    ! The first element in the requested block-column
    integer :: iovlap    ! Block overlap counter
    integer :: loc_iblk  ! Block-column counter local to this proc
    integer :: ptr       ! Pointer to data array
    integer :: seg       ! Segment index
    integer :: seg_type  ! Segment type for this segment
    integer :: seg_start ! Start position of this segment in the index
    type(PARAL_INFO), pointer :: row_par ! Parallel strategy for rows
    type(PARAL_INFO), pointer :: col_par ! Parallel strategy for columns

    ! Get library entry
    ilib = mat%lib
    row_blks = library(ilib)%row_blks
    col_blks = library(ilib)%col_blks

    ! rc2013: assign the pointers for par based on library structure
    row_par => library(ilib)%row_par
    col_par => library(ilib)%col_par

    ! Check arguments
    call utils_assert(.not. mat%iscmplx, 'Error in sparse_get_element_real: &
         &real matrices only.')

    if(pub_debug) then
       call utils_assert(jrow >= 1 .and. jrow <= library(ilib)%nrows, &
            'Error in sparse_get_element_real: invalid row index.')
       call utils_assert(icol >= 1 .and. icol <= library(ilib)%mcols, &
            'Error in sparse_get_element_real: invalid column index.')
    end if
    ! Check this column is local to this processor
    call utils_assert(icol>=sparse_first_elem_on_proc(pub_my_proc_id,mat,'C').and.&
         icol<=sparse_last_elem_on_proc(pub_my_proc_id,mat,'C'), 'Error in &
         &sparse_get_element_real: requested column is not local to proc ', &
         pub_my_proc_id)

    el = 0.0_DP
    nrows = 0

    ! Obtain information about the relevant block & segment for this element
    seg = row_par%proc_of_elem(jrow,row_blks)
    seg_type = library(ilib)%seg_info(s_type,seg)
    if (seg_type==SEG_BLANK) return
    seg_start = library(ilib)%seg_info(s_idx,seg)
    jblk = row_par%atom_of_elem(jrow,row_blks)
    iblk = col_par%atom_of_elem(icol,col_blks)
    jelem0 = row_par%first_elem_on_atom(jblk,row_blks)
    ielem0 = col_par%first_elem_on_atom(iblk,col_blks)
    if (seg_type==SEG_DENSE) then
       nrows = row_par%num_elems_on_proc(seg,row_blks)
    else if (seg_type==SEG_SPARSE) then
       nrows = row_par%num_elems_on_atom(jblk,row_blks)
    end if

    ! Find local index for the block-column
    loc_iblk = iblk - col_par%my_first_blk + 1

    do iovlap=library(ilib)%blk_idx(seg_start+loc_iblk-1), &
         library(ilib)%blk_idx(seg_start+loc_iblk)-1
       kblk = library(ilib)%blk_idx(iovlap)   ! block-row
       if (kblk < jblk) cycle
       if (kblk > jblk) exit
       ptr = library(ilib)%blk_ptr(iovlap)
       el = mat%dmtx(ptr + (icol - ielem0)*nrows + (jrow - jelem0))
    end do

  end subroutine sparse_get_element_real


  !============================================================================!
  ! This subroutine solves the linear equation BC = A for C using the          !
  ! algorithm of Gibson, Haydock and LaFemina.                                 !
  ! DOI:https://doi.org/10.1103/PhysRevB.47.9229 (see section "An Approximate  !
  ! Hopping Matrix"). This algorithm uses submatrices formed from the first-   !
  ! order connectivity in the sparsity pattern to form lower dimensioned       !
  ! linear problems which can be solved independently, in a serial fashion on  !
  ! multiple procs. The method should become competitive once the matrices are !
  ! large and sparse enough.                                                   !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  ! cmat (output)    : The solution matrix                                     !
  ! amat (input)     : The right hand side matrix                              !
  ! bmat (input)     : The left hand side matrix                               !
  ! aplusb_sparsity (optional) : Whether to use the union of A,B and C         !
  !                              sparsity pattern in computing the elements of !
  !                              C. Solves only for NZ blocks of C sparsity    !
  !                              pattern by default.                           !
  !----------------------------------------------------------------------------!
  ! Notes:                                                                     !
  ! This subspace solve routine is mainly useful when B is very sparse, for    !
  ! example when B is the overlap matrix and A is the Hamiltonin matrix, in    !
  ! which case C will be the contra-covariant Hamiltonian matrix known as the  !
  ! "hopping matrix" in LaFemina et al. If B is not very sparse, the routine   !
  ! become very inefficient and an alternative should be used.                 !
  !----------------------------------------------------------------------------!
  ! Written by Jolyon Aarons, May 2017.                                        !
  !============================================================================!
  subroutine sparse_solve2(cmat,amat,bmat,aplusb_sparsity)
    use comms, only: pub_my_proc_id, comms_reduce, pub_total_num_procs, comms_send, &
         comms_recv, comms_barrier, pub_on_root, comms_wait, pub_my_proc_id
    use constants, only: cmplx_0, cmplx_1
    use parallel_strategy, only: PARAL_INFO
    use timer, only: timer_clock, wrappers_etime
    use utils, only: utils_dealloc_check, utils_alloc_check, utils_assert, utils_unit
    implicit none

    type(SPAM3), intent(inout)    :: cmat
    type(SPAM3), intent(in)       :: amat
    type(SPAM3), intent(in)       :: bmat
    logical, optional, intent(in) :: aplusb_sparsity

    integer :: ierr
    integer :: loc_col
    integer :: num_col_blks
    integer :: ient
    integer :: ilib
    integer :: jcol,krow
    real(kind=dp), dimension(:,:), allocatable :: hopping_sub
    real(kind=dp), dimension(:,:), allocatable :: subspace_ham, subspace_olap
    real(kind=dp), dimension(:,:), allocatable :: subspace_buffer_ham, subspace_buffer_olap

    complex(kind=dp), dimension(:,:), allocatable :: zhopping_sub
    complex(kind=dp), dimension(:,:), allocatable :: zsubspace_ham, zsubspace_olap
    complex(kind=dp), dimension(:,:), allocatable :: zsubspace_buffer_ham, zsubspace_buffer_olap

    integer, dimension(:), allocatable :: pivot
    integer, dimension(:), allocatable :: entriesa, entriesb
    integer, dimension(:), allocatable :: recentries
    integer :: info
    type(PARAL_INFO), pointer :: col_par, row_par ! rc2013: parallel strategies
    integer :: colproc
    integer :: icol
    integer :: nblkcols
    integer :: colcount, rowcount
    integer :: colelems, rowelems
    integer :: nelems_on_proc, nentries_on_proc
    integer :: nrecentries
    integer :: first_entry_on_proc, firstcol_thatproc, firstcol_thisproc
    integer :: ncols_thatproc, ncols_thisproc
    integer :: iproc
    integer :: send_proc, recv_proc

    integer :: max_subspace
    type(SPAM3) :: tmpmat, tmpmat2

    integer :: send_handle
    logical :: add_to_stack

    integer :: col_blks, row_blks
    integer :: offset, offset2
    logical :: need_data, send_data

    logical :: c_sparsity

    integer :: outunit
    character(len=32) :: filename

    real(kind=dp) :: solve_time, start_time
    logical :: iscmplx

    ! LAPACK subroutines
    external :: zgesv, dgesv


    add_to_stack=.false.
    c_sparsity=.true.
    if(present(aplusb_sparsity)) then
       c_sparsity=.not.aplusb_sparsity
    end if

    ilib = amat%lib
    col_blks = library(ilib)%col_blks
    row_blks = library(ilib)%row_blks
    row_par => library(ilib)%row_par
    col_par => library(ilib)%col_par

    call utils_assert((amat%iscmplx.eqv.bmat%iscmplx).and. &
         (amat%iscmplx.eqv.cmat%iscmplx), &
         "Error in sparse_solve: A,B and C must all be either real or complex.")

    iscmplx=amat%iscmplx

    num_col_blks=col_par%my_last_blk
    if(pub_total_num_procs>1) then
       call comms_reduce('MAX',num_col_blks)
    end if

    if(allocated(subspace_ham)) then
       deallocate(subspace_ham,stat=ierr)
       call utils_dealloc_check('sparse_solve','subspace_ham',ierr)
    end if
    if(allocated(subspace_olap)) then
       deallocate(subspace_olap,stat=ierr)
       call utils_dealloc_check('sparse_solve','subspace_olap',ierr)
    end if
    if(allocated(hopping_sub)) then
       deallocate(hopping_sub,stat=ierr)
       call utils_dealloc_check('sparse_solve','hopping_sub',ierr)
    end if

    if(allocated(zsubspace_ham)) then
       deallocate(zsubspace_ham,stat=ierr)
       call utils_dealloc_check('sparse_solve','zsubspace_ham',ierr)
    end if
    if(allocated(zsubspace_olap)) then
       deallocate(zsubspace_olap,stat=ierr)
       call utils_dealloc_check('sparse_solve','zsubspace_olap',ierr)
    end if
    if(allocated(zhopping_sub)) then
       deallocate(zhopping_sub,stat=ierr)
       call utils_dealloc_check('sparse_solve','zhopping_sub',ierr)
    end if


    if(allocated(pivot)) then
       deallocate(pivot,stat=ierr)
       call utils_dealloc_check('sparse_solve','pivot',ierr)
    end if

    if(allocated(entriesa)) then
       deallocate(entriesa,stat=ierr)
       call utils_dealloc_check('sparse_solve','entriesa',ierr)
    end if
    if(allocated(entriesb)) then
       deallocate(entriesb,stat=ierr)
       call utils_dealloc_check('sparse_solve','entriesb',ierr)
    end if


    max_subspace=1
    do icol=col_par%my_first_blk,col_par%my_last_blk
       if(c_sparsity) then
          call get_entries(icol,cmat,entriesa)
       else
          call get_entries(icol,amat,entriesa,bmat,entriesb,.true.)
       end if
       max_subspace=max(max_subspace,size(entriesa))
       if(allocated(entriesa)) then
          deallocate(entriesa,stat=ierr)
          call utils_dealloc_check('sparse_solve','entriesa',ierr)
       end if

       if(allocated(entriesb)) then
          deallocate(entriesb,stat=ierr)
          call utils_dealloc_check('sparse_solve','entriesb',ierr)
       end if
    end do
    if(pub_total_num_procs>1) then
       call comms_reduce('MAX',max_subspace)
    end if

    ! Work out the max number of column blocks on a proc
    nblkcols=col_par%my_last_blk-col_par%my_first_blk+1
    if(pub_total_num_procs>1) then
       call comms_reduce('MAX',nblkcols)
    end if

    ! We loop over the max number of column blocks so that the procs which have
    ! finished their own data are still ready to meet send requests.x
    do icol=1,nblkcols
       loc_col=col_par%my_first_blk+icol-1

       need_data = (loc_col>=col_par%my_first_blk .and.loc_col<=col_par%my_last_blk)

       if(need_data) then
          if(c_sparsity) then
             call get_entries(loc_col,cmat,entriesa)
          else
             call get_entries(loc_col,amat,entriesa,bmat,entriesb,.true.)
          end if
          call sparse_sizeof_subspace_on_proc(amat,entriesa,rowelems,colelems, &
               ncols_thisproc,nelems_on_proc,first_entry_on_proc)

          if(iscmplx) then
             allocate(zsubspace_ham(rowelems,colelems),stat=ierr)
             call utils_alloc_check('sparse_solve','zsubspace_ham',ierr)

             allocate(zsubspace_olap(rowelems,colelems),stat=ierr)
             call utils_alloc_check('sparse_solve','zsubspace_olap',ierr)
          else
             allocate(subspace_ham(rowelems,colelems),stat=ierr)
             call utils_alloc_check('sparse_solve','subspace_ham',ierr)

             allocate(subspace_olap(rowelems,colelems),stat=ierr)
             call utils_alloc_check('sparse_solve','subspace_olap',ierr)
          end if
       end if

       if(pub_total_num_procs>1) then
          call comms_barrier()
       end if

       do iproc=1,pub_total_num_procs

          send_proc = mod(pub_my_proc_id+iproc,pub_total_num_procs)                      ! sending DATA to this proc
          recv_proc = mod(pub_total_num_procs+pub_my_proc_id-iproc,pub_total_num_procs)  ! recving DATA from this proc

          !when send_proc==recv_proc -> handle this proc!
          if(pub_total_num_procs==1) then
             send_data=need_data
          else
             call comms_send(recv_proc,need_data,return_handle=send_handle, &
                  add_to_stack=add_to_stack)
             call comms_recv(send_proc,send_data)
             call comms_wait(send_handle,add_to_stack)
          end if

          if(pub_total_num_procs>1) then
             call comms_barrier()
          end if

          if(need_data.and.(pub_total_num_procs>1)) then
             !send entries size to recv_proc
             call comms_send(recv_proc,size(entriesa), &
                  return_handle=send_handle,add_to_stack=add_to_stack)
          end if

          if(send_data) then
             if(pub_total_num_procs>1) then
                !recv entries size from send_proc
                call comms_recv(send_proc,nrecentries)
             else
                nrecentries=size(entriesa)
             end if

             !make space to receive entries
             allocate(recentries(nrecentries),stat=ierr)
             call utils_alloc_check('sparse_solve','recentries',ierr)

          end if
          if(need_data.and.(pub_total_num_procs>1)) then
             call comms_wait(send_handle,add_to_stack)
          end if

          if(pub_total_num_procs>1) then
             call comms_barrier()
          end if

          if(need_data.and.(pub_total_num_procs>1)) then
             !send and recv entries
             call comms_send(recv_proc,entriesa(1),size(entriesa), &
                  return_handle=send_handle,add_to_stack=add_to_stack)
          end if
          if(send_data) then
             if(pub_total_num_procs>1) then
                call comms_recv(send_proc,recentries(1),nrecentries)
             else
                recentries(1:nrecentries)=entriesa(1:nrecentries)
             end if

             ! get local subspace data
             call sparse_sizeof_subspace_on_proc(amat,recentries,rowelems, &
                  colelems,ncols_thisproc,nelems_on_proc,first_entry_on_proc)

             if(iscmplx) then
                allocate(zsubspace_buffer_olap(rowelems,nelems_on_proc),stat=ierr)
                call utils_alloc_check('sparse_solve','zsubspace_buffer_olap',ierr)
                zsubspace_buffer_olap=cmplx_0
             else
                allocate(subspace_buffer_olap(rowelems,nelems_on_proc),stat=ierr)
                call utils_alloc_check('sparse_solve','subspace_buffer_olap',ierr)
                subspace_buffer_olap=0.0_dp
             end if

             rowcount=0
             colcount=0

             if(iscmplx) then
                do jcol=first_entry_on_proc,nrecentries
                   if(recentries(jcol)<col_par%my_first_blk) cycle
                   if(recentries(jcol)>col_par%my_last_blk) exit
                   rowcount=0
                   do krow=1,nrecentries
                      call sparse_get_block(zsubspace_buffer_olap(rowcount+1: &
                           rowcount+row_par%num_elems_on_atom(recentries(krow),&
                           row_blks),colcount+1:&
                           colcount+col_par%num_elems_on_atom(recentries(jcol),&
                           col_blks)),bmat,recentries(krow),recentries(jcol))

                      rowcount=rowcount+ &
                           row_par%num_elems_on_atom(recentries(krow),row_blks)
                   end do
                   colcount=colcount+ &
                        col_par%num_elems_on_atom(recentries(jcol),col_blks)
                end do
             else
                do jcol=first_entry_on_proc,nrecentries
                   if(recentries(jcol)<col_par%my_first_blk) cycle
                   if(recentries(jcol)>col_par%my_last_blk) exit
                   rowcount=0
                   do krow=1,nrecentries
                      call sparse_get_block(subspace_buffer_olap(rowcount+1: &
                           rowcount+row_par%num_elems_on_atom(recentries(krow),&
                           row_blks),colcount+1:&
                           colcount+col_par%num_elems_on_atom(recentries(jcol),&
                           col_blks)),bmat,recentries(krow),recentries(jcol))

                      rowcount=rowcount+ &
                           row_par%num_elems_on_atom(recentries(krow),row_blks)
                   end do
                   colcount=colcount+ &
                        col_par%num_elems_on_atom(recentries(jcol),col_blks)
                end do
             end if

             if(iscmplx) then
                allocate(zsubspace_buffer_ham(rowelems,nelems_on_proc),stat=ierr)
                call utils_alloc_check('sparse_solve','zsubspace_buffer_ham',ierr)
                zsubspace_buffer_ham=cmplx_0
             else
                allocate(subspace_buffer_ham(rowelems,nelems_on_proc),stat=ierr)
                call utils_alloc_check('sparse_solve','subspace_buffer_ham',ierr)
                subspace_buffer_ham=0.0_dp
             end if

             rowcount=0
             colcount=0

             if(iscmplx) then
                do jcol=first_entry_on_proc,nrecentries
                   if(recentries(jcol)<col_par%my_first_blk) cycle
                   if(recentries(jcol)>col_par%my_last_blk) exit
                   rowcount=0
                   do krow=1,nrecentries
                      call sparse_get_block(zsubspace_buffer_ham(rowcount+1: &
                           rowcount+row_par%num_elems_on_atom(recentries(krow),&
                           row_blks),colcount+1: &
                           colcount+col_par%num_elems_on_atom(recentries(jcol),&
                           col_blks)),amat,recentries(krow),recentries(jcol))

                      rowcount = rowcount + &
                           row_par%num_elems_on_atom(recentries(krow),row_blks)
                   end do
                   colcount = colcount + &
                        col_par%num_elems_on_atom(recentries(jcol),col_blks)
                end do
             else
                do jcol=first_entry_on_proc,nrecentries
                   if(recentries(jcol)<col_par%my_first_blk) cycle
                   if(recentries(jcol)>col_par%my_last_blk) exit
                   rowcount=0
                   do krow=1,nrecentries
                      call sparse_get_block(subspace_buffer_ham(rowcount+1: &
                           rowcount+row_par%num_elems_on_atom(recentries(krow),&
                           row_blks),colcount+1: &
                           colcount+col_par%num_elems_on_atom(recentries(jcol),&
                           col_blks)),amat,recentries(krow),recentries(jcol))

                      rowcount = rowcount + &
                           row_par%num_elems_on_atom(recentries(krow),row_blks)
                   end do
                   colcount = colcount + &
                        col_par%num_elems_on_atom(recentries(jcol),col_blks)
                end do
             end if

          end if

          if(need_data.and.(pub_total_num_procs>1)) then
             call comms_wait(send_handle,add_to_stack)
          end if

          if(pub_total_num_procs>1) then
             call comms_barrier()
          end if

          if(send_data.and.(pub_total_num_procs>1)) then
             !send sizes to send_proc
             call comms_send(send_proc,nelems_on_proc, &
                  return_handle=send_handle,add_to_stack=add_to_stack)
          end if

          if(need_data) then
             !recv sizes from recv_proc
             if(pub_total_num_procs>1) then
                call comms_recv(recv_proc,ncols_thatproc)
             else
                ncols_thatproc=nelems_on_proc
             end if
          end if

          if(send_data.and.(pub_total_num_procs>1)) then
             call comms_wait(send_handle,add_to_stack)
          end if

          if(pub_total_num_procs>1) then
             call comms_barrier()
          end if

          if(send_data) then
             firstcol_thisproc=1
             do jcol=1,nrecentries
                if(recentries(jcol)<col_par%my_first_blk) then
                   firstcol_thisproc=firstcol_thisproc + &
                        col_par%num_elems_on_atom(recentries(jcol),col_blks)
                else
                   exit
                end if
             end do

             if(pub_total_num_procs>1) then
                call comms_send(send_proc,firstcol_thisproc, &
                     return_handle=send_handle,add_to_stack=add_to_stack)
             end if

             !deallocate entries
             deallocate(recentries,stat=ierr)
             call utils_dealloc_check('sparse_solve','recentries',ierr)
          end if

          if(need_data) then
             if(pub_total_num_procs>1) then
                call comms_recv(recv_proc,firstcol_thatproc)
             else
                firstcol_thatproc=firstcol_thisproc
             end if
          end if

          if(send_data.and.(pub_total_num_procs>1)) then
             call comms_wait(send_handle,add_to_stack)
          end if

          if(pub_total_num_procs>1) then
             call comms_barrier()
          end if

          if(send_data.and.(pub_total_num_procs>1)) then
             !send and recv olap data
             if(rowelems*nelems_on_proc>0) then
                if(iscmplx) then
                   call comms_send(send_proc,zsubspace_buffer_olap(1,1), &
                        rowelems*nelems_on_proc,&
                        return_handle=send_handle,add_to_stack=add_to_stack)
                else
                   call comms_send(send_proc,subspace_buffer_olap(1,1), &
                        rowelems*nelems_on_proc,&
                        return_handle=send_handle,add_to_stack=add_to_stack)
                end if
             end if
          end if
          if(need_data) then
             if(iscmplx) then
                if(size(zsubspace_olap,1)*ncols_thatproc>0) then
                   if(pub_total_num_procs>1) then
                      call comms_recv(recv_proc, &
                           zsubspace_olap(1,firstcol_thatproc), &
                           size(zsubspace_olap,1)*ncols_thatproc)
                   else
                      zsubspace_olap(1:size(zsubspace_olap,1),firstcol_thatproc:firstcol_thatproc+ncols_thatproc-1)=&
                           zsubspace_buffer_olap(1:rowelems,1:nelems_on_proc)
                   end if
                end if
             else
                if(size(subspace_olap,1)*ncols_thatproc>0) then
                   if(pub_total_num_procs>1) then
                      call comms_recv(recv_proc, &
                           subspace_olap(1,firstcol_thatproc), &
                           size(subspace_olap,1)*ncols_thatproc)
                   else
                      subspace_olap(1:size(subspace_olap,1),firstcol_thatproc:firstcol_thatproc+ncols_thatproc-1)=&
                           subspace_buffer_olap(1:rowelems,1:nelems_on_proc)
                   end if
                end if
             end if
          end if

          if(send_data.and.(pub_total_num_procs>1)) then
             call comms_wait(send_handle,add_to_stack)
          end if

          if(pub_total_num_procs>1) then
             call comms_barrier()
          end if

          if(send_data.and.(pub_total_num_procs>1)) then
             !send and recv ham data
             if(rowelems*nelems_on_proc>0) then
                if(iscmplx) then
                   call comms_send(send_proc,zsubspace_buffer_ham(1,1), &
                        rowelems*nelems_on_proc,&
                        return_handle=send_handle,add_to_stack=add_to_stack)
                else
                   call comms_send(send_proc,subspace_buffer_ham(1,1), &
                        rowelems*nelems_on_proc,&
                        return_handle=send_handle,add_to_stack=add_to_stack)
                end if
             end if
          end if
          if(need_data) then
             if(iscmplx) then
                if(size(zsubspace_ham,1)*ncols_thatproc>0) then
                   if(pub_total_num_procs>1) then
                      call comms_recv(recv_proc, &
                           zsubspace_ham(1,firstcol_thatproc), &
                           size(zsubspace_ham,1)*ncols_thatproc)
                   else
                      zsubspace_ham(1:size(zsubspace_ham,1),firstcol_thatproc:firstcol_thatproc+ncols_thatproc-1)=&
                           zsubspace_buffer_ham(1:rowelems,1:nelems_on_proc)
                   end if
                end if
             else
                if(size(subspace_ham,1)*ncols_thatproc>0) then
                   if(pub_total_num_procs>1) then
                      call comms_recv(recv_proc, &
                           subspace_ham(1,firstcol_thatproc), &
                           size(subspace_ham,1)*ncols_thatproc)
                   else
                      subspace_ham(1:size(subspace_ham,1),firstcol_thatproc:firstcol_thatproc+ncols_thatproc-1)=&
                           subspace_buffer_ham(1:rowelems,1:nelems_on_proc)
                   end if
                end if
             end if
          end if

          if(send_data.and.(pub_total_num_procs>1)) then
             call comms_wait(send_handle,add_to_stack)
          end if

          if(pub_total_num_procs>1) then
             call comms_barrier()
          end if

          if(send_data) then
             if(iscmplx) then
                !deallocate buffer
                deallocate(zsubspace_buffer_olap,stat=ierr)
                call utils_dealloc_check('sparse_solve', &
                     'zsubspace_buffer_olap',ierr)

                !deallocate buffer
                deallocate(zsubspace_buffer_ham,stat=ierr)
                call utils_dealloc_check('sparse_solve', &
                     'zsubspace_buffer_ham',ierr)
             else
                !deallocate buffer
                deallocate(subspace_buffer_olap,stat=ierr)
                call utils_dealloc_check('sparse_solve', &
                     'subspace_buffer_olap',ierr)

                !deallocate buffer
                deallocate(subspace_buffer_ham,stat=ierr)
                call utils_dealloc_check('sparse_solve', &
                     'subspace_buffer_ham',ierr)
             end if

          end if


       end do

       if(pub_total_num_procs>1) then
          call comms_barrier()
       end if

       if(need_data) then

          if(iscmplx) then
             allocate(zhopping_sub(size(zsubspace_ham,1), &
                  size(zsubspace_ham,2)),stat=ierr)
             call utils_alloc_check('sparse_solve','zhopping_sub',ierr)
             zhopping_sub(:,:) = zsubspace_ham
             allocate(pivot(size(zsubspace_ham,1)), stat=ierr)
             call utils_alloc_check('sparse_solve','pivot',ierr)
             pivot=0

             call zgesv(size(zsubspace_olap,1), size(zsubspace_ham,1), &
                  zsubspace_olap, size(zsubspace_olap,1), &
                  pivot, zhopping_sub, size(zsubspace_olap,1), info)
             call utils_assert(info==0, &
                  'Error in sparse_solve: ZGESV returned info = ', info)
          else
             allocate(hopping_sub(size(subspace_ham,1), &
                  size(subspace_ham,2)),stat=ierr)
             call utils_alloc_check('sparse_solve','hopping_sub',ierr)
             hopping_sub(:,:) = subspace_ham
             allocate(pivot(size(subspace_ham,1)), stat=ierr)
             call utils_alloc_check('sparse_solve','pivot',ierr)
             pivot=0

             call dgesv(size(subspace_olap,1), size(subspace_ham,1), &
                  subspace_olap, size(subspace_olap,1),&
                  pivot, hopping_sub, size(subspace_olap,1), info)
             call utils_assert(info==0, &
                  'Error in sparse_solve: DGESV returned info = ', info)
          end if


          offset=0
          do ient=1,size(entriesa)
             if(entriesa(ient)>loc_col-1) cycle
             offset=offset+row_par%num_elems_on_atom(entriesa(ient),row_blks)
          end do



          offset2=0
          if(iscmplx) then
             do ient=1,size(entriesa)
                call sparse_put_block(zhopping_sub(offset2+1: &
                     offset2+row_par%num_elems_on_atom(entriesa(ient),row_blks),&
                     offset+1: &
                     offset+col_par%num_elems_on_atom(entriesa(ient),col_blks)), &
                     cmat,entriesa(ient),loc_col)
                offset2=offset2+row_par%num_elems_on_atom(entriesa(ient),row_blks)
             end do
          else
             do ient=1,size(entriesa)
                call sparse_put_block(hopping_sub(offset2+1: &
                     offset2+row_par%num_elems_on_atom(entriesa(ient),row_blks), &
                     offset+1: &
                     offset+col_par%num_elems_on_atom(entriesa(ient),col_blks)), &
                     cmat,entriesa(ient),loc_col)
                offset2=offset2+row_par%num_elems_on_atom(entriesa(ient),row_blks)
             end do
          end if


          if(iscmplx) then
             deallocate(zsubspace_ham,stat=ierr)
             call utils_dealloc_check('sparse_solve','zsubspace_ham',ierr)

             deallocate(zsubspace_olap,stat=ierr)
             call utils_dealloc_check('sparse_solve','zsubspace_olap',ierr)

             deallocate(zhopping_sub,stat=ierr)
             call utils_dealloc_check('sparse_solve','zhopping_sub',ierr)
          else
             deallocate(subspace_ham,stat=ierr)
             call utils_dealloc_check('sparse_solve','subspace_ham',ierr)

             deallocate(subspace_olap,stat=ierr)
             call utils_dealloc_check('sparse_solve','subspace_olap',ierr)

             deallocate(hopping_sub,stat=ierr)
             call utils_dealloc_check('sparse_solve','hopping_sub',ierr)

          end if


          deallocate(pivot,stat=ierr)
          call utils_dealloc_check('sparse_solve','pivot',ierr)

          if(allocated(entriesa)) then
             deallocate(entriesa,stat=ierr)
             call utils_dealloc_check('sparse_solve','entriesa',ierr)
          end if

          if(allocated(entriesb)) then
             deallocate(entriesb,stat=ierr)
             call utils_dealloc_check('sparse_solve','entriesb',ierr)
          end if
       end if

    end do

    if(pub_total_num_procs>1) then
       call comms_barrier()
    end if

  contains

    !============================================================================!
    ! This contained subroutine fetches the indices of the nonzero blocks of a   !
    ! column of a SPAM3 matrix.                                                  !
    !----------------------------------------------------------------------------!
    ! Arguments:                                                                 !
    ! icol (input)     : The desired column from the matrix                      !
    ! mat  (input)     : The SPAM3 matrix                                        !
    ! entries(output)  : The indices of the blocks in the column                 !
    !----------------------------------------------------------------------------!
    ! Written by Jolyon Aarons, May 2017.                                        !
    !============================================================================!
    subroutine get_entries(icol,mat,entries,pmat,entriesp,combine_spaces)
      use comms, only: pub_my_proc_id,comms_reduce,comms_bcast, &
           comms_allgather, comms_send,comms_recv, comms_barrier, &
           pub_total_num_procs
      use constants, only: LONG
      use parallel_strategy, only: PARAL_INFO
      use rundat, only: pub_debug
      use utils, only: utils_assert, utils_alloc_check, &
           utils_dealloc_check, utils_heapsort
      implicit none
      integer,     intent(in) :: icol
      type(SPAM3), intent(in) :: mat
      integer, dimension(:), allocatable, intent(inout) :: entries
      type(SPAM3), optional, intent(in) :: pmat
      integer, dimension(:), allocatable, optional, intent(inout) :: entriesp
      logical, optional, intent(in) :: combine_spaces

      integer(kind=LONG), allocatable, dimension(:) :: comentries
      integer, allocatable, dimension(:) :: pivot

      integer :: jcol,jstart,pstart
      integer :: krow

      integer :: nsubspace, psubspace

      integer :: ierr
      logical :: have_p
      logical :: combine_mats


      combine_mats=.false.
      if(present(combine_spaces)) then
         combine_mats=combine_spaces
      end if
      have_p=.false.
      if(present(pmat)) then
         have_p=.true.
      end if

      if(have_p) then
         call utils_assert(present(entriesp),"Error: if P mat provided, &
              &P subspace must also be provided!")
      end if

      call sparse_get_col_block_count(nsubspace,mat,icol)
      if(have_p) then
         call sparse_get_col_block_count(psubspace,pmat,icol)
      end if


      allocate(entries(nsubspace),stat=ierr)
      call utils_alloc_check('get_entries','entries',ierr)
      entries=-1

      if(have_p) then
         allocate(entriesp(psubspace),stat=ierr)
         call utils_alloc_check('get_entries','entriesp',ierr)
         entriesp=-1
      end if

      call sparse_get_col_block_indices(entries,mat,icol)
      if(have_p) then
         call sparse_get_col_block_indices(entriesp,pmat,icol)
      end if


      if(have_p) then
         if(combine_mats) then
            allocate(comentries(nsubspace+psubspace),stat=ierr)
            call utils_alloc_check('get_entries','comentries',ierr)
            comentries=-1
            comentries(1:nsubspace)=entries
            comentries(1+nsubspace:psubspace+nsubspace)=entriesp
            allocate(pivot(nsubspace+psubspace),stat=ierr)
            call utils_alloc_check('get_entries','pivot',ierr)
            call utils_heapsort(nsubspace+psubspace,comentries,pivot)

            jcol=1
            pivot=0
            pivot(jcol)=comentries(jcol)
            do krow=2,nsubspace+psubspace
               if(comentries(krow)/=pivot(jcol)) then
                  jcol=jcol+1
                  pivot(jcol)=comentries(krow)
               end if
            end do
            deallocate(comentries,stat=ierr)
            call utils_dealloc_check('get_entries','comentries',ierr)
            nsubspace=jcol
            psubspace=jcol

            deallocate(entries,stat=ierr)
            call utils_dealloc_check('get_entries','comentries',ierr)
            deallocate(entriesp,stat=ierr)
            call utils_dealloc_check('get_entries','comentries',ierr)
            allocate(entries(nsubspace),stat=ierr)
            call utils_alloc_check('get_entries','entries',ierr)
            allocate(entriesp(psubspace),stat=ierr)
            call utils_alloc_check('get_entries','entriesp',ierr)

            entries=pivot(1:nsubspace)
            entriesp=pivot(1:psubspace)

            deallocate(pivot,stat=ierr)
            call utils_dealloc_check('get_entries','pivot',ierr)

         end if
      end if

    end subroutine get_entries

  end subroutine sparse_solve2


  !============================================================================!
  ! This subroutine fetches the indices of the nonzero blocks of a column of a !
  ! SPAM3 matrix. The number of blocks must be known beforehand and indices    !
  ! allocated to an appropriate length.                                        !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  ! indices(output)  : The indices of the blocks in the column                 !
  ! mat  (input)     : The SPAM3 matrix                                        !
  ! iblk (input)     : The desired block column (atom index) from the matrix   !
  !----------------------------------------------------------------------------!
  ! Written by Jolyon Aarons, May 2017.                                        !
  !============================================================================!
  subroutine sparse_get_col_block_indices(indices,mat,iblk)

    use comms, only: pub_my_proc_id, pub_total_num_procs
    use parallel_strategy, only: PARAL_INFO
    use rundat, only: pub_debug
    use utils, only: utils_assert

    implicit none

    ! Arguments
    integer, intent(out) :: indices(:)    ! The required col
    type(SPAM3), intent(in) :: mat            ! The matrix
    integer, intent(in) :: iblk               ! The atom index of the column

    ! Local variables
    integer :: ilib      ! Library entry for mat
    integer :: row_blks  ! Identifier for column blocking scheme
    integer :: col_blks
    integer :: iovlap    ! Block overlap counter
    integer :: loc_iblk  ! Block-column counter local to this proc
    integer :: seg       ! Segment index
    integer :: seg_type  ! Segment type for this segment
    integer :: seg_start ! Start position of this segment in the index
    type(PARAL_INFO), pointer :: col_par ! parallel strategy for columns
    integer :: counter

    ! Get library entry
    ilib = mat%lib
    row_blks = library(ilib)%row_blks
    col_blks = library(ilib)%col_blks
    ! rc2013: assign the pointers for par based on library structure
    col_par => library(ilib)%col_par

    ! Check arguments
    call utils_assert(ilib /= 0, &
         'Error in sparse_get_col_block_indices: Matrix does not exist.')
         ! ^ jd: Otherwise uninitialized matrices segfault below

    if(pub_debug) then
       call utils_assert(iblk >= 1 .and. iblk <= library(ilib)%nblk, &
          'Error in sparse_get_col_block_indices: invalid column index.')
    end if

    ! Check this block-column is local to this processor
    call utils_assert(iblk >= col_par%my_first_blk .and. &
         iblk <= col_par%my_last_blk, 'Error in sparse_get_col_block_indices: &
         &requested block is not local to proc ', &
         pub_my_proc_id)

    counter=0

    ! Find local index for the block-column
    loc_iblk = iblk - col_par%my_first_blk + 1
    do seg=0,pub_total_num_procs-1
       ! Segment is blank
       if (library(ilib)%seg_info(s_type,seg) == SEG_BLANK) cycle

!       ! Obtain information about the relevant block-column for this elem

       seg_start = library(ilib)%seg_info(s_idx,seg)
       do iovlap=library(ilib)%blk_idx(seg_start+loc_iblk-1), &
            library(ilib)%blk_idx(seg_start+loc_iblk)-1
          counter=counter+1
          indices(counter) = library(ilib)%blk_idx(iovlap)   ! block-row
       end do

    end do  ! seg

  end subroutine sparse_get_col_block_indices


  !============================================================================!
  ! This subroutine counts the number of nonzero blocks in a column of a       !
  ! SPAM3 matrix. The number of blocks can then be used to allocate a vector   !
  ! of an appropriate length to call sparse_get_col_block_indices.             !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  ! count(output)    : The count of blocks in the column                       !
  ! mat  (input)     : The SPAM3 matrix                                        !
  ! iblk (input)     : The desired block column (atom index) from the matrix   !
  !----------------------------------------------------------------------------!
  ! Written by Jolyon Aarons, May 2017.                                        !
  !============================================================================!
  subroutine sparse_get_col_block_count(count,mat,iblk)

    use comms, only: pub_my_proc_id, pub_total_num_procs
    use parallel_strategy, only: PARAL_INFO
    use rundat, only: pub_debug
    use utils, only: utils_assert

    implicit none

    ! Arguments
    integer, intent(out) :: count    ! The required columncount
    type(SPAM3), intent(in) :: mat            ! The matrix
    integer, intent(in) :: iblk               ! The atom index of the column

    ! Local variables
    integer :: ilib      ! Library entry for mat
    integer :: row_blks  ! Identifier for column blocking scheme
    integer :: col_blks
    integer :: iovlap    ! Block overlap counter
    integer :: loc_iblk  ! Block-column counter local to this proc
    integer :: seg       ! Segment index
    integer :: seg_type  ! Segment type for this segment
    integer :: seg_start ! Start position of this segment in the index
    type(PARAL_INFO), pointer :: col_par ! parallel strategy for columns
    integer :: counter
    integer :: jblk
    real(kind=dp) :: tmpblock(6,6)

    ! Get library entry
    ilib = mat%lib
    row_blks = library(ilib)%row_blks
    col_blks = library(ilib)%col_blks
    ! rc2013: assign the pointers for par based on library structure
    col_par => library(ilib)%col_par

    ! Check arguments
    call utils_assert(ilib /= 0, &
         'Error in sparse_get_col_block_count: Matrix does not exist.')
         ! ^ jd: Otherwise uninitialized matrices segfault below

    if(pub_debug) then
       call utils_assert(iblk >= 1 .and. iblk <= library(ilib)%nblk, &
          'Error in sparse_get_col_block_count: invalid column index.')
    end if

    ! Check this block-column is local to this processor
    call utils_assert(iblk >= col_par%my_first_blk .and. &
         iblk <= col_par%my_last_blk, 'Error in sparse_get_col_block_count: &
         &requested block is not local to proc ', &
         pub_my_proc_id)

    counter=0

    ! Find local index for the block-column
    loc_iblk = iblk - col_par%my_first_blk + 1
    do seg=0,pub_total_num_procs-1
       ! Segment is blank
       if (library(ilib)%seg_info(s_type,seg) == SEG_BLANK) cycle

!       ! Obtain information about the relevant block-column for this elem
       seg_start = library(ilib)%seg_info(s_idx,seg)
       do iovlap=library(ilib)%blk_idx(seg_start+loc_iblk-1), &
            library(ilib)%blk_idx(seg_start+loc_iblk)-1
          jblk = library(ilib)%blk_idx(iovlap)   ! block-row
          counter=counter+1
       end do

    end do  ! seg

    count=counter

  end subroutine sparse_get_col_block_count

  !============================================================================!
  ! This subroutine returns information about the elements and blocks in a     !
  ! SPAM3 matrix on a proc when masked by the pattern in entries.              !
  ! The idea here is that entries contains a number of block indices, if a     !
  ! block row-column pair, i,j are both present in entries, then that block is !
  ! part of the subspace and to be considered in the returned information.     !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  ! mat  (input)     : The SPAM3 matrix                                        !
  ! entries(input)   : The list of blocks in the subspace                      !
  ! rowelems(output) : The total number of row elements in this subspace       !
  ! colelems(output) : The total number of column elements in this subspace    !
  ! nblkcols_on_proc (output) : The number of nonzero column blocks in the     !
  !                             subspace on this proc.                         !
  ! nelems_on_proc (output) : The number / elements in subspace on this proc   !
  ! first_entry_on_proc(output) : The index of the first block on this proc    !
  !----------------------------------------------------------------------------!
  ! Written by Jolyon Aarons, May 2017.                                        !
  !============================================================================!
  subroutine sparse_sizeof_subspace_on_proc(mat,entries,rowelems,colelems, &
       nblkcols_on_proc,nelems_on_proc,first_entry_on_proc)
    use comms, only: pub_my_proc_id
    implicit none
    type(SPAM3), intent(in)    :: mat
    integer,  dimension(:),   intent(in)    :: entries
    integer,     intent(out)   :: rowelems
    integer,     intent(out)   :: colelems
    integer,     intent(out)   :: nblkcols_on_proc
    integer,     intent(out)   :: nelems_on_proc
    integer,     intent(out)   :: first_entry_on_proc

    integer :: ilib
    integer :: row_blks, col_blks
    type(PARAL_INFO), pointer :: col_par
    integer :: nsubspace

    integer :: jcol
    integer :: elemcounter, blkcounter

    ! library entry
    ilib = mat%lib
    row_blks = library(ilib)%row_blks
    col_blks = library(ilib)%col_blks

    ! column parallel strategy
    col_par => library(ilib)%col_par

    nsubspace=size(entries)

    elemcounter=0
    do jcol=1,nsubspace
       ! rc2013: PAR_CHECK!
       elemcounter=elemcounter+col_par%num_elems_on_atom(entries(jcol),row_blks)
    end do
    rowelems=elemcounter

    elemcounter=0
    do jcol=1,nsubspace
       ! rc2013: PAR_CHECK!
       elemcounter=elemcounter+col_par%num_elems_on_atom(entries(jcol),col_blks)
    end do
    colelems=elemcounter

    first_entry_on_proc=-1
    blkcounter=0
    elemcounter=0
    do jcol=1,nsubspace
       if(entries(jcol)<col_par%my_first_blk) cycle
       if(entries(jcol)>col_par%my_last_blk) exit
       if(first_entry_on_proc<0) first_entry_on_proc=jcol
       blkcounter=blkcounter+1
       elemcounter=elemcounter+col_par%num_elems_on_atom(entries(jcol),col_blks)
    end do
    nblkcols_on_proc=blkcounter
    nelems_on_proc=elemcounter
    first_entry_on_proc=max(1,first_entry_on_proc)

  end subroutine sparse_sizeof_subspace_on_proc




  !============================================================================!
  ! This subroutine returns a block of a given block sparse matrix local to    !
  ! this processor.                                                            !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   blk (output)   : The required block of the real matrix                   !
  !   mat (input)    : The real matrix whose element is required               !
  !   jblk (input)   : The atom index of the row required                      !
  !   iblk (input)   : The atom index of the column required                   !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, May 2009 based on SPAM2 version.                 !
  !----------------------------------------------------------------------------!
  ! SPAM2 version written by Peter Haynes, in July 2006.                       !
  !============================================================================!

  subroutine sparse_get_block_real(blk,mat,jblk,iblk)

    use comms, only: pub_my_proc_id
    use parallel_strategy, only: PARAL_INFO
    use rundat, only: pub_debug
    use utils, only: utils_assert

    implicit none

    ! Arguments
    real(kind=DP), intent(out) :: blk(:,:)    ! The required block
    type(SPAM3), intent(in) :: mat            ! The matrix
    integer, intent(in) :: jblk               ! The atom index of the row
    integer, intent(in) :: iblk               ! The atom index of the column

    ! Local variables
    integer :: ilib      ! Library entry for mat
    integer :: row_blks  ! Identifier for column blocking scheme
    integer :: col_blks  ! Identifier for row blocking scheme
    integer :: kblk      ! Block counter
    integer :: ielem     ! Element row counter
    integer :: jelem     ! Element column counter
    integer :: iovlap    ! Block overlap counter
    integer :: loc_iblk  ! Block-column counter local to this proc
    integer :: ptr       ! Pointer to data array
    integer :: seg       ! Segment index
    integer :: seg_type  ! Segment type for this segment
    integer :: seg_start ! Start position of this segment in the index
    integer :: nrows     ! Total number of rows in this block/segment
    type(PARAL_INFO), pointer :: row_par ! Parallel strategy for rows
    type(PARAL_INFO), pointer :: col_par ! parallel strategy for columns

    ! Get library entry
    ilib = mat%lib
    row_blks = library(ilib)%row_blks
    col_blks = library(ilib)%col_blks

    ! rc2013: assign the pointers for par based on library structure
    row_par => library(ilib)%row_par
    col_par => library(ilib)%col_par

    ! Check arguments
    call utils_assert(ilib /= 0, &
         'Error in sparse_get_block_real: Matrix does not exist.')
         ! ^ jd: Otherwise uninitialized matrices segfault below

    call utils_assert(.not. mat%iscmplx, 'Error in sparse_get_block_real: &
         &real matrices only.')
    if(pub_debug) then
       call utils_assert(jblk >= 1 .and. jblk <= library(ilib)%nblk, &
            'Error in sparse_get_block_real: invalid row index.')
       call utils_assert(iblk >= 1 .and. iblk <= library(ilib)%mblk, &
          'Error in sparse_get_block_real: invalid column index.')
       call utils_assert(size(blk,1) >= row_par%num_elems_on_atom(jblk,row_blks), &
          'Error in sparse_get_block_real: too few rows in blk.')
       call utils_assert(size(blk,2) >= col_par%num_elems_on_atom(iblk,col_blks), &
          'Error in sparse_get_block_real: too few columns in blk.')
    end if

    ! Check this block-column is local to this processor
    call utils_assert(iblk >= col_par%my_first_blk .and. iblk <= col_par%my_last_blk, 'Error &
         &in sparse_get_block_real: requested block is not local to proc ', &
         pub_my_proc_id)

    ! Set block to zero
    blk = 0.0_DP
    nrows = 0

    ! Find information about this segment
    seg = row_par%proc_of_elem( &
          row_par%first_elem_on_atom(jblk,row_blks),row_blks)
    seg_type = library(ilib)%seg_info(s_type,seg)
    if (seg_type == SEG_BLANK) return

    ! Find number of rows in this block/segment
    if (seg_type == SEG_DENSE) then
       nrows = row_par%num_elems_on_proc(seg,row_blks)
    else if (seg_type == SEG_SPARSE) then
       nrows = row_par%num_elems_on_atom(jblk,row_blks)
    end if

    ! Find local index for the block-column
    loc_iblk = iblk - col_par%my_first_blk + 1
    seg_start = library(ilib)%seg_info(s_idx,seg)
    do iovlap=library(ilib)%blk_idx(seg_start+loc_iblk-1), &
         library(ilib)%blk_idx(seg_start+loc_iblk)-1
       kblk = library(ilib)%blk_idx(iovlap)   ! block-row
       if (kblk < jblk) cycle
       if (kblk > jblk) exit
       ptr = library(ilib)%blk_ptr(iovlap)
       do ielem=1,col_par%num_elems_on_atom(iblk,col_blks)
          do jelem=1,row_par%num_elems_on_atom(jblk,row_blks)
             blk(jelem,ielem) = mat%dmtx(ptr+(ielem-1)*nrows+jelem-1)
          end do
       end do
    end do

  end subroutine sparse_get_block_real

  !============================================================================!
  ! This subroutine inserts a block of a given block sparse matrix local to    !
  ! this processor.                                                            !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   blk (input)    : The block to insert into the real matrix                !
  !   mat (output)   : The real matrix into which to insert the block          !
  !   jblk (input)   : The atom index of the row to insert                     !
  !   iblk (input)   : The atom index of the column to insert                  !
  !   accum (in, opt): If passed and .true., the block is added to existing    !
  !                    data in the matrix rather than replacing it.            !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, May 2009 based on SPAM2 version.                 !
  ! Trivial extension ('accum') by Jacek Dziedzic, March 2019.                 !
  !----------------------------------------------------------------------------!
  ! SPAM2 version written by Peter Haynes, in July 2006.                       !
  !============================================================================!

  subroutine sparse_put_block_real(blk,mat,jblk,iblk,accum)

    use comms, only: pub_my_proc_id
    use parallel_strategy, only: PARAL_INFO
    use rundat, only: pub_debug
    use utils, only: utils_assert

    implicit none

    ! Arguments
    real(kind=DP), intent(in) :: blk(:,:)     ! The block to insert
    type(SPAM3), intent(inout) :: mat         ! The matrix
    integer, intent(in) :: jblk               ! The atom index of the row
    integer, intent(in) :: iblk               ! The atom index of the column
    logical, intent(in), optional :: accum    ! store mode or accumulate mode

    ! Local variables
    integer :: ilib      ! Library entry for mat
    integer :: row_blks  ! Identifier for column blocking scheme
    integer :: col_blks  ! Identifier for row blocking scheme
    integer :: kblk      ! Block counter
    integer :: ielem     ! Element row counter
    integer :: jelem     ! Element column counter
    integer :: iovlap    ! Block overlap counter
    integer :: loc_iblk  ! Block-column counter local to this proc
    integer :: ptr       ! Pointer to data array
    integer :: seg       ! Segment index
    integer :: seg_type  ! Segment type for this segment
    integer :: seg_start ! Start position of this segment in the index
    integer :: nrows     ! Total number of rows in this block/segment
    logical :: loc_accum
    type(PARAL_INFO), pointer :: row_par ! Parallel strategy for rows
    type(PARAL_INFO), pointer :: col_par ! parallel_strategy for columns

    if(present(accum)) then
       loc_accum = accum
    else
       loc_accum = .false.
    end if

    ! Get library entry
    ilib = mat%lib
    row_blks = library(ilib)%row_blks
    col_blks = library(ilib)%col_blks
    nrows = 0

    ! rc2013: assign the pointers for par based on library structure
    row_par => library(ilib)%row_par
    col_par => library(ilib)%col_par

    ! Check arguments
    if(pub_debug) then
       call utils_assert(jblk >= 1 .and. jblk <= library(ilib)%nblk, &
            'Error in sparse_put_block_real: invalid row index.')
       call utils_assert(iblk >= 1 .and. iblk <= library(ilib)%mblk, &
            'Error in sparse_put_block_real: invalid column index.')
       call utils_assert(size(blk,1) >= row_par%num_elems_on_atom(jblk,row_blks), &
            'Error in sparse_put_block_real: too few rows in blk.')
       call utils_assert(size(blk,2) >= col_par%num_elems_on_atom(iblk,col_blks), &
            'Error in sparse_put_block_real: too few columns in blk.')
    end if

    ! Check this block-column is local to this processor
    call utils_assert(iblk >= col_par%my_first_blk .and. iblk <= col_par%my_last_blk, 'Error &
         &in sparse_put_block_real: requested block is not local to proc ', &
         pub_my_proc_id)

    ! Find information about this segment
    seg = row_par%proc_of_elem(row_par%first_elem_on_atom(jblk,row_blks),row_blks)
    seg_type = library(ilib)%seg_info(s_type,seg)
    if (seg_type == SEG_BLANK) return

    ! Find number of rows in this block/segment
    if (seg_type == SEG_DENSE) then
       nrows = row_par%num_elems_on_proc(seg,row_blks)
    else if (seg_type == SEG_SPARSE) then
       nrows = row_par%num_elems_on_atom(jblk,row_blks)
    end if

    ! Find local index for the block-column
    loc_iblk = iblk - col_par%my_first_blk + 1
    seg_start = library(ilib)%seg_info(s_idx,seg)
    do iovlap=library(ilib)%blk_idx(seg_start+loc_iblk-1), &
         library(ilib)%blk_idx(seg_start+loc_iblk)-1
       kblk = library(ilib)%blk_idx(iovlap)   ! block-row
       if (kblk < jblk) cycle
       if (kblk > jblk) exit
       ptr = library(ilib)%blk_ptr(iovlap)
       if(.not. loc_accum) then
          ! jd: Default mode -- store block
          do ielem=1,col_par%num_elems_on_atom(iblk,col_blks)
             do jelem=1,row_par%num_elems_on_atom(jblk,row_blks)
                if (mat%iscmplx) then
                   mat%zmtx(ptr+(ielem-1)*nrows+jelem-1) = &
                              cmplx(blk(jelem,ielem), 0.0_DP, kind=DP)
                else
                   mat%dmtx(ptr+(ielem-1)*nrows+jelem-1) = blk(jelem,ielem)
                end if
             end do
          end do
       else
          ! jd: Accumulate mode -- add block to existing data
          do ielem=1,col_par%num_elems_on_atom(iblk,col_blks)
             do jelem=1,row_par%num_elems_on_atom(jblk,row_blks)
                if (mat%iscmplx) then
                   mat%zmtx(ptr+(ielem-1)*nrows+jelem-1) = &
                        mat%zmtx(ptr+(ielem-1)*nrows+jelem-1) + &
                        cmplx(blk(jelem,ielem), 0.0_DP, kind=DP)
                else
                   mat%dmtx(ptr+(ielem-1)*nrows+jelem-1) = &
                        mat%dmtx(ptr+(ielem-1)*nrows+jelem-1) + &
                        blk(jelem,ielem)
                end if
             end do
          end do
       end if
    end do

  end subroutine sparse_put_block_real

  !============================================================================!
  ! This subroutine returns a column of a given block sparse matrix local to   !
  ! this processor.                                                            !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   data (output)  : The required column of the real matrix                  !
  !   mat (input)    : The real matrix whose column is required                !
  !   elem (input)   : The element index of the column required                !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, May 2009 based on SPAM2 version.                 !
  !----------------------------------------------------------------------------!
  ! SPAM2 version written by Peter Haynes, in May 2004.                        !
  !============================================================================!

  subroutine sparse_get_col_real(data,mat,elem)

    use comms, only: pub_my_proc_id, pub_total_num_procs
    use parallel_strategy, only: PARAL_INFO
    use utils, only: utils_assert

    implicit none

    ! Arguments
    real(kind=DP), intent(inout) :: data(:)    ! The required column
    type(SPAM3), intent(in) :: mat             ! The matrix
    integer, intent(in) :: elem                ! The elem index of the col

    ! Local variables
    integer :: ilib      ! Library entry for mat
    integer :: row_blks  ! Identifier for column blocking scheme
    integer :: col_blks  ! Identifier for row blocking scheme
    integer :: iblk      ! The block-column of the requested column
    integer :: ielem0    ! The first elem in the requested block-column
    integer :: iovlap    ! Block overlap counter
    integer :: loc_iblk  ! Block-column counter local to this proc
    integer :: jblk      ! Block-row counter
    integer :: jelems    ! Number of rows in block-row
    integer :: jelem0    ! First row in block-row
    integer :: ptr       ! Pointer to data array
    integer :: seg       ! Segment index
    integer :: seg_start ! Start position of this segment in the index
    type(PARAL_INFO), pointer :: row_par ! Parallel strategy for rows
    type(PARAL_INFO), pointer :: col_par ! parallel strategy for columns

    ! Get library entry
    ilib = mat%lib
    row_blks = library(ilib)%row_blks
    col_blks = library(ilib)%col_blks

    ! rc2013: assign the pointers for par based on library structure
    row_par => library(ilib)%row_par
    col_par => library(ilib)%col_par

    ! Check arguments
    call utils_assert(.not. mat%iscmplx, 'Error in sparse_get_col_real: &
         &real matrices only.')
    call utils_assert(size(data) >= library(ilib)%nrows, &
          'Error in sparse_get_col_real: insufficient space in output array.')
    call utils_assert(elem >= 1 .and. elem <= library(ilib)%mcols, &
         'Error in sparse_get_col_real: invalid column index.')

    ! Check this column is local to this processor
    call utils_assert(elem>=sparse_first_elem_on_proc(pub_my_proc_id,mat,'C') &
         .and. elem<=sparse_last_elem_on_proc(pub_my_proc_id,mat,'C'), &
         'Error in sparse_get_col_real: requested column is not local &
         &to proc ', pub_my_proc_id)

    ! Loop over proc-proc segments
    do seg=0,pub_total_num_procs-1


       ! Segment is dense
       if (library(ilib)%seg_info(s_type,seg) == SEG_DENSE) then

          jelems = row_par%num_elems_on_proc(seg,row_blks)
          ptr = library(ilib)%seg_info(s_ptr,seg) + &
               (elem - col_par%first_elem_on_proc(pub_my_proc_id,col_blks)) * &
               row_par%num_elems_on_proc(seg,row_blks)

          data(sparse_first_elem_on_proc(seg,mat,'R'):sparse_last_elem_on_proc(seg, &
               mat,'R')) = mat%dmtx(ptr:ptr+jelems-1)

       ! Segment is sparse
       else if (library(ilib)%seg_info(s_type,seg) == SEG_SPARSE) then

          ! Obtain information about the relevant block-column for this elem
          iblk = col_par%atom_of_elem(elem,col_blks)
          ielem0 = col_par%first_elem_on_atom(iblk,col_blks)

          ! Find local index for the block-column
          loc_iblk = iblk - col_par%my_first_blk + 1
          seg_start = library(ilib)%seg_info(s_idx,seg)
          do iovlap=library(ilib)%blk_idx(seg_start+loc_iblk-1), &
               library(ilib)%blk_idx(seg_start+loc_iblk)-1
             jblk = library(ilib)%blk_idx(iovlap)   ! block-row
             jelems = row_par%num_elems_on_atom(jblk,row_blks)
             jelem0 = row_par%first_elem_on_atom(jblk,row_blks)
             ptr = library(ilib)%blk_ptr(iovlap) + (elem - ielem0) * jelems
             data(jelem0:jelem0+jelems-1) = mat%dmtx(ptr:ptr+jelems-1)
          end do

       end if

    end do  ! seg

  end subroutine sparse_get_col_real

  !============================================================================!
  ! This subroutine inserts a column of a given block sparse matrix local to   !
  ! this processor.                                                            !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   data (input)   : The column of the real matrix to insert                 !
  !   mat  (output)  : The real matrix whose column is required                !
  !   elem (input)   : The element index of the column to insert               !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, May 2009 based on SPAM2 version.                 !
  !----------------------------------------------------------------------------!
  ! SPAM2 version written by Peter Haynes, in May 2004.                        !
  !============================================================================!

  subroutine sparse_put_col_real(data,mat,elem)

    use comms, only: pub_my_proc_id, pub_total_num_procs
    use parallel_strategy, only: PARAL_INFO
    use utils, only: utils_assert

    implicit none

    ! Arguments
    real(kind=DP), intent(in) :: data(:)      ! The required column
    type(SPAM3), intent(inout) :: mat         ! The matrix
    integer, intent(in) :: elem               ! The elem index of the col

    ! Local variables
    integer :: ilib      ! Library entry for mat
    integer :: row_blks  ! Identifier for column blocking scheme
    integer :: col_blks  ! Identifier for row blocking scheme
    integer :: iblk      ! The block-column of the requested column
    integer :: ielem0    ! The first elem in the requested block-column
    integer :: iovlap    ! Block overlap counter
    integer :: loc_iblk  ! Block-column counter local to this proc
    integer :: jblk      ! Block-row counter
    integer :: jelems    ! Number of rows in block-row
    integer :: jelem0    ! First row in block-row
    integer :: ptr       ! Pointer to data array
    integer :: seg       ! Segment index
    integer :: seg_start ! Start position of this segment in the index
    type(PARAL_INFO), pointer :: row_par ! Parallel strategy for rows
    type(PARAL_INFO), pointer :: col_par ! parallel strategy for columns

    ! get library entry
    ilib = mat%lib
    row_blks = library(ilib)%row_blks
    col_blks = library(ilib)%col_blks

    ! rc2013: assign the pointers for par based on library structure
    row_par => library(ilib)%row_par
    col_par => library(ilib)%col_par

    ! Check arguments
    call utils_assert(size(data) >= library(ilib)%nrows, &
         'Error in sparse_put_col_real: insufficient space in output array.')
    call utils_assert (elem >= 1 .and. elem <= library(ilib)%mcols, &
         'Error in sparse_put_col_real: invalid column index.')

    ! Check this column is local to this processor
    call utils_assert(elem>=sparse_first_elem_on_proc(pub_my_proc_id,mat,'C').and.&
         elem<=sparse_last_elem_on_proc(pub_my_proc_id,mat,'C'), 'Error in &
         &sparse_put_col_real: requested column is not local to proc ', &
         pub_my_proc_id)

    ! Loop over proc-proc segments
    do seg=0,pub_total_num_procs-1

       ! Segment is dense
       if (library(ilib)%seg_info(s_type,seg) == SEG_DENSE) then

          jelems = row_par%num_elems_on_proc(seg,row_blks)
          ptr = library(ilib)%seg_info(s_ptr,seg) + &
               (elem - col_par%first_elem_on_proc(pub_my_proc_id,col_blks)) * &
               row_par%num_elems_on_proc(seg,row_blks)

          if (mat%iscmplx) then
             mat%zmtx(ptr:ptr+jelems-1) = cmplx(data( &
                  row_par%first_elem_on_proc(seg,row_blks): &
                  row_par%first_elem_on_proc(seg+1,row_blks)-1), 0.0_DP,kind=DP)
          else
             mat%dmtx(ptr:ptr+jelems-1) = data( &
                  row_par%first_elem_on_proc(seg,row_blks): &
                  row_par%first_elem_on_proc(seg+1,row_blks)-1)
          end if

       ! Segment is sparse
       else if (library(ilib)%seg_info(s_type,seg) == SEG_SPARSE) then

          ! Obtain information about the relevant block-column for this elem
          iblk = col_par%atom_of_elem(elem,col_blks)
          ielem0 = col_par%first_elem_on_atom(iblk,col_blks)

          ! Find local index for the block-column
          loc_iblk = iblk - col_par%my_first_blk + 1
          seg_start = library(ilib)%seg_info(s_idx,seg)
          do iovlap=library(ilib)%blk_idx(seg_start+loc_iblk-1), &
               library(ilib)%blk_idx(seg_start+loc_iblk)-1
             jblk = library(ilib)%blk_idx(iovlap)   ! block-row
             jelems = row_par%num_elems_on_atom(jblk,row_blks)
             jelem0 = row_par%first_elem_on_atom(jblk,row_blks)
             ptr = library(ilib)%blk_ptr(iovlap) + (elem - ielem0) * jelems
             if (mat%iscmplx) then
                mat%zmtx(ptr:ptr+jelems-1) = &
                           cmplx(data(jelem0:jelem0+jelems-1), 0.0_DP, kind=DP)
             else
                mat%dmtx(ptr:ptr+jelems-1) = data(jelem0:jelem0+jelems-1)
             end if
          end do

       ! Nothing to do if segment is blank

       end if

    end do  ! seg

  end subroutine sparse_put_col_real

  !============================================================================!
  ! This subroutine clears an array according to the sparsity pattern of a     !
  ! column of a given block sparse matrix local to this processor.             !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   data (output)  : The real array to clear                                 !
  !   mat  (output)  : The real matrix whose sparsity pattern is to be used    !
  !   elem (input)   : The element index of the column to use                  !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, May 2009 based on SPAM2 version.                 !
  !----------------------------------------------------------------------------!
  ! SPAM2 version written by Peter Haynes, in May 2004.                        !
  !============================================================================!

  subroutine sparse_clr_col_real(data,mat,elem)

    use comms, only: pub_total_num_procs, pub_my_proc_id
    use parallel_strategy, only: PARAL_INFO
    use utils, only: utils_assert

    implicit none

    ! Arguments
    real(kind=DP), intent(inout) :: data(:)    ! The column to insert
    type(SPAM3), intent(in) :: mat             ! The matrix
    integer, intent(in) :: elem                ! The elem index of the col

    ! Local variables
    integer :: ilib      ! Library entry for mat
    integer :: row_blks  ! Identifier for column blocking scheme
    integer :: col_blks  ! Identifier for row blocking scheme
    integer :: iblk      ! The block-column of the requested column
    integer :: iovlap    ! Block overlap counter
    integer :: loc_iblk  ! Block-column counter local to this proc
    integer :: jblk      ! Block-row counter
    integer :: jelems    ! Number of rows in block-row
    integer :: jelem0    ! First row in block-row
    integer :: seg       ! Segment index
    integer :: seg_start ! Start position of this segment in the index
    type(PARAL_INFO), pointer :: row_par ! Parallel strategy for rows
    type(PARAL_INFO), pointer :: col_par ! parallel_strategy for columns

    ! Get library entry
    ilib = mat%lib
    row_blks = library(ilib)%row_blks
    col_blks = library(ilib)%col_blks

    ! rc2013: assign the pointers for par based on library structure
    row_par => library(ilib)%row_par
    col_par => library(ilib)%col_par

    ! Check arguments
    call utils_assert(.not. mat%iscmplx, &
         'Error in sparse_clr_col_real: real matrices only.')
    call utils_assert(size(data) >= library(ilib)%nrows, &
         'Error in sparse_clr_col_real: insufficient space in input array.')
    call utils_assert(elem >= 1 .and. elem <= library(ilib)%mcols, &
         'Error in sparse_clr_col_real: invalid column index.')

    ! Obtain information about the relevant block-column for this elem
    iblk = col_par%atom_of_elem(elem,col_blks)

    ! Check this block-column is local to this processor
    call utils_assert(iblk >= col_par%my_first_blk .and. &
         iblk <= col_par%my_last_blk, 'Error in sparse_clr_col_real: &
         &requested column is not local to proc ', pub_my_proc_id)

    ! Find local index for the block-column
    loc_iblk = iblk - col_par%my_first_blk + 1

    do seg=0,pub_total_num_procs-1

       ! Nothing to do if segment is blank
       if (library(ilib)%seg_info(s_type,seg)==SEG_BLANK) cycle

       ! Loop over the nonzero blocks of this segment
       seg_start = library(ilib)%seg_info(s_idx,seg)
       do iovlap=library(ilib)%blk_idx(seg_start+loc_iblk-1), &
            library(ilib)%blk_idx(seg_start+loc_iblk)-1
          jblk = library(ilib)%blk_idx(iovlap)   ! block-row

          ! Set the nonzero elements of this col to zero
          jelems = row_par%num_elems_on_atom(jblk,row_blks)
          jelem0 = row_par%first_elem_on_atom(jblk,row_blks)
          data(jelem0:jelem0+jelems-1) = 0.0_DP

       end do  ! iovlap

    end do  ! seg

  end subroutine sparse_clr_col_real


  !============================================================================!
  ! This subroutine converts a block sparse matrix in a SPAM3 into a dense     !
  ! full square real array, or collates a dense matrix across procs to one     !
  ! single full square array                                                   !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   dest   (output) : The destination full square matrix                     !
  !   src    (input)  : The source block sparse matrix                         !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, June 2009.                                       !
  !============================================================================!

  subroutine sparse_spam3tofull_real(dest,src)

    use comms, only: comms_bcast, pub_total_num_procs
    use parallel_strategy, only: PARAL_INFO
    use utils, only: utils_assert

    implicit none

    ! Arguments
    real(kind=DP), intent(out) :: dest(:,:)     ! The full square matrix
    type(SPAM3), intent(in) :: src              ! The block sparse matrix

    ! Local variables
    integer :: srclib        ! Library entry for src
    integer :: row_blks      ! Identifier for column blocking scheme
    integer :: col_blks      ! Identifier for row blocking scheme
    integer :: iblk,jblk     ! Block counters
    integer :: loc_iblk      ! Local counter for iblk
    integer :: idx           ! Index
    integer :: ptr           ! Pointer
    integer :: ielemonat     ! elem in atom block iblk counter
    integer :: jelemonat     ! elem in atom block jblk counter
    integer :: ielem,jelem   ! elem counters
    integer :: ielems,jelems ! Numbers of element rows/cols in this block
    integer :: nrows         ! Total number of rows in this block/segment
    integer :: proc          ! Proc loop counter
    integer :: seg           ! Segment index
    integer :: seg_type      ! Start position of this segment in the index
    integer :: seg_start     ! Start position of this segment in the index
    type(PARAL_INFO), pointer :: row_par ! Parallel strategy for rows
    type(PARAL_INFO), pointer :: col_par ! parallel strategy for columns

    ! Get library entry for src
    srclib = src%lib
    row_blks = library(srclib)%row_blks
    col_blks = library(srclib)%col_blks
    nrows = 0

    ! rc2013: assign the pointers for par based on library structure
    row_par => library(srclib)%row_par
    col_par => library(srclib)%col_par

    ! Check arguments
    call utils_assert(size(dest,1) == library(srclib)%nrows, 'Error in &
         &sparse_spam3tofull_real: incompatible numbers of rows in arguments.')
    call utils_assert(size(dest,2) == library(srclib)%mcols, 'Error in &
         &sparse_spam3tofull_real: incompatible numbers of cols in arguments.')
    call utils_assert(.not. src%iscmplx, &
         'Error in sparse_spam3tofull_real: real matrices only.')

    ! Zero full square matrix
    dest(:,:) = 0.0_DP

    ! Loop over the segments of src on this proc
    do seg=0,pub_total_num_procs-1

       ! Cycle if this segment is blank
       seg_type = library(srclib)%seg_info(s_type,seg)
       if (seg_type == SEG_BLANK) cycle
       if (seg_type == SEG_DENSE) nrows = row_par%num_elems_on_proc(seg,row_blks)

       ! Loop over block-columns of src on this proc
       seg_start = library(srclib)%seg_info(s_idx,seg)
       loc_iblk = 0
       do iblk=col_par%my_first_blk,col_par%my_last_blk
          loc_iblk = loc_iblk + 1
          ielem = col_par%first_elem_on_atom(iblk,col_blks)
          ielems = col_par%num_elems_on_atom(iblk,col_blks)

          ! Loop over block-rows jblk in column iblk of src
          do idx=library(srclib)%blk_idx(seg_start+loc_iblk-1), &
               library(srclib)%blk_idx(seg_start+loc_iblk)-1
             jblk = library(srclib)%blk_idx(idx)
             jelem = row_par%first_elem_on_atom(jblk,row_blks)
             jelems = row_par%num_elems_on_atom(jblk,row_blks)

             if (seg_type == SEG_SPARSE) nrows = jelems

             ! Copy block
             ptr = library(srclib)%blk_ptr(idx)
             do ielemonat=0,ielems-1      ! Element col in block
                do jelemonat=0,jelems-1   ! Element row in block
                   dest(jelem+jelemonat,ielem+ielemonat) = &
                        src%dmtx(ptr+ielemonat*nrows+jelemonat)
                end do
             end do

          end do  ! Loop over block-rows in column iblk of src

       end do  ! Loop over block-columns of src

    end do  ! seg

    ! Broadcast result across all procs
    do proc=0,pub_total_num_procs-1

       ! Get first column and number of columns
       ielem = col_par%first_elem_on_proc(proc,col_blks)
       ielems = col_par%num_elems_on_proc(proc,col_blks)

       ! Broadcast this part of the matrix
       ! rc2013: only if there's something here...
       if (ielem .ne. 0 .and. ielems .gt. 0) &
            call comms_bcast(proc,dest(1,ielem),ielems*library(srclib)%nrows)
       !call comms_bcast(proc,dest(1,ielem),ielems*library(srclib)%nrows)

    end do  ! Loop over procs

  end subroutine sparse_spam3tofull_real

  !============================================================================!
  ! This subroutine converts a block sparse matrix in a SPAM3 into a dense     !
  ! full square real array, or collates a dense matrix across procs to one     !
  ! single full square array                                                   !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   dest   (output) : The destination full square matrix                     !
  !   src    (input)  : The source block sparse matrix                         !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, June 2009.                                       !
  !============================================================================!

  subroutine sparse_fulltospam3_real(dest,src)

    use comms, only: pub_total_num_procs
    use parallel_strategy, only: PARAL_INFO
    use utils, only: utils_assert

    implicit none

    ! Arguments
    type(SPAM3), intent(inout) :: dest       ! The block sparse matrix
    real(kind=DP), intent(in) :: src(:,:)    ! The full square matrix

    ! Local variables
    integer :: destlib       ! Library entry for dest
    integer :: row_blks      ! Identifier for column blocking scheme
    integer :: col_blks      ! Identifier for row blocking scheme
    integer :: iblk,jblk     ! Block counters
    integer :: loc_iblk      ! Local counter for iblk
    integer :: idx           ! Index
    integer :: ptr           ! Pointer
    integer :: ielemonat     ! Element in atom block iblk counter
    integer :: jelemonat     ! Element in atom block jblk counter
    integer :: ielem,jelem   ! Element start positions of a block
    integer :: ielems,jelems ! Numbers of element rows/cols in this block
    integer :: nrows         ! Total number of rows in the block/segment
    integer :: seg           ! Segment index
    integer :: seg_type      ! Segment type for this segment
    integer :: seg_start     ! Start position of this segment in the index
    type(PARAL_INFO), pointer :: row_par ! Parallel strategy for rows
    type(PARAL_INFO), pointer :: col_par ! parallel strategy for columns

    ! Get library entry for src
    destlib = dest%lib
    row_blks = library(destlib)%row_blks
    col_blks = library(destlib)%col_blks
    nrows = 0

    ! rc2013: assign the pointers for par based on library structure
    row_par => library(destlib)%row_par
    col_par => library(destlib)%col_par

    ! Check arguments
    call utils_assert(size(src,1) == library(destlib)%nrows, 'Error in &
         &sparse_fulltospam3_real: incompatible numbers of rows in arguments.')
    call utils_assert(size(src,2) == library(destlib)%mcols, 'Error in &
         &sparse_fulltospam3_real: incompatible numbers of cols in arguments.')
    call utils_assert(.not. dest%iscmplx, &
         'Error in sparse_fulltospam3_real: real matrices only.')

    ! Loop over the segments of dest on this proc
    do seg=0,pub_total_num_procs-1

       ! Cycle if this segment is blank
       seg_type = library(destlib)%seg_info(s_type,seg)
       if (seg_type == SEG_BLANK) cycle
       if (seg_type == SEG_DENSE) nrows = row_par%num_elems_on_proc(seg,row_blks)

       ! Loop over block-columns of dest on this proc
       seg_start = library(destlib)%seg_info(s_idx,seg)
       loc_iblk = 0
       do iblk=col_par%my_first_blk,col_par%my_last_blk
          loc_iblk = loc_iblk + 1
          ielems = col_par%num_elems_on_atom(iblk,col_blks)

          ! Loop over block-rows jblk in column iblk of dest
          do idx=library(destlib)%blk_idx(seg_start+loc_iblk-1), &
               library(destlib)%blk_idx(seg_start+loc_iblk)-1
             jblk = library(destlib)%blk_idx(idx)
             jelems = row_par%num_elems_on_atom(jblk,row_blks)
             if (seg_type == SEG_SPARSE) nrows = jelems

             ! Copy block
             ielem = col_par%first_elem_on_atom(iblk,col_blks)
             ptr = library(destlib)%blk_ptr(idx)
             do ielemonat=0,ielems-1      ! Element col in block
                jelem = row_par%first_elem_on_atom(jblk,row_blks)
                do jelemonat=0,jelems-1   ! Element row in block
                   dest%dmtx(ptr+ielemonat*nrows+jelemonat) = &
                        src(jelem+jelemonat,ielem+ielemonat)
                end do
             end do

          end do  ! Loop over block-rows in column iblk of dest

       end do  ! Loop over block-columns of dest

    end do  ! seg

  end subroutine sparse_fulltospam3_real

  !============================================================================!
  ! This subroutine scales and optionally shifts the eigenvalues of a matrix.  !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   mat    (inout) : The matrix to be rescaled and shifted                   !
  !   alpha  (input) : The real scaling parameter                              !
  !   beta   (input) : The optional real shift parameter                       !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, May 2009 based on SPAM2 version.                 !
  !----------------------------------------------------------------------------!
  ! SPAM2 version written by Peter Haynes, March 2004.                         !
  ! Revised for distributed data, June 2006.                                   !
  ! Dense matrix version by Nicholas Hine, Dec 2007.                           !
  !============================================================================!

  subroutine sparse_scale_real(mat,alpha,beta)

    use comms, only: pub_my_proc_id
    use parallel_strategy, only: PARAL_INFO
    use utils, only: utils_assert

    implicit none

    ! Arguments
    type(SPAM3), intent(inout)          :: mat    ! The matrix to be operated on
    real(kind=DP), intent(in)           :: alpha  ! The scaling parameter
    real(kind=DP), optional, intent(in) :: beta   ! The shift parameter

    ! Local variables
    integer :: ilib      ! Library entry for mat
    integer :: blks      ! Identifier for blocking scheme
    integer :: iblk      ! Block-column loop counter
    integer :: loc_iblk  ! Loop counter for iblk on this proc
    integer :: ielems    ! Number of columns in block
    integer :: ptr       ! Pointer to diagonal entries
    integer :: ielem     ! Element in block-column loop counter
    integer :: seg       ! Segment index
    integer :: seg_start ! Segment type for this segment
    type(PARAL_INFO), pointer :: row_par ! Parallel strategy for rows
    type(PARAL_INFO), pointer :: col_par ! Parallel strategy fotr col_par


    ! Rescale the elements if required
    if (alpha /= 1.0_DP) then
       if (mat%iscmplx) then
          mat%zmtx = alpha * mat%zmtx
       else
          mat%dmtx = alpha * mat%dmtx
       end if
    end if

    ! Shift the eigenvalues if required
    if (present(beta)) then

       ! Get library entry for mat
       ilib = mat%lib
       blks = library(ilib)%row_blks
       seg = pub_my_proc_id

       ! rc2013: assign the pointers for par based on library structure
       row_par => library(ilib)%row_par
       col_par => library(ilib)%col_par

       ! Can only shift eigenvalues of square matrix
       call utils_assert(library(ilib)%row_blks == library(ilib)%col_blks, &
            'Error in sparse_scale: cannot shift eigenvalues of non-square &
            &matrices.')

       if (library(ilib)%seg_info(s_type,seg)==SEG_DENSE) then

          ! Set pointer to first diagonal element
          ptr = library(ilib)%seg_info(s_ptr,seg)

          ! Get number of columns in this segment
          ielems = col_par%num_elems_on_proc(seg,blks)

          ! Add identity matrix scaled by beta to all elements
          if (mat%iscmplx) then
             ! Loop over columns on my proc
             do ielem=1,ielems
                mat%zmtx(ptr) = mat%zmtx(ptr) + cmplx(beta,0.0_DP,kind=DP)
                ptr = ptr + ielems + 1
             end do
          else
             ! Loop over columns on my proc
             do ielem=1,col_par%num_elems_on_proc(seg,blks)
                mat%dmtx(ptr) = mat%dmtx(ptr) + beta
                ptr = ptr + ielems + 1
             end do
          end if

       else if (library(ilib)%seg_info(s_type,seg)==SEG_SPARSE) then

          seg_start = library(ilib)%seg_info(s_idx,seg)

          ! Loop over atom-blocks iblk on my proc
          loc_iblk = 0
          do iblk=col_par%my_first_blk,col_par%my_last_blk
             loc_iblk = loc_iblk + 1

             ! Get number of elements on this atom
             ielems = col_par%num_elems_on_atom(iblk,blks)

             ! Find pointer to start of diagonal block in this column
             ptr = library(ilib)%blk_ptr(seg_start+loc_iblk-1)

             ! Loop over diagonal elements in this diagonal block and add
             ! identity matrix scaled by beta
             if (mat%iscmplx) then
                do ielem=1,ielems
                   mat%zmtx(ptr) = mat%zmtx(ptr) + cmplx(beta,0.0_DP,kind=DP)
                   ptr = ptr + ielems + 1
                end do
             else
                do ielem=1,ielems
                   mat%dmtx(ptr) = mat%dmtx(ptr) + beta
                   ptr = ptr + ielems + 1
                end do
             end if

          end do  ! Loop over atom-blocks iblk on this proc

       endif

    end if


  end subroutine sparse_scale_real

  !============================================================================!
  ! This subroutine sets a SPAM3 matrix to garbage.
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   matrix (in/out): The matrix to operate on.                               !
  !----------------------------------------------------------------------------!
  ! Written by Jacek Dziedzic on 2016.10.14.                                   !
  !============================================================================!

  subroutine sparse_set_to_garbage(matrix)

    use constants, only: garbage_real, garbage_complex
    use timer, only: timer_clock

    implicit none

    type(SPAM3), intent(inout) :: matrix

    ! -------------------------------------------------------------------------

    ! jd: To make sure the overhead of invalidating matrices in negligible
    call timer_clock('sparse_set_to_garbage',1)

    if(matrix%iscmplx) then
       matrix%zmtx = garbage_complex
    else
       matrix%dmtx = garbage_real
       ! rc2013: EMBED_FIX!
       !matrix%dmtx = 0.0_DP
    end if

    ! jd: Zap all elements that "don't exist", but are stored in dense blocks
    if (any(library(matrix%lib)%seg_info(s_type,:)==SEG_DENSE)) then
       call sparse_enforce_sparsity(matrix)
    end if

    call timer_clock('sparse_set_to_garbage',2)

  end subroutine sparse_set_to_garbage

  !============================================================================!
  ! This subroutine performs an axpy operation on the matrices:                !
  !   y := y + alpha * x                                                       !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   ymat  (inout) : The sparse matrix y                                      !
  !   xmat  (input) : The sparse matrix x                                      !
  !   alpha (input) : The real parameter alpha                                 !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, May 2009 based on SPAM2 version.                 !
  !----------------------------------------------------------------------------!
  ! SPAM2 version written by Peter Haynes, March 2004.                         !
  ! Revised for distributed data, June 2006.                                   !
  ! Dense matrix version by Nicholas Hine, Dec 2007.                           !
  ! Modified by Andrea Greco to allow safe conversion from complex xmat to     !
  ! real ymat, September 2016.                                                 !
  !============================================================================!

  subroutine sparse_axpy_real(ymat,xmat,alpha,cmplx_to_real)

    use comms, only: pub_total_num_procs
    use parallel_strategy, only: PARAL_INFO
!$  use rundat, only: pub_threads_max
    use utils, only: utils_assert, utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(SPAM3), intent(inout) :: ymat    ! The sparse matrix y
    type(SPAM3), intent(in) :: xmat       ! The sparse matrix x
    real(kind=DP), intent(in) :: alpha    ! The parameter alpha
    ! agrecocmplx
    logical, optional, intent(in) :: cmplx_to_real

    ! Local variables
    integer :: xlib,ylib    ! Library pointers for x and y
    integer :: row_blks     ! Identifier for column blocking scheme
    integer :: col_blks     ! Identifier for row blocking scheme
    integer :: iblk,jblk    ! Atom loop counters
    integer :: icol         ! Atom loop counters
    integer :: loc_iblk     ! Atom loop counter on this proc
    integer :: idx          ! Index loop counter
    integer :: iel          ! Element loop counter
    integer :: ifound       ! Loop counter over columns found
    integer :: nfound       ! Number of columns found
    integer :: xstart,xend  ! Start and end of block data of xmat
    integer :: ystart,yend  ! Start and end of block data of ymat
    integer :: ystride      ! Stride to get to next column of src
    integer :: xstride      ! Stride to get to next column of dest
    integer :: seg          ! Segment index
    integer :: x_seg_type   ! Segment type for segments of x
    integer :: y_seg_type   ! Segment type for segments of y
    integer :: x_seg_start  ! Start position of segments of x in the index
    integer :: y_seg_start  ! Start position of segments of y in the index
    integer :: ierr         ! Error flag
    integer, allocatable :: tfound_idx(:), tfound_ptr(:)
    logical, allocatable :: tfound_blk(:)
    type(PARAL_INFO), pointer :: row_par ! Parallel strategy for rows
    type(PARAL_INFO), pointer :: col_par ! Parallel strategy for cols
    logical :: loc_cmplx_to_real ! agrecocmplx

    ! agrecocmplx
    loc_cmplx_to_real = .false.
    if (present(cmplx_to_real)) loc_cmplx_to_real = cmplx_to_real

    ! agrecocmplx: protection  against unsafe conversion from complex to real
    call utils_assert(.not.(xmat%iscmplx).or.(ymat%iscmplx).or.loc_cmplx_to_real, &
         'ERROR: unsafe conversion of complex matrix to real. Call routine &
         &sparse_axpy_real with cmplx_to_real=.true. to force the conversion.')

    ! Extract library entries
    xlib = xmat%lib
    ylib = ymat%lib
    row_blks = library(xlib)%row_blks
    col_blks = library(xlib)%col_blks
    xstride = 0; ystride = 0

    ! rc2013: assign the pointers for par based on library structure
    ! rc2013: x and y must have the same parallel strategies
    row_par => library(xlib)%row_par
    col_par => library(xlib)%col_par

    ! Check matrices share same blocking scheme
    call utils_assert( (library(ylib)%row_blks == row_blks) .and. &
         (library(ylib)%col_blks == col_blks), 'Error in sparse_axpy_real: &
         &x and y matrices have different blocking schemes.')

    ! rc2013: check for consistency of parallel strategies
    call sparse_par_check(ylib, xlib, routine='sparse_axpy')

    ! Check for alpha = 0 (no operation required)
    if (alpha == 0.0_DP) then
       return
    end if

    ! Check whether the matrices share the same structure
    if ( xlib == ylib ) then

       ! Trivial axpy of whole data arrays
       if (xmat%iscmplx) then
          if (ymat%iscmplx) then
             ymat%zmtx = ymat%zmtx + alpha * xmat%zmtx
          else
             do iel=1,size(ymat%dmtx)
                ymat%dmtx(iel) = ymat%dmtx(iel) + &
                     alpha * real(xmat%zmtx(iel),kind=DP)
             end do
          end if
       else
          if (ymat%iscmplx) then
             do iel=1,size(ymat%zmtx)
                ymat%zmtx(iel) = ymat%zmtx(iel) + &
                     cmplx(alpha * xmat%dmtx(iel),0.0_DP,kind=DP)
             end do
          else
             ymat%dmtx = ymat%dmtx + alpha * xmat%dmtx
          end if
       end if

    ! Structures differ, so proceed segment-by-segment
    else

!$OMP PARALLEL NUM_THREADS(pub_threads_max) DEFAULT(none) &
!$OMP PRIVATE(seg,x_seg_type,y_seg_type,loc_iblk,iblk,idx,jblk,nfound, &
!$OMP      ifound,x_seg_start,y_seg_start,xstart,xend,ystart,yend,xstride, &
!$OMP      ystride,icol,iel,tfound_blk,tfound_ptr,tfound_idx,ierr) &
!$OMP SHARED(xmat,ymat,xlib,ylib,library,row_blks,col_blks,alpha, &
!$OMP      row_par,col_par,pub_total_num_procs,pub_threads_max)

       allocate(tfound_idx(library(xlib)%nblk),stat=ierr)
       call utils_alloc_check('sparse_axpy_real','tfound_idx',ierr)
       allocate(tfound_ptr(library(xlib)%nblk),stat=ierr)
       call utils_alloc_check('sparse_axpy_real','tfound_ptr',ierr)
       allocate(tfound_blk(library(xlib)%nblk),stat=ierr)
       call utils_alloc_check('sparse_axpy_real','tfound_blk',ierr)
       tfound_blk = .false.

!$OMP DO
       do seg=0,pub_total_num_procs-1

          x_seg_type = library(xlib)%seg_info(s_type,seg)
          y_seg_type = library(ylib)%seg_info(s_type,seg)

          ! Nothing to add to if this segment in y is blank
          if (y_seg_type==SEG_BLANK) cycle

          ! Nothing to add if this segment in x is blank
          if (x_seg_type == SEG_BLANK) cycle

          ! Find information about these segments
          x_seg_start = library(xlib)%seg_info(s_idx,seg)
          y_seg_start = library(ylib)%seg_info(s_idx,seg)
          if (y_seg_type == SEG_DENSE) ystride = &
               row_par%num_elems_on_proc(seg,row_blks)
          if (x_seg_type == SEG_DENSE) xstride = &
               row_par%num_elems_on_proc(seg,row_blks)

          ! Loop over atom blocks (block-columns) on this proc
          loc_iblk = 0
          do iblk=col_par%my_first_blk,col_par%my_last_blk
             loc_iblk = loc_iblk + 1

             ! Reset counters of number of rows found
             nfound = 0

             ! Loop over block-rows in x for this segment
             do idx=library(xlib)%blk_idx(x_seg_start+loc_iblk-1), &
                  library(xlib)%blk_idx(x_seg_start+loc_iblk)-1
                jblk = library(xlib)%blk_idx(idx)

                ! Mark flag for this row and add to list
                nfound = nfound + 1
                tfound_blk(jblk) = .true.
                tfound_idx(nfound) = jblk
                tfound_ptr(jblk) = library(xlib)%blk_ptr(idx)

             end do

             ! Loop over block-rows in y
             do idx=library(ylib)%blk_idx(y_seg_start+loc_iblk-1), &
                  library(ylib)%blk_idx(y_seg_start+loc_iblk)-1
                jblk = library(ylib)%blk_idx(idx)

                ! Check whether this block-row exists in source
                if (tfound_blk(jblk)) then

                   ! Find pointers to start of dest and src data
                   ystart = library(ylib)%blk_ptr(idx)
                   yend = ystart + row_par%num_elems_on_atom(jblk,row_blks) - 1
                   xstart = tfound_ptr(jblk)
                   xend = xstart + row_par%num_elems_on_atom(jblk,row_blks) - 1

                   ! Find stride to next column of dest data
                   if (y_seg_type == SEG_SPARSE) then
                      ystride = row_par%num_elems_on_atom(jblk,row_blks)
                   end if

                   ! Find stride to next column of source data
                   if (x_seg_type == SEG_SPARSE) then
                      xstride = row_par%num_elems_on_atom(jblk,row_blks)
                   end if

                   ! Copy columns from source to destination
                   do icol=1,col_par%num_elems_on_atom(iblk,col_blks)

                      if (xmat%iscmplx) then
                         if (ymat%iscmplx) then
                            ymat%zmtx(ystart:yend) = &
                                 ymat%zmtx(ystart:yend) + &
                                 alpha * xmat%zmtx(xstart:xend)
                         else
                            do iel=ystart,yend
                               ymat%dmtx(iel) = ymat%dmtx(iel) + alpha * &
                                    real(xmat%zmtx(iel-ystart+xstart), &
                                    kind=DP)
                            end do
                         end if
                      else
                         if (ymat%iscmplx) then
                            do iel=ystart,yend
                               ymat%zmtx(iel) = ymat%zmtx(iel) + cmplx( &
                                    alpha * xmat%dmtx(iel-ystart+xstart), &
                                    0.0_DP,kind=DP)
                            end do
                         else
                            ymat%dmtx(ystart:yend) = &
                                 ymat%dmtx(ystart:yend) + &
                                 alpha * xmat%dmtx(xstart:xend)
                         end if
                      end if

                      ystart = ystart + ystride
                      yend = yend + ystride
                      xstart = xstart + xstride
                      xend = xend + xstride

                   end do

                end if

             end do  ! Loop over destination block-rows

             ! Reset flags
             do ifound=1,nfound
                jblk = tfound_idx(ifound)
                tfound_blk(jblk) = .false.
             end do

          end do  ! Loop over atom block-columns

       end do  ! seg
!$OMP END DO

      deallocate(tfound_blk,stat=ierr)
      call utils_dealloc_check('sparse_axpy_real','tfound_blk',ierr)
      deallocate(tfound_ptr,stat=ierr)
      call utils_dealloc_check('sparse_axpy_real','tfound_ptr',ierr)
      deallocate(tfound_idx,stat=ierr)
      call utils_dealloc_check('sparse_axpy_real','tfound_idx',ierr)

!$OMP END PARALLEL

    end if


  end subroutine sparse_axpy_real


  !============================================================================!
  ! This subroutine inserts an element into a given block sparse matrix local  !
  ! to this processor.                                                         !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   el (input)     : The element of the real matrix to insert                !
  !   mat (inout)    : The real matrix                                         !
  !   jrow (input)   : The element index of the row required                   !
  !   icol (input)   : The element index of the column required                !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, May 2009 based on SPAM2 version.                 !
  !----------------------------------------------------------------------------!
  ! SPAM2 version written by Peter Haynes, in June 2006.                       !
  !============================================================================!

  subroutine sparse_put_element_complex(el,mat,jrow,icol)

    use comms, only: pub_my_proc_id
    use parallel_strategy, only: PARAL_INFO
    use rundat, only: pub_debug
    use utils, only: utils_assert

    implicit none

    ! Arguments
    complex(kind=DP), intent(in) :: el        ! The element to insert
    type(SPAM3), intent(inout) :: mat         ! The matrix
    integer, intent(in) :: jrow               ! The index of the row
    integer, intent(in) :: icol               ! The index of the column

    ! Local variables
    integer :: ilib      ! Library entry for mat
    integer :: row_blks  ! Identifier for column blocking scheme
    integer :: col_blks  ! Identifier for row blocking scheme
    integer :: iblk      ! The block-row to which the element belongs
    integer :: jblk      ! The block-col to which the element belongs
    integer :: kblk      ! Block counter
    integer :: nrows     ! Total number of rows in this block/segment
    integer :: jelem0    ! The first element in the requested block-row
    integer :: ielem0    ! The first element in the requested block-column
    integer :: iovlap    ! Block overlap counter
    integer :: loc_iblk  ! Block-column counter local to this proc
    integer :: ptr       ! Pointer to data array
    integer :: seg       ! Segment index
    integer :: seg_type  ! Segment type for this segment
    integer :: seg_start ! Start position of this segment in the index
    type(PARAL_INFO), pointer :: row_par ! parallel strategy for rows
    type(PARAL_INFO), pointer :: col_par ! parallel strategy for columns

    ! Get library entry
    ilib = mat%lib
    row_blks = library(ilib)%row_blks
    col_blks = library(ilib)%col_blks
    nrows = 0

    ! rc2013: assign the pointers for par based on library structure
    row_par => library(ilib)%row_par
    col_par => library(ilib)%col_par

    ! Check arguments
    call utils_assert(mat%iscmplx, &
         'Error in sparse_put_element_complex: complex matrices only.')
    if(pub_debug) then
       call utils_assert(jrow >= 1 .and. jrow <= library(ilib)%nrows, &
            'Error in sparse_put_element_complex: invalid row index.')
       call utils_assert(icol >= 1 .and. icol <= library(ilib)%mcols, &
            'Error in sparse_put_element_complex: invalid column index.')
    end if
    ! Check this column is local to this processor
    call utils_assert(icol>=col_par%first_elem_on_proc(pub_my_proc_id,col_blks).and.&
         icol<col_par%first_elem_on_proc(pub_my_proc_id+1,col_blks), 'Error in &
         &sparse_put_element_complex: requested column is not local to proc ', &
         pub_my_proc_id)

    ! Obtain information about the relevant block & segment for this element
    seg = row_par%proc_of_elem(jrow,row_blks)
    seg_type = library(ilib)%seg_info(s_type,seg)
    if (seg_type==SEG_BLANK) return
    seg_start = library(ilib)%seg_info(s_idx,seg)
    jblk = row_par%atom_of_elem(jrow,row_blks)
    iblk = col_par%atom_of_elem(icol,col_blks)
    jelem0 = row_par%first_elem_on_atom(jblk,row_blks)
    ielem0 = col_par%first_elem_on_atom(iblk,col_blks)
    if (seg_type==SEG_DENSE) then
       nrows = row_par%num_elems_on_proc(seg,row_blks)
    else if (seg_type==SEG_SPARSE) then
       nrows = row_par%num_elems_on_atom(jblk,row_blks)
    end if

    ! Find local index for the block-column
    loc_iblk = iblk - col_par%my_first_blk + 1

    do iovlap=library(ilib)%blk_idx(seg_start+loc_iblk-1), &
         library(ilib)%blk_idx(seg_start+loc_iblk)-1
       kblk = library(ilib)%blk_idx(iovlap)   ! block-row
       if (kblk < jblk) cycle
       if (kblk > jblk) exit
       ptr = library(ilib)%blk_ptr(iovlap)
       mat%zmtx(ptr + (icol - ielem0)*nrows + (jrow - jelem0)) = el
    end do

  end subroutine sparse_put_element_complex

  !============================================================================!
  ! This subroutine retrieves an element into a given block sparse matrix      !
  ! local to this processor.                                                   !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   el (output)    : The element of the real matrix retrieved                !
  !   mat (inout)    : The real matrix                                         !
  !   jrow (input)   : The element index of the row required                   !
  !   icol (input)   : The element index of the column required                !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, May 2009 based on SPAM2 version.                 !
  !----------------------------------------------------------------------------!
  ! SPAM2 version written by Peter Haynes, in June 2006.                       !
  !============================================================================!

  subroutine sparse_get_element_complex(el,mat,jrow,icol)

    use comms, only: pub_my_proc_id
    use constants, only: cmplx_0
    use parallel_strategy, only: PARAL_INFO
    use rundat, only: pub_debug
    use utils, only: utils_assert

    implicit none

    ! Arguments
    complex(kind=DP), intent(out) :: el       ! The element to insert
    type(SPAM3), intent(in) :: mat            ! The matrix
    integer, intent(in) :: jrow               ! The index of the row
    integer, intent(in) :: icol               ! The index of the column

    ! Local variables
    integer :: ilib      ! Library entry for mat
    integer :: row_blks  ! Identifier for column blocking scheme
    integer :: col_blks  ! Identifier for row blocking scheme
    integer :: iblk      ! The block-row to which the element belongs
    integer :: jblk      ! The block-col to which the element belongs
    integer :: kblk      ! Block counter
    integer :: nrows     ! Total number of rows in this block/segment
    integer :: jelem0    ! The first element in the requested block-row
    integer :: ielem0    ! The first element in the requested block-column
    integer :: iovlap    ! Block overlap counter
    integer :: loc_iblk  ! Block-column counter local to this proc
    integer :: ptr       ! Pointer to data array
    integer :: seg       ! Segment index
    integer :: seg_type  ! Segment type for this segment
    integer :: seg_start ! Start position of this segment in the index
    type(PARAL_INFO), pointer :: row_par ! parallel strategy for rows
    type(PARAL_INFO), pointer :: col_par ! parallel strategy for columns

    ! Get library entry
    ilib = mat%lib
    row_blks = library(ilib)%row_blks
    col_blks = library(ilib)%col_blks

    ! rc2013: assign the pointers for par based on library structure
    row_par => library(ilib)%row_par
    col_par => library(ilib)%col_par

    ! Check arguments
    call utils_assert(mat%iscmplx, &
         'Error in sparse_get_element_complex: complex matrices only.')
    if(pub_debug) then
       call utils_assert(jrow >= 1 .and. jrow <= library(ilib)%nrows, &
            'Error in sparse_get_element_complex: invalid row index.')
       call utils_assert(icol >= 1 .and. icol <= library(ilib)%mcols, &
            'Error in sparse_get_element_complex: invalid column index.')
    end if
    ! Check this column is local to this processor
    call utils_assert(icol>=sparse_first_elem_on_proc(pub_my_proc_id,mat,'C').and.&
         icol<=sparse_last_elem_on_proc(pub_my_proc_id,mat,'C'), 'Error in &
         &sparse_get_element_complex: requested column is not local to proc ', &
         pub_my_proc_id)

    nrows = 0
    el = cmplx_0

    ! Obtain information about the relevant block & segment for this element
    seg = row_par%proc_of_elem(jrow,row_blks)
    seg_type = library(ilib)%seg_info(s_type,seg)
    if (seg_type==SEG_BLANK) return
    seg_start = library(ilib)%seg_info(s_idx,seg)
    jblk = row_par%atom_of_elem(jrow,row_blks)
    iblk = col_par%atom_of_elem(icol,col_blks)
    jelem0 = row_par%first_elem_on_atom(jblk,row_blks)
    ielem0 = col_par%first_elem_on_atom(iblk,col_blks)
    if (seg_type==SEG_DENSE) then
       nrows = row_par%num_elems_on_proc(seg,row_blks)
    else if (seg_type==SEG_SPARSE) then
       nrows = row_par%num_elems_on_atom(jblk,row_blks)
    end if

    ! Find local index for the block-column
    loc_iblk = iblk - col_par%my_first_blk + 1

    do iovlap=library(ilib)%blk_idx(seg_start+loc_iblk-1), &
         library(ilib)%blk_idx(seg_start+loc_iblk)-1
       kblk = library(ilib)%blk_idx(iovlap)   ! block-row
       if (kblk < jblk) cycle
       if (kblk > jblk) exit
       ptr = library(ilib)%blk_ptr(iovlap)
       el = mat%zmtx(ptr + (icol - ielem0)*nrows + (jrow - jelem0))
    end do

  end subroutine sparse_get_element_complex

  !============================================================================!
  ! This subroutine returns a block of a given block sparse matrix local to    !
  ! this processor.                                                            !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   blk (output)   : The required block of the real matrix                   !
  !   mat (input)    : The real matrix whose element is required               !
  !   jblk (input)   : The atom index of the row required                      !
  !   iblk (input)   : The atom index of the column required                   !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, May 2009 based on SPAM2 version.                 !
  !----------------------------------------------------------------------------!
  ! SPAM2 version written by Peter Haynes, in July 2006.                       !
  !============================================================================!

  subroutine sparse_get_block_complex(blk,mat,jblk,iblk)

    use comms, only: pub_my_proc_id
    use constants, only: cmplx_0
    use parallel_strategy, only: PARAL_INFO
    use rundat, only: pub_debug
    use utils, only: utils_assert

    implicit none

    ! Arguments
    complex(kind=DP), intent(out) :: blk(:,:) ! The required block
    type(SPAM3), intent(in) :: mat            ! The matrix
    integer, intent(in) :: jblk               ! The atom index of the row
    integer, intent(in) :: iblk               ! The atom index of the column

    ! Local variables
    integer :: ilib      ! Library entry for mat
    integer :: row_blks  ! Identifier for column blocking scheme
    integer :: col_blks  ! Identifier for row blocking scheme
    integer :: kblk      ! Block counter
    integer :: ielem     ! Element row counter
    integer :: jelem     ! Element column counter
    integer :: iovlap    ! Block overlap counter
    integer :: loc_iblk  ! Block-column counter local to this proc
    integer :: ptr       ! Pointer to data array
    integer :: seg       ! Segment index
    integer :: seg_type  ! Segment type for this segment
    integer :: seg_start ! Start position of this segment in the index
    integer :: nrows     ! Total number of rows in this block/segment
    type(PARAL_INFO), pointer :: row_par
    type(PARAL_INFO), pointer :: col_par

    ! Get library entry
    ilib = mat%lib
    row_blks = library(ilib)%row_blks
    col_blks = library(ilib)%col_blks

    ! rc2013: assign the pointers for par based on library structure
    row_par => library(ilib)%row_par
    col_par => library(ilib)%col_par

    ! Check arguments
    call utils_assert(mat%iscmplx, &
         'Error in sparse_get_block_complex: complex matrices only.')
    if(pub_debug) then
       call utils_assert(jblk >= 1 .and. jblk <= library(ilib)%nblk, &
            'Error in sparse_get_block_complex: invalid row index.')
       call utils_assert(iblk >= 1 .and. iblk <= library(ilib)%nblk, &
            'Error in sparse_get_block_complex: invalid column index.')
       call utils_assert(size(blk,1) >= row_par%num_elems_on_atom(jblk,row_blks), &
            'Error in sparse_get_block_complex: too few rows in blk.')
       call utils_assert(size(blk,2) >= col_par%num_elems_on_atom(iblk,col_blks), &
            'Error in sparse_get_block_complex: too few columns in blk.')
    end if

    ! Check this block-column is local to this processor
    call utils_assert(iblk >= col_par%my_first_blk .and. iblk <= col_par%my_last_blk, 'Error &
         &in sparse_get_block_complex: requested block is not local to proc ', &
         pub_my_proc_id)

    ! Set block to zero
    blk = cmplx_0
    nrows = 0

    ! Find information about this segment
    seg = row_par%proc_of_elem(row_par%first_elem_on_atom(jblk,row_blks),row_blks)
    seg_type = library(ilib)%seg_info(s_type,seg)
    if (seg_type == SEG_BLANK) return

    ! Find number of rows in this block/segment
    if (seg_type == SEG_DENSE) then
       nrows = row_par%num_elems_on_proc(seg,row_blks)
    else if (seg_type == SEG_SPARSE) then
       nrows = row_par%num_elems_on_atom(jblk,row_blks)
    end if

    ! Find local index for the block-column
    loc_iblk = iblk - col_par%my_first_blk + 1
    seg_start = library(ilib)%seg_info(s_idx,seg)
    do iovlap=library(ilib)%blk_idx(seg_start+loc_iblk-1), &
         library(ilib)%blk_idx(seg_start+loc_iblk)-1
       kblk = library(ilib)%blk_idx(iovlap)   ! block-row
       if (kblk < jblk) cycle
       if (kblk > jblk) exit
       ptr = library(ilib)%blk_ptr(iovlap)
       do ielem=1,col_par%num_elems_on_atom(iblk,col_blks)
          do jelem=1,row_par%num_elems_on_atom(jblk,row_blks)
             blk(jelem,ielem) = mat%zmtx(ptr+(ielem-1)*nrows+jelem-1)
          end do
       end do
    end do

  end subroutine sparse_get_block_complex

  !============================================================================!
  ! This subroutine inserts a block of a given block sparse matrix local to    !
  ! this processor.                                                            !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   blk (input)    : The block to insert into the real matrix                !
  !   mat (output)   : The real matrix into which to insert the block          !
  !   jblk (input)   : The atom index of the row to insert                     !
  !   iblk (input)   : The atom index of the column to insert                  !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, May 2009 based on SPAM2 version.                 !
  !----------------------------------------------------------------------------!
  ! SPAM2 version written by Peter Haynes, in July 2006.                       !
  !============================================================================!

  subroutine sparse_put_block_complex(blk,mat,jblk,iblk)

    use comms, only: pub_my_proc_id
    use parallel_strategy, only: PARAL_INFO
    use rundat, only: pub_debug
    use utils, only: utils_assert

    implicit none

    ! Arguments
    complex(kind=DP), intent(in) :: blk(:,:)  ! The block to insert
    type(SPAM3), intent(inout) :: mat         ! The matrix
    integer, intent(in) :: jblk               ! The atom index of the row
    integer, intent(in) :: iblk               ! The atom index of the column

    ! Local variables
    integer :: ilib      ! Library entry for mat
    integer :: row_blks  ! Identifier for column blocking scheme
    integer :: col_blks  ! Identifier for row blocking scheme
    integer :: kblk      ! Block counter
    integer :: ielem     ! Element row counter
    integer :: jelem     ! Element column counter
    integer :: iovlap    ! Block overlap counter
    integer :: loc_iblk  ! Block-column counter local to this proc
    integer :: ptr       ! Pointer to data array
    integer :: seg       ! Segment index
    integer :: seg_type  ! Segment type for this segment
    integer :: seg_start ! Start position of this segment in the index
    integer :: nrows     ! Total number of rows in this block/segment
    type(PARAL_INFO), pointer :: row_par
    type(PARAL_INFO), pointer :: col_par

    ! Get library entry
    ilib = mat%lib
    row_blks = library(ilib)%row_blks
    col_blks = library(ilib)%col_blks
    nrows = 0

    ! rc2013: assign the pointers for par based on library structure
    row_par => library(ilib)%row_par
    col_par => library(ilib)%col_par

    ! Check arguments
    call utils_assert(mat%iscmplx, &
         'Error in sparse_put_block_complex: complex matrices only.')
    if(pub_debug) then
       call utils_assert(jblk >= 1 .and. jblk <= library(ilib)%nblk, &
            'Error in sparse_put_block_complex: invalid row index.')
       call utils_assert(iblk >= 1 .and. iblk <= library(ilib)%nblk, &
            'Error in sparse_put_block_complex: invalid column index.')
       call utils_assert(size(blk,1) >= row_par%num_elems_on_atom(jblk,row_blks), &
            'Error in sparse_put_block_complex: too few rows in blk.')
       call utils_assert(size(blk,2) >= col_par%num_elems_on_atom(iblk,col_blks), &
            'Error in sparse_get_block_complex: too few columns in blk.')
    end if

    ! Check this block-column is local to this processor
    call utils_assert(iblk >= col_par%my_first_blk .and. iblk <= col_par%my_last_blk, 'Error &
         &in sparse_put_block_complex: requested block is not local to proc ', &
         pub_my_proc_id)

    ! Find information about this segment
    seg = row_par%proc_of_elem(row_par%first_elem_on_atom(jblk,row_blks),row_blks)
    seg_type = library(ilib)%seg_info(s_type,seg)
    if (seg_type == SEG_BLANK) return

    ! Find number of rows in this block/segment
    if (seg_type == SEG_DENSE) then
       nrows = row_par%num_elems_on_proc(seg,row_blks)
    else if (seg_type == SEG_SPARSE) then
       nrows = row_par%num_elems_on_atom(jblk,row_blks)
    end if

    ! Find local index for the block-column
    loc_iblk = iblk - col_par%my_first_blk + 1
    seg_start = library(ilib)%seg_info(s_idx,seg)
    do iovlap=library(ilib)%blk_idx(seg_start+loc_iblk-1), &
         library(ilib)%blk_idx(seg_start+loc_iblk)-1
       kblk = library(ilib)%blk_idx(iovlap)   ! block-row
       if (kblk < jblk) cycle
       if (kblk > jblk) exit
       ptr = library(ilib)%blk_ptr(iovlap)
       do ielem=1,col_par%num_elems_on_atom(iblk,col_blks)
          do jelem=1,row_par%num_elems_on_atom(jblk,row_blks)
             mat%zmtx(ptr+(ielem-1)*nrows+jelem-1) = blk(jelem,ielem)
          end do
       end do
    end do

  end subroutine sparse_put_block_complex

  !============================================================================!
  ! This subroutine returns a column of a given block sparse matrix local to   !
  ! this processor.                                                            !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   data (output)  : The required column of the real matrix                  !
  !   mat (input)    : The real matrix whose column is required                !
  !   elem (input)   : The element index of the column required                !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, May 2009 based on SPAM2 version.                 !
  !----------------------------------------------------------------------------!
  ! SPAM2 version written by Peter Haynes, in May 2004.                        !
  !============================================================================!

  subroutine sparse_get_col_complex(data,mat,elem)

    use comms, only: pub_my_proc_id, pub_total_num_procs
    use parallel_strategy, only: PARAL_INFO
    use utils, only: utils_assert

    implicit none

    ! Arguments
    complex(kind=DP), intent(inout) :: data(:) ! The required column
    type(SPAM3), intent(in) :: mat             ! The matrix
    integer, intent(in) :: elem                ! The elem index of the col

    ! Local variables
    integer :: ilib      ! Library entry for mat
    integer :: row_blks  ! Identifier for column blocking scheme
    integer :: col_blks  ! Identifier for row blocking scheme
    integer :: iblk      ! The block-column of the requested column
    integer :: ielem0    ! The first elem in the requested block-column
    integer :: iovlap    ! Block overlap counter
    integer :: loc_iblk  ! Block-column counter local to this proc
    integer :: jblk      ! Block-row counter
    integer :: jelems    ! Number of rows in block-row
    integer :: jelem0    ! First row in block-row
    integer :: ptr       ! Pointer to data array
    integer :: seg       ! Segment index
    integer :: seg_start ! Start position of this segment in the index
    type(PARAL_INFO), pointer :: row_par ! Parallel strategy for rows
    type(PARAL_INFO), pointer :: col_par ! Parallel strategy for columns

    ! Get library entry
    ilib = mat%lib
    row_blks = library(ilib)%row_blks
    col_blks = library(ilib)%col_blks

    ! rc2013: assign the pointers for par based on library structure
    row_par => library(ilib)%row_par
    col_par => library(ilib)%col_par

    ! Check arguments
    call utils_assert(mat%iscmplx, &
         'Error in sparse_get_col_complex: complex matrices only.')
    call utils_assert(size(data) >= library(ilib)%nrows, &
         'Error in sparse_get_col_complex: insufficient space in output array.')
    call utils_assert(elem >= 1 .and. elem <= library(ilib)%mcols, &
         'Error in sparse_get_col_complex: invalid column index.')

    ! Check this column is local to this processor
    call utils_assert(elem>=sparse_first_elem_on_proc(pub_my_proc_id,mat,'C').and.&
         elem<=sparse_last_elem_on_proc(pub_my_proc_id,mat,'C'), 'Error in &
         &sparse_get_col_complex: requested column is not local to proc ', &
         pub_my_proc_id)

    ! Loop over proc-proc segments
    do seg=0,pub_total_num_procs-1

       ! Segment is dense
       if (library(ilib)%seg_info(s_type,seg) == SEG_DENSE) then

          jelems = row_par%num_elems_on_proc(seg,row_blks)
          ptr = library(ilib)%seg_info(s_ptr,seg) + &
               (elem - col_par%first_elem_on_proc(pub_my_proc_id,col_blks)) * &
               row_par%num_elems_on_proc(seg,row_blks)

          data(row_par%first_elem_on_proc(seg,row_blks):row_par%first_elem_on_proc(seg+1, &
               row_blks)-1) = mat%zmtx(ptr:ptr+jelems-1)

       ! Segment is sparse
       else if (library(ilib)%seg_info(s_type,seg) == SEG_SPARSE) then

          ! Obtain information about the relevant block-column for this elem
          iblk = col_par%atom_of_elem(elem,col_blks)
          ielem0 = col_par%first_elem_on_atom(iblk,col_blks)

          ! Find local index for the block-column
          loc_iblk = iblk - col_par%my_first_blk + 1
          seg_start = library(ilib)%seg_info(s_idx,seg)
          do iovlap=library(ilib)%blk_idx(seg_start+loc_iblk-1), &
               library(ilib)%blk_idx(seg_start+loc_iblk)-1
             jblk = library(ilib)%blk_idx(iovlap)   ! block-row
             jelems = row_par%num_elems_on_atom(jblk,row_blks)
             jelem0 = row_par%first_elem_on_atom(jblk,row_blks)
             ptr = library(ilib)%blk_ptr(iovlap) + (elem - ielem0) * jelems
             data(jelem0:jelem0+jelems-1) = mat%zmtx(ptr:ptr+jelems-1)
          end do

       end if

    end do  ! seg

  end subroutine sparse_get_col_complex

  !============================================================================!
  ! This subroutine inserts a column of a given block sparse matrix local to   !
  ! this processor.                                                            !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   data (input)   : The column of the real matrix to insert                 !
  !   mat  (output)  : The real matrix whose column is required                !
  !   elem (input)   : The element index of the column to insert               !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, May 2009 based on SPAM2 version.                 !
  !----------------------------------------------------------------------------!
  ! SPAM2 version written by Peter Haynes, in May 2004.                        !
  !============================================================================!

  subroutine sparse_put_col_complex(data,mat,elem)

    use comms, only: pub_my_proc_id, pub_total_num_procs
    use parallel_strategy, only: PARAL_INFO
    use utils, only: utils_assert

    implicit none

    ! Arguments
    complex(kind=DP), intent(in) :: data(:)   ! The required column
    type(SPAM3), intent(inout) :: mat         ! The matrix
    integer, intent(in) :: elem               ! The elem index of the col

    ! Local variables
    integer :: ilib      ! Library entry for mat
    integer :: row_blks  ! Identifier for column blocking scheme
    integer :: col_blks  ! Identifier for row blocking scheme
    integer :: iblk      ! The block-column of the requested column
    integer :: ielem0    ! The first elem in the requested block-column
    integer :: iovlap    ! Block overlap counter
    integer :: loc_iblk  ! Block-column counter local to this proc
    integer :: jblk      ! Block-row counter
    integer :: jelems    ! Number of rows in block-row
    integer :: jelem0    ! First row in block-row
    integer :: ptr       ! Pointer to data array
    integer :: seg       ! Segment index
    integer :: seg_start ! Start position of this segment in the index
    type(PARAL_INFO), pointer :: row_par ! Parallel strategy for rows
    type(PARAL_INFO), pointer :: col_par ! Parallel strategy for columns

    ! Get library entry
    ilib = mat%lib
    row_blks = library(ilib)%row_blks
    col_blks = library(ilib)%col_blks

    ! rc2013: assign the pointers for par based on library structure
    row_par => library(ilib)%row_par
    col_par => library(ilib)%col_par

    ! Check arguments
    call utils_assert(mat%iscmplx, &
         'Error in sparse_put_col_complex: complex matrices only.')
    call utils_assert(size(data) >= library(ilib)%nrows, &
         'Error in sparse_put_col_complex: insufficient space in output array.')
    call utils_assert(elem >= 1 .and. elem <= library(ilib)%mcols, &
         'Error in sparse_put_col_complex: invalid column index.')

    ! Check this column is local to this processor
    call utils_assert(elem>=sparse_first_elem_on_proc(pub_my_proc_id,mat,'C').and.&
         elem<=sparse_last_elem_on_proc(pub_my_proc_id,mat,'C'), 'Error in &
         &sparse_put_col_complex: requested column is not local to proc ', &
         pub_my_proc_id)

    ! Loop over proc-proc segments
    do seg=0,pub_total_num_procs-1

       ! Segment is dense
       if (library(ilib)%seg_info(s_type,seg) == SEG_DENSE) then

          jelems = row_par%num_elems_on_proc(seg,row_blks)
          ptr = library(ilib)%seg_info(s_ptr,seg) + &
               (elem - col_par%first_elem_on_proc(pub_my_proc_id,col_blks)) * &
               row_par%num_elems_on_proc(seg,row_blks)

          mat%zmtx(ptr:ptr+jelems-1) = data(row_par%first_elem_on_proc(seg,row_blks): &
               row_par%first_elem_on_proc(seg+1,row_blks)-1)

       ! Segment is sparse
       else if (library(ilib)%seg_info(s_type,seg) == SEG_SPARSE) then

          ! Obtain information about the relevant block-column for this elem
          iblk = col_par%atom_of_elem(elem,col_blks)
          ielem0 = col_par%first_elem_on_atom(iblk,col_blks)

          ! Find local index for the block-column
          loc_iblk = iblk - col_par%my_first_blk + 1
          seg_start = library(ilib)%seg_info(s_idx,seg)
          do iovlap=library(ilib)%blk_idx(seg_start+loc_iblk-1), &
               library(ilib)%blk_idx(seg_start+loc_iblk)-1
             jblk = library(ilib)%blk_idx(iovlap)   ! block-row
             jelems = row_par%num_elems_on_atom(jblk,row_blks)
             jelem0 = row_par%first_elem_on_atom(jblk,row_blks)
             ptr = library(ilib)%blk_ptr(iovlap) + (elem - ielem0) * jelems
             mat%zmtx(ptr:ptr+jelems-1) = data(jelem0:jelem0+jelems-1)
          end do

       ! Nothing to do if segment is blank

       end if

    end do  ! seg

  end subroutine sparse_put_col_complex

  !============================================================================!
  ! This subroutine clears an array according to the sparsity pattern of a     !
  ! column of a given block sparse matrix local to this processor.             !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   data (output)  : The real array to clear                                 !
  !   mat  (output)  : The real matrix whose sparsity pattern is to be used    !
  !   elem (input)   : The element index of the column to use                  !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, May 2009 based on SPAM2 version.                 !
  !----------------------------------------------------------------------------!
  ! SPAM2 version written by Peter Haynes, in May 2004.                        !
  !============================================================================!

  subroutine sparse_clr_col_complex(data,mat,elem)

    use comms, only: pub_total_num_procs, pub_my_proc_id
    use parallel_strategy, only: PARAL_INFO
    use utils, only: utils_assert

    implicit none

    ! Arguments
    complex(kind=DP), intent(inout) :: data(:) ! The column to insert
    type(SPAM3), intent(in) :: mat             ! The matrix
    integer, intent(in) :: elem                ! The elem index of the col

    ! Local variables
    integer :: ilib      ! Library entry for mat
    integer :: row_blks  ! Identifier for column blocking scheme
    integer :: col_blks  ! Identifier for row blocking scheme
    integer :: iblk      ! The block-column of the requested column
    integer :: iovlap    ! Block overlap counter
    integer :: loc_iblk  ! Block-column counter local to this proc
    integer :: jblk      ! Block-row counter
    integer :: jelems    ! Number of rows in block-row
    integer :: jelem0    ! First row in block-row
    integer :: seg       ! Segment index
    integer :: seg_start ! Start position of this segment in the index
    type(PARAL_INFO), pointer :: row_par ! Parallel strategy for rows
    type(PARAL_INFO), pointer :: col_par ! Parallel strategy for columns

    ! Get library entry
    ilib = mat%lib
    row_blks = library(ilib)%row_blks
    col_blks = library(ilib)%col_blks

    ! rc2013: assign the pointers for par based on library structure
    row_par => library(ilib)%row_par
    col_par => library(ilib)%col_par

    ! Check arguments
    call utils_assert(mat%iscmplx, &
         'Error in sparse_clr_col_complex: complex matrices only.')
    call utils_assert(size(data) >= library(ilib)%nrows, &
         'Error in sparse_clr_col_complex: insufficient space in input array.')
    call utils_assert(elem >= 1 .and. elem <= library(ilib)%mcols, &
         'Error in sparse_clr_col_complex: invalid column index.')

    ! Obtain information about the relevant block-column for this elem
    iblk = col_par%atom_of_elem(elem,col_blks)

    ! Check this block-column is local to this processor
    call utils_assert(iblk >= col_par%my_first_blk .and. iblk <= col_par%my_last_blk, 'Error &
         &in sparse_clr_col_complex: requested column is not local to proc ', &
         pub_my_proc_id)

    ! Find local index for the block-column
    loc_iblk = iblk - col_par%my_first_blk + 1

    do seg=0,pub_total_num_procs-1

       ! Nothing to do if segment is blank
       if (library(ilib)%seg_info(s_type,seg)==SEG_BLANK) cycle

       ! Loop over the nonzero blocks of this segment
       seg_start = library(ilib)%seg_info(s_idx,seg)
       do iovlap=library(ilib)%blk_idx(seg_start+loc_iblk-1), &
            library(ilib)%blk_idx(seg_start+loc_iblk)-1
          jblk = library(ilib)%blk_idx(iovlap)   ! block-row

          ! Set the nonzero elements of this col to zero
          jelems = row_par%num_elems_on_atom(jblk,row_blks)
          jelem0 = row_par%first_elem_on_atom(jblk,row_blks)
          data(jelem0:jelem0+jelems-1) = cmplx(0.0_DP,0.0_DP,kind=DP)

       end do  ! iovlap

    end do  ! seg

  end subroutine sparse_clr_col_complex

  !============================================================================!
  ! This subroutine is a wrapper to insert an element of type COEF into a      !
  ! sparse matrix                                                              !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   el (input)     : The element of the matrix to insert                     !
  !   mat (inout)    : The matrix (real or complex)                            !
  !   jrow (input)   : The element index of the row required                   !
  !   icol (input)   : The element index of the column required                !
  !----------------------------------------------------------------------------!
  ! Written by Andrea Greco, February 2016                                     !
  !============================================================================!

  subroutine sparse_put_element_coef(el, mat, jrow, icol)

     use datatypes, only: COEF

     implicit none

     ! Arguments
     type(COEF), intent(in)     :: el        ! The element to insert
     type(SPAM3), intent(inout) :: mat       ! The matrix
     integer, intent(in)        :: jrow      ! The index of the row
     integer, intent(in)        :: icol      ! The index of the column

     ! jme: No need to check that element and matrix are of the same type,
     ! since the check will be done inside sparse_put_element_real/complex.

     if (el%iscmplx) then
        call sparse_put_element_complex(el%z, mat, jrow, icol)
     ! real case
     else
        call sparse_put_element_real(el%d, mat, jrow, icol)
     end if

  end subroutine sparse_put_element_coef


  !============================================================================!
  ! This subroutine is a wrapper to extract an element of type COEF from a     !
  ! sparse matrix                                                              !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   el (input)     : The element of the matrix to extract                    !
  !   mat (inout)    : The matrix (real or complex)                            !
  !   jrow (input)   : The element index of the row required                   !
  !   icol (input)   : The element index of the column required                !
  !----------------------------------------------------------------------------!
  ! Written by Andrea Greco, February 2016                                     !
  !============================================================================!

  subroutine sparse_get_element_coef(el, mat, jrow, icol)

     use datatypes, only: COEF

     implicit none

     ! Arguments
     type(COEF), intent(inout) :: el               ! The element to extract
     type(SPAM3), intent(in) :: mat                ! The matrix
     integer, intent(in)     :: jrow               ! The index of the row
     integer, intent(in)     :: icol               ! The index of the column

     ! jme: No need to check that element and matrix are of the same type,
     ! since the check will be done inside sparse_get_element_real/complex.

     if (el%iscmplx) then
        call sparse_get_element_complex(el%z, mat, jrow, icol)
     ! real case
     else
        call sparse_get_element_real(el%d, mat, jrow, icol)
     end if

  end subroutine sparse_get_element_coef


  !============================================================================!
  ! Returns the 1-(induced) norm of a matrix.                                  !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   norm (output)  : The calculated 1-norm of the matrix                     !
  !   mat (input)    : The real matrix whose 1-norm is required                !
  !----------------------------------------------------------------------------!
  ! Written by Jolyon Aarons, March 2015 based on sparse_get_col_real.         !
  !----------------------------------------------------------------------------!
  ! Which was witten by Nicholas Hine, May 2009 based on SPAM2 version.        !
  ! SPAM2 version written by Peter Haynes, in May 2004.                        !
  !============================================================================!

  subroutine sparse_1norm(norm,mat)

    use comms, only: pub_my_proc_id, pub_total_num_procs, comms_reduce
    use parallel_strategy, only: PARAL_INFO
    use utils, only: utils_assert

    implicit none

    ! Arguments
    real(kind=dp), intent(out) :: norm         ! The 1-norm to compute
    type(SPAM3), intent(in) :: mat             ! The matrix

    ! Local variables
    integer :: ilib      ! Library entry for mat
    integer :: row_blks  ! Identifier for column blocking scheme
    integer :: col_blks  ! Identifier for row blocking scheme
    integer :: iblk      ! The block-column of the requested column
    integer :: ielem0    ! The first elem in the requested block-column
    integer :: iovlap    ! Block overlap counter
    integer :: loc_iblk  ! Block-column counter local to this proc
    integer :: jblk      ! Block-row counter
    integer :: jelems    ! Number of rows in block-row
    integer :: jelem0    ! First row in block-row
    integer :: ptr       ! Pointer to data array
    integer :: seg       ! Segment index
    integer :: seg_start ! Start position of this segment in the index
    real(kind=dp) :: locnorm   ! Local contribution to norm of matrix
    real(kind=dp) :: colsum    ! Sum of absolute elements of column
    integer       :: elem      ! Column index
    type(PARAL_INFO), pointer :: row_par
    type(PARAL_INFO), pointer :: col_par


    ! Get library entry
    ilib = mat%lib
    row_blks = library(ilib)%row_blks
    col_blks = library(ilib)%col_blks

    ! rc2013: assign the pointers for par based on library structure
    row_par => library(ilib)%row_par
    col_par => library(ilib)%col_par

    locnorm=0.0_dp
    norm=0.0_dp

    if(mat%iscmplx) then
       ! Loop over columns on this processor
       colloopz: do elem=col_par%first_elem_on_proc(pub_my_proc_id,col_blks),col_par%first_elem_on_proc(pub_my_proc_id+1,col_blks)-1
          colsum=0.0_dp
          ! Loop over proc-proc segments
          segloopz: do seg=0,pub_total_num_procs-1

             ! Segment is dense
             if (library(ilib)%seg_info(s_type,seg) == SEG_DENSE) then

                jelems = row_par%num_elems_on_proc(seg,row_blks)
                ptr = library(ilib)%seg_info(s_ptr,seg) + &
                     (elem - col_par%first_elem_on_proc(pub_my_proc_id,col_blks)) * &
                     row_par%num_elems_on_proc(seg,row_blks)

                colsum = colsum + sum(abs(mat%zmtx(ptr:ptr+jelems-1)))

                ! Segment is sparse
             else if (library(ilib)%seg_info(s_type,seg) == SEG_SPARSE) then

                ! Obtain information about the relevant block-column for this elem
                iblk = col_par%atom_of_elem(elem,col_blks)
                ielem0 = col_par%first_elem_on_atom(iblk,col_blks)

                ! Find local index for the block-column
                loc_iblk = iblk - col_par%my_first_blk + 1
                seg_start = library(ilib)%seg_info(s_idx,seg)
                do iovlap=library(ilib)%blk_idx(seg_start+loc_iblk-1), &
                     library(ilib)%blk_idx(seg_start+loc_iblk)-1
                   jblk = library(ilib)%blk_idx(iovlap)   ! block-row
                   jelems = row_par%num_elems_on_atom(jblk,row_blks)
                   jelem0 = row_par%first_elem_on_atom(jblk,row_blks)
                   ptr = library(ilib)%blk_ptr(iovlap) + (elem - ielem0) * jelems
                   colsum = colsum + sum(abs(mat%zmtx(ptr:ptr+jelems-1)))
                end do

             end if

          end do segloopz
          if(colsum>locnorm) locnorm=colsum
       end do colloopz
    else
       ! Loop over columns on this processor
       colloop: do elem=col_par%first_elem_on_proc(pub_my_proc_id,col_blks),col_par%first_elem_on_proc(pub_my_proc_id+1,col_blks)-1
          colsum=0.0_dp
          ! Loop over proc-proc segments
          segloop: do seg=0,pub_total_num_procs-1

             ! Segment is dense
             if (library(ilib)%seg_info(s_type,seg) == SEG_DENSE) then

                jelems = row_par%num_elems_on_proc(seg,row_blks)
                ptr = library(ilib)%seg_info(s_ptr,seg) + &
                     (elem - col_par%first_elem_on_proc(pub_my_proc_id,col_blks)) * &
                     row_par%num_elems_on_proc(seg,row_blks)

                colsum = colsum + sum(abs(mat%dmtx(ptr:ptr+jelems-1)))

                ! Segment is sparse
             else if (library(ilib)%seg_info(s_type,seg) == SEG_SPARSE) then

                ! Obtain information about the relevant block-column for this elem
                iblk = col_par%atom_of_elem(elem,col_blks)
                ielem0 = col_par%first_elem_on_atom(iblk,col_blks)

                ! Find local index for the block-column
                loc_iblk = iblk - col_par%my_first_blk + 1
                seg_start = library(ilib)%seg_info(s_idx,seg)
                do iovlap=library(ilib)%blk_idx(seg_start+loc_iblk-1), &
                     library(ilib)%blk_idx(seg_start+loc_iblk)-1
                   jblk = library(ilib)%blk_idx(iovlap)   ! block-row
                   jelems = row_par%num_elems_on_atom(jblk,row_blks)
                   jelem0 = row_par%first_elem_on_atom(jblk,row_blks)
                   ptr = library(ilib)%blk_ptr(iovlap) + (elem - ielem0) * jelems
                   colsum = colsum + sum(abs(mat%dmtx(ptr:ptr+jelems-1)))
                end do

             end if

          end do segloop
          if(colsum>locnorm) locnorm=colsum
       end do colloop
    end if

    call comms_reduce('MAX',locnorm)
    norm=locnorm

  end subroutine sparse_1norm

  !============================================================================!
  ! This function calculates the entrywise norm of the data in a given sparse  !
  ! matrix.                                                                    !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   mat   (input) : The sparse matrix to check                               !
  !   ord   (input) : order to raise entries to                                !
  !----------------------------------------------------------------------------!
  ! Written by Jolyon Aarons, June 2016.                                       !
  !============================================================================!

  function sparse_entrywise_norm(mat,ord) result(norm)

    use comms, only: comms_reduce, comms_bcast, pub_root_proc_id

    implicit none

    ! Arguments
    type(SPAM3), intent(in)     :: mat
    integer,     intent(in)     :: ord
    real(kind=DP)               :: norm

    norm = 0.0_dp
    if(mat%iscmplx) then
       norm=sum(abs(mat%zmtx)**ord)
    else
       norm=sum(abs(mat%dmtx)**ord)
    end if

    ! Sum over all procs
    call comms_reduce("SUM",norm)

    norm = norm**(1.0_dp/real(ord,dp))

    ! Synchronise to ensure all procs hold same data to avoid parallel desync
    call comms_bcast(pub_root_proc_id,norm)

  end function sparse_entrywise_norm

  !============================================================================!
  ! This function checks if any entries in a given sparse matrix are listed    !
  ! as "Not a Number" (NaN) for debugging purposes.                            !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   mat   (input) : The sparse matrix to check                               !
  !----------------------------------------------------------------------------!
  ! Written by Jolyon Aarons, June 2016.                                       !
  !============================================================================!

  function sparse_any_isnan(mat) result(isnan)

    use comms, only: comms_reduce
    use utils, only: utils_isnan

    implicit none

    ! Arguments
    type(SPAM3), intent(in) :: mat
    logical                 :: isnan

    ! Local Variables
    integer                 :: ii
    logical                 :: loc_isnan

    loc_isnan=.false.

    ! Check if local data contains NaN
    if (mat%iscmplx) then
       do ii=1,library(mat%lib)%my_nze
          loc_isnan=loc_isnan.or.utils_isnan(aimag(mat%zmtx(ii)))
          loc_isnan=loc_isnan.or.utils_isnan(real(mat%zmtx(ii),kind=DP))
       end do
    else
       do ii=1,library(mat%lib)%my_nze
          loc_isnan=loc_isnan.or.utils_isnan(mat%dmtx(ii))
       end do
    end if

    ! Check if any proc has found a NaN
    call comms_reduce('OR',loc_isnan)

    isnan = loc_isnan

  end function sparse_any_isnan


  !============================================================================!
  ! This subroutine converts a block sparse matrix in a SPAM3 into a dense     !
  ! full square real array, or collates a dense matrix across procs to one     !
  ! single full square array                                                   !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   dest   (output) : The destination full square matrix                     !
  !   src    (input)  : The source block sparse matrix                         !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, June 2009.                                       !
  !============================================================================!

  subroutine sparse_spam3tofull_complex(dest,src)

    use comms, only: comms_bcast, pub_total_num_procs
    use constants, only: cmplx_0
    use parallel_strategy, only: PARAL_INFO
    use utils, only: utils_assert

    implicit none

    ! Arguments
    complex(kind=DP), intent(out) :: dest(:,:)  ! The full square matrix
    type(SPAM3), intent(in) :: src              ! The block sparse matrix

    ! Local variables
    integer :: srclib        ! Library entry for src
    integer :: row_blks      ! Identifier for column blocking scheme
    integer :: col_blks      ! Identifier for row blocking scheme
    integer :: iblk,jblk     ! Block counters
    integer :: loc_iblk      ! Local counter for iblk
    integer :: idx           ! Index
    integer :: ptr           ! Pointer
    integer :: ielemonat     ! elem in atom block iblk counter
    integer :: jelemonat     ! elem in atom block jblk counter
    integer :: ielem,jelem   ! elem counters
    integer :: ielems,jelems ! Numbers of element rows/cols in this block
    integer :: nrows         ! Total number of rows in this block/segment
    integer :: proc          ! Proc loop counter
    integer :: seg           ! Segment index
    integer :: seg_type      ! Start position of this segment in the index
    integer :: seg_start     ! Start position of this segment in the index
    type(PARAL_INFO), pointer :: row_par
    type(PARAL_INFO), pointer :: col_par ! Parallel strategy for columns

    ! Get library entry for src
    srclib = src%lib
    row_blks = library(srclib)%row_blks
    col_blks = library(srclib)%col_blks

    ! rc2013: assign the pointers for par based on library structure
    row_par => library(srclib)%row_par
    col_par => library(srclib)%col_par

    ! Check arguments
    call utils_assert(size(dest,1) == library(srclib)%nrows, 'Error in sparse_&
         &spam3tofull_complex: incompatible numbers of rows in arguments.')
    call utils_assert(size(dest,2) == library(srclib)%mcols, 'Error in sparse_&
         &spam3tofull_complex: incompatible numbers of cols in arguments.')

    ! Zero full square matrix
    dest = cmplx_0
    nrows = 0

    ! Loop over the segments of src on this proc
    do seg=0,pub_total_num_procs-1

       ! Cycle if this segment is blank
       seg_type = library(srclib)%seg_info(s_type,seg)
       if (seg_type == SEG_BLANK) cycle
       if (seg_type == SEG_DENSE) nrows = row_par%num_elems_on_proc(seg,row_blks)

       ! Loop over block-columns of src on this proc
       seg_start = library(srclib)%seg_info(s_idx,seg)
       loc_iblk = 0
       do iblk=col_par%my_first_blk,col_par%my_last_blk
          loc_iblk = loc_iblk + 1
          ielem = col_par%first_elem_on_atom(iblk,col_blks)
          ielems = col_par%num_elems_on_atom(iblk,col_blks)

          ! Loop over block-rows jblk in column iblk of src
          do idx=library(srclib)%blk_idx(seg_start+loc_iblk-1), &
               library(srclib)%blk_idx(seg_start+loc_iblk)-1
             jblk = library(srclib)%blk_idx(idx)
             jelem = row_par%first_elem_on_atom(jblk,row_blks)
             jelems = row_par%num_elems_on_atom(jblk,row_blks)

             if (seg_type == SEG_SPARSE) nrows = jelems

             ! Copy block
             ptr = library(srclib)%blk_ptr(idx)
             if (src%iscmplx) then
                do ielemonat=0,ielems-1      ! Element col in block
                   do jelemonat=0,jelems-1   ! Element row in block
                      dest(jelem+jelemonat,ielem+ielemonat) = &
                           src%zmtx(ptr+ielemonat*nrows+jelemonat)
                   end do
                end do
             else
                do ielemonat=0,ielems-1      ! Element col in block
                   do jelemonat=0,jelems-1   ! Element row in block
                      dest(jelem+jelemonat,ielem+ielemonat) = &
                           cmplx(src%dmtx(ptr+ielemonat*nrows+jelemonat), &
                           0.0_DP,kind=DP)
                   end do
                end do
             end if

          end do  ! Loop over block-rows in column iblk of src

       end do  ! Loop over block-columns of src

    end do  ! seg

    ! Broadcast result across all procs
    do proc=0,pub_total_num_procs-1

       ! Get first column and number of columns
       ielem = col_par%first_elem_on_proc(proc,col_blks)
       ielems = col_par%num_elems_on_proc(proc,col_blks)

       ! Broadcast this part of the matrix
       call comms_bcast(proc,dest(1,ielem),ielems*library(srclib)%nrows)

    end do  ! Loop over procs

  end subroutine sparse_spam3tofull_complex

  !============================================================================!
  ! This subroutine converts a block sparse matrix in a SPAM3 into a dense     !
  ! full square real array, or collates a dense matrix across procs to one     !
  ! single full square array                                                   !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   dest   (output) : The destination full square matrix                     !
  !   src    (input)  : The source block sparse matrix                         !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, June 2009.                                       !
  !============================================================================!

  subroutine sparse_fulltospam3_complex(dest,src)

    use comms, only: pub_total_num_procs
    use parallel_strategy, only: PARAL_INFO
    use utils, only: utils_assert

    implicit none

    ! Arguments
    type(SPAM3), intent(inout) :: dest       ! The block sparse matrix
    complex(kind=DP), intent(in) :: src(:,:) ! The full square matrix

    ! Local variables
    integer :: destlib       ! Library entry for dest
    integer :: row_blks      ! Identifier for column blocking scheme
    integer :: col_blks      ! Identifier for row blocking scheme
    integer :: iblk,jblk     ! Block counters
    integer :: loc_iblk      ! Local counter for iblk
    integer :: idx           ! Index
    integer :: ptr           ! Pointer
    integer :: ielemonat     ! Element in atom block iblk counter
    integer :: jelemonat     ! Element in atom block jblk counter
    integer :: ielem,jelem   ! Element start positions of a block
    integer :: ielems,jelems ! Numbers of element rows/cols in this block
    integer :: nrows         ! Total number of rows in the block/segment
    integer :: seg           ! Segment index
    integer :: seg_type      ! Segment type for this segment
    integer :: seg_start     ! Start position of this segment in the index
    type(PARAL_INFO), pointer :: row_par ! Parallel strategy for columns
    type(PARAL_INFO), pointer :: col_par ! Parallel strategy for columns

    ! Get library entry for src
    destlib = dest%lib
    row_blks = library(destlib)%row_blks
    col_blks = library(destlib)%col_blks
    nrows = 0

    ! rc2013: assign the pointers for par based on library structure
    row_par => library(destlib)%row_par
    col_par => library(destlib)%col_par

    ! Check arguments
    call utils_assert(size(src,1) == library(destlib)%nrows, 'Error in sparse_&
         &fulltospam3_complex: incompatible numbers of rows in arguments.')
    call utils_assert(size(src,2) == library(destlib)%mcols, 'Error in sparse_&
         &fulltospam3_complex: incompatible numbers of cols in arguments.')
    call utils_assert(dest%iscmplx, &
         'Error in sparse_fulltospam3_complex: complex matrices only.')

    ! Loop over the segments of dest on this proc
    do seg=0,pub_total_num_procs-1

       ! Cycle if this segment is blank
       seg_type = library(destlib)%seg_info(s_type,seg)
       if (seg_type == SEG_BLANK) cycle
       if (seg_type == SEG_DENSE) nrows=row_par%num_elems_on_proc(seg,row_blks)

       ! Loop over block-columns of dest on this proc
       seg_start = library(destlib)%seg_info(s_idx,seg)
       loc_iblk = 0
       do iblk=col_par%my_first_blk,col_par%my_last_blk
          loc_iblk = loc_iblk + 1
          ielems = col_par%num_elems_on_atom(iblk,col_blks)

          ! Loop over block-rows jblk in column iblk of dest
          do idx=library(destlib)%blk_idx(seg_start+loc_iblk-1), &
               library(destlib)%blk_idx(seg_start+loc_iblk)-1
             jblk = library(destlib)%blk_idx(idx)
             jelems = row_par%num_elems_on_atom(jblk,row_blks)
             if (seg_type == SEG_SPARSE) nrows = jelems

             ! Copy block
             ielem = col_par%first_elem_on_atom(iblk,col_blks)
             ptr = library(destlib)%blk_ptr(idx)
             do ielemonat=0,ielems-1      ! Element col in block
                jelem = row_par%first_elem_on_atom(jblk,row_blks)
                do jelemonat=0,jelems-1   ! Element row in block
                   dest%zmtx(ptr+ielemonat*nrows+jelemonat) = &
                        src(jelem+jelemonat,ielem+ielemonat)
                end do
             end do

          end do  ! Loop over block-rows in column iblk of dest

       end do  ! Loop over block-columns of dest

    end do  ! seg

  end subroutine sparse_fulltospam3_complex

  !============================================================================!
  ! This subroutine scales and optionally shifts the eigenvalues of a matrix.  !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   mat    (inout) : The matrix to be rescaled and shifted                   !
  !   alpha  (input) : The complex scaling parameter                           !
  !   beta   (input) : The optional real shift parameter                       !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, May 2009 based on SPAM2 version.                 !
  !----------------------------------------------------------------------------!
  ! SPAM2 version written by Peter Haynes, March 2004.                         !
  ! Revised for distributed data, June 2006.                                   !
  ! Dense matrix version by Nicholas Hine, Dec 2007.                           !
  !============================================================================!

  subroutine sparse_scale_complex(mat,alpha,beta)

    use comms, only: pub_my_proc_id
    use parallel_strategy, only: PARAL_INFO
    use utils, only: utils_assert

    implicit none

    ! Arguments
    type(SPAM3), intent(inout)          :: mat    ! The matrix to be operated on
    complex(kind=DP), intent(in)        :: alpha  ! The scaling parameter
    real(kind=DP), optional, intent(in) :: beta   ! The shift parameter

    ! Local variables
    integer :: ilib      ! Library entry for mat
    integer :: blks      ! Identifier for blocking scheme
    integer :: iblk      ! Block-column loop counter
    integer :: loc_iblk  ! Loop counter for iblk on this proc
    integer :: ielems    ! Number of columns in block
    integer :: ptr       ! Pointer to diagonal entries
    integer :: ielem     ! Element in block-column loop counte
    integer :: seg       ! Segment index
    integer :: seg_start ! Segment type for this segment
    type(PARAL_INFO), pointer :: row_par ! Parallel strategy for rows
    type(PARAL_INFO), pointer :: col_par ! Parallel strategy for columns


    ! This only makes sense if mat is complex...
    call utils_assert(mat%iscmplx, &
         'Error in sparse_scale_complex: mat must be complex.')

    ! Rescale the elements if required
    if (alpha /= (1.0_DP,0.0_DP)) then
       if (mat%iscmplx) then
          mat%zmtx = alpha * mat%zmtx
       else
          mat%dmtx = real(alpha * mat%dmtx,kind=DP)
       end if
    end if

    ! Shift the eigenvalues if required
    if (present(beta)) then

       ! Get library entry for mat
       ilib = mat%lib
       blks = library(ilib)%row_blks
       seg = pub_my_proc_id

       ! rc2013: assign the pointers for par based on library structure
       row_par => library(ilib)%row_par
       col_par => library(ilib)%col_par

       ! Can only shift eigenvalues of square matrix
       call utils_assert(library(ilib)%row_blks==library(ilib)%col_blks, &
            'Error in sparse_scale: cannot shift eigenvalues of non-square &
            &matrices.')

       if (library(ilib)%seg_info(s_type,seg)==SEG_DENSE) then

          ! Set pointer to first diagonal element
          ptr = library(ilib)%seg_info(s_ptr,seg)

          ! Get number of columns in this segment
          ielems = col_par%num_elems_on_proc(seg,blks)

          ! Add identity matrix scaled by beta to all elements
          if (mat%iscmplx) then
             ! Loop over columns on my proc
             do ielem=1,ielems
                mat%zmtx(ptr) = mat%zmtx(ptr) + cmplx(beta,0.0_DP,kind=DP)
                ptr = ptr + ielems + 1
             end do
          else
             ! Loop over columns on my proc
             do ielem=1,col_par%num_elems_on_proc(seg,blks)
                mat%dmtx(ptr) = mat%dmtx(ptr) + beta
                ptr = ptr + ielems + 1
             end do
          end if

       else if (library(ilib)%seg_info(s_type,seg)==SEG_SPARSE) then

          seg_start = library(ilib)%seg_info(s_idx,seg)

          ! Loop over atom-blocks iblk on my proc
          loc_iblk = 0
          do iblk=col_par%my_first_blk,col_par%my_last_blk
             loc_iblk = loc_iblk + 1

             ! Get number of elements on this atom
             ielems = col_par%num_elems_on_atom(iblk,blks)

             ! Find pointer to start of diagonal block in this column
             ptr = library(ilib)%blk_ptr(seg_start+loc_iblk-1)

             ! Loop over diagonal elements in this diagonal block and add
             ! identity matrix scaled by beta
             if (mat%iscmplx) then
                do ielem=1,ielems
                   mat%zmtx(ptr) = mat%zmtx(ptr) + cmplx(beta,0.0_DP,kind=DP)
                   ptr = ptr + ielems + 1
                end do
             else
                do ielem=1,ielems
                   mat%dmtx(ptr) = mat%dmtx(ptr) + beta
                   ptr = ptr + ielems + 1
                end do
             end if

          end do  ! Loop over atom-blocks iblk on this proc

       endif

    end if


  end subroutine sparse_scale_complex

  !============================================================================!
  ! This subroutine scales and optionally shifts the eigenvalues of a matrix.  !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   mat    (inout) : The matrix to be rescaled and shifted                   !
  !   alpha  (input) : The complex scaling parameter                           !
  !   beta   (input) : The optional complex shift parameter                    !
  !----------------------------------------------------------------------------!
  ! Modified from sparse_scale_complex by Jolyon Aarons, May 2017.             !
  !----------------------------------------------------------------------------!
  ! SPAM2 version written by Peter Haynes, March 2004.                         !
  ! Revised for distributed data, June 2006.                                   !
  ! Dense matrix version by Nicholas Hine, Dec 2007.                           !
  !============================================================================!

  subroutine sparse_scale_complex2(mat,alpha,beta)
    use constants, only: cmplx_1
    use comms, only: pub_my_proc_id
    use parallel_strategy, only: PARAL_INFO
    use utils, only: utils_assert

    implicit none

    ! Arguments
    type(SPAM3), intent(inout)          :: mat    ! The matrix to be operated on
    complex(kind=DP), intent(in)        :: alpha  ! The scaling parameter
    complex(kind=DP), intent(in) :: beta   ! The shift parameter

    ! Local variables
    integer :: ilib      ! Library entry for mat
    integer :: blks      ! Identifier for blocking scheme
    integer :: iblk      ! Block-column loop counter
    integer :: loc_iblk  ! Loop counter for iblk on this proc
    integer :: ielems    ! Number of columns in block
    integer :: ptr       ! Pointer to diagonal entries
    integer :: ielem     ! Element in block-column loop counte
    integer :: seg       ! Segment index
    integer :: seg_start ! Segment type for this segment
    type(PARAL_INFO), pointer :: col_par ! Parallel strategy for columns


    ! This only makes sense if mat is complex...
    call utils_assert(mat%iscmplx, &
         'Error in sparse_scale_complex2: mat must be complex.')

    ! Rescale the elements if required
    if (alpha /= cmplx_1) then
       if (mat%iscmplx) then
          mat%zmtx = alpha * mat%zmtx
       else
          mat%dmtx = real(alpha * mat%dmtx,kind=DP)
       end if
    end if

    ! Shift the eigenvalues if required
!    if (present(beta)) then

       ! Get library entry for mat
       ilib = mat%lib
       blks = library(ilib)%row_blks
       seg = pub_my_proc_id

       ! rc2013: assign the pointers for par based on library structure
       col_par => library(ilib)%col_par

       ! Can only shift eigenvalues of square matrix
       call utils_assert(library(ilib)%row_blks==library(ilib)%col_blks, &
            'Error in sparse_scale: cannot shift eigenvalues of non-square &
            &matrices.')

       if (library(ilib)%seg_info(s_type,seg)==SEG_DENSE) then

          ! Set pointer to first diagonal element
          ptr = library(ilib)%seg_info(s_ptr,seg)

          ! Get number of columns in this segment
          ielems = col_par%num_elems_on_proc(seg,blks)

          ! Add identity matrix scaled by beta to all elements
          if (mat%iscmplx) then
             ! Loop over columns on my proc
             do ielem=1,ielems
                mat%zmtx(ptr) = mat%zmtx(ptr) + beta
                ptr = ptr + ielems + 1
             end do
          else
             ! Loop over columns on my proc
             do ielem=1,col_par%num_elems_on_proc(seg,blks)
                mat%dmtx(ptr) = mat%dmtx(ptr) + real(beta,dp)
                ptr = ptr + ielems + 1
             end do
          end if

       else if (library(ilib)%seg_info(s_type,seg)==SEG_SPARSE) then

          seg_start = library(ilib)%seg_info(s_idx,seg)

          ! Loop over atom-blocks iblk on my proc
          loc_iblk = 0
          do iblk=col_par%my_first_blk,col_par%my_last_blk
             loc_iblk = loc_iblk + 1

             ! Get number of elements on this atom
             ielems = col_par%num_elems_on_atom(iblk,blks)

             ! Find pointer to start of diagonal block in this column
             ptr = library(ilib)%blk_ptr(seg_start+loc_iblk-1)

             ! Loop over diagonal elements in this diagonal block and add
             ! identity matrix scaled by beta
             if (mat%iscmplx) then
                do ielem=1,ielems
                   mat%zmtx(ptr) = mat%zmtx(ptr) + beta
                   ptr = ptr + ielems + 1
                end do
             else
                do ielem=1,ielems
                   mat%dmtx(ptr) = mat%dmtx(ptr) + real(beta,dp)
                   ptr = ptr + ielems + 1
                end do
             end if

          end do  ! Loop over atom-blocks iblk on this proc

       endif

!    end if


  end subroutine sparse_scale_complex2

  !============================================================================!
  ! This subroutine transforms a complex matrix into its conjugate.            !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   amat    (inout) : The matrix to be complex conjugated                    !
  !   bmat    (inout) : if present, stores the conjugate of amat               !
  !----------------------------------------------------------------------------!
  ! Written by Andrea Greco, March 2016.                                       !
  !----------------------------------------------------------------------------!
  !============================================================================!

  subroutine sparse_conjugate(amat, bmat)

    use utils, only: utils_assert

    implicit none

    ! Arguments
    type(SPAM3), intent(inout)           :: amat    ! The matrix to be operated on
    type(SPAM3), optional, intent(inout) :: bmat    ! Conjugate stored here if present

    ! This only makes sense if matrices are complex...
    call utils_assert(amat%iscmplx, &
         'Error in sparse_conjugate: amat must be complex.')

    ! store the conjugate of amat in bmat
    if (present(bmat)) then
       call utils_assert(bmat%iscmplx, &
            'Error in sparse_conjugate: bmat must be complex.')
       ! rc2013: check for consistency of parallel strategies
       call sparse_par_check(amat%lib, bmat%lib, routine='sparse_conjugate')
       bmat%zmtx = conjg(amat%zmtx)
    ! modify amat directly
    else
       amat%zmtx = conjg(amat%zmtx)
    end if

  end subroutine sparse_conjugate


  !============================================================================!
  ! This subroutine performs an axpy operation on the matrices:                !
  !   y := y + alpha * x                                                       !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   ymat  (inout) : The sparse matrix y                                      !
  !   xmat  (input) : The sparse matrix x                                      !
  !   alpha (input) : The complex parameter alpha                              !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, May 2009 based on SPAM2 version.                 !
  !----------------------------------------------------------------------------!
  ! SPAM2 version written by Peter Haynes, March 2004.                         !
  ! Revised for distributed data, June 2006.                                   !
  ! Dense matrix version by Nicholas Hine, Dec 2007.                           !
  ! Modified by Andrea Greco to protect safe conversion from complex xmat      !
  ! to real ymat, September 2016.                                              !
  !============================================================================!

  subroutine sparse_axpy_complex(ymat,xmat,alpha,cmplx_to_real)

    use comms, only: pub_total_num_procs
    use parallel_strategy, only: PARAL_INFO
!$  use rundat, only: pub_threads_max
    use utils, only: utils_assert, utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(SPAM3), intent(inout) :: ymat    ! The sparse matrix y
    type(SPAM3), intent(in) :: xmat       ! The sparse matrix x
    complex(kind=DP), intent(in) :: alpha ! The parameter alpha
    ! agrecocmplx
    logical, optional, intent(in) :: cmplx_to_real

    ! Local variables
    integer :: xlib,ylib    ! Library pointers for x and y
    integer :: row_blks     ! Identifier for column blocking scheme
    integer :: col_blks     ! Identifier for row blocking scheme
    integer :: iblk,jblk    ! Atom loop counters
    integer :: icol         ! Atom loop counters
    integer :: loc_iblk     ! Atom loop counter on this proc
    integer :: idx          ! Index loop counter
    integer :: iel          ! Element loop counter
    integer :: ifound       ! Loop counter over columns found
    integer :: nfound       ! Number of columns found
    integer :: xstart,xend  ! Start and end of block data of xmat
    integer :: ystart,yend  ! Start and end of block data of ymat
    integer :: ystride      ! Stride to get to next column of src
    integer :: xstride      ! Stride to get to next column of dest
    integer :: seg          ! Segment index
    integer :: x_seg_type   ! Segment type for segments of x
    integer :: y_seg_type   ! Segment type for segments of y
    integer :: x_seg_start  ! Start position of segments of x in the index
    integer :: y_seg_start  ! Start position of segments of y in the index
    integer :: ierr         ! Error flag
    integer, allocatable :: tfound_idx(:), tfound_ptr(:)
    logical, allocatable :: tfound_blk(:)
    logical :: loc_cmplx_to_real         ! agrecocmplx
    type(PARAL_INFO), pointer :: row_par ! Parallel strategy for rows
    type(PARAL_INFO), pointer :: col_par ! Parallel strategy for columns

    ! agrecocmplx
    loc_cmplx_to_real = .false.
    if (present(cmplx_to_real)) loc_cmplx_to_real = cmplx_to_real

    ! agrecocmplx: protection against unsafe conversion from complex to real
    call utils_assert( ((aimag(alpha)==0.0_DP) .and. (.not.(xmat%iscmplx))) &
         .or. (ymat%iscmplx) .or. loc_cmplx_to_real, &
         'ERROR: unsafe conversion of complex matrix to real. Call routine &
         &sparse_axpy_complex with cmplx_to_real=.true. to force the conversion.')

    ! Extract library entries
    xlib = xmat%lib
    ylib = ymat%lib
    row_blks = library(xlib)%row_blks
    col_blks = library(xlib)%col_blks
    xstride = 0; ystride = 0

    ! rc2013: assign the pointers for par based on library structure
    ! rc2013: ax+y only possible if x and y have same parallel strategies
    row_par => library(xlib)%row_par
    col_par => library(xlib)%col_par

    ! Check matrices share same blocking scheme
    call utils_assert ( (library(xlib)%row_blks == row_blks) .and. &
         (library(xlib)%col_blks == col_blks), 'Error in sparse_axpy_complex: &
         &x and y matrices have different blocking schemes.')

    ! rc2013: check for consistency of parallel strategies
    call sparse_par_check(ylib, xlib, routine='sparse_axpy')

    ! Check for alpha = 0 (no operation required)
    if (alpha == (0.0_DP,0.0_DP)) then
       return
    end if

    ! Check whether the matrices share the same structure
    if (xlib == ylib) then

       ! Trivial axpy of whole data arrays
       if (xmat%iscmplx) then
          if (ymat%iscmplx) then
             ymat%zmtx = ymat%zmtx + alpha * xmat%zmtx
          else
             do iel=1,size(ymat%dmtx)
                ! jme: bug (incorrect algebra) fixed on 11/05/2017.
                ymat%dmtx(iel) = ymat%dmtx(iel) + &
                     real(alpha * xmat%zmtx(iel), kind=DP)
             end do
          end if
       else
          if (ymat%iscmplx) then
             do iel=1,size(ymat%zmtx)
                ymat%zmtx(iel) = ymat%zmtx(iel) + &
                     alpha * cmplx(xmat%dmtx(iel),0.0_DP,kind=DP)
             end do
          else
             ymat%dmtx = ymat%dmtx + real(alpha,kind=DP) * xmat%dmtx
          end if
       end if

    ! Structures differ, so proceed segment-by-segment
    else

!$OMP PARALLEL NUM_THREADS(pub_threads_max) DEFAULT(none) &
!$OMP PRIVATE(seg,x_seg_type,y_seg_type,loc_iblk,iblk,idx,jblk,nfound, &
!$OMP      ifound,x_seg_start,y_seg_start,xstart,xend,ystart,yend,xstride, &
!$OMP      ystride,icol,iel,tfound_blk,tfound_idx,tfound_ptr,ierr) &
!$OMP SHARED(xmat,ymat,xlib,ylib,library,row_blks,col_blks,alpha, &
!$OMP      col_par,row_par,pub_total_num_procs,pub_threads_max)

       allocate(tfound_idx(library(xlib)%nblk),stat=ierr)
       call utils_alloc_check('sparse_axpy_complex','tfound_idx',ierr)
       allocate(tfound_ptr(library(xlib)%nblk),stat=ierr)
       call utils_alloc_check('sparse_axpy_complex','tfound_ptr',ierr)
       allocate(tfound_blk(library(xlib)%nblk),stat=ierr)
       call utils_alloc_check('sparse_axpy_complex','tfound_blk',ierr)
       tfound_blk = .false.

!$OMP DO
       do seg=0,pub_total_num_procs-1

          x_seg_type = library(xlib)%seg_info(s_type,seg)
          y_seg_type = library(ylib)%seg_info(s_type,seg)

          ! Nothing to add to if this segment in y is blank
          if (y_seg_type==SEG_BLANK) cycle

          ! Nothing to add if this segment in x is blank
          if (x_seg_type == SEG_BLANK) cycle

          ! Find information about these segments
          x_seg_start = library(xlib)%seg_info(s_idx,seg)
          y_seg_start = library(ylib)%seg_info(s_idx,seg)
          if (y_seg_type == SEG_DENSE) ystride = &
               row_par%num_elems_on_proc(seg,row_blks)
          if (x_seg_type == SEG_DENSE) xstride = &
               row_par%num_elems_on_proc(seg,row_blks)

          ! Loop over atom blocks (block-columns) on this proc
          loc_iblk = 0
          do iblk=col_par%my_first_blk,col_par%my_last_blk
             loc_iblk = loc_iblk + 1

             ! Reset counters of number of rows found
             nfound = 0

             ! Loop over block-rows in x for this segment
             do idx=library(xlib)%blk_idx(x_seg_start+loc_iblk-1), &
                  library(xlib)%blk_idx(x_seg_start+loc_iblk)-1
                jblk = library(xlib)%blk_idx(idx)

                ! Mark flag for this row and add to list
                nfound = nfound + 1
                tfound_blk(jblk) = .true.
                tfound_idx(nfound) = jblk
                tfound_ptr(jblk) = library(xlib)%blk_ptr(idx)

             end do

             ! Loop over block-rows in y
             do idx=library(ylib)%blk_idx(y_seg_start+loc_iblk-1), &
                  library(ylib)%blk_idx(y_seg_start+loc_iblk)-1
                jblk = library(ylib)%blk_idx(idx)

                ! Check whether this block-row exists in source
                if (tfound_blk(jblk)) then

                   ! Find pointers to start of dest and src data
                   ystart = library(ylib)%blk_ptr(idx)
                   yend = ystart + row_par%num_elems_on_atom(jblk,row_blks) - 1
                   xstart = tfound_ptr(jblk)
                   xend = xstart + row_par%num_elems_on_atom(jblk,row_blks) - 1

                   ! Find stride to next column of dest data
                   if (y_seg_type == SEG_SPARSE) then
                      ystride = row_par%num_elems_on_atom(jblk,row_blks)
                   end if

                   ! Find stride to next column of source data
                   if (x_seg_type == SEG_SPARSE) then
                      xstride = row_par%num_elems_on_atom(jblk,row_blks)
                   end if

                   ! Copy columns from source to destination
                   do icol=1,col_par%num_elems_on_atom(iblk,col_blks)

                      if (xmat%iscmplx) then
                         if (ymat%iscmplx) then
                            ymat%zmtx(ystart:yend) = &
                                 ymat%zmtx(ystart:yend) + &
                                 alpha * xmat%zmtx(xstart:xend)
                         else
                            do iel=ystart,yend
                               ymat%dmtx(iel) = ymat%dmtx(iel) + real( &
                                    alpha * xmat%zmtx(iel-ystart+xstart), &
                                    kind=DP)
                            end do
                         end if
                      else
                         if (ymat%iscmplx) then
                            do iel=ystart,yend
                               ymat%zmtx(iel) = ymat%zmtx(iel) + alpha * &
                                    cmplx(xmat%dmtx(iel-ystart+xstart), &
                                    0.0_DP,kind=DP)
                            end do
                         else
                            ymat%dmtx(ystart:yend) = &
                                 ymat%dmtx(ystart:yend) + &
                                 real(alpha,kind=DP) * xmat%dmtx(xstart:xend)
                         end if
                      end if

                      ystart = ystart + ystride
                      yend = yend + ystride
                      xstart = xstart + xstride
                      xend = xend + xstride

                   end do

                end if

             end do  ! Loop over destination block-rows

             ! Reset flags
             do ifound=1,nfound
                jblk = tfound_idx(ifound)
                tfound_blk(jblk) = .false.
             end do

          end do  ! Loop over atom block-columns

       end do  ! seg
!$OMP END DO

      deallocate(tfound_blk,stat=ierr)
      call utils_dealloc_check('sparse_axpy_complex','tfound_blk',ierr)
      deallocate(tfound_ptr,stat=ierr)
      call utils_dealloc_check('sparse_axpy_complex','tfound_ptr',ierr)
      deallocate(tfound_idx,stat=ierr)
      call utils_dealloc_check('sparse_axpy_complex','tfound_idx',ierr)

!$OMP END PARALLEL

    end if


  end subroutine sparse_axpy_complex


  !============================================================================!
  ! This subroutine copies the data from one sparse matrix to another, taking  !
  ! different structures into account if necessary.                            !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   dest (inout) : The destination sparse matrix                             !
  !   src  (input) : The source sparse matrix                                  !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, May 2009 based on SPAM2 version.                 !
  !----------------------------------------------------------------------------!
  ! SPAM2 version written by Peter Haynes, March 2004.                         !
  ! Revised for distributed data, June 2006.                                   !
  ! Minor Modification for dense matrices by Nicholas Hine, Dec 2007.          !
  ! Added argument to check complex2real conversion by Andrea Greco, Sept 2016.!
  !============================================================================!

  subroutine sparse_copy(dest,src,cmplx_to_real)

    use comms, only: pub_total_num_procs
    use parallel_strategy, only: PARAL_INFO
!$  use rundat, only: pub_threads_max
    use utils, only: utils_assert, utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(SPAM3), intent(inout) :: dest   ! The destination sparse matrix
    type(SPAM3), intent(in) :: src       ! The source sparse matrix
    ! agrecocmplx
    logical, optional, intent(in) :: cmplx_to_real ! check complex2real

    ! Local variables
    integer :: destlib        ! Library entry for dest
    integer :: srclib         ! Library entry for src
    integer :: row_blks       ! Identifier for column blocking scheme
    integer :: col_blks       ! Identifier for row blocking scheme
    integer :: iblk,jblk      ! Atom loop counters
    integer :: loc_iblk       ! Atom loop counter on this proc
    integer :: idx            ! Index loop counter
    integer :: iel            ! Element loop counter
    integer :: ifound         ! Loop counter over columns found
    integer :: nfound         ! Number of columns found
    integer :: icol           ! Column element loop counter
    integer :: dstart,dend    ! Start and end of block data of dest
    integer :: sstart,send    ! Start and end of block data of src
    integer :: dstride        ! Stride to get to next column of src
    integer :: sstride        ! Stride to get to next column of dest
    integer :: seg            ! Segment index
    integer :: src_seg_type   ! Segment type for segments of src
    integer :: dest_seg_type  ! Segment type for segments of dest
    integer :: src_seg_start  ! Start position of segments of src in the index
    integer :: dest_seg_start ! Start position of segments of dest in the index
    integer :: ierr           ! Error flag
    integer, allocatable :: tfound_idx(:), tfound_ptr(:)
    logical, allocatable :: tfound_blk(:)
    type(PARAL_INFO), pointer :: row_par ! Parallel strategy for rows
    type(PARAL_INFO), pointer :: col_par ! Parallel strategy for columns
    logical :: loc_cmplx_to_real         ! agrecocmplx

    loc_cmplx_to_real = .false.
    if (present(cmplx_to_real)) loc_cmplx_to_real = cmplx_to_real

    ! agrecocmplx: protection  against unsafe conversion from complex to real
    call utils_assert(.not.(src%iscmplx).or.(dest%iscmplx).or.loc_cmplx_to_real, &
         'ERROR: unsafe conversion of complex matrix to real. Call routine &
         &sparse_copy with cmplx_to_real=.true. to force the conversion.')

    ! Extract library entries
    destlib = dest%lib
    srclib = src%lib
    row_blks = library(srclib)%row_blks
    col_blks = library(srclib)%col_blks

    ! rc2013: assign the pointers for par based on library structure
    ! rc2013: src and dest must have the same parallel strategies
    row_par => library(destlib)%row_par
    col_par => library(destlib)%col_par

    ! Check matrices share same blocking scheme
    call utils_assert( (library(destlib)%row_blks == row_blks) .and. &
         (library(destlib)%col_blks == col_blks), 'Error in sparse_copy: &
         &src and dest matrices have different blocking schemes.')

    ! rc2013: check for consistency of parallel strategies
    call sparse_par_check(destlib, srclib, routine='sparse_copy')

    ! Check whether the matrices share the same structure
    if ( srclib == destlib ) then

       ! Trivial copy of data arrays
       if (src%iscmplx) then
          if (dest%iscmplx) then
             dest%zmtx = src%zmtx
          else
!$OMP PARALLEL DO NUM_THREADS(pub_threads_max) &
!$OMP PRIVATE(iel) SHARED(src,dest,pub_threads_max)
             do iel=1,size(dest%dmtx)
                dest%dmtx(iel) = real(src%zmtx(iel),kind=DP)
             end do
!$OMP END PARALLEL DO
          end if
       else
          if (dest%iscmplx) then
!$OMP PARALLEL DO NUM_THREADS(pub_threads_max) &
!$OMP PRIVATE(iel) SHARED(src,dest,pub_threads_max)
             do iel=1,size(dest%zmtx)
                dest%zmtx(iel) = cmplx(src%dmtx(iel),0.0_DP,kind=DP)
             end do
!$OMP END PARALLEL DO
          else
!             dest%dmtx = src%dmtx
             do iel=1,size(dest%dmtx) ! ndmh and jhl52
                dest%dmtx(iel) = src%dmtx(iel)
             end do
          end if
       end if

    else

!$OMP PARALLEL NUM_THREADS(pub_threads_max) DEFAULT(none) &
!$OMP PRIVATE(seg,src_seg_type,dest_seg_type,sstart,send,dstart, &
!$OMP      dend,src_seg_start,dest_seg_start,loc_iblk,iblk,nfound,idx,jblk, &
!$OMP      sstride,dstride,icol,ifound,tfound_blk,tfound_idx,tfound_ptr,ierr) &
!$OMP SHARED(src,dest,srclib,destlib,pub_total_num_procs,row_blks,col_blks, &
!$OMP      col_par,row_par, &
!num_elems_on_atom,num_elems_on_proc, &
!$OMP      library,pub_threads_max)

       ! rc2013: allocate with number of block-rows
       allocate(tfound_idx(library(srclib)%nblk),stat=ierr)
       call utils_alloc_check('sparse_copy','tfound_idx',ierr)
       allocate(tfound_ptr(library(srclib)%nblk),stat=ierr)
       call utils_alloc_check('sparse_copy','tfound_ptr',ierr)
       allocate(tfound_blk(library(srclib)%nblk),stat=ierr)
       call utils_alloc_check('sparse_copy','tfound_blk',ierr)
       tfound_blk = .false.

!$OMP DO
       do seg=0,pub_total_num_procs-1

          src_seg_type = library(srclib)%seg_info(s_type,seg)
          dest_seg_type = library(destlib)%seg_info(s_type,seg)
          sstart = library(srclib)%seg_info(s_ptr,seg)
          send = library(srclib)%seg_info(s_ptr,seg+1) - 1
          dstart = library(destlib)%seg_info(s_ptr,seg)
          dend = library(destlib)%seg_info(s_ptr,seg+1) - 1

          ! Nothing to copy if destination segment is blank
          if (dest_seg_type==SEG_BLANK) cycle

          ! Fill dest data with zeros if source segment is blank
          if (src_seg_type == SEG_BLANK) then

             if (dest%iscmplx) then
                dest%zmtx(dstart:dend) = cmplx(0.0_DP,0.0_DP,kind=DP)
             else
                dest%dmtx(dstart:dend) = 0.0_DP
             end if

             cycle

          end if

          src_seg_start = library(srclib)%seg_info(s_idx,seg)
          dest_seg_start = library(destlib)%seg_info(s_idx,seg)

          ! Loop over block-columns on this proc
          loc_iblk = 0
          do iblk=col_par%my_first_blk,col_par%my_last_blk
             loc_iblk = loc_iblk + 1

             ! Reset counters of number of rows found
             nfound = 0

             ! Loop over block-rows in source for this segment
             do idx=library(srclib)%blk_idx(src_seg_start+loc_iblk-1), &
                  library(srclib)%blk_idx(src_seg_start+loc_iblk)-1
                jblk = library(srclib)%blk_idx(idx)

                ! Mark flag for this row and add to list
                nfound = nfound + 1
                tfound_blk(jblk) = .true.
                tfound_idx(nfound) = jblk
                tfound_ptr(jblk) = library(srclib)%blk_ptr(idx)

             end do  ! idx

             ! Loop over block-rows in destination
             do idx=library(destlib)%blk_idx(dest_seg_start+loc_iblk-1), &
                  library(destlib)%blk_idx(dest_seg_start+loc_iblk)-1
                jblk = library(destlib)%blk_idx(idx)

                ! Check whether this block-row exists in source
                if (tfound_blk(jblk)) then

                   ! Find pointers to start of dest and src data
                   dstart = library(destlib)%blk_ptr(idx)
                   dend = dstart + row_par%num_elems_on_atom(jblk,row_blks) - 1
                   sstart = tfound_ptr(jblk)
                   send = sstart + row_par%num_elems_on_atom(jblk,row_blks) - 1

                   ! Find stride to next column of dest data
                   if (dest_seg_type == SEG_DENSE) then
                      dstride = row_par%num_elems_on_proc(seg,row_blks)
                   else
                      dstride = row_par%num_elems_on_atom(jblk,row_blks)
                   end if

                   ! Find stride to next column of source data
                   if (src_seg_type == SEG_DENSE) then
                      sstride = row_par%num_elems_on_proc(seg,row_blks)
                   else
                      sstride = row_par%num_elems_on_atom(jblk,row_blks)
                   end if

                   ! Copy columns from source to destination
                   do icol=1,col_par%num_elems_on_atom(iblk,col_blks)

                      if (src%iscmplx) then
                         if (dest%iscmplx) then
                            dest%zmtx(dstart:dend) = src%zmtx(sstart:send)
                         else
                            do iel=dstart,dend
                               dest%dmtx(iel) = &
                                    real(src%zmtx(iel-dstart+sstart),kind=DP)
                            end do
                         end if
                      else
                         if (dest%iscmplx) then
                            do iel=dstart,dend
                               dest%zmtx(iel) = cmplx(src%dmtx(iel-&
                                    dstart+sstart),0.0_DP,kind=DP)
                            end do
                         else
                            dest%dmtx(dstart:dend) = src%dmtx(sstart:send)
                         end if
                      end if

                      dstart = dstart + dstride
                      dend = dend + dstride
                      sstart = sstart + sstride
                      send = send + sstride

                   end do

                else

                   ! Set data for this block to zero
                   dstart = library(destlib)%blk_ptr(idx)
                   dend = dstart + row_par%num_elems_on_atom(jblk,row_blks) - 1

                   if (dest_seg_type == SEG_DENSE) then
                      dstride = row_par%num_elems_on_proc(seg,row_blks)
                   else
                      dstride = row_par%num_elems_on_atom(jblk,row_blks)
                   end if

                   do icol=1,col_par%num_elems_on_atom(iblk,col_blks)

                      if (dest%iscmplx) then
                         dest%zmtx(dstart:dend) = (0.0_DP,0.0_DP)
                      else
                         dest%dmtx(dstart:dend) = 0.0_DP
                      end if

                      dstart = dstart + dstride
                      dend = dend + dstride

                   end do  ! icol

                end if

             end do  ! Loop over destination block-rows

             ! Reset flags
             do ifound=1,nfound
                jblk = tfound_idx(ifound)
                tfound_blk(jblk) = .false.
             end do

          end do  ! Loop over atom block-columns

       end do  ! seg
!$OMP END DO

      deallocate(tfound_blk,stat=ierr)
      call utils_dealloc_check('sparse_copy','tfound_blk',ierr)
      deallocate(tfound_ptr,stat=ierr)
      call utils_dealloc_check('sparse_copy','tfound_ptr',ierr)
      deallocate(tfound_idx,stat=ierr)
      call utils_dealloc_check('sparse_copy','tfound_idx',ierr)

!$OMP END PARALLEL

    end if

  end subroutine sparse_copy

  !============================================================================!
  ! This function checks if a complex sparse matrix is hermitian, by computing !
  ! the RMS of the difference mat_dagger-mat                                   !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   mat   (input)  : The sparse complex matrix to check                      !
  !----------------------------------------------------------------------------!
  ! Written by Andrea Greco, October 2016                                      !
  !============================================================================!

  real(kind=DP) function sparse_check_hermitian(mat)

    use utils, only: utils_assert

    implicit none

    ! Arguments
    type(SPAM3), intent(in)       :: mat           ! input complex matrix

    ! Local variables
    type(SPAM3)                   :: mat_buffer    ! mat_dagger - mat

    ! agrecocmplx: use only with complex matrices
    call utils_assert(mat%iscmplx, 'ERROR in routine sparse_check_hermitian: &
         &matrix must be complex')

    ! agrecocmplx: compute transpose conjugate of mat
    call sparse_create(mat_buffer, mat)
    call sparse_transpose(mat_buffer, mat)

    ! agrecocmplx: compute mat_dagger - mat
    call sparse_axpy(mat_buffer, mat, -1.0_DP)

    ! agrecocmplx: compute RMS of mat_dagger - mat
    sparse_check_hermitian = sparse_rms_element(mat_buffer)

    call sparse_destroy(mat_buffer)

  end function sparse_check_hermitian

  !============================================================================!
  ! This function calculates the largest imaginary part of a sparse matrix     !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   mat   (input)  : The sparse complex matrix to check                      !
  !----------------------------------------------------------------------------!
  ! Written by Andrea Greco, September 2016                                    !
  !============================================================================!

  real(kind=DP) function sparse_max_imag_part(mat)

    use comms, only: comms_reduce
    use utils, only: utils_assert

    implicit none

    ! Arguments
    type(SPAM3), intent(in)       :: mat           ! input complex matrix

    ! Local variables
    integer                       :: my_nze        ! number of nonzero elements
    real(kind=DP)                 :: loc_max_imag  ! local max imaginary part
    real(kind=DP)                 :: glob_max_imag ! global max imaginary part

    call utils_assert(mat%iscmplx, 'ERROR in routine sparse_max_imag_part: &
         &matrix must be complex')

    ! jme: Get number of nonzero elements in the proc
    my_nze = library(mat%lib)%my_nze

    ! get max value of imaginary part on this proc
    if (my_nze > 0) then
       loc_max_imag = maxval(abs(aimag(mat%zmtx(1:my_nze))))
    else
       loc_max_imag = 0.0_DP
    end if

    !  across processors
    glob_max_imag = loc_max_imag
    call comms_reduce('MAX',glob_max_imag)

    sparse_max_imag_part = glob_max_imag

  end function sparse_max_imag_part

  !============================================================================!
  ! This subroutine performs a safe conversion of a complex matrix to a real   !
  ! one. The imaginary part of the complex source matrix is checked against    !
  ! a specified threshold.                                                     !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   dest  (output) : The resulting sparse real matrix                        !
  !   src   (input)  : The original sparse complex matrix                      !
  !   threshold      : max value allowed for the imaginary part of src         !
  !----------------------------------------------------------------------------!
  ! Written by Andrea Greco, September 2016                                    !
  !============================================================================!

  subroutine sparse_take_real_part(dest,src,threshold)

    use utils, only: utils_assert, utils_abort

    implicit none

    ! Arguments
    type(SPAM3), intent(in)       :: src           ! input complex matrix
    type(SPAM3), intent(inout)    :: dest          ! output real matrix
    real(kind=DP), intent(in)     :: threshold     ! threshold to accept input matrix

    ! Local variables
    real(kind=DP)                 :: max_imag      ! global max imaginary part
    character(len=100)            :: message       ! error message

    call utils_assert(src%iscmplx, 'ERROR in routine sparse_take_real_part: &
         &source matrix must be complex')
    call utils_assert(.not.(dest%iscmplx), 'ERROR in routine &
         &sparse_take_real_part: destination matrix must be real')

    ! calculates largest imaginary part of source matrix
    max_imag = sparse_max_imag_part(src)

    ! agrecocmplx: abort if imaginary part of complex source matrix
    ! is too large
    if (max_imag>threshold) then
       write(message,'(a,f15.8)') 'Error in routine sparse_take_real_part: &
          &largest imaginary part of source is ', max_imag
       call utils_abort(message)
    ! agrecocmplx: call sparse_copy with complx_to_real=.true. since
    ! it is safe to throw away the imaginary part
    else
       call sparse_copy(dest,src,cmplx_to_real=.true.)
    end if

  end subroutine sparse_take_real_part

  !============================================================================!
  ! This subroutine performs a safe axpy operation involving a complex matrix  !
  ! ymat and a real one xmat. The imaginary part of xmat is checked against a  !
  ! specified threshold. The alpha parameter is real.                          !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   ymat  (output) : The resulting sparse real matrix                        !
  !   xmat   (input)  : The original sparse complex matrix                     !
  !   alpha  (input)  : y = y + alpha * x                                      !
  !   threshold      : max value allowed for the imaginary part of src         !
  !----------------------------------------------------------------------------!
  ! Written by Andrea Greco, September 2016                                    !
  !============================================================================!

  subroutine sparse_safe_axpy_real(ymat,xmat,alpha,threshold)

    use utils, only: utils_assert, utils_abort

    implicit none

    ! Arguments
    type(SPAM3), intent(in)       :: xmat          ! input complex matrix
    type(SPAM3), intent(inout)    :: ymat          ! output real matrix
    real(kind=DP), intent(in)     :: alpha
    real(kind=DP), intent(in)     :: threshold     ! threshold to accept input matrix

    ! Local variables
    real(kind=DP)                 :: max_imag      ! global max imaginary part
    character(len=100)            :: message       ! error message

    call utils_assert(xmat%iscmplx, 'ERROR in routine sparse_safe_axpy: &
         &x matrix must be complex')
    call utils_assert(.not.(ymat%iscmplx), 'ERROR in routine &
         &sparse_safe_axpy: y matrix must be real')

    ! calculates largest imaginary part of source matrix
    max_imag = sparse_max_imag_part(xmat)

    ! agrecocmplx: abort if imaginary part of complex source matrix
    ! is too large
    if (max_imag>threshold) then
       write(message,'(a,f15.8)') 'Error in routine sparse_safe_axpy: &
          &largest imaginary part of x is ', max_imag
       call utils_abort(message)
    ! agrecocmplx: call sparse_axpy with complx_to_real=.true. since
    ! it is safe to throw away the imaginary part
    else
       call sparse_axpy_real(ymat,xmat,alpha,cmplx_to_real=.true.)
    end if

  end subroutine sparse_safe_axpy_real

  !============================================================================!
  ! This subroutine performs a safe axpy operation involving a complex matrix  !
  ! ymat and a real one xmat. The imaginary part of xmat is checked against a  !
  ! specified threshold. The alpha parameter is complex.                       !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   ymat  (output) : The resulting sparse real matrix                        !
  !   xmat   (input)  : The original sparse complex matrix                     !
  !   alpha  (input)  : y = y + alpha * x                                      !
  !   threshold      : max value allowed for the imaginary part of src         !
  !----------------------------------------------------------------------------!
  ! Written by Andrea Greco, September 2016                                    !
  !============================================================================!

  subroutine sparse_safe_axpy_complex(ymat,xmat,alpha,threshold)

    use utils, only: utils_assert, utils_abort

    implicit none

    ! Arguments
    type(SPAM3), intent(in)       :: xmat          ! input complex matrix
    type(SPAM3), intent(inout)    :: ymat          ! output real matrix
    complex(kind=DP), intent(in)  :: alpha
    real(kind=DP), intent(in)     :: threshold     ! threshold to accept input matrix

    ! Local variables
    real(kind=DP)                 :: max_imag      ! global max imaginary part
    character(len=100)            :: message       ! error message

    call utils_assert(xmat%iscmplx, 'ERROR in routine sparse_safe_axpy: &
         &x matrix must be complex')
    call utils_assert(.not.(ymat%iscmplx), 'ERROR in routine &
         &sparse_safe_axpy: y matrix must be real')

    ! calculates largest imaginary part of source matrix
    ! agrecocmplx: should we check the imaginary part of only xmat,
    ! or alpha*xmat? since alpha is complex, is it possible that
    ! xmat has non-zero imag part but alpha*xmat does?
    ! currently in sparse_axpy, alpha is multiplied by the real part
    ! of xmat, instead of taking the real part of alpha*xmat
    max_imag = sparse_max_imag_part(xmat)

    ! agrecocmplx: abort if imaginary part of complex source matrix
    ! is too large
    if (max_imag>threshold) then
       write(message,'(a,f15.8)') 'Error in routine sparse_safe_axpy: &
          &largest imaginary part of x is ', max_imag
       call utils_abort(message)
    ! agrecocmplx: call sparse_axpy with complx_to_real=.true. since
    ! it is safe to throw away the imaginary part
    else
       call sparse_axpy_complex(ymat,xmat,alpha,cmplx_to_real=.true.)
    end if

  end subroutine sparse_safe_axpy_complex


  !============================================================================!
  ! This subroutine performs a matmul operation on the matrices:               !
  !   C := A . B                                                               !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   cmat  (inout) : The sparse matrix C                                      !
  !   amat  (input) : The sparse matrix A                                      !
  !   bmat  (input) : The sparse matrix B                                      !
  !----------------------------------------------------------------------------!
  ! SPAM3 version written by Nicholas Hine, June 2009.                         !
  ! Allows mixed sparse-dense matrices divided segment-by-segment.             !
  ! Comms system revised by Nicholas Hine in February 2012.                    !
  ! OpenMP parallelised by Nicholas Hine in May 2013.                          !
  !----------------------------------------------------------------------------!
  ! SPAM2 version written by Peter Haynes, March 2004.                         !
  ! Revised for distributed data, June 2006.                                   !
  ! Modified by Nicholas Hine, November 2007 to loop-unroll block matrix       !
  ! multiplication in the inner loop.                                          !
  ! Modified by Nicholas Hine, January 2008 to remove blocking comms.          !
  ! Re-written by Nicholas Hine, November 2008 to improve comms efficiency on  !
  ! very large numbers of procs.                                               !
  !============================================================================!

  subroutine sparse_product(cmat,amat,bmat,allow_mix_types)

    use comms, only: comms_send, comms_barrier, &
         pub_my_proc_id, pub_total_num_procs, &
         pub_comms_group_size, pub_num_comms_groups, pub_first_proc_in_group
    use linalg, only: linalg_dgemm_serial, linalg_zgemm_serial
    use rundat, only: pub_timings_level
!$  use rundat, only: pub_threads_max
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_assert, &
         utils_abort
#ifdef ITC_TRACE
    use vt
#endif

    implicit none

    ! Arguments
    type(SPAM3), intent(inout) :: cmat   ! The sparse matrix C
    type(SPAM3), intent(in) :: amat      ! The sparse matrix A
    type(SPAM3), intent(in) :: bmat      ! The sparse matrix B
    logical, intent(in), optional :: allow_mix_types

    ! Local variables
    integer :: alib,blib,clib   ! Pointers to library entries
    integer :: arow_blks        ! Identifier for A row, C row blocking scheme
    integer :: bcol_blks        ! Identifier for B col, C col blocking scheme
    integer :: acol_blks        ! Identifier for A col, B row blocking scheme
    integer :: iblk,jblk,kblk   ! Atom block loop counters
    integer :: recvproc         ! Proc index of recv origin
    integer :: loc_iblk         ! Atom block loop counter local to this proc
    integer :: loc_jblk         ! Atom block loop counter local to this proc
    integer :: iidx,jidx        ! Indices
    integer :: step             ! Proc step counter
    integer :: nrows            ! Number of rows of A and C in block/segment
    integer :: ncols            ! Number of cols of B and C in block/segment
    integer :: nsum             ! Number of cols of A and rows of B
    integer :: aptr,bptr,cptr   ! Pointer to blocks/segs in matrices A, B and C
    integer :: ielem,jelem      ! Element row/col counters
    integer :: ielem0,jelem0    ! Element row/col block start positions
    integer :: ielem1,jelem1    ! Element row/col block start positions
    integer :: ac_row_seg       ! Segment index of row segment in A and C
    integer :: bc_col_seg       ! Segment index of col segment in B and C
    integer :: bc_col_seg_proc  ! Proc corresponding to bc_col_seg
    integer :: segprod          ! Loop counter for ac and bc segs
    integer :: aseg_start       ! Start position in A index of current segment
    integer :: bseg_start       ! Start position in B index of current segment
    integer :: cseg_start       ! Start position in C index of current segment
    integer :: lda              ! Leading dimension of array storing block of A
    integer :: ldb              ! Leading dimension of array storing block of B
    integer :: ldc              ! Leading dimension of array storing block of C
    integer :: ifound,nfound    ! Index and number of found blocks of C
    logical :: a_dense,a_sparse ! Shorthand flags for type of segment of A
    logical :: b_dense,b_sparse ! Shorthand flags for type of segment of B
    logical :: c_dense,c_sparse ! Shorthand flags for type of segment of C
    integer :: reqproc          ! Proc to request from
    integer :: nextproc         ! Proc being received by master thread
    integer :: next_reqproc     ! Next proc to be requested from
    integer :: reqstep          ! Loop used to find next step to request data
    integer :: ierr             ! Error flag

    ! Internal arrays for commonly-used 'reverse row/column look-up'
    integer, allocatable :: tfound_idx(:)
    integer, allocatable :: tfound_ptr(:)
    logical, allocatable :: tfound_blk(:)

    type(COM3) :: acom
    type(SHARE3) :: bshare
    type(SHARE3) :: cshare

    ! rc2013: local pointers to row and column parallel strategies
    type(PARAL_INFO), pointer :: arow_par, acol_par
    type(PARAL_INFO), pointer :: brow_par, bcol_par
    type(PARAL_INFO), pointer :: crow_par, ccol_par

    character(len=*), parameter :: myself = 'sparse_product'
    logical :: loc_allow_mix_types

!$  integer, external :: omp_get_thread_num
    ! BLAS subroutines
    external :: sparse_zgemm_149, sparse_dgemm_149

    ! Extract library entries
    alib = amat%lib
    blib = bmat%lib
    clib = cmat%lib

    ! jd: Sanity check against garbage lib indices. This helps detects errors
    !     (see e.g. #1895) that only manifest in release mode and so the
    !     RTL's bound-check doesn't help.
    call utils_assert(alib >= lbound(library,1) .and. &
         alib <= ubound(library,1), myself//': alib out of range', alib)
    call utils_assert(blib >= lbound(library,1) .and. &
         blib <= ubound(library,1), myself//': blib out of range', blib)
    call utils_assert(clib >= lbound(library,1) .and. &
         clib <= ubound(library,1), myself//': clib out of range', clib)

    arow_blks = library(alib)%row_blks
    bcol_blks = library(blib)%col_blks
    acol_blks = library(alib)%col_blks

    ! rc2013: assign the pointers for par based on library structure
    ! In fact, need crow=arow, ccol=bcol, acol=brow (sparse_par_check)
    ! Assign them all anyway for clarity later
    arow_par => library(alib)%row_par
    brow_par => library(blib)%row_par
    crow_par => library(clib)%row_par
    acol_par => library(alib)%col_par
    bcol_par => library(blib)%col_par
    ccol_par => library(clib)%col_par

    ! rc2013: check that the parallel strategies match
    call sparse_par_check(clib,alib,blib,'sparse_product')

    ! Special timer level for measuring individual sparse_product calls
    if (iand(pub_timings_level,8) /= 0) then
       !call timer_clock('sparse_product_' &
       !     //trim(library(alib)%structure)//'_' &
       !     //trim(library(blib)%structure)//'_' &
       !     //trim(library(clib)%structure),1)
    end if
    call timer_clock('sparse_product',1)
#ifdef ITC_TRACE
    call VTBEGIN(vt_sparse_product,vt_err)
#endif

    ! Check consistency of matrix types and blocking schemes
    call utils_assert( (acol_blks == library(blib)%row_blks).and.&
         (arow_blks == library(clib)%row_blks).and.&
         (bcol_blks == library(clib)%col_blks), &
         'Error in sparse_product: incompatible blocking schemes.')

    ! jme: check types
    loc_allow_mix_types = .false.
    if (present(allow_mix_types)) loc_allow_mix_types = allow_mix_types

    if ( (.not. cmat%iscmplx) .and. &
         (amat%iscmplx .or. bmat%iscmplx) ) then
       ! real result but one of the inputs is complex
       call utils_abort('Error in '//myself//': incompatible argument types: &
            &iscmplx components of cmat (result), amat, bmat are: ', &
            opt_logical_to_print1=cmat%iscmplx, &
            opt_logical_to_print2=amat%iscmplx, &
            opt_logical_to_print3=bmat%iscmplx)
    else if ( ( (cmat%iscmplx .neqv. bmat%iscmplx) .or. &
                (cmat%iscmplx .neqv. amat%iscmplx) ) .and. &
              (.not. loc_allow_mix_types) ) then
       ! types are mixed [C is complex; A or B (or both) are real]
       ! but this has not been enabled
       call utils_abort('Error in '//myself//': incompatible argument types: &
            &mixing of types not enabled, but iscmplx components of &
            &cmat (result), amat, bmat are: ', &
            opt_logical_to_print1=cmat%iscmplx, &
            opt_logical_to_print2=amat%iscmplx, &
            opt_logical_to_print3=bmat%iscmplx)
    end if

    ! Allocate arrays and initialise comms
    call sparse_com_allocate(acom,amat,3, &
         alloc_mtx=.true.,cropped=.true.,seg=.false.,groupseg=.false., &
         cmplx_duplicate=(cmat%iscmplx.and.(.not.amat%iscmplx)))

    ! Set up which other procs will be requesting data
    do recvproc=0,pub_total_num_procs-1
       if (library(blib)%idx_seg_lens(recvproc)/=0) then
          acom%index_reqs(recvproc) = index_not_sent
          acom%data_reqs(recvproc) = data_not_sent
       end if
    end do

    ! Initialise B, C data sharing
    if (iand(pub_timings_level,8) /= 0) call timer_clock('sparse_product_init_share',1)
    call sparse_init_share(bshare,bmat,.true.,.false., &
         cmplx_duplicate=(cmat%iscmplx.and.(.not.bmat%iscmplx)))
    call sparse_init_share(cshare,cmat,.false.,.false.)
    if (iand(pub_timings_level,8) /= 0) call timer_clock('sparse_product_init_share',2)

    ! Initialise amat comms
    call sparse_init_comms(acom,grouped=.true.)

    recvproc = -1; reqstep = -1

    ! The product to be formed is:
    !  C   = A   B
    !   ki    kj  ji

!$OMP PARALLEL NUM_THREADS(pub_threads_max) DEFAULT(none) &
!$OMP PRIVATE(bc_col_seg,bc_col_seg_proc,b_sparse,b_dense,nextproc,next_reqproc, &
!$OMP      ac_row_seg,a_sparse,a_dense,c_sparse,c_dense,nrows,ncols,nsum, &
!$OMP      aptr,bptr,cptr,cseg_start,loc_iblk,iblk,ielem0,ielem1,iidx, &
!$OMP      jblk,jidx,jelem0,jelem1,ielem,jelem,aseg_start,bseg_start,reqstep, &
!$OMP      lda,ldb,ldc,nfound,loc_jblk,kblk,ifound,reqproc,segprod,step, &
!$OMP      tfound_blk,tfound_idx,tfound_ptr,ierr) &
!$OMP SHARED(pub_timings_level,pub_first_proc_in_group,recvproc, &
!$OMP      pub_total_num_procs,pub_comms_group_size,amat,acol_blks,acom, &
!$OMP      bcol_blks,bshare,cshare,pub_num_comms_groups,pub_my_proc_id, &
!$OMP      arow_par,brow_par,crow_par,acol_par,bcol_par,ccol_par, &
!$OMP      arow_blks,library,clib,blib, &
!$OMP      pub_threads_max)

    allocate(tfound_idx(library(clib)%nblk),stat=ierr)
    call utils_alloc_check('sparse_product','tfound_idx',ierr)
    allocate(tfound_ptr(library(clib)%nblk),stat=ierr)
    call utils_alloc_check('sparse_product','tfound_ptr',ierr)
    allocate(tfound_blk(library(clib)%nblk),stat=ierr)
    call utils_alloc_check('sparse_product','tfound_blk',ierr)
    tfound_blk = .false.

    reqproc = -1; nextproc = -1; next_reqproc = -1

    ! Loop over processors
    do step=0,pub_num_comms_groups-1

!$     if (omp_get_thread_num()==0) then
       ! Receive index, pointers and data for this step if required
       nextproc = modulo(pub_my_proc_id-step*pub_comms_group_size &
            +pub_total_num_procs,pub_total_num_procs)
       if (any(bshare%seginfobuf(s_type,nextproc,0:pub_comms_group_size-1) &
            /=SEG_BLANK)) then
          if (.not.acom%got_data) then
             call sparse_get_crop_step_data(acom,amat,nextproc,.false., &
                  maxval(library(blib)%idx_lens), &
                  maxval(library(clib)%idx_lens), &
                  bshare%idxbuf,cshare%idxbuf, &
                  bshare%seginfobuf,cshare%seginfobuf, &
                  tfound_blk,tfound_idx,bcol_par)
           end if
       end if

       ! Request index and data for next required step
       do reqstep=step+1,pub_num_comms_groups-1

          next_reqproc = modulo(pub_my_proc_id-(reqstep)*pub_comms_group_size &
               +pub_total_num_procs,pub_total_num_procs)

          if (any(bshare%seginfobuf(s_type,next_reqproc,0:pub_comms_group_size-1) &
               /=SEG_BLANK).and.(next_reqproc/=pub_my_proc_id)) then

             if (next_reqproc==reqproc) exit ! already requested

             reqproc = next_reqproc

             ! Send request for index and pointers
             call comms_send(reqproc,index_needed,1,tag=IDX_REQ_TAG)

             ! Start asynchronous receives for index and pointers
             call sparse_recv_index(acom,reqproc,1,async=.true.)
             call sparse_recv_pointers(acom,reqproc,1,async=.true.)
             acom%got_index = .false.
             acom%got_reqs = .false.
             acom%got_data = .false.

             exit
          end if

       end do

!$     end if

!$     if (iand(pub_timings_level,8) /= 0) call timer_clock('sparse_product_omp_barrier',1)
!$OMP BARRIER
!$     if (iand(pub_timings_level,8) /= 0) call timer_clock('sparse_product_omp_barrier',2)

!$     if (omp_get_thread_num()==0) then
       recvproc = nextproc
       ! If there has been data received from next proc, swap buffer marked
       ! as "used" buffer, so next proc's data does not overwrite it.
       if (any(bshare%seginfobuf(s_type,nextproc,0:pub_comms_group_size-1) &
            /=SEG_BLANK)) then
          acom%ubuf = acom%cbuf
          acom%cbuf = acom%cbuf + 1
          if (acom%cbuf>acom%nbuf) acom%cbuf = 2

          if (cshare%iscmplx.and.(.not.amat%iscmplx)) then
             ! jme: copy real recvbuf into complex recvbuf
             acom%zmtxrecvbuf(:,acom%ubuf) = &
                  cmplx(acom%dmtxrecvbuf(:,acom%ubuf), 0.0_DP, kind=DP)
          end if
       end if

!$     end if

!$     if (iand(pub_timings_level,8) /= 0) call timer_clock('sparse_product_omp_barrier',1)
!$OMP FLUSH(acom,recvproc)
!$OMP BARRIER
!$     if (iand(pub_timings_level,8) /= 0) call timer_clock('sparse_product_omp_barrier',2)

       ! Loop over col segments of B and C in this group of procs
!$OMP DO SCHEDULE(DYNAMIC)
       do segprod=0,pub_comms_group_size*pub_total_num_procs-1

          bc_col_seg = segprod / pub_total_num_procs
          ac_row_seg = modulo(segprod,pub_total_num_procs)

          ! Find the proc index in the full list for this bc_col_seg
          bc_col_seg_proc = bc_col_seg + pub_first_proc_in_group

          ! No work for this recvproc if corresponding segment of B is blank
          if (bshare%seginfobuf(s_type,recvproc,bc_col_seg)==SEG_BLANK) cycle

          ! Find out if corresponding section of B is sparse or dense
          b_sparse = (bshare%seginfobuf(s_type,recvproc,bc_col_seg)==SEG_SPARSE)
          b_dense = (bshare%seginfobuf(s_type,recvproc,bc_col_seg)==SEG_DENSE)

          ! Skip this segment if A or C are blank
          if ((acom%seginfobuf(s_type,ac_row_seg,acom%ubuf)==SEG_BLANK).or.&
               (cshare%seginfobuf(s_type,ac_row_seg,bc_col_seg)==SEG_BLANK)) cycle

          ! Find out information about these segments of A and C
          a_sparse = (acom%seginfobuf(s_type,ac_row_seg,acom%ubuf)==SEG_SPARSE)
          a_dense = (acom%seginfobuf(s_type,ac_row_seg,acom%ubuf)==SEG_DENSE)
          c_sparse = (cshare%seginfobuf(s_type,ac_row_seg,bc_col_seg)==SEG_SPARSE)
          c_dense = (cshare%seginfobuf(s_type,ac_row_seg,bc_col_seg)==SEG_DENSE)
          nrows = arow_par%num_elems_on_proc(ac_row_seg,arow_blks) ! rows of A and C
          ncols = bcol_par%num_elems_on_proc(bc_col_seg_proc,bcol_blks) ! cols of C
          nsum = acol_par%num_elems_on_proc(recvproc,acol_blks) ! cols of A
          aptr = acom%seginfobuf(s_ptr,ac_row_seg,acom%ubuf)
          bptr = bshare%seginfobuf(s_ptr,recvproc,bc_col_seg)
          cptr = cshare%seginfobuf(s_ptr,ac_row_seg,bc_col_seg)

          if (a_dense.and.b_dense.and.c_dense) then
             if (iand(pub_timings_level,8) /= 0) call timer_clock('sparse_product_dense_dgemm',1)
#ifdef ITC_TRACE
             call VTBEGIN(vt_sparse_product_dense,vt_err)
#endif
             ! Multiply segments of A and B with Lapack call
             if (cshare%iscmplx) then
                call linalg_zgemm_serial(cshare%zmtxbuf(cptr,bc_col_seg), &
                     acom%zmtxrecvbuf(aptr,acom%ubuf), &
                     bshare%zmtxbuf(bptr,bc_col_seg), &
                     nrows,ncols,nsum,nrows,nsum,nrows)
             else
                call linalg_dgemm_serial(cshare%dmtxbuf(cptr,bc_col_seg), &
                     acom%dmtxrecvbuf(aptr,acom%ubuf), &
                     bshare%dmtxbuf(bptr,bc_col_seg), &
                     nrows,ncols,nsum,nrows,nsum,nrows)
             end if
             if (iand(pub_timings_level,8) /= 0) call timer_clock('sparse_product_dense_dgemm',2)
#ifdef ITC_TRACE
             call VTEND(vt_sparse_product_dense,vt_err)
#endif
          else if (a_dense.and.b_dense.and.c_sparse) then
             if (iand(pub_timings_level,8) /= 0) call timer_clock('sparse_product_picks',1)
#ifdef ITC_TRACE
             call VTBEGIN(vt_sparse_product_picks,vt_err)
#endif
             ! Loop over elements of C only, calculating product explicitly
             cseg_start = cshare%seginfobuf(s_idx,ac_row_seg,bc_col_seg)

             ! Loop over col-blocks of C
             loc_iblk = 0
             do iblk=ccol_par%first_atom_on_proc(bc_col_seg_proc), &
                  ccol_par%first_atom_on_proc(bc_col_seg_proc+1)-1
                loc_iblk = loc_iblk + 1

                ! Find range of col offsets in this block of C
                ielem0 = bcol_par%first_elem_on_atom(iblk,bcol_blks) - &
                     bcol_par%first_elem_on_proc(bc_col_seg_proc,bcol_blks)
                ielem1 = ielem0 + bcol_par%num_elems_on_atom(iblk,bcol_blks) - 1

                ! Loop over block-rows jblk in block-column iblk of C
                do iidx=cshare%idxbuf(cseg_start+loc_iblk-1,bc_col_seg), &
                     cshare%idxbuf(cseg_start+loc_iblk,bc_col_seg)-1
                   jblk = cshare%idxbuf(iidx,bc_col_seg)
                   cptr = cshare%ptrbuf(iidx,bc_col_seg)

                   ! Find range of row offsets in this block of C
                   jelem0 = arow_par%first_elem_on_atom(jblk,arow_blks) - &
                        arow_par%first_elem_on_proc(ac_row_seg,arow_blks)
                   jelem1 = jelem0 + arow_par%num_elems_on_atom(jblk,arow_blks) - 1

                   ! Loop over elements of this block of C
                   if (cshare%iscmplx) then
                      do ielem=ielem0,ielem1
                         do jelem=jelem0,jelem1
                            cshare%zmtxbuf(cptr,bc_col_seg) = &
                                 cshare%zmtxbuf(cptr,bc_col_seg) + sum( &
                                 bshare%zmtxbuf(bptr+ielem*nsum:bptr+(ielem+1)*nsum-1:1,bc_col_seg) * &
                                 acom%zmtxrecvbuf(aptr+jelem:aptr+jelem+nrows*(nsum-1):nrows,acom%ubuf))
                            cptr = cptr + 1
                         end do
                      end do
                   else
                     do ielem=ielem0,ielem1
                         do jelem=jelem0,jelem1

                            cshare%dmtxbuf(cptr,bc_col_seg) = cshare%dmtxbuf(cptr,bc_col_seg) + &
                                 sum(bshare%dmtxbuf(bptr+ielem*nsum:bptr+(ielem+1)*nsum-1:1,bc_col_seg) * &
                                 acom%dmtxrecvbuf(aptr+jelem:aptr+jelem+nrows*(nsum-1):nrows,acom%ubuf))
                            cptr = cptr + 1
                         end do
                      end do
                   end if

                end do  ! iidx
             end do  ! iblk
             if (iand(pub_timings_level,8) /= 0) call timer_clock('sparse_product_picks',2)
#ifdef ITC_TRACE
             call VTEND(vt_sparse_product_picks,vt_err)
#endif
          else ! any other combination of sparse and dense segments
             if (iand(pub_timings_level,8) /= 0) call timer_clock('sparse_product_sparse',1)
#ifdef ITC_TRACE
             call VTBEGIN(vt_sparse_product_sparse,vt_err)
#endif
             ! Set up shorthand variables for these segments
             aseg_start = acom%seginfobuf(s_idx,ac_row_seg,acom%ubuf)
             bseg_start = bshare%seginfobuf(s_idx,recvproc,bc_col_seg)
             cseg_start = cshare%seginfobuf(s_idx,ac_row_seg,bc_col_seg)
             if (a_dense) lda = arow_par%num_elems_on_proc(ac_row_seg,arow_blks)
             if (b_dense) ldb = acol_par%num_elems_on_proc(recvproc,acol_blks)
             if (c_dense) ldc = arow_par%num_elems_on_proc(ac_row_seg,arow_blks)

             ! Loop over block-columns iblk of B and C on this proc
             ! rc2013: use col_par
             loc_iblk = 0
             do iblk=ccol_par%first_atom_on_proc(bc_col_seg_proc), &
                  ccol_par%first_atom_on_proc(bc_col_seg_proc+1)-1
                loc_iblk = loc_iblk + 1

                ! Set number of columns in this block-column of C
                ncols = bcol_par%num_elems_on_atom(iblk,bcol_blks)

                ! Only proceed if there are nonzero block-rows in this column
                ! in this segment of B
                if (bshare%idxbuf(bseg_start+loc_iblk-1,bc_col_seg)== &
                     bshare%idxbuf(bseg_start+loc_iblk,bc_col_seg)) cycle

                ! Set number of blocks found to zero
                nfound = 0

                ! Loop over block-rows kblk of C in column iblk
                do iidx=cshare%idxbuf(cseg_start+loc_iblk-1,bc_col_seg), &
                     cshare%idxbuf(cseg_start+loc_iblk,bc_col_seg)-1
                   kblk = cshare%idxbuf(iidx,bc_col_seg)

                   ! Mark this block-row as found and note row and pointer
                   tfound_blk(kblk) = .true.
                   nfound = nfound + 1
                   tfound_idx(nfound) = kblk
                   tfound_ptr(kblk) = cshare%ptrbuf(iidx,bc_col_seg)

                end do  ! Loop over block-rows kblk of C in column iblk

                ! Don't check B index if C has no nze's in this col of
                ! this A-row segment
                if (nfound==0) cycle

                ! Loop over block-rows jblk of B in column iblk
                do iidx=bshare%idxbuf(bseg_start+loc_iblk-1,bc_col_seg), &
                     bshare%idxbuf(bseg_start+loc_iblk,bc_col_seg)-1
                   jblk = bshare%idxbuf(iidx,bc_col_seg)

                   ! Find number of this atom on recvproc
                   loc_jblk = jblk - brow_par%first_atom_on_proc(recvproc) + 1
                   !loc_jblk = jblk - par%first_atom_on_proc(recvproc) + 1

                   ! Set number of rows in this block-row of B
                   nsum = acol_par%num_elems_on_atom(jblk,acol_blks)
                   if (b_sparse) ldb = nsum

                   ! Copy pointer to this block in B
                   bptr = bshare%ptrbuf(iidx,bc_col_seg)

                   ! Loop over block-rows kblk of A in column jblk
                   do jidx=acom%idxbuf(aseg_start+loc_jblk-1,acom%ubuf), &
                        acom%idxbuf(aseg_start+loc_jblk,acom%ubuf)-1
                      kblk = acom%idxbuf(jidx,acom%ubuf)
                      ! If this row occurs in C do block multiply
                      if (tfound_blk(kblk)) then
                         nrows = arow_par%num_elems_on_atom(kblk,arow_blks)
                         if (a_sparse) lda = nrows
                         if (c_sparse) ldc = nrows
                         cptr = tfound_ptr(kblk)
                         aptr = acom%ptrbuf(jidx,acom%ubuf)

                         if (cshare%iscmplx) then
                            if ((nrows==1.or.nrows==4.or.nrows==9).and. &
                                (ncols==1.or.ncols==4.or.ncols==9).and. &
                                (nsum ==1.or.nsum ==4.or.nsum ==9)) then
                               call sparse_zgemm_149(nrows,ncols,nsum, &
                                    acom%zmtxrecvbuf(aptr,acom%ubuf),lda, &
                                    bshare%zmtxbuf(bptr,bc_col_seg),ldb, &
                                    cshare%zmtxbuf(cptr,bc_col_seg),ldc)
                            else
                               call linalg_zgemm_serial(&
                                    cshare%zmtxbuf(cptr,bc_col_seg), &
                                    acom%zmtxrecvbuf(aptr,acom%ubuf), &
                                    bshare%zmtxbuf(bptr,bc_col_seg), &
                                    nrows,ncols,nsum,lda,ldb,ldc)
                            end if
                         else
                            if ((nrows==1.or.nrows==4.or.nrows==9).and. &
                                 (ncols==1.or.ncols==4.or.ncols==9).and. &
                                 (nsum ==1.or.nsum ==4.or.nsum ==9)) then
                               call sparse_dgemm_149(nrows,ncols,nsum, &
                                    acom%dmtxrecvbuf(aptr,acom%ubuf),lda, &
                                    bshare%dmtxbuf(bptr,bc_col_seg),ldb, &
                                    cshare%dmtxbuf(cptr,bc_col_seg),ldc)
                            else
                               call linalg_dgemm_serial(&
                                    cshare%dmtxbuf(cptr,bc_col_seg), &
                                    acom%dmtxrecvbuf(aptr,acom%ubuf), &
                                    bshare%dmtxbuf(bptr,bc_col_seg), &
                                    nrows,ncols,nsum,lda,ldb,ldc)
                            end if
                         end if

                      end if

                   end do  ! Column blocks kblk of A in row jblk

                end do  ! Block-rows jblk of B in column iblk

                ! Reset found block-rows flags
                do ifound=1,nfound
                   iidx = tfound_idx(ifound)
                   tfound_blk(iidx) = .false.
                end do

             end do  ! Loop over block-columns iblk of B and C
             if (iand(pub_timings_level,8) /= 0) call timer_clock('sparse_product_sparse',2)
#ifdef ITC_TRACE
             call VTEND(vt_sparse_product_sparse,vt_err)
#endif
          end if

          ! Check incoming requests once per segment (Master thread only)
!$        if (omp_get_thread_num()==0) then
          call sparse_check_send_requests(acom,amat)
          ! Check ongoing receive of next step's data
          if (reqproc>=0) then
             if (any(bshare%seginfobuf(s_type,reqproc,0:pub_comms_group_size-1) &
                  /=SEG_BLANK).and.(reqproc/=pub_my_proc_id)) then
                if (.not.acom%got_index) then
                   call sparse_check_next_recv(acom,reqproc, &
                        maxval(library(blib)%idx_lens), &
                        maxval(library(clib)%idx_lens), &
                        bshare%idxbuf,cshare%idxbuf,bshare%seginfobuf, &
                        cshare%seginfobuf,tfound_blk,tfound_idx,bcol_par,'L')
                end if
                if (acom%got_index.and.(.not.acom%got_data)) then
                   call sparse_get_crop_step_data(acom,amat,reqproc,.true., &
                        maxval(library(blib)%idx_lens), &
                        maxval(library(clib)%idx_lens), &
                        bshare%idxbuf,cshare%idxbuf, &
                        bshare%seginfobuf,cshare%seginfobuf, &
                        tfound_blk,tfound_idx,bcol_par)
                end if

             end if
          end if
!$        end if

       end do  ! segprod
!$OMP END DO NOWAIT

    end do  ! step

    deallocate(tfound_blk,stat=ierr)
    call utils_dealloc_check('sparse_exit','tfound_blk',ierr)
    deallocate(tfound_ptr,stat=ierr)
    call utils_dealloc_check('sparse_exit','tfound_ptr',ierr)
    deallocate(tfound_idx,stat=ierr)
    call utils_dealloc_check('sparse_exit','tfound_idx',ierr)

!$OMP END PARALLEL

    ! Wait for final sends and deallocate communication buffers
    if (iand(pub_timings_level,8) /= 0) call timer_clock('sparse_product_exit_comms',1)
    call sparse_exit_comms(acom,amat)
    call sparse_com_deallocate(acom,dealloc_mtx=.true.,cropped=.true.)
    if (iand(pub_timings_level,8) /= 0) call timer_clock('sparse_product_exit_comms',2)

    if (iand(pub_timings_level,8) /= 0) call timer_clock('sparse_product_exit_share',1)
    call sparse_exit_share(bshare,.false.,.false.)
    call sparse_exit_share(cshare,.true.,.false.,cmat)
    if (iand(pub_timings_level,8) /= 0) call timer_clock('sparse_product_exit_share',2)

    call comms_barrier

    if (any(library(clib)%seg_info(s_type,:)==SEG_DENSE)) then
       call sparse_enforce_sparsity(cmat)
    end if

    ! ndmh: special timer level for measuring individual sparse_product calls
    call timer_clock('sparse_product',2)
    if (iand(pub_timings_level,8) /= 0) then
       !call timer_clock('sparse_product_' &
       !     //trim(library(alib)%structure)//'_' &
       !     //trim(library(blib)%structure)//'_' &
       !     //trim(library(clib)%structure),2)
    end if
#ifdef ITC_TRACE
    call VTEND(vt_sparse_product,vt_err)
#endif

  end subroutine sparse_product

  subroutine sparse_print_index(lib)

   !========================================================================!
   ! This subroutine prints the cropped indices of A, B and C for debugging !
   !------------------------------------------------------------------------!
   ! Arguments:                                                             !
   !   None                                                                 !
   !------------------------------------------------------------------------!
   ! Written by Nicholas Hine, Nov 2008.                                    !
   !========================================================================!

   use comms, only: pub_my_proc_id

   ! Arguments
   integer, intent(in) :: lib

   ! Local Variables
   integer :: nprint
   integer, parameter :: pp=10

   !write(pub_my_proc_id+60,'(a,i3,a,i3)')'PROC ',pub_my_proc_id, &
   !     ' from ',recvproc
   !write(pub_my_proc_id+60,'(a,a)') 'amat%structure=', amat%structure
   !nprint = library(alib)%idx_lens(recvproc)-1
   !call internal_fmt_int_array('acom%idxrecvbuf(:,1)',nprint,pp, &
   !     acom%idxbuf(1:nprint,1))
   !nprint = library(alib)%idx_lens(recvproc)-1
   !call internal_fmt_int_array('acom%idxrecvbuf(:,2)',nprint,pp, &
   !     acom%idxbuf(1:nprint,2))
   !nprint = library(alib)%idx_lens(recvproc)-1
   !call internal_fmt_int_array('acom%ptrrecvbuf(:,1)',nprint,pp, &
   !     acom%ptrbuf(1:nprint,1))
   !nprint = library(alib)%idx_lens(recvproc)-1
   !call internal_fmt_int_array('acom%ptrrecvbuf(:,2)',nprint,pp, &
   !     acom%ptrbuf(1:nprint,2))

   write(pub_my_proc_id+60,'(a,a)') 'lib%structure=', library(lib)%structure
   nprint = library(lib)%idx_lens(pub_my_proc_id)-1
   call internal_fmt_int_array('library(lib)%blk_idx',nprint,pp, &
        library(lib)%blk_idx(1:nprint))
   call internal_fmt_int_array('library(lib)%blk_idx',nprint,pp, &
        library(lib)%blk_ptr(1:nprint+1))

contains

    subroutine internal_fmt_int_array(label,nprint,pp,array)

      !========================================================================!
      ! This subroutine prints an array of integers in a nicely-spaced format  !
      !------------------------------------------------------------------------!
      ! Arguments:                                                             !
      !   None                                                                 !
      !------------------------------------------------------------------------!
      ! Written by Nicholas Hine, Nov 2008.                                    !
      !========================================================================!

       ! Arguments
       character(*),intent(in) :: label
       integer,intent(in) :: nprint
       integer,intent(in) :: pp
       integer,intent(in) :: array(nprint)

       ! Locals
       integer :: nlines,nrem
       character(len=20) :: fmt, tmp1, tmp2, tmp3

       write(pub_my_proc_id+60,'(a)') label
       nlines = nprint / pp
       nrem = modulo(nprint,pp)
       write(tmp1,'(i10)') nlines
       write(tmp2,'(i10)') pp
       write(tmp3,'(i10)') nrem
       if (nrem>0) then
          write(fmt,'(7a)')'(',trim(adjustl(tmp1)),'(',trim(adjustl(tmp2)), &
               'i8/)',trim(adjustl(tmp3)),'i8)'
       else
          write(fmt,'(5a)')'(',trim(adjustl(tmp1)),'(',trim(adjustl(tmp2)), &
               'i8/))'
       end if
       write(pub_my_proc_id+60,fmt) array(1:nprint)

    end subroutine internal_fmt_int_array

  end subroutine sparse_print_index

  !============================================================================!
  ! This subroutine crops the received index of A so that only the nonzero     !
  ! blocks which contribute to C on this proc are present, and forms           !
  ! the pointer request buffer for the required sections of the data of A      !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   com  (input) : The sparse matrix communicator of the matrix (COM3 type)  !
  !   recvproc (input) : The remote proc whose index is being processed.       !
  !   bbuflen,cbuflen,bidxbuf,cidxbuf,bseginfobuf,cseginfobuf (input) : Matrix !
  !    index sizes and indices for the B and C matrices in this sparse product !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine in July 2009                                      !
  ! Modified in July 2013 to send full segments if cropping is not worthwhile  !
  ! Modified by Robert Charlton to use correct parallel strategies, 21/04/18.  !
  !============================================================================!

  subroutine sparse_crop_index(com,recvproc, &
       bbuflen,cbuflen,bidxbuf,cidxbuf,bseginfobuf,cseginfobuf, &
       tfound_blk,tfound_idx,bcol_par)

    use comms, only: pub_comms_group_size, pub_first_proc_in_group, &
         pub_total_num_procs
#ifdef ITC_TRACE
    use vt
#endif

    implicit none

    ! Arguments
    type(COM3), intent(inout) :: com    ! rc2013: COM3 for A
    integer, intent(in) :: recvproc
    integer, intent(in) :: bbuflen, cbuflen
    integer, intent(in) :: bidxbuf(1:bbuflen,0:pub_comms_group_size-1)
    integer, intent(in) :: cidxbuf(1:cbuflen,0:pub_comms_group_size-1)
    integer, intent(in) :: bseginfobuf(1:3,0:pub_total_num_procs, &
           0:pub_comms_group_size-1)
    integer, intent(in) :: cseginfobuf(1:3,0:pub_total_num_procs, &
           0:pub_comms_group_size-1)
    integer, intent(inout) :: tfound_idx(:)
    logical, intent(inout) :: tfound_blk(:)
    type(PARAL_INFO), pointer, intent(in) :: bcol_par ! Parallel strategy for B cols

    ! Local Variables
    integer :: jidx_crop
    integer :: aseg_start_crop
    integer :: nstride
    logical :: any_b_dense
    logical :: a_dense,a_sparse ! Shorthand flags for type of segment of A
    integer :: nrows, ncols, ifound, nfound
    integer :: ac_row_seg
    integer :: iidx, loc_iblk, kblk
    integer :: aseg_start, bseg_start, cseg_start
    integer :: bc_col_seg, bc_col_seg_proc
    integer :: jidx, jblk, loc_jblk
    integer :: send_all_thresh

#ifdef ITC_TRACE
    call VTBEGIN(vt_sparse_product_crop,vt_err)
#endif

    ! Reset the cropped index for A
    com%idxbuf(1:library(com%lib)%idx_lens(recvproc),com%cbuf) = 0
    com%ptrbuf(1:library(com%lib)%idx_lens(recvproc),com%cbuf) = 0
    com%seginfobuf(:,:,com%cbuf) = 0

    ! Initialisations
    com%nreqs = 0
    nstride = -1
    com%reqdatlen = 1
    send_all_thresh = 1

    ! Check if any of the row-segments corresponding to recvproc
    ! in this group of col-segments of B is dense
    any_b_dense = any(bseginfobuf(s_type,recvproc, &
         0:pub_comms_group_size-1)==SEG_DENSE)

    ! Loop over (row) segments of index of A and C
    do ac_row_seg=0,pub_total_num_procs-1

       ! Nothing to do here if A is blank for this segment
       if (com%seginfobuf(s_type,ac_row_seg,1)==SEG_BLANK) cycle
       ! Or if this row of C is blank for all the segments in this group
       if (all(cseginfobuf(s_type,ac_row_seg, &
            0:pub_comms_group_size-1)==SEG_BLANK)) cycle

       ! Find index start position for this segment of A
       aseg_start = com%seginfobuf(s_idx,ac_row_seg,1)

       ! Find out if this segment of A is sparse or dense and set nstride
       a_sparse = (com%seginfobuf(s_type,ac_row_seg,1)==SEG_SPARSE)
       a_dense = (com%seginfobuf(s_type,ac_row_seg,1)==SEG_DENSE)
       if (a_dense) nstride = com%row_par%num_elems_on_proc(ac_row_seg, &
            library(com%lib)%row_blks)

       ! Mark where this segment will start in received data buffer
       com%seginfobuf(s_ptr,ac_row_seg,com%cbuf) = com%reqdatlen

       ! If both A and any B segs are dense, request this whole segment of A
       if (a_dense.and.any_b_dense) then
          ! Find number of rows and columns of this segment
          ncols = com%col_par%num_elems_on_proc(recvproc, &
               library(com%lib)%col_blks)
          nrows = com%row_par%num_elems_on_proc(ac_row_seg, &
               library(com%lib)%row_blks)
          ! Record request pointer for A data on other proc
          com%nreqs = com%nreqs + 1
          com%ptrreqsendbuf(1,com%nreqs) = com%seginfobuf(s_ptr,ac_row_seg,1)
          com%ptrreqsendbuf(2,com%nreqs) = ncols
          com%ptrreqsendbuf(3,com%nreqs) = nrows
          com%ptrreqsendbuf(4,com%nreqs) = nstride

          ! Add this whole segment to the cropped index of A
          ! Loop over atoms on recvproc (cols of segment of A)
          do loc_iblk=1,com%col_par%num_atoms_on_proc(recvproc)

             ! Loop over block-rows kblk of A in column iblk
             do jidx=com%idxbuf(aseg_start+loc_iblk-1,1), &
                  com%idxbuf(aseg_start+loc_iblk,1)-1

                ! Add this block to the cropped index, adjusting pointer
                ! to correct offset
                com%idxbuf(jidx,com%cbuf) = com%idxbuf(jidx,1)
                com%ptrbuf(jidx,com%cbuf) = com%ptrbuf(jidx,1) &
                     - com%seginfobuf(s_ptr,ac_row_seg,1) + com%reqdatlen
             end do

          end do

          ! Advance the request pointer past this whole segment
          com%reqdatlen = com%reqdatlen + nrows*ncols

       else

       ! None of the B segments which will multiply this A segment in this
       ! group are dense, so we have to go through all the B segments and
       ! mark which blocks need to be requested from the proc that holds A.
       do bc_col_seg=0,pub_comms_group_size-1

          if (bseginfobuf(s_type,recvproc,bc_col_seg)==SEG_BLANK) cycle
          if (cseginfobuf(s_type,ac_row_seg,bc_col_seg)==SEG_BLANK) cycle

          bc_col_seg_proc = bc_col_seg + pub_first_proc_in_group

          cseg_start = cseginfobuf(s_idx,ac_row_seg,bc_col_seg)
          bseg_start = bseginfobuf(s_idx,recvproc,bc_col_seg)

          ! Loop over block columns of B and C
          do loc_iblk=1,bcol_par%num_atoms_on_proc(bc_col_seg_proc)

             ! Set number of blocks found to zero
             nfound = 0

             ! Loop over block-rows kblk of C in column iblk
             do iidx=cidxbuf(cseg_start+loc_iblk-1,bc_col_seg), &
                  cidxbuf(cseg_start+loc_iblk,bc_col_seg)-1
                kblk = cidxbuf(iidx,bc_col_seg)

                ! Mark this block-row as found and note index and pointer
                tfound_blk(kblk) = .true.
                nfound = nfound + 1
                tfound_idx(nfound) = kblk

             end do  ! Loop over block-rows kblk of C in column iblk

             ! Loop over block-rows jblk of B in column iblk
             do iidx=bidxbuf(bseg_start+loc_iblk-1,bc_col_seg), &
                  bidxbuf(bseg_start+loc_iblk,bc_col_seg)-1
                jblk = bidxbuf(iidx,bc_col_seg)

                ! Find local number of this atom
                loc_jblk = jblk - com%col_par%first_atom_on_proc(recvproc) + 1

                ! Set number of cols in this block-row
                ncols = com%col_par%num_elems_on_atom(jblk,library(com%lib)%col_blks)

                ! Loop over block-rows kblk of A in column jblk
                do jidx=com%idxbuf(aseg_start+loc_jblk-1,1), &
                     com%idxbuf(aseg_start+loc_jblk,1)-1
                   kblk = com%idxbuf(jidx,1)

                   ! If this row occurs in C, add this block to the
                   ! cropped index of A and request it from recvproc
                   if (tfound_blk(kblk)) then

                      nrows = com%row_par%num_elems_on_atom(kblk, &
                           library(com%lib)%row_blks)
                      if (a_sparse) nstride = nrows

                      ! If this entry has not yet been recorded in the
                      ! cropped index aidxrecvbuf(:,com%cbuf)
                      if (com%idxbuf(jidx,com%cbuf)/=kblk) then
                         com%nreqs = com%nreqs + 1
                         ! Record request pointer for A data on other proc
                         com%ptrreqsendbuf(1,com%nreqs) = com%ptrbuf(jidx,1)
                         com%ptrreqsendbuf(2,com%nreqs) = ncols
                         com%ptrreqsendbuf(3,com%nreqs) = nrows
                         com%ptrreqsendbuf(4,com%nreqs) = nstride
                         ! Add this block to the cropped index of A
                         com%idxbuf(jidx,com%cbuf) = kblk
                         com%ptrbuf(jidx,com%cbuf) = com%reqdatlen
                         com%reqdatlen = com%reqdatlen + nrows*ncols
                      end if

                   end if

                end do  ! Column blocks kblk of A in row jblk

             end do  ! Block-rows jblk of B in column iblk

             ! Reset found block-rows flags
             do ifound=1,nfound
                iidx = tfound_idx(ifound)
                tfound_blk(iidx) = .false.
             end do

          end do  ! Loop over block-columns iblk of B and C

       end do ! Col-Segments bc_col_seg of B,C in this group

       end if

    end do ! Row-Segments ac_row_seg of remote A cols

    ! Request pointer finishes at 1 beyond end of array, so set to end
    com%reqdatlen = com%reqdatlen - 1

    ! Check if it is worth sending everything, rather than the cropped buffer
    if (com%reqdatlen>=send_all_thresh * &
         com%seginfobuf(s_ptr,pub_total_num_procs,1)-1) then

       ! Copy full index into cbuf
       com%idxbuf(1:library(com%lib)%idx_lens(recvproc),com%cbuf) = &
            com%idxbuf(1:library(com%lib)%idx_lens(recvproc),1)
       com%ptrbuf(1:library(com%lib)%idx_lens(recvproc),com%cbuf) = &
            com%ptrbuf(1:library(com%lib)%idx_lens(recvproc),1)
       com%seginfobuf(:,:,com%cbuf) = com%seginfobuf(:,:,1)

       ! Ask for whole of data
       com%reqdatlen = com%seginfobuf(s_ptr,pub_total_num_procs,1)-1
       com%nreqs = -1

       return

    end if

    ! Shift the cropped indices of A down
    ! Find the start of the block list for first column of A
    jidx_crop = 1
    aseg_start_crop = 1

    do ac_row_seg=0,pub_total_num_procs-1

       if (com%seginfobuf(s_type,ac_row_seg,1)==SEG_DENSE) then
          a_dense = .true.
       else
          a_dense = .false.
       end if

       ! Mark segment start in cropped segment info
       com%seginfobuf(s_idx,ac_row_seg,com%cbuf) = aseg_start_crop

       ! Nothing to do here if A or C are blank for this segment
       if ((com%seginfobuf(s_type,ac_row_seg,1)==SEG_BLANK).or.&
            (all(cseginfobuf(s_type,ac_row_seg,0:pub_comms_group_size-1) &
            ==SEG_BLANK))) then
          com%seginfobuf(s_type,ac_row_seg,com%cbuf) = SEG_BLANK
          cycle
       end if

       ! Find segment start positions
       aseg_start = com%seginfobuf(s_idx,ac_row_seg,1)

       ! Leave space for column start positions
       jidx_crop = jidx_crop + com%col_par%num_atoms_on_proc(recvproc) + 1

       ! Loop over block-columns jblk of A
       do loc_jblk=1,com%col_par%num_atoms_on_proc(recvproc)

          ! Mark the column start in the cropped index of A
          com%idxbuf(aseg_start_crop+loc_jblk-1,com%cbuf) = jidx_crop

          ! Loop over block-rows kblk of A in column jblk
          do jidx=com%idxbuf(aseg_start+loc_jblk-1,1), &
               com%idxbuf(aseg_start+loc_jblk,1)-1
             ! Check if we recorded this block in first pass
             if (com%idxbuf(jidx,com%cbuf)>0) then
                ! Add it to the cropped index of A
                com%idxbuf(jidx_crop,com%cbuf) = com%idxbuf(jidx,com%cbuf)
                com%ptrbuf(jidx_crop,com%cbuf) = com%ptrbuf(jidx,com%cbuf)
                jidx_crop = jidx_crop + 1
             end if
          end do

       end do

       ! Mark the last column end in the cropped index of A
       com%idxbuf(aseg_start_crop+com%col_par%num_atoms_on_proc(recvproc),com%cbuf) &
            = jidx_crop

       ! If there are no blocks in this segment of the cropped index, mark
       ! it as blank so it is skipped
       if (aseg_start_crop == jidx_crop) then
          com%seginfobuf(s_type,ac_row_seg,com%cbuf) = SEG_BLANK
       else
          ! If these segments of A and B are both dense, mark cropped segment
          ! as dense, otherwise sparse
          if (a_dense.and.any_b_dense) then
             com%seginfobuf(s_type,ac_row_seg,com%cbuf) = SEG_DENSE
          else
             com%seginfobuf(s_type,ac_row_seg,com%cbuf) = SEG_SPARSE
          end if
          aseg_start_crop = jidx_crop
       end if

    end do

    ! Record end of last segment
    com%seginfobuf(s_idx,pub_total_num_procs,com%cbuf) = jidx_crop
    com%seginfobuf(s_ptr,pub_total_num_procs,com%cbuf) = com%reqdatlen + 1

#ifdef ITC_TRACE
    call VTEND(vt_sparse_product_crop,vt_err)
#endif

  end subroutine sparse_crop_index

  !============================================================================!
  ! This subroutine receives the index, segment info, pointers and data of A   !
  ! from recvproc, while checking for any new requests arriving                !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   com  (input) : The sparse matrix communicator of the matrix (COM3 type)  !
  !   mat  (input) : The sparse matrix itself (SPAM3 type)                     !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine in July 2009, later made standalone.              !
  !============================================================================!

  subroutine sparse_get_crop_step_data(com,mat,recvproc,oneshot, &
       bbuflen,cbuflen,bidxbuf,cidxbuf,bseginfobuf,cseginfobuf, &
       tfound_blk,tfound_idx,bcol_par)

    use comms, only: comms_send, comms_wait, comms_waitany, comms_test, &
         pub_my_proc_id, pub_null_handle, pub_total_num_procs, &
         pub_comms_group_size
    use rundat, only: pub_timings_level
    use timer, only: timer_clock

    implicit none

    ! Arguments
    type(COM3), intent(inout) :: com
    type(SPAM3), intent(in) :: mat
    integer, intent(in) :: recvproc
    logical, intent(in) :: oneshot
    integer, intent(in) :: bbuflen, cbuflen
    integer, intent(in) :: bidxbuf(1:bbuflen,0:pub_comms_group_size-1)
    integer, intent(in) :: cidxbuf(1:cbuflen,0:pub_comms_group_size-1)
    integer, intent(in) :: bseginfobuf(1:3,0:pub_total_num_procs, &
           0:pub_comms_group_size-1)
    integer, intent(in) :: cseginfobuf(1:3,0:pub_total_num_procs, &
           0:pub_comms_group_size-1)
    integer, intent(inout) :: tfound_idx(:)
    logical, intent(inout) :: tfound_blk(:)
    type(PARAL_INFO), pointer, intent(in) :: bcol_par

    ! Local Variables
    logical :: commstest

    if (iand(pub_timings_level,8) /= 0) &
         call timer_clock('sparse_product_get_crop_step_data',1)

    ! Local proc
    if (recvproc==pub_my_proc_id) then
       ! Call 'comms' routines to copy local data only
       call sparse_recv_index(com,recvproc,com%cbuf,.true.)
       call sparse_recv_pointers(com,recvproc,com%cbuf,.true.)
       call sparse_recv_data(com,mat,recvproc,com%cbuf,.true.)
    ! Remote proc
    else
       ! Loop until all required components have been received
       do

         ! Test if any data needs to be sent
         call sparse_check_send_requests(com,mat)

         ! Test if we have received the index of recvproc
         if (.not.com%got_index) then
            call sparse_check_next_recv(com,recvproc, &
                 bbuflen,cbuflen,bidxbuf,cidxbuf,bseginfobuf,cseginfobuf, &
                 tfound_blk,tfound_idx,bcol_par,'B')
            if (com%got_index) cycle
         end if
         if (com%got_index.and.com%nreqs==0) then
            com%got_data = .true.
            exit
         end if

         ! Check if the request send has been completed
         if (com%got_index .and. (.not.com%got_reqs) .and. &
              (com%handles(com%req_send_hdl)==pub_null_handle)) then
             call comms_wait(com%handles(com%req_send_hdl))
             com%got_reqs = .true.
             cycle
         end if

         ! Check if the data receive has been completed
         if (com%got_index .and. com%got_reqs .and. .not.(com%got_data)) then
            call comms_test(commstest,com%handles(com%data_hdl))
            if (commstest.or.(com%handles(com%data_hdl)==pub_null_handle)) then
               call comms_wait(com%handles(com%data_hdl))
               ! All done, so leave loop
               com%got_data = .true.
               exit
            end if
         end if

         ! Leave now if the current main loop step is not yet finished, as
         ! this thread still has work to do
         if (oneshot) exit

         ! This point is only reached if there are no outstanding requests
         ! and we are still waiting for the data before continuing to the next
         ! main loop step
         call comms_waitany(com%num_handles,com%handles)

       end do
    end if

    if (iand(pub_timings_level,8) /= 0) &
        call timer_clock('sparse_product_get_crop_step_data',2)

  end subroutine sparse_get_crop_step_data

  !============================================================================!
  ! This subroutine checks whether the index of a given communicator has been  !
  ! received, and if so, it proceeds to crop it so as to determine what data   !
  ! to request from the remote proc. Then it requests said data from the       !
  ! remote proc into the next available buffer.                                !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   com      (input) : The sparse matrix communicator of the matrix (COM3)   !
  !   recvproc (input) : The remote proc to receive from.                      !
  !   bbuflen,cbuflen,bidxbuf,cidxbuf,bseginfobuf,cseginfobuf (input) : Matrix !
  !    index sizes and indices for the B and C matrices in this sparse product !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine in June 2013, based on part of get_crop_step_data !
  !============================================================================!

  subroutine sparse_check_next_recv(com,recvproc, &
       bbuflen,cbuflen,bidxbuf,cidxbuf,bseginfobuf,cseginfobuf, &
       tfound_blk,tfound_idx,bcol_par,whichcall)

    use comms, only: comms_send, comms_irecv, comms_wait, comms_test, &
         pub_null_handle, pub_comms_group_size, pub_total_num_procs
    use rundat, only: pub_timings_level
    use timer, only: timer_clock

    implicit none

    ! Arguments
    type(COM3), intent(inout) :: com
    integer, intent(in) :: recvproc
    integer, intent(in) :: bbuflen, cbuflen
    integer, intent(in) :: bidxbuf(1:bbuflen,0:pub_comms_group_size-1)
    integer, intent(in) :: cidxbuf(1:cbuflen,0:pub_comms_group_size-1)
    integer, intent(in) :: bseginfobuf(1:3,0:pub_total_num_procs, &
           0:pub_comms_group_size-1)
    integer, intent(in) :: cseginfobuf(1:3,0:pub_total_num_procs, &
           0:pub_comms_group_size-1)
    integer, intent(inout) :: tfound_idx(:)
    logical, intent(inout) :: tfound_blk(:)
    type(PARAL_INFO), pointer, intent(in) :: bcol_par
    character, intent(in) :: whichcall

    ! Local Variables
    logical :: commstest

    ! Check if the index receive has been completed
    call comms_test(commstest,com%handles(com%idx_hdl))
    if (commstest.or.(com%handles(com%idx_hdl)==pub_null_handle)) then
       call comms_wait(com%handles(com%idx_hdl))
    else
       return ! no sign of it yet...
    end if

    ! Check if the info receive has been completed
    call comms_test(commstest,com%handles(com%info_hdl))
    if (commstest.or.(com%handles(com%info_hdl)==pub_null_handle)) then
       call comms_wait(com%handles(com%info_hdl))
    else
       return ! no sign of it yet...
    end if

    ! Check if the ptr receive has been completed
    call comms_test(commstest,com%handles(com%ptr_hdl))
    if (commstest.or.(com%handles(com%ptr_hdl)==pub_null_handle)) then
       call comms_wait(com%handles(com%ptr_hdl))
    else
       return ! no sign of it yet...
    end if

    com%got_index = .true.
    ! Crop the index of A received from recvproc
    ! and form the list of required sections of data
    if (iand(pub_timings_level,8) /= 0) call timer_clock('sparse_product_crop_index'//whichcall,1)
    call sparse_crop_index(com,recvproc, &
         bbuflen,cbuflen,bidxbuf,cidxbuf,bseginfobuf,cseginfobuf,tfound_blk, &
         tfound_idx, bcol_par)
    if (iand(pub_timings_level,8) /= 0) call timer_clock('sparse_product_crop_index'//whichcall,2)

    ! If there are no blocks required from recvproc, send 0 to
    ! recvproc and move on to next proc (all segs of A will be
    ! blank in the cropped index)
    if (com%nreqs==0) then
       call comms_send(recvproc,0,1,tag=DATA_REQ_TAG)
       return
    end if

    ! Send number of blocks required to recvproc as request
    call comms_send(recvproc,com%nreqs,1,tag=DATA_REQ_TAG)

    ! Start asynchronous receives for the matrix data
    if (com%iscmplx) then
       call comms_irecv(recvproc,com%zmtxrecvbuf(:,com%cbuf),com%reqdatlen, &
            tag=DATA_TAG, handle=com%handles(com%data_hdl))
    else
       call comms_irecv(recvproc,com%dmtxrecvbuf(:,com%cbuf),com%reqdatlen, &
            tag=DATA_TAG, handle=com%handles(com%data_hdl))
    end if

    ! Send the request list to recvproc (skip if nreqs==-1)
    if (com%nreqs>0) then
       call comms_send(recvproc,com%ptrreqsendbuf,com%nreqs*4, &
            tag=PTR_REQ_TAG,return_handle=com%handles(com%req_send_hdl), &
            add_to_stack=.false.)
    end if

  end subroutine sparse_check_next_recv

  !============================================================================!
  ! This subroutine receives the index, segment info, pointers and data of A   !
  ! from recvproc, while checking for any new requests arriving                !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   com  (input) : The sparse matrix communicator of the matrix (COM3 type)  !
  !   mat  (input) : The sparse matrix itself (SPAM3 type)                     !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine in July 2009, later made standalone.              !
  !============================================================================!

  subroutine sparse_get_step_data(com,mat,recvproc)

    use comms, only: comms_send, comms_irecv, comms_wait, comms_waitany, &
         pub_my_proc_id, pub_null_handle, pub_comms_group_size, &
         pub_my_rank_in_group

    implicit none

    ! Arguments
    type(COM3), intent(inout) :: com
    type(SPAM3), intent(in) :: mat
    integer, intent(in) :: recvproc

    ! Local Variables
    integer :: seg_start, seg_full_start, seg_ptr_start
    integer :: iblk, iidx, buflen
    integer :: shiftproc

    ! Local proc
    if (recvproc==pub_my_proc_id) then
       ! Call 'comms' routines to copy local data only
       call sparse_recv_index(com,recvproc,2,.true.)
       call sparse_recv_pointers(com,recvproc,2,.true.)
       call sparse_recv_data(com,mat,recvproc,2,.true.)

    ! Remote proc
    else
       ! Loop until all required components have been received
       do

         ! Test if any data needs to be sent
         call sparse_check_send_requests(com,mat)

         ! Check if the index receives have been completed
         if ( (com%handles(com%idx_hdl)==pub_null_handle) .and. &
              (com%handles(com%info_hdl)==pub_null_handle) .and. &
              (com%handles(com%ptr_hdl)==pub_null_handle).and. &
              (com%handles(com%data_hdl)==pub_null_handle)) then
             call comms_wait(com%handles(com%idx_hdl))
             call comms_wait(com%handles(com%info_hdl))
             call comms_wait(com%handles(com%ptr_hdl))
             call comms_wait(com%handles(com%data_hdl))

             ! Shift to other buffer
             com%seginfobuf(:,:,2) = com%seginfobuf(:,:,1)
             com%idxbuf(:,2) = com%idxbuf(:,1)
             com%ptrbuf(:,2) = com%ptrbuf(:,1)
             if (com%iscmplx) then
                com%zmtxrecvbuf(:,2) = com%zmtxrecvbuf(:,1)
             else
                com%dmtxrecvbuf(:,2) = com%dmtxrecvbuf(:,1)
             end if

            ! All done, so leave loop
            exit
         end if

         ! This point is only reached if there are no outstanding requests
         ! and we are still waiting for the data
         call comms_waitany(com%num_handles,com%handles)

       end do
    end if

    if (com%segonly) then
       ! Index counts from start position in full list on recvproc
       ! so shift it down by seg_start, so that the nzb's start at
       ! index position par%num_atoms_on_proc(recvproc) + 2
       com%seginfobuf(s_idx,0,2) = 1
       seg_full_start = com%idxbuf(1,2)
       do iblk=1,com%col_par%num_atoms_on_proc(recvproc)+1
          com%idxbuf(iblk,2) = com%idxbuf(iblk,2) - seg_full_start + &
               com%col_par%num_atoms_on_proc(recvproc) + 2
       end do

       ! Pointers count from first element of segment in full array
       ! so shift them down so that they start from 1 in the received array
       buflen = library(com%lib)%idx_seg_lens(recvproc)
       seg_ptr_start = com%seginfobuf(s_ptr,0,2)
       com%seginfobuf(s_ptr,0,2) = com%seginfobuf(s_ptr,0,2) - seg_ptr_start + 1
       do iblk=com%col_par%num_atoms_on_proc(recvproc)+2,buflen
          com%ptrbuf(iblk,2) = com%ptrbuf(iblk,2) - seg_ptr_start + 1
       end do
       if (recvproc == pub_my_proc_id) then
          do iblk=1,com%col_par%num_atoms_on_proc(recvproc)
             com%ptrbuf(iblk,2) = com%ptrbuf(iblk,2) - seg_ptr_start + 1
          end do
       end if

    else if (com%groupsegonly) then

       ! Index counts from start position in full list on recvproc
       ! so shift it down by seg_start, so that the nzb's start at
       ! index position com%col_par%num_atoms_on_proc(recvproc) + 2
       seg_full_start = com%seginfobuf(s_idx,0,2)
       seg_ptr_start = com%seginfobuf(s_ptr,0,2)

       ! Shift info buffer entries to be relative to start of group of segments
       do shiftproc=0,pub_comms_group_size-1
          com%seginfobuf(s_idx,shiftproc,2) = &
               com%seginfobuf(s_idx,shiftproc,2) - seg_full_start + 1
          com%seginfobuf(s_ptr,shiftproc,2) = &
               com%seginfobuf(s_ptr,shiftproc,2) - seg_ptr_start + 1
       end do

       ! Shift whole list of received block start positions so that they count
       ! from the start of the received buffer
       do shiftproc=0,pub_comms_group_size-1
          ! Ignore this segment if blank
          if (com%seginfobuf(s_type,shiftproc,2)==SEG_BLANK) cycle
          ! Find start position of this segment in the original full buffer
          seg_full_start = com%idxbuf(com%seginfobuf(s_idx,shiftproc,2),2)
          ! Find start position of this segment in the received buffer from
          ! com%seginfobuf, which has already been shifted
          seg_start = com%seginfobuf(s_idx,shiftproc,2)

          ! Now shift indices of start positions of each block column in
          ! this segment
          do iblk=1,com%col_par%num_atoms_on_proc(recvproc)+1
             com%idxbuf(seg_start+iblk-1,2) = &
                  com%idxbuf(seg_start+iblk-1,2) - seg_full_start + &
                  com%col_par%num_atoms_on_proc(recvproc) + seg_start + 1
          end do
       end do

       ! The received pointers count from the first element of the segment in
       ! the original full array, so we need to so shift them so that they
       ! are counted relative to the start of the received array in the
       ! communicator, for each of the received segments in the group
       do shiftproc=0,pub_comms_group_size-1
          ! Skip this segment if it is blank
          if (com%seginfobuf(s_type,shiftproc,2)==SEG_BLANK) cycle
          ! Find start position of this segment in received buffer
          seg_start = com%seginfobuf(s_idx,shiftproc,2)
          ! Loop over block-columns in received segment
          do iblk=1,com%col_par%num_atoms_on_proc(recvproc)
             ! Loop over nonzero blocks in received segment, and shift pointers
             do iidx=com%idxbuf(seg_start+iblk-1,2), &
                  com%idxbuf(seg_start+iblk,2)-1
                com%ptrbuf(iidx,2) = com%ptrbuf(iidx,2) - seg_ptr_start + 1
             end do
          end do
       end do

       ! If we have the local data, adjust pointers to diagonal blocks so that
       ! they count an offset from the start position in the buffer array in
       ! the communicator, rather than the full original data array.
       if (recvproc == pub_my_proc_id) then
          ! Skip this if the segment is blank
          if (com%seginfobuf(s_type,pub_my_rank_in_group,2)/=SEG_BLANK) then
             ! Shift pointers to diagonal blocks
             seg_start = com%seginfobuf(s_idx,pub_my_rank_in_group,2)
             do iblk=1,com%col_par%num_atoms_on_proc(recvproc)
                com%ptrbuf(seg_start+iblk-1,2) = &
                     com%ptrbuf(seg_start+iblk-1,2) - seg_ptr_start + 1
             end do
          end if
       end if

    end if

  end subroutine sparse_get_step_data


  !============================================================================!
  ! This subroutine receives the index and segment info of a matrix            !
  ! from recvproc, while checking for any new requests arriving                !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   com  (input) : The sparse matrix communicator of the matrix (COM3 type)  !
  !   mat  (input) : The sparse matrix itself (SPAM3 type)                     !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine in July 2009, later made standalone.              !
  !============================================================================!

  subroutine sparse_get_step_index(com,recvproc)

    use comms, only: comms_send, comms_irecv, comms_wait, comms_waitany, &
         pub_my_proc_id, pub_null_handle

    implicit none

    ! Arguments
    type(COM3), intent(inout) :: com
    integer, intent(in) :: recvproc

    ! Local proc
    if (recvproc==pub_my_proc_id) then
       ! Call 'comms' routines to copy local data only
       call sparse_recv_index(com,recvproc,2,.true.)

    ! Remote proc
    else
       ! Loop until all required components have been received
       do

         ! Test if any data needs to be sent
         call sparse_check_send_requests(com)

         ! Check if the index receives have been completed
         if ( (com%handles(com%idx_hdl)==pub_null_handle) .and. &
              (com%handles(com%info_hdl)==pub_null_handle)) then
             call comms_wait(com%handles(com%idx_hdl))
             call comms_wait(com%handles(com%info_hdl))

             ! Shift to other buffer
             com%seginfobuf(:,:,2) = com%seginfobuf(:,:,1)
             com%idxbuf(:,2) = com%idxbuf(:,1)

            ! All done, so leave loop
            exit
         end if

         ! This point is only reached if there are no outstanding requests
         ! and we are still waiting for the data
         call comms_waitany(com%num_handles,com%handles)

       end do
    end if

  end subroutine sparse_get_step_index

  !============================================================================!
  ! This subroutine checks the index and data requests from all procs and      !
  ! sends off indices, pointers and data as required. Does not exit until      !
  ! either there are no more outstanding requests, or there are outstanding    !
  ! data requests only but the data send buffer is still in use.               !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   com  (input) : The sparse matrix communicator of the matrix (COM3 type)  !
  !   mat  (input) : The sparse matrix itself (SPAM3 type)                     !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine in July 2009, later made standalone.              !
  !============================================================================!

  subroutine sparse_check_send_requests(com,mat)

    use comms, only: comms_probe, comms_recv, comms_send, &
         comms_test, comms_wait, pub_my_proc_id, &
         pub_total_num_procs, pub_comms_group_size
    use rundat, only: pub_timings_level
    use timer, only: timer_clock
#ifdef ITC_TRACE
    use vt
#endif

    implicit none

    ! Arguments
    type(COM3), intent(inout) :: com
    type(SPAM3), intent(in), optional :: mat

    ! Local Variables
    integer :: sendproc
    integer :: jreq,sendreqs
    integer :: sendptr,origptr
    integer :: nrows,ncols,nstride,icol
    integer :: procstride
    logical :: can_exit
    logical :: probe_test

    ! Start Timer
    if ((iand(pub_timings_level,8) /= 0).and.(com%cropped)) &
         call timer_clock('sparse_product_check_send_requests',1)

    ! Test whether previous data send has completed yet
    ! and ensure that ongoing mpi calls complete
    if (.not.com%send_buffer_free(1)) then
       call comms_test(com%send_buffer_free(1),com%handles(com%data_send_hdl))
    else
       call comms_probe(probe_test,modulo(pub_my_proc_id+1, &
            pub_total_num_procs))
    end if

    ! Loop over procs until no more messages are left
    sendproc = pub_my_proc_id
    can_exit = .true.
    procstride = 1
    if (com%cropped) procstride = pub_comms_group_size
    do

       ! Check if index request has been received
       if (com%index_reqs(sendproc)==index_needed) then

          ! Make sure request has finished being received
          call comms_wait(com%handles(sendproc+1))

          ! Different response depending on whether this routine is
          ! called from sparse_product (hence cropped=true), or
          ! sparse_trace/transpose/expand (hence seg_only=true), or
          ! sparse_create (hence send only index).
          if (com%cropped) then
             ! Send the index and pointers only
             call sparse_send_index(com,sendproc)
             call sparse_send_pointers(com,sendproc)
          else if (allocated(com%dmtxrecvbuf)) then
             ! Send segment index, pointers and data
             call sparse_send_index(com,sendproc)
             call sparse_send_pointers(com,sendproc)
             call sparse_send_data(com,mat,sendproc)
          else
             ! Send segment index only
             call sparse_send_index(com,sendproc)
          end if

          ! Mark completion of index send to sendproc
          com%index_reqs(sendproc) = index_sent

          ! Need to restart check, since a request came in
          can_exit = .false.

       end if

       ! Check if the send buffer has become free recently
       if ((com%data_reqs(sendproc)>=0).and.(.not.com%send_buffer_free(1))) then
          call comms_test(com%send_buffer_free(1),com%handles(com%data_send_hdl))
       end if

       if (com%data_reqs(sendproc)>=-1) then
#ifdef ITC_TRACE
          call VTBEGIN(vt_sparse_product_data_reply,vt_err)
#endif
          ! Make sure request has finished being received
          call comms_wait(com%handles(sendproc+pub_total_num_procs+1))

          if (com%data_reqs(sendproc)==0) then
             ! No data to send, so just mark this proc as done and move on
             com%data_reqs(sendproc) = data_sent
          else if (com%data_reqs(sendproc)==-1) then
             ! Send all data (no request buffer will have been sent)
             if (mat%iscmplx) then
                call comms_send(sendproc,mat%zmtx(:),library(mat%lib)%my_nze,tag=DATA_TAG)
             else
                call comms_send(sendproc,mat%dmtx(:),library(mat%lib)%my_nze,tag=DATA_TAG)
             end if

             ! Mark completion of data send to sendproc
             com%data_reqs(sendproc) = data_sent
          else
             ! Receive the pointer request buffer
             sendreqs = 4*com%data_reqs(sendproc)
             call comms_recv(sendproc,com%ptrreqrecvbuf,sendreqs,tag=PTR_REQ_TAG)
             sendreqs = com%data_reqs(sendproc)
             sendptr = 1

             ! Do not re-fill the send buffer unless the previous data send
             ! is done. The corresponding receive was posted before the
             ! send started so there is no chance of lockup.
             call comms_wait(com%handles(com%data_send_hdl))

             ! Fill com%(d/z)mtxsendbuf with the requested data
             do jreq=1,sendreqs
                origptr = com%ptrreqrecvbuf(1,jreq)
                ncols = com%ptrreqrecvbuf(2,jreq)
                nrows = com%ptrreqrecvbuf(3,jreq)
                nstride = com%ptrreqrecvbuf(4,jreq)
                if (mat%iscmplx) then
                   do icol=1,ncols
                      com%zmtxsendbuf(sendptr:sendptr+nrows-1,1) = &
                           mat%zmtx(origptr:origptr+nrows-1)
                      origptr = origptr + nstride
                      sendptr = sendptr + nrows
                   end do
                else
                   do icol=1,ncols
                      com%dmtxsendbuf(sendptr:sendptr+nrows-1,1) = &
                           mat%dmtx(origptr:origptr+nrows-1)
                      origptr = origptr + nstride
                      sendptr = sendptr + nrows
                   end do
                end if
             end do

             ! Set sendptr to the length of the data
             sendptr = sendptr - 1

             ! Send the data requested by sendproc
             ! Do not add these to the send stack - the handles are
             ! managed explicitly to prevent buffer-overwrite
             if (mat%iscmplx) then
                call comms_send(sendproc,com%zmtxsendbuf(:,1),sendptr,tag=DATA_TAG, &
                     return_handle=com%handles(com%data_send_hdl), &
                     add_to_stack=.false.)
             else
                call comms_send(sendproc,com%dmtxsendbuf(:,1),sendptr,tag=DATA_TAG, &
                     return_handle=com%handles(com%data_send_hdl), &
                     add_to_stack=.false.)
             end if

             ! Mark completion of data send to sendproc
             com%data_reqs(sendproc) = data_sent

             ! Prevent send buffer being re-used until this send is finished
             com%send_buffer_free(1) = .false.

          end if

          ! Need to restart check, since a new request may have come in
          can_exit = .false.
#ifdef ITC_TRACE
          call VTEND(vt_sparse_product_data_reply,vt_err)
#endif
       end if

       ! Advance to next proc there might be a request from
       ! Only exit if nothing happened this time round
       sendproc = sendproc + procstride
       if (sendproc>=pub_total_num_procs) &
            sendproc = sendproc - pub_total_num_procs
       if (sendproc==pub_my_proc_id) then
          if (can_exit) then
             exit
          else
             sendproc = pub_my_proc_id
             can_exit = .true.
          end if
       end if

    end do

    ! Stop Timer
    if ((iand(pub_timings_level,8) /= 0).and.(com%cropped)) &
         call timer_clock('sparse_product_check_send_requests',2)

  end subroutine sparse_check_send_requests


  !============================================================================!
  ! This subroutine initialises storage for sharing matrix data between the    !
  ! members of a comms group.                                                  !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   share (input) : The sparse matrix sharer for the matrix (SHARE3 type)    !
  !   mat   (input) : The sparse matrix itself (SPAM3 type)                    !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine in May 2013 based on code from July 2011.         !
  !============================================================================!

  subroutine sparse_init_share(share,mat,share_data,idx_only,cmplx_duplicate)

    use comms, only: comms_allgather, comms_free, pub_total_num_procs, &
         pub_comms_group_size, pub_my_rank_in_group, pub_first_proc_in_group, &
         pub_last_proc_in_group, pub_my_proc_id, pub_group_comm
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_assert

    implicit none

    ! Arguments
    type(SHARE3), intent(inout) :: share
    type(SPAM3), intent(in) :: mat
    logical, intent(in) :: share_data
    logical, intent(in) :: idx_only
    logical, intent(in), optional :: cmplx_duplicate

    ! Local Variables
    integer :: sendproc
    integer :: buflen, datlen
    integer :: ierr
    integer, allocatable :: displs(:), lengths(:)
    logical :: loc_cmplx_duplicate

    ! jme: will we create a complex duplicate of the local data?
    ! (if mat%iscmplx this will not be used at all)
    loc_cmplx_duplicate = .false.
    if (present(cmplx_duplicate)) loc_cmplx_duplicate = cmplx_duplicate

    share%lib = mat%lib
    share%iscmplx = mat%iscmplx

    ! jd: Sanity check against garbage lib indices. This helps detects errors
    !     (see e.g. #1895) that only manifest in release mode and so the
    !     RTL's bound-check doesn't help.
    call utils_assert(mat%lib >= lbound(library,1) .and. &
         mat%lib <= ubound(library,1), &
         'sparse_init_share: mat%lib out of range', mat%lib)

    ! rc2013: point the par strategies to the library stuctures.
    share%row_par=>library(mat%lib)%row_par
    share%col_par=>library(mat%lib)%col_par

    ! Allocate buffers for shared data
    allocate(share%seginfobuf(1:3,0:pub_total_num_procs, &
         0:pub_comms_group_size-1),stat=ierr)
    call utils_alloc_check('sparse_init_share','share%seginfobuf',ierr)
    buflen = maxval(library(share%lib)%idx_lens)
    allocate(share%idxbuf(1:buflen,0:pub_comms_group_size-1),stat=ierr)
    call utils_alloc_check('sparse_init_share','share%idxbuf',ierr)
    if (.not.idx_only) then
       allocate(share%ptrbuf(1:buflen,0:pub_comms_group_size-1),stat=ierr)
       call utils_alloc_check('sparse_init_share','share%ptrbuf',ierr)
       datlen = library(share%lib)%max_nze
       if (share%iscmplx) then
          allocate(share%zmtxbuf(datlen,0:pub_comms_group_size-1),stat=ierr)
          call utils_alloc_check('sparse_init_share','share%zmtxbuf',ierr)
          allocate(share%dmtxbuf(1,1),stat=ierr)
          call utils_alloc_check('sparse_init_share','share%dmtxbuf',ierr)
       else
          allocate(share%dmtxbuf(datlen,0:pub_comms_group_size-1),stat=ierr)
          call utils_alloc_check('sparse_init_share','share%dmtxbuf',ierr)
          if (loc_cmplx_duplicate) then
             allocate(share%zmtxbuf(datlen,0:pub_comms_group_size-1),stat=ierr)
             call utils_alloc_check('sparse_init_share','share%zmtxbuf',ierr)
          else
             allocate(share%zmtxbuf(1,1),stat=ierr)
             call utils_alloc_check('sparse_init_share','share%zmtxbuf',ierr)
          end if
       end if
    end if

    ! Length and displacement buffers for allgatherv
    allocate(displs(0:pub_comms_group_size-1),stat=ierr)
    call utils_alloc_check('sparse_init_share','displs',ierr)
    allocate(lengths(0:pub_comms_group_size-1),stat=ierr)
    call utils_alloc_check('sparse_init_share','lengths',ierr)

    ! Share segment info, index and pointers for mat over all procs in group
    call comms_allgather(share%seginfobuf(1:3,0:pub_total_num_procs,0), &
         library(share%lib)%seg_info(1:3,0:pub_total_num_procs), &
         comm=pub_group_comm)
    lengths(0:pub_comms_group_size-1) = library(share%lib)%idx_lens( &
         pub_first_proc_in_group:pub_last_proc_in_group)
    do sendproc=0,pub_comms_group_size-1
       displs(sendproc) = sendproc*size(share%idxbuf,1)
    end do
    buflen = max(library(share%lib)%idx_lens(pub_my_proc_id),1)
    call comms_allgather(share%idxbuf(:,0),library(share%lib)%blk_idx(1:buflen), &
         length_src=lengths(pub_my_rank_in_group), &
         lengths_dest=lengths,displs_dest=displs,comm=pub_group_comm)
    if (.not.idx_only) then
       call comms_allgather(share%ptrbuf(:,0),library(share%lib)%blk_ptr(1:buflen), &
            length_src=lengths(pub_my_rank_in_group), &
            lengths_dest=lengths,displs_dest=displs,comm=pub_group_comm)
    end if

    if ((.not.share_data).or.(idx_only)) then

       if (.not.idx_only) then
          ! Zero result matrix C
          if (share%iscmplx) then
             share%zmtxbuf = (0.0_DP,0.0_DP)
          else
             share%dmtxbuf = 0.0_DP
          end if
       end if

    else

       lengths(0:pub_comms_group_size-1) = &
            share%seginfobuf(s_ptr,pub_total_num_procs,0:pub_comms_group_size-1) - 1
       do sendproc=0,pub_comms_group_size-1
          if (share%iscmplx) then
             displs(sendproc) = sendproc*size(share%zmtxbuf,1)
          else
             displs(sendproc) = sendproc*size(share%dmtxbuf,1)
          end if
       end do

       datlen = max(library(share%lib)%my_nze,1)
       if (share%iscmplx) then
          call comms_allgather(share%zmtxbuf(:,0),mat%zmtx(1:datlen), &
               length_src=lengths(pub_my_rank_in_group), &
               lengths_dest=lengths,displs_dest=displs,comm=pub_group_comm)
       else
          call comms_allgather(share%dmtxbuf(:,0),mat%dmtx(1:datlen), &
               length_src=lengths(pub_my_rank_in_group), &
               lengths_dest=lengths,displs_dest=displs,comm=pub_group_comm)
       end if

    end if

    ! Free up memory in comms_mod
    call comms_free

    ! jme: copy data (duplicate)
    if (loc_cmplx_duplicate .and. (.not. share%iscmplx)) then
       share%zmtxbuf(:,:) = cmplx(share%dmtxbuf, 0.0_DP, kind=DP)
    end if

    ! Deallocate length and displacement buffers for allgatherv
    deallocate(displs,stat=ierr)
    call utils_dealloc_check('sparse_init_share','displs',ierr)
    deallocate(lengths,stat=ierr)
    call utils_dealloc_check('sparse_init_share','lengths',ierr)

  end subroutine sparse_init_share


  !============================================================================!
  ! This subroutine deallocates storage for sharing matrix data between the    !
  ! members of a comms group, and performs a sum over the data if required.    !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   share (input) : The sparse matrix sharer for the matrix (SHARE3 type)    !
  !   mat   (input) : The sparse matrix itself (SPAM3 type)                    !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine in May 2013 based on code from July 2011.         !
  !============================================================================!

  subroutine sparse_exit_share(share,sum_data,idx_only,mat)

    use comms, only: comms_reduce, comms_free,pub_comms_group_size, &
         pub_total_num_procs,pub_group_comm
    use utils, only: utils_dealloc_check, utils_assert

    implicit none

    ! Arguments
    type(SHARE3), intent(inout) :: share
    logical, intent(in) :: sum_data
    logical, intent(in) :: idx_only
    type(SPAM3), intent(inout), optional :: mat

    ! Local Variables
    integer :: proc
    integer :: datlen
    integer :: ierr

    ! Sum resulting data over all procs in group
    if (sum_data) then
       call utils_assert(present(mat), 'Error in sparse_exit_share: &
            &matrix not provided, but sum_data=.true.')
       do proc=0,pub_comms_group_size-1
          datlen = share%seginfobuf(s_ptr,pub_total_num_procs,proc) - 1
          if (mat%iscmplx) then
             call comms_reduce('SUM',mat%zmtx(:), &
                  length=datlen,comm=pub_group_comm,root=proc, &
                  z_array_src=share%zmtxbuf(1:datlen,proc))
          else
             call comms_reduce('SUM',mat%dmtx(:), &
                  length=datlen,comm=pub_group_comm,root=proc, &
                  d_array_src=share%dmtxbuf(1:datlen,proc))
          end if
       end do
       call comms_free
    end if

    ! Deallocate communication buffers
    if (.not.idx_only) then
       deallocate(share%dmtxbuf,stat=ierr)
       call utils_dealloc_check('sparse_exit_share','share%dmtxbuf',ierr)
       deallocate(share%zmtxbuf,stat=ierr)
       call utils_dealloc_check('sparse_exit_share','share%zmtxbuf',ierr)
       deallocate(share%ptrbuf,stat=ierr)
       call utils_dealloc_check('sparse_exit_share','share%ptrbuf',ierr)
    end if
    deallocate(share%idxbuf,stat=ierr)
    call utils_dealloc_check('sparse_exit_share','share%idxbuf',ierr)
    deallocate(share%seginfobuf,stat=ierr)
    call utils_dealloc_check('sparse_exit_share','share%seginfobuf',ierr)

  end subroutine sparse_exit_share


  !============================================================================!
  ! This subroutine initialises the requests and handles required for the      !
  ! communications operations of sparse_product                                !
  ! The communicator must already have had its %index_reqs and %data_reqs      !
  ! set up prior to entry to this routine.                                     !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   com  (input) : The sparse matrix communicator of the matrix (COM3 type)  !
  !   mat  (input) : The sparse matrix itself (SPAM3 type)                     !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine in July 2009, later made standalone.              !
  !============================================================================!

  subroutine sparse_init_comms(com,grouped)

    use comms, only: comms_irecv, pub_my_proc_id, pub_total_num_procs, &
         pub_my_rank_in_group, pub_comms_group_size

    implicit none

    ! Arguments
    type(COM3), intent(inout) :: com
    logical, intent(in) :: grouped

    ! Local Variables
    integer :: reqproc, recvproc
    integer :: rank_in_group,group_size
    integer :: first_in_req_group,last_in_req_group

    ! Set up handle indices for various messages
    if (com%cropped) then
       com%num_handles = 2*pub_total_num_procs + 6
       com%data_send_hdl = com%num_handles - 5   ! send operation for data
       com%info_hdl = com%num_handles - 4 ! recv operation for com%seginfobuf
       com%idx_hdl = com%num_handles - 3  ! recv operation for index
       com%ptr_hdl = com%num_handles - 2  ! recv operation for pointers
       com%data_hdl = com%num_handles - 1 ! recv operation for data
       com%req_send_hdl = com%num_handles ! send operation for requests
    else
       com%num_handles = 2*pub_total_num_procs + 4
       com%info_hdl = com%num_handles - 3 ! recv operation for com%seginfobuf
       com%idx_hdl = com%num_handles - 2  ! recv operation for index
       com%ptr_hdl = com%num_handles - 1  ! recv operation for pointers
       com%data_hdl = com%num_handles     ! recv operation for data
    end if

    ! Initialise arrays and variables
    com%send_buffer_free(:) = .true.
    com%index_reqs(pub_my_proc_id) = index_sent
    com%data_reqs(pub_my_proc_id) = data_sent
    com%cbuf = 2
    com%got_index = .true.
    com%got_reqs = .true.
    com%got_data = .false. ! so as to copy in first step data

    if (grouped) then
       rank_in_group = pub_my_rank_in_group
       group_size = pub_comms_group_size
    else
       rank_in_group = 0
       group_size = 1
    end if

    ! Start asynchronous request receive for all those procs which
    ! will be requesting anything from this proc, by dint of having the
    ! same rank within their group as this proc
    do reqproc=rank_in_group,pub_total_num_procs-1,group_size
       ! Find range of procs belonging to the same group as reqprocs
       first_in_req_group = reqproc - rank_in_group
       last_in_req_group = first_in_req_group + group_size - 1
       ! Proc never needs to send to itself, so override for pub_my_proc_id
       if (reqproc==pub_my_proc_id) then
          com%index_reqs(first_in_req_group:last_in_req_group) = index_sent
          if (grouped) com%data_reqs(first_in_req_group:last_in_req_group) &
               = data_sent
       end if
       ! Check if any of the group of procs reqproc is in have nonzero
       ! index_reqs
       if (any(com%index_reqs(first_in_req_group:last_in_req_group)/=index_sent)) then
          ! Record the fact that we are waiting for an index and data request
          ! from this proc, and clear the index_reqs and data_reqs for the
          ! other procs in the same group as reqproc, since they will not
          ! request anything from this proc
          com%index_reqs(reqproc) = index_not_sent
          if (com%cropped) com%data_reqs(reqproc) = data_not_sent
          do recvproc=first_in_req_group,last_in_req_group
             if (recvproc/=reqproc) then
                com%index_reqs(recvproc) = index_sent
                if (com%cropped) com%data_reqs(recvproc) = data_sent
             end if
          end do
          ! Start receive for index request
          call comms_irecv(reqproc,com%index_reqs(reqproc),1,tag=IDX_REQ_TAG, &
               handle=com%handles(reqproc+1))
          ! Start receive for data request flag
          if (com%cropped) call comms_irecv(reqproc,com%data_reqs(reqproc),1, &
               tag=DATA_REQ_TAG,handle= &
               com%handles(reqproc+pub_total_num_procs+1))
       end if
    end do

  end subroutine sparse_init_comms

  !============================================================================!
  ! This subroutine ensures all the requests and ongoing send operations       !
  ! have completed then deallocates temporary arrays, and in the case of       !
  ! comms groups being used, sums the result matrix over the comms group       !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   com  (input) : The sparse matrix communicator of the matrix (COM3 type)  !
  !   mat  (input) : The sparse matrix itself (SPAM3 type)                     !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine in July 2009, later made standalone.              !
  !============================================================================!

  subroutine sparse_exit_comms(com,mat)

    use comms, only: comms_free, comms_probe, comms_wait, comms_waitany, &
         pub_null_handle, pub_total_num_procs

    implicit none

    ! Arguments
    type(COM3), intent(inout) :: com
    type(SPAM3), intent(in), optional :: mat

    ! Ensure all requests are dealt with before deallocating
    do
       ! If there is something left to send, check incoming requests
       call sparse_check_send_requests(com,mat)

       ! Check all requests have been dealt with and that all request handles
       ! are null before exiting
       if ((all(com%index_reqs(:) == index_sent)) .and. &
            (all(com%data_reqs(:) == data_sent)) .and. &
            (all(com%handles(1:2*pub_total_num_procs+1)==pub_null_handle))) exit

       ! Wait until something new happens: check only the requests and the
       ! data send handle (as this may be blocking the final outgoing
       ! sends from starting), as the index and pointer and data receive
       ! handles will already be null.
       call comms_waitany(com%num_handles-5,com%handles)
    end do

    ! Complete all outgoing sends before deallocating
    if (com%cropped) call comms_wait(com%handles(com%data_send_hdl))
    call comms_free

  end subroutine sparse_exit_comms


  !============================================================================!
  ! This function performs a trace on a sparse matrix or the product of two    !
  ! sparse matrices (the whole product is not calculated).                     !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   amat  (input) : The sparse matrix A                                      !
  !   bmat  (input) : The (optional) sparse matrix B                           !
  !   @docme: opA
  !----------------------------------------------------------------------------!
  ! SPAM3 version written by Nicholas Hine, June 2009.                         !
  ! Made to use asynchronous comms by Nicholas Hine, Feb 2012                  !
  ! Revised to use bmat sharing and group comms by Nicholas Hine, June 2013    !
  !----------------------------------------------------------------------------!
  ! Original SPAM2 version written by Peter Haynes, March 2004.                !
  ! Revised for distributed data, June 2006.                                   !
  ! Dense matrix version by Nicholas Hine, Dec 2007.                           !
  ! Modified by Andrea Greco to allow operation on amat, May 2016.             !
  ! Currently only conjugacy of amat is supported.                             !
  !============================================================================!

  real(kind=DP) function sparse_trace(amat,bmat,opA)

    use comms, only: comms_reduce, comms_send, pub_my_proc_id, &
         pub_total_num_procs, pub_comms_group_size, pub_num_comms_groups, &
         pub_first_proc_in_group
!$  use rundat, only: pub_threads_max
    use timer, only: timer_clock
    use utils, only: utils_assert
#ifdef ITC_TRACE
    use vt
#endif

    implicit none

    ! Arguments
    type(SPAM3), intent(in) :: amat                ! The sparse matrix A
    type(SPAM3), optional, intent(in) :: bmat      ! The sparse matrix B
    character(len=1), optional, intent(in) :: opA  ! operation to be applied to A

    ! Local variables
    integer :: alib,blib       ! Pointers to library entries
    integer :: iblk,jblk,kblk  ! Atom block loop counters
    integer :: recvproc        ! Proc index of recv origin
    integer :: reqproc         ! Proc index of next recv origin
    integer :: groupfirstproc  ! First proc in recvproc's group
    integer :: grouplastproc   ! Last proc in recvproc's group
    integer :: step            ! Proc loop counter
    integer :: loc_iblk        ! Atom block loop counter local to this proc
    integer :: loc_jblk        ! Atom block loop counter local to this proc
    integer :: iidx,jidx       ! Indices
    integer :: ielems,jelems   ! Numbers of element rows/cols in this block
    integer :: ielem           ! Loop counter over elems on atom
    integer :: jelem           ! Loop counter over elems on atom
    integer :: aptr,bptr       ! Pointers to blocks in matrices A and B
    integer :: bseg_start      ! Start position of this segment in index of B
    integer :: seg_type        ! Segment type for this segment
    integer :: arow_seg        ! Loop counter for received segments of A
    integer :: bcol_seg        ! Loop counter for shared segments of B
    integer :: bcol_seg_proc   ! Proc in full list corresponding to bcol_seg
    integer :: aseg_start      ! Start position of this segment in received index
    integer :: arow_blks       ! Identifier for A row, B col row blocking scheme
    integer :: brow_blks       ! Identifier for B row, A col row blocking scheme
    real(kind=DP) :: local_tr  ! Local trace on this proc
    type(COM3) :: acom
    type(SHARE3) :: bshare(1)
    real(kind=DP), external :: ddot
    ! agrecocmplx
    character(len=1) :: loc_opA
    ! rc2013: parallel strategies associated with the matrices
    type(PARAL_INFO), pointer :: arow_par
    type(PARAL_INFO), pointer :: acol_par
    type(PARAL_INFO), pointer :: brow_par
    type(PARAL_INFO), pointer :: bcol_par

!$  integer, external :: omp_get_thread_num

#ifdef ITC_TRACE
    call VTBEGIN(vt_sparse_trace,vt_err)
#endif

    call timer_clock('sparse_trace',1)

    ! agrecocmplx: only 'C' and 'N' are currently supported for opA
    if (present(opA)) then
        loc_opA = opA
    else
        loc_opA = 'N'
    end if

    ! qoh: Initialise variables
    blib = -1
    brow_blks = -1

    ! Extract library entries
    alib = amat%lib
    if (present(bmat)) blib = bmat%lib

    ! jd: Sanity check against garbage lib indices. This helps detects errors
    !     (see e.g. #1895) that only manifest in release mode and so the
    !     RTL's bound-check doesn't help.
    call utils_assert(alib >= lbound(library,1) .and. &
         alib <= ubound(library,1), 'sparse_trace: alib out of range', alib)
    if(present(bmat)) then
       call utils_assert(blib >= lbound(library,1) .and. &
            blib <= ubound(library,1), 'sparse_trace: blib out of range', blib)
    end if

    ! Check consistency of arguments
    if (.not.present(bmat)) then
       call utils_assert(library(alib)%row_blks==library(alib)%col_blks, &
            'Error in sparse_trace: matrix is not square.')
       ! rc2013: row/col strategies of A must match for trace
       call sparse_par_check(alib, alib, routine='sparse_trace', trace=.true.)
       arow_blks = library(alib)%row_blks
    else
       call utils_assert(amat%iscmplx .eqv. bmat%iscmplx, &
            'Error in sparse_trace: incompatible argument types.')
       call utils_assert(library(alib)%row_blks==library(blib)%col_blks, &
            'Error in sparse_trace: incompatible matrix blocking schemes.')
       call utils_assert(library(alib)%col_blks==library(blib)%row_blks, &
            'Error in sparse_trace: incompatible matrix blocking schemes.')
       arow_blks = library(alib)%row_blks
       brow_blks = library(blib)%row_blks

       ! rc2013: check for consistency of parallel strategies
       call sparse_par_check(alib, blib, routine='sparse_trace', trace=.true.)
       brow_par => library(blib)%row_par
       bcol_par => library(blib)%col_par
    end if
    arow_par => library(alib)%row_par
    acol_par => library(alib)%col_par

    ! Zero trace accumulator for this proc
    local_tr = 0.0_DP
    ielems = -1

    ! If only one argument is given, calculate the trace of matrix A
    if (.not. present(bmat)) then

       seg_type = library(alib)%seg_info(s_type,pub_my_proc_id)

       if (seg_type/=SEG_BLANK) then

          ! Loop over atom-blocks iblk on this proc
          loc_iblk = 0
          aseg_start = library(alib)%seg_info(s_idx,pub_my_proc_id)
          ! rc2013: arow_par=acol_par here.
          do iblk=acol_par%my_first_blk,acol_par%my_last_blk
             loc_iblk = loc_iblk + 1

             ! Get number of elems in this block
             if (seg_type==SEG_DENSE) then
                ielems = arow_par%num_elems_on_proc(pub_my_proc_id,arow_blks)
             else if (seg_type==SEG_SPARSE) then
                ielems = arow_par%num_elems_on_atom(iblk,arow_blks)
             end if

             ! Find pointer to start of diagonal block in this column
             aptr = library(alib)%blk_ptr(aseg_start+loc_iblk-1)

             ! Accumulate trace from this diagonal block
             if (amat%iscmplx) then
                do ielem=1,arow_par%num_elems_on_atom(iblk,arow_blks)
                   ! agrecocmplx: take conjugate of amat if specified
                   if (loc_opA=='C') then
                      local_tr = local_tr + real(conjg(amat%zmtx(aptr)),kind=DP)
                   else
                      local_tr = local_tr + real(amat%zmtx(aptr),kind=DP)
                   end if
                   aptr = aptr + ielems + 1
                end do
             else
                do ielem=1,arow_par%num_elems_on_atom(iblk,arow_blks)
                   local_tr = local_tr + amat%dmtx(aptr)
                   aptr = aptr + ielems + 1
                end do
             end if

          end do  ! Loop over atom-blocks iblk on this proc

       end if

    else  ! Calculate the trace of the product of A and B

       ! Initialise bmat data sharing
       call sparse_init_share(bshare(1),bmat,.true.,.false.)

       ! Allocate communication buffers
       call sparse_com_allocate(acom,amat,2, &
            alloc_mtx=.true.,cropped=.false.,seg=.false.,groupseg=.true.)

       ! Set up which other procs will be requesting data
       do step=0,pub_num_comms_groups-1
          recvproc = modulo(pub_my_proc_id-step*pub_comms_group_size &
               +pub_total_num_procs,pub_total_num_procs)
          groupfirstproc = int(recvproc/pub_comms_group_size) &
               *pub_comms_group_size
          grouplastproc = (int(recvproc/pub_comms_group_size)+1) &
               *pub_comms_group_size - 1
          if (any(library(alib)%seg_info(s_type,groupfirstproc:grouplastproc) &
               /=SEG_BLANK)) then
             acom%index_reqs(recvproc) = index_not_sent
          end if
       end do

       ! Initialise amat comms
       call sparse_init_comms(acom,grouped=.false.)

!$OMP PARALLEL NUM_THREADS(pub_threads_max) DEFAULT(none) &
!$OMP PRIVATE(step,arow_seg,bcol_seg,bcol_seg_proc, &
!$OMP      groupfirstproc,ielems,jelems,ielem,jelem,aptr,bptr,iblk,jblk,kblk, &
!$OMP      loc_iblk,loc_jblk,iidx,jidx,bseg_start,seg_type,aseg_start) &
!$OMP SHARED (arow_blks,brow_blks,acom,pub_num_comms_groups,pub_my_proc_id, &
!$OMP      pub_total_num_procs,pub_comms_group_size,library,alib,amat, &
!$OMP      pub_first_proc_in_group,bshare, &
!$OMP      arow_par,brow_par,bcol_par,recvproc,reqproc,loc_opA) &
!$OMP REDUCTION (+:local_tr)

       ! Loop over processors
       do step=0,pub_num_comms_groups-1

!$OMP BARRIER
!$        if (omp_get_thread_num()==0) then
          ! Receive data for this step (master thread only)
          recvproc = modulo(pub_my_proc_id-step*pub_comms_group_size &
               + pub_total_num_procs,pub_total_num_procs)
          if (library(alib)%idx_groupseg_lens(recvproc)>0) &
               call sparse_get_step_data(acom,amat,recvproc)

          ! Request index and data for next step if required
          reqproc = modulo(pub_my_proc_id-(step+1)*pub_comms_group_size &
               + pub_total_num_procs,pub_total_num_procs)
          if ((library(alib)%idx_groupseg_lens(reqproc)/=0).and. &
               (reqproc/=pub_my_proc_id)) then
             ! Send request for index, pointers and data
             call comms_send(reqproc,index_needed,1,tag=IDX_REQ_TAG)
             ! Start asynchronous receives for index and pointers
             call sparse_recv_index(acom,reqproc,1,async=.true.)
             call sparse_recv_pointers(acom,reqproc,1,async=.true.)
             call sparse_recv_data(acom,amat,reqproc,1,async=.true.)
          end if
!$        end if
!$OMP FLUSH(acom,recvproc)
!$OMP BARRIER
          if (library(alib)%idx_groupseg_lens(recvproc)==0) cycle

!$OMP DO
          do arow_seg=0,pub_comms_group_size-1

          bcol_seg = arow_seg
          groupfirstproc = int(recvproc/pub_comms_group_size) &
               *pub_comms_group_size
          bcol_seg_proc = bcol_seg + pub_first_proc_in_group

          ! Find type and start position of the corresponding segment in B
          seg_type = bshare(1)%seginfobuf(s_type,recvproc,bcol_seg)
          bseg_start = bshare(1)%seginfobuf(s_idx,recvproc,bcol_seg)

          ! Override for dense-dense trace (skips indexing)
          if ((acom%seginfobuf(s_type,arow_seg,2) == SEG_DENSE) .and. &
               (seg_type == SEG_DENSE)) then

             ! Find range of rows on A and cols on B
             ielems = arow_par%num_elems_on_proc(bcol_seg_proc,arow_blks)
             ! Find range of cols on A and rows on B
             jelems = brow_par%num_elems_on_proc(recvproc,brow_blks)

             do jelem=0,jelems-1

                ! Find first elem in src needed for dest
                aptr = acom%seginfobuf(s_ptr,arow_seg,2)+jelem*ielems
                bptr = bshare(1)%seginfobuf(s_ptr,recvproc,bcol_seg) + jelem

                if (amat%iscmplx) then
                   do ielem=0,ielems-1
                      ! agrecocmplx: take conjugate of amat if specified
                      if (loc_opA=='C') then
                         local_tr = local_tr + real(conjg(acom%zmtxrecvbuf(aptr,2)) * &
                              bshare(1)%zmtxbuf(bptr,bcol_seg),kind=DP)
                      else
                         local_tr = local_tr + real(acom%zmtxrecvbuf(aptr,2) * &
                              bshare(1)%zmtxbuf(bptr,bcol_seg),kind=DP)
                      end if
                      aptr = aptr + 1
                      bptr = bptr + jelems
                   end do
                else
                   ! call dot-product BLAS routine
                   local_tr = local_tr + ddot(ielems, &
                        acom%dmtxrecvbuf(aptr,2),1, &
                        bshare(1)%dmtxbuf(bptr,bcol_seg),jelems)
                   !do ielem=0,ielems-1
                   !   local_tr = local_tr + acom%dmtxrecvbuf(aptr,2) &
                   !        * bshare(1)%dmtxbuf(bptr,bcol_seg)
                   !   aptr = aptr + 1
                   !   bptr = bptr + jelems
                   !end do
                endif

             end do

          ! one or other segment is sparse, so go through block-by-block
          else if (.not.((seg_type == SEG_BLANK).or. &
               (acom%seginfobuf(s_type,arow_seg,2) == SEG_BLANK))) then

             ! Calculate contribution to trace
             ! Loop over block-columns iblk of B on proc of col seg
             loc_iblk = 0
             ! rc2013: loop with bcol_par (same as arow_par)
             do iblk=bcol_par%first_atom_on_proc(bcol_seg_proc), &
                  bcol_par%first_atom_on_proc(bcol_seg_proc+1)-1
                loc_iblk = loc_iblk + 1

                ! Number of elems in this block
                if (acom%seginfobuf(s_type,arow_seg,2) == SEG_DENSE) then
                   ielems = arow_par%num_elems_on_proc(bcol_seg_proc,arow_blks)
                else
                   ielems = arow_par%num_elems_on_atom(iblk,arow_blks)
                end if

                ! Loop over block-rows jblk in block-column iblk of B
                do iidx=bshare(1)%idxbuf(bseg_start+loc_iblk-1,bcol_seg), &
                     bshare(1)%idxbuf(bseg_start+loc_iblk,bcol_seg)-1
                   jblk = bshare(1)%idxbuf(iidx,bcol_seg)

                   ! Find local number of this atom
                   loc_jblk = jblk - brow_par%first_atom_on_proc(recvproc) + 1

                   ! Find start position of this segment in idxbuf
                   aseg_start = acom%seginfobuf(s_idx,arow_seg,2)

                   ! Copy pointer to this block in B
                   bptr = bshare(1)%ptrbuf(iidx,bcol_seg)
                   if (seg_type == SEG_DENSE) then
                      jelems = brow_par%num_elems_on_proc(recvproc,brow_blks)
                   else
                      jelems = brow_par%num_elems_on_atom(jblk,brow_blks)
                   end if

                   ! Loop over block-rows kblk in block-column jblk of A
                   do jidx=acom%idxbuf(aseg_start+loc_jblk-1,2), &
                        acom%idxbuf(aseg_start+loc_jblk,2)-1
                      kblk = acom%idxbuf(jidx,2)

                      ! If this is a diagonal block, it will contribute
                      if (kblk == iblk) then

                         ! Copy pointer to this block in A
                         aptr = acom%ptrbuf(jidx,2)

                         ! Perform block trace
                         if (amat%iscmplx) then
                            do ielem=0,arow_par%num_elems_on_atom(iblk,arow_blks)-1
                               do jelem=0,brow_par%num_elems_on_atom(jblk,brow_blks)-1
                                  ! agrecocmplx: take the conjugate of amat if specified
                                  if (loc_opA=='C') then
                                     local_tr = local_tr + real( &
                                          conjg(acom%zmtxrecvbuf(aptr+ielem+jelem*ielems,2)) * &
                                          bshare(1)%zmtxbuf(bptr+jelem+ielem*jelems,bcol_seg),kind=DP)
                                  else
                                     local_tr = local_tr + real( &
                                          acom%zmtxrecvbuf(aptr+ielem+jelem*ielems,2) * &
                                          bshare(1)%zmtxbuf(bptr+jelem+ielem*jelems,bcol_seg),kind=DP)
                                  end if
                               end do
                            end do
                         else
                            do ielem=0,arow_par%num_elems_on_atom(iblk,arow_blks)-1
                               do jelem=0,brow_par%num_elems_on_atom(jblk,brow_blks)-1
                                  local_tr = local_tr + &
                                       acom%dmtxrecvbuf(aptr+ielem+jelem*ielems,2) * &
                                       bshare(1)%dmtxbuf(bptr+jelem+ielem*jelems,bcol_seg)
                               end do
                            end do
                         end if
                      end if  ! Diagonal block

                   end do  ! Block-rows kblk in block-column jblk of A on proc

                end do  ! Block-rows jblk in block-column iblk of B

             end do  ! Loop over atom-blocks iblk on proc of col seg

          end if

          ! Check incoming requests after each non-skipped segment
          ! (master thread only)
!$        if (omp_get_thread_num()==0) then
          call sparse_check_send_requests(acom,amat)
!$        endif

          end do
!$OMP END DO

       end do  ! Loop over processors

       ! Wait for final sends and deallocate communication buffers
!$     if (omp_get_thread_num()==0) then
       call sparse_exit_comms(acom,amat)
!$     endif

!$OMP END PARALLEL

       call sparse_com_deallocate(acom,dealloc_mtx=.true.,cropped=.false.)

       ! Deallocate data sharing structure for bmat
       call sparse_exit_share(bshare(1),.false.,.false.)

    end if

    ! Sum up contributions over all procs
    call comms_reduce('SUM',local_tr)
    sparse_trace = local_tr

    call timer_clock('sparse_trace',2)
#ifdef ITC_TRACE
    call VTEND(vt_sparse_trace,vt_err)
#endif

  end function sparse_trace


  !============================================================================!
  ! This subroutine forms the transpose of a sparse matrix.                    !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   dest  (output) : The resulting transpose                                 !
  !   src   (input)  : The matrix to be transposed                             !
  !----------------------------------------------------------------------------!
  ! SPAM3 version written by Nicholas Hine, June 2009.                         !
  !----------------------------------------------------------------------------!
  ! SPAM2 version written by Peter Haynes, March 2004.                         !
  ! Revised for distributed data, June 2006.                                   !
  ! Dense matrix version by Nicholas Hine, Dec 2007.                           !
  !============================================================================!

  subroutine sparse_transpose(dest,src)

    use comms, only: comms_barrier, comms_free, comms_send, &
         pub_my_proc_id, pub_total_num_procs, &
         pub_comms_group_size, pub_num_comms_groups, pub_first_proc_in_group
!$  use rundat, only: pub_threads_max
    use utils, only: utils_assert

    implicit none

    ! Arguments
    type(SPAM3), intent(inout) :: dest    ! The resulting transpose matrix
    type(SPAM3), intent(in)    :: src     ! The matrix to be transposed

    ! Local variables
    integer :: destlib         ! Library entry for dest
    integer :: srclib          ! Library entry for src
    integer :: drow_scol_blks  ! Identifier for D row, S col blocking scheme
    integer :: srow_dcol_blks  ! Identifier for S row, D col blocking scheme
    integer :: recvproc        ! Proc index of recv origin
    integer :: reqproc         ! Proc index of send destination
    integer :: groupfirstproc  ! First proc in recvproc's group
    integer :: grouplastproc   ! Last proc in recvproc's group
    integer :: step            ! Proc loop counter
    integer :: iblk,jblk,kblk  ! Block loop counters
    integer :: loc_iblk        ! Loop counter for blocks on proc
    integer :: loc_jblk        ! Loop counter for blocks on proc
    integer :: iidx,jidx       ! Indices
    integer :: ielems,jelems   ! Numbers of elements in blocks/segments
    integer :: ielem,jelem     ! Element counter in block/segment
    integer :: sseg_start      ! Start position of segment of src in index
    integer :: srcptr          ! Pointer to data in src
    integer :: dseg_start      ! Start position of segment of dest in index
    integer :: destptr         ! Pointer to start of data in dest
    integer :: srow_seg        ! Loop counter for received segments of src
    integer :: dcol_seg        ! Loop counter for shared segments of dest
    integer :: dcol_seg_proc   ! Proc in full list corresponding to dcol_seg
    logical :: src_dense       ! Whether this source matrix segment is dense
    logical :: src_sparse      ! Whether this source matrix segment is sparse
    logical :: dest_dense      ! Whether this dest matrix segment is dense
    logical :: dest_sparse     ! Whether this dest matrix segment is sparse
    type(COM3) :: srccom
    type(SHARE3) :: destshare
    ! rc2013: row and column parallel strategies.
    type(PARAL_INFO), pointer :: drow_par
    type(PARAL_INFO), pointer :: dcol_par
!$  integer, external :: omp_get_thread_num

    ! Get library entries for matrices
    destlib = dest%lib
    srclib = src%lib
    drow_scol_blks = library(srclib)%col_blks
    srow_dcol_blks = library(srclib)%row_blks
    ielems = 0; jelems = 0

    ! rc2013: assign parallel strategies from library structures
    ! rc2013: drow_par=scol_par and dcol_par=srow_par
    drow_par => library(srclib)%col_par
    dcol_par => library(srclib)%row_par

    ! Check arguments
    call utils_assert(drow_scol_blks == library(destlib)%row_blks, &
         'Error in sparse_transpose: source cols and dest rows have &
         &different blocking schemes.')
    call utils_assert(srow_dcol_blks == library(destlib)%col_blks, &
         'Error in sparse_transpose: source rows and dest cols have &
         &different blocking schemes.')
    call utils_assert(dest%iscmplx .eqv. src%iscmplx, &
         'Error in sparse_transpose: mixed real/complex argument types.')

    ! rc2013: check for consistency of parallel strategies
    ! Include flag for transpose (swap row<->col check)
    call sparse_par_check(destlib, srclib, &
         routine='sparse_transpose', trace=.true.)

    ! Initialise dest data sharing
    call sparse_init_share(destshare,dest,.false.,.false.)

    ! Allocate communication buffers
    call sparse_com_allocate(srccom,src,2, &
         alloc_mtx=.true.,cropped=.false.,seg=.false.,groupseg=.true.)

    ! Set up which other procs will be requesting data
    do step=0,pub_num_comms_groups-1
       recvproc = modulo(pub_my_proc_id-step*pub_comms_group_size &
            +pub_total_num_procs,pub_total_num_procs)
       groupfirstproc = int(recvproc/pub_comms_group_size) &
            *pub_comms_group_size
       grouplastproc = (int(recvproc/pub_comms_group_size)+1) &
            *pub_comms_group_size - 1
       if (any(library(srclib)%seg_info(s_type,groupfirstproc:grouplastproc) &
            /=SEG_BLANK)) then
          srccom%index_reqs(recvproc) = index_not_sent
       end if
    end do

    ! Initialise amat comms
    call sparse_init_comms(srccom,grouped=.false.)

!$OMP PARALLEL NUM_THREADS(pub_threads_max) DEFAULT(none) &
!$OMP PRIVATE(step,srow_seg,dcol_seg,dcol_seg_proc, &
!$OMP      groupfirstproc,ielems,jelems,ielem,jelem,srcptr,destptr, &
!$OMP      iblk,jblk,kblk,loc_iblk,loc_jblk,iidx,jidx,dseg_start, &
!$OMP      sseg_start,dest_dense,src_dense,dest_sparse,src_sparse) &
!$OMP SHARED (src,srccom,destshare,drow_scol_blks,srow_dcol_blks, &
!$OMP      pub_num_comms_groups,pub_my_proc_id,pub_total_num_procs, &
!$OMP      pub_comms_group_size,library,srclib,pub_threads_max, &
!$OMP      pub_first_proc_in_group, &
!$OMP      drow_par,dcol_par,recvproc,reqproc)

    ! Loop over processors
    do step=0,pub_num_comms_groups-1

!$OMP BARRIER
!$     if (omp_get_thread_num()==0) then
       ! Receive data for this step (master thread only)
       recvproc = modulo(pub_my_proc_id-step*pub_comms_group_size &
            + pub_total_num_procs,pub_total_num_procs)
       if (library(srclib)%idx_groupseg_lens(recvproc)>0) &
            call sparse_get_step_data(srccom,src,recvproc)

       ! Request index and data for next step if required
       reqproc = modulo(pub_my_proc_id-(step+1)*pub_comms_group_size &
            + pub_total_num_procs,pub_total_num_procs)
       if ((library(srclib)%idx_groupseg_lens(reqproc)/=0).and. &
            (reqproc/=pub_my_proc_id)) then
          ! Send request for index, pointers and data
          call comms_send(reqproc,index_needed,1,tag=IDX_REQ_TAG)
          ! Start asynchronous receives for index and pointers
          call sparse_recv_index(srccom,reqproc,1,async=.true.)
          call sparse_recv_pointers(srccom,reqproc,1,async=.true.)
          call sparse_recv_data(srccom,src,reqproc,1,async=.true.)
       end if
!$     end if
!$OMP FLUSH(srccom,recvproc)
!$OMP BARRIER
       if (library(srclib)%idx_groupseg_lens(recvproc)==0) cycle

!$OMP DO
       do srow_seg=0,pub_comms_group_size-1

       dcol_seg = srow_seg
       groupfirstproc = int(recvproc/pub_comms_group_size) &
            *pub_comms_group_size
       dcol_seg_proc = dcol_seg + pub_first_proc_in_group

       ! Check if there is any work to do for this segment pair
       if (destshare%seginfobuf(s_type,recvproc,dcol_seg)==SEG_BLANK) cycle
       if (srccom%seginfobuf(s_type,srow_seg,2)==SEG_BLANK) cycle

       ! Find segment types of this pair
       dest_dense = (destshare%seginfobuf(s_type,recvproc,dcol_seg)==SEG_DENSE)
       dest_sparse = (destshare%seginfobuf(s_type,recvproc,dcol_seg)==SEG_SPARSE)
       src_dense = (srccom%seginfobuf(s_type,srow_seg,2)==SEG_DENSE)
       src_sparse = (srccom%seginfobuf(s_type,srow_seg,2)==SEG_SPARSE)

       ! Override for dense-dense transpose (skips indexing)
       if (src_dense.and.dest_dense) then

          ! Find range of cols on src and rows on dest
          jelems = drow_par%num_elems_on_proc(recvproc,drow_scol_blks)
          ! Find range of rows on src and cols on dest
          ielems = dcol_par%num_elems_on_proc(dcol_seg_proc,srow_dcol_blks)

          do jelem=0,jelems-1

             ! Find first elem in src needed for dest
             srcptr = srccom%seginfobuf(s_ptr,srow_seg,2) &
                  + jelem*ielems
             destptr = destshare%seginfobuf(s_ptr,recvproc,dcol_seg) + jelem

             ! Copy columns of proc to rows of my_proc
             if (src%iscmplx) then
                do ielem=0,ielems-1
                   destshare%zmtxbuf(destptr,dcol_seg) &
                        = conjg(srccom%zmtxrecvbuf(srcptr,2))
                   srcptr = srcptr + 1
                   destptr = destptr + jelems
                end do
             else
                do ielem=0,ielems-1
                   destshare%dmtxbuf(destptr,dcol_seg) &
                        = srccom%dmtxrecvbuf(srcptr,2)
                   srcptr = srcptr + 1
                   destptr = destptr + jelems
                end do
             end if

          end do

       ! one or other segment is sparse, so go through block-by-block
       else

          if (dest_dense) jelems = &
               drow_par%num_elems_on_proc(recvproc,drow_scol_blks)
          if (src_dense) ielems = &
               dcol_par%num_elems_on_proc(dcol_seg_proc,srow_dcol_blks)

          ! Loop over block-columns iblk of dest on proc of col seg
          loc_iblk = 0
          dseg_start = destshare%seginfobuf(s_idx,recvproc,dcol_seg)

          do iblk=dcol_par%first_atom_on_proc(dcol_seg_proc), &
               dcol_par%first_atom_on_proc(dcol_seg_proc+1)-1
             loc_iblk = loc_iblk + 1

             ! Number of elems in this block
             if (src_sparse) ielems = dcol_par%num_elems_on_atom(iblk,srow_dcol_blks)

             ! Loop over block-rows jblk in block-column iblk of dest
             do iidx=destshare%idxbuf(dseg_start+loc_iblk-1,dcol_seg), &
                  destshare%idxbuf(dseg_start+loc_iblk,dcol_seg)-1
                jblk = destshare%idxbuf(iidx,dcol_seg)

                ! Find local number of this atom on recvproc
                loc_jblk = jblk - drow_par%first_atom_on_proc(recvproc) + 1

                ! Find start position of this segment in idxbuf
                sseg_start = srccom%seginfobuf(s_idx,srow_seg,2)

                ! Find pointer to this block in dest, and stride length
                destptr = destshare%ptrbuf(iidx,dcol_seg)
                if (dest_sparse) jelems &
                     = drow_par%num_elems_on_atom(jblk,drow_scol_blks)

                ! Now find this block in src
                do jidx=srccom%idxbuf(sseg_start+loc_jblk-1,2), &
                     srccom%idxbuf(sseg_start+loc_jblk,2)-1
                   kblk = srccom%idxbuf(jidx,2)
                   if (kblk == iblk) then

                      ! Find pointer to first element in block and stride length
                      srcptr = srccom%ptrbuf(jidx,2)

                      ! Copy elements while transposing block
                      if (src%iscmplx) then
                         do jelem=0,drow_par%num_elems_on_atom(jblk,drow_scol_blks)-1
                            do ielem=0,dcol_par%num_elems_on_atom(iblk,srow_dcol_blks)-1
                               destshare%zmtxbuf(destptr+ielem*jelems+jelem,dcol_seg) = &
                                    conjg(srccom%zmtxrecvbuf(srcptr+jelem*ielems+ielem,2))
                            end do
                         end do
                      else
                         do jelem=0,drow_par%num_elems_on_atom(jblk,drow_scol_blks)-1
                            do ielem=0,dcol_par%num_elems_on_atom(iblk,srow_dcol_blks)-1
                               destshare%dmtxbuf(destptr+ielem*jelems+jelem,dcol_seg) = &
                                    srccom%dmtxrecvbuf(srcptr+jelem*ielems+ielem,2)
                            end do
                         end do
                      end if

                      exit
                   end if

                end do  ! Loop finding block in src

             end do  ! Loop over block-rows in col iblk of dest

          end do  ! Loop over block-colks iblk of dest

       end if

       ! Check incoming requests after each non-skipped segment
       ! (master thread only)
!$     if (omp_get_thread_num()==0) then
       call sparse_check_send_requests(srccom,src)
!$     endif

       end do
!$OMP END DO

    end do  ! step

    ! Wait for final sends and deallocate communication buffers
!$  if (omp_get_thread_num()==0) then
    call sparse_exit_comms(srccom,src)
!$  endif

!$OMP END PARALLEL

    call sparse_com_deallocate(srccom,dealloc_mtx=.true.,cropped=.false.)

    ! Deallocate data sharing for destmat and sum over procs in group
    call sparse_exit_share(destshare,.true.,.false.,dest)

    call comms_free
    call comms_barrier

    if (any(library(destlib)%seg_info(s_type,:)==SEG_DENSE)) then
       call sparse_enforce_sparsity(dest)
    end if

  end subroutine sparse_transpose

  !============================================================================!
  ! This subroutine expands a set of elements related by symmetry to fill a    !
  ! square sparse matrix. Choice of elements depends on variable 'pattern'     !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   mat  (inout)  : The matrix to be expanded                                !
  !   pattern (in)  : PATTERN_LOWER means lower triangle elements filled in    !
  !                   from upper triangle ones                                 !
  !                   PATTERN_ALTERNATE means alternating elements along a row !
  !                   or down a column come from upper and lower triangles     !
  !   sym  (in)     : Whether the matrix is symmetric (sym=.false. means       !
  !                   signs of elements are flipped on transposition).         !
  !----------------------------------------------------------------------------!
  ! SPAM3 version written by Nicholas Hine, May 2009 based on SPAM2 version.   !
  !----------------------------------------------------------------------------!
  ! SPAM2 version written by Peter Haynes, March 2004.                         !
  ! Revised for distributed data, June 2006.                                   !
  ! Dense matrix version by Nicholas Hine, Dec 2007.                           !
  ! Facility for anti-symmetric (-Hermitian) matrices added PDH May 2008.      !
  !============================================================================!

  subroutine sparse_expand(mat,pattern,sym)

    use comms, only: comms_barrier, comms_free, comms_send, &
         pub_my_proc_id, pub_total_num_procs
    use utils, only: utils_assert

    implicit none

    ! Arguments
    type(SPAM3), intent(inout) :: mat    ! The matrix
    integer,intent(in) :: pattern
    logical,optional,intent(in) :: sym

    ! Local variables
    integer :: ilib            ! Library entry for matrix
    integer :: col_blks        ! Identifier for column blocking scheme
    integer :: row_blks        ! Identifier for row blocking scheme
    integer :: reqproc         ! Proc index of next recv origin
    integer :: recvproc        ! Proc index of recv origin
    integer :: step            ! Proc loop counter
    integer :: iblk,jblk,kblk  ! Block loop counters
    integer :: loc_iblk        ! Loop counter for blocks on proc
    integer :: loc_jblk        ! Loop counter for blocks on proc
    integer :: iidx,jidx       ! Indices
    integer :: ielems,jelems   ! Numbers of elements in blocks/segments
    integer :: ielem,jelem     ! Element counter in block/segment
    integer :: ielem0,jelem0   ! Element row/col counter to top left of block
    integer :: seg_type        ! Segment index
    integer :: seg_start       ! Start position of this segment in the index
    integer :: srcptr          ! Pointer to source data
    integer :: destptr         ! Pointer to dest data
    real(kind=DP) :: sgn       ! Sign to multiply by while transposing block
    type(COM3) :: matcom
    type(PARAL_INFO), pointer :: row_par ! Parallel strategy for rows.
    type(PARAL_INFO), pointer :: col_par ! Parallel strategy for columns.

    ! Deal with optional argument
    sgn = 1.0_DP
    if (present(sym)) then
       if (.not. sym) sgn = -1.0_DP
    end if

    ! Get library entries for matrices
    ilib = mat%lib
    row_blks = library(ilib)%row_blks
    col_blks = library(ilib)%col_blks

    ! rc2013: assign par pointers
    row_par => library(ilib)%row_par
    col_par => library(ilib)%col_par

    ! Check this matrix can be expanded
    call utils_assert(row_blks == col_blks, &
         'Error in sparse_expand: Non-square matrices can not be expanded.')

    ! rc2013: parallel strategies must match for square matrices
    !call sparse_par_check(ilib, ilib, routine='sparse_expand', trace=.true.)

    ! Allocate communication buffers
    call sparse_com_allocate(matcom,mat,2, &
         alloc_mtx=.true.,cropped=.false.,seg=.true.,groupseg=.false.)

    ! Set up which other procs will be requesting data
    do recvproc=0,pub_total_num_procs-1
       if (library(ilib)%seg_info(s_type,recvproc)/=SEG_BLANK) then
          matcom%index_reqs(recvproc) = index_not_sent
       end if
    end do

    ! Initialise amat comms
    call sparse_init_comms(matcom,grouped=.false.)

    ! Loop over procs
    do step=0,pub_total_num_procs-1

       ! Receive data for this step
       recvproc = modulo(pub_my_proc_id-step+pub_total_num_procs, &
            pub_total_num_procs)
       if (library(ilib)%idx_seg_lens(recvproc)>0) &
            call sparse_get_step_data(matcom,mat,recvproc)

       ! Request index and data for next step if required
       reqproc = modulo(pub_my_proc_id-step-1+pub_total_num_procs, &
            pub_total_num_procs)
       if ((library(ilib)%idx_seg_lens(reqproc)/=0).and. &
            (reqproc/=pub_my_proc_id)) then

          ! Send request for index, pointers and data
          call comms_send(reqproc,index_needed,1,tag=IDX_REQ_TAG)

          ! Start asynchronous receives for index and pointers
          call sparse_recv_index(matcom,reqproc,1,async=.true.)
          call sparse_recv_pointers(matcom,reqproc,1,async=.true.)
          call sparse_recv_data(matcom,mat,reqproc,1,async=.true.)
       end if
       if (library(ilib)%idx_seg_lens(recvproc)==0) cycle

       seg_type = library(ilib)%seg_info(s_type,recvproc)

       ! Override for dense-dense expand (skips indexing)
       if ((matcom%seginfobuf(s_type,0,2) == SEG_DENSE) .and. &
            (seg_type == SEG_DENSE)) then

          ! Find range of cols in this segment
          jelems = col_par%num_elems_on_proc(recvproc,col_blks)
          jelem0 = col_par%first_elem_on_proc(recvproc,col_blks)
          ! Find range of rows in this segment
          ielems = row_par%num_elems_on_proc(pub_my_proc_id,row_blks)
          ielem0 = row_par%first_elem_on_proc(pub_my_proc_id,row_blks)

          do jelem=0,jelems-1

             ! Find first elem in src needed for dest
             srcptr = jelem*ielems + 1
             destptr = library(ilib)%seg_info(s_ptr,recvproc) + jelem

             ! Copy columns of proc to rows of my_proc
             ! agrecocmplx: fixed so we take the conjugate in
             ! the complex case (matrix needs to be Hermitian)
             if (mat%iscmplx) then
                do ielem=0,ielems-1
                   if (pattern_test()) mat%zmtx(destptr) = &
                        sgn*conjg(matcom%zmtxrecvbuf(srcptr,2))
                   srcptr = srcptr + 1
                   destptr = destptr + jelems
                end do
             else
                do ielem=0,ielems-1
                   if (pattern_test()) mat%dmtx(destptr) = &
                        sgn*matcom%dmtxrecvbuf(srcptr,2)
                   srcptr = srcptr + 1
                   destptr = destptr + jelems
                end do
             end if

          end do

       ! One or other segment is sparse, so go through block-by-block
       else if (.not.((seg_type == SEG_BLANK).or. &
            (matcom%seginfobuf(s_type,0,2) == SEG_BLANK))) then

          ! Loop over block-columns iblk of dest on this proc
          loc_iblk = 0
          seg_start = library(ilib)%seg_info(s_idx,recvproc)
          do iblk=col_par%my_first_blk,col_par%my_last_blk
             loc_iblk = loc_iblk + 1

             ! Number of elems in this block/segment, and origin column
             if (matcom%seginfobuf(s_type,0,2) == SEG_DENSE) then
                ielems = row_par%num_elems_on_proc(pub_my_proc_id,row_blks)
             else
                ielems = row_par%num_elems_on_atom(iblk,row_blks)
             end if
             ielem0 = col_par%first_elem_on_atom(iblk,col_blks)

             ! Loop over block-rows jblk in block-column iblk of dest
             do iidx=library(ilib)%blk_idx(seg_start+loc_iblk-1), &
                     library(ilib)%blk_idx(seg_start+loc_iblk)-1
                jblk = library(ilib)%blk_idx(iidx)
                jelem0 = row_par%first_elem_on_atom(jblk,row_blks)

                ! Find local number of this atom on recvproc
                loc_jblk = jblk - row_par%first_atom_on_proc(recvproc) + 1

                ! Find pointer to this block in dest, and stride length
                destptr = library(ilib)%blk_ptr(iidx)
                if (seg_type == SEG_DENSE) then
                   jelems = col_par%num_elems_on_proc(recvproc,col_blks)
                else
                   jelems = col_par%num_elems_on_atom(jblk,col_blks)
                end if

                ! Now find this block in src
                do jidx=matcom%idxbuf(loc_jblk,2),matcom%idxbuf(loc_jblk+1,2)-1
                   kblk = matcom%idxbuf(jidx,2)
                   if (kblk == iblk) then

                      ! Find pointer to first element in block and stride length
                      srcptr = matcom%ptrbuf(jidx,2)

                      ! Copy elements while expanding block
                      if (mat%iscmplx) then
                         do jelem=0,col_par%num_elems_on_atom(jblk,col_blks)-1
                            do ielem=0,row_par%num_elems_on_atom(iblk,row_blks)-1
                               if (pattern_test()) mat%zmtx(destptr+ &
                                    ielem*jelems+jelem) = sgn*conjg( &
                                    matcom%zmtxrecvbuf(srcptr+jelem*ielems+ielem,2))
                            end do
                         end do
                      else
                         do jelem=0,col_par%num_elems_on_atom(jblk,col_blks)-1
                            do ielem=0,row_par%num_elems_on_atom(iblk,row_blks)-1
                               if (pattern_test()) mat%dmtx(destptr+ &
                                    ielem*jelems+jelem) = sgn*matcom%dmtxrecvbuf( &
                                    srcptr+jelem*ielems+ielem,2)
                            end do
                         end do
                      end if

                      exit
                   end if

                end do  ! Loop finding block in src

             end do  ! Loop over block-rows in col iblk of dest

          end do  ! Loop over block-colks iblk of dest

       end if

       ! Check incoming requests after each non-skipped segment
       call sparse_check_send_requests(matcom,mat)

    end do  ! step

    ! Finish comms and deallocate communication buffers
    call sparse_exit_comms(matcom,mat)
    call sparse_com_deallocate(matcom,dealloc_mtx=.true.,cropped=.false.)

    call comms_free
    call comms_barrier

contains

    logical function pattern_test()

       ! Locals
       integer :: row,col   ! row,col of src (also col,row of dest)

       pattern_test = .false.

       row = ielem0 + ielem
       col = jelem0 + jelem

       ! Never include diagonal elements
       if (row == col) return

       ! Check for i>j for lower triangle expansion
       if (pattern == PATTERN_LOWER) then
          if (row > col) pattern_test = .true.
       end if

       if (pattern == PATTERN_ALTERNATE) then
          if (((col<row).and.(mod(col+row,2)==1)).or. &
              ((col>row).and.(mod(col+row,2)==0))) pattern_test = .true.
       end if

       return

    end function pattern_test

  end subroutine sparse_expand

  !============================================================================!
  ! This subroutine obtains the maximum eigenvalue of a given sparse matrix by !
  ! an iterative conjugate gradients procedure.                                !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   mat (input)    : The matrix whose eigenvalue is desired                  !
  !   met (input)    : The metric to use                                       !
  !   eval (output)  : The eigenvalue estimate                                 !
  !   tol (input)    : Optional tolerance for eigenvalue estimate              !
  !   min_val (input): Optional flag to converge minimum eigenvalue instead    !
  !----------------------------------------------------------------------------!
  ! Written by Peter Haynes, March 2004.                                       !
  ! Revised for distributed data, June 2006.                                   !
  ! Modified by Andrea Greco to add compatibility with complex matrices,       !
  ! February 2016.                                                             !
  !============================================================================!

  subroutine sparse_extremal_eigenvalue(mat,met,eval,tol,min_val,evec, &
       allow_warnings)

    ! agrecocmplx: add dependencies to use new types
    use comms, only: comms_bcast, comms_reduce, pub_on_root, &
         pub_my_proc_id, pub_root_proc_id, pub_total_num_procs
    use constants, only: SAFE_DIV_EPS
    use datatypes, only: FUNCTIONS, COEF, data_functions_dot, &
         data_functions_scale, data_functions_axpy, &
         data_set_to_zero, data_functions_copy
!$  use rundat, only: pub_threads_max
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_assert, &
         utils_isnan

    implicit none

    ! Arguments
    type(SPAM3), intent(in) :: mat     ! The matrix whose eigenvalue is desired
    type(SPAM3), intent(in) :: met     ! The metric to be used
    real(kind=DP), intent(out) :: eval ! The eigenvalue
    real(kind=DP), optional, intent(in) :: tol  ! Optional tolerance
    logical, optional, intent(in) :: min_val
    ! agrecocmplx
    type(FUNCTIONS), optional, intent(inout) :: evec
    logical, intent(in), optional :: allow_warnings

    ! Local variables
    integer, parameter :: maxiter = 200  ! Maximum number of iterations
    integer, parameter :: cgreset = 5    ! Conjugate gradients reset frequency
    integer :: n                         ! Vector length
    integer :: iter                      ! Iteration counter
    integer :: conv                      ! Converged iteration counter
    logical :: local_min_val             ! local copy of min_val
    logical :: loc_warnings              ! local copy of allow_warnings
    real(kind=DP) :: local_tol           ! Local copy of tolerance
    real(kind=DP) :: norm                ! Normalisation
    real(kind=DP) :: normfac             ! Renormalisation factor
    real(kind=DP) :: gdotg               ! Norm of gradient
    real(kind=DP) :: gamma               ! Conjugate gradients gamma
    real(kind=DP) :: xdotd,ddotd,ddoth,ddoty    ! Dot products for line search
    real(kind=DP) :: a,b,c,disc,q,step,oldeval  ! Variables for line search
    ! agrecocmplx
    type(FUNCTIONS) :: xvec  ! Trial eigenvector
    type(FUNCTIONS) :: yvec  ! Product of matrix and search dirn
    type(FUNCTIONS) :: grad  ! Gradient
    type(FUNCTIONS) :: dirn  ! Search direction
    type(FUNCTIONS) :: matx  ! Product of matrix and trial vector
    type(FUNCTIONS) :: metx  ! Product of metric and trial vector
    type(COEF) :: temp_coef
    logical :: loc_cmplx


    ! This routine will only function for square matrices
    call utils_assert(library(mat%lib)%row_blks==library(mat%lib)%col_blks, &
         'Error in sparse_extremal_eigenvalue: square matrices only.')

    ! Either both real or both complex matrices
    call utils_assert((mat%iscmplx .eqv. met%iscmplx), &
         'Error in sparse_extremal_eigenvalue: incompatible argument types.')

    ! Deal with optional argument
    if (present(tol)) then
       local_tol = tol
    else
       local_tol = 0.01_DP
    end if
    if (present(allow_warnings)) then
       loc_warnings = allow_warnings
    else
       loc_warnings = .false.
    end if


    if (present(min_val)) then
       local_min_val = min_val
    else
       local_min_val = .false.
    endif

    ! Find vector length
    n = library(mat%lib)%nrows

    ! agrecocmplx
    loc_cmplx = mat%iscmplx

    ! Allocate workspace
    call internal_alloc_work

    ! Initialise eigenvector to random guess on one proc and broadcast
    ! ndmh_2012_11_01: workaround for BG/Q compilation issue.
    ! agrecocmplx: take care of complex case as well
    if (pub_on_root) then
        do iter=1,n
           call random_number(norm)
           if (loc_cmplx) then
              !call random_number(gamma)
              !xvec%z(iter)=cmplx(norm,gamma,kind=DP)
              xvec%z(iter)=cmplx(norm,0.0_DP,kind=DP)
           else
              xvec%d(iter)=norm
           end if
        end do
    end if

    ! agrecocmplx
    if (loc_cmplx) then
       call comms_bcast(pub_root_proc_id,xvec%z)
    else
       call comms_bcast(pub_root_proc_id,xvec%d)
    end if

    ! Normalise
    ! agrecocmplx
    call internal_matvec(metx,met,xvec)
    temp_coef = data_functions_dot(xvec,metx)
    if (temp_coef%iscmplx) then
       norm = real(temp_coef%z,kind=DP)
    else
       norm = temp_coef%d
    end if
    ! jd: Makes debugging easier
    call utils_assert(norm > SAFE_DIV_EPS, 'sparse_extremal_eigenvalue(): &
         &Normalisation factor not positive. This almost definitely indicates &
         &a failure to update the overlap and/or inv_overlap of an NGWF_REP &
         &after NGWFs changed.')
    normfac = 1.0_DP / sqrt(norm)

    ! agrecocmplx
    call data_functions_scale(xvec,normfac)
    call data_functions_scale(metx,normfac)

    ! Initialise CG
    gdotg = 1.0_DP
    gamma = 0.0_DP
    ! agrecocmplx
    call data_set_to_zero(dirn)

    ! Obtain initial estimate of eigenvalue and matx := mat . x
    call internal_matvec(matx,mat,xvec)
    temp_coef = data_functions_dot(metx,matx)
    if (temp_coef%iscmplx) then
       eval = real(temp_coef%z,kind=DP)
    else
       eval = temp_coef%d
    end if

    ! ndmh: suggested by Victor Milman to avoid possible WIN32 parallel desync
    call comms_bcast(pub_root_proc_id,eval)

    ! Start of conjugate gradients loop
    conv = 0
    do iter=1,maxiter

       ! Obtain gradient grad
       ! agrecocmplx
       call data_functions_copy(grad,matx)
       call data_functions_axpy(grad,xvec,-eval)

       ! ndmh: avoid parallel desynchronisation of grad, which
       ! ndmh: really should not happen but does on DARWIN (??)
       ! agrecocmplx
       if (loc_cmplx) then
          call comms_bcast(pub_root_proc_id,grad%z)
       else
          call comms_bcast(pub_root_proc_id,grad%d)
       end if

       ! Obtain conjugate direction dirn
       gamma = 1.0_DP / gdotg
       call internal_matvec(yvec,met,grad)
       ! agrecocmplx
       temp_coef = data_functions_dot(grad,yvec)
       if (temp_coef%iscmplx) then
          gdotg = real(temp_coef%z,kind=DP)
       else
          gdotg = temp_coef%d
       end if

       if (mod(iter,cgreset) == 0) then
          gamma = 0.0_DP
       else
          gamma = gamma * gdotg
       end if

       ! agrecocmplx
       call data_functions_scale(dirn,gamma)
       call data_functions_axpy(dirn,grad,-1.0_DP)

       ! Orthogonalise to xvec
       ! agrecocmplx
       temp_coef = data_functions_dot(metx,dirn)
       if (temp_coef%iscmplx) then
          xdotd = real(temp_coef%z,kind=DP)
       else
          xdotd = temp_coef%d
       end if

       ! agrecocmplx
       call data_functions_axpy(dirn,xvec,-xdotd)

       ! Obtain yvec := mat . dirn
       call internal_matvec(yvec,mat,dirn)

       ! Obtain metx := met .dirn
       call internal_matvec(metx,met,dirn)

       ! Obtain line search coefficients
       ! agrecocmplx
       temp_coef = data_functions_dot(metx,dirn)
       if (temp_coef%iscmplx) then
          ddotd = real(temp_coef%z,kind=DP)
       else
          ddotd = temp_coef%d
       end if
       temp_coef = data_functions_dot(metx,matx)
       if (temp_coef%iscmplx) then
          ddoth = real(temp_coef%z,kind=DP)
       else
          ddoth = temp_coef%d
       end if
       temp_coef = data_functions_dot(metx,yvec)
       if (temp_coef%iscmplx) then
          ddoty = real(temp_coef%z,kind=DP)
       else
          ddoty = temp_coef%d
       end if

       a = ddotd * ddoth
       call comms_bcast(pub_root_proc_id,a)
       if (abs(a) < epsilon(1.0_DP)) exit

       b = eval * ddotd - ddoty
       call comms_bcast(pub_root_proc_id,b)
       if (abs(b) < epsilon(1.0_DP)) exit

       c = -ddoth
       disc = b*b - 4.0_DP*a*c
       call comms_bcast(pub_root_proc_id,c)
       call comms_bcast(pub_root_proc_id,disc)
       if (disc < 0.0_DP) exit

       !tjz21: Check if min or max eval is requested
       if(local_min_val) then
          q = -0.5_DP * (b - sign(sqrt(disc),b))
       else
          q = -0.5_DP * (b + sign(sqrt(disc),b))
       endif

       if (b < 0.0_DP) then
          step = q / a
       else
          step = c / q
       end if

       ! Update vector xvec and matx
       ! agrecocmplx
       call data_functions_axpy(xvec,dirn,step)
       call data_functions_axpy(matx,yvec,step)

       ! ndmh: avoid parallel desynchronisation of xvec and matx
       ! agrecocmplx
       if (loc_cmplx) then
          call comms_bcast(pub_root_proc_id,xvec%z)
          call comms_bcast(pub_root_proc_id,matx%z)
       else
          call comms_bcast(pub_root_proc_id,xvec%d)
          call comms_bcast(pub_root_proc_id,matx%d)
       end if

       ! Renormalise
       call internal_matvec(metx,met,xvec)
       ! agrecocmplx
       temp_coef = data_functions_dot(xvec,metx)
       if (temp_coef%iscmplx) then
          norm = real(temp_coef%z,kind=DP)
       else
          norm = temp_coef%d
       end if
       normfac = 1.0_DP / sqrt(norm)

       call comms_bcast(pub_root_proc_id,normfac)

       ! agrecocmplx
       call data_functions_scale(xvec,normfac)
       call data_functions_scale(matx,normfac)
       call data_functions_scale(metx,normfac)

       ! Re-estimate eigenvalue
       oldeval = eval
       temp_coef = data_functions_dot(metx,matx)
       if (temp_coef%iscmplx) then
          eval = real(temp_coef%z,kind=DP)
       else
          eval = temp_coef%d
       end if

       ! ndmh: prevent WIN32 desync as above
       call comms_bcast(pub_root_proc_id,eval)

       ! Check convergence (to better than tol)
       if (abs(eval - oldeval) < local_tol) then
          conv = conv + 1
       else
          conv = 0
       end if
       ! ndmh: prevent WIN32 desync as above
       call comms_bcast(pub_root_proc_id,conv)
       if (conv > cgreset .or. abs(gdotg) < epsilon(1.0_DP)) exit

    end do

    if (loc_warnings.and.((iter==maxiter).or.(utils_isnan(eval)))) then
       if (pub_on_root) then
          write(stdout,'(a)') 'WARNING in sparse_extremal_eigenvalue: &
               &convergence did not reach set threshold:'
          write(stdout,'(a,i5,3f20.14)') 'iter = ',iter,' eval = ',eval, &
               ' oldeval = ',oldeval,' step = ',step,' gdotg = ',gdotg
       end if
    end if

    if(present(evec)) then
       call data_functions_copy(evec,xvec)
       if (loc_cmplx) then
          call comms_bcast(pub_root_proc_id,evec%z)
       else
          call comms_bcast(pub_root_proc_id,evec%d)
       end if
    endif

    ! Deallocate workspace
    call internal_dealloc_work


  contains

    !==========================================================================!
    ! This subroutine performs a matrix-vector product between a block-sparse  !
    ! matrix and a single dense vector.                                        !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   yvec (output) : the result vector                                      !
    !   amat (input)  : the block sparse matrix                                !
    !   xvec (input)  : the original vector                                    !
    !--------------------------------------------------------------------------!
    ! Written by Peter Haynes, March 2004.                                     !
    ! Revised for distributed data, June 2006.                                 !
    ! Modified for SPAM3 by Nicholas Hine, June 2009.                          !
    ! Modified by Andrea Greco to allow compatibility with complex case,       !
    ! February 2016.                                                           !
    !==========================================================================!

    subroutine internal_matvec(yvec,amat,xvec)

      implicit none

      ! Arguments
      ! agrecocmplx
      type(FUNCTIONS), intent(inout) :: yvec  ! result vector
      type(SPAM3), intent(in) :: amat         ! block sparse matrix
      type(FUNCTIONS), intent(in) :: xvec     ! input vector

      ! Local variables
      integer :: alib        ! Library entry for amat
      integer :: blks        ! Identifier for row blocking scheme
      integer :: iblk        ! Block-column counter
      integer :: loc_iblk    ! Block-column counter local to this proc
      integer :: jblk        ! Block-row counter
      integer :: idx         ! Index
      integer :: ielems      ! Number of elems in block-column
      integer :: jelems      ! Number of elems in block-row
      integer :: ielem       ! Row elem counter for input vector
      integer :: jelem       ! Row elem counter for result vector
      integer :: seg         ! Segment index
      integer :: seg_start   ! Start position of this segment in the index
      type(PARAL_INFO), pointer :: col_par ! Parallel strategy for columns.

      ! BLAS subroutines
      external :: zgemv, dgemv

      ! Get library entry
      alib = amat%lib
      blks = library(alib)%row_blks

      ! rc2013: assign parallel strategy from library structure
      col_par => library(alib)%col_par

      ! Set result vector to zero
      ! agrecocmplx
      call data_set_to_zero(yvec)

      ! The product to be formed is:
      !  y  = A   x
      !   j    ji  i

      ! Loop over segments of the matrix on this proc
!$OMP PARALLEL DO NUM_THREADS(pub_threads_max) DEFAULT (none) &
!$OMP PRIVATE(seg,ielems,ielem,jelems,jelem,seg_start,loc_iblk,iblk, &
!$OMP      idx,jblk) &
!$OMP SHARED(alib,amat,blks,library,pub_total_num_procs,xvec,yvec, &
!$OMP      col_par,pub_my_proc_id, &
!$OMP      pub_threads_max)
      do seg=0,pub_total_num_procs-1

         if (library(alib)%seg_info(s_type,seg)==SEG_DENSE) then

            ! Do matrix-vector multiplication
            ielems = col_par%num_elems_on_proc(pub_my_proc_id,blks)
            ielem = col_par%first_elem_on_proc(pub_my_proc_id,blks)
            jelems = col_par%num_elems_on_proc(seg,blks)
            jelem = col_par%first_elem_on_proc(seg,blks)

            ! agrecocmplx: complex case
            if (amat%iscmplx) then
               call zgemv('N',jelems,ielems,(1.0_DP,0.0_DP), &
                    amat%zmtx(library(alib)%seg_info(s_ptr,seg)), &
                    jelems, xvec%z(ielem), 1, &
                    (1.0_DP,0.0_DP), yvec%z(jelem), 1)
            ! real case
            else
               call dgemv('N',jelems,ielems,1.0_DP, &
                    amat%dmtx(library(alib)%seg_info(s_ptr,seg)), &
                    jelems, xvec%d(ielem), 1, &
                    1.0_DP, yvec%d(jelem), 1)

            end if

         else if (library(alib)%seg_info(s_type,seg)==SEG_SPARSE) then

            ! Loop over block-columns of matrix and input vector on this proc
            seg_start = library(alib)%seg_info(s_idx,seg)
            loc_iblk = 0
            do iblk=col_par%my_first_blk,col_par%my_last_blk
               loc_iblk = loc_iblk + 1

               ! Number of elems and first elem in this block-column
               ielems = col_par%num_elems_on_atom(iblk,blks)
               ielem = col_par%first_elem_on_atom(iblk,blks)

               ! Loop over block-rows jblk in block-col iblk of matrix
               do idx=library(alib)%blk_idx(seg_start+loc_iblk-1), &
                    library(alib)%blk_idx(seg_start+loc_iblk)-1
                  jblk = library(alib)%blk_idx(idx)

                  ! Number of elems and first elem in this block-row
                  jelems = col_par%num_elems_on_atom(jblk,blks)
                  jelem = col_par%first_elem_on_atom(jblk,blks)

                  ! Do block multiplication
                  ! agrecocmplx: complex case
                  if (amat%iscmplx) then
                     call zgemv('N',jelems,ielems,(1.0_DP,0.0_DP), &
                          amat%zmtx(library(alib)%blk_ptr(idx)), &
                          jelems, xvec%z(ielem), 1, &
                          (1.0_DP,0.0_DP), yvec%z(jelem), 1)
                  ! real case
                  else
                     call dgemv('N',jelems,ielems,1.0_DP, &
                          amat%dmtx(library(alib)%blk_ptr(idx)), &
                          jelems, xvec%d(ielem), 1, &
                          1.0_DP, yvec%d(jelem), 1)
                  end if

               end do  ! Loop over block-rows jblk in block-col iblk of matrix

            end do  ! Loop over block-columns of matrix and input vector

         end if  ! seg_info(s_type,seg)

      end do  ! seg
!$OMP END PARALLEL DO

      ! Sum contributions from all procs
      ! agrecocmplx
      if (loc_cmplx) then
         call comms_reduce('SUM',yvec%z)
      else
         call comms_reduce('SUM',yvec%d)
      end if

    end subroutine internal_matvec


    !==========================================================================!
    ! This subroutine allocates workspace for the parent subroutine.           !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   None                                                                   !
    !--------------------------------------------------------------------------!
    ! Written by Peter Haynes, March 2004.                                     !
    ! Revised for distributed data, June 2006.                                 !
    ! Modified by Andrea Greco to allow compatibility with complex case,       !
    ! February 2016.                                                           !
    !==========================================================================!

    subroutine internal_alloc_work

      ! agrecocmplx
      use datatypes, only: data_functions_alloc

      implicit none

      ! agrecocmplx
      ! Allocate workspace using appropriate functions
      call data_functions_alloc(xvec,n,iscmplx=loc_cmplx)

      call data_functions_alloc(yvec,n,iscmplx=loc_cmplx)

      call data_functions_alloc(grad,n,iscmplx=loc_cmplx)

      call data_functions_alloc(dirn,n,iscmplx=loc_cmplx)

      call data_functions_alloc(matx,n,iscmplx=loc_cmplx)

      call data_functions_alloc(metx,n,iscmplx=loc_cmplx)

    end subroutine internal_alloc_work


    !==========================================================================!
    ! This subroutine deallocates workspace for the parent subroutine.         !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   None                                                                   !
    !--------------------------------------------------------------------------!
    ! Written by Peter Haynes, March 2004.                                     !
    ! Revised for distributed data, June 2006.                                 !
    ! Modified by Andrea Greco to allow compatibility with complex case,       !
    ! February 2016.                                                           !
    !==========================================================================!

    subroutine internal_dealloc_work

      ! agrecocmplx
      use datatypes, only: data_functions_dealloc

      implicit none

      ! agrecocmplx
      ! Deallocate workspace using appropriate routines
      call data_functions_dealloc(metx)

      call data_functions_dealloc(matx)

      call data_functions_dealloc(dirn)

      call data_functions_dealloc(grad)

      call data_functions_dealloc(yvec)

      call data_functions_dealloc(xvec)

    end subroutine internal_dealloc_work

  end subroutine sparse_extremal_eigenvalue


  !============================================================================!
  ! This subroutine enforces the sparsity pattern of a given sparse matrix, in !
  ! the cases where some or all of its segments are declared 'dense'.          !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   mat  (inout) : The sparse matrix                                         !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, June 2009.                                       !
  !============================================================================!

  subroutine sparse_enforce_sparsity(mat)

    use comms, only: pub_my_proc_id, pub_total_num_procs
!$  use rundat, only: pub_threads_max
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(SPAM3), intent(inout) :: mat   ! The sparse matrix

    ! Local variables
    integer :: lib             ! Library pointers for x and y
    integer :: row_blks        ! Identifier for column blocking scheme
    integer :: col_blks        ! Identifier for row blocking scheme
    integer :: seg             ! Segment index
    integer :: seg_start       ! Start position of this segment in the index
    integer :: iblk            ! Block-column counter
    integer :: loc_iblk        ! Block-column counter local to this proc
    integer :: jblk            ! Block-row counter
    integer :: iidx            ! Index, Pointer
    integer :: ptr,bptr        ! Matrix and block pointers
    integer :: col_nzb         ! Number of nonzero blocks in this column
    integer :: ielems,jelems   ! Numbers of elements in each block
    integer :: ielem           ! ELement in block loop counter
    integer :: nrows           ! Total number of rows in this segment
    integer :: ierr            ! Error flag
    ! rc2013: local row and column parallel strategy pointers
    type(PARAL_INFO), pointer :: row_par
    type(PARAL_INFO), pointer :: col_par

    real(kind=DP), allocatable :: dblk_col(:)
    complex(kind=DP), allocatable :: zblk_col(:)

    ! Get library entries
    lib = mat%lib
    row_blks = library(lib)%row_blks
    col_blks = library(lib)%col_blks
    ! rc2013: get pointers to row and column parallel strategies
    row_par => library(lib)%row_par
    col_par => library(lib)%col_par

    ! Loop over segments
!$OMP PARALLEL NUM_THREADS(pub_threads_max) DEFAULT(none) &
!$OMP PRIVATE(seg,nrows,loc_iblk,iblk,ielems,col_nzb,iidx,jblk,jelems, &
!$OMP      ptr,bptr,ielem,seg_start,zblk_col,dblk_col,ierr) &
!$OMP SHARED(mat,library,lib,row_blks,col_blks,pub_total_num_procs, &
!$OMP      col_par,row_par, pub_my_proc_id,pub_threads_max)

    ! Allocate temporary storage for block column
    if (mat%iscmplx) then
       allocate(zblk_col(maxval(row_par%num_elems_on_proc(:,row_blks)) * &
            maxval(col_par%num_elems_on_atom(:,col_blks))),stat=ierr)
       call utils_alloc_check('sparse_enforce_sparsity','zblk_col',ierr)
       zblk_col(:) = (0.0_DP,0.0_DP)
    else
       allocate(dblk_col(maxval(row_par%num_elems_on_proc(:,row_blks)) * &
            maxval(col_par%num_elems_on_atom(:,col_blks))),stat=ierr)
       call utils_alloc_check('sparse_enforce_sparsity','dblk_col',ierr)
       dblk_col(:) = 0.0_DP
    end if

!$OMP DO
    do seg=0,pub_total_num_procs-1

       ! Only anything to do if this is a dense segment
       if (library(lib)%seg_info(s_type,seg)/=SEG_DENSE) cycle

       seg_start = library(lib)%seg_info(s_idx,seg)
       nrows = row_par%num_elems_on_proc(seg,row_blks)

       ! Loop over block columns on this proc
       loc_iblk = 0
       do iblk=col_par%my_first_blk,col_par%my_last_blk
          loc_iblk = loc_iblk + 1

          ielems = col_par%num_elems_on_atom(iblk,col_blks)
          if (ielems==0) cycle

          ! Find number of nonzero blocks in this column of this segment
          col_nzb = library(lib)%blk_idx(seg_start+loc_iblk) - &
               library(lib)%blk_idx(seg_start+loc_iblk-1)
          ! Move on if all the blocks are nonzero (nothing to do)
          if (col_nzb == col_par%num_atoms_on_proc(seg) .and. col_nzb.gt.0) cycle

          ! Otherwise, go through and copy the blocks in the index to blk_col
          do iidx=library(lib)%blk_idx(seg_start+loc_iblk-1), &
               library(lib)%blk_idx(seg_start+loc_iblk)-1
             jblk = library(lib)%blk_idx(iidx)
             jelems = row_par%num_elems_on_atom(jblk,row_blks)
             ptr = library(lib)%blk_ptr(iidx)
             bptr = row_par%first_elem_on_atom(jblk,row_blks) - &
                  row_par%first_elem_on_proc(seg,row_blks) + 1

             if (mat%iscmplx) then
                do ielem=1,ielems
                   zblk_col(bptr:bptr+jelems-1) = mat%zmtx(ptr:ptr+jelems-1)
                   ptr = ptr + nrows
                   bptr = bptr + nrows
                end do
             else
                do ielem=1,ielems
                   dblk_col(bptr:bptr+jelems-1) = mat%dmtx(ptr:ptr+jelems-1)
                   ptr = ptr + nrows
                   bptr = bptr + nrows
                end do
             end if

          end do

          ! Copy blk_col back to mtx
          ptr = library(lib)%seg_info(s_ptr,seg) + nrows * &
               (col_par%first_elem_on_atom(iblk,col_blks) - &
               col_par%first_elem_on_proc(pub_my_proc_id,col_blks))
          bptr = 1
          if (mat%iscmplx) then
             do ielem=1,ielems
                mat%zmtx(ptr:ptr+nrows-1) = zblk_col(bptr:bptr+nrows-1)
                ptr = ptr + nrows
                bptr = bptr + nrows
             end do
          else
             do ielem=1,ielems
                mat%dmtx(ptr:ptr+nrows-1) = dblk_col(bptr:bptr+nrows-1)
                ptr = ptr + nrows
                bptr = bptr + nrows
             end do
          end if

          ! Reset nonzero blocks of blk_col to zero
          do iidx=library(lib)%blk_idx(seg_start+loc_iblk-1), &
               library(lib)%blk_idx(seg_start+loc_iblk)-1
             jblk = library(lib)%blk_idx(iidx)
             jelems = row_par%num_elems_on_atom(jblk,row_blks)
             bptr = row_par%first_elem_on_atom(jblk,row_blks) - &
                  row_par%first_elem_on_proc(seg,row_blks) + 1

             if (mat%iscmplx) then
                do ielem=1,ielems
                   zblk_col(bptr:bptr+jelems-1) = (0.0_DP,0.0_DP)
                   bptr = bptr + nrows
                end do
             else
                do ielem=1,ielems
                   dblk_col(bptr:bptr+jelems-1) = 0.0_DP
                   bptr = bptr + nrows
                end do
             end if
          end do

       end do  ! iblk

    end do  ! seg
!$OMP END DO

    if (mat%iscmplx) then
       deallocate(zblk_col,stat=ierr)
       call utils_dealloc_check('sparse_enforce_sparsity','zblk_col',ierr)
    else
       deallocate(dblk_col,stat=ierr)
       call utils_dealloc_check('sparse_enforce_sparsity','dblk_col',ierr)
    end if

!$OMP END PARALLEL

  end subroutine sparse_enforce_sparsity

  !============================================================================!
  ! This subroutine calculates the outer product matrix uv of two vectors. The !
  ! input matrix sparsity pattern is respected in the final answer, so lots of !
  ! nonzero terms of the outer product may get skipped. It is the user's       !
  ! responsibility to ensure that a suitable sparsity pattern is passed in     !
  ! that will not result in unwanted truncation.                               !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !  mat (inout) : The sparse matrix                                           !
  !  uvec (in)   : The real row vector u                                       !
  !  vvec (in)   : The real column vector v                                    !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, August 2015.                                     !
  !============================================================================!

  subroutine sparse_outer_product_real(mat,uvec,vvec)

    use comms, only: pub_my_proc_id, pub_total_num_procs
    use rundat, only: pub_threads_max

    implicit none

    ! Arguments
    type(SPAM3), intent(inout) :: mat   ! The sparse matrix
    real(kind=DP), intent(in) :: uvec(:),vvec(:)

    ! Local variables
    integer :: lib             ! Library pointer
    integer :: row_blks        ! Identifier for column blocking scheme
    integer :: col_blks        ! Identifier for row blocking scheme
    integer :: seg             ! Segment index
    integer :: seg_type        ! Segment type for this segment
    integer :: seg_start       ! Start position of this segment in the index
    integer :: iblk            ! Block-column counter
    integer :: loc_iblk        ! Block-column counter local to this proc
    integer :: jblk            ! Block-row counter
    integer :: iidx            ! Index counter
    integer :: ptr             ! Matrix data pointer
    integer :: nrows           ! Number of rows to get to start of next col
    integer :: ielems,jelems   ! Numbers of row/col elements in each block
    integer :: ielem0,jelem0   ! First row/col elements in block
    integer :: ielem,jelem     ! ELement in row/col block loop counter
    ! rc2013: local pointers to row and column parallel strategies
    type(PARAL_INFO), pointer :: row_par
    type(PARAL_INFO), pointer :: col_par

    ! Get library entries
    lib = mat%lib
    row_blks = library(lib)%row_blks
    col_blks = library(lib)%col_blks

    ! rc2013: specify row and column par pointers
    row_par => library(lib)%row_par
    col_par => library(lib)%col_par

    ! Loop over segments
!$OMP PARALLEL NUM_THREADS(pub_threads_max) DEFAULT(none) &
!$OMP PRIVATE(seg,loc_iblk,iblk,ielems,ielem0,iidx,jblk,jelems,jelem0, &
!$OMP      ptr,ielem,jelem,seg_start,seg_type,nrows) &
!$OMP SHARED(mat,library,lib,row_blks,col_blks,pub_total_num_procs, &
!$OMP      col_par,row_par,uvec,vvec, pub_my_proc_id,pub_threads_max)

!$OMP DO
    do seg=0,pub_total_num_procs-1

       seg_start = library(lib)%seg_info(s_idx,seg)
       seg_type = library(lib)%seg_info(s_type,seg)
       if (seg_type == SEG_DENSE) nrows = row_par%num_elems_on_proc(seg,row_blks)
       if (seg_type == SEG_BLANK) cycle

       ! Loop over block columns on this proc
       loc_iblk = 0
       do iblk=col_par%my_first_blk,col_par%my_last_blk
          loc_iblk = loc_iblk + 1

          ! Find range of colunns in this block
          ielems = col_par%num_elems_on_atom(iblk,col_blks)
          ielem0 = col_par%first_elem_on_atom(iblk,col_blks)
          if (ielems==0) cycle

          ! Cycle through nonzero blocks in this block column
          do iidx=library(lib)%blk_idx(seg_start+loc_iblk-1), &
               library(lib)%blk_idx(seg_start+loc_iblk)-1
             jblk = library(lib)%blk_idx(iidx)
             jelems = row_par%num_elems_on_atom(jblk,row_blks)
             jelem0 = row_par%first_elem_on_atom(jblk,row_blks)

             ! Set jump at end of column to number of rows in block if sparse
             if (seg_type == SEG_SPARSE) nrows = jelems

             ! Now actually form the outer product elements in this block
             if (mat%iscmplx) then
                do ielem=1,ielems
                   ptr = library(lib)%blk_ptr(iidx) + nrows*(ielem-1)
                   do jelem=1,jelems
                      mat%zmtx(ptr) = cmplx(uvec(jelem0+jelem-1) * &
                           vvec(ielem0+ielem-1),0.0_DP,kind=DP)
                      ptr = ptr + 1
                   end do  ! jelem
                end do   ! ielem
             else
                do ielem=1,ielems
                   ptr = library(lib)%blk_ptr(iidx) + nrows*(ielem-1)
                   do jelem=1,jelems
                      mat%dmtx(ptr) = uvec(jelem0+jelem-1) * &
                           vvec(ielem0+ielem-1)
                      ptr = ptr + 1
                   end do  ! jelem
                end do   ! ielem
              end if  ! iscmplx
          end do  ! iidx

       end do  ! iblk

    end do  ! seg
!$OMP END DO

!$OMP END PARALLEL

  end subroutine sparse_outer_product_real


  !============================================================================!
  ! This subroutine calculates the outer product matrix uv of two vectors. The !
  ! input matrix sparsity pattern is respected in the final answer, so lots of !
  ! nonzero terms of the outer product may get skipped. It is the user's       !
  ! responsibility to ensure that a suitable sparsity pattern is passed in     !
  ! that will not result in unwanted truncation.                               !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !  mat (inout) : The sparse matrix                                           !
  !  uvec (in)   : The complex row vector u                                    !
  !  vvec (in)   : The complex column vector v                                 !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, August 2015.                                     !
  !============================================================================!

  subroutine sparse_outer_product_complex(mat,uvec,vvec)

    use comms, only: pub_my_proc_id, pub_total_num_procs
    use parallel_strategy, only: PARAL_INFO
    use rundat, only: pub_threads_max
    use utils, only: utils_assert

    implicit none

    ! Arguments
    type(SPAM3), intent(inout) :: mat   ! The sparse matrix
    complex(kind=DP), intent(in) :: uvec(:),vvec(:)

    ! Local variables
    integer :: lib             ! Library pointer
    integer :: row_blks        ! Identifier for column blocking scheme
    integer :: col_blks        ! Identifier for row blocking scheme
    integer :: seg             ! Segment index
    integer :: seg_type        ! Segment type for this segment
    integer :: seg_start       ! Start position of this segment in the index
    integer :: iblk            ! Block-column counter
    integer :: loc_iblk        ! Block-column counter local to this proc
    integer :: jblk            ! Block-row counter
    integer :: iidx            ! Index counter
    integer :: ptr             ! Matrix data pointer
    integer :: nrows           ! Number of rows to get to start of next col
    integer :: ielems,jelems   ! Numbers of row/col elements in each block
    integer :: ielem0,jelem0   ! First row/col elements in block
    integer :: ielem,jelem     ! ELement in row/col block loop counter
    type(PARAL_INFO), pointer :: row_par ! Parallel strategy for rows.
    type(PARAL_INFO), pointer :: col_par ! Parallel strategy for columns.

    call utils_assert(mat%iscmplx, 'Error in sparse_outer_product_complex: &
         &routine is only valid for complex matrices ')

    ! Get library entries
    lib = mat%lib
    row_blks = library(lib)%row_blks
    col_blks = library(lib)%col_blks

    ! rc2013: specify row and column par pointers
    row_par => library(lib)%row_par
    col_par => library(lib)%col_par

    ! Loop over segments
!$OMP PARALLEL NUM_THREADS(pub_threads_max) DEFAULT(none) &
!$OMP PRIVATE(seg,loc_iblk,iblk,ielems,ielem0,iidx,jblk,jelems,jelem0, &
!$OMP      ptr,ielem,jelem,seg_start,seg_type,nrows) &
!$OMP SHARED(mat,library,lib,row_blks,col_blks,pub_total_num_procs, &
!$OMP      col_par,row_par,uvec,vvec, pub_my_proc_id,pub_threads_max)

!$OMP DO
    do seg=0,pub_total_num_procs-1

       seg_start = library(lib)%seg_info(s_idx,seg)
       seg_type = library(lib)%seg_info(s_type,seg)
       if (seg_type == SEG_DENSE) nrows = row_par%num_elems_on_proc(seg,row_blks)
       if (seg_type == SEG_BLANK) cycle

       ! Loop over block columns on this proc
       loc_iblk = 0
       do iblk=col_par%my_first_blk,col_par%my_last_blk
          loc_iblk = loc_iblk + 1

          ! Find range of colunns in this block
          ielems = col_par%num_elems_on_atom(iblk,col_blks)
          ielem0 = col_par%first_elem_on_atom(iblk,col_blks)
          if (ielems==0) cycle

          ! Cycle through nonzero blocks in this block column
          do iidx=library(lib)%blk_idx(seg_start+loc_iblk-1), &
               library(lib)%blk_idx(seg_start+loc_iblk)-1
             jblk = library(lib)%blk_idx(iidx)
             jelems = row_par%num_elems_on_atom(jblk,row_blks)
             jelem0 = row_par%first_elem_on_atom(jblk,row_blks)

             ! Set jump at end of column to number of rows in block if sparse
             if (seg_type == SEG_SPARSE) nrows = jelems

             ! Now actually form the outer product elements in this block
             if (mat%iscmplx) then
                do ielem=1,ielems
                   ptr = library(lib)%blk_ptr(iidx) + nrows*(ielem-1)
                   do jelem=1,jelems
                      mat%zmtx(ptr) = uvec(jelem0+jelem-1) * &
                           conjg(vvec(ielem0+ielem-1))
                      ptr = ptr + 1
                   end do  ! jelem
                end do   ! ielem
             else
                ! Not sure this would make sense. Feel free to add if it does.
             end if  ! iscmplx
          end do  ! iidx

       end do  ! iblk

    end do  ! seg
!$OMP END DO

!$OMP END PARALLEL

  end subroutine sparse_outer_product_complex


  !==========================================================================!
  ! This subroutine initialises the inverse of a sparse matrix so that it is !
  ! ready to be passed to sparse_hotelling_invert and inverted using the     !
  ! Hotelling Algorithm.                                                     !
  ! Based on the theory in T. Ozaki, Phys. Rev. B., vol 64, page 195110.     !
  !--------------------------------------------------------------------------!
  ! Arguments:                                                               !
  !   inv_mat  (input) : The inverse to initialise                           !
  !   mat      (input) : The matrix to be inverted                           !
  !--------------------------------------------------------------------------!
  ! Originally written by Chris-Kriton Skylaris on 31/7/2003.                !
  ! Modified to support SPAM 2 by Peter Haynes on 21/7/2004.                 !
  ! Adapted for SPAM3 matrices by Nicholas Hine, June 2009.                  !
  ! Moved to sparse_mod by Nicholas Hine, September 2009.                    !
  !==========================================================================!

  subroutine sparse_hotelling_init(inv_mat,mat)

    use comms, only: pub_on_root, comms_reduce, pub_my_proc_id
    use constants, only: DP, stdout, cmplx_0
    use utils, only: utils_alloc_check, utils_dealloc_check, &
        utils_assert

    implicit none

    ! Arguments
    type(SPAM3), intent(inout) :: inv_mat
    type(SPAM3), intent(in)  :: mat

    ! Local Variables
    real(kind=DP) :: sigma, sum_col
    real(kind=DP), allocatable, dimension(:) :: drow
    complex(kind=DP), allocatable, dimension(:) :: zrow
    integer :: nn, row, col
    integer :: ierr
    integer :: ilib, col_blks
    type(PARAL_INFO), pointer :: col_par

    ! agrecocmplx: either both real or both complex
    call utils_assert((mat%iscmplx .eqv. inv_mat%iscmplx), &
         'Error in sparse_hotelling_init: incompatible argument types.')

    ! Allocate workspace
    ilib = mat%lib
    nn = library(ilib)%nrows
    col_blks = library(ilib)%col_blks
    col_par => library(ilib)%col_par

    if (library(ilib)%row_blks /= col_blks) then
       ! Error
    end if
    if (mat%iscmplx) then
       allocate(zrow(nn),stat=ierr)
       call utils_alloc_check('sparse_hotelling_init','zrow',ierr)
       zrow = cmplx_0
    else
       allocate(drow(nn),stat=ierr)
       call utils_alloc_check('sparse_hotelling_init','drow',ierr)
       drow = 0.0_DP
    end if

    sigma = 0.0_DP

    ! Loop over all cols of S on this processor
    do col=col_par%first_elem_on_proc(pub_my_proc_id,col_blks), &
         col_par%first_elem_on_proc(pub_my_proc_id+1,col_blks)-1

       if (mat%iscmplx) then

          ! Get column of S
          call sparse_get_col(zrow,mat,col)

          ! Sum column of S
          sum_col = 0.0_DP
          do row=1,nn
             ! jd: Ignore 'inaccessible' elements in dense blocks, whose values
             !     are undefined.
             if(sparse_element_exists(mat,row,col)) then
                sum_col = sum_col + abs(zrow(row))
             end if
          end do
          sigma = max(sigma,sum_col)

          ! Clear column of S
          call sparse_clr_col(zrow,mat,col)

       else

          ! Get column of S
          call sparse_get_col(drow,mat,col)

          ! Sum column of S
          sum_col = 0.0_DP
          do row=1,nn
             sum_col = sum_col + abs(drow(row))
          end do
          sigma = max(sigma,sum_col)

          ! Clear column of S
          call sparse_clr_col(drow,mat,col)

       end if

    end do

    call comms_reduce('MAX',sigma)

    ! Deallocate workspace
    if (mat%iscmplx) then
       deallocate(zrow,stat=ierr)
       call utils_dealloc_check('sparse_hotelling_init','zrow',ierr)
    else
       deallocate(drow,stat=ierr)
       call utils_dealloc_check('sparse_hotelling_init','drow',ierr)
    end if

    ! Copy mat into inv_mat
    call sparse_copy(inv_mat,mat)

    if (sigma > epsilon(1.0_DP)) then
       sigma = 1.0_DP / (sigma*sigma)
    else
       sigma = 0.001_DP
       if (pub_on_root) write(stdout,'(a)') &
            'WARNING in sparse_hotelling_init: zero overlap matrix'
    end if

    ! cks: scale by sigma
    call sparse_scale(inv_mat,sigma)

  end subroutine sparse_hotelling_init

  !==========================================================================!
  ! This subroutine applies successive quadratically convergent Hotelling    !
  ! iterations to improve an approximate inverse overlap matrix, based on    !
  ! the theory in the paper by T. Ozaki, Phys. Rev. B., vol 64, page 195110. !
  ! The iterations continue (unless they exceed the num_iter input           !
  ! parameter) until convergence to machine precision is reached - but       !
  ! taking into account limitations arising from the truncation of the       !
  ! inverse overlap matrix.                                                  !
  !--------------------------------------------------------------------------!
  ! Arguments:                                                               !
  !   inv_mat  (input) : The inverse to calculate                            !
  !   mat      (input) : The matrix to be inverted                           !
  !--------------------------------------------------------------------------!
  ! Originally written by Chris-Kriton Skylaris on 31/7/2003                 !
  ! Modified to support SPAM 2 by Peter Haynes on 21/7/2004.                 !
  ! Modified to test for convergence by Chris-Kriton Skylaris on 4/10/2004.  !
  ! Minor modifications for parallel SPAM 2 by Peter Haynes                  !
  ! Adapted to use SPAM3 matrices by Nicholas Hine, June 2009.               !
  ! Moved to sparse_mod by Nicholas Hine, September 2009.                    !
  !==========================================================================!

  subroutine sparse_hotelling_invert(inv_mat,mat,show_output, &
       max_resid_converged, num_iter, final_max_resid)

    use comms, only: pub_on_root
    use constants, only: DP, stdout
    use timer, only: timer_clock
    use utils, only: utils_assert

    implicit none

    ! Arguments
    type(SPAM3), intent(inout) :: inv_mat
    type(SPAM3), intent(in) :: mat
    real(kind=DP), intent(in) :: max_resid_converged
    logical, intent(in) :: show_output
    integer, intent(in) :: num_iter
    real(kind=DP), intent(out), optional :: final_max_resid

    ! Local Variables
    integer :: iter
    logical :: quit_early  ! pdh: flag
    real(kind=DP) :: max_resid  ! maximum value of residual
    real(kind=DP) :: frob_norm
    real(kind=DP) :: previous_frob_norm
    type(SPAM3) :: mim,inv_tmp

    ! Start timer
    call timer_clock('sparse_hotelling_invert',1)

    ! agrecocmplx: either both real or both complex
    call utils_assert((mat%iscmplx .eqv. inv_mat%iscmplx), &
         'Error in sparse_hotelling_invert: incompatible argument types.')

    ! Initialisations
    previous_frob_norm = huge(1.0_DP)
    quit_early = .false.

    ! Allocate workspace
    call sparse_create(inv_tmp,inv_mat)
    call sparse_create(mim,mat,inv_mat)

    ! Start of Hotelling iteration loop
    hotel_sparse_loop: do iter=1,num_iter

       ! MIM := M*(M^-1)
       call sparse_product(mim,mat,inv_mat)

       ! cks: MIM := I-M*(M^-1)
       call sparse_scale(mim,-1.0_DP,1.0_DP)

       ! cks: Calculate Frobenius norm of MIM -->0
       frob_norm = sparse_rms_element(mim) * &
            sqrt(sparse_num_element(mim))

       ! cks: Maximum element of residual I-M*M_n^-1
       max_resid = sparse_max_abs_element(mim)

       if (show_output.and.pub_on_root) then
          write(stdout,'(t12,i5,tr5,e16.8,tr3,e16.8)') iter, frob_norm, &
               max_resid
       end if

       ! cks: Test for convergence due to machine precision
       ! cks: or inverse overlap truncation
       if (frob_norm >= previous_frob_norm .or. &
            max_resid <= max_resid_converged) then
          quit_early = .true.
          exit
       else
          previous_frob_norm = frob_norm
       endif

       ! cks: MIM := I + MIM := 2I -M*(M^-1)
       call sparse_scale(mim,1.0_DP,1.0_DP)

       ! (M^-1) := (M^-1)*(2*I-M*(M^-1))
       call sparse_copy(inv_tmp,inv_mat)
       call sparse_product(inv_mat,inv_tmp,mim)

    end do hotel_sparse_loop

    if (present(final_max_resid)) final_max_resid = max_resid

    if (pub_on_root .and. .not. quit_early) write(stdout,*) &
         'WARNING: max Hotelling iterations exceeded. Last norm=',frob_norm

    ! Deallocate workspace
    call sparse_destroy(mim)
    call sparse_destroy(inv_tmp)

    ! Stop timer
    call timer_clock('sparse_hotelling_invert',2)

  end subroutine sparse_hotelling_invert


! These routines only need to exist if the code is being compiled with ScaLAPACK
#ifdef SCALAPACK

  !============================================================================!
  ! This subroutine converts a block sparse matrix in a SPAM3 into a dense     !
  ! array in the BLACS distributed format suitable for use with ScaLAPACK      !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   dest   (output) : The destination BLACS matrix (must be real)            !
  !   src    (input)  : The source SPAM3 matrix (must be real)                 !
  !   ld     (input)  : Leading dimension of the array dest                    !
  !   nc     (input)  : Number of columns of the array dest                    !
  !   desc   (input)  : BLACS descriptor for the array dest                    !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, October 2009.                                    !
  !============================================================================!

  subroutine sparse_spam3toblacs_real(dest,src,ld,nc,desc)

    use comms, only: comms_bcast, pub_my_proc_id, pub_on_root, &
         pub_total_num_procs
    use parallel_strategy, only: PARAL_INFO
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_assert

    implicit none

    ! Arguments
    integer,intent(in) :: ld,nc
    real(kind=DP),intent(out) :: dest(ld,nc)
    type(SPAM3),intent(in) :: src
    integer,intent(in) :: desc(50)

    ! LAPACK subroutine
    external :: pdelset

    ! Local Variables
    integer :: ierr
    integer :: srclib        ! Library entry for src
    integer :: row_blks      ! Identifier for column blocking scheme
    integer :: col_blks      ! Identifier for row blocking scheme
    integer :: iblk,jblk     ! Block counters
    integer :: loc_iblk      ! Local counter for iblk
    integer :: idx           ! Index
    integer :: ptr           ! Pointer
    integer :: ielemonat     ! elem in atom block iblk counter
    integer :: jelemonat     ! elem in atom block jblk counter
    integer :: ielem,jelem   ! elem counters
    integer :: ielems,jelems ! Numbers of element rows/cols in this block
    integer :: nrows         ! Total number of rows in this block/segment
    integer :: proc          ! Proc loop counter
    integer :: seg           ! Segment index
    integer :: seg_type      ! Start position of this segment in the index
    integer :: seg_start     ! Start position of this segment in the index
    real(kind=DP), allocatable :: buffer(:,:),loc_buffer(:,:)
    type(PARAL_INFO), pointer :: col_par ! Parallel strategy for columns.
    type(PARAL_INFO), pointer :: row_par ! Parallel strategy for columns.

    ! Get library entry for src
    srclib = src%lib
    row_blks = library(srclib)%row_blks
    col_blks = library(srclib)%col_blks

    ! rc2013: specify row and column par pointers
    row_par => library(srclib)%row_par
    col_par => library(srclib)%col_par

    ! Check arguments
    call utils_assert(.not. src%iscmplx, &
         'Error in sparse_spam3toblacs_real: real matrices only.')

    ! Allocate buffers
    nrows = library(srclib)%nrows
    allocate(buffer(nrows,maxval(col_par%num_elems_on_proc(:,col_blks))),stat=ierr)
    call utils_alloc_check('sparse_spam3toblacs_real','buffer',ierr)
    allocate(loc_buffer(nrows,col_par%num_elems_on_proc(pub_my_proc_id,col_blks)), &
         stat=ierr)
    call utils_alloc_check('sparse_spam3toblacs_real','loc_buffer',ierr)

    ! Zero destination matrix and buffer
    dest(:,:) = 0.0_DP
    buffer(:,:) = 0.0_DP
    loc_buffer(:,:) = 0.0_DP

    ! Loop over the segments of src on this proc
    do seg=0,pub_total_num_procs-1

       ! Cycle if this segment is blank
       seg_type = library(srclib)%seg_info(s_type,seg)
       if (seg_type == SEG_BLANK) cycle
       if (seg_type == SEG_DENSE) nrows = row_par%num_elems_on_proc(seg,row_blks)

       ! Loop over block-columns of src on this proc
       seg_start = library(srclib)%seg_info(s_idx,seg)
       loc_iblk = 0
       do iblk=col_par%my_first_blk,col_par%my_last_blk
          loc_iblk = loc_iblk + 1
          ielem = col_par%first_elem_on_atom(iblk,col_blks) - &
               col_par%first_elem_on_proc(pub_my_proc_id,col_blks) + 1
          ielems = col_par%num_elems_on_atom(iblk,col_blks)

          ! Loop over block-rows jblk in column iblk of src
          do idx=library(srclib)%blk_idx(seg_start+loc_iblk-1), &
               library(srclib)%blk_idx(seg_start+loc_iblk)-1
             jblk = library(srclib)%blk_idx(idx)
             jelem = row_par%first_elem_on_atom(jblk,row_blks)
             jelems = row_par%num_elems_on_atom(jblk,row_blks)

             if (seg_type == SEG_SPARSE) nrows = jelems

             ! Copy block
             ptr = library(srclib)%blk_ptr(idx)
             do ielemonat=0,ielems-1      ! Element col in block
                do jelemonat=0,jelems-1   ! Element row in block
                   loc_buffer(jelem+jelemonat,ielem+ielemonat) = &
                        src%dmtx(ptr+ielemonat*nrows+jelemonat)
                end do
             end do

          end do  ! Loop over block-rows in column iblk of src

       end do  ! Loop over block-columns of src

    end do  ! seg

    ! Broadcast results across all procs
    nrows = library(srclib)%nrows
    do proc=0,pub_total_num_procs-1

       ! Copy local buffer to buffer if proc is local proc
       if (proc==pub_my_proc_id) then
          ielems = col_par%num_elems_on_proc(pub_my_proc_id,col_blks)
          buffer(1:nrows,1:ielems) = loc_buffer(1:nrows,1:ielems)
       end if

       ! Broadcast this part of the matrix
       call comms_bcast(proc,buffer(1,1), &
            col_par%num_elems_on_proc(proc,col_blks)*library(srclib)%nrows)

       ! Loop over rows of dest, setting values
       do ielem=col_par%first_elem_on_proc(proc,col_blks), &
            col_par%first_elem_on_proc(proc+1,col_blks)-1
          ielems = ielem - col_par%first_elem_on_proc(proc,col_blks) + 1
          do jelem=1,nrows
             call pdelset(dest,jelem,ielem,desc,buffer(jelem,ielems))
          end do
       end do

    end do  ! Loop over procs

    ! Dellocate buffers
    deallocate(buffer,stat=ierr)
    call utils_dealloc_check('sparse_spam3toblacs_real','buffer',ierr)
    deallocate(loc_buffer,stat=ierr)
    call utils_dealloc_check('sparse_spam3toblacs_real','loc_buffer',ierr)

  end subroutine sparse_spam3toblacs_real

  !============================================================================!
  ! This subroutine converts a block sparse matrix in a SPAM3 into a dense     !
  ! array in the BLACS distributed format suitable for use with ScaLAPACK      !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   dest   (output) : The destination BLACS matrix (must be complex)         !
  !   src    (input)  : The source SPAM3 matrix (must be complex)              !
  !   ld     (input)  : Leading dimension of the array dest                    !
  !   nc     (input)  : Number of columns of the array dest                    !
  !   desc   (input)  : BLACS descriptor for the array dest                    !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, October 2009.                                    !
  !============================================================================!

  subroutine sparse_spam3toblacs_complex(dest,src,ld,nc,desc)

    use comms, only: comms_bcast, pub_my_proc_id, pub_total_num_procs
    use parallel_strategy, only: PARAL_INFO
    use utils, only: utils_assert, utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    integer,intent(in) :: ld,nc
    complex(kind=DP),intent(out) :: dest(ld,nc)
    type(SPAM3),intent(in) :: src
    integer,intent(in) :: desc(50)

    ! LAPACK subroutine
    external :: pzelset

    ! Local Variables
    integer :: ierr
    integer :: srclib        ! Library entry for src
    integer :: row_blks      ! Identifier for column blocking scheme
    integer :: col_blks      ! Identifier for row blocking scheme
    integer :: iblk,jblk     ! Block counters
    integer :: loc_iblk      ! Local counter for iblk
    integer :: idx           ! Index
    integer :: ptr           ! Pointer
    integer :: ielemonat     ! elem in atom block iblk counter
    integer :: jelemonat     ! elem in atom block jblk counter
    integer :: ielem,jelem   ! elem counters
    integer :: ielems,jelems ! Numbers of element rows/cols in this block
    integer :: nrows         ! Total number of rows in this block/segment
    integer :: proc          ! Proc loop counter
    integer :: seg           ! Segment index
    integer :: seg_type      ! Start position of this segment in the index
    integer :: seg_start     ! Start position of this segment in the index
    complex(kind=DP), allocatable :: buffer(:,:),loc_buffer(:,:)
    type(PARAL_INFO), pointer :: row_par ! Parallel strategy for rows.
    type(PARAL_INFO), pointer :: col_par ! Parallel strategy for columns.

    ! Get library entry for src
    srclib = src%lib
    row_blks = library(srclib)%row_blks
    col_blks = library(srclib)%col_blks

    ! rc2013: specify row and column par pointers
    row_par => library(srclib)%row_par
    col_par => library(srclib)%col_par

    ! Allocate buffers
    nrows = library(srclib)%nrows
    allocate(buffer(nrows,maxval(col_par%num_elems_on_proc(:,col_blks))),stat=ierr)
    call utils_alloc_check('sparse_spam3toblacs_complex','buffer',ierr)
    allocate(loc_buffer(nrows,col_par%num_elems_on_proc(pub_my_proc_id,col_blks)), &
         stat=ierr)
    call utils_alloc_check('sparse_spam3toblacs_complex','loc_buffer',ierr)

    ! Zero destination matrix and buffer
    dest(:,:) = cmplx(0.0_DP,0.0_DP,kind=DP)
    buffer(:,:) = cmplx(0.0_DP,0.0_DP,kind=DP)
    loc_buffer(:,:) = cmplx(0.0_DP,0.0_DP,kind=DP)

    ! Loop over the segments of src on this proc
    do seg=0,pub_total_num_procs-1

       ! Cycle if this segment is blank
       seg_type = library(srclib)%seg_info(s_type,seg)
       if (seg_type == SEG_BLANK) cycle
       if (seg_type == SEG_DENSE) nrows = row_par%num_elems_on_proc(seg,row_blks)

       ! Loop over block-columns of src on this proc
       seg_start = library(srclib)%seg_info(s_idx,seg)
       loc_iblk = 0
       do iblk=col_par%my_first_blk,col_par%my_last_blk
          loc_iblk = loc_iblk + 1
          ielem = col_par%first_elem_on_atom(iblk,col_blks) - &
               col_par%first_elem_on_proc(pub_my_proc_id,col_blks) + 1
          ielems = col_par%num_elems_on_atom(iblk,col_blks)

          ! Loop over block-rows jblk in column iblk of src
          do idx=library(srclib)%blk_idx(seg_start+loc_iblk-1), &
               library(srclib)%blk_idx(seg_start+loc_iblk)-1
             jblk = library(srclib)%blk_idx(idx)
             jelem = row_par%first_elem_on_atom(jblk,row_blks)
             jelems = row_par%num_elems_on_atom(jblk,row_blks)

             if (seg_type == SEG_SPARSE) nrows = jelems

             ! Copy block
             ptr = library(srclib)%blk_ptr(idx)
             if (src%iscmplx) then
                do ielemonat=0,ielems-1      ! Element col in block
                   do jelemonat=0,jelems-1   ! Element row in block
                      loc_buffer(jelem+jelemonat,ielem+ielemonat) = &
                           src%zmtx(ptr+ielemonat*nrows+jelemonat)
                   end do
                end do
             else
                do ielemonat=0,ielems-1      ! Element col in block
                   do jelemonat=0,jelems-1   ! Element row in block
                      loc_buffer(jelem+jelemonat,ielem+ielemonat) = &
                           cmplx(src%dmtx(ptr+ielemonat*nrows+jelemonat), &
                           0.0_DP,kind=DP)
                   end do
                end do
             end if

          end do  ! Loop over block-rows in column iblk of src

       end do  ! Loop over block-columns of src

    end do  ! seg

    ! Broadcast results across all procs
    nrows = library(srclib)%nrows
    do proc=0,pub_total_num_procs-1

       ! Copy local buffer to buffer if proc is local proc
       if (proc==pub_my_proc_id) then
          ielems = col_par%num_elems_on_proc(pub_my_proc_id,col_blks)
          buffer(1:nrows,1:ielems) = loc_buffer(1:nrows,1:ielems)
       end if

       ! Broadcast this part of the matrix
       call comms_bcast(proc,buffer(1,1), &
            col_par%num_elems_on_proc(proc,col_blks)*library(srclib)%nrows)

       ! Loop over rows of dest, setting values
       do ielem=col_par%first_elem_on_proc(proc,col_blks), &
            col_par%first_elem_on_proc(proc+1,col_blks)-1
          ielems = ielem - col_par%first_elem_on_proc(proc,col_blks) + 1
          do jelem=1,nrows
             call pzelset(dest,jelem,ielem,desc,buffer(jelem,ielems))
          end do
       end do

    end do  ! Loop over procs

    ! Dellocate buffers
    deallocate(buffer,stat=ierr)
    call utils_dealloc_check('sparse_spam3toblacs_complex','buffer',ierr)
    deallocate(loc_buffer,stat=ierr)
    call utils_dealloc_check('sparse_spam3toblacs_complex','loc_buffer',ierr)

  end subroutine sparse_spam3toblacs_complex

  !============================================================================!
  ! This subroutine converts a dense array in the BLACS distributed format     !
  ! suitable for use with ScaLAPACK into a SPAM3 block sparse matrix           !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   dest   (output) : The destination SPAM3 matrix (must be real)            !
  !   src    (input)  : The source BLACS dense matrix (real)                   !
  !   ld     (input)  : Leading dimension of the array dest                    !
  !   nc     (input)  : Number of columns of the array dest                    !
  !   desc   (input)  : BLACS descriptor for the array dest                    !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, February 2010.                                   !
  !============================================================================!

  subroutine sparse_blacstospam3_real(dest,src,ld,nc,desc)

    use comms, only: comms_reduce, pub_my_proc_id, pub_total_num_procs
    use parallel_strategy, only: PARAL_INFO
    use utils, only: utils_assert, utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    integer,intent(in) :: ld,nc
    type(SPAM3),intent(inout) :: dest
    real(kind=DP),intent(in) :: src(ld,nc)
    integer,intent(in) :: desc(50)

    ! ScaLAPACK subroutine
    external :: pdelget

    ! Local Variables
    integer :: ierr
    integer :: destlib       ! Library entry for dest
    integer :: row_blks      ! Identifier for column blocking scheme
    integer :: col_blks      ! Identifier for row blocking scheme
    integer :: iblk,jblk     ! Block counters
    integer :: loc_iblk      ! Local counter for iblk
    integer :: idx           ! Index
    integer :: ptr           ! Pointer
    integer :: ielemonat     ! elem in atom block iblk counter
    integer :: jelemonat     ! elem in atom block jblk counter
    integer :: ielem,jelem   ! elem counters
    integer :: ielems,jelems ! Numbers of element rows/cols in this block
    integer :: nrows         ! Total number of rows in this block/segment
    integer :: proc          ! Proc loop counter
    integer :: seg           ! Segment index
    integer :: seg_type      ! Start position of this segment in the index
    integer :: seg_start     ! Start position of this segment in the index
    real(kind=DP), allocatable :: buffer(:,:),loc_buffer(:,:)
    type(PARAL_INFO), pointer :: row_par ! Parallel strategy for rows.
    type(PARAL_INFO), pointer :: col_par ! Parallel strategy for columns.

    ! Get library entry for src
    destlib = dest%lib
    row_blks = library(destlib)%row_blks
    col_blks = library(destlib)%col_blks
    nrows = library(destlib)%nrows

    ! rc2013: specify row and column par pointers
    row_par => library(destlib)%row_par
    col_par => library(destlib)%col_par

    ! Check arguments
    call utils_assert(.not. dest%iscmplx, &
         'Error in sparse_blacstospam3_real: real matrices only.')

    ! Allocate buffers
    allocate(buffer(nrows,maxval(col_par%num_elems_on_proc(:,col_blks))),stat=ierr)
    call utils_alloc_check('sparse_blacstospam3_real','buffer',ierr)
    allocate(loc_buffer(nrows,col_par%num_elems_on_proc(pub_my_proc_id,col_blks)), &
         stat=ierr)
    call utils_alloc_check('sparse_blacstospam3_real','loc_buffer',ierr)

    ! Broadcast results across all procs
    nrows = library(destlib)%nrows
    do proc=0,pub_total_num_procs-1

       ! Zero buffer
       buffer(:,:) = 0.0_DP

       ! Loop over rows of src, getting values
       do ielem=col_par%first_elem_on_proc(proc,col_blks), &
            col_par%first_elem_on_proc(proc+1,col_blks)-1
          ielems = ielem - col_par%first_elem_on_proc(proc,col_blks) + 1
          do jelem=1,nrows
             call pdelget('O','O',buffer(jelem,ielems),src,jelem,ielem,desc)
          end do
       end do

       ! Sum this part of the matrix over procs
       call comms_reduce('SUM',buffer, &
            col_par%num_elems_on_proc(proc,col_blks)*library(destlib)%nrows)

       ! Copy buffer to local buffer if proc is local proc
       if (proc==pub_my_proc_id) then
          ielems = col_par%num_elems_on_proc(pub_my_proc_id,col_blks)
          loc_buffer(1:nrows,1:ielems) = buffer(1:nrows,1:ielems)
       end if

    end do  ! Loop over procs

    ! Loop over the segments of dest on this proc
    do seg=0,pub_total_num_procs-1

       ! Cycle if this segment is blank
       seg_type = library(destlib)%seg_info(s_type,seg)
       if (seg_type == SEG_BLANK) cycle
       if (seg_type == SEG_DENSE) nrows = row_par%num_elems_on_proc(seg,row_blks)

       ! Loop over block-columns of dest on this proc
       seg_start = library(destlib)%seg_info(s_idx,seg)
       loc_iblk = 0
       do iblk=col_par%my_first_blk,col_par%my_last_blk
          loc_iblk = loc_iblk + 1
          ielem = col_par%first_elem_on_atom(iblk,col_blks) - &
               col_par%first_elem_on_proc(pub_my_proc_id,col_blks) + 1
          ielems = col_par%num_elems_on_atom(iblk,col_blks)

          ! Loop over block-rows jblk in column iblk of dest
          do idx=library(destlib)%blk_idx(seg_start+loc_iblk-1), &
               library(destlib)%blk_idx(seg_start+loc_iblk)-1
             jblk = library(destlib)%blk_idx(idx)
             jelem = row_par%first_elem_on_atom(jblk,row_blks)
             jelems = row_par%num_elems_on_atom(jblk,row_blks)

             if (seg_type == SEG_SPARSE) nrows = jelems

             ! Copy block
             ptr = library(destlib)%blk_ptr(idx)
             do ielemonat=0,ielems-1      ! Element col in block
                do jelemonat=0,jelems-1   ! Element row in block
                   dest%dmtx(ptr+ielemonat*nrows+jelemonat) = &
                        loc_buffer(jelem+jelemonat,ielem+ielemonat)
                end do
             end do

          end do  ! Loop over block-rows in column iblk of src

       end do  ! Loop over block-columns of src

    end do  ! seg

    ! Dellocate buffers
    deallocate(buffer,stat=ierr)
    call utils_dealloc_check('sparse_blacstospam3_real','buffer',ierr)
    deallocate(loc_buffer,stat=ierr)
    call utils_dealloc_check('sparse_blacstospam3_real','loc_buffer',ierr)

  end subroutine sparse_blacstospam3_real

  !============================================================================!
  ! This subroutine converts a dense array in the BLACS distributed format     !
  ! suitable for use with ScaLAPACK into a SPAM3 block sparse matrix           !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   dest   (output) : The destination SPAM3 matrix (must be complex)         !
  !   src    (input)  : The source BLACS dense matrix (complex)                !
  !   ld     (input)  : Leading dimension of the array dest                    !
  !   nc     (input)  : Number of columns of the array dest                    !
  !   desc   (input)  : BLACS descriptor for the array dest                    !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, February 2010.                                   !
  !============================================================================!

  subroutine sparse_blacstospam3_complex(dest,src,ld,nc,desc)

    use comms, only: comms_reduce, pub_my_proc_id, pub_total_num_procs
    use parallel_strategy, only: PARAL_INFO
    use utils, only: utils_assert, utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    integer,intent(in) :: ld,nc
    type(SPAM3),intent(inout) :: dest
    complex(kind=DP),intent(in) :: src(ld,nc)
    integer,intent(in) :: desc(50)

    ! ScaLAPACK subroutine
    external :: pzelget

    ! Local Variables
    integer :: ierr
    integer :: destlib       ! Library entry for dest
    integer :: row_blks      ! Identifier for column blocking scheme
    integer :: col_blks      ! Identifier for row blocking scheme
    integer :: iblk,jblk     ! Block counters
    integer :: loc_iblk      ! Local counter for iblk
    integer :: idx           ! Index
    integer :: ptr           ! Pointer
    integer :: ielemonat     ! elem in atom block iblk counter
    integer :: jelemonat     ! elem in atom block jblk counter
    integer :: ielem,jelem   ! elem counters
    integer :: ielems,jelems ! Numbers of element rows/cols in this block
    integer :: nrows         ! Total number of rows in this block/segment
    integer :: proc          ! Proc loop counter
    integer :: seg           ! Segment index
    integer :: seg_type      ! Start position of this segment in the index
    integer :: seg_start     ! Start position of this segment in the index
    complex(kind=DP), allocatable :: buffer(:,:),loc_buffer(:,:)
    type(PARAL_INFO), pointer :: row_par ! Parallel strategy for rows.
    type(PARAL_INFO), pointer :: col_par ! Parallel strategy for columns.

    ! Get library entry for src
    destlib = dest%lib
    row_blks = library(destlib)%row_blks
    col_blks = library(destlib)%col_blks
    nrows = library(destlib)%nrows

    ! rc2013: specify row and column par pointers
    row_par => library(destlib)%col_par
    col_par => library(destlib)%col_par

    ! Check arguments
    call utils_assert(dest%iscmplx, &
         'Error in sparse_blacstospam3_real: complex matrices only.')

    ! Allocate buffers
    allocate(buffer(nrows,maxval(col_par%num_elems_on_proc(:,col_blks))),stat=ierr)
    call utils_alloc_check('sparse_blacstospam3_complex','buffer',ierr)
    allocate(loc_buffer(nrows,col_par%num_elems_on_proc(pub_my_proc_id,col_blks)), &
         stat=ierr)
    call utils_alloc_check('sparse_blacstospam3_complex','loc_buffer',ierr)

    ! Broadcast results across all procs
    nrows = library(destlib)%nrows
    do proc=0,pub_total_num_procs-1

       ! Zero buffer
       buffer(:,:) = cmplx(0.0_DP,0.0_DP,kind=DP)

       ! Loop over rows of src, getting values
       do ielem=col_par%first_elem_on_proc(proc,col_blks), &
            col_par%first_elem_on_proc(proc+1,col_blks)-1
          ielems = ielem - col_par%first_elem_on_proc(proc,col_blks) + 1
          do jelem=1,nrows
             call pzelget('O','O',buffer(jelem,ielems),src,jelem,ielem,desc)
          end do
       end do

       ! Sum this part of the matrix over procs
       call comms_reduce('SUM',buffer, &
            col_par%num_elems_on_proc(proc,col_blks)*library(destlib)%nrows)

       ! Copy buffer to local buffer if proc is local proc
       if (proc==pub_my_proc_id) then
          ielems = col_par%num_elems_on_proc(pub_my_proc_id,col_blks)
          loc_buffer(1:nrows,1:ielems) = buffer(1:nrows,1:ielems)
       end if

    end do  ! Loop over procs

    ! Loop over the segments of dest on this proc
    do seg=0,pub_total_num_procs-1

       ! Cycle if this segment is blank
       seg_type = library(destlib)%seg_info(s_type,seg)
       if (seg_type == SEG_BLANK) cycle
       if (seg_type == SEG_DENSE) nrows = row_par%num_elems_on_proc(seg,row_blks)

       ! Loop over block-columns of dest on this proc
       seg_start = library(destlib)%seg_info(s_idx,seg)
       loc_iblk = 0
       do iblk=col_par%my_first_blk,col_par%my_last_blk
          loc_iblk = loc_iblk + 1
          ielem = col_par%first_elem_on_atom(iblk,col_blks) - &
               col_par%first_elem_on_proc(pub_my_proc_id,col_blks) + 1
          ielems = col_par%num_elems_on_atom(iblk,col_blks)

          ! Loop over block-rows jblk in column iblk of dest
          do idx=library(destlib)%blk_idx(seg_start+loc_iblk-1), &
               library(destlib)%blk_idx(seg_start+loc_iblk)-1
             jblk = library(destlib)%blk_idx(idx)
             jelem = row_par%first_elem_on_atom(jblk,row_blks)
             jelems = row_par%num_elems_on_atom(jblk,row_blks)

             if (seg_type == SEG_SPARSE) nrows = jelems

             ! Copy block
             ptr = library(destlib)%blk_ptr(idx)
             do ielemonat=0,ielems-1      ! Element col in block
                do jelemonat=0,jelems-1   ! Element row in block
                   dest%zmtx(ptr+ielemonat*nrows+jelemonat) = &
                        loc_buffer(jelem+jelemonat,ielem+ielemonat)
                end do
             end do

          end do  ! Loop over block-rows in column iblk of src

       end do  ! Loop over block-columns of src

    end do  ! seg

    ! Dellocate buffers
    deallocate(buffer,stat=ierr)
    call utils_dealloc_check('sparse_blacstospam3_complex','buffer',ierr)
    deallocate(loc_buffer,stat=ierr)
    call utils_dealloc_check('sparse_blacstospam3_complex','loc_buffer',ierr)

  end subroutine sparse_blacstospam3_complex

#endif

!------------------------------------------------------------------------------

  !============================================================================!
  ! This subroutine is a wrapper to write a single sparse matrix to file by    !
  ! calling sparse_write_matrix with a single-element 2D array.                !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   mat      (input) : The sparse matrix to be written                       !
  !   filename (input) : The filename to use                                   !
  !   unit     (input) : The (optional) unit number to use if file is open     !
  !   append   (input) : (optional) flag to activate append (default: replace) !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, June 2009.                                       !
  ! Append option added by J.M. Escartin, Summer 2015.                         !
  ! Modified to allow writing 2D array of SPAM3 by F. Corsetti, November 2015. !
  !============================================================================!

  subroutine sparse_write_scalar(mat,filename,unit,append)

    implicit none

    ! Arguments
    type(SPAM3), intent(in) :: mat            ! Matrix to write
    character(len=*), intent(in) :: filename  ! Filename to use
    integer, optional, intent(in) :: unit     ! I/O unit
    logical, optional, intent(in) :: append   ! append (avoid replace)

    ! Local Variables
    type(SPAM3) :: mat_vec(1,1)

    call sparse_create(mat_vec(1,1),mat)
    call sparse_copy(mat_vec(1,1),mat)
    call sparse_write_matrix(mat_vec,filename,unit,append)
    call sparse_destroy(mat_vec(1,1))

  end subroutine sparse_write_scalar

!------------------------------------------------------------------------------

  !============================================================================!
  ! This subroutine is a wrapper to write a vector of sparse matrices to file  !
  ! by calling sparse_write_matrix with a (single x multi)-element 2D array.   !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   mat      (input) : The sparse matrix to be written                       !
  !   filename (input) : The filename to use                                   !
  !   unit     (input) : The (optional) unit number to use if file is open     !
  !   append   (input) : (optional) flag to activate append (default: replace) !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, June 2009.                                       !
  ! Append option added by J.M. Escartin, Summer 2015.                         !
  ! Modified to allow writing 2D array of SPAM3 by F. Corsetti, November 2015. !
  !============================================================================!

  subroutine sparse_write_vector(mat,filename,unit,append)

    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(SPAM3), intent(in) :: mat(:)         ! Matrix to write
    character(len=*), intent(in) :: filename  ! Filename to use
    integer, optional, intent(in) :: unit     ! I/O unit
    logical, optional, intent(in) :: append   ! append (avoid replace)

    ! Local Variables
    integer :: ierr
    integer :: i
    type(SPAM3), allocatable :: mat_vec(:,:)

    allocate(mat_vec(1,size(mat)),stat=ierr)
    call utils_alloc_check('sparse_write_vector','mat_vec',ierr)
    do i=1,size(mat)
       call sparse_create(mat_vec(1,i),mat(i))
       call sparse_copy(mat_vec(1,i),mat(i))
    end do
    call sparse_write_matrix(mat_vec,filename,unit,append)
    do i=1,size(mat)
       call sparse_destroy(mat_vec(1,i))
    end do
    deallocate(mat_vec,stat=ierr)
    call utils_dealloc_check('sparse_write_vector','mat_vec',ierr)

  end subroutine sparse_write_vector

!------------------------------------------------------------------------------

  !============================================================================!
  ! This subroutine writes a matrix of SPAM3 sparse matrice to a file.         !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   mat      (input) : The sparse matrix to be written                       !
  !   filename (input) : The filename to use                                   !
  !   unit     (input) : The (optional) unit number to use if file is open     !
  !   append   (input) : (optional) flag to activate append (default: replace) !
  !----------------------------------------------------------------------------!
  ! Written by Peter Haynes, January 2007.                                     !
  ! Modification for dense matrices by Nicholas Hine, Dec 2007.                !
  ! Significantly modified for SPAM3 by Nicholas Hine, June 2009.              !
  ! Append option added by J.M. Escartin, Summer 2015.                         !
  ! Modified to allow writing 2D array of SPAM3 by F. Corsetti, November 2015. !
  ! Modified for multiple parallel strategies by Robert Charlton, June 2018.   !
  !============================================================================!

  subroutine sparse_write_matrix(mat,filename,unit,append)

    use comms, only: comms_barrier, comms_bcast, comms_send, comms_free, &
         comms_recv, comms_reduce, pub_my_proc_id, pub_root_proc_id, &
         pub_total_num_procs, pub_on_root
#ifdef HDF5
    use hdf5
#endif
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_unit, &
         utils_assert

    implicit none

    ! Arguments
    type(SPAM3), intent(in) :: mat(:,:)       ! Matrix to write
    character(len=*), intent(in) :: filename  ! Filename to use
    integer, optional, intent(in) :: unit     ! I/O unit
    logical, optional, intent(in) :: append   ! append (avoid replace)

    ! Local variables
    integer :: ncomps(2)            ! Number of components for each dimension
    integer :: tot_ncomps           ! Total number of components
    integer :: icomp1, icomp2       ! Loop counters for components
    integer :: row_blks             ! Identifier for row blocking scheme
    integer :: col_blks             ! Identifier for col blocking scheme
    integer :: ierr                 ! Error flag
    integer :: ilib                 ! Library entry for mat
    integer :: iunit                ! Unit number for writing file
    integer :: iblk                 ! Atom loop counter
    integer :: iat_orig             ! Atom iblk in original order
    integer :: jblk                 ! Second atom
    integer :: jat_orig             ! Atom jblk in original order
    integer :: loc_iblk             ! Local atom counter for iat
    integer :: iidx                 ! Index loop counter
    integer :: jidx                 ! Index info counter
    integer :: idxlen               ! Index info length
    integer :: proc                 ! Processor loop counter
    integer :: datlen               ! Buffer length for data
    integer :: blk_start,blk_end    ! First and last blocks on proc
    integer :: ielems               ! Number of elements on atom iat
    integer, allocatable :: block_sizes(:)  ! Block sizes
    integer, allocatable :: idxbuf(:)       ! Comm buffer for index
    integer, allocatable :: idxinfo(:)      ! Index info for writing
    character(len=11) :: tag        ! Field tag
    real(kind=DP), allocatable :: dmtxbuf(:)    ! Comm buffer for real data
    complex(kind=DP), allocatable :: zmtxbuf(:) ! Comm buffer for complex data
    logical :: square               ! Is the matrix square?
    logical :: loc_append           ! append (avoid replace)
    type(PARAL_INFO), pointer :: row_par ! Parallel strategy for rows.
    type(PARAL_INFO), pointer :: col_par ! Parallel strategy for columns.
    ! rc2013: extra block sizes for row + col structures
    integer, allocatable :: nblock_sizes(:)  ! Block sizes
    integer, allocatable :: mblock_sizes(:)  ! Block sizes
#ifdef HDF5
    character(50) :: dsetname, str1, str2

    integer(hid_t) :: filesystem_id
    integer(hid_t) :: file_id
    integer(hid_t) :: group_id
    integer(hid_t) :: attr_id
    integer(hid_t) :: aspace_id1, aspace_id2, aspace_idb
    integer(hid_t) :: dset_id
    integer(hid_t) :: dspace_id
    integer(hsize_t) :: dims(1), adims(1), adimsb(2)

    integer :: dkn_dims(2)
    ! rc2013: in general block_sizesb will need different array lengths
    ! rc2013: for rows and columns
    integer, allocatable :: block_sizesb(:,:)
#endif

    ! Obtain number of components
    ncomps = shape(mat)
    tot_ncomps=ncomps(1)*ncomps(2)

    ! Obtain library entry for mat
    ilib = mat(1,1)%lib
    row_blks = library(ilib)%row_blks
    col_blks = library(ilib)%col_blks

    ! rc2013: check that the parallel strategies match i.e. the matrix is square
    ! Treat this like a trace (rows/cols must match)
!    call sparse_par_check(ilib, ilib, routine='sparse_write_matrix', trace=.true.)

    ! rc2013: specify row and column par pointers
    row_par => library(ilib)%row_par
    col_par => library(ilib)%col_par

    ! Check that the matrix is square
    square = .true.
    if (row_blks/=col_blks .or. row_par%par_index/=col_par%par_index) &
         square = .false.

    ! Determine whether the open call will be append or not (=replace).
    if (present(append)) then
       loc_append = append
    else
       loc_append = .false.
    end if

    ! Preliminaries on root proc only
    if (pub_on_root) then

#ifdef HDF5
       call utils_assert(.not.present(append), &
            'ERROR: append functionality not available yet with HDF5.')

       iunit = -1

       call h5open_f(ierr)
       call h5fcreate_f(trim(filename),h5f_acc_trunc_f,filesystem_id,ierr)
       call h5gcreate_f(filesystem_id,'sparse',file_id,ierr)

       adims=(/1/)
       call h5screate_simple_f(1,adims,aspace_id1,ierr)
       call h5acreate_f(file_id,'version',h5t_native_integer, &
            aspace_id1,attr_id,ierr)
       call h5awrite_f(attr_id,h5t_native_integer,nint(file_version*100), &
            adims,ierr)
       call h5aclose_f(attr_id,ierr)
       call h5acreate_f(file_id,'num_atoms',h5t_native_integer, &
            aspace_id1,attr_id,ierr)
       call h5awrite_f(attr_id,h5t_native_integer,library(ilib)%nblk,adims,ierr)
       call h5aclose_f(attr_id,ierr)
       call h5acreate_f(file_id,'type',h5t_native_integer, &
            aspace_id1,attr_id,ierr)
       if (mat(1,1)%iscmplx) then
          call h5awrite_f(attr_id,h5t_native_integer,-1,adims,ierr)
       else
          call h5awrite_f(attr_id,h5t_native_integer,0,adims,ierr)
       end if
       call h5aclose_f(attr_id,ierr)
       call h5acreate_f(file_id,'num_comps',h5t_native_integer, &
            aspace_id1,attr_id,ierr)
       call h5awrite_f(attr_id,h5t_native_integer,tot_ncomps,adims,ierr)
       call h5aclose_f(attr_id,ierr)
       call h5acreate_f(file_id,'num_chunks',h5t_native_integer, &
            aspace_id1,attr_id,ierr)
       call h5awrite_f(attr_id,h5t_native_integer,pub_total_num_procs,adims, &
            ierr)
       call h5aclose_f(attr_id,ierr)
       if (square) then
          call h5acreate_f(file_id,'num_ngwfs',h5t_native_integer, &
               aspace_id1,attr_id,ierr)
          call h5awrite_f(attr_id,h5t_native_integer,library(ilib)%nrows, &
               adims,ierr)
          call h5aclose_f(attr_id,ierr)
          adims=(/library(ilib)%nblk/)
          call h5screate_simple_f(1,adims,aspace_idb,ierr)
          call h5acreate_f(file_id,'block_sizes',h5t_native_integer, &
               aspace_idb,attr_id,ierr)
          allocate(block_sizes(library(ilib)%nblk),stat=ierr)
          call utils_alloc_check('sparse_write_matrix','block_sizes',ierr)
          do iblk=1,library(ilib)%nblk
             iat_orig = row_par%orig_atom(iblk)
             block_sizes(iat_orig) = row_par%num_elems_on_atom(iblk,row_blks)
          end do
          call h5awrite_f(attr_id,h5t_native_integer,block_sizes,adims,ierr)
          deallocate(block_sizes,stat=ierr)
          call utils_dealloc_check('sparse_write_matrix','block_sizes',ierr)
          call h5aclose_f(attr_id,ierr)
          adims=(/1/)
       else
          adims=(/2/)
          call h5screate_simple_f(1,adims,aspace_id2,ierr)
          call h5acreate_f(file_id,'dkn_dims',h5t_native_integer, &
               aspace_id2,attr_id,ierr)
          dkn_dims=(/library(ilib)%nrows,library(ilib)%mcols/)
          call h5awrite_f(attr_id,h5t_native_integer,dkn_dims,adims,ierr)
          call h5aclose_f(attr_id,ierr)
          adims=(/1/)
          adimsb=(/2,library(ilib)%nblk/)
          call h5screate_simple_f(2,adimsb,aspace_idb,ierr)
          call h5acreate_f(file_id,'block_dims',h5t_native_integer, &
               aspace_idb,attr_id,ierr)
          allocate(block_sizesb(2,library(ilib)%nblk),stat=ierr)
          call utils_alloc_check('sparse_write_matrix','block_sizesb',ierr)
          ! rc2013: split up this do loop for rows and columns separately
          do iblk=1,library(ilib)%nblk
             iat_orig = row_par%orig_atom(iblk)
             block_sizesb(1,iat_orig) = row_par%num_elems_on_atom(iblk,row_blks)
          end do
          ! This will only work if row_par=col_par i.e. nblk=mblk
          do iblk=1,library(ilib)%mblk
             iat_orig = col_par%orig_atom(iblk)
             block_sizesb(2,iat_orig) = col_par%num_elems_on_atom(iblk,col_blks)
          end do
          call h5awrite_f(attr_id,h5t_native_integer,block_sizesb,adimsb,ierr)
          deallocate(block_sizesb,stat=ierr)
          call utils_dealloc_check('sparse_write_matrix','block_sizesb',ierr)
          call h5aclose_f(attr_id,ierr)
       end if
#else
       if (present(unit)) then
          iunit = unit
          call utils_assert(.not.present(append), &
               'ERROR: Optional append provided in call to &
               &sparse_write_matrix, but unused - unit also specified.')
       else
          iunit = utils_unit()

          if (loc_append) then
             open(unit=iunit,file=trim(filename),form='unformatted', &
                  action='write',iostat=ierr,position='append',status='old')
          else
             open(unit=iunit,file=trim(filename),form='unformatted', &
                  action='write',iostat=ierr,status='replace')
          end if

          call utils_assert(ierr == 0, 'Error in sparse_write_matrix: &
               &opening file "'//trim(filename)//'" failed with code ',ierr)
       end if

       ! Version - version*100 so 150 means version 1.5
       tag = 'VERSION'
       write(iunit) tag//'I'
       write(iunit) 1
       write(iunit) nint(file_version*100)

       ! Number of atoms/blocks
       if(row_par%par_index == col_par%par_index) then
          tag = 'NUM_ATOMS'
          write(iunit) tag//'I'
          write(iunit) 1
          write(iunit) library(ilib)%nblk
       else
          ! rc2013: print nblk and mblk if they're not the same
          tag = 'NUM_ATOMS_R'
          write(iunit) tag//'I'
          write(iunit) 1
          write(iunit) library(ilib)%nblk

          tag = 'NUM_ATOMS_C'
          write(iunit) tag//'I'
          write(iunit) 1
          write(iunit) library(ilib)%mblk
       end if

       ! Number of NGWFs
       if (square) then
          tag = 'NUM_NGWFS'
          write(iunit) tag//'I'
          write(iunit) 1
          write(iunit) library(ilib)%nrows
       else
          tag = 'NUM_ROWS'
          write(iunit) tag//'I'
          write(iunit) 1
          write(iunit) library(ilib)%nrows

          tag = 'NUM_COLS'
          write(iunit) tag//'I'
          write(iunit) 1
          write(iunit) library(ilib)%mcols
       end if

       ! Data type: real or complex
       tag = 'TYPE'
       write(iunit) tag//'I'
       write(iunit) 1
       if (mat(1,1)%iscmplx) then
          write(iunit) -1
       else
          write(iunit) 0
       end if

       ! Number of components
       tag = 'NUM_COMPS'
       write(iunit) tag//'I'
       write(iunit) 1
       write(iunit) tot_ncomps

       ! List of block sizes (in original order)
       ! rc2013: nblk =/= mblk if there are different par's involved
       allocate(nblock_sizes(library(ilib)%nblk),stat=ierr)
       call utils_alloc_check('sparse_write_matrix','nblock_sizes',ierr)
       allocate(mblock_sizes(library(ilib)%mblk),stat=ierr)
       call utils_alloc_check('sparse_write_matrix','mblock_sizes',ierr)
       if (square) then
          tag = 'BLOCK_SIZES'
          write(iunit) tag//'I'
          write(iunit) library(ilib)%nblk
          do iblk=1,library(ilib)%nblk
             iat_orig = row_par%orig_atom(iblk)
             nblock_sizes(iat_orig) = row_par%num_elems_on_atom(iblk,row_blks)
          end do
          write(iunit) nblock_sizes
       else
          tag = 'BLOCK_ROWS'
          write(iunit) tag//'I'
          write(iunit) library(ilib)%nblk
          do iblk=1,library(ilib)%nblk
             iat_orig = row_par%orig_atom(iblk)
             nblock_sizes(iat_orig) = row_par%num_elems_on_atom(iblk,row_blks)
          end do
          write(iunit) nblock_sizes

          tag = 'BLOCK_COLS'
          write(iunit) tag//'I'
          write(iunit) library(ilib)%mblk
          do iblk=1,library(ilib)%mblk
             iat_orig = col_par%orig_atom(iblk)
             mblock_sizes(iat_orig) = col_par%num_elems_on_atom(iblk,col_blks)
          end do
          write(iunit) mblock_sizes
       end if
       deallocate(nblock_sizes,stat=ierr)
       call utils_dealloc_check('sparse_write_matrix','nblock_sizes',ierr)
       deallocate(mblock_sizes,stat=ierr)
       call utils_dealloc_check('sparse_write_matrix','mblock_sizes',ierr)
#endif


       ! End of root proc preliminaries
    else
       iunit = -1 ! qoh: Initialise to avoid compiler warning
    end if

    ! Allocate buffers to receive data communicated from procs
    if (mat(1,1)%iscmplx) then
       allocate(zmtxbuf(library(ilib)%max_nze*tot_ncomps),stat=ierr)
       call utils_alloc_check('sparse_write_matrix','zmtxbuf',ierr)
    else
       allocate(dmtxbuf(library(ilib)%max_nze*tot_ncomps),stat=ierr)
       call utils_alloc_check('sparse_write_matrix','dmtxbuf',ierr)
    end if

    ! Global preliminaries
    idxlen = sparse_index_length(mat(1,1)) + 1
    call comms_reduce('MAX',idxlen)
    allocate(idxinfo(idxlen*3),stat=ierr)
    call utils_alloc_check('sparse_write_matrix','idxinfo',ierr)
    allocate(idxbuf(idxlen),stat=ierr)
    call utils_alloc_check('sparse_write_matrix','idxbuf',ierr)

    ! Get the index for this proc in non-segmented format
    call sparse_generate_index(idxbuf,mat(1,1))

    ! Read segmented data into unsegmented buffers
    call internal_unsegment

    call comms_barrier

    ! Loop over processors and communicate with root proc
    do proc=0,pub_total_num_procs-1
       blk_start = col_par%first_atom_on_proc(proc)
       blk_end = col_par%first_atom_on_proc(proc+1) - 1
       ! rc2013: skip this proc if there's nothing here
       if(blk_start == 0) cycle

       ! Send index of mat stored on proc to root proc
       if ((pub_my_proc_id == proc).and.(.not.pub_on_root)) then
          idxbuf(idxlen) = datlen
          call comms_send(pub_root_proc_id,idxbuf,idxlen)
       end if
       if (proc /= pub_root_proc_id .and. pub_on_root) then
          call comms_recv(proc,idxbuf,idxlen)
       end if

       ! Send data of mat stored on proc to root proc
       if (pub_my_proc_id == proc) then

          if (.not.pub_on_root) then
             if (mat(1,1)%iscmplx) then
                call comms_send(pub_root_proc_id,zmtxbuf,datlen*tot_ncomps)
             else
                call comms_send(pub_root_proc_id,dmtxbuf,datlen*tot_ncomps)
             end if
          end if

       end if

       ! Receive data of mat on root proc from other proc
       if (pub_on_root .and. proc /= pub_root_proc_id) then
          datlen = idxbuf(idxlen)
          if (mat(1,1)%iscmplx) then
             call comms_recv(proc,zmtxbuf,datlen*tot_ncomps)
          else
             call comms_recv(proc,dmtxbuf,datlen*tot_ncomps)
          end if
       end if

       if (pub_on_root) then
          ! Now write out indexing information for this data
          jidx = 0
          loc_iblk = 0
          do iblk=blk_start,blk_end
          !do iblk=col_par%my_first_blk,col_par%my_last_blk
             loc_iblk = loc_iblk + 1
             iat_orig = col_par%orig_atom(iblk)
             ielems = col_par%num_elems_on_atom(iblk,col_blks)

             ! Loop over block-rows in column iat
             do iidx=idxbuf(loc_iblk),idxbuf(loc_iblk+1)-1
                jblk = idxbuf(iidx)
                jat_orig = row_par%orig_atom(jblk)

                ! Generate index info: row, column and number of elements
                idxinfo(jidx+1) = jat_orig
                idxinfo(jidx+2) = iat_orig
                idxinfo(jidx+3) = ielems * row_par%num_elems_on_atom(jblk,row_blks)
                jidx = jidx + 3
             end do
          end do
#ifdef HDF5
          write(str1,*) proc+1
          dsetname="chunk_"//trim(adjustl(str1))
          call h5gcreate_f(file_id,adjustl(dsetname),group_id,ierr)
          call h5acreate_f(group_id,'index_dim',h5t_native_integer, &
               aspace_id1,attr_id,ierr)
          call h5awrite_f(attr_id,h5t_native_integer,jidx,adims,ierr)
          call h5aclose_f(attr_id,ierr)
          call h5acreate_f(group_id,'data_dim',h5t_native_integer, &
               aspace_id1,attr_id,ierr)
          call h5awrite_f(attr_id,h5t_native_integer,datlen*tot_ncomps,adims,ierr)
          call h5aclose_f(attr_id,ierr)
          dims=(/jidx/)
          call h5screate_simple_f(1,dims,dspace_id,ierr)
          call h5dcreate_f(group_id,'index',h5t_native_integer,dspace_id, &
               dset_id,ierr)
          call h5dwrite_f(dset_id,h5t_native_integer,idxinfo(1:jidx),dims,ierr)
          call h5dclose_f(dset_id,ierr)
          call h5sclose_f(dspace_id,ierr)
          dims=(/datlen*tot_ncomps/)
          call h5screate_simple_f(1,dims,dspace_id,ierr)
          if (mat(1,1)%iscmplx) then
             call h5dcreate_f(group_id,'data_real',h5t_native_double, &
                  dspace_id,dset_id,ierr)
             call h5dwrite_f(dset_id,h5t_native_double, &
                  real(zmtxbuf(1:datlen*tot_ncomps),DP),dims,ierr)
             call h5dclose_f(dset_id,ierr)
             call h5dcreate_f(group_id,'data_imag',h5t_native_double, &
                  dspace_id,dset_id,ierr)
             call h5dwrite_f(dset_id,h5t_native_double, &
                  aimag(zmtxbuf(1:datlen*tot_ncomps)),dims,ierr)
             call h5dclose_f(dset_id,ierr)
          else
             call h5dcreate_f(group_id,'data',h5t_native_double,dspace_id, &
                  dset_id,ierr)
             call h5dwrite_f(dset_id,h5t_native_double, &
                  dmtxbuf(1:datlen*tot_ncomps),dims,ierr)
             call h5dclose_f(dset_id,ierr)
          endif
          call h5sclose_f(dspace_id,ierr)
          call h5gclose_f(group_id,ierr)
#else
          tag = 'INDEX'
          write(iunit) tag//'I'
          write(iunit) jidx
          write(iunit) idxinfo(1:jidx)

          ! Now the data itself
          tag = 'DATA'
          if (mat(1,1)%iscmplx) then
             write(iunit) tag//'Z'
             write(iunit) datlen*tot_ncomps
             write(iunit) zmtxbuf(1:datlen*tot_ncomps)
          else
             write(iunit) tag//'D'
             write(iunit) datlen*tot_ncomps
             write(iunit) dmtxbuf(1:datlen*tot_ncomps)
          end if
#endif

       end if

       ! Loop over procs
    end do

    if (pub_on_root) then
#ifdef HDF5
       call h5gclose_f(file_id,ierr)
       call h5fclose_f(filesystem_id,ierr)
       call h5close_f(ierr)
#else
       ! Write a terminating tag
       tag = 'ENDFILE'
       write(iunit) tag//'I'
       write(iunit) 1
       write(iunit) 0

       ! Close file for output
       if (.not. present(unit)) then
          close(unit=iunit,iostat=ierr)
          call utils_assert(ierr == 0, 'Error in sparse_write_matrix: &
               &closing file "'//trim(filename)//'" failed with code ',ierr)
       end if
#endif
    end if

    ! Re-sync procs
    call comms_free
    call comms_barrier

    ! Deallocate workspace
    deallocate(idxbuf,stat=ierr)
    call utils_dealloc_check('sparse_write_matrix','idxbuf',ierr)
    deallocate(idxinfo,stat=ierr)
    call utils_dealloc_check('sparse_write_matrix','idxinfo',ierr)

    if (mat(1,1)%iscmplx) then
       deallocate(zmtxbuf,stat=ierr)
       call utils_dealloc_check('sparse_write_matrix','zmtxbuf',ierr)
    else
       deallocate(dmtxbuf,stat=ierr)
       call utils_dealloc_check('sparse_write_matrix','dmtxbuf',ierr)
    end if

  contains

    subroutine internal_unsegment

      implicit none

      ! Locals
      integer :: iidx,ptr
      integer :: iblk,loc_iblk
      integer :: jblk
      integer :: jelems
      integer :: ielem
      integer :: max_rows_on_atom
      integer :: max_cols_on_atom
      real (kind=DP), allocatable :: dblk(:,:)
      complex (kind=DP), allocatable :: zblk(:,:)

      ! Allocate temporary array to store block
      max_rows_on_atom = maxval(row_par%num_elems_on_atom(:,row_blks))
      max_cols_on_atom = maxval(col_par%num_elems_on_atom(:,col_blks))
      if (mat(1,1)%iscmplx) then
         allocate(zblk(max_rows_on_atom,max_cols_on_atom),stat=ierr)
         call utils_alloc_check('internal_unsegment (sparse_write_matrix)', &
              'zblk',ierr)
      else
         allocate(dblk(max_rows_on_atom,max_cols_on_atom),stat=ierr)
         call utils_alloc_check('internal_unsegment (sparse_write_matrix)', &
              'dblk',ierr)
      end if

      ptr = 1

      ! Loop over components of the matrix array
      do icomp1=1,ncomps(1)
         do icomp2=1,ncomps(2)

            loc_iblk = 0
            ! Loop over block-cols on this proc
            do iblk=col_par%my_first_blk,col_par%my_last_blk
               loc_iblk = loc_iblk + 1

               ! Loop over nonzero block-rows in this block-column
               do iidx=idxbuf(loc_iblk),idxbuf(loc_iblk+1)-1
                  jblk = idxbuf(iidx)
                  jelems = row_par%num_elems_on_atom(jblk,row_blks)

                  ! Get this block and write it into buffer
                  if (mat(1,1)%iscmplx) then
                     call sparse_get_block(zblk,mat(icomp1,icomp2),jblk,iblk)
                     do ielem=1,col_par%num_elems_on_atom(iblk,col_blks)
                        zmtxbuf(ptr:ptr+jelems-1) = zblk(1:jelems,ielem)
                        ptr = ptr + jelems
                     end do
                  else
                     call sparse_get_block(dblk,mat(icomp1,icomp2),jblk,iblk)
                     do ielem=1,col_par%num_elems_on_atom(iblk,col_blks)
                        dmtxbuf(ptr:ptr+jelems-1) = dblk(1:jelems,ielem)
                        ptr = ptr + jelems
                     end do
                  end if

               end do

            end do

         end do
      end do

      ! Record total length of data in each matrix
      datlen = (ptr - 1) / tot_ncomps

      ! Deallocate temporary array for block
      if (mat(1,1)%iscmplx) then
         deallocate(zblk,stat=ierr)
         call utils_dealloc_check('internal_unsegment (sparse_write_matrix)', &
              'zblk',ierr)
      else
         deallocate(dblk,stat=ierr)
         call utils_dealloc_check('internal_unsegment (sparse_write_matrix)', &
              'dblk',ierr)
      end if

    end subroutine internal_unsegment

  end subroutine sparse_write_matrix

!------------------------------------------------------------------------------

  !============================================================================!
  ! This subroutine is a wrapper to read a single sparse matrix from a file by !
  ! calling sparse_read_matrix with a single-element 2D array.                 !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   mat      (inout) : The sparse matrix to be read                          !
  !   filename (input) : The filename to use                                   !
  !   unit     (inout) : The (optional) unit number to use if file is open     !
  !   get_unit (input) : (optional) unit is unknown at entry, valid at output  !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, June 2009                                        !
  ! get_unit option added by J.M. Escartin, Summer 2015.                       !
  ! Modified to allow reading 2D array of SPAM3 by F. Corsetti, November 2015. !
  !============================================================================!

  subroutine sparse_read_scalar(mat, filename, unit, get_unit)

    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(SPAM3), intent(inout) :: mat          ! Matrix to write
    character(len=*), intent(in) :: filename   ! Filename to use
    integer, optional, intent(inout) :: unit   ! I/O unit
    logical, optional, intent(in) :: get_unit  ! get unit number from this call

    ! Locals
    integer :: ierr
    type(SPAM3),allocatable :: mat_copy(:,:)

    ! Allocate and create a temporary SPAM3 2D array of length 1x1 to read into
    allocate(mat_copy(1,1),stat=ierr)
    call utils_alloc_check('sparse_read_scalar','mat_copy',ierr)
    call sparse_create(mat_copy(1,1),mat)

    ! Call matrix version of sparse_read
    call sparse_read_matrix(mat_copy, filename, unit, get_unit)

    ! Copy temporary array version into argument
    call sparse_copy(mat,mat_copy(1,1))

    ! Destroy and deallocate temporary SPAM3 array
    call sparse_destroy(mat_copy(1,1))
    deallocate(mat_copy,stat=ierr)
    call utils_dealloc_check('sparse_read_scalar','mat_copy',ierr)

  end subroutine sparse_read_scalar

!------------------------------------------------------------------------------

  !============================================================================!
  ! This subroutine is a wrapper to read a vector of sparse matrices from a    !
  ! file by calling sparse_read_matrix with a (single x multi)-element 2D      !
  ! array                                                                      !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   mat      (inout) : The sparse matrix to be read                          !
  !   filename (input) : The filename to use                                   !
  !   unit     (inout) : The (optional) unit number to use if file is open     !
  !   get_unit (input) : (optional) unit is unknown at entry, valid at output  !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, June 2009                                        !
  ! get_unit option added by J.M. Escartin, Summer 2015.                       !
  ! Modified to allow reading 2D array of SPAM3 by F. Corsetti, November 2015. !
  !============================================================================!

  subroutine sparse_read_vector(mat, filename, unit, get_unit)

    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(SPAM3), intent(inout) :: mat(:)       ! Matrix to write
    character(len=*), intent(in) :: filename   ! Filename to use
    integer, optional, intent(inout) :: unit   ! I/O unit
    logical, optional, intent(in) :: get_unit  ! get unit number from this call

    ! Locals
    integer :: ierr
    integer :: i
    type(SPAM3),allocatable :: mat_copy(:,:)

    ! Allocate and create a temporary SPAM3 2D array of length 1xsize(mat) to read into
    allocate(mat_copy(1,size(mat)),stat=ierr)
    call utils_alloc_check('sparse_read_vector','mat_copy',ierr)
    do i=1,size(mat)
       call sparse_create(mat_copy(1,i),mat(i))
    end do

    ! Call matrix version of sparse_read
    call sparse_read_matrix(mat_copy, filename, unit, get_unit)

    ! Copy temporary array version into argument
    do i=1,size(mat)
       call sparse_copy(mat(i),mat_copy(1,i))
    end do

    ! Destroy and deallocate temporary SPAM3 array
    do i=1,size(mat)
       call sparse_destroy(mat_copy(1,i))
    end do
    deallocate(mat_copy,stat=ierr)
    call utils_dealloc_check('sparse_read_vector','mat_copy',ierr)

  end subroutine sparse_read_vector

!------------------------------------------------------------------------------

  !============================================================================!
  ! This subroutine reads a matrix of sparse matrices from a file.             !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   mat     (output) : The sparse matrix to be read                          !
  !   filename (input) : The filename                                          !
  !   unit     (inout) : The (optional) unit number to use if file is open     !
  !   get_unit (input) : (optional) unit is unknown at entry, valid at output  !
  !----------------------------------------------------------------------------!
  ! Written by Peter Haynes, February 2007.                                    !
  ! Modification for SPAM3 by Nicholas Hine, June 2009                         !
  ! get_unit option added by J.M. Escartin, Summer 2015.                       !
  ! Modified to allow reading 2D array of SPAM3 by F. Corsetti, November 2015. !
  ! Modified for multiple parallel strategies by Robert Charlton, June 2018.   !
  !============================================================================!

  subroutine sparse_read_matrix(mat, filename, unit, get_unit)

    use comms, only: comms_barrier, comms_bcast, comms_send, &
         comms_recv, comms_free, pub_my_proc_id, pub_on_root, &
         pub_root_proc_id, pub_total_num_procs
#ifdef HDF5
    use hdf5
#endif
    use rundat, only: pub_devel_code, pub_use_swx
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_unit, &
         utils_assert, utils_abort, utils_devel_code, utils_feature_not_supported

    implicit none

    ! Arguments
    type(SPAM3), intent(inout) :: mat(:,:)     ! Matrix to read
    character(len=*), intent(in) :: filename   ! Filename to use
    integer, optional, intent(inout) :: unit   ! I/O unit
    logical, optional, intent(in) :: get_unit  ! get unit number from this call

    ! Local variables
    integer :: ncomps(2)            ! Number of components for each dimension
    integer :: tot_ncomps           ! Total number of components
    integer :: icomp1, icomp2       ! Loop counters for components
    integer :: ierr                 ! Error flag
    integer :: ilib                 ! Library entry for mat
    integer :: row_blks             ! Identifier for blocking scheme
    integer :: col_blks             ! Identifier for blocking scheme
    integer :: iunit                ! Unit number for writing file
    integer :: idum                 ! Dummy variable
    integer :: file_nat             ! Number of atoms from file
    integer :: file_num             ! Number of NGWFs from file
    integer :: iblk                 ! Block loop counter
    integer :: iat_orig             ! Atom iblk in original order
    integer :: jblk                 ! Second block
    integer :: jat_orig             ! Atom jblk in original order
    integer :: iidx                 ! Index loop counter
    integer :: ielem,jelems         ! Element loop counter
    integer :: iptr                 ! Pointer
    integer :: idxlen               ! Index info length
    integer :: proc                 ! Processor loop counter
    integer :: idxptr               ! Index pointer
    integer :: datptr               ! Data pointer
    integer :: idxptr0              ! Index pointer start for a proc
    integer :: datptr0              ! Data pointer start for a proc
    integer :: comp0                ! Offset between components in data
    integer :: blksize              ! Block size
    integer :: datlen               ! Buffer length for data
    integer :: size_idxinfo         ! Workaround for apparent bug in PGI F90
    logical :: endfile_found        ! Workaround for bug in ifort 2013
    logical :: square               ! Is the matrix square?
    logical :: read_bigendian       ! Whether to try to read big-endian data
    integer, allocatable :: block_sizes(:)  ! Block sizes
    integer, allocatable :: ibuf(:)         ! Integer buffer
    real(kind=DP), allocatable :: dbuf(:)   ! Double precision real buffer
    complex(kind=DP), allocatable :: zbuf(:)! Complex buffer
    integer, allocatable :: idxbuf(:)       ! Comm buffer for index
    integer, allocatable :: idxinfo(:)      ! Index info for reading
    character(len=12) :: tag        ! Field tag
    character(len=12) :: cbuf       ! Character buffer
    real(kind=DP), allocatable :: dmtxbuf(:)    ! Comm buffer for real data
    complex(kind=DP), allocatable :: zmtxbuf(:) ! Comm buffer for complex data
    integer :: max_rows_on_atom                ! Maximum number of elements
    integer :: max_cols_on_atom                ! Maximum number of elements
    real (kind=DP), allocatable :: dblk(:,:)    ! Buffer for real block
    complex (kind=DP), allocatable :: zblk(:,:) ! Buffer for complex block
    type(PARAL_INFO), pointer :: row_par ! Parallel strategy for rows.
    type(PARAL_INFO), pointer :: col_par ! Parallel strategy for columns.
    integer :: file_natr             ! Number of row-atoms from file
    integer :: file_natc             ! Number of column-atoms from file
#ifdef HDF5
    character(50) :: dsetname, str1, str2

    integer(hid_t) :: filesystem_id
    integer(hid_t) :: file_id
    integer(hid_t) :: group_id
    integer(hid_t) :: attr_id
    integer(hid_t) :: dset_id
    integer(hsize_t) :: dims(1), adims(1)

    integer :: iel
    integer :: ichunk, nchunks
#endif

    logical :: loc_get_unit

    ! Obtain number of components
    ncomps = shape(mat)
    tot_ncomps=ncomps(1)*ncomps(2)

    ! Obtain library entry for mat
    ilib = mat(1,1)%lib
    row_blks = library(ilib)%row_blks
    col_blks = library(ilib)%col_blks

    ! rc2013: specify row and column par pointers
    row_par => library(ilib)%row_par
    col_par => library(ilib)%col_par

    ! Check if the matrix is square
    square = .true.
    if (row_blks/=col_blks .or. row_par%par_index/=col_par%par_index) &
         square = .false.

    ! Determine whether we want to get the unit number from this call.
    if (present(get_unit)) then
       loc_get_unit = get_unit
    else
       loc_get_unit = .false.
    end if

    ! Check if we want to read big-endian formatted data
    read_bigendian = utils_devel_code(.false., 'RESTART','BIGENDIAN', &
       pub_devel_code)

    ! Preliminaries on root proc only
    if (pub_on_root) then

#ifdef HDF5
       !call utils_assert(.not.present(get_unit), &
       !     'ERROR: append functionality not available yet with HDF5.')

       iunit = -1

       call h5open_f(ierr)
       call h5fopen_f(trim(filename),h5f_acc_rdonly_f,filesystem_id,ierr)
       call utils_assert(ierr == 0, 'Error in sparse_read_matrix: &
            &opening file "'//trim(filename)//'" failed with code ',ierr)
       call h5gopen_f(filesystem_id,'sparse',file_id,ierr)
       call utils_assert(ierr == 0, 'Error in sparse_read_matrix: &
            &opening root group in file "'//trim(filename)//'" failed with code ',ierr)

       adims=(/1/)
       call h5aopen_f(file_id,'version',attr_id,ierr)
       call h5aread_f(attr_id,h5t_native_integer,idum,adims,ierr)
       if (idum > nint(file_version*100)) then
          write(stdout,'(a,f5.2,3a,f5.2,a)') &
               'WARNING in sparse_read_matrix: version', &
               real(idum,kind=DP)*0.01_DP,' of file "',trim(filename), &
               '" is more recent than implemented version ', &
               file_version,': continuing'
       end if
       call h5aclose_f(attr_id,ierr)
       call h5aopen_f(file_id,'num_atoms',attr_id,ierr)
       call h5aread_f(attr_id,h5t_native_integer,file_nat,adims,ierr)
       call utils_assert(file_nat == library(ilib)%nblk, &
            'Error in sparse_read_matrix: atom number mismatch - file &
            &"'//trim(filename)//'" specifies a different number of atoms &
            &than expected by the library', file_nat, library(ilib)%nblk)
       call h5aclose_f(attr_id,ierr)
       call h5aopen_f(file_id,'type',attr_id,ierr)
       call h5aread_f(attr_id,h5t_native_integer,idum,adims,ierr)
       call utils_assert(idum /= -1 .or. mat(1,1)%iscmplx, &
            'Error in sparse_read_matrix: data type mismatch - file "'//&
            trim(filename)//'" specifies complex data but expecting real &
            &data.')
       call utils_assert(idum /= 0 .or. .not. mat(1,1)%iscmplx, &
            'Error in sparse_read_matrix: data type mismatch - file "'//&
            trim(filename)//'" specifies real data but expecting complex &
            &data.')
       call h5aclose_f(attr_id,ierr)
       call h5aopen_f(file_id,'num_comps',attr_id,ierr)
       call h5aread_f(attr_id,h5t_native_integer,idum,adims,ierr)
       call utils_assert(idum == tot_ncomps, 'Error in sparse_read_vector: &
            &component number mismatch - file "'//trim(filename)//&
            '" specifies idum components but expecting ncomps. &
            &Values for idum and ncomps follow. ', idum, tot_ncomps)
       call h5aclose_f(attr_id,ierr)
       call h5aopen_f(file_id,'num_ngwfs',attr_id,ierr)
       call h5aread_f(attr_id,h5t_native_integer,file_num,adims,ierr)
       if(pub_use_swx) then
          call utils_assert(file_num == library(ilib)%nrows, &
               'Error in sparse_read_matrix: NGWF (or other basis function) &
               &number mismatch - file &
               &"'//trim(filename)//'" specifies this many NGWFs &
               &(or other basis functions)', file_num)
       else
          call utils_assert(file_num == library(ilib)%nrows, &
               'Error in sparse_read_matrix: NGWF number mismatch - file &
               &"'//trim(filename)//'" specifies this many NGWFs ', file_num)
       end if
       call h5aclose_f(attr_id,ierr)
       adims=(/library(ilib)%nblk/)
       call h5aopen_f(file_id,'block_sizes',attr_id,ierr)
       allocate(block_sizes(library(ilib)%nblk),stat=ierr)
       call utils_alloc_check('sparse_read_matrix','block_sizes',ierr)
       call h5aread_f(attr_id,h5t_native_integer,block_sizes,adims,ierr)
       do iblk=1,library(ilib)%nblk
          iat_orig = row_par%orig_atom(iblk)
          call utils_assert(block_sizes(iat_orig) == &
               row_par%num_elems_on_atom(iblk, row_blks), 'Error in sparse_read_&
               &vector: incompatible structure in file "'//trim(filename)&
               //'".')
       end do
       deallocate(block_sizes,stat=ierr)
       call utils_dealloc_check('sparse_read_matrix','block_sizes',ierr)
       call h5aclose_f(attr_id,ierr)
#else
       ! We only use input unit if it is present and its value is valid.
       if ( present(unit) .and. (.not. loc_get_unit) ) then
          iunit = unit
       else
          iunit = utils_unit()
          if (.not.read_bigendian) then
             open(unit=iunit,file=trim(filename),form='unformatted', &
                  action='read',iostat=ierr,status='old')
          else

#ifndef F2008
             ! jd: 'convert' is non-standard Fortran.
             write(stdout,*) 'WARNING: Over-riding data format to big-endian &
                  &when reading matrix'
             open(unit=iunit,file=trim(filename),form='unformatted', &
                  action='read',iostat=ierr,status='old',convert='big_endian')
#else
             call utils_feature_not_supported('sparse_read_matrix:'//&
                  'Big-endian input has not been compiled in because it &
                  &needs "convert", which is not part of the F2008 standard, &
                  &and you asked for an F2008-compliant binary with -D2008.')
#endif

          end if
          call utils_assert(ierr == 0, 'Error in sparse_read_matrix: &
               &opening file "'//trim(filename)//'" failed with code ',ierr)
       end if

       ! Initialise some variables which should be set within the file
       file_nat  = -1
       file_natr = -1
       file_natc = -1
       file_num  = -1
#endif
    else
       iunit = -1 ! qoh: Initialise to avoid compiler warning
    end if

    ! If required, assign value of return value of unit.
    if (loc_get_unit) then
       if (present(unit)) then
          unit = iunit
       else
          call utils_abort('Error in sparse_read_matrix: &
               &asked to retrieve unit but unit argument missing.')
       end if
    end if

    ! Global preliminaries
    ! Allocate temporary array to store block
    max_rows_on_atom = maxval(row_par%num_elems_on_atom(:,row_blks))
    max_cols_on_atom = maxval(col_par%num_elems_on_atom(:,col_blks))
    if (mat(1,1)%iscmplx) then
       allocate(zblk(max_rows_on_atom,max_cols_on_atom),stat=ierr)
       call utils_alloc_check('sparse_read_matrix', &
            'zblk',ierr)
    else
       allocate(dblk(max_rows_on_atom,max_cols_on_atom),stat=ierr)
       call utils_alloc_check('sparse_read_matrix', &
            'dblk',ierr)
    end if
    allocate(ibuf(2*pub_total_num_procs),stat=ierr)
    call utils_alloc_check('sparse_read_matrix','ibuf',ierr)

    do icomp1=1,ncomps(1)
       do icomp2=1,ncomps(2)
          if (mat(1,1)%iscmplx) then
             mat(icomp1,icomp2)%zmtx = (0.0_DP,0.0_DP)
          else
             mat(icomp1,icomp2)%dmtx = 0.0_DP
          end if
       end do
    end do

    endfile_found = .false.

    ! This is what the root proc does: reads from file and sends data to
    ! procs as appropriate
    if (pub_on_root) then

       size_idxinfo = 0 ! qoh: Initialise to prevent compiler warning
#ifdef HDF5
       adims=(/1/)
       call h5aopen_f(file_id,'num_chunks',attr_id,ierr)
       call h5aread_f(attr_id,h5t_native_integer,nchunks,adims,ierr)
       call h5aclose_f(attr_id,ierr)

       do ichunk=1,nchunks
          write(str1,*) ichunk
          dsetname="chunk_"//trim(adjustl(str1))
          call h5gopen_f(file_id,adjustl(dsetname),group_id,ierr)

          call h5aopen_f(group_id,'index_dim',attr_id,ierr)
          call h5aread_f(attr_id,h5t_native_integer,idxlen,adims,ierr)
          if (allocated(idxinfo)) then
             if (size_idxinfo < idxlen) then
                deallocate(idxinfo,stat=ierr)
                call utils_dealloc_check('sparse_read_matrix','idxinfo',ierr)
             end if
          end if
          if (.not. allocated(idxinfo)) then
             allocate(idxinfo(idxlen),stat=ierr)
             call utils_alloc_check('sparse_read_matrix','idxinfo',ierr)
             size_idxinfo = idxlen
          end if
          call h5aclose_f(attr_id,ierr)
          dims=(/idxlen/)
          call h5dopen_f(group_id,'index',dset_id,ierr)
          call h5dread_f(dset_id,h5t_native_integer,idxinfo(1:idxlen),dims,ierr)
          call h5dclose_f(dset_id,ierr)

          call h5aopen_f(group_id,'data_dim',attr_id,ierr)
          call h5aread_f(attr_id,h5t_native_integer,datlen,adims,ierr)
          call utils_assert(mod(datlen,tot_ncomps) == 0, 'Error in sparse_read_&
                &vector: length of DATA record in file "'//trim(filename)//&
                '" is incompatible with number of components.')
          datlen = datlen / tot_ncomps
          call h5aclose_f(attr_id,ierr)
          if (mat(1,1)%iscmplx) then
             if (allocated(zbuf)) then
                if (size(zbuf) < datlen*tot_ncomps) then
                   deallocate(zbuf,stat=ierr)
                   call utils_dealloc_check('sparse_read_matrix','zbuf',ierr)
                end if
             end if
             if (.not. allocated(zbuf)) then
                allocate(zbuf(datlen*tot_ncomps),stat=ierr)
                call utils_alloc_check('sparse_read_matrix','zbuf',ierr)
             end if
             dims=(/datlen*tot_ncomps/)
             call h5dopen_f(group_id,'data_real',dset_id,ierr)
             allocate(dbuf(2*datlen*tot_ncomps),stat=ierr)
             call utils_alloc_check('sparse_read_matrix','dbuf',ierr)
             call h5dread_f(dset_id,h5t_native_double,dbuf(1:datlen*tot_ncomps), &
                  dims,ierr)
             call h5dclose_f(dset_id,ierr)
             call h5dopen_f(group_id,'data_imag',dset_id,ierr)
             call h5dread_f(dset_id,h5t_native_double, &
                  dbuf(datlen*tot_ncomps+1:2*datlen*tot_ncomps),dims,ierr)
             do iel=1,datlen*tot_ncomps
                zbuf(iel)=cmplx(dbuf(iel),dbuf(datlen*tot_ncomps+iel),DP)
             end do
             deallocate(dbuf,stat=ierr)
             call utils_dealloc_check('sparse_read_matrix','dbuf',ierr)
             call h5dclose_f(dset_id,ierr)
          else
             if (allocated(dbuf)) then
                if (size(dbuf) < datlen*tot_ncomps) then
                   deallocate(dbuf,stat=ierr)
                   call utils_dealloc_check('sparse_read_matrix','dbuf',ierr)
                end if
             end if
             if (.not. allocated(dbuf)) then
                allocate(dbuf(datlen*tot_ncomps),stat=ierr)
                call utils_alloc_check('sparse_read_matrix','dbuf',ierr)
             end if
             dims=(/datlen*tot_ncomps/)
             call h5dopen_f(group_id,'data',dset_id,ierr)
             call h5dread_f(dset_id,h5t_native_double,dbuf(1:datlen*tot_ncomps), &
                  dims,ierr)
             call h5dclose_f(dset_id,ierr)
          end if

          call h5gclose_f(group_id,ierr)
#else
       ! Loop until end of file reached
       do

          ! Read tag and work out what to do with record
          read(iunit) tag

          select case (trim(tag(1:11)))

             ! File version information
          case ('VERSION')
             read(iunit) idum
             call utils_assert(idum == 1, 'Error in sparse_read_matrix: &
                  &malformed VERSION record in file "'//trim(filename)//'".')
             read(iunit) idum
             if (idum > nint(file_version*100)) then
                write(stdout,'(a,f5.2,3a,f5.2,a)') &
                     'WARNING in sparse_read_matrix: version', &
                     real(idum,kind=DP)*0.01_DP,' of file "',trim(filename), &
                     '" is more recent than implemented version ', &
                     file_version,': continuing'
             end if

             ! Number of atoms (square matrices only)
          case ('NUM_ATOMS')
             read(iunit) idum
             call utils_assert(idum == 1, 'Error in sparse_read_matrix: &
                  &malformed NUM_ATOMS record in file "'//trim(filename)//'".')
             read(iunit) file_nat
             call utils_assert(file_nat == library(ilib)%nblk, &
                  'Error in sparse_read_matrix: atom number mismatch - file &
                  &"'//trim(filename)//'" specifies a different number of atoms &
                  &than expected by the library', file_nat, library(ilib)%nblk)

              ! Number of atoms (rows, if the matrix isn't square)
          case ('NUM_ATOMS_R')
             read(iunit) idum
             call utils_assert(idum == 1, 'Error in sparse_read_matrix: &
                  &malformed NUM_ATOMS_R record in file "'//trim(filename)//'".')
             read(iunit) file_natr
             call utils_assert(file_natr == library(ilib)%nblk, &
                  'Error in sparse_read_matrix: atom number mismatch - file &
                  &"'//trim(filename)//'" specifies nat=', file_natr)

              ! Number of atoms (cols, if the matrix isn't square)
          case ('NUM_ATOMS_C')
             read(iunit) idum
             call utils_assert(idum == 1, 'Error in sparse_read_matrix: &
                  &malformed NUM_ATOMS_C record in file "'//trim(filename)//'".')
             read(iunit) file_natc
             call utils_assert(file_natc == library(ilib)%mblk, &
                  'Error in sparse_read_matrix: atom number mismatch - file &
                  &"'//trim(filename)//'" specifies nat=', file_natc)

             ! Number of NGWFs
          case ('NUM_NGWFS')
             read(iunit) idum
             call utils_assert(idum == 1, 'Error in sparse_read_matrix: &
                  &malformed NUM_NGWFS record in file "'//trim(filename)//'".')
             read(iunit) file_num
             if(pub_use_swx) then
                call utils_assert(file_num == library(ilib)%nrows, &
                     'Error in sparse_read_matrix: NGWF (or other basis &
                     &function) number mismatch - file &
                     &"'//trim(filename)//'" specifies this many NGWFs &
                     &(or other basis functions)', file_num)
             else
                call utils_assert(file_num == library(ilib)%nrows, &
                     'Error in sparse_read_matrix: NGWF number mismatch - file &
                     &"'//trim(filename)//'" specifies this many NGWFs ', file_num)
             end if
             ! Number of Rows
          case ('NUM_ROWS')
             read(iunit) idum
             call utils_assert(idum == 1, 'Error in sparse_read_vector: &
                  &malformed NUM_ROWS record in file "'//trim(filename)//'".')
             read(iunit) file_num
             call utils_assert(file_num == library(ilib)%nrows, &
                  'Error in sparse_read_vector: Row number mismatch - file &
                  &"'//trim(filename)//'" specifies this many Rows ', file_num)

             ! Number of Cols
          case ('NUM_COLS')
             read(iunit) idum
             call utils_assert(idum == 1, 'Error in sparse_read_vector: &
                  &malformed NUM_COLS record in file "'//trim(filename)//'".')
             read(iunit) file_num
             call utils_assert(file_num == library(ilib)%mcols, &
                  'Error in sparse_read_vector: Col number mismatch - file &
                  &"'//trim(filename)//'" specifies this many Cols ', file_num)

             ! Data type of matrix
          case ('TYPE')
             read(iunit) idum
             call utils_assert(idum == 1, 'Error in sparse_read_matrix: &
                  &malformed TYPE record in file "'//trim(filename)//'"')
             read(iunit) idum
             call utils_assert(idum /= -1 .or. mat(1,1)%iscmplx, &
                  'Error in sparse_read_matrix: data type mismatch - file "'//&
                  trim(filename)//'" specifies complex data but expecting real &
                  &data.')
             call utils_assert(idum /= 0 .or. .not. mat(1,1)%iscmplx, &
                  'Error in sparse_read_matrix: data type mismatch - file "'//&
                  trim(filename)//'" specifies real data but expecting complex &
                  &data.')

             ! Number of components
          case ('NUM_COMPS')
             read(iunit) idum
             call utils_assert(idum == 1, 'Error in sparse_read_matrix: &
                  &malformed NUM_COMPS record in file "'//trim(filename)//'".')
             read(iunit) idum
             call utils_assert(idum == tot_ncomps, 'Error in sparse_read_vector: &
                  &component number mismatch - file "'//trim(filename)//&
                  '" specifies idum components but expecting ncomps. &
                  &Values for idum and ncomps follow. ', idum, tot_ncomps)

             ! Block sizes (number of NGWFs on each atom)
          case ('BLOCK_SIZES')
             ! Error since this should only be found for a square matrix
             call utils_assert(square, 'Error in sparse_read_matrix: &
                  &found BLOCK_SIZES record but matrix is not square.')
             read(iunit) idum
             if (file_nat < 0) file_nat = idum
             call utils_assert(file_nat == idum, 'Error in sparse_read_matrix: &
                  &inconsistent atom number data in file "'//trim(filename)//&
                  '".')
             call utils_assert(file_nat == library(ilib)%nblk, &
                  'Error in sparse_read_matrix: atom number mismatch - file &
                  &"'//trim(filename)//'" specifies this many atoms ', file_nat)
             allocate(block_sizes(library(ilib)%nblk),stat=ierr)
             call utils_alloc_check('sparse_read_matrix','block_sizes',ierr)
             read(iunit) block_sizes
             do iblk=1,library(ilib)%nblk
                iat_orig = row_par%orig_atom(iblk)
                call utils_assert(block_sizes(iat_orig) == &
                     row_par%num_elems_on_atom(iblk, row_blks), 'Error in sparse_read_&
                     &vector: incompatible structure in file "'//trim(filename)&
                     //'".')
             end do
             deallocate(block_sizes,stat=ierr)
             call utils_dealloc_check('sparse_read_matrix','block_sizes',ierr)

          case ('BLOCK_ROWS')
                ! Error since this should not be found for a square matrix
             call utils_assert(.not. square, 'Error in sparse_read_matrix: &
                  &found BLOCK_ROWS record but matrix is square')
             read(iunit) idum
             if (file_natr < 0) file_natr = idum
             call utils_assert(file_natr == idum, 'Error in sparse_read_matrix: &
                  &inconsistent atom number data in file "'//trim(filename)//&
                  '".')
             call utils_assert(file_natr == library(ilib)%nblk, &
                  'Error in sparse_read_matrix: atom number mismatch - file &
                  &"'//trim(filename)//'" specifies this many atoms ', file_natr)
             allocate(block_sizes(library(ilib)%nblk),stat=ierr)
             call utils_alloc_check('sparse_read_matrix','block_sizes',ierr)
             read(iunit) block_sizes
             do iblk=1,library(ilib)%nblk
                iat_orig = row_par%orig_atom(iblk)
                call utils_assert(block_sizes(iat_orig) == &
                     row_par%num_elems_on_atom(iblk, row_blks), 'Error in sparse_read_&
                     &vector: incompatible structure in file "'//trim(filename)&
                     //'".')
             end do
             deallocate(block_sizes,stat=ierr)
             call utils_dealloc_check('sparse_read_matrix','block_sizes',ierr)

          case ('BLOCK_COLS')
                ! Error since this should not be found for a square matrix
             call utils_assert(.not. square, 'Error in sparse_read_matrix: &
                  &found BLOCK_COLS record but matrix is square.')
             read(iunit) idum
             if (file_natc < 0) file_natc = idum
             call utils_assert(file_natc == idum, 'Error in sparse_read_matrix: &
                  &inconsistent atom number data in file "'//trim(filename)//&
                  '".')
             call utils_assert(file_natc == library(ilib)%mblk, &
                  'Error in sparse_read_matrix: atom number mismatch - file &
                  &"'//trim(filename)//'" specifies this many atoms ', file_natc)
             allocate(block_sizes(library(ilib)%mblk),stat=ierr)
             call utils_alloc_check('sparse_read_matrix','block_sizes',ierr)
             read(iunit) block_sizes
             do iblk=1,library(ilib)%mblk
                iat_orig = col_par%orig_atom(iblk)
                call utils_assert(block_sizes(iat_orig) == &
                     col_par%num_elems_on_atom(iblk, col_blks), 'Error in sparse_read_&
                     &vector: incompatible structure in file "'&
                     //trim(filename)//'".')
             end do
             deallocate(block_sizes,stat=ierr)
             call utils_dealloc_check('sparse_read_matrix','block_sizes',ierr)

             ! Block data
          case ('INDEX')

             ! Do we have all relevant info?
             if (file_nat < 0) then
                write(stdout,'(3a)') 'WARNING in sparse_read_matrix: &
                     &file "',trim(filename), &
                     '" does not specify number of atoms: continuing'
                file_nat = library(ilib)%nblk
             end if
             if (file_num < 0) then
                write(stdout,'(3a)') 'WARNING in sparse_read_matrix: &
                     &file "',trim(filename), &
                     '" does not specify number of NGWFs: continuing'
                file_num = library(ilib)%nrows
             end if

             ! Read in index information for data
             read(iunit) idxlen
             if (allocated(idxinfo)) then
                if (size_idxinfo < idxlen) then
                   deallocate(idxinfo,stat=ierr)
                   call utils_dealloc_check('sparse_read_matrix','idxinfo',ierr)
                end if
             end if
             if (.not. allocated(idxinfo)) then
                allocate(idxinfo(idxlen),stat=ierr)
                call utils_alloc_check('sparse_read_matrix','idxinfo',ierr)
                size_idxinfo = idxlen
             end if
             read(iunit) idxinfo(1:idxlen)

             ! Read in data which follows in next record
             read(iunit) cbuf
             call utils_assert(cbuf(1:4) == 'DATA', 'Error in sparse_read_&
                  &vector: no DATA record following INDEX record in file "'//&
                  trim(filename)//'".')
             call utils_assert((cbuf(12:12) /= 'D' .or. .not. mat(1,1)%iscmplx) &
                  .and. (cbuf(12:12) /= 'Z' .or. mat(1,1)%iscmplx), 'Error in &
                  &sparse_read_matrix: incompatible data type in DATA block in &
                  &file "'//trim(filename)//'".')
             read(iunit) datlen
             call utils_assert(mod(datlen,tot_ncomps) == 0, 'Error in sparse_read_&
                  &vector: length of DATA record in file "'//trim(filename)//&
                  '" is incompatible with number of components.')
             datlen = datlen / tot_ncomps
             if (mat(1,1)%iscmplx) then
                if (allocated(zbuf)) then
                   if (size(zbuf) < datlen*tot_ncomps) then
                      deallocate(zbuf,stat=ierr)
                      call utils_dealloc_check('sparse_read_matrix','zbuf',ierr)
                   end if
                end if
                if (.not. allocated(zbuf)) then
                   allocate(zbuf(datlen*tot_ncomps),stat=ierr)
                   call utils_alloc_check('sparse_read_matrix','zbuf',ierr)
                end if
                read(iunit) zbuf(1:datlen*tot_ncomps)
             else
                if (allocated(dbuf)) then
                   if (size(dbuf) < datlen*tot_ncomps) then
                      deallocate(dbuf,stat=ierr)
                      call utils_dealloc_check('sparse_read_matrix','dbuf',ierr)
                   end if
                end if
                if (.not. allocated(dbuf)) then
                   allocate(dbuf(datlen*tot_ncomps),stat=ierr)
                   call utils_alloc_check('sparse_read_matrix','dbuf',ierr)
                end if
                read(iunit) dbuf(1:datlen*tot_ncomps)
             end if
#endif

             ! Work out how much data must be sent to each proc
             ibuf(1:2*pub_total_num_procs) = 0
             do iidx=1,idxlen,3
                iat_orig = idxinfo(iidx+1)
                proc = col_par%proc_of_atom(iat_orig)
                ibuf(2*proc+1) = ibuf(2*proc+1) + 3
                ibuf(2*proc+2) = ibuf(2*proc+2) + idxinfo(iidx+2)
             end do
             call comms_bcast(pub_root_proc_id,ibuf,2*pub_total_num_procs)
             call comms_barrier

             ! Allocate buffers for communication
             if (allocated(idxbuf)) then
                if (size(idxbuf) < idxlen) then
                   deallocate(idxbuf,stat=ierr)
                   call utils_dealloc_check('sparse_read_matrix','idxbuf',ierr)
                end if
             end if
             if (.not. allocated(idxbuf)) then
                allocate(idxbuf(idxlen),stat=ierr)
                call utils_alloc_check('sparse_read_matrix','idxbuf',ierr)
             end if
             if (mat(1,1)%iscmplx) then
                if (allocated(zmtxbuf)) then
                   if (size(zmtxbuf) < datlen*tot_ncomps) then
                      deallocate(zmtxbuf,stat=ierr)
                      call utils_dealloc_check('sparse_read_matrix', &
                           'zmtxbuf',ierr)
                   end if
                end if
                if (.not. allocated(zmtxbuf)) then
                   allocate(zmtxbuf(datlen*tot_ncomps),stat=ierr)
                   call utils_alloc_check('sparse_read_matrix','zmtxbuf',ierr)
                end if
             else
                if (allocated(dmtxbuf)) then
                   if (size(dmtxbuf) < datlen*tot_ncomps) then
                      deallocate(dmtxbuf,stat=ierr)
                      call utils_dealloc_check('sparse_read_matrix', &
                           'dmtxbuf',ierr)
                   end if
                end if
                if (.not. allocated(dmtxbuf)) then
                   allocate(dmtxbuf(datlen*tot_ncomps),stat=ierr)
                   call utils_alloc_check('sparse_read_matrix','dmtxbuf',ierr)
                end if
             end if

             ! Deposit data which is stored on root proc (no comms needed)
             iptr = 1
             do iidx=1,idxlen,3
                iat_orig = idxinfo(iidx+1)
                blksize = idxinfo(iidx+2)
                if (col_par%proc_of_atom(iat_orig) == pub_root_proc_id) then
                   iblk = col_par%distr_atom(iat_orig)
                   jat_orig = idxinfo(iidx)
                   jblk = row_par%distr_atom(jat_orig)
                   jelems = row_par%num_elems_on_atom(jblk,row_blks)
                   comp0 = 0
                   do icomp1=1,ncomps(1)
                      do icomp2=1,ncomps(2)
                         if (mat(1,1)%iscmplx) then
                            do ielem=1,col_par%num_elems_on_atom(iblk,col_blks)
                               zblk(1:jelems,ielem) = zbuf(comp0+iptr+&
                                    (ielem-1)*jelems:comp0+iptr+ielem*jelems-1)
                            end do
                            call sparse_put_block(zblk,mat(icomp1,icomp2),jblk,iblk)
                         else
                            do ielem=1,col_par%num_elems_on_atom(iblk,col_blks)
                               dblk(1:jelems,ielem) = dbuf(comp0+iptr+&
                                    (ielem-1)*jelems:comp0+iptr+ielem*jelems-1)
                            end do
                            call sparse_put_block(dblk,mat(icomp1,icomp2),jblk,iblk)
                         end if
                         comp0 = comp0 + datlen
                      end do
                   end do
                end if
                iptr = iptr + blksize
             end do

             ! Now deal with remaining procs (comms required)
             idxptr = 1
             datptr = 1
             do proc=1,pub_total_num_procs-1
                if (ibuf(2*proc+1) == 0) cycle
                iptr = 1
                idxptr0 = idxptr
                datptr0 = datptr
                do iidx=1,idxlen,3
                   iat_orig = idxinfo(iidx+1)
                   blksize = idxinfo(iidx+2)
                   if (proc == col_par%proc_of_atom(iat_orig)) then
                      idxbuf(idxptr) = idxinfo(iidx)
                      idxbuf(idxptr+1) = iat_orig
                      idxbuf(idxptr+2) = blksize
                      comp0 = 0
                      do icomp1=1,ncomps(1)
                         do icomp2=1,ncomps(2)
                            if (mat(1,1)%iscmplx) then
                               zmtxbuf(datptr:datptr+blksize-1) = &
                                    zbuf(comp0+iptr:comp0+iptr+blksize-1)
                            else
                               dmtxbuf(datptr:datptr+blksize-1) = &
                                    dbuf(comp0+iptr:comp0+iptr+blksize-1)
                            end if
                            comp0 = comp0 + datlen
                            datptr = datptr + blksize
                         end do
                      end do
                      idxptr = idxptr + 3
                   end if
                   iptr = iptr + blksize
                end do
                call comms_send(proc,idxbuf(idxptr0),idxptr-idxptr0)
                if (mat(1,1)%iscmplx) then
                   call comms_send(proc,zmtxbuf(datptr0),datptr-datptr0)
                else
                   call comms_send(proc,dmtxbuf(datptr0),datptr-datptr0)
                end if
             end do

#ifdef HDF5
       end do
#else
             ! End of file
          case ('ENDFILE')
#endif
             ibuf(1:2*pub_total_num_procs) = -1
             call comms_bcast(pub_root_proc_id,ibuf,2*pub_total_num_procs)
             endfile_found = .true.

#ifndef HDF5
             ! Unknown tag
          case default
             if (tag(12:12) == 'I' .or. tag(12:12) == 'D' .or. &
                  tag(12:12) == 'Z') then
                write(stdout,'(5a)') 'WARNING in sparse_read_matrix: &
                     &unknown record tag ',trim(tag(1:11)),' found in file "', &
                     trim(filename),'": continuing'
             else
                call utils_abort('Error in sparse_read_matrix: &
                     &malformed record tag '//trim(tag(1:11))//&
                     ' found in file "'//trim(filename)//'".')
             end if
             read(iunit) idum
             call utils_assert(idum > 0, 'Error in sparse_read_matrix: &
                  &malformed record length '//trim(tag(1:11))//&
                     ' found in file "'//trim(filename)//'".')

             if (tag(12:12) == 'I') then
                if (allocated(ibuf)) then
                   if (size(ibuf) < idum) then
                      deallocate(ibuf,stat=ierr)
                      call utils_dealloc_check('sparse_read_matrix','ibuf',ierr)
                   end if
                end if
                if (.not. allocated(ibuf)) then
                   allocate(ibuf(idum),stat=ierr)
                   call utils_alloc_check('sparse_read_matrix','ibuf',ierr)
                end if
                read(iunit) ibuf(1:idum)
             else if (tag(12:12) == 'D') then
                if (allocated(dbuf)) then
                   if (size(dbuf) < idum) then
                      deallocate(dbuf,stat=ierr)
                      call utils_dealloc_check('sparse_read_matrix','dbuf',ierr)
                   end if
                end if
                if (.not. allocated(dbuf)) then
                   allocate(dbuf(idum),stat=ierr)
                   call utils_alloc_check('sparse_read_matrix','dbuf',ierr)
                end if
                read(iunit) dbuf(1:idum)
             else if (tag(12:12) == 'Z') then
                if (allocated(zbuf)) then
                   if (size(zbuf) < idum) then
                      deallocate(zbuf,stat=ierr)
                      call utils_dealloc_check('sparse_read_matrix','zbuf',ierr)
                   end if
                end if
                if (.not. allocated(zbuf)) then
                   allocate(zbuf(idum),stat=ierr)
                   call utils_alloc_check('sparse_read_matrix','zbuf',ierr)
                end if
                read(iunit) zbuf(1:idum)
             end if
          end select

          if (endfile_found) exit
       end do
#endif

       ! Other procs wait for information from root proc
    else

       ! Loop until end-of-file signal is received
       do
          call comms_bcast(pub_root_proc_id,ibuf,2*pub_total_num_procs)
          idxlen = ibuf(2*pub_my_proc_id+1)
          datlen = ibuf(2*pub_my_proc_id+2)
          if (idxlen == -1) exit  ! signal that end of file reached
          call comms_barrier
          if (idxlen > 0) then
             if (allocated(idxbuf)) then
                if (size(idxbuf) < idxlen) then
                   deallocate(idxbuf,stat=ierr)
                   call utils_dealloc_check('sparse_read_matrix','idxbuf',ierr)
                end if
             end if
             if (.not. allocated(idxbuf)) then
                allocate(idxbuf(idxlen),stat=ierr)
                call utils_alloc_check('sparse_read_matrix','idxbuf',ierr)
             end if
             call comms_recv(pub_root_proc_id,idxbuf,idxlen)
             if (mat(1,1)%iscmplx) then
                if (allocated(zmtxbuf)) then
                   if (size(zmtxbuf) < datlen*tot_ncomps) then
                      deallocate(zmtxbuf,stat=ierr)
                      call utils_dealloc_check('sparse_read_matrix', &
                           'zmtxbuf',ierr)
                   end if
                end if
                if (.not. allocated(zmtxbuf)) then
                   allocate(zmtxbuf(datlen*tot_ncomps),stat=ierr)
                   call utils_alloc_check('sparse_read_matrix','zmtxbuf',ierr)
                end if
                call comms_recv(pub_root_proc_id,zmtxbuf,datlen*tot_ncomps)
             else
                if (allocated(dmtxbuf)) then
                   if (size(dmtxbuf) < datlen*tot_ncomps) then
                      deallocate(dmtxbuf,stat=ierr)
                      call utils_dealloc_check('sparse_read_matrix', &
                           'dmtxbuf',ierr)
                   end if
                end if
                if (.not. allocated(dmtxbuf)) then
                   allocate(dmtxbuf(datlen*tot_ncomps),stat=ierr)
                   call utils_alloc_check('sparse_read_matrix','dmtxbuf',ierr)
                end if
                call comms_recv(pub_root_proc_id,dmtxbuf,datlen*tot_ncomps)
             end if

             iptr = 1
             do iidx=1,idxlen,3
                jat_orig = idxbuf(iidx)
                jblk = row_par%distr_atom(jat_orig)
                iat_orig = idxbuf(iidx+1)
                iblk = col_par%distr_atom(iat_orig)
                blksize = idxbuf(iidx+2)
                jelems = row_par%num_elems_on_atom(jblk,row_blks)
                do icomp1=1,ncomps(1)
                   do icomp2=1,ncomps(2)
                      if (mat(1,1)%iscmplx) then
                         do ielem=1,col_par%num_elems_on_atom(iblk,col_blks)
                            zblk(1:jelems,ielem) = &
                                 zmtxbuf(iptr+(ielem-1)*jelems:iptr+ielem*jelems-1)
                         end do
                         call sparse_put_block(zblk,mat(icomp1,icomp2),jblk,iblk)
                      else
                         do ielem=1,col_par%num_elems_on_atom(iblk,col_blks)
                            dblk(1:jelems,ielem) = &
                                 dmtxbuf(iptr+(ielem-1)*jelems:iptr+ielem*jelems-1)
                         end do
                         call sparse_put_block(dblk,mat(icomp1,icomp2),jblk,iblk)
                      end if
                      iptr = iptr + blksize
                   end do
                end do

             end do

          end if

       end do
    end if

    ! Close file for input
    if (pub_on_root .and. (.not. present(unit))) then
#ifdef HDF5
       call h5gclose_f(file_id,ierr)
       call h5fclose_f(filesystem_id,ierr)
       call h5close_f(ierr)
#else
       close(unit=iunit,iostat=ierr)
       call utils_assert(ierr == 0, &
            'Error in sparse_read_matrix: closing file "'//trim(filename)//&
            '" failed with code ',ierr)
#endif
    end if

    ! Re-sync before continuing
    call comms_free
    call comms_barrier

    ! Global finalisation
    if (allocated(zmtxbuf)) then
       deallocate(zmtxbuf,stat=ierr)
       call utils_dealloc_check('sparse_read_matrix','zmtxbuf',ierr)
    end if
    if (allocated(dmtxbuf)) then
       deallocate(dmtxbuf,stat=ierr)
       call utils_dealloc_check('sparse_read_matrix','dmtxbuf',ierr)
    end if
    if (allocated(idxbuf)) then
       deallocate(idxbuf,stat=ierr)
       call utils_dealloc_check('sparse_read_matrix','idxbuf',ierr)
    end if
    if (allocated(idxinfo)) then
       deallocate(idxinfo,stat=ierr)
       call utils_dealloc_check('sparse_read_matrix','idxinfo',ierr)
    end if
    if (allocated(zbuf)) then
       deallocate(zbuf,stat=ierr)
       call utils_dealloc_check('sparse_read_matrix','zbuf',ierr)
    end if
    if (allocated(dbuf)) then
       deallocate(dbuf,stat=ierr)
       call utils_dealloc_check('sparse_read_matrix','dbuf',ierr)
    end if
    if (mat(1,1)%iscmplx) then
       deallocate(zblk,stat=ierr)
       call utils_dealloc_check('sparse_read_matrix','zblk',ierr)
    else
       deallocate(dblk,stat=ierr)
       call utils_dealloc_check('sparse_read_matrix','dblk',ierr)
    end if
    deallocate(ibuf,stat=ierr)
    call utils_dealloc_check('sparse_read_matrix','ibuf',ierr)

  end subroutine sparse_read_matrix

!------------------------------------------------------------------------------

  !============================================================================!
  ! This subroutine allocates and initialises the segment type array for a     !
  ! SPAM3 hierarchical-sparsity matrix                                         !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   str      (inout) : The matrix structure to be allocated                  !
  !   my_nze   (inout) : Number of nonzero elements on proc                    !
  !   seg_nzb  (in)    : Number of nonzero blocks in each segment on proc      !
  !   seg_nze  (in)    : Number of nonzero elements in each segment on proc    !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, May 2009                                         !
  !============================================================================!

  subroutine sparse_segments_alloc(str,my_nze,seg_nze,seg_nzb,row_par,col_par)

    use comms, only: comms_alltoall, pub_my_proc_id, pub_total_num_procs, &
         pub_num_comms_groups, pub_comms_group_size
    use rundat, only: pub_dense_threshold
    use utils, only: utils_alloc_check

    implicit none

    ! Arguments
    type(STRUC3), intent(inout) :: str ! The matrix structure
    integer, intent(inout) :: my_nze   ! Number of nonzero elements on this proc
    integer, intent(inout) :: seg_nze(0:pub_total_num_procs-1)
    integer, intent(in) :: seg_nzb(0:pub_total_num_procs-1)
    type(PARAL_INFO), pointer, intent(in) :: row_par ! Parallel strategy for rows.
    type(PARAL_INFO), pointer, intent(in) :: col_par ! Parallel strategy for columns.

    ! Local variables
    integer :: ierr                    ! Error flag
    integer :: proc                    ! Proc counter
    integer :: groupfirstproc,grouplastproc
    integer :: group                   ! Comms group counter
    real(kind=DP) :: seg_density       ! Density of the current segment
    integer :: seg_blks                ! Total number of blocks in this segment

    ! Allocate the segment info array
    allocate(str%seg_info(3,0:pub_total_num_procs),stat=ierr)
    call utils_alloc_check('sparse_segments_alloc','str%seg_info('// &
         trim(str%structure)//')',ierr)

    ! Decide on which segments are dense, sparse and blank
    do proc=0,pub_total_num_procs-1

       seg_blks = row_par%num_atoms_on_proc(proc) * &
            col_par%num_atoms_on_proc(pub_my_proc_id)
       ! rc2013: EMBED_FIX!
       if (seg_blks == 0) then
          str%seg_info(s_type,proc) = SEG_BLANK
          cycle ! rc2013: nothing on this proc, move along...
       end if

       ! Calculate the density of nonzero blocks in this segment on this proc
       seg_density = real(seg_nzb(proc),kind=DP)/real(seg_blks,kind=DP)

       ! Determine what type of segment this should be
       if (seg_density >= pub_dense_threshold) then

          ! Segment is dense
          str%seg_info(s_type,proc) = SEG_DENSE

          ! Remove the number of counted nonzero elements from total
          ! and add on the total number in the segment
          my_nze = my_nze - seg_nze(proc)
          seg_nze(proc) = col_par%num_elems_on_proc(pub_my_proc_id,str%col_blks) * &
               row_par%num_elems_on_proc(proc,str%row_blks)
          my_nze = my_nze + seg_nze(proc)

       else if (seg_density > 0.0_DP) then

          ! Segment sparse-indexed
          str%seg_info(s_type,proc) = SEG_SPARSE

       else

          ! Segment contains no elements so can always be skipped
          str%seg_info(s_type,proc) = SEG_BLANK

       end if

    end do

    ! Zero final seg_type (unused)
    str%seg_info(s_type,pub_total_num_procs) = 0

    ! Allocate list of lengths of local proc's segments of index on other procs
    allocate(str%idx_seg_lens(0:pub_total_num_procs-1),stat=ierr)
    call utils_alloc_check('sparse_segments_alloc','str%idx_seg_lens('// &
         trim(str%structure)//')',ierr)

    ! Allocate list of lengths of local proc's segments of data on other procs
    allocate(str%mtx_seg_lens(0:pub_total_num_procs-1),stat=ierr)
    call utils_alloc_check('sparse_segments_alloc','str%mtx_seg_lens('// &
         trim(str%structure)//')',ierr)

    ! Allocate list of lengths of local proc's groups of segments of index
    allocate(str%idx_groupseg_lens(0:pub_total_num_procs-1),stat=ierr)
    call utils_alloc_check('sparse_segments_alloc','str%idx_groupseg_lens('// &
         trim(str%structure)//')',ierr)

    ! Allocate list of lengths of local proc's groups of segments of data
    allocate(str%mtx_groupseg_lens(0:pub_total_num_procs-1),stat=ierr)
    call utils_alloc_check('sparse_segments_alloc','str%mtx_groupseg_lens('// &
         trim(str%structure)//')',ierr)

    ! Find lengths of other proc's segments on local procs
    do proc=0,pub_total_num_procs-1
       if (str%seg_info(s_type,proc)==SEG_BLANK) then
          str%idx_seg_lens(proc) = 0
       else
          str%idx_seg_lens(proc) = col_par%num_atoms_on_proc(pub_my_proc_id) &
               + 1 + seg_nzb(proc)
       end if
       str%mtx_seg_lens(proc) = seg_nze(proc)
    end do

    str%idx_groupseg_lens = 0
    str%mtx_groupseg_lens = 0
    do group=0,pub_num_comms_groups-1
       groupfirstproc = group*pub_comms_group_size
       grouplastproc = (group+1)*pub_comms_group_size - 1
       do proc=groupfirstproc,grouplastproc
          str%idx_groupseg_lens(proc) = sum(str%idx_seg_lens(groupfirstproc: &
               grouplastproc))
          str%mtx_groupseg_lens(proc) = sum(str%mtx_seg_lens(groupfirstproc: &
               grouplastproc))
       end do
    end do

    ! Share with other procs (only need to retain other procs' segment lengths
    ! corresponding to this proc)
    call comms_alltoall(str%idx_seg_lens,1)
    call comms_alltoall(str%mtx_seg_lens,1)
    call comms_alltoall(str%idx_groupseg_lens,1)
    call comms_alltoall(str%mtx_groupseg_lens,1)

  end subroutine sparse_segments_alloc

  !============================================================================!
  ! This subroutine allocates the internal arrays of a block sparse matrix     !
  ! structure.                                                                 !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   str     (inout) : The matrix structure to be allocated                   !
  !   nze     (input) : Number of nonzero elements                             !
  !   nzb     (input) : Number of nonzero blocks                               !
  !   my_nze  (input) : Number of nonzero elements on this processor           !
  !   my_nzb  (input) : Number of nonzero blocks on this processor             !
  !----------------------------------------------------------------------------!
  ! Written by Peter Haynes, March 2004.                                       !
  ! Revised for distributed data, June 2006.                                   !
  ! Revised for SPAM3 by Nicholas Hine, May 2009.                              !
  ! Revised for multiple parallel strategies by Robert Charlton, December 2017.!
  !============================================================================!

  subroutine sparse_struc_alloc(str,nze,nzb,my_nze,my_nzb,row_par,col_par)

    use comms, only: comms_reduce, pub_my_proc_id, pub_total_num_procs
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(STRUC3), intent(inout) :: str ! The matrix structure
    integer(kind=LONG), intent(in) :: nze ! Number of nonzero elements
    integer, intent(in) :: nzb         ! Number of nonzero blocks
    integer, intent(in) :: my_nze      ! Number of nonzero elements on this proc
    integer, intent(in) :: my_nzb      ! Number of nonzero blocks on this proc
    type(PARAL_INFO), pointer :: row_par ! Parallel strategy for rows
    type(PARAL_INFO), pointer :: col_par ! Parallel strategy for columns

    ! Local variables
    integer :: ierr                    ! Error flag
    integer :: nindices                ! Number of non-blank segment indices
    integer :: my_idx_len              ! Index length on this proc
    integer, allocatable :: ibuf(:)    ! Buffer for communication

    ! Store sizes
    str%nrows = sum(row_par%num_elems_on_proc(:,str%row_blks))
    str%mcols = sum(col_par%num_elems_on_proc(:,str%col_blks))
    str%nze = nze
    ! rc2013: no. of blks in strucutre. This could differ for row/col par's
    str%nblk = row_par%nat
    str%mblk = col_par%nat
    str%nzb = nzb
    str%my_nze = my_nze
    str%my_nzb = my_nzb
    str%my_nblks = col_par%my_last_blk - col_par%my_first_blk + 1
    ! rc2013: include parallel strategies for rows and columns in matrix structure
    str%row_par => row_par
    str%col_par => col_par

    nindices = pub_total_num_procs - count(str%seg_info(s_type,:)==SEG_BLANK)

    my_idx_len = (str%my_nblks + 1)*nindices + str%my_nzb

    ! Ensure arrays are allocated to at least size 1
    my_idx_len = max(my_idx_len,1)

    ! Allocate list of index lengths
    allocate(str%idx_lens(0:pub_total_num_procs-1),stat=ierr)
    call utils_alloc_check('sparse_struc_alloc','str%idx_lens('// &
         trim(str%structure)//')',ierr)

    ! Allocate block index
    allocate(str%blk_idx(my_idx_len),stat=ierr)
    call utils_alloc_check('sparse_struc_alloc','str%blk_idx('// &
         trim(str%structure)//')',ierr)

    ! Allocate block pointers
    allocate(str%blk_ptr(my_idx_len+1),stat=ierr)
    call utils_alloc_check('sparse_struc_alloc','str%blk_ptr('// &
         trim(str%structure)//')',ierr)

    ! Allocate buffer for communication
    allocate(ibuf(0:max(2,pub_total_num_procs-1)),stat=ierr)
    call utils_alloc_check('sparse_struc_alloc','ibuf',ierr)

    ! Find maximum sizes
    ibuf = 0
    ibuf(0) = my_nze ; ibuf(1) = str%my_nblks ; ibuf(2) = my_nzb
    call comms_reduce('MAX',ibuf(0:2))
    str%max_nze = ibuf(0)
    str%max_nblks = ibuf(1)
    str%max_nzb = ibuf(2)

    ! Set up list of index lengths
    ibuf = 0
    ibuf(pub_my_proc_id) = my_idx_len
    call comms_reduce('SUM',ibuf)
    str%idx_lens(0:pub_total_num_procs-1) = ibuf(0:pub_total_num_procs-1)

    ! Deallocate buffer for communication
    deallocate(ibuf,stat=ierr)
    call utils_dealloc_check('sparse_struc_alloc','ibuf',ierr)

  end subroutine sparse_struc_alloc

  !============================================================================!
  ! This subroutine copies the internal arrays of a block sparse matrix        !
  ! structure from one structure to another.                                   !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   str_dest (inout) : The matrix structure to be copied into                !
  !   str_src  (input) : Number of nonzero elements                            !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, April 2011.                                      !
  !============================================================================!

  subroutine sparse_struc_copy(str_dest,str_src)

    use comms, only: comms_reduce, pub_total_num_procs
    use utils, only: utils_alloc_check

    implicit none

    ! Arguments
    type(STRUC3), intent(inout) :: str_dest ! The matrix structure
    type(STRUC3), intent(in) :: str_src  ! The matrix structure

    ! Local variables
    integer :: ierr                    ! Error flag
    integer :: my_idx_len              ! Index length on this proc

    ! Store basic info
    str_dest%structure = str_src%structure
    str_dest%nrows = str_src%nrows
    str_dest%mcols = str_src%mcols
    str_dest%nze = str_src%nze
    str_dest%nblk = str_src%nblk
    str_dest%mblk = str_src%mblk
    str_dest%nzb = str_src%nzb
    str_dest%my_nze = str_src%my_nze
    str_dest%my_nzb = str_src%my_nzb
    str_dest%my_nblks = str_src%my_nblks
    ! rc2013: pointers for parallel strategies
    str_dest%row_par=>str_src%row_par
    str_dest%col_par=>str_src%col_par

    ! Allocate the segment info array
    allocate(str_dest%seg_info(3,0:pub_total_num_procs),stat=ierr)
    call utils_alloc_check('sparse_struc_copy','str%seg_info('// &
         trim(str_dest%structure)//')',ierr)
    str_dest%seg_info(:,:) = str_src%seg_info(:,:)

    ! Allocate list of index lengths
    allocate(str_dest%idx_lens(0:pub_total_num_procs-1),stat=ierr)
    call utils_alloc_check('sparse_struc_alloc','str%idx_lens('// &
         trim(str_dest%structure)//')',ierr)
    str_dest%idx_lens = str_src%idx_lens

    ! Allocate block index
    my_idx_len = size(str_src%blk_idx)
    allocate(str_dest%blk_idx(my_idx_len),stat=ierr)
    call utils_alloc_check('sparse_struc_alloc','str%blk_idx('// &
         trim(str_dest%structure)//')',ierr)
    str_dest%blk_idx = str_src%blk_idx

    ! Allocate block pointers
    allocate(str_dest%blk_ptr(my_idx_len+1),stat=ierr)
    call utils_alloc_check('sparse_struc_alloc','str%blk_ptr('// &
         trim(str_dest%structure)//')',ierr)
    str_dest%blk_ptr = str_src%blk_ptr

    str_dest%max_nze = str_src%max_nze
    str_dest%max_nblks = str_src%max_nblks
    str_dest%max_nzb = str_src%max_nzb

  end subroutine sparse_struc_copy

  !============================================================================!
  ! This subroutine allocates the data arrays of a block sparse matrix.        !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   mat     (inout)           : The matrix structure to be allocated         !
  !   iscmplx (optional, input) : TRUE if data is complex, otherwise FALSE     !
  !----------------------------------------------------------------------------!
  ! Written by Peter Haynes, March 2004.                                       !
  ! Revised for distributed data, June 2006.                                   !
  !============================================================================!

  subroutine sparse_data_alloc(mat,iscmplx)

    use constants, only: real_size, cmplx_size
    use utils, only: utils_alloc_check

    implicit none

    ! Arguments
    type(SPAM3), intent(inout) :: mat        ! The matrix structure
    logical, optional, intent(in) :: iscmplx ! Flag to indicate complex data

    ! Local variables
    integer :: ierr          ! Error flag
    logical :: loc_iscmplx   ! Local copy of iscmplx input flag
    integer :: my_nze        ! Length of data arrays for this matrix

    ! Make copy of optional input argument
    if (present(iscmplx)) then
       loc_iscmplx = iscmplx
    else
       loc_iscmplx = .false.
    end if

    ! Find length of data arrays from library
    my_nze = max(library(mat%lib)%my_nze,1)

    ! Allocate space for matrix elements
    if (loc_iscmplx) then

       allocate(mat%zmtx(my_nze),stat=ierr)
       call utils_alloc_check('sparse_data_alloc','mat%zmtx('// &
            trim(library(mat%lib)%structure)//')',ierr)
       mat%zmtx = (0.0_DP,0.0_DP)

    else

       allocate(mat%dmtx(my_nze),stat=ierr)
       call utils_alloc_check('sparse_data_alloc','mat%dmtx('// &
            trim(library(mat%lib)%structure)//')',ierr)
       mat%dmtx = 0.0_DP

    end if

    ! Fill in remaining details
    mat%iscmplx = loc_iscmplx
    mat%structure = library(mat%lib)%structure

    ! Record matrix memory usage
    if (mat%iscmplx) then
       local_mat_mem = local_mat_mem + int(library(mat%lib)%my_nze,kind=LONG)* &
            cmplx_size
       global_mat_mem = global_mat_mem + library(mat%lib)%nze*cmplx_size
    else
       local_mat_mem = local_mat_mem + int(library(mat%lib)%my_nze,kind=LONG)* &
            real_size
       global_mat_mem = global_mat_mem + library(mat%lib)%nze*real_size
    end if

  end subroutine sparse_data_alloc

  !============================================================================!
  ! This subroutine deallocates a block sparse matrix structure.               !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   str     (inout) : The matrix structure to be deallocated                 !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, May 2009.                                        !
  !============================================================================!

  subroutine sparse_segments_dealloc(str)

    use utils, only: utils_dealloc_check

    implicit none

    ! Arguments
    type(STRUC3), intent(inout) :: str        ! The matrix structure

    ! Local variables
    integer :: ierr          ! Error flag

    ! Deallocate list of lengths of local proc's segment of data on other procs
    if (allocated(str%mtx_groupseg_lens)) then
       deallocate(str%mtx_groupseg_lens,stat=ierr)
       call utils_dealloc_check('sparse_segments_dealloc', &
            'str%mtx_groupseg_lens('//trim(str%structure)//')',ierr)
    end if

    ! Deallocate list of local proc's segment of index lengths on other procs
    if (allocated(str%idx_groupseg_lens)) then
       deallocate(str%idx_groupseg_lens,stat=ierr)
       call utils_dealloc_check('sparse_segments_dealloc', &
            'str%idx_groupseg_lens('//trim(str%structure)//')',ierr)
    end if

    ! Deallocate list of lengths of local proc's segment of data on other procs
    if (allocated(str%mtx_seg_lens)) then
       deallocate(str%mtx_seg_lens,stat=ierr)
       call utils_dealloc_check('sparse_segments_dealloc','str%mtx_seg_lens(' &
            //trim(str%structure)//')',ierr)
    end if

    ! Deallocate list of local proc's segment of index lengths on other procs
    if (allocated(str%idx_seg_lens)) then
       deallocate(str%idx_seg_lens,stat=ierr)
       call utils_dealloc_check('sparse_segments_dealloc','str%idx_seg_lens(' &
            //trim(str%structure)//')',ierr)
    end if

    ! Deallocate segment information
    if (allocated(str%seg_info)) then
       deallocate(str%seg_info,stat=ierr)
       call utils_dealloc_check('sparse_segments_dealloc','str%seg_info(' &
            //trim(str%structure)//')',ierr)
    end if

  end subroutine sparse_segments_dealloc


  !============================================================================!
  ! This subroutine deallocates a block sparse matrix structure.               !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   str     (inout) : The matrix structure to be deallocated                 !
  !----------------------------------------------------------------------------!
  ! Written by Peter Haynes, March 2004.                                       !
  ! Revised for distributed data, June 2006.                                   !
  ! Adapted for SPAM3 by Nicholas Hine, May 2009.                              !
  !============================================================================!

  subroutine sparse_struc_dealloc(str)

    use utils, only: utils_dealloc_check

    implicit none

    ! Arguments
    type(STRUC3), intent(inout) :: str        ! The matrix structure

    ! Local variables
    integer :: ierr          ! Error flag

    ! Deallocate block pointers
    if (allocated(str%blk_ptr)) then
       deallocate(str%blk_ptr,stat=ierr)
       call utils_dealloc_check('sparse_struc_dealloc','str%blk_ptr('// &
            trim(str%structure)//')',ierr)
    end if

    ! Deallocate block index
    if (allocated(str%blk_idx)) then
       deallocate(str%blk_idx,stat=ierr)
       call utils_dealloc_check('sparse_struc_dealloc','str%blk_idx('// &
            trim(str%structure)//')',ierr)
    end if

    ! Deallocate list of index lengths
    if (allocated(str%idx_lens)) then
       deallocate(str%idx_lens,stat=ierr)
       call utils_dealloc_check('sparse_struc_dealloc','str%idx_lens('// &
            trim(str%structure)//')',ierr)
    end if

    ! Deallocate segment related arrays
    call sparse_segments_dealloc(str)

  end subroutine sparse_struc_dealloc


  !============================================================================!
  ! This subroutine deallocates the data arrays of a block sparse matrix.      !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   mat     (inout) : The matrix structure to be deallocated                 !
  !----------------------------------------------------------------------------!
  ! Written by Peter Haynes, March 2004.                                       !
  ! Revised for distributed data, June 2006.                                   !
  !============================================================================!

  subroutine sparse_data_dealloc(mat)

    use constants, only: real_size, cmplx_size
    use utils, only: utils_dealloc_check

    implicit none

    ! Arguments
    type(SPAM3), intent(inout) :: mat        ! The matrix structure

    ! Local variables
    integer :: ierr          ! Error flag

    ! Deallocate data arrays
    if (allocated(mat%dmtx)) then
       deallocate(mat%dmtx,stat=ierr)
       call utils_dealloc_check('sparse_data_dealloc','mat%dmtx('// &
            trim(library(mat%lib)%structure)//')',ierr)
    end if

    if (allocated(mat%zmtx)) then
       deallocate(mat%zmtx,stat=ierr)
       call utils_dealloc_check('sparse_data_dealloc','mat%zmtx('// &
            trim(library(mat%lib)%structure)//')',ierr)
    end if

    ! Record matrix memory usage
    if (mat%iscmplx) then
       local_mat_mem = local_mat_mem - int(library(mat%lib)%my_nze,kind=LONG)* &
            cmplx_size
       global_mat_mem = global_mat_mem - library(mat%lib)%nze*cmplx_size
    else
       local_mat_mem = local_mat_mem - int(library(mat%lib)%my_nze,kind=LONG)* &
            real_size
       global_mat_mem = global_mat_mem - library(mat%lib)%nze*real_size
    end if

  end subroutine sparse_data_dealloc


  !==========================================================================!
  ! This subroutine searches the library for an entry which matches the      !
  ! given structure identification string.                                   !
  !--------------------------------------------------------------------------!
  ! Arguments:                                                               !
  !   name (inout) : Structure to search for in library                      !
  !   idx (output) : Entry number in library if found, otherwise 0           !
  !--------------------------------------------------------------------------!
  ! Written by Peter Haynes, March 2004.                                     !
  ! Moved to private module subroutine, Nicholas Hine, Dec 2007              !
  !==========================================================================!

  subroutine sparse_search_library(name,idx)

    implicit none

    ! Arguments
    character(len=*), intent(in) :: name     ! Structure identifier to search
    integer, intent(out) :: idx              ! Index in library if found

    ! Local variables
    integer :: ilib    ! Counter over library entries

    ! Set index to zero
    idx = 0

    ! Loop over library entries
    do ilib=1,num_library

       ! Compare structure identification strings
       if (trim(library(ilib)%structure) == trim(name)) then

          ! This library entry matches
          idx = ilib

          exit
       end if

    end do

  end subroutine sparse_search_library


  !============================================================================!
  ! This routine lists the library in the event that it fills up.              !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   name   (input) : Identifying routine name                                !
  !----------------------------------------------------------------------------!
  ! Written by Peter Haynes, March 2006.                                       !
  !============================================================================!

  subroutine sparse_library_full(name)

    use comms, only: pub_on_root
    use utils, only: utils_abort

    implicit none

    ! Argument
    character(len=*), intent(in) :: name  ! Identifying name

    ! Local variables
    integer :: ilib      ! Library entry

    ! Report error and library information
    if (pub_on_root) then
       write(stdout,'(3a)') 'Error in ',trim(name),': library full'
       write(stdout,'(a,i3)') '  Library size: ',max_library
       write(stdout,'(a,i3)') 'Current number: ',num_library
       write(stdout,'(a)') '  Library contents:'
       do ilib=1,max_library
          write(stdout,'(4x,i3,2a)') ilib,': ', &
               adjustl(library(ilib)%structure)
       end do
    end if

    ! And quit
    call utils_abort('Intentional exit in sparse_library_full().')

  end subroutine sparse_library_full


  !==========================================================================!
  ! This subroutine searches the library for an entry which matches the      !
  ! given structure identification string.                                   !
  !--------------------------------------------------------------------------!
  ! Arguments:                                                               !
  !     name (in)     : Structure to search for in library                   !
  !   exists (output) : Returns 'T' if 'name' is in the library              !
  !--------------------------------------------------------------------------!
  ! Written by Robert Charlton, 01/07/2019.                                  !
  !==========================================================================!

  subroutine sparse_mat_in_library(name,exists)

    implicit none

    ! Arguments
    character(len=*), intent(in) :: name     ! Structure identifier to search
    logical, intent(out) :: exists             ! Index in library if found

    ! Local variables
    integer :: ilib    ! Counter over library entries

    ! Set index to zero
    exists = .false.

    ! Loop over library entries
    do ilib=1,num_library

       ! Compare structure identification strings
       if (trim(library(ilib)%structure) == trim(name)) then

          ! This library entry matches
          exists = .true.

          exit
       end if

    end do

  end subroutine sparse_mat_in_library


  !============================================================================!
  ! This routine allocates a COM3 communicator. Takes multiple options for     !
  ! allocating different arrays, depending on which routine called it.         !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   com   (inout) : COM3 type with communicator arrays.                      !
  !   mat   (input) : SPAM3 matrix.                                            !
  !   nbuf  (input) : number of buffers to allocate.                           !
  !   alloc_mtx (in): whether to allocate the data parts of the communicator.  !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, February 2012.                                   !
  !============================================================================!

  subroutine sparse_com_allocate(com,mat,nbuf,alloc_mtx,cropped,seg,groupseg, &
       cmplx_duplicate)

    use comms, only: pub_total_num_procs, pub_null_handle
    use utils, only: utils_alloc_check

    implicit none

    ! Arguments
    type(COM3), intent(inout) :: com
    type(SPAM3), intent(in) :: mat
    integer, intent(in) :: nbuf
    logical, intent(in) :: alloc_mtx
    logical, intent(in) :: cropped
    logical, intent(in) :: seg
    logical, intent(in) :: groupseg
    logical, intent(in), optional :: cmplx_duplicate

    ! Local Variables
    integer :: ierr
    logical :: loc_cmplx_duplicate

    ! Set up basic info
    com%lib = mat%lib
    com%iscmplx = mat%iscmplx
    com%cropped = cropped
    com%segonly = seg
    com%groupsegonly = groupseg
    com%nbuf = nbuf
    ! rc2013: obtain the parallel strategies info using library strucures.
    com%row_par=>library(mat%lib)%row_par
    com%col_par=>library(mat%lib)%col_par

    if (seg) then
       com%buflen = maxval(library(mat%lib)%idx_seg_lens)
       com%datlen = maxval(library(mat%lib)%mtx_seg_lens)
    else if (groupseg) then
       com%buflen = maxval(library(mat%lib)%idx_groupseg_lens)
       com%datlen = maxval(library(mat%lib)%mtx_groupseg_lens)
    else
       com%buflen = maxval(library(mat%lib)%idx_lens)
       com%datlen = library(mat%lib)%max_nze
    end if

    ! jme: will we create a complex buffer duplicate of the real one?
    ! (if mat%iscmplx this will not be used at all)
    loc_cmplx_duplicate = .false.
    if (present(cmplx_duplicate)) loc_cmplx_duplicate = cmplx_duplicate

    ! Allocate communication buffers
    allocate(com%seginfobuf(3,0:pub_total_num_procs,1:nbuf),stat=ierr)
    call utils_alloc_check('sparse_com_allocate','com%seginfobuf',ierr)
    allocate(com%idxbuf(1:com%buflen,1:nbuf),stat=ierr)
    call utils_alloc_check('sparse_com_allocate','com%idxbuf',ierr)
    allocate(com%ptrbuf(1:com%buflen,1:nbuf),stat=ierr)
    call utils_alloc_check('sparse_com_allocate','com%ptrbuf',ierr)
    if (alloc_mtx) then
       if (com%iscmplx) then
          allocate(com%zmtxrecvbuf(com%datlen,1:nbuf),stat=ierr)
          call utils_alloc_check('sparse_com_allocate','com%zmtxrecvbuf',ierr)
          allocate(com%dmtxrecvbuf(1,2:2),stat=ierr)
          call utils_alloc_check('sparse_com_allocate','com%dmtxrecvbuf',ierr)
       else
          allocate(com%dmtxrecvbuf(com%datlen,1:nbuf),stat=ierr)
          call utils_alloc_check('sparse_com_allocate','com%dmtxrecvbuf',ierr)
          if (loc_cmplx_duplicate) then
             allocate(com%zmtxrecvbuf(com%datlen,1:nbuf),stat=ierr)
             call utils_alloc_check('sparse_com_allocate','com%zmtxrecvbuf',ierr)
          else
             allocate(com%zmtxrecvbuf(1,2:2),stat=ierr)
             call utils_alloc_check('sparse_com_allocate','com%zmtxrecvbuf',ierr)
          end if
       end if
    end if

    ! Allocate comms buffers for cropped send/recv
    if (cropped) then
       allocate(com%ptrreqrecvbuf(4,1:com%buflen),stat=ierr)
       call utils_alloc_check('sparse_com_allocate','com%ptrreqrecvbuf',ierr)
       allocate(com%ptrreqsendbuf(4,1:com%buflen),stat=ierr)
       call utils_alloc_check('sparse_com_allocate','com%ptrreqsendbuf',ierr)
       if (alloc_mtx) then
          if (com%iscmplx) then
             allocate(com%zmtxsendbuf(com%datlen,1:nbuf),stat=ierr)
             call utils_alloc_check('sparse_com_allocate','com%zmtxsendbuf',ierr)
             allocate(com%dmtxsendbuf(1,2:2),stat=ierr)
             call utils_alloc_check('sparse_com_allocate','com%dmtxsendbuf',ierr)
          else
             allocate(com%dmtxsendbuf(com%datlen,1:nbuf),stat=ierr)
             call utils_alloc_check('sparse_com_allocate','com%dmtxsendbuf',ierr)
             allocate(com%zmtxsendbuf(1,2:2),stat=ierr)
             call utils_alloc_check('sparse_com_allocate','com%zmtxsendbuf',ierr)
          end if
       end if
    end if

    ! Allocate request and handle buffers
    com%num_handles = 2*pub_total_num_procs + 6
    allocate(com%index_reqs(0:pub_total_num_procs-1),stat=ierr)
    call utils_alloc_check('sparse_com_allocate','com%index_reqs',ierr)
    allocate(com%data_reqs(0:pub_total_num_procs-1),stat=ierr)
    call utils_alloc_check('sparse_com_allocate','com%data_reqs',ierr)
    allocate(com%handles(1:com%num_handles),stat=ierr)
    call utils_alloc_check('sparse_com_allocate','com%handles',ierr)
    allocate(com%send_buffer_free(1:nbuf),stat=ierr)
    call utils_alloc_check('sparse_com_allocate','com%send_buffer_free',ierr)
    com%index_reqs(:) = index_sent
    com%data_reqs(:) = data_sent
    com%handles(:) = pub_null_handle

  end subroutine sparse_com_allocate


  !============================================================================!
  ! This routine deallocates a COM3 communicator.                              !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   com   (inout) : COM3 type with communicator arrays.                      !
  !   dealloc_mtx (in): whether to deallocate the data parts.                  !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, February 2012.                                   !
  !============================================================================!

  subroutine sparse_com_deallocate(com,dealloc_mtx,cropped)

    use utils, only: utils_dealloc_check

    implicit none

    ! Arguments
    type(COM3), intent(inout) :: com
    logical, intent(in) :: dealloc_mtx
    logical, intent(in) :: cropped

    ! Local Variables
    integer :: ierr

    ! Deallocate request and handle buffers
    deallocate(com%send_buffer_free,stat=ierr)
    call utils_dealloc_check('sparse_com_deallocate','com%send_buffer_free', &
         ierr)
    deallocate(com%handles,stat=ierr)
    call utils_dealloc_check('sparse_com_deallocate','com%handles',ierr)
    deallocate(com%data_reqs,stat=ierr)
    call utils_dealloc_check('sparse_com_deallocate','com%data_reqs',ierr)
    deallocate(com%index_reqs,stat=ierr)
    call utils_dealloc_check('sparse_com_deallocate','com%index_reqs',ierr)

    ! Allocate comms buffers for cropped send/recv
    if (cropped) then
       if (dealloc_mtx) then
          deallocate(com%dmtxsendbuf,stat=ierr)
          call utils_dealloc_check('sparse_com_deallocate','com%dmtxsendbuf',ierr)
          deallocate(com%zmtxsendbuf,stat=ierr)
          call utils_dealloc_check('sparse_com_deallocate','com%zmtxsendbuf',ierr)
       end if
       deallocate(com%ptrreqsendbuf,stat=ierr)
       call utils_dealloc_check('sparse_com_deallocate','com%ptrreqsendbuf',ierr)
       deallocate(com%ptrreqrecvbuf,stat=ierr)
       call utils_dealloc_check('sparse_com_deallocate','com%ptrreqrecvbuf',ierr)
    end if

    ! Deallocate communication buffers
    if (dealloc_mtx) then
       deallocate(com%zmtxrecvbuf,stat=ierr)
       call utils_dealloc_check('sparse_com_deallocate','com%zmtxrecvbuf',ierr)
       deallocate(com%dmtxrecvbuf,stat=ierr)
       call utils_dealloc_check('sparse_com_deallocate','com%dmtxrecvbuf',ierr)
    end if

    deallocate(com%ptrbuf,stat=ierr)
    call utils_dealloc_check('sparse_com_deallocate','com%ptrbuf',ierr)
    deallocate(com%idxbuf,stat=ierr)
    call utils_dealloc_check('sparse_com_deallocate','com%idxbuf',ierr)
    deallocate(com%seginfobuf,stat=ierr)
    call utils_dealloc_check('sparse_com_deallocate','com%seginfobuf',ierr)

  end subroutine sparse_com_deallocate


  !==========================================================================!
  ! This subroutine sends either a segment of the index or the whole of the  !
  ! index of a given sparse matrix structure to another proc                 !
  !--------------------------------------------------------------------------!
  ! Arguments:                                                               !
  !   lib                                                                    !
  !   destproc                                                               !
  !--------------------------------------------------------------------------!
  ! Written by Nicholas Hine, May 2009.                                      !
  !==========================================================================!

  subroutine sparse_send_index(com,destproc,idx_handle,info_handle)

    use comms, only: comms_send, pub_my_proc_id, pub_total_num_procs, &
         pub_null_handle, pub_comms_group_size

    implicit none

    ! Arguments
    type(COM3), intent(in) :: com
    integer, intent(in) :: destproc
    integer, intent(inout), optional :: idx_handle
    integer, intent(inout), optional :: info_handle

    ! Locals
    integer :: buflen, seg_start
    integer :: destgroupfirstproc,destgrouplastproc

    if (pub_my_proc_id==destproc) return

    if (com%segonly) then
       seg_start = library(com%lib)%seg_info(s_idx,destproc)
       buflen = library(com%lib)%seg_info(s_idx,destproc+1) - seg_start
       if (buflen < 1) return
       if (present(idx_handle).and.present(info_handle)) then
          call comms_send(destproc,library(com%lib)%seg_info(:,destproc), &
               3,tag=SEGINFO_TAG,return_handle=info_handle)
          call comms_send(destproc,library(com%lib)%blk_idx(seg_start), &
               buflen,tag=BLKIDX_TAG,return_handle=idx_handle)
       else
          call comms_send(destproc,library(com%lib)%seg_info(:,destproc), &
               3,tag=SEGINFO_TAG)
          call comms_send(destproc,library(com%lib)%blk_idx(seg_start), &
               buflen,tag=BLKIDX_TAG)
       end if
    else if (com%groupsegonly) then
       ! Send info and index of all local segments corresponding to the
       ! procs of destproc's comms group
       destgroupfirstproc = &
            int(destproc/pub_comms_group_size)*pub_comms_group_size
       seg_start = library(com%lib)%seg_info(s_idx,destgroupfirstproc)
       destgrouplastproc = (int(destproc/pub_comms_group_size)+1) * &
            pub_comms_group_size - 1
       buflen = library(com%lib)%seg_info(s_idx,destgrouplastproc+1) - seg_start
       if (buflen < 1) return
       if (present(idx_handle).and.present(info_handle)) then
          call comms_send(destproc,library(com%lib)%seg_info( &
               :,destgroupfirstproc:destgrouplastproc), &
               3*pub_comms_group_size,tag=SEGINFO_TAG,return_handle=info_handle)
          call comms_send(destproc,library(com%lib)%blk_idx(seg_start), &
               buflen,tag=BLKIDX_TAG,return_handle=idx_handle)
       else
          call comms_send(destproc,library(com%lib)%seg_info( &
               :,destgroupfirstproc:destgrouplastproc), &
               3*pub_comms_group_size,tag=SEGINFO_TAG)
          call comms_send(destproc,library(com%lib)%blk_idx(seg_start), &
               buflen,tag=BLKIDX_TAG)
       end if
    else
       buflen = library(com%lib)%idx_lens(pub_my_proc_id)
       if (buflen < 1) then
          if (present(idx_handle).and.present(info_handle)) then
             idx_handle = pub_null_handle
             info_handle = pub_null_handle
          end if
          return
       end if
       if (present(idx_handle).and.present(info_handle)) then
          call comms_send(destproc,library(com%lib)%seg_info(:,:), &
               3*pub_total_num_procs+3,tag=SEGINFO_TAG, &
               return_handle=info_handle)
          call comms_send(destproc,library(com%lib)%blk_idx, &
               buflen,tag=BLKIDX_TAG,return_handle=idx_handle)
       else ! synchronous receives - no need for handles
          call comms_send(destproc,library(com%lib)%seg_info(:,:), &
               3*pub_total_num_procs+3,tag=SEGINFO_TAG)
          call comms_send(destproc,library(com%lib)%blk_idx, &
               buflen,tag=BLKIDX_TAG)
       end if
    end if

  end subroutine sparse_send_index


  !==========================================================================!
  ! This subroutine receives either the whole index or a segment of the      !
  ! index of a given sparse matrix structure from another proc               !
  !--------------------------------------------------------------------------!
  ! Arguments:                                                               !
  !   lib                                                                    !
  !   srcproc                                                                !
  !   idxbuf                                                                 !
  !   seginfobuf                                                             !
  !--------------------------------------------------------------------------!
  ! Written by Nicholas Hine, May 2009.                                      !
  !==========================================================================!

  subroutine sparse_recv_index(com,srcproc,ibuf,async)

    use comms, only: comms_recv, comms_irecv, pub_my_proc_id, &
         pub_total_num_procs, pub_comms_group_size, pub_first_proc_in_group, &
         pub_last_proc_in_group

    implicit none

    ! Arguments
    type(COM3), intent(inout) :: com      ! Sparse Matrix Communicator object
    integer, intent(in) :: srcproc        ! Source proc from which to receive
    integer, intent(in) :: ibuf           ! Which buffer to use
    logical, intent(in) :: async          ! Whether to send asynchronously

    ! Locals
    integer :: buflen      ! Buffer Length
    integer :: seg_start   ! Segment start position

    ! If segonly==.true., we are receiving only the segment of the index
    ! corresponding to the local proc
    if (com%segonly) then

       buflen = library(com%lib)%idx_seg_lens(srcproc)
       if (buflen == 0) then
          com%seginfobuf(s_type,0,ibuf) = SEG_BLANK
          return
       end if
       if (pub_my_proc_id /= srcproc) then
          if (async) then
             call comms_irecv(srcproc,com%seginfobuf(:,0,ibuf),3, &
                  tag=SEGINFO_TAG,handle=com%handles(com%info_hdl))
             call comms_irecv(srcproc,com%idxbuf(:,ibuf),buflen, &
                  tag=BLKIDX_TAG,handle=com%handles(com%idx_hdl))
          else
             call comms_recv(srcproc,com%seginfobuf(:,0,ibuf),3,tag=SEGINFO_TAG)
             call comms_recv(srcproc,com%idxbuf(:,ibuf),buflen,tag=BLKIDX_TAG)
          end if
       else
          seg_start = library(com%lib)%seg_info(s_idx,srcproc)
          com%seginfobuf(:,0,ibuf) = library(com%lib)%seg_info(:,srcproc)
          com%idxbuf(1:buflen,ibuf) = &
               library(com%lib)%blk_idx(seg_start:seg_start+buflen-1)
       end if

    ! If groupsegonly==.true., we are receiving only the segments of the
    ! index corresponding to the local group's procs
    else if (com%groupsegonly) then

       ! Recv group segments of all procs in group of srcproc
       buflen = library(com%lib)%idx_groupseg_lens(srcproc)
       if (buflen == 0) then
          com%seginfobuf(s_type,0,ibuf) = SEG_BLANK
          return
       end if
       if (pub_my_proc_id /= srcproc) then
          if (async) then
             call comms_irecv(srcproc,com%seginfobuf(:, &
                  0:pub_comms_group_size-1,ibuf),3*pub_comms_group_size, &
                  tag=SEGINFO_TAG,handle=com%handles(com%info_hdl))
             call comms_irecv(srcproc,com%idxbuf(:,ibuf),buflen, &
                  tag=BLKIDX_TAG,handle=com%handles(com%idx_hdl))
          else
             call comms_recv(srcproc,com%seginfobuf(:, &
                  0:pub_comms_group_size-1,ibuf),3*pub_comms_group_size, &
                  tag=SEGINFO_TAG)
             call comms_recv(srcproc,com%idxbuf(:,ibuf),buflen,tag=BLKIDX_TAG)
          end if
       else
          seg_start = library(com%lib)%seg_info(s_idx,pub_first_proc_in_group)
         com%seginfobuf(:,0:pub_comms_group_size-1,ibuf) = &
               library(com%lib)%seg_info(:, &
               pub_first_proc_in_group:pub_last_proc_in_group)
          com%idxbuf(1:buflen,ibuf) = &
               library(com%lib)%blk_idx(seg_start:seg_start+buflen-1)
       end if

    else ! else we are receiving the whole index

       buflen = library(com%lib)%idx_lens(srcproc)
       if (buflen == 0) then
          com%seginfobuf(s_type,:,ibuf) = SEG_BLANK
          return
       end if
       if (pub_my_proc_id /= srcproc) then
          if (async) then
             call comms_irecv(srcproc,com%seginfobuf(1,0,ibuf), &
                  3*pub_total_num_procs+3,tag=SEGINFO_TAG, &
                  handle=com%handles(com%info_hdl))
             call comms_irecv(srcproc,com%idxbuf(1,ibuf),buflen, &
                  tag=BLKIDX_TAG,handle=com%handles(com%idx_hdl))
          else
             call comms_recv(srcproc,com%seginfobuf(1,0,ibuf), &
                  3*pub_total_num_procs+3,tag=SEGINFO_TAG)
             call comms_recv(srcproc,com%idxbuf(1,ibuf),buflen,tag=BLKIDX_TAG)
          end if
       else
          com%seginfobuf(:,:,ibuf) = library(com%lib)%seg_info(:,:)
          com%idxbuf(1:buflen,ibuf) = library(com%lib)%blk_idx(1:buflen)
       end if

    end if

  end subroutine sparse_recv_index


  !==========================================================================!
  ! This subroutine sends either a segment of the pointers or the whole set  !
  ! of pointers for a given sparse matrix structure to another proc          !
  !--------------------------------------------------------------------------!
  ! Arguments:                                                               !
  !   lib                                                                    !
  !   destproc                                                               !
  !   ptr_handle                                                             !
  !--------------------------------------------------------------------------!
  ! Written by Nicholas Hine, May 2009.                                      !
  !==========================================================================!

  subroutine sparse_send_pointers(com,destproc,ptr_handle)

    use comms, only: comms_send, pub_my_proc_id, pub_comms_group_size

    implicit none

    ! Arguments
    type(COM3), intent(in) :: com     ! Matrix Communicator
    integer, intent(in) :: destproc   ! Destination proc to which to send
    integer, intent(inout), optional :: ptr_handle

    ! Locals
    integer :: buflen, seg_start
    integer :: destgroupfirstproc,destgrouplastproc

    if (pub_my_proc_id == destproc) return

    if (com%segonly) then
       seg_start = library(com%lib)%seg_info(s_idx,destproc)
       buflen = library(com%lib)%seg_info(s_idx,destproc+1) - seg_start
       call comms_send(destproc,library(com%lib)%blk_ptr(seg_start), &
            buflen,tag=BLKPTR_TAG)
    else if (com%groupsegonly) then
       ! Send pointers of all local segments corresponding to the
       ! procs of destproc's comms group
       destgroupfirstproc = &
            int(destproc/pub_comms_group_size)*pub_comms_group_size
       seg_start = library(com%lib)%seg_info(s_idx,destgroupfirstproc)
       destgrouplastproc = (int(destproc/pub_comms_group_size)+1) * &
            pub_comms_group_size - 1
       buflen = library(com%lib)%seg_info(s_idx,destgrouplastproc+1) - seg_start
       if (buflen < 1) return
       call comms_send(destproc,library(com%lib)%blk_ptr(seg_start), &
            buflen,tag=BLKPTR_TAG)
    else
       buflen = library(com%lib)%idx_lens(pub_my_proc_id)
       if (present(ptr_handle)) then
          call comms_send(destproc,library(com%lib)%blk_ptr,buflen, &
               tag=BLKPTR_TAG,return_handle=ptr_handle)
       else
          call comms_send(destproc,library(com%lib)%blk_ptr,buflen, &
               tag=BLKPTR_TAG)
       end if
    end if

  end subroutine sparse_send_pointers


  !==========================================================================!
  ! This subroutine receives either the whole set of pointers or a segment   !
  ! of the pointers for a given sparse matrix structure from another proc.   !
  !--------------------------------------------------------------------------!
  ! Arguments:                                                               !
  !   lib                                                                    !
  !   srcproc                                                                !
  !   ibuf                                                                   !
  !   async                                                                  !
  !--------------------------------------------------------------------------!
  ! Written by Nicholas Hine, May 2009.                                      !
  !==========================================================================!

  subroutine sparse_recv_pointers(com,srcproc,ibuf,async)

    use comms, only: comms_recv, comms_irecv, pub_my_proc_id, &
         pub_first_proc_in_group

    implicit none

    ! Arguments
    type(COM3), intent(inout) :: com      ! Sparse Matrix Communicator object
    integer, intent(in) :: srcproc        ! Source proc from which to receive
    integer, intent(in) :: ibuf           ! Which buffer to use
    logical, intent(in) :: async          ! Whether to send asynchronously

    ! Locals
    integer :: buflen
    integer :: seg_start

    ! If segonly==.true. we are receiving this proc's segment of the index only
    if (com%segonly) then
       buflen = library(com%lib)%idx_seg_lens(srcproc)
       seg_start = library(com%lib)%seg_info(s_idx,srcproc)
    else if (com%groupsegonly) then
       ! Recv group segments of all procs in group of srcproc
       buflen = library(com%lib)%idx_groupseg_lens(srcproc)
       if (pub_my_proc_id == srcproc) seg_start &
            = library(com%lib)%seg_info(s_idx,pub_first_proc_in_group)
    else ! else we are receiving the whole set of pointers
       buflen = library(com%lib)%idx_lens(srcproc)
       seg_start = 1
    end if

    if (pub_my_proc_id /= srcproc) then
       if (async) then
          call comms_irecv(srcproc,com%ptrbuf(:,ibuf),buflen,tag=BLKPTR_TAG, &
               handle=com%handles(com%ptr_hdl))
       else
          call comms_recv(srcproc,com%ptrbuf(:,ibuf),buflen,tag=BLKPTR_TAG)
       end if
    else
       com%ptrbuf(1:buflen,ibuf) = &
            library(com%lib)%blk_ptr(seg_start:seg_start+buflen-1)
    end if

  end subroutine sparse_recv_pointers


  !==========================================================================!
  ! This subroutine sends either the whole of the matrix data of a given     !
  ! proc or the segment of it corresponding to the local proc's rows         !
  !--------------------------------------------------------------------------!
  ! Arguments:                                                               !
  !   mat                                                                    !
  !   destproc                                                               !
  !--------------------------------------------------------------------------!
  ! Written by Nicholas Hine, May 2009.                                      !
  !==========================================================================!

  subroutine sparse_send_data(com,mat,destproc)

    use comms, only: comms_send, pub_my_proc_id, pub_total_num_procs, &
         pub_comms_group_size

    implicit none

    ! Arguments
    type(COM3), intent(in) :: com
    type(SPAM3),intent(in) :: mat
    integer, intent(in) :: destproc

    ! Locals
    integer :: datlen
    integer :: seg_start
    integer :: destgroupfirstproc,destgrouplastproc

    if (pub_my_proc_id == destproc) return

    ! If seg==.true., we are sending the segment of the data only
    if (com%segonly) then
       seg_start = library(mat%lib)%seg_info(s_ptr,destproc)
       datlen = library(mat%lib)%seg_info(s_ptr,destproc+1) - seg_start
    else if (com%groupsegonly) then
       ! Send data of all local segments corresponding to the
       ! procs of destproc's comms group
       destgroupfirstproc = &
            int(destproc/pub_comms_group_size)*pub_comms_group_size
       seg_start = library(com%lib)%seg_info(s_ptr,destgroupfirstproc)
       destgrouplastproc = (int(destproc/pub_comms_group_size)+1) * &
            pub_comms_group_size - 1
       datlen = library(com%lib)%seg_info(s_ptr,destgrouplastproc+1) - seg_start
    else
       seg_start = 1
       datlen = library(mat%lib)%seg_info(s_ptr,pub_total_num_procs) - 1
    end if

    if (mat%iscmplx) then
       call comms_send(destproc,mat%zmtx(seg_start:seg_start+datlen-1), &
            datlen,tag=DATA_TAG)
    else
       call comms_send(destproc,mat%dmtx(seg_start:seg_start+datlen-1), &
            datlen,tag=DATA_TAG)
    end if

  end subroutine sparse_send_data


  !==========================================================================!
  ! This subroutine receives either the whole set of data or a segment       !
  ! of the data for a given sparse matrix structure from another proc.       !
  !--------------------------------------------------------------------------!
  ! Arguments:                                                               !
  !   mat                                                                    !
  !   srcproc                                                                !
  !   dbuf                                                                   !
  !   zbuf                                                                   !
  !   totlen                                                                 !
  !--------------------------------------------------------------------------!
  ! Written by Nicholas Hine, May 2009.                                      !
  !==========================================================================!

  subroutine sparse_recv_data(com,mat,srcproc,ibuf,async)

    use comms, only: comms_recv, comms_irecv, pub_my_proc_id, &
         pub_total_num_procs,pub_first_proc_in_group

    implicit none

    ! Arguments
    type(COM3), intent(inout) :: com
    type(SPAM3),intent(in) :: mat
    integer, intent(in) :: srcproc
    integer, intent(in) :: ibuf
    logical, intent(in) :: async          ! Whether to send asynchronously

    ! Locals
    integer :: seg_start
    integer :: datlen

    if (com%segonly) then
       datlen = library(mat%lib)%mtx_seg_lens(srcproc)
    else if (com%groupsegonly) then
       ! Recv data of all procs in group of srcproc
       datlen = library(mat%lib)%mtx_groupseg_lens(srcproc)
    else
       datlen = com%seginfobuf(s_ptr,pub_total_num_procs,2)-1
    end if

    if (pub_my_proc_id /= srcproc) then
       if (mat%iscmplx) then
          if (async) then
             call comms_irecv(srcproc,com%zmtxrecvbuf(:,ibuf),datlen, &
                  tag=DATA_TAG,handle=com%handles(com%data_hdl))
          else
             call comms_recv(srcproc,com%zmtxrecvbuf(:,ibuf),datlen, &
                  tag=DATA_TAG)
          end if
       else
          if (async) then
             call comms_irecv(srcproc,com%dmtxrecvbuf(:,ibuf),datlen, &
                  tag=DATA_TAG,handle=com%handles(com%data_hdl))
          else
             call comms_recv(srcproc,com%dmtxrecvbuf(:,ibuf),datlen, &
                  tag=DATA_TAG)
          end if
       end if
    else
       if (com%segonly) then
          seg_start = library(mat%lib)%seg_info(s_ptr,srcproc)
       else if (com%groupsegonly) then
          seg_start = library(mat%lib)%seg_info(s_ptr,pub_first_proc_in_group)
       else
          seg_start = 1
       end if
       if (mat%iscmplx) then
          com%zmtxrecvbuf(1:datlen,ibuf) = &
               mat%zmtx(seg_start:seg_start+datlen-1)
       else
          com%dmtxrecvbuf(1:datlen,ibuf) = &
               mat%dmtx(seg_start:seg_start+datlen-1)
       end if
    end if

  end subroutine sparse_recv_data


  !============================================================================!
  ! This subroutine converts a SPAM3 matrix into a column-indexed,             !
  ! non-segmented form                                                         !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   mat      (input)  : The sparse matrix to be converted                    !
  !   idxinfo  (output) : The block-index in a non-segmented form              !
  !   dmtx     (output) : The actual non-zero matrix elements                  !
  !----------------------------------------------------------------------------!
  ! Written by Simon M.-M. Dubois (January 2010)                               !
  !    * largely inspired by the sparse_write_vector subroutine written        !
  !      by Nicholas Hine                                                      !
  ! Modified for non-square matrices by Robert Charlton, 12/09/2018.           !
  !============================================================================!

  subroutine sparse_convert_unsegment_real(mat, idxlen, dmtxlen, idxinfo, dmtx)

    use comms, only: comms_barrier, comms_bcast, comms_send, &
         comms_recv, comms_reduce
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_unit, &
         utils_assert

    implicit none

    ! Arguments
    type(SPAM3), intent(in)               :: mat(:)
    integer, intent(out)                  :: idxlen
    integer, intent(out)                  :: dmtxlen
    integer, intent(inout), optional      :: idxinfo(:)
    real(kind=DP), intent(inout), optional:: dmtx(:)

    ! Local variables
    integer :: ncomps               ! Number of components
    integer :: icomp                ! Loop counter for components
    integer :: row_blks, col_blks   ! Identifier for blocking scheme
    integer :: my_nblks             ! Number of blocks stored on this proc
    integer :: my_nzb               ! Number of non-zero blocks on this proc
    integer :: ierr                 ! Error flag
    integer :: ilib                 ! Library entry for mat
    integer :: iblk                 ! Atom loop counter
    integer :: iat_orig             ! Atom iblk in original order
    integer :: jblk                 ! Second atom
    integer :: jat_orig             ! Atom jblk in original order
    integer :: loc_iblk             ! Local atom counter for iat
    integer :: iidx                 ! Index loop counter
    integer :: jidx                 ! Index info counter
    integer :: ielems               ! Number of elements on atom iat
    integer, allocatable :: idxbuf(:)       ! Comm buffer for index
    ! rc2013: local pointers to parallel strategies
    type(PARAL_INFO), pointer :: row_par
    type(PARAL_INFO), pointer :: col_par

    ! Check arguments
    call utils_assert(.not. mat(1)%iscmplx, 'Error in &
            &sparse_convert_unsegment_real: real matrices only.')
    call utils_assert(present(idxinfo).eqv.present(dmtx), &
         'Error in sparse_convert_unsegment_real: &
         &only all/none optional arguments allowed.')

    ! Obtain number of components
    ncomps = size(mat)

    ! Obtain library entry for mat
    ilib = mat(1)%lib
    row_blks = library(ilib)%row_blks
    col_blks = library(ilib)%col_blks
    my_nblks = library(ilib)%my_nblks
    my_nzb = library(ilib)%my_nzb

    ! rc2013: assign parallel strategies
    row_par => library(ilib)%row_par
    col_par => library(ilib)%col_par

    ! Check that the matrix is square
    ! rc2013: leave this out for testing non-square structures
    !call utils_assert(blks == library(ilib)%col_blks, &
    !     'Error in sparse_convert_unsegment_real: &
    !     &row and column blocking schemes do not match.')

    ! rc2013: check that the row/col parallel strategies match
    call sparse_par_check(ilib, ilib, &
        routine='sparse_convert_unsegment_real', trace=.true.)

    ! Determine the space needed to store dmtx and idxinfo
    dmtxlen = library(ilib)%my_nze*ncomps
    idxlen  = my_nzb*3 + 2
    if (.not.present(dmtx)) return

    ! Check that dmtx and idxinfo have the convenient dimension
    call utils_assert(dmtxlen .eq. size(dmtx), &
       'Error in sparse_convert_unsegment_real: &
       &storage space for data incorrectly allocated.')

    call utils_assert(idxlen .eq. size(idxinfo), &
         'Error in sparse_convert_unsegment_real: &
         &storage space for index incorrectly allocated.')

    ! Get the index for this proc in non-segmented format
    allocate(idxbuf(my_nblks+my_nzb+2),stat=ierr)
    call utils_alloc_check('sparse_convert_unsegment_real','idxbuf',ierr)
    call sparse_generate_index(idxbuf,mat(1))

    ! Read segmented data into unsegmented buffers
    call internal_unsegment

    ! Write indexing information for the unsegment data into idxinfo :
    ! idxinfo(1)           > is set to idxlen
    ! idxinfo(2)           > is set to dmtxlen
    ! idxinfo(n*3+3:n*3+5) > contains the information regarding the nth
    !                        block of mat written in a column-indexed
    !                        non-segmented form

    idxinfo(1)=idxlen
    idxinfo(2)=dmtxlen


    jidx = 2
    loc_iblk = 0
    do iblk=col_par%my_first_blk,col_par%my_last_blk
       loc_iblk = loc_iblk + 1
       iat_orig = col_par%orig_atom(iblk)
       ielems = col_par%num_elems_on_atom(iblk,col_blks)

       ! Loop over block-rows in column iat
       do iidx=idxbuf(loc_iblk),idxbuf(loc_iblk+1)-1
          jblk = idxbuf(iidx)
          ! rc2013: row/col_par should match
          jat_orig = row_par%orig_atom(jblk)

          ! Generate index info: row, column and number of elements
          idxinfo(jidx+1) = jat_orig
          idxinfo(jidx+2) = iat_orig
          idxinfo(jidx+3) = ielems * row_par%num_elems_on_atom(jblk,row_blks)

          jidx = jidx + 3
       end do
    end do


    ! Re-sync procs
    call comms_barrier

    ! Deallocate workspace
    deallocate(idxbuf,stat=ierr)
    call utils_dealloc_check('sparse_convert_unsegment_real','idxbuf',ierr)

  contains

    subroutine internal_unsegment

      implicit none

      ! Locals
      integer :: iidx,ptr
      integer :: iblk,loc_iblk
      integer :: jblk
      integer :: jelems
      integer :: ielem
      integer :: max_elems_on_atom
      real (kind=DP), allocatable :: dblk(:,:)

      ! Allocate temporary array to store block
      max_elems_on_atom = maxval(col_par%num_elems_on_atom(:,col_blks))
      allocate(dblk(max_elems_on_atom,max_elems_on_atom),stat=ierr)
      call utils_alloc_check('internal_unsegment  &
           &(sparse_convert_unsegment_real)', 'dblk',ierr)

      ptr = 1

      ! Loop over components of the matrix array
      do icomp=1,ncomps

         loc_iblk = 0

         ! Loop over block-cols on this proc
         do iblk=col_par%my_first_blk,col_par%my_last_blk
            loc_iblk = loc_iblk + 1

            ! Loop over nonzero block-rows in this block-column
            do iidx=idxbuf(loc_iblk),idxbuf(loc_iblk+1)-1
               jblk = idxbuf(iidx)
               jelems = row_par%num_elems_on_atom(jblk,col_blks)

               ! Get this block and write it into buffer
               call sparse_get_block(dblk,mat(icomp),jblk,iblk)
               do ielem=1,col_par%num_elems_on_atom(iblk,col_blks)
                  dmtx(ptr:ptr+jelems-1) = dblk(1:jelems,ielem)
                  ptr = ptr + jelems
               end do

            end do

         end do

      end do

      ! Deallocate temporary array for block
      deallocate(dblk,stat=ierr)
      call utils_dealloc_check('internal_unsegment  &
           &(sparse_convert_unsegment_real)','dblk',ierr)

    end subroutine internal_unsegment

  end subroutine sparse_convert_unsegment_real


  !============================================================================!
  ! This subroutine converts a SPAM3 matrix into a column-indexed,             !
  ! non-segmented form                                                         !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   mat      (input)  : The sparse matrix to be converted                    !
  !   idxinfo  (output) : The block-index in a non-segmented form              !
  !   zmtx     (output) : The actual non-zero matrix elements                  !
  !----------------------------------------------------------------------------!
  ! Written by J.M. Escartin based on sparse_convert_unsegment_real.           !
  !============================================================================!

  subroutine sparse_convert_unsegment_complex(mat, idxlen, zmtxlen, &
       idxinfo, zmtx)

    use comms, only: comms_barrier, comms_bcast, comms_send, &
         comms_recv, comms_reduce
    use parallel_strategy, only: PARAL_INFO
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_unit, &
         utils_assert

    implicit none

    ! Arguments
    type(SPAM3), intent(in)               :: mat(:)
    integer, intent(out)                  :: idxlen
    integer, intent(out)                  :: zmtxlen
    integer, intent(out),optional         :: idxinfo(:)
    complex(kind=DP), intent(out),optional:: zmtx(:)

    ! Local variables
    integer :: ncomps               ! Number of components
    integer :: icomp                ! Loop counter for components
    integer :: row_blks, col_blks   ! Identifier for blocking scheme
    integer :: my_nblks             ! Number of blocks stored on this proc
    integer :: my_nzb               ! Number of non-zero blocks on this proc
    integer :: ierr                 ! Error flag
    integer :: ilib                 ! Library entry for mat
    integer :: iblk                 ! Atom loop counter
    integer :: iat_orig             ! Atom iblk in original order
    integer :: jblk                 ! Second atom
    integer :: jat_orig             ! Atom jblk in original order
    integer :: loc_iblk             ! Local atom counter for iat
    integer :: iidx                 ! Index loop counter
    integer :: jidx                 ! Index info counter
    integer :: ielems               ! Number of elements on atom iat
    integer, allocatable :: idxbuf(:)       ! Comm buffer for index
    character(len=*), parameter :: myself = 'sparse_convert_unsegment_complex'
    ! rc2013: use local pointers to parallel strategies for rows and columns
    type(PARAL_INFO), pointer :: row_par
    type(PARAL_INFO), pointer :: col_par


    ! Check arguments
    call utils_assert(mat(1)%iscmplx, 'Error in '//myself// &
            ': complex matrices only.')
    call utils_assert(present(idxinfo).eqv.present(zmtx), &
         'Error in sparse_convert_unsegment_complex: &
         &only all/none optional arguments allowed.')

    ! Obtain number of components
    ncomps = size(mat)

    ! Obtain library entry for mat
    ilib = mat(1)%lib
    row_blks = library(ilib)%row_blks
    col_blks = library(ilib)%col_blks
    my_nblks = library(ilib)%my_nblks
    my_nzb = library(ilib)%my_nzb

    ! rc2013: assign the parallel strategies from library structure
    row_par => library(ilib)%row_par
    col_par => library(ilib)%col_par

    ! Check that the matrix is square
    ! rc2013: leave this out for testing non-square structures
    !call utils_assert(blks == library(ilib)%col_blks, &
    !     'Error in '//myself//': row and column blocking schemes do not match.')

    ! Determine the space needed to store zmtx and idxinfo
    zmtxlen = library(ilib)%my_nze*ncomps
    idxlen  = my_nzb*3 + 2
    if (.not.present(zmtx)) return

    ! Check that zmtx and idxinfo have the convenient dimension
    call utils_assert(zmtxlen .eq. size(zmtx), &
         'Error in'//myself//': storage space for data incorrectly allocated.')

    call utils_assert(idxlen .eq. size(idxinfo), &
         'Error in'//myself//': storage space for index incorrectly allocated.')

    ! Get the index for this proc in non-segmented format
    allocate(idxbuf(my_nblks+my_nzb+2),stat=ierr)
    call utils_alloc_check(myself, 'idxbuf', ierr)
    call sparse_generate_index(idxbuf,mat(1))

    ! Read segmented data into unsegmented buffers
    call internal_unsegment

    ! Write indexing information for the unsegment data into idxinfo :
    ! idxinfo(1)           > is set to idxlen
    ! idxinfo(2)           > is set to zmtxlen
    ! idxinfo(n*3+3:n*3+5) > contains the information regarding the nth
    !                        block of mat written in a column-indexed
    !                        non-segmented form

    idxinfo(1)=idxlen
    idxinfo(2)=zmtxlen


    jidx = 2
    loc_iblk = 0
    do iblk=col_par%my_first_blk,col_par%my_last_blk
       loc_iblk = loc_iblk + 1
       iat_orig = col_par%orig_atom(iblk)
       ielems = col_par%num_elems_on_atom(iblk,col_blks)

       ! Loop over block-rows in column iat
       do iidx=idxbuf(loc_iblk),idxbuf(loc_iblk+1)-1
          jblk = idxbuf(iidx)
          jat_orig = row_par%orig_atom(jblk)

          ! Generate index info: row, column and number of elements
          idxinfo(jidx+1) = jat_orig
          idxinfo(jidx+2) = iat_orig
          idxinfo(jidx+3) = ielems * row_par%num_elems_on_atom(jblk,row_blks)

          jidx = jidx + 3
       end do
    end do


    ! Re-sync procs
    call comms_barrier

    ! Deallocate workspace
    deallocate(idxbuf,stat=ierr)
    call utils_dealloc_check(myself, 'idxbuf', ierr)

  contains

    subroutine internal_unsegment

      implicit none

      ! Locals
      integer :: iidx,ptr
      integer :: iblk,loc_iblk
      integer :: jblk
      integer :: jelems
      integer :: ielem
      integer :: max_elems_on_atom
      ! jme: bug (wrong type) fixed on 11/05/2017.
      complex (kind=DP), allocatable :: zblk(:,:)

      ! Allocate temporary array to store block
      max_elems_on_atom = maxval(col_par%num_elems_on_atom(:,col_blks))
      allocate(zblk(max_elems_on_atom,max_elems_on_atom),stat=ierr)
      call utils_alloc_check('internal_unsegment  &
           &(sparse_convert_unsegment_complex)', 'zblk',ierr)

      ptr = 1

      ! Loop over components of the matrix array
      do icomp=1,ncomps

         loc_iblk = 0

         ! Loop over block-cols on this proc
         do iblk=col_par%my_first_blk,col_par%my_last_blk
            loc_iblk = loc_iblk + 1

            ! Loop over nonzero block-rows in this block-column
            do iidx=idxbuf(loc_iblk),idxbuf(loc_iblk+1)-1
               jblk = idxbuf(iidx)
               jelems = row_par%num_elems_on_atom(jblk,row_blks)

               ! Get this block and write it into buffer
               call sparse_get_block(zblk,mat(icomp),jblk,iblk)
               do ielem=1,col_par%num_elems_on_atom(iblk,col_blks)
                  zmtx(ptr:ptr+jelems-1) = zblk(1:jelems,ielem)
                  ptr = ptr + jelems
               end do

            end do

         end do

      end do

      ! Deallocate temporary array for block
      deallocate(zblk,stat=ierr)
      call utils_dealloc_check('internal_unsegment  &
           &(sparse_convert_unsegment_complex)','zblk',ierr)

    end subroutine internal_unsegment

  end subroutine sparse_convert_unsegment_complex


  !============================================================================!
  ! This subroutine converts a block sparse matrix stored in a column-indexed, !
  ! non-segmented form into a SPAM3                                            !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   mat      (output)  : The sparse matrix to be converted                   !
  !   idxinfo  (input)   : The block-index in a non-segmented form             !
  !   dmtx     (input)   : The actual non-zero matrix elements                 !
  !----------------------------------------------------------------------------!
  ! Written by Simon M.-M. Dubois (January 2010)                               !
  !    * largely inspired by the sparse_read_vector subroutine written         !
  !      by Nicholas Hine                                                      !
  !============================================================================!

  subroutine sparse_convert_segment_real(mat,idxinfo,dmtx)

    use comms, only: comms_alltoall, comms_barrier, &
         comms_bcast, comms_free, comms_send, comms_recv, pub_my_proc_id, &
         pub_total_num_procs, comms_irecv, comms_probe, comms_test, &
         comms_waitany
    use parallel_strategy, only: PARAL_INFO
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_unit, &
         utils_flush, utils_assert

    implicit none

    ! Arguments
    type(SPAM3), intent(inout) :: mat(:)        ! Matrix to fill
    integer, intent(in)        :: idxinfo(:)
    real(kind=DP), intent(in)  :: dmtx(:)

    ! Local variables
    integer :: ncomps               ! Number of components
    integer :: icomp                ! Loop counter for components
    integer :: ierr                 ! Error flag
    integer :: ilib                 ! Library entry for mat
    integer :: row_blks, col_blks   ! Identifier for blocking scheme
    integer :: iblk                 ! Block loop counter
    integer :: iat_orig             ! Atom iblk in original order
    integer :: jblk                 ! Second block
    integer :: jat_orig             ! Atom jblk in original order
    integer :: iidx                 ! Index loop counter
    integer :: ielem,jelems         ! Element loop counter
    integer :: iptr                 ! Pointer
    integer :: idxlen               ! Index info length
    integer :: proc_s               ! Processor loop counters
                                    ! (sender proc in p2p commications)
    integer :: proc_r               ! Processor loop counters
                                    ! (receiver proc in p2p commications)
    integer :: idxptr               ! Index pointer
    integer :: datptr               ! Data pointer
    integer :: comp0                ! Offset between components in data
    integer :: blksize              ! Block size
    integer :: datlen               ! Buffer length for data
    integer :: max_elems_on_atom                ! Maximum number of elements
    integer :: send_num(2)
    integer :: receive_num(2)
    integer :: send_ptr(2)
    integer :: rcv_ptr(2)
    integer  :: ibuf_s(2*pub_total_num_procs)           ! Integer buffers
    integer  :: ibuf_r(2*pub_total_num_procs)           ! Integer buffers
    integer, allocatable :: idxbuf_s(:)         ! Comm buffers for index
    integer, allocatable :: idxbuf_r(:)         ! Comm buffers for index
    real(kind=DP), allocatable :: dmtxbuf_s(:)  ! Comm buffer for real data
    real(kind=DP), allocatable :: dmtxbuf_r(:)  ! Comm buffer for real data
    real (kind=DP), allocatable :: dblk(:,:)    ! Buffer for real block
    ! rc2013: local pointers for parallel strategies
    type(PARAL_INFO), pointer :: row_par
    type(PARAL_INFO), pointer :: col_par

    ! Obtain number of components
    ncomps = size(mat)
    icomp = 1

    ! Obtain library entry for mat
    ilib = mat(1)%lib
    row_blks = library(ilib)%row_blks
    col_blks = library(ilib)%row_blks

    ! rc2013: assign parallel strategies from library structures
    row_par => library(ilib)%row_par
    col_par => library(ilib)%col_par

    ! Check that the matrix is square
    ! rc2013: leave this out for testing non-square structures
    !call utils_assert(blks == library(ilib)%col_blks, &
    !     'Error in sparse_convert_segment_real: row and column blocking &
    !     &schemes do not match.')

    ! rc2013: check that the row/col parallel strategies match
    call sparse_par_check(ilib, ilib, &
        routine='sparse_convert_unsegment_real', trace=.true.)

    ! Global preliminaries
    ! Allocate temporary array to store block
    max_elems_on_atom = maxval(col_par%num_elems_on_atom(:,col_blks))
    allocate(dblk(max_elems_on_atom,max_elems_on_atom),stat=ierr)
    call utils_alloc_check('sparse_convert_segment_real', &
         'dblk',ierr)

    do icomp=1,ncomps
       mat(icomp)%dmtx = 0.0_DP
    end do

    idxlen = idxinfo(1)
    datlen = idxinfo(2)/ncomps

    ! Work out how much data must be sent to other procs
    ibuf_r(:) = 0
    ibuf_s(:) = 0
    do iidx=3,idxlen,3
       iat_orig = idxinfo(iidx+1)
       proc_r = col_par%proc_of_atom(iat_orig)
       if (proc_r .ne. pub_my_proc_id) then
          ibuf_s(2*proc_r+1) =  ibuf_s(2*proc_r+1) + 3
          ibuf_s(2*proc_r+2) =  ibuf_s(2*proc_r+2) +  idxinfo(iidx+2)
       endif
    enddo

    ibuf_r(:) = ibuf_s(:)             ! data to be sent to other procs
    call comms_alltoall(ibuf_r,2)     ! data to be received by other procs

    receive_num = 0
    send_num = 0
    do proc_r=0,pub_total_num_procs-1
       send_num(:) = send_num(:) + ibuf_s(2*proc_r+1:2*proc_r+2)
       receive_num(:) = receive_num(:) + ibuf_r(2*proc_r+1:2*proc_r+2)
    enddo

    ! Allocate buffers for data to be send
    if (send_num(1).ne.0) then
       allocate(idxbuf_s(send_num(1)),stat=ierr)
       call utils_alloc_check('sparse_convert_segment_real', &
            'idxbuf_s',ierr)
       allocate(dmtxbuf_s(send_num(2)*ncomps),stat=ierr)
       call utils_alloc_check('sparse_convert_segment_real', &
            'dmtxbuf_s',ierr)
       idxbuf_s = 0
       dmtxbuf_s = 0.0_DP
    endif

    ! Allocate buffers for data to be received
    if (receive_num(1).ne.0) then
       allocate(idxbuf_r(receive_num(1)),stat=ierr)
       call utils_alloc_check('sparse_convert_segment_real', &
            'idxbuf_r',ierr)
       allocate(dmtxbuf_r(receive_num(2)*ncomps),stat=ierr)
       call utils_alloc_check('sparse_convert_segment_real', &
            'dmtxbuf_r',ierr)
       idxbuf_r = 0
       dmtxbuf_r = 0.0_DP
    endif


    ! If required, fill the buffer arrays for outgoing communications
    if (send_num(1).ne.0) then

       idxptr = 1
       datptr = 1
       do proc_r=0,pub_total_num_procs-1

          if (ibuf_s(2*proc_r+1) == 0) then
             cycle
          endif

          iptr = 1
          do iidx=3,idxlen,3
             iat_orig = idxinfo(iidx+1)
             blksize = idxinfo(iidx+2)
             if (proc_r == col_par%proc_of_atom(iat_orig)) then
                idxbuf_s(idxptr) = idxinfo(iidx)
                idxbuf_s(idxptr+1) = iat_orig
                idxbuf_s(idxptr+2) = blksize
                comp0 = 0
                do icomp=1,ncomps
                   dmtxbuf_s(datptr:datptr+blksize-1) = &
                        dmtx(comp0+iptr:comp0+iptr+blksize-1)
                   comp0 = comp0 + datlen
                   datptr = datptr + blksize
                enddo
                idxptr = idxptr + 3
             endif
             iptr = iptr + blksize
          enddo

       enddo
    endif

    ! Synchronous loop over all procs to send/receive data
    rcv_ptr(:) = 1
    do proc_s=0,pub_total_num_procs-1

       if (pub_my_proc_id == proc_s) then
          send_ptr(:) = 1
          do proc_r=0,pub_total_num_procs-1
             if (ibuf_s(2*proc_r+1) == 0) then
                cycle
             endif
             call comms_send(proc_r,idxbuf_s(send_ptr(1):send_ptr(1)+ibuf_s(2*proc_r+1)-1), &
                        length=ibuf_s(2*proc_r+1),tag=1)
             call comms_send(proc_r,dmtxbuf_s(send_ptr(2):send_ptr(2)+ibuf_s(2*proc_r+2)*ncomps-1), &
                        length=ibuf_s(2*proc_r+2)*ncomps,tag=2)
             send_ptr(1) = send_ptr(1) + ibuf_s(2*proc_r+1)
             send_ptr(2) = send_ptr(2) + ibuf_s(2*proc_r+2)*ncomps
          enddo
       else
          if (ibuf_r(2*proc_s+1) .ne. 0) then
             call comms_recv(proc_s,idxbuf_r(rcv_ptr(1):rcv_ptr(1)+ibuf_r(2*proc_s+1)-1), &
                        length=ibuf_r(2*proc_s+1),tag=1)
             call comms_recv(proc_s,dmtxbuf_r(rcv_ptr(2):rcv_ptr(2)+ibuf_r(2*proc_s+2)*ncomps-1), &
                        length=ibuf_r(2*proc_s+2)*ncomps,tag=2)
             rcv_ptr(1) = rcv_ptr(1) + ibuf_r(2*proc_s+1)
             rcv_ptr(2) = rcv_ptr(2) + ibuf_r(2*proc_s+2)*ncomps
          endif
       endif

       call comms_barrier
    enddo

    ! Deposit data that are stored locally
    iptr = 1
    do iidx=3,idxlen,3
       iat_orig = idxinfo(iidx+1)
       blksize = idxinfo(iidx+2)
       ! rc2013: row/col schemes should match
       if (col_par%proc_of_atom(iat_orig) == pub_my_proc_id) then

          iblk = col_par%distr_atom(iat_orig)
          jat_orig = idxinfo(iidx)
          jblk = row_par%distr_atom(jat_orig)
          jelems = row_par%num_elems_on_atom(jblk,row_blks)
          comp0 = 0
          do icomp=1,ncomps
             do ielem=1,col_par%num_elems_on_atom(iblk,col_blks)
                dblk(1:jelems,ielem) = dmtx(comp0+iptr+&
                     (ielem-1)*jelems:comp0+iptr+ielem*jelems-1)
             end do
             call sparse_put_block(dblk,mat(icomp),jblk,iblk)
             comp0 = comp0 + datlen
          end do
       end if
       iptr = iptr + blksize
    end do

    ! If required, empty the buffer arrays for incoming communications
    if (receive_num(1).ne.0) then

       iptr = 1
       do iidx=1,receive_num(1),3
          jat_orig = idxbuf_r(iidx)
          jblk = col_par%distr_atom(jat_orig)
          iat_orig = idxbuf_r(iidx+1)
          iblk = row_par%distr_atom(iat_orig)
          blksize = idxbuf_r(iidx+2)
          jelems = col_par%num_elems_on_atom(jblk,col_blks)
          do icomp=1,ncomps
             do ielem=1,row_par%num_elems_on_atom(iblk,row_blks)
                dblk(1:jelems,ielem) = &
                     dmtxbuf_r(iptr+(ielem-1)*jelems:iptr+ielem*jelems-1)
             end do
             call sparse_put_block(dblk,mat(icomp),jblk,iblk)
             iptr = iptr + blksize
          end do
       end do

    endif

    call comms_free

    ! Deallocate communication buffers
    if (allocated(dmtxbuf_s)) then
       deallocate(dmtxbuf_s,stat=ierr)
       call utils_dealloc_check('sparse_convert_segment_real', &
            'dmtxbuf_s',ierr)
    end if
    if (allocated(dmtxbuf_r)) then
       deallocate(dmtxbuf_r,stat=ierr)
       call utils_dealloc_check('sparse_convert_segment_real', &
            'dmtxbuf_r',ierr)
    end if
    if (allocated(idxbuf_s)) then
       deallocate(idxbuf_s,stat=ierr)
       call utils_dealloc_check('sparse_convert_segment_real', &
            'idxbuf_s',ierr)
    end if
    if (allocated(idxbuf_r)) then
       deallocate(idxbuf_r,stat=ierr)
       call utils_dealloc_check('sparse_convert_segment_real', &
            'idxbuf_r',ierr)
    end if

    deallocate(dblk,stat=ierr)
    call utils_dealloc_check('sparse_convert_segment_real','dblk',ierr)

    ! Re-sync before continuing
    call comms_barrier

  end subroutine sparse_convert_segment_real


  !============================================================================!
  ! This subroutine converts a block sparse matrix stored in a column-indexed, !
  ! non-segmented form into a SPAM3                                            !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   mat      (output)  : The sparse matrix to be converted                   !
  !   idxinfo  (input)   : The block-index in a non-segmented form             !
  !   zmtx     (input)   : The actual non-zero matrix elements                 !
  !----------------------------------------------------------------------------!
  ! Adapted by J.M. Escartin from sparse_convert_segment_real                  !
  !============================================================================!

  subroutine sparse_convert_segment_complex(mat, idxinfo, zmtx)

    use constants, only: cmplx_0
    use comms, only: comms_alltoall, comms_barrier, &
         comms_bcast, comms_free, comms_send, comms_recv, pub_my_proc_id, &
         pub_total_num_procs, comms_irecv, comms_probe, comms_test, &
         comms_waitany
    use parallel_strategy, only: PARAL_INFO
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_unit, &
         utils_flush, utils_assert

    implicit none

    ! Arguments
    type(SPAM3), intent(inout) :: mat(:)        ! Matrix to fill
    integer, intent(in)        :: idxinfo(:)
    complex(kind=DP), intent(in)  :: zmtx(:)

    ! Local variables
    integer :: ncomps               ! Number of components
    integer :: icomp                ! Loop counter for components
    integer :: ierr                 ! Error flag
    integer :: ilib                 ! Library entry for mat
    integer :: row_blks, col_blks   ! Identifier for blocking scheme
    integer :: iblk                 ! Block loop counter
    integer :: iat_orig             ! Atom iblk in original order
    integer :: jblk                 ! Second block
    integer :: jat_orig             ! Atom jblk in original order
    integer :: iidx                 ! Index loop counter
    integer :: ielem,jelems         ! Element loop counter
    integer :: iptr                 ! Pointer
    integer :: idxlen               ! Index info length
    integer :: proc_s               ! Processor loop counters
                                    ! (sender proc in p2p commications)
    integer :: proc_r               ! Processor loop counters
                                    ! (receiver proc in p2p commications)
    integer :: idxptr               ! Index pointer
    integer :: datptr               ! Data pointer
    integer :: comp0                ! Offset between components in data
    integer :: blksize              ! Block size
    integer :: datlen               ! Buffer length for data
    integer :: max_elems_on_atom                ! Maximum number of elements
    integer :: send_num(2)
    integer :: receive_num(2)
    integer :: send_ptr(2)
    integer :: rcv_ptr(2)
    integer  :: ibuf_s(2*pub_total_num_procs)           ! Integer buffers
    integer  :: ibuf_r(2*pub_total_num_procs)           ! Integer buffers
    integer, allocatable :: idxbuf_s(:)         ! Comm buffers for index
    integer, allocatable :: idxbuf_r(:)         ! Comm buffers for index
    complex(kind=DP), allocatable :: zmtxbuf_s(:)  ! Comm buffer for complex data
    complex(kind=DP), allocatable :: zmtxbuf_r(:)  ! Comm buffer for complex data
    complex (kind=DP), allocatable :: zblk(:,:)    ! Buffer for complex block
    character(len=*), parameter :: myself = 'sparse_convert_segment_complex'
    type(PARAL_INFO), pointer :: row_par ! Parallel strategy for rows.
    type(PARAL_INFO), pointer :: col_par ! Parallel strategy for columns.

    ! Obtain number of components
    ncomps = size(mat)
    icomp = 1

    ! Obtain library entry for mat
    ilib = mat(1)%lib
    row_blks = library(ilib)%row_blks
    col_blks = library(ilib)%col_blks

    ! rc2013: assign parallel strategies using library structures.
    row_par => library(ilib)%row_par
    col_par => library(ilib)%col_par

    ! Check that the matrix is square
    ! rc2013: leave this out for testing non-square structures
    !call utils_assert(blks == library(ilib)%col_blks, &
    !     'Error in '//myself//': row and column blocking schemes do not match.')

    ! rc2013: check that the row/col parallel strategies match
    call sparse_par_check(ilib, ilib, &
        routine=myself, trace=.true.)

    ! Global preliminaries
    ! Allocate temporary array to store block
    max_elems_on_atom = maxval(col_par%num_elems_on_atom(:,col_blks))
    allocate(zblk(max_elems_on_atom,max_elems_on_atom),stat=ierr)
    call utils_alloc_check(myself, 'zblk', ierr)

    do icomp=1,ncomps
       mat(icomp)%zmtx = cmplx_0
    end do

    idxlen = idxinfo(1)
    datlen = idxinfo(2)/ncomps

    ! Work out how much data must be sent to other procs
    ibuf_r(:) = 0
    ibuf_s(:) = 0
    do iidx=3,idxlen,3
       iat_orig = idxinfo(iidx+1)
       proc_r = col_par%proc_of_atom(iat_orig)
       if (proc_r .ne. pub_my_proc_id) then
          ibuf_s(2*proc_r+1) =  ibuf_s(2*proc_r+1) + 3
          ibuf_s(2*proc_r+2) =  ibuf_s(2*proc_r+2) +  idxinfo(iidx+2)
       endif
    enddo

    ibuf_r(:) = ibuf_s(:)             ! data to be sent to other procs
    call comms_alltoall(ibuf_r,2)     ! data to be received by other procs

    receive_num = 0
    send_num = 0
    do proc_r=0,pub_total_num_procs-1
       send_num(:) = send_num(:) + ibuf_s(2*proc_r+1:2*proc_r+2)
       receive_num(:) = receive_num(:) + ibuf_r(2*proc_r+1:2*proc_r+2)
    enddo

    ! Allocate buffers for data to be send
    if (send_num(1).ne.0) then
       allocate(idxbuf_s(send_num(1)),stat=ierr)
       call utils_alloc_check(myself, 'idxbuf_s', ierr)
       allocate(zmtxbuf_s(send_num(2)*ncomps),stat=ierr)
       call utils_alloc_check(myself, 'zmtxbuf_s', ierr)
       idxbuf_s = 0
       zmtxbuf_s = cmplx_0
    endif

    ! Allocate buffers for data to be received
    if (receive_num(1).ne.0) then
       allocate(idxbuf_r(receive_num(1)),stat=ierr)
       call utils_alloc_check(myself, 'idxbuf_r', ierr)
       allocate(zmtxbuf_r(receive_num(2)*ncomps),stat=ierr)
       call utils_alloc_check(myself, 'zmtxbuf_r', ierr)
       idxbuf_r = 0
       zmtxbuf_r = cmplx_0
    endif


    ! If required, fill the buffer arrays for outgoing communications
    if (send_num(1).ne.0) then

       idxptr = 1
       datptr = 1
       do proc_r=0,pub_total_num_procs-1

          if (ibuf_s(2*proc_r+1) == 0) then
             cycle
          endif

          iptr = 1
          do iidx=3,idxlen,3
             iat_orig = idxinfo(iidx+1)
             blksize = idxinfo(iidx+2)
             if (proc_r == col_par%proc_of_atom(iat_orig)) then
                idxbuf_s(idxptr) = idxinfo(iidx)
                idxbuf_s(idxptr+1) = iat_orig
                idxbuf_s(idxptr+2) = blksize
                comp0 = 0
                do icomp=1,ncomps
                   zmtxbuf_s(datptr:datptr+blksize-1) = &
                        zmtx(comp0+iptr:comp0+iptr+blksize-1)
                   comp0 = comp0 + datlen
                   datptr = datptr + blksize
                enddo
                idxptr = idxptr + 3
             endif
             iptr = iptr + blksize
          enddo

       enddo
    endif

    ! Synchronous loop over all procs to send/receive data
    rcv_ptr(:) = 1
    do proc_s=0,pub_total_num_procs-1

       if (pub_my_proc_id == proc_s) then
          send_ptr(:) = 1
          do proc_r=0,pub_total_num_procs-1
             if (ibuf_s(2*proc_r+1) == 0) then
                cycle
             endif
             call comms_send(proc_r,idxbuf_s(send_ptr(1):send_ptr(1)+ibuf_s(2*proc_r+1)-1), &
                        length=ibuf_s(2*proc_r+1),tag=1)
             call comms_send(proc_r,zmtxbuf_s(send_ptr(2):send_ptr(2)+ibuf_s(2*proc_r+2)*ncomps-1), &
                        length=ibuf_s(2*proc_r+2)*ncomps,tag=2)
             send_ptr(1) = send_ptr(1) + ibuf_s(2*proc_r+1)
             send_ptr(2) = send_ptr(2) + ibuf_s(2*proc_r+2)*ncomps
          enddo
       else
          if (ibuf_r(2*proc_s+1) .ne. 0) then
             call comms_recv(proc_s,idxbuf_r(rcv_ptr(1):rcv_ptr(1)+ibuf_r(2*proc_s+1)-1), &
                        length=ibuf_r(2*proc_s+1),tag=1)
             call comms_recv(proc_s,zmtxbuf_r(rcv_ptr(2):rcv_ptr(2)+ibuf_r(2*proc_s+2)*ncomps-1), &
                        length=ibuf_r(2*proc_s+2)*ncomps,tag=2)
             rcv_ptr(1) = rcv_ptr(1) + ibuf_r(2*proc_s+1)
             rcv_ptr(2) = rcv_ptr(2) + ibuf_r(2*proc_s+2)*ncomps
          endif
       endif

       call comms_barrier
    enddo

    ! Deposit data that are stored locally
    iptr = 1
    do iidx=3,idxlen,3
       iat_orig = idxinfo(iidx+1)
       blksize = idxinfo(iidx+2)
       if (col_par%proc_of_atom(iat_orig) == pub_my_proc_id) then

          iblk = col_par%distr_atom(iat_orig)
          jat_orig = idxinfo(iidx)
          jblk = row_par%distr_atom(jat_orig)
          jelems = row_par%num_elems_on_atom(jblk,row_blks)
          comp0 = 0
          do icomp=1,ncomps
             do ielem=1,col_par%num_elems_on_atom(iblk,col_blks)
                zblk(1:jelems,ielem) = zmtx(comp0+iptr+&
                     (ielem-1)*jelems:comp0+iptr+ielem*jelems-1)
             end do
             call sparse_put_block(zblk,mat(icomp),jblk,iblk)
             comp0 = comp0 + datlen
          end do
       end if
       iptr = iptr + blksize
    end do

    ! If required, empty the buffer arrays for incoming communications
    if (receive_num(1).ne.0) then

       iptr = 1
       do iidx=1,receive_num(1),3
          jat_orig = idxbuf_r(iidx)
          jblk = col_par%distr_atom(jat_orig)
          iat_orig = idxbuf_r(iidx+1)
          iblk = row_par%distr_atom(iat_orig)
          blksize = idxbuf_r(iidx+2)
          jelems = col_par%num_elems_on_atom(jblk,col_blks)
          do icomp=1,ncomps
             do ielem=1,row_par%num_elems_on_atom(iblk,row_blks)
                zblk(1:jelems,ielem) = &
                     zmtxbuf_r(iptr+(ielem-1)*jelems:iptr+ielem*jelems-1)
             end do
             call sparse_put_block(zblk,mat(icomp),jblk,iblk)
             iptr = iptr + blksize
          end do
       end do

    endif

    call comms_free

    ! Deallocate communication buffers
    if (allocated(zmtxbuf_s)) then
       deallocate(zmtxbuf_s, stat=ierr)
       call utils_dealloc_check(myself, 'zmtxbuf_s', ierr)
    end if
    if (allocated(zmtxbuf_r)) then
       deallocate(zmtxbuf_r, stat=ierr)
       call utils_dealloc_check(myself, 'zmtxbuf_r', ierr)
    end if
    if (allocated(idxbuf_s)) then
       deallocate(idxbuf_s, stat=ierr)
       call utils_dealloc_check(myself, 'idxbuf_s', ierr)
    end if
    if (allocated(idxbuf_r)) then
       deallocate(idxbuf_r, stat=ierr)
       call utils_dealloc_check(myself, 'idxbuf_r', ierr)
    end if

    deallocate(zblk, stat=ierr)
    call utils_dealloc_check(myself, 'zblk', ierr)

    ! Re-sync before continuing
    call comms_barrier

  end subroutine sparse_convert_segment_complex


  !============================================================================!
  ! This subroutine returns the parallel strategy associated with a given      !
  ! sparse matrix. Can select row or column par if requested.                  !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   par      (output) : The parallel strategy pointer to return.             !
  !   mat      (input)  : The sparse matrix.                                   !
  !   rowcol   (input)  : Identifier for row/column parallel strategy.         !
  !----------------------------------------------------------------------------!
  ! Created by Robert Charlton, 21/02/2017.                                    !
  !============================================================================!

  subroutine sparse_get_par(par,mat,rowcol)

    use parallel_strategy, only: PARAL_INFO
    use utils, only: utils_abort

    implicit none

    ! Arguments
    type(SPAM3), intent(in)                  :: mat    ! The sparse matrix
    type(PARAL_INFO), pointer, intent(inout) :: par    ! Parallel strategy
    character(len=1), intent(in), optional   :: rowcol ! Row or column par?

    ! rc2013: nullify pointer if associated
    ! if(associated(par)) nullify(par)

    ! jd: ^ The above approach generates heaps of valgrind warnings that 'par'
    !     has been used without being defined. It is not OK to test an
    !     uninitialised pointer for association "Using an undefined pointer as
    !     an argument to the ASSOCIATED statement is an error because an
    !     undefined pointer cannot be referenced anywhere!" and "You can use
    !     NULLIFY to initialize an undefined pointer, giving it disassociated
    !     status. Then the pointer can be tested using the intrinsic function
    !     ASSOCIATED. I don't think there is any harm in nullifying a pointer
    !     that has not been associated, so I'll just drop the offending 'if'.
    nullify(par)

    ! rc2013: Point to row or column parallel strategy depending on input
    if(present(rowcol)) then
       if((rowcol=='R') .or. (rowcol=='r')) then
          par => library(mat%lib)%row_par
       else if((rowcol=='C') .or. (rowcol=='c')) then
          par => library(mat%lib)%col_par
       else
          call utils_abort('Error in sparse_get_par: &
          & invalid parallel strategy identifier.')
       endif
    else
       par => library(mat%lib)%row_par
    endif

  end subroutine sparse_get_par


  !============================================================================!
  ! This subroutine prints out the blocking information for an embedding       !
  ! matrix. Useful as diagnostic tool when combining matrices.                 !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   mat (inout) : The embedding sparse matrix.                               !
  !----------------------------------------------------------------------------!
  ! Written by Robert Charlton, 13/08/2017.                                    !
  !============================================================================!

  subroutine sparse_diagnose(mat)

    use comms, only: pub_on_root

    implicit none

    ! Arguments
    type(SPAM3), intent(in)      :: mat   ! The destination sparse matrix

    ! Local variables
    integer :: lib

    lib = mat%lib
    if(pub_on_root) then
        write(stdout,'(a12,a30)') 'structure = ', library(lib)%structure
        write(stdout,'(a18,i4)') '      nblk      = ', library(lib)%nblk
        write(stdout,'(a18,i4)') '      mblk      = ', library(lib)%mblk
    endif


  end subroutine sparse_diagnose

  !============================================================================!
  ! This function calculates the RMS element value of a 2D array of sparse     !
  ! matrices.                                                                  !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   mat   (input)  : The matrix to be assessed                               !
  !----------------------------------------------------------------------------!
  ! Written by Robert Charlton, 14/08/2017.                                    !
  ! Based on sparse_rms_element by Peter Haynes.                               !
  !============================================================================!

  real(kind=DP) function sparse_rms_element_array(mat)

    use comms, only: comms_reduce

    implicit none

    ! Arguments
    type(SPAM3), intent(in) :: mat(:,:)   ! The sparse matrices

    ! Local variables
    integer :: ilib            ! Library entry for mat
    integer :: iel             ! Element loop counter
    real(kind=DP) :: loc_sum   ! Sum accumulated on this processor
    real(kind=DP) :: glob_sum  ! Sum accumulated across all processors
    integer :: loc_nze         ! Total nze across all matrix blocks
    integer :: isub, jsub

    ! rc2013: Sum number of non-zero elements across all matrix blocks
    loc_nze = 0
    ! Sum modulus squared of each element on this processor
    loc_sum = 0.0_DP
    do jsub=1,size(mat,dim=1)
       do isub=1,size(mat,dim=1)
          ! Get library entry for mat
          ilib = mat(isub,jsub)%lib
          if (mat(isub,jsub)%iscmplx) then
             do iel=1,library(ilib)%my_nze
                 loc_sum = loc_sum + real(mat(isub,jsub)%zmtx(iel)* &
                     conjg(mat(isub,jsub)%zmtx(iel)),kind=DP)
             end do
          else
             do iel=1,library(ilib)%my_nze
                loc_sum = loc_sum + &
                     mat(isub,jsub)%dmtx(iel)*mat(isub,jsub)%dmtx(iel)
             end do
          end if
          loc_nze = loc_nze + library(ilib)%nze
       end do
    end do

    ! Sum up across processors
    glob_sum = loc_sum
    call comms_reduce('SUM',glob_sum)

    ! Calculate RMS value
    sparse_rms_element_array = sqrt(glob_sum / real(loc_nze,kind=DP))

  end function sparse_rms_element_array

  !============================================================================!
  ! This subroutine obtains the maximum eigenvalue of a given 2D array of      !
  ! sparse matrices by an iterative conjugate gradients procedure.             !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   mat (input)    : The matrix whose eigenvalue is desired                  !
  !   met (input)    : The metric to use                                       !
  !   eval (output)  : The eigenvalue estimate                                 !
  !   tol (input)    : Optional tolerance for eigenvalue estimate              !
  !   min_val (input): Optional flag to converge minimum eigenvalue instead    !
  !----------------------------------------------------------------------------!
  ! Written by Robert Charlton, 11/09/2017.                                    !
  ! Based on sparse_extremal_eigenvalue by Peter Haynes and Andrea Greco.      !
  !============================================================================!

  subroutine sparse_extremal_eigenvalue_array(mat,met,eval,tol,min_val,evec, &
       allow_warnings)

    ! agrecocmplx: add dependencies to use new types
    use comms, only: comms_bcast, comms_reduce, pub_on_root, &
         pub_my_proc_id, pub_root_proc_id, pub_total_num_procs
    use datatypes, only: FUNCTIONS, COEF, data_functions_dot, &
         data_functions_scale, data_functions_axpy, &
         data_set_to_zero, data_functions_copy
!$  use rundat, only: pub_threads_max
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_assert, &
         utils_isnan

    implicit none

    ! Arguments
    type(SPAM3), intent(in) :: mat(:,:) ! The matrix whose eigenvalue is desired
    type(SPAM3), intent(in) :: met(:,:) ! The metric to be used
    real(kind=DP), intent(out) :: eval  ! The eigenvalue
    real(kind=DP), optional, intent(in) :: tol  ! Optional tolerance
    logical, optional, intent(in) :: min_val
    ! agrecocmplx
    type(FUNCTIONS), optional, intent(inout) :: evec(:)
    logical, intent(in), optional :: allow_warnings

    ! Local variables
    integer, parameter :: maxiter = 200  ! Maximum number of iterations
    integer, parameter :: cgreset = 5    ! Conjugate gradients reset frequency
    integer, allocatable :: n(:)         ! Vector sub-block lengths
    integer :: iter                      ! Iteration counter
    integer :: conv                      ! Converged iteration counter
    logical :: local_min_val             ! local copy of min_val
    logical :: loc_warnings              ! local copy of allow_warnings
    real(kind=DP) :: local_tol           ! Local copy of tolerance
    real(kind=DP) :: norm                ! Normalisation
    real(kind=DP) :: normfac             ! Renormalisation factor
    real(kind=DP) :: gdotg               ! Norm of gradient
    real(kind=DP) :: gamma               ! Conjugate gradients gamma
    real(kind=DP) :: xdotd,ddotd,ddoth,ddoty    ! Dot products for line search
    real(kind=DP) :: a,b,c,disc,q,step,oldeval  ! Variables for line search
    ! rc2013: all these vectors have components for each subsystem
    ! agrecocmplx
    type(FUNCTIONS), allocatable :: xvec(:)  ! Trial eigenvector
    type(FUNCTIONS), allocatable :: yvec(:)  ! Product of matrix and search dirn
    type(FUNCTIONS), allocatable :: grad(:)  ! Gradient
    type(FUNCTIONS), allocatable :: dirn(:)  ! Search direction
    type(FUNCTIONS), allocatable :: matx(:)  ! Product of matrix and trial vector
    type(FUNCTIONS), allocatable :: metx(:)  ! Product of metric and trial vector
    !type(COEF) :: temp_coef
    logical :: loc_cmplx
    integer :: mrows, ncols, isub, jsub
    integer :: ierr

    ! rc2013: set size of matrix array
    mrows = size(mat,dim=1)
    ncols = size(mat,dim=2)

    ! This routine will only function for square matrices
    call utils_assert(mrows==ncols, &
         'Error in sparse_extremal_eigenvalue_array: square matrices only.')

    allocate(n(mrows),stat=ierr)
    call utils_alloc_check('sparse_extremal_eigenvalue_array', 'n', ierr)

    ! Deal with optional argument
    if (present(tol)) then
       local_tol = tol
    else
       local_tol = 0.01_DP
    end if
    if (present(allow_warnings)) then
       loc_warnings = allow_warnings
    else
       loc_warnings = .false.
    end if


    if (present(min_val)) then
       local_min_val = min_val
    else
       local_min_val = .false.
    endif

    ! Find vector length
    ! rc2013: nrows and mcols may differ (non-square off-diagonal matrices)
    do isub=1,mrows
        n(isub) = library(mat(1,isub)%lib)%mcols
    end do

    ! agrecocmplx
    ! rc2013: in principle matrices in array are either all real or all complex
    loc_cmplx = mat(1,1)%iscmplx

    ! Allocate workspace
    call internal_alloc_work

    ! Initialise eigenvector to random guess on one proc and broadcast
    ! ndmh_2012_11_01: workaround for BG/Q compilation issue.
    ! agrecocmplx: take care of complex case as well
    if (pub_on_root) then
        ! rc2013: randomise each vector sub-block
        do isub=1,mrows
            do iter=1,n(isub)
                call random_number(norm)
                if (loc_cmplx) then
                    !call random_number(gamma)
                    !xvec%z(iter)=cmplx(norm,gamma,kind=DP)
                    xvec(isub)%z(iter)=cmplx(norm,0.0_DP,kind=DP)
                else
                    xvec(isub)%d(iter)=norm
                end if
            end do
        end do
    end if

    ! agrecocmplx
    do isub=1,mrows
        if (loc_cmplx) then
            call comms_bcast(pub_root_proc_id,xvec(isub)%z)
        else
            call comms_bcast(pub_root_proc_id,xvec(isub)%d)
        end if
    end do

    ! Normalise
    ! agrecocmplx
    call internal_matvec(metx,met,xvec)
    norm = internal_get_coeff(xvec, metx)
    normfac = 1.0_DP / sqrt(norm)

    do isub=1,mrows
       ! agrecocmplx
       call data_functions_scale(xvec(isub),normfac)
       call data_functions_scale(metx(isub),normfac)

       ! agrecocmplx
       call data_set_to_zero(dirn(isub))
    end do

    ! Initialise CG
    gdotg = 1.0_DP
    gamma = 0.0_DP

    ! Obtain initial estimate of eigenvalue and matx := mat . x
    call internal_matvec(matx,mat,xvec)
    eval = internal_get_coeff(metx, matx)

    ! ndmh: suggested by Victor Milman to avoid possible WIN32 parallel desync
    call comms_bcast(pub_root_proc_id,eval)

    ! Start of conjugate gradients loop
    conv = 0
    do iter=1,maxiter

       ! Obtain gradient grad
       ! agrecocmplx
       do isub=1,mrows
          call data_functions_copy(grad(isub),matx(isub))
          call data_functions_axpy(grad(isub),xvec(isub),-eval)

          ! ndmh: avoid parallel desynchronisation of grad, which
          ! ndmh: really should not happen but does on DARWIN (??)
          ! agrecocmplx
          if (loc_cmplx) then
              call comms_bcast(pub_root_proc_id,grad(isub)%z)
          else
              call comms_bcast(pub_root_proc_id,grad(isub)%d)
          end if
       end do

       ! Obtain conjugate direction dirn
       gamma = 1.0_DP / gdotg
       call internal_matvec(yvec,met,grad)
       ! agrecocmplx
       gdotg = internal_get_coeff(grad, yvec)

       if (mod(iter,cgreset) == 0) then
          gamma = 0.0_DP
       else
          gamma = gamma * gdotg
       end if

       ! agrecocmplx
       do isub=1,mrows
          call data_functions_scale(dirn(isub),gamma)
          call data_functions_axpy(dirn(isub),grad(isub),-1.0_DP)
       end do

       ! Orthogonalise to xvec
       ! agrecocmplx
       xdotd = internal_get_coeff(metx, dirn)

       ! agrecocmplx
       do isub=1,mrows
          call data_functions_axpy(dirn(isub),xvec(isub),-xdotd)
       end do

       ! Obtain yvec := mat . dirn
       call internal_matvec(yvec,mat,dirn)

       ! Obtain metx := met .dirn
       call internal_matvec(metx,met,dirn)

       ! Obtain line search coefficients
       ! agrecocmplx
       ddotd = internal_get_coeff(metx, dirn)
       ddoth = internal_get_coeff(metx, matx)
       ddoty = internal_get_coeff(metx, yvec)

       a = ddotd * ddoth
       call comms_bcast(pub_root_proc_id,a)
       if (abs(a) < epsilon(1.0_DP)) exit

       b = eval * ddotd - ddoty
       call comms_bcast(pub_root_proc_id,b)
       if (abs(b) < epsilon(1.0_DP)) exit

       c = -ddoth
       disc = b*b - 4.0_DP*a*c
       call comms_bcast(pub_root_proc_id,c)
       call comms_bcast(pub_root_proc_id,disc)
       if (disc < 0.0_DP) exit

       !tjz21: Check if min or max eval is requested
       if(local_min_val) then
          q = -0.5_DP * (b - sign(sqrt(disc),b))
       else
          q = -0.5_DP * (b + sign(sqrt(disc),b))
       endif

       if (b < 0.0_DP) then
          step = q / a
       else
          step = c / q
       end if

       do isub=1,mrows
          ! Update vector xvec and matx
          ! agrecocmplx
          call data_functions_axpy(xvec(isub),dirn(isub),step)
          call data_functions_axpy(matx(isub),yvec(isub),step)

          ! ndmh: avoid parallel desynchronisation of xvec and matx
          ! agrecocmplx
          if (loc_cmplx) then
              call comms_bcast(pub_root_proc_id,xvec(isub)%z)
              call comms_bcast(pub_root_proc_id,matx(isub)%z)
          else
              call comms_bcast(pub_root_proc_id,xvec(isub)%d)
              call comms_bcast(pub_root_proc_id,matx(isub)%d)
          end if
       end do

       ! Renormalise
       call internal_matvec(metx,met,xvec)
       ! agrecocmplx
       norm = internal_get_coeff(xvec, metx)
       normfac = 1.0_DP / sqrt(norm)

       call comms_bcast(pub_root_proc_id,normfac)

       ! agrecocmplx
       do isub=1,mrows
          call data_functions_scale(xvec(isub),normfac)
          call data_functions_scale(matx(isub),normfac)
          call data_functions_scale(metx(isub),normfac)
       end do

       ! Re-estimate eigenvalue
       oldeval = eval
       eval = internal_get_coeff(metx, matx)

       ! ndmh: prevent WIN32 desync as above
       call comms_bcast(pub_root_proc_id,eval)

       ! Check convergence (to better than tol)
       if (abs(eval - oldeval) < local_tol) then
          conv = conv + 1
       else
          conv = 0
       end if
       ! ndmh: prevent WIN32 desync as above
       call comms_bcast(pub_root_proc_id,conv)
       if (conv > cgreset .or. abs(gdotg) < epsilon(1.0_DP)) exit

    end do

    if (loc_warnings.and.((iter==maxiter).or.(utils_isnan(eval)))) then
       if (pub_on_root) then
          write(stdout,'(a)') 'WARNING in sparse_extremal_eigenvalue_array: &
               &convergence did not reach set threshold:'
          write(stdout,'(a,i5,3f20.14)') 'iter = ',iter,' eval = ',eval, &
               ' oldeval = ',oldeval,' step = ',step,' gdotg = ',gdotg
       end if
    end if

    if(present(evec)) then
       do isub=1,mrows
          call data_functions_copy(evec(isub),xvec(isub))
          if (loc_cmplx) then
              call comms_bcast(pub_root_proc_id,evec(isub)%z)
          else
              call comms_bcast(pub_root_proc_id,evec(isub)%d)
          end if
       end do
    endif

    deallocate(n,stat=ierr)
    call utils_dealloc_check('sparse_extremal_eigenvalue_array', 'n', ierr)

    ! Deallocate workspace
    call internal_dealloc_work


  contains

    !==========================================================================!
    ! This subroutine performs a matrix-vector product between a block-sparse  !
    ! matrix and a single dense vector.                                        !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   yvec (output) : the result vector                                      !
    !   amat (input)  : the block sparse matrix                                !
    !   xvec (input)  : the original vector                                    !
    !--------------------------------------------------------------------------!
    ! Written by Peter Haynes, March 2004.                                     !
    ! Revised for distributed data, June 2006.                                 !
    ! Modified for SPAM3 by Nicholas Hine, June 2009.                          !
    ! Modified by Andrea Greco to allow compatibility with complex case,       !
    ! February 2016.                                                           !
    ! Modified for multiple subsystems by Robert Charlton, 11/09/2017.         !
    !==========================================================================!

    subroutine internal_matvec(yvec,amat,xvec)

      use datatypes, only: data_functions_alloc, data_functions_dealloc
      use parallel_strategy, only: PARAL_INFO

      implicit none

      ! Arguments
      ! agrecocmplx
      type(FUNCTIONS), intent(inout) :: yvec(:)  ! result vector
      type(SPAM3), intent(in) :: amat(:,:)       ! block sparse matrix
      type(FUNCTIONS), intent(in) :: xvec(:)     ! input vector

      ! Local variables
      integer :: alib        ! Library entry for amat
      integer :: row_blks    ! Identifier for row blocking scheme
      integer :: col_blks    ! Identifier for row blocking scheme
      integer :: iblk        ! Block-column counter
      integer :: loc_iblk    ! Block-column counter local to this proc
      integer :: jblk        ! Block-row counter
      integer :: idx         ! Index
      integer :: ielems      ! Number of elems in block-column
      integer :: jelems      ! Number of elems in block-row
      integer :: ielem       ! Row elem counter for input vector
      integer :: jelem       ! Row elem counter for result vector
      integer :: seg         ! Segment index
      integer :: seg_start   ! Start position of this segment in the index
      type(PARAL_INFO), pointer :: col_par ! Parallel strategy for columns.
      type(PARAL_INFO), pointer :: row_par ! Parallel strategy for rows.
      type(FUNCTIONS), allocatable :: tmp_y(:,:) ! rc2013: placeholder for matrix multiplication
      type(FUNCTIONS) :: tmp_y2 ! rc2013: placeholder for matrix multiplication
      type(COEF)    :: temp, tmp_dot
      real(kind=DP) :: tmp_sum, total_sum
      integer :: kk

      allocate(tmp_y(mrows,ncols),stat=ierr)
      call utils_alloc_check('internal_matvec','tmp_y',ierr)

      do jsub=1,mrows
         ! Set result vector to zero
         ! agrecocmplx
         call data_set_to_zero(yvec(jsub))

         ! rc2013: allocate data functions
         call data_functions_alloc(tmp_y2,n(jsub),iscmplx=loc_cmplx)

         do isub=1,ncols

            call data_functions_alloc(tmp_y(jsub,isub),n(jsub),iscmplx=loc_cmplx)
            ! Get library entry
            alib = amat(jsub,isub)%lib
            row_blks = library(alib)%row_blks
            col_blks = library(alib)%col_blks

            ! rc2013: assign parallel strategy from library structure
            col_par => library(alib)%col_par
            row_par => library(alib)%row_par

            ! The product to be formed is:
            !  y  = A   x
            !   j    ji  i
            call data_set_to_zero(tmp_y(jsub,isub))
            call data_set_to_zero(tmp_y2)

            ! Loop over segments of the matrix on this proc
!$OMP PARALLEL DO NUM_THREADS(pub_threads_max) DEFAULT (none) &
!$OMP PRIVATE(seg,ielems,ielem,jelems,jelem,seg_start,loc_iblk,iblk, &
!$OMP      idx,jblk) &
!$OMP SHARED(alib,amat,library,pub_total_num_procs,xvec,yvec,tmp_y, &
!$OMP      row_par,col_par,pub_my_proc_id,row_blks,col_blks, &
!$OMP      pub_threads_max,isub,jsub)
            do seg=0,pub_total_num_procs-1

               if (library(alib)%seg_info(s_type,seg)==SEG_DENSE) then

                  ! Do matrix-vector multiplication
                  ielems = col_par%num_elems_on_proc(pub_my_proc_id,col_blks)
                  ielem  = col_par%first_elem_on_proc(pub_my_proc_id,col_blks)
                  jelems = row_par%num_elems_on_proc(seg,row_blks)
                  jelem  = row_par%first_elem_on_proc(seg,row_blks)

                  ! agrecocmplx: complex case
                  if (amat(jsub,isub)%iscmplx) then
                     call zgemv('N',jelems,ielems,(1.0_DP,0.0_DP), &
                          amat(jsub,isub)%zmtx(library(alib)%seg_info(s_ptr, &
                          seg)), jelems,xvec(isub)%z(ielem),1,(1.0_DP,0.0_DP), &
                          tmp_y(jsub,isub)%z(jelem),1)
                  ! real case
                  else
                     call dgemv('N',jelems,ielems,1.0_DP, &
                          amat(jsub,isub)%dmtx(library(alib)%seg_info(s_ptr, &
                          seg)), jelems,xvec(isub)%d(ielem),1,1.0_DP, &
                          tmp_y(jsub,isub)%d(jelem),1)
                  end if

               else if (library(alib)%seg_info(s_type,seg)==SEG_SPARSE) then

                  ! Loop over block-columns of matrix and input vector on this proc
                  seg_start = library(alib)%seg_info(s_idx,seg)
                  loc_iblk = 0
                  do iblk=col_par%my_first_blk,col_par%my_last_blk
                     loc_iblk = loc_iblk + 1

                     ! Number of elems and first elem in this block-column
                     ielems = col_par%num_elems_on_atom(iblk,col_blks)
                     ielem = col_par%first_elem_on_atom(iblk,col_blks)

                     ! Loop over block-rows jblk in block-col iblk of matrix
                     do idx=library(alib)%blk_idx(seg_start+loc_iblk-1), &
                          library(alib)%blk_idx(seg_start+loc_iblk)-1

                        jblk = library(alib)%blk_idx(idx)

                        ! Number of elems and first elem in this block-row
                        jelems = row_par%num_elems_on_atom(jblk,row_blks)
                        jelem = row_par%first_elem_on_atom(jblk,row_blks)

                        ! Do block multiplication
                        ! agrecocmplx: complex case
                        if (amat(jsub,isub)%iscmplx) then
                           call zgemv('N',jelems,ielems,(1.0_DP,0.0_DP), &
                                amat(jsub,isub)%zmtx(library(alib)%blk_ptr(idx)), &
                                jelems,xvec(isub)%z(ielem),1, &
                                cmplx(1.0_DP,0.0_DP,kind=DP), &
                                tmp_y(jsub,isub)%z(jelem),1)
                        ! real case
                        else
                           call dgemv('N',jelems,ielems,1.0_DP, &
                                amat(jsub,isub)%dmtx(library(alib)%blk_ptr(idx)), &
                                jelems,xvec(isub)%d(ielem),1,1.0_DP, &
                                tmp_y(jsub,isub)%d(jelem),1)
                        end if

                     end do  ! Loop over block-rows jblk in block-col iblk of matrix

                  end do  ! Loop over block-columns of matrix and input vector

               end if  ! seg_info(s_type,seg)
            end do  ! seg
!$OMP END PARALLEL DO
            call data_functions_axpy(yvec(jsub), tmp_y(jsub,isub), 1.0_DP)
         end do ! cols loop
         ! rc2013: free up memory
         call data_functions_dealloc(tmp_y2)
         ! Sum contributions from all procs
         ! agrecocmplx
         if (loc_cmplx) then
            call comms_reduce('SUM',yvec(jsub)%z)
         else
            call comms_reduce('SUM',yvec(jsub)%d)
         end if

      end do ! rows loop

      ! rc2013: deallocate
      do jsub=1,ncols
         do isub=1,mrows
            call data_functions_dealloc(tmp_y(isub,jsub))
         end do
      end do

      deallocate(tmp_y,stat=ierr)
      call utils_dealloc_check('internal_matvec','tmp_y',ierr)

    end subroutine internal_matvec


    !==========================================================================!
    ! This subroutine allocates workspace for the parent subroutine.           !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   None                                                                   !
    !--------------------------------------------------------------------------!
    ! Written by Peter Haynes, March 2004.                                     !
    ! Revised for distributed data, June 2006.                                 !
    ! Modified by Andrea Greco to allow compatibility with complex case,       !
    ! February 2016.                                                           !
    ! Modified for multiple subsystems by Robert Charlton, 11/09/2017.         !
    !==========================================================================!

    subroutine internal_alloc_work

      ! agrecocmplx
      use datatypes, only: data_functions_alloc

      implicit none

      ! rc2013: allocate memory for all arrays
      allocate(xvec(mrows),stat=ierr)
      call utils_alloc_check('sparse_extremal_eigenvalue_array','xvec',ierr)
      allocate(yvec(mrows),stat=ierr)
      call utils_alloc_check('sparse_extremal_eigenvalue_array','yvec',ierr)
      allocate(grad(mrows),stat=ierr)
      call utils_alloc_check('sparse_extremal_eigenvalue_array','grad',ierr)
      allocate(dirn(mrows),stat=ierr)
      call utils_alloc_check('sparse_extremal_eigenvalue_array','dirn',ierr)
      allocate(matx(mrows),stat=ierr)
      call utils_alloc_check('sparse_extremal_eigenvalue_array','matx',ierr)
      allocate(metx(mrows),stat=ierr)
      call utils_alloc_check('sparse_extremal_eigenvalue_array','metx',ierr)

      ! agrecocmplx
      ! Allocate workspace using appropriate functions
      do isub=1,mrows
         call data_functions_alloc(xvec(isub),n(isub),iscmplx=loc_cmplx)
         call data_functions_alloc(yvec(isub),n(isub),iscmplx=loc_cmplx)
         call data_functions_alloc(grad(isub),n(isub),iscmplx=loc_cmplx)
         call data_functions_alloc(dirn(isub),n(isub),iscmplx=loc_cmplx)
         call data_functions_alloc(matx(isub),n(isub),iscmplx=loc_cmplx)
         call data_functions_alloc(metx(isub),n(isub),iscmplx=loc_cmplx)
      end do

    end subroutine internal_alloc_work


    !==========================================================================!
    ! This subroutine deallocates workspace for the parent subroutine.         !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   None                                                                   !
    !--------------------------------------------------------------------------!
    ! Written by Peter Haynes, March 2004.                                     !
    ! Revised for distributed data, June 2006.                                 !
    ! Modified by Andrea Greco to allow compatibility with complex case,       !
    ! February 2016.                                                           !
    ! Modified for multiple subsystems by Robert Charlton, 11/09/2017.         !
    !==========================================================================!

    subroutine internal_dealloc_work

      ! agrecocmplx
      use datatypes, only: data_functions_dealloc

      implicit none

      ! agrecocmplx
      ! Deallocate workspace using appropriate routines
      do isub=1,mrows
         call data_functions_dealloc(metx(isub))
         call data_functions_dealloc(matx(isub))
         call data_functions_dealloc(dirn(isub))
         call data_functions_dealloc(grad(isub))
         call data_functions_dealloc(yvec(isub))
         call data_functions_dealloc(xvec(isub))
      end do

      deallocate(metx,stat=ierr)
      call utils_dealloc_check('sparse_extremal_eigenvalue_array','metx',ierr)
      deallocate(matx,stat=ierr)
      call utils_dealloc_check('sparse_extremal_eigenvalue_array','matx',ierr)
      deallocate(dirn,stat=ierr)
      call utils_dealloc_check('sparse_extremal_eigenvalue_array','dirn',ierr)
      deallocate(grad,stat=ierr)
      call utils_dealloc_check('sparse_extremal_eigenvalue_array','grad',ierr)
      deallocate(yvec,stat=ierr)
      call utils_dealloc_check('sparse_extremal_eigenvalue_array','yvec',ierr)
      deallocate(xvec,stat=ierr)
      call utils_dealloc_check('sparse_extremal_eigenvalue_array','xvec',ierr)

    end subroutine internal_dealloc_work

    !==========================================================================!
    ! This function calculates the dot product of two vectors using real or    !
    ! complex data functions.                                                  !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   xvec (input)  ::  Row vector.                                          !
    !   yvec (input)  ::  Column vector.                                       !
    !--------------------------------------------------------------------------!
    ! Written by Robert Charlton, 11/09/2017.                                  !
    !==========================================================================!

    real(kind=DP) function internal_get_coeff(xvec, yvec) result(dotp)

      use datatypes, only: FUNCTIONS, COEF, data_functions_dot

      implicit none

      ! Arguments
      type(FUNCTIONS), intent(in)  :: xvec(:)     ! input vector
      type(FUNCTIONS), intent(in)  :: yvec(:)     ! input vector

      ! Local variables
      type(COEF)    :: temp_coef

      dotp = 0.0_DP
      do isub=1,mrows
         temp_coef = data_functions_dot(xvec(isub),yvec(isub))
         if (temp_coef%iscmplx) then
            dotp = dotp + real(temp_coef%z,kind=DP)
         else
            dotp = dotp + temp_coef%d
         endif
      end do

    end function internal_get_coeff

  end subroutine sparse_extremal_eigenvalue_array

  !============================================================================!
  ! This subroutine performs consistency checks of the parallel strategies of  !
  ! 2 (e.g. sparse_copy) or 3 (e.g. sparse_product) matrices.                  !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   clib  (inout)   : The library index for sparse matrix C.                 !
  !   alib  (input)   : The library index for sparse matrix A.                 !
  !   blib  (input)   : The library index for (optional) sparse matrix B.      !
  !   routine (input) : Name of routine calling us (for debugging).            !
  !   trace  (input)  : Treat this as trace of a matrix product C_ij*A_ji.     !
  !----------------------------------------------------------------------------!
  ! Written by Robert Charlton, 23/09/2017.                                    !
  !============================================================================!

  subroutine sparse_par_check(clib, alib, blib, routine, trace)

    use parallel_strategy, only: PARAL_INFO
    use utils, only: utils_assert

    implicit none

    ! Arguments
    integer, intent(in)           :: clib    ! The sparse matrix C library index
    integer, intent(in)           :: alib    ! The sparse matrix A library index
    integer, intent(in), optional :: blib    ! The sparse matrix B library index
    character(len=*), intent(in)  :: routine ! Identifier for the routine calling us
    logical, intent(in), optional :: trace    ! Treat C and A as trace of a product?

    ! Local variables
    ! rc2013: local pointers to row and column parallel strategies
    type(PARAL_INFO), pointer :: arow_par, acol_par
    type(PARAL_INFO), pointer :: brow_par, bcol_par
    type(PARAL_INFO), pointer :: crow_par, ccol_par

    ! jd: Sanity check against garbage lib indices. This helps detects errors
    !     (see e.g. #1895) that only manifest in release mode and so the
    !     RTL's bound-check doesn't help.

    call utils_assert(alib >= lbound(library,1) .and. &
         alib <= ubound(library,1), 'sparse_par_check: alib out of range', alib)
    call utils_assert(clib >= lbound(library,1) .and. &
         clib <= ubound(library,1), 'sparse_par_check: clib out of range', clib)

    ! rc2013: assign the pointers for par based on library structure
    arow_par => library(alib)%row_par
    crow_par => library(clib)%row_par
    acol_par => library(alib)%col_par
    ccol_par => library(clib)%col_par
    if(present(blib)) then
       call utils_assert(blib >= lbound(library,1) .and. &
            blib <= ubound(library,1), 'sparse_par_check: blib out of range', &
            blib)
       brow_par => library(blib)%row_par
       bcol_par => library(blib)%col_par
    endif

    ! rc2013: check consistency based on arguments passed
    if(present(blib)) then
       ! rc2013: check that the parallel strategies match
       call utils_assert((acol_par%par_index==brow_par%par_index), &
            'Error in sparse_par_check: acol and brow parallel strategies &
            &do not match in '//trim(routine))
       call utils_assert((ccol_par%par_index==bcol_par%par_index), &
            'Error in sparse_par_check: ccol and bcol parallel strategies &
            &do not match in '//trim(routine))
    else
       if(present(trace)) then
          ! rc2013: when calculating trace of product, need C_{ij}*A_{ji}
          call utils_assert((ccol_par%par_index==arow_par%par_index), &
               'Error in sparse_par_check: ccol and arow parallel strategies &
               &do not match (trace) in '//trim(routine))
          call utils_assert((crow_par%par_index==acol_par%par_index), &
               'Error in sparse_par_check: crow and acol parallel strategies &
               &do not match (trace) in '//trim(routine))
       else
          call utils_assert((crow_par%par_index==arow_par%par_index), &
               'Error in sparse_par_check: crow and arow parallel strategies &
               &do not match in '//trim(routine))
          call utils_assert((ccol_par%par_index==acol_par%par_index), &
               'Error in sparse_par_check: ccol and acol parallel strategies &
               &do not match in '//trim(routine))
       endif
    endif

  end subroutine sparse_par_check

  !============================================================================!
  ! Returns the index of the first atom in a given comms group of procs.       !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !    par        (in) : Parallel strategy that we wish to probe.              !
  !    first_proc (in) : The first proc in this comms group.                   !
  !    last_proc  (in) : The last proc in this comms group.                    !
  !    group_size (in) : The (optional) size of this group.                    !
  !----------------------------------------------------------------------------!
  ! Written by Robert Charlton, 11/07/2018.                                    !
  !============================================================================!

  integer function sparse_first_atom_in_group(par, first_proc, &
       last_proc, group_size)

    implicit none

    ! Arguments
    type(PARAL_INFO), intent(in), pointer :: par
    integer, intent(in) :: first_proc
    integer, intent(in), optional :: last_proc
    integer, intent(in), optional :: group_size

    ! Local variables
    integer :: iproc, loc_last_proc, iat

    if(present(last_proc)) then
       loc_last_proc = last_proc
    else if(present(group_size)) then
       loc_last_proc = first_proc + group_size - 1
    end if

    ! rc2013: loop over all the procs of this group
    do iproc=first_proc,loc_last_proc
       ! rc2013: exit once we find a proc with atoms
       iat = par%first_atom_on_proc(iproc)
       if(iat .gt. 0) exit
    end do
    sparse_first_atom_in_group = iat

  end function sparse_first_atom_in_group

  !============================================================================!
  ! Returns the index of the last atom in a given comms group of procs.        !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !    par        (in) : Parallel strategy that we wish to probe.              !
  !    first_proc (in) : The first proc in this comms group.                   !
  !    last_proc  (in) : The last proc in this comms group.                    !
  !    group_size (in) : The (optional) size of this group.                    !
  !----------------------------------------------------------------------------!
  ! Written by Robert Charlton, 11/07/2018.                                    !
  !============================================================================!

  integer function sparse_last_atom_in_group(par, first_proc, last_proc)

    implicit none

    ! Arguments
    type(PARAL_INFO), intent(in), pointer :: par
    integer, intent(in) :: first_proc
    integer, intent(in), optional :: last_proc

    ! Local variables
    integer :: iproc, iat

    ! rc2013: loop over all the procs of this group
    do iproc=last_proc,first_proc,-1
       ! rc2013: exit once we find a proc with atoms
       iat = par%first_atom_on_proc(iproc) + par%num_atoms_on_proc(iproc) - 1
       if(iat .gt. 0) exit
    end do
    sparse_last_atom_in_group = iat

  end function sparse_last_atom_in_group


end module sparse

!==========================================================================!
! These subroutines override d/zgemm when nrows, ncols, nsum are all in    !
! the set {1,4,9}                                                          !
! They perform the operation c = c + a.b.                                  !
!--------------------------------------------------------------------------!
! Arguments:                                                               !
!   ncol,nrow : Dimensions of output matrix c,                             !
!                also number of rows of b and cols of a respectively       !
!   nsum      : Number of cols of a and rows of b                          !
!   a (input) : The block matrix a                                         !
!   a (input) : The block matrix b                                         !
!   c (inout) : The block matrix c                                         !
!--------------------------------------------------------------------------!
! This table shows the 27 types of block multiply that might be happening: !
! type nrow ncol nsum  matrix multiply                                     !
!    1    1    1    1  (1,1)*(1,1)->(1,1)                                  !
!    2    1    1    4  (1,4)*(4,1)->(1,1)                                  !
!    3    1    1    9  (1,9)*(9,1)->(1,1)                                  !
!    4    1    4    1  (1,1)*(1,4)->(1,4)                                  !
!    5    1    4    4  (1,4)*(4,4)->(1,4)                                  !
!    6    1    4    9  (1,9)*(9,4)->(1,4)                                  !
!    7    1    9    1  (1,1)*(1,9)->(1,9)                                  !
!    8    1    9    4  (1,4)*(4,9)->(1,9)                                  !
!    9    1    9    9  (1,9)*(9,9)->(1,9)                                  !
!   10    4    1    1  (4,1)*(1,1)->(4,1)                                  !
!   11    4    1    4  (4,4)*(4,1)->(4,1)                                  !
!   12    4    1    9  (4,9)*(9,1)->(4,1)                                  !
!   13    4    4    1  (4,1)*(1,4)->(4,4)                                  !
!   14    4    4    4  (4,4)*(4,4)->(4,4)                                  !
!   15    4    4    9  (4,9)*(9,4)->(4,4)                                  !
!   16    4    9    1  (4,1)*(1,9)->(4,9)                                  !
!   17    4    9    4  (4,4)*(4,9)->(4,9)                                  !
!   18    4    9    9  (4,9)*(9,9)->(4,9)                                  !
!   19    9    1    1  (9,1)*(1,1)->(9,1)                                  !
!   20    9    1    4  (9,4)*(4,1)->(9,1)                                  !
!   21    9    1    9  (9,9)*(9,1)->(9,1)                                  !
!   22    9    4    1  (9,1)*(1,4)->(9,4)                                  !
!   23    9    4    4  (9,4)*(4,4)->(9,4)                                  !
!   24    9    4    9  (9,9)*(9,4)->(9,4)                                  !
!   25    9    9    1  (9,1)*(1,9)->(9,9)                                  !
!   26    9    9    4  (9,4)*(4,9)->(9,9)                                  !
!   27    9    9    9  (9,9)*(9,9)->(9,9)                                  !
!--------------------------------------------------------------------------!
! Written by Nicholas Hine, November 2008                                  !
!==========================================================================!
subroutine sparse_dgemm_149(nrow,ncol,nsum,a,lda,b,ldb,c,ldc)

  use constants, only: DP

  implicit none

  ! Arguments
  integer,intent(in) :: nrow,ncol,nsum
  integer,intent(in) :: lda,ldb,ldc
  real(kind=DP),intent(in) :: a(lda,9),b(ldb,9)
  real(kind=DP),intent(inout) :: c(ldc,9)

  ! Locals
  integer :: i,j

  if      (nrow==1.and.ncol==1.and.nsum==1) then ! 1
     c(1,1)=c(1,1)+a(1,1)*b(1,1)
  else if (nrow==1.and.ncol==1.and.nsum==4) then ! 2
     c(1,1)=c(1,1)+sum(a(1,1:4)*b(1:4,1))
  else if (nrow==1.and.ncol==1.and.nsum==9) then ! 3
     c(1,1)=c(1,1)+sum(a(1,1:9)*b(1:9,1))
  else if (nrow==1.and.ncol==4.and.nsum==1) then ! 4
     c(1,1:4)=c(1,1:4)+a(1,1)*b(1,1:4)
  else if (nrow==1.and.ncol==4.and.nsum==4) then ! 5
     c(1,1)=c(1,1)+sum(a(1,1:4)*b(1:4,1))
     c(1,2)=c(1,2)+sum(a(1,1:4)*b(1:4,2))
     c(1,3)=c(1,3)+sum(a(1,1:4)*b(1:4,3))
     c(1,4)=c(1,4)+sum(a(1,1:4)*b(1:4,4))
  else if (nrow==1.and.ncol==4.and.nsum==9) then ! 6
     c(1,1)=c(1,1)+sum(a(1,1:9)*b(1:9,1))
     c(1,2)=c(1,2)+sum(a(1,1:9)*b(1:9,2))
     c(1,3)=c(1,3)+sum(a(1,1:9)*b(1:9,3))
     c(1,4)=c(1,4)+sum(a(1,1:9)*b(1:9,4))
  else if (nrow==1.and.ncol==9.and.nsum==1) then ! 7
     c(1,1:9)=c(1,1:9)+a(1,1)*b(1,1:9)
  else if (nrow==4.and.ncol==1.and.nsum==1) then ! 10
     c(1:4,1)=c(1:4,1)+a(1:4,1)*b(1,1)
  else if (nrow==4.and.ncol==1.and.nsum==4) then ! 11
     c(1,1)=c(1,1)+sum(a(1,1:4)*b(1:4,1))
     c(2,1)=c(2,1)+sum(a(2,1:4)*b(1:4,1))
     c(3,1)=c(3,1)+sum(a(3,1:4)*b(1:4,1))
     c(4,1)=c(4,1)+sum(a(4,1:4)*b(1:4,1))
  else if (nrow==4.and.ncol==1.and.nsum==9) then ! 12
     c(1,1)=c(1,1)+sum(a(1,1:9)*b(1:9,1))
     c(2,1)=c(2,1)+sum(a(2,1:9)*b(1:9,1))
     c(3,1)=c(3,1)+sum(a(3,1:9)*b(1:9,1))
     c(4,1)=c(4,1)+sum(a(4,1:9)*b(1:9,1))
  else if (nrow==4.and.ncol==4.and.nsum==1) then ! 13
     c(1:4,1)=c(1:4,1)+a(1:4,1)*b(1,1)
     c(1:4,2)=c(1:4,2)+a(1:4,1)*b(1,2)
     c(1:4,3)=c(1:4,3)+a(1:4,1)*b(1,3)
     c(1:4,4)=c(1:4,4)+a(1:4,1)*b(1,4)
  else if (nrow==4.and.ncol==4.and.nsum==4) then ! 14
     c(1,1)=c(1,1)+sum(a(1,1:4)*b(1:4,1))
     c(2,1)=c(2,1)+sum(a(2,1:4)*b(1:4,1))
     c(3,1)=c(3,1)+sum(a(3,1:4)*b(1:4,1))
     c(4,1)=c(4,1)+sum(a(4,1:4)*b(1:4,1))
     c(1,2)=c(1,2)+sum(a(1,1:4)*b(1:4,2))
     c(2,2)=c(2,2)+sum(a(2,1:4)*b(1:4,2))
     c(3,2)=c(3,2)+sum(a(3,1:4)*b(1:4,2))
     c(4,2)=c(4,2)+sum(a(4,1:4)*b(1:4,2))
     c(1,3)=c(1,3)+sum(a(1,1:4)*b(1:4,3))
     c(2,3)=c(2,3)+sum(a(2,1:4)*b(1:4,3))
     c(3,3)=c(3,3)+sum(a(3,1:4)*b(1:4,3))
     c(4,3)=c(4,3)+sum(a(4,1:4)*b(1:4,3))
     c(1,4)=c(1,4)+sum(a(1,1:4)*b(1:4,4))
     c(2,4)=c(2,4)+sum(a(2,1:4)*b(1:4,4))
     c(3,4)=c(3,4)+sum(a(3,1:4)*b(1:4,4))
     c(4,4)=c(4,4)+sum(a(4,1:4)*b(1:4,4))
  else if (nrow==4.and.ncol==4.and.nsum==9) then ! 15
     c(1,1)=c(1,1)+sum(a(1,1:9)*b(1:9,1))
     c(2,1)=c(2,1)+sum(a(2,1:9)*b(1:9,1))
     c(3,1)=c(3,1)+sum(a(3,1:9)*b(1:9,1))
     c(4,1)=c(4,1)+sum(a(4,1:9)*b(1:9,1))
     c(1,2)=c(1,2)+sum(a(1,1:9)*b(1:9,2))
     c(2,2)=c(2,2)+sum(a(2,1:9)*b(1:9,2))
     c(3,2)=c(3,2)+sum(a(3,1:9)*b(1:9,2))
     c(4,2)=c(4,2)+sum(a(4,1:9)*b(1:9,2))
     c(1,3)=c(1,3)+sum(a(1,1:9)*b(1:9,3))
     c(2,3)=c(2,3)+sum(a(2,1:9)*b(1:9,3))
     c(3,3)=c(3,3)+sum(a(3,1:9)*b(1:9,3))
     c(4,3)=c(4,3)+sum(a(4,1:9)*b(1:9,3))
     c(1,4)=c(1,4)+sum(a(1,1:9)*b(1:9,4))
     c(2,4)=c(2,4)+sum(a(2,1:9)*b(1:9,4))
     c(3,4)=c(3,4)+sum(a(3,1:9)*b(1:9,4))
     c(4,4)=c(4,4)+sum(a(4,1:9)*b(1:9,4))
  else if (nrow==9.and.ncol==1.and.nsum==1) then ! 19
     c(1:9,1)=c(1:9,1)+a(1:9,1)*b(1,1)
  else !  any of the rest
     do i=1,ncol
        do j=1,nrow
           c(j,i) = c(j,i) + sum(a(j,1:nsum)*b(1:nsum,i))
        end do
     end do
  end if

end subroutine sparse_dgemm_149

subroutine sparse_zgemm_149(nrow,ncol,nsum,a,lda,b,ldb,c,ldc)

  use constants, only: DP

  implicit none

  ! Arguments
  integer,intent(in) :: nrow,ncol,nsum
  integer,intent(in) :: lda,ldb,ldc
  complex(kind=DP),intent(in) :: a(lda,9),b(ldb,9)
  complex(kind=DP),intent(inout) :: c(ldc,9)

  ! Locals
  integer :: i,j

  if      (nrow==1.and.ncol==1.and.nsum==1) then ! 1
     c(1,1)=c(1,1)+a(1,1)*b(1,1)
  else if (nrow==1.and.ncol==1.and.nsum==4) then ! 2
     c(1,1)=c(1,1)+sum(a(1,1:4)*b(1:4,1))
  else if (nrow==1.and.ncol==1.and.nsum==9) then ! 3
     c(1,1)=c(1,1)+sum(a(1,1:9)*b(1:9,1))
  else if (nrow==1.and.ncol==4.and.nsum==1) then ! 4
     c(1,1:4)=c(1,1:4)+a(1,1)*b(1,1:4)
  else if (nrow==1.and.ncol==4.and.nsum==4) then ! 5
     c(1,1)=c(1,1)+sum(a(1,1:4)*b(1:4,1))
     c(1,2)=c(1,2)+sum(a(1,1:4)*b(1:4,2))
     c(1,3)=c(1,3)+sum(a(1,1:4)*b(1:4,3))
     c(1,4)=c(1,4)+sum(a(1,1:4)*b(1:4,4))
  else if (nrow==1.and.ncol==4.and.nsum==9) then ! 6
     c(1,1)=c(1,1)+sum(a(1,1:9)*b(1:9,1))
     c(1,2)=c(1,2)+sum(a(1,1:9)*b(1:9,2))
     c(1,3)=c(1,3)+sum(a(1,1:9)*b(1:9,3))
     c(1,4)=c(1,4)+sum(a(1,1:9)*b(1:9,4))
  else if (nrow==1.and.ncol==9.and.nsum==1) then ! 7
     c(1,1:9)=c(1,1:9)+a(1,1)*b(1,1:9)
  else if (nrow==4.and.ncol==1.and.nsum==1) then ! 10
     c(1:4,1)=c(1:4,1)+a(1:4,1)*b(1,1)
  else if (nrow==4.and.ncol==1.and.nsum==4) then ! 11
     c(1,1)=c(1,1)+sum(a(1,1:4)*b(1:4,1))
     c(2,1)=c(2,1)+sum(a(2,1:4)*b(1:4,1))
     c(3,1)=c(3,1)+sum(a(3,1:4)*b(1:4,1))
     c(4,1)=c(4,1)+sum(a(4,1:4)*b(1:4,1))
  else if (nrow==4.and.ncol==1.and.nsum==9) then ! 12
     c(1,1)=c(1,1)+sum(a(1,1:9)*b(1:9,1))
     c(2,1)=c(2,1)+sum(a(2,1:9)*b(1:9,1))
     c(3,1)=c(3,1)+sum(a(3,1:9)*b(1:9,1))
     c(4,1)=c(4,1)+sum(a(4,1:9)*b(1:9,1))
  else if (nrow==4.and.ncol==4.and.nsum==1) then ! 13
     c(1:4,1)=c(1:4,1)+a(1:4,1)*b(1,1)
     c(1:4,2)=c(1:4,2)+a(1:4,1)*b(1,2)
     c(1:4,3)=c(1:4,3)+a(1:4,1)*b(1,3)
     c(1:4,4)=c(1:4,4)+a(1:4,1)*b(1,4)
  else if (nrow==4.and.ncol==4.and.nsum==4) then ! 14
     c(1,1)=c(1,1)+sum(a(1,1:4)*b(1:4,1))
     c(2,1)=c(2,1)+sum(a(2,1:4)*b(1:4,1))
     c(3,1)=c(3,1)+sum(a(3,1:4)*b(1:4,1))
     c(4,1)=c(4,1)+sum(a(4,1:4)*b(1:4,1))
     c(1,2)=c(1,2)+sum(a(1,1:4)*b(1:4,2))
     c(2,2)=c(2,2)+sum(a(2,1:4)*b(1:4,2))
     c(3,2)=c(3,2)+sum(a(3,1:4)*b(1:4,2))
     c(4,2)=c(4,2)+sum(a(4,1:4)*b(1:4,2))
     c(1,3)=c(1,3)+sum(a(1,1:4)*b(1:4,3))
     c(2,3)=c(2,3)+sum(a(2,1:4)*b(1:4,3))
     c(3,3)=c(3,3)+sum(a(3,1:4)*b(1:4,3))
     c(4,3)=c(4,3)+sum(a(4,1:4)*b(1:4,3))
     c(1,4)=c(1,4)+sum(a(1,1:4)*b(1:4,4))
     c(2,4)=c(2,4)+sum(a(2,1:4)*b(1:4,4))
     c(3,4)=c(3,4)+sum(a(3,1:4)*b(1:4,4))
     c(4,4)=c(4,4)+sum(a(4,1:4)*b(1:4,4))
  else if (nrow==4.and.ncol==4.and.nsum==9) then ! 15
     c(1,1)=c(1,1)+sum(a(1,1:9)*b(1:9,1))
     c(2,1)=c(2,1)+sum(a(2,1:9)*b(1:9,1))
     c(3,1)=c(3,1)+sum(a(3,1:9)*b(1:9,1))
     c(4,1)=c(4,1)+sum(a(4,1:9)*b(1:9,1))
     c(1,2)=c(1,2)+sum(a(1,1:9)*b(1:9,2))
     c(2,2)=c(2,2)+sum(a(2,1:9)*b(1:9,2))
     c(3,2)=c(3,2)+sum(a(3,1:9)*b(1:9,2))
     c(4,2)=c(4,2)+sum(a(4,1:9)*b(1:9,2))
     c(1,3)=c(1,3)+sum(a(1,1:9)*b(1:9,3))
     c(2,3)=c(2,3)+sum(a(2,1:9)*b(1:9,3))
     c(3,3)=c(3,3)+sum(a(3,1:9)*b(1:9,3))
     c(4,3)=c(4,3)+sum(a(4,1:9)*b(1:9,3))
     c(1,4)=c(1,4)+sum(a(1,1:9)*b(1:9,4))
     c(2,4)=c(2,4)+sum(a(2,1:9)*b(1:9,4))
     c(3,4)=c(3,4)+sum(a(3,1:9)*b(1:9,4))
     c(4,4)=c(4,4)+sum(a(4,1:9)*b(1:9,4))
  else if (nrow==9.and.ncol==1.and.nsum==1) then ! 19
     c(1:9,1)=c(1:9,1)+a(1:9,1)*b(1,1)
  else !  any of the rest
     do i=1,ncol
        do j=1,nrow
           c(j,i) = c(j,i) + sum(a(j,1:nsum)*b(1:nsum,i))
        end do
     end do
  end if

end subroutine sparse_zgemm_149

