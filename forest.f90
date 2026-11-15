!! Part 2 and 3 of the tutorial: a forest of quadtrees on a 4 x 2 grid of
!! unit root cells, and the global arrays vtkHyperTreeGrid (XML v2 / vtkHDF)
!! expects for it.
program forest
  use quadtree_mod
  use field_mod, only: bump_forest
  implicit none

  integer, parameter :: nx = 4, ny = 2            ! root cells per direction
  integer, parameter :: ntrees = nx*ny
  integer, parameter :: m = 2**MAX_LEVEL          ! finest cells per root edge
  real(dp), parameter :: tol = 0.10_dp

  ! Coarse grid: Dimensions counts points, so cells = points - 1 per axis.
  real(dp) :: xcoords(nx + 1), ycoords(ny + 1)
  type(qnode)    :: trees(0:ntrees - 1)
  type(htree_bf) :: bf(0:ntrees - 1)
  integer :: i, j, t
  integer, parameter :: cw = 2                   ! characters per finest cell
  character :: canvas(0:cw*nx*m, 0:2*ny*m)
  integer, allocatable :: tree_ids(:), depth_per_tree(:), nverts(:), global_start(:)
  logical, allocatable :: descriptors(:)
  real(dp), allocatable :: cell_data(:)

  xcoords = [(real(i, dp), i = 0, nx)]
  ycoords = [(real(j, dp), j = 0, ny)]

  ! 1. One tree per root cell, each built depth first.
  !    Tree index t = i + nx*j: i (x) fastest, VTK's default root indexing.
  do j = 0, ny - 1
    do i = 0, nx - 1
      t = i + nx*j
      trees(t) = make_root(xcoords(i + 1), ycoords(j + 1), xcoords(i + 2) - xcoords(i + 1))
      call refine(trees(t), bump_forest, tol)
    end do
  end do

  canvas = ' '
  do j = 0, ny - 1
    do i = 0, nx - 1
      call draw_leaf_boxes(trees(i + nx*j), canvas, cw, i*m, j*m)
    end do
  end do
  print '(a)', '== Leaf cells of the forest, labelled by level (y up, x right) =='
  do j = 2*ny*m, 0, -1
    print '(*(a))', canvas(:, j)
  end do

  ! 2. Each tree in breadth-first form.
  print '(/a)', '== Per tree =='
  print '(a)', ' tree (i,j)  depth  cells/depth          descriptor'
  do t = 0, ntrees - 1
    bf(t) = to_breadth_first(trees(t))
    if (.not. check_bf(bf(t))) error stop 'inconsistent tree'
    print '(i5, " (", i0, ",", i0, ")", i6, 3x, a20, 1x, a)', t, mod(t, nx), t/nx, bf(t)%depth, &
      ints(bf(t)%nverts_per_depth), by_depth(bf(t))
  end do

  ! 3. Concatenate into the global arrays, trees in TreeIds order.
  allocate (tree_ids(0), depth_per_tree(0), nverts(0), global_start(0), &
            descriptors(0), cell_data(0))
  do t = 0, ntrees - 1
    tree_ids       = [tree_ids, t]
    depth_per_tree = [depth_per_tree, bf(t)%depth]
    nverts         = [nverts, bf(t)%nverts_per_depth]
    global_start   = [global_start, size(cell_data)]    ! 0-based offset
    descriptors    = [descriptors, bf(t)%descriptor]
    cell_data      = [cell_data, bf(t)%values]
  end do

  print '(/a)', '== Global arrays for the whole grid =='
  print '(a, 3(i0, :, 1x))', 'Dimensions (points)        : ', nx + 1, ny + 1, 1
  print '(a, *(f4.1, :, 1x))', 'XCoordinates               : ', xcoords
  print '(a, *(f4.1, :, 1x))', 'YCoordinates               : ', ycoords
  print '(a)', 'BranchFactor               : 2'
  print '(a, *(i0, :, 1x))', 'TreeIds                    : ', tree_ids
  print '(a, *(i0, :, 1x))', 'DepthPerTree               : ', depth_per_tree
  print '(a, *(i0, :, 1x))', 'NumberOfVerticesPerDepth   : ', nverts
  print '(a, i0, a)', 'Descriptors                : ', size(descriptors), ' bits'
  print '(a, a)', '  ', bits_to_string(descriptors)
  print '(a, *(z2.2, :, 1x))', '  packed, MSB first (hex)  : ', pack_bits_msb_first(descriptors)
  print '(a, *(i0, :, 1x))', 'GlobalIndexStart per tree  : ', global_start
  print '(a, i0)', 'NumberOfCells (cell data)  : ', size(cell_data)
  print '(a, i0)', 'leaves                     : ', size(cell_data) - count(descriptors)

  print '(/a)', '== Locating one cell: the value at the centre of tree 1 =='
  ! Tree 1 is root cell (1,0). Its root is cell 0 of that tree, so its
  ! value sits at cell_data(global_start + 0), 1-based in Fortran.
  print '(a, f8.5)', 'cell_data(GlobalIndexStart(1) + 1) = ', cell_data(global_start(2) + 1)
  print '(a, f8.5)', 'f at the centre of root cell (1,0) = ', bump_forest(1.5_dp, 0.5_dp)

contains

  function ints(a) result(s)
    integer, intent(in) :: a(:)
    character(len=20) :: s
    integer :: k
    s = ''
    do k = 1, size(a)
      write (s(len_trim(s) + 2:), '(i0)') a(k)
    end do
    s = adjustl(s)
  end function ints

  !> The descriptor split by depth, with spaces between depths.
  function by_depth(b) result(s)
    type(htree_bf), intent(in) :: b
    character(len=:), allocatable :: s
    integer :: dd, f1, l1
    s = ''
    f1 = 1
    do dd = 1, b%depth - 1
      l1 = f1 + b%nverts_per_depth(dd) - 1
      s = s // bits_to_string(b%descriptor(f1:l1)) // ' '
      f1 = l1 + 1
    end do
    if (len(s) == 0) s = '(none: root is a leaf)'
  end function by_depth

end program forest
