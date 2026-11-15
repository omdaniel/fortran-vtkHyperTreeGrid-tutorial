!! Part 1 and 3 of the tutorial: one quadtree over the unit square,
!! built depth first, then laid out breadth first as vtkHyperTreeGrid does.
program single_tree
  use, intrinsic :: iso_fortran_env, only: output_unit
  use quadtree_mod
  use field_mod, only: bump_unit
  implicit none

  real(dp), parameter :: tol = 0.10_dp
  integer, parameter :: n = 2**MAX_LEVEL
  type(qnode), target :: root
  type(htree_bf) :: bf, bfq
  character :: canvas(0:4*n, 0:2*n)
  integer :: r, d, first, last
  integer, allocatable :: bytes(:)

  ! 1. Build depth first.
  root = make_root(0.0_dp, 0.0_dp, 1.0_dp)
  call refine(root, bump_unit, tol)

  print '(a)', '== Depth-first (build) order =='
  call print_depth_first(root, output_unit, '')

  ! 2. Where the leaves ended up: each leaf drawn at its true size.
  canvas = ' '
  call draw_leaf_boxes(root, canvas, 4, 0, 0)
  print '(/a)', '== Leaf cells, labelled by level (y up, x right) =='
  do r = 2*n, 0, -1
    print '(*(a))', canvas(:, r)
  end do

  ! 3. Breadth-first layout.
  bf = to_breadth_first(root)
  bfq = to_breadth_first_queue(root)

  print '(/a)', '== Breadth-first (vtkHyperTreeGrid) layout =='
  print '(a, i0)', 'depth (levels present)    : ', bf%depth
  print '(a, *(i0, :, 1x))', 'NumberOfVerticesPerDepth  : ', bf%nverts_per_depth
  print '(a, i0)', 'total cells               : ', size(bf%values)
  print '(a, i0)', 'leaves                    : ', &
    size(bf%values) - count(bf%descriptor)
  print '(a, i0, a)', 'descriptor bits           : ', size(bf%descriptor), &
    '  (deepest level omitted)'
  first = 1
  do d = 1, bf%depth - 1
    last = first + bf%nverts_per_depth(d) - 1
    print '(a, i0, a, a)', '  depth ', d - 1, ': ', bits_to_string(bf%descriptor(first:last))
    first = last + 1
  end do
  bytes = pack_bits_msb_first(bf%descriptor)
  print '(a, *(z2.2, :, 1x))', 'descriptor bytes (hex)    : ', bytes

  print '(a)', 'cell values, breadth first, by depth:'
  first = 1
  do d = 1, bf%depth
    last = first + bf%nverts_per_depth(d) - 1
    print '(a, i0, a, *(f6.3, :, 1x))', '  depth ', d - 1, ': ', bf%values(first:min(last, first + 7))
    if (last > first + 7) print '(a, i0, a)', '           ... (', last - first - 7, ' more)'
    first = last + 1
  end do

  ! 4. Cross-checks.
  print '(/a)', '== Checks =='
  print '(a, l1)', 'layout is self-consistent         : ', check_bf(bf)
  print '(a, l1)', 'bucket walk == FIFO queue (bits)  : ', &
    all(bf%descriptor .eqv. bfq%descriptor) .and. size(bf%descriptor) == size(bfq%descriptor)
  print '(a, l1)', 'bucket walk == FIFO queue (values): ', &
    all(bf%values == bfq%values) .and. size(bf%values) == size(bfq%values)
end program single_tree
