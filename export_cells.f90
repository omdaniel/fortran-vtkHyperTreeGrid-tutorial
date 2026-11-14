!! export_cells: write every cell of both examples to plain text, in
!! depth-first (build) order, so that make_figures.py can draw them.
!!
!! Columns: tree  path  level  refined  x0  y0  h  value
!! "path" is the string of child indices from the root ("r" for the root).
program export_cells
  use quadtree_mod
  use field_mod, only: bump_unit, bump_forest
  implicit none

  real(dp), parameter :: tol = 0.10_dp
  integer, parameter :: nx = 4, ny = 2
  type(qnode) :: root, tree
  integer :: u, i, j

  ! Single tree on the unit square (same setup as single_tree.f90).
  root = make_root(0.0_dp, 0.0_dp, 1.0_dp)
  call refine(root, bump_unit, tol)
  open (newunit=u, file='cells_single.txt', status='replace', action='write')
  write (u, '(a)') '# tree path level refined x0 y0 h value'
  call write_cells(root, 0, '')
  close (u)

  ! Forest on [0,4] x [0,2] (same setup as forest.f90), tree t = i + nx*j.
  open (newunit=u, file='cells_forest.txt', status='replace', action='write')
  write (u, '(a)') '# tree path level refined x0 y0 h value'
  do j = 0, ny - 1
    do i = 0, nx - 1
      tree = make_root(real(i, dp), real(j, dp), 1.0_dp)
      call refine(tree, bump_forest, tol)
      call write_cells(tree, i + nx*j, '')
    end do
  end do
  close (u)

contains

  recursive subroutine write_cells(node, t, path)
    type(qnode), intent(in) :: node
    integer, intent(in) :: t
    character(*), intent(in) :: path
    integer :: ic
    character(len=1) :: digit
    write (u, '(i0, 1x, a, 1x, i0, 1x, i0, 4(1x, es23.15e3))') t, 'r' // path, node%level, &
      merge(1, 0, allocated(node%child)), node%x0, node%y0, node%h, node%value
    if (allocated(node%child)) then
      do ic = 0, NCHILD - 1
        write (digit, '(i1)') ic
        call write_cells(node%child(ic), t, path // digit)
      end do
    end if
  end subroutine write_cells

end program export_cells
