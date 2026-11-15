!! quadtree_mod: a minimal adaptive quadtree, built depth first, and its
!! breadth-first (vtkHyperTreeGrid) serialisation.
!!
!! Conventions follow VTK's vtkHyperTree / vtkHyperTreeGrid:
!!   * child index ic = ix + 2*iy, with ix, iy in {0,1}  (x is the fastest bit)
!!   * every cell (refined or leaf) carries a value
!!   * the descriptor holds one bit per cell, breadth first, level by level,
!!     1 = refined, 0 = leaf, and the deepest level is omitted
!!   * NumberOfVerticesPerDepth holds the cell count at each depth
!!
!! Requires a Fortran 2008 compiler (recursive allocatable components).
module quadtree_mod
  use, intrinsic :: iso_fortran_env, only: dp => real64
  implicit none
  private

  public :: dp, field2d
  public :: qnode, make_root, refine, print_depth_first, draw_leaf_boxes
  public :: htree_bf, to_breadth_first, to_breadth_first_queue, check_bf
  public :: bits_to_string, pack_bits_msb_first

  !> Deepest refinement level a leaf may reach (the root is level 0).
  integer, parameter, public :: MAX_LEVEL = 3
  integer, parameter, public :: NCHILD = 4

  abstract interface
    pure function field2d(x, y) result(f)
      import :: dp
      real(dp), intent(in) :: x, y
      real(dp) :: f
    end function field2d
  end interface

  !> One cell of the tree. A leaf has no children allocated.
  type :: qnode
    integer  :: level = 0
    real(dp) :: x0 = 0, y0 = 0      ! lower-left corner
    real(dp) :: h = 1               ! edge length (square cells)
    real(dp) :: value = 0           ! f sampled at the cell centre
    type(qnode), allocatable :: child(:)   ! child(0:3) when refined
  end type qnode

  !> One tree in breadth-first form, as vtkHyperTreeGrid stores it.
  type :: htree_bf
    integer :: depth = 0                         ! number of levels present
    integer,  allocatable :: nverts_per_depth(:) ! NumberOfVerticesPerDepth
    logical,  allocatable :: descriptor(:)       ! levels 0 .. depth-2
    real(dp), allocatable :: values(:)           ! cell data, breadth first
  end type htree_bf

  ! Work space for the depth-first to breadth-first conversion.
  type :: depth_bucket
    logical,  allocatable :: refined(:)
    real(dp), allocatable :: value(:)
  end type depth_bucket

  ! A pointer wrapper, so that a queue can hold references to nodes.
  type :: node_ptr
    type(qnode), pointer :: p => null()
  end type node_ptr

contains

  !-------------------------------------------------------------------------
  ! Building the tree, depth first
  !-------------------------------------------------------------------------

  pure function make_root(x0, y0, h) result(root)
    real(dp), intent(in) :: x0, y0, h
    type(qnode) :: root
    root%level = 0
    root%x0 = x0
    root%y0 = y0
    root%h = h
  end function make_root

  !> Refinement test: sample f at the four corners and the centre and refine
  !> when the spread of those five samples exceeds tol. Five samples per cell
  !> is exactly what makes this test blind to features narrower than the cell.
  logical function needs_refinement(node, f, tol)
    type(qnode), intent(in) :: node
    procedure(field2d) :: f
    real(dp), intent(in) :: tol
    real(dp) :: s(5)
    associate (x0 => node%x0, y0 => node%y0, h => node%h)
      s = [ f(x0, y0), f(x0 + h, y0), f(x0, y0 + h), f(x0 + h, y0 + h), &
            f(x0 + 0.5_dp*h, y0 + 0.5_dp*h) ]
    end associate
    needs_refinement = (maxval(s) - minval(s)) > tol
  end function needs_refinement

  !> Depth-first construction: a child's whole subtree is finished before
  !> its next sibling is started.
  recursive subroutine refine(node, f, tol)
    type(qnode), intent(inout) :: node
    procedure(field2d) :: f
    real(dp), intent(in) :: tol
    integer :: ic, ix, iy

    node%value = f(node%x0 + 0.5_dp*node%h, node%y0 + 0.5_dp*node%h)
    if (node%level >= MAX_LEVEL) return
    if (.not. needs_refinement(node, f, tol)) return

    allocate (node%child(0:NCHILD - 1))
    do ic = 0, NCHILD - 1
      ix = mod(ic, 2)   ! bit 0: -x (0) or +x (1)
      iy = ic / 2       ! bit 1: -y (0) or +y (1)
      associate (c => node%child(ic))
        c%level = node%level + 1
        c%h  = 0.5_dp*node%h
        c%x0 = node%x0 + ix*c%h
        c%y0 = node%y0 + iy*c%h
        call refine(c, f, tol)
      end associate
    end do
  end subroutine refine

  !> Print the tree in depth-first pre-order: the order it was built in.
  recursive subroutine print_depth_first(node, unit, path)
    type(qnode), intent(in) :: node
    integer, intent(in) :: unit
    character(*), intent(in) :: path
    integer :: ic
    character(len=1) :: digit
    character(len=:), allocatable :: label

    label = path
    if (len(path) == 0) label = 'root'
    write (unit, '(a, a, a, i0, a, f8.5, a)') repeat('  ', node%level), &
      merge('+ ', '- ', allocated(node%child)), label // '  level ', &
      node%level, '  f(centre) = ', node%value, &
      merge('  refined', '         ', allocated(node%child))
    if (allocated(node%child)) then
      do ic = 0, NCHILD - 1
        write (digit, '(i1)') ic
        call print_depth_first(node%child(ic), unit, path // digit)
      end do
    end if
  end subroutine print_depth_first

  !> Draw every leaf as a box at its true size, with its level in the middle.
  !> canvas(c, r) is the picture: column c runs along x, row r along y
  !> (r = 0 at the bottom). Each finest cell is cw characters wide and 2 rows
  !> tall, and (i0, j0) is the offset of this cell in finest cells.
  !> Edges are drawn only over blanks and corners always, so a corner that
  !> lies on a larger neighbour's edge still comes out as +.
  recursive subroutine draw_leaf_boxes(node, canvas, cw, i0, j0)
    type(qnode), intent(in) :: node
    character, intent(inout) :: canvas(0:, 0:)
    integer, intent(in) :: cw, i0, j0
    integer :: ic, w, c0, c1, r0, r1, k
    w = 2**(MAX_LEVEL - node%level)     ! width of this cell in finest cells
    if (allocated(node%child)) then
      do ic = 0, NCHILD - 1
        call draw_leaf_boxes(node%child(ic), canvas, cw, &
                             i0 + mod(ic, 2)*w/2, j0 + (ic/2)*w/2)
      end do
      return
    end if
    c0 = cw*i0
    c1 = cw*(i0 + w)
    r0 = 2*j0
    r1 = 2*(j0 + w)
    do k = c0 + 1, c1 - 1
      if (canvas(k, r0) == ' ') canvas(k, r0) = '-'
      if (canvas(k, r1) == ' ') canvas(k, r1) = '-'
    end do
    do k = r0 + 1, r1 - 1
      if (canvas(c0, k) == ' ') canvas(c0, k) = '|'
      if (canvas(c1, k) == ' ') canvas(c1, k) = '|'
    end do
    canvas([c0, c1], r0) = '+'
    canvas([c0, c1], r1) = '+'
    canvas((c0 + c1)/2, (r0 + r1)/2) = achar(iachar('0') + node%level)
  end subroutine draw_leaf_boxes

  !-------------------------------------------------------------------------
  ! Depth first to breadth first
  !-------------------------------------------------------------------------

  !> The conversion VTK itself uses (vtkHyperTree::ComputeBreadthFirstOrderDescriptor):
  !> walk the tree depth first, but append each cell to a list for its depth.
  !> Within one depth, a depth-first walk meets cells left to right, which
  !> is exactly breadth-first order, so concatenating the lists gives the
  !> breadth-first layout with no queue.
  function to_breadth_first(root) result(bf)
    type(qnode), intent(in) :: root
    type(htree_bf) :: bf
    type(depth_bucket) :: b(0:MAX_LEVEL)
    integer :: d

    do d = 0, MAX_LEVEL
      allocate (b(d)%refined(0), b(d)%value(0))
    end do
    call collect(root, b)

    bf%depth = count([(size(b(d)%value) > 0, d = 0, MAX_LEVEL)])
    allocate (bf%nverts_per_depth(0), bf%descriptor(0), bf%values(0))
    do d = 0, bf%depth - 1
      bf%nverts_per_depth = [bf%nverts_per_depth, size(b(d)%value)]
      bf%values = [bf%values, b(d)%value]
      ! The deepest level is never refined, so it is left out.
      if (d < bf%depth - 1) bf%descriptor = [bf%descriptor, b(d)%refined]
    end do
  end function to_breadth_first

  recursive subroutine collect(node, b)
    type(qnode), intent(in) :: node
    type(depth_bucket), intent(inout) :: b(0:)
    integer :: ic
    b(node%level)%value   = [b(node%level)%value, node%value]
    b(node%level)%refined = [b(node%level)%refined, allocated(node%child)]
    if (allocated(node%child)) then
      do ic = 0, NCHILD - 1
        call collect(node%child(ic), b)
      end do
    end if
  end subroutine collect

  !> The textbook alternative, used here only as a cross-check: a FIFO queue.
  !> The queue holds pointers into the tree, so root needs the target
  !> attribute here and in the caller. The tree is only read through them.
  function to_breadth_first_queue(root) result(bf)
    type(qnode), intent(in), target :: root
    type(htree_bf) :: bf
    type(node_ptr), allocatable :: queue(:)
    type(qnode), pointer :: n
    logical, allocatable :: refined(:)
    integer, allocatable :: level(:)
    integer :: head, ic, d, nlast

    allocate (queue(1), refined(0), level(0), bf%values(0))
    queue(1)%p => root
    head = 1
    do while (head <= size(queue))
      n => queue(head)%p          ! dequeue
      bf%values = [bf%values, n%value]
      refined = [refined, allocated(n%child)]
      level = [level, n%level]
      if (allocated(n%child)) then
        do ic = 0, NCHILD - 1     ! enqueue the children, in VTK child order
          queue = [queue, node_ptr(n%child(ic))]
        end do
      end if
      head = head + 1
    end do

    bf%depth = maxval(level) + 1
    allocate (bf%nverts_per_depth(bf%depth))
    do d = 0, bf%depth - 1
      bf%nverts_per_depth(d + 1) = count(level == d)
    end do
    nlast = bf%nverts_per_depth(bf%depth)
    allocate (bf%descriptor, source=refined(1:size(refined) - nlast))
  end function to_breadth_first_queue

  !> Consistency of a breadth-first tree: each refined cell at depth d
  !> contributes 4 vertices at depth d+1, and the sizes add up.
  logical function check_bf(bf) result(ok)
    type(htree_bf), intent(in) :: bf
    integer :: d, first, last
    ok = bf%nverts_per_depth(1) == 1
    ok = ok .and. size(bf%values) == sum(bf%nverts_per_depth)
    ok = ok .and. size(bf%descriptor) == sum(bf%nverts_per_depth(1:bf%depth - 1))
    first = 1
    do d = 1, bf%depth - 1
      last = first + bf%nverts_per_depth(d) - 1
      ok = ok .and. NCHILD*count(bf%descriptor(first:last)) == bf%nverts_per_depth(d + 1)
      first = last + 1
    end do
  end function check_bf

  !-------------------------------------------------------------------------
  ! Bit helpers
  !-------------------------------------------------------------------------

  pure function bits_to_string(bits) result(s)
    logical, intent(in) :: bits(:)
    character(len=size(bits)) :: s
    integer :: i
    do i = 1, size(bits)
      s(i:i) = merge('1', '0', bits(i))
    end do
  end function bits_to_string

  !> Pack bits into bytes the way vtkBitArray does: bit id k lives in byte
  !> k/8 under the mask 0x80 >> mod(k, 8), i.e. most significant bit first.
  pure function pack_bits_msb_first(bits) result(bytes)
    logical, intent(in) :: bits(:)
    integer, allocatable :: bytes(:)
    integer :: k
    allocate (bytes((size(bits) + 7)/8), source=0)
    do k = 0, size(bits) - 1
      if (bits(k + 1)) bytes(k/8 + 1) = ior(bytes(k/8 + 1), ishft(128, -mod(k, 8)))
    end do
  end function pack_bits_msb_first

end module quadtree_mod
