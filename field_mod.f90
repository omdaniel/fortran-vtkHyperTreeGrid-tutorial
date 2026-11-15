!! field_mod: the cherry-picked function the tutorial maps.
!!
!! A single smooth Gaussian bump. Its width sigma is comparable to the root
!! cell size, so the five-point refinement test in quadtree_mod always "sees"
!! the bump from the coarsest level down. That is a deliberate choice: a bump
!! much narrower than a cell can fall between the samples and never trigger
!! refinement (aliasing). Nothing here guards against that.
module field_mod
  use quadtree_mod, only: dp
  implicit none
  private
  public :: bump_unit, bump_forest

contains

  !> Bump on the unit square, for the single-tree example.
  pure function bump_unit(x, y) result(f)
    real(dp), intent(in) :: x, y
    real(dp) :: f
    real(dp), parameter :: xc = 0.20_dp, yc = 0.90_dp, sigma = 0.15_dp
    f = exp(-((x - xc)**2 + (y - yc)**2) / (2*sigma**2))
  end function bump_unit

  !> Bump on [0,4] x [0,2], for the forest example (unit root cells).
  pure function bump_forest(x, y) result(f)
    real(dp), intent(in) :: x, y
    real(dp) :: f
    real(dp), parameter :: xc = 1.30_dp, yc = 0.85_dp, sigma = 0.35_dp
    f = exp(-((x - xc)**2 + (y - yc)**2) / (2*sigma**2))
  end function bump_forest

end module field_mod
