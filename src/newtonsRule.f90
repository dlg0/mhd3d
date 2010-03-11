function newtonsRule ( r0, u1, u3, theta0, rI )
  implicit none

	integer, parameter ::  DBL = selected_real_kind ( p = 13, r = 200 )
	real(kind=DBL) :: newtonsRule
  real(kind=DBL), intent ( in ) :: r0, u1, u3, theta0, rI
  real(kind=DBL) :: tol, errV, fR, dFdR, ans, r
  integer, parameter :: maxN = 100
  real(kind=DBL), dimension ( maxN ) :: rStorage
  integer :: i

  tol = 1e-3
  errV = 1.0
  r = r0

  i = 1
  do
    
    rStorage(i) = r

    fR  = u3 ** 2 * r ** 4 * cos ( theta0 ) ** 2 / rI ** 3 &
      - u1 * r - rI

    dFdR  = 4.0 * u3 ** 2 * r ** 3 * cos ( theta0 ) ** 2 / rI ** 3 &
      - u1;

    if ( dFdR == 0. ) then

      newtonsRule = r
      return

    end if

    ans = r - fR / dFdR
    errV = ans - r

    if ( errV < 0. ) errV = -errV

    r  = ans

    i = i + 1

    if ( i >= maxN .or. errV < tol ) exit

  end do 

  newtonsRule = ans
  return

end function newtonsRule
