module constants
  implicit none
  save

	integer, parameter :: SGL	= selected_real_kind ( p = 6, r = 37 )
	integer, parameter :: DBL	= selected_real_kind ( p = 13, r = 200 )

  real(kind=DBL), parameter :: rE = 6356.75e3
  real(kind=DBL), parameter :: rI = 1.02
  real(kind=DBL), parameter :: pi = 3.141593
  real(kind=DBL), parameter :: degToRad = pi / 180.0
  real(kind=DBL), parameter :: radToDeg = 180.0 / pi
  real(kind=DBL), parameter :: u0_ = 1.2566e-6;
  real(kind=DBL), parameter :: e0_ = 8.8541878e-12;
  real(kind=DBL), parameter :: c_  = 3e8;
  real(kind=DBL), parameter :: u0  = u0_ / rE;
  real(kind=DBL), parameter :: e0  = e0_ * rE ** 3;
  real(kind=DBL), parameter :: c = c_ / rE;

end module constants

module model
	use constants
	implicit none 

	integer, parameter :: nU1 = 40
  integer, parameter :: nU2 = 42
  integer, parameter :: nU3 = 41
	integer, parameter :: maxK = 6
	integer, parameter :: maxM = 4
  
  real(kind=DBL) :: minCoLat = 20d0
  real(kind=DBL) :: maxCoLat = 70d0

end module model

module mhd_grid
  use constants
	use model
  implicit none
 
  real(kind=DBL) :: minU2
  real(kind=DBL) :: maxU2

  real(kind=DBL) :: dTheta0 
  real(kind=DBL), dimension ( nU1 ) :: theta0, theta0_e2
  real(kind=DBL) :: dU2
  
  ! u1, u2, u3 coords
  
  real(kind=DBL), dimension ( nU1 ) :: u1
  real(kind=DBL), dimension ( nU2 ) :: u2, u2_e1
  real(kind=DBL), dimension ( nU3 ) :: u3
	real(kind=DBL) :: u2_wrap ( nU2 + 2 ), u2_e1_wrap ( nU2 + 2 )

  real(kind=DBL), dimension ( nU1, nU2, nU3 ) :: r, theta, phi
  real(kind=DBL), dimension ( nU1, nU2, nU3 ) :: x, y, z

	!	Fit coordinates

	real(kind=DBL), dimension ( nU1, nU2 ) :: fitTheta0_h3, fitTheta0_e2, fitTheta0_e1, &
		fitPhi_h3, fitPhi_e2, fitPhi_e1

  contains
    subroutine make_grid
			implicit none 
      
      integer :: i, j, k
      real(kind=DBL) :: u3SpaceFac, A
      real(kind=DBL), dimension ( nU3 / 2 + 1 ) :: u3Y, u3_half
      real(kind=DBL) :: newtonsRule

      dTheta0 = ( maxCoLat - minCoLat ) / ( nU1 - 1 )
      theta0  = (/ ( i * dTheta0 + minCoLat, i = 0, nU1 - 1 ) /) * degToRad
      
      minU2 = 0.0
      maxU2 = 2.0 * pi - 2.0 * pi / nU2
      dU2 = ( maxU2 - minU2 ) / ( nU2 - 1.0 )
     
      u1 = cos ( theta0 ) ** 2 - 1.0
      u2 = (/ ( j * dU2 + minU2, j = 0, nU2 - 1 ) /)
			u2_wrap (2:nU2+1)	= u2
			u2_wrap (nU2+2)	= u2(1) + 2.0 * pi
  		u2_wrap (1)	= u2(nU2) - 2.0 * pi
  
      ! Design u3 with variable spacing according to a gaussian shape
    
      u3SpaceFac = 0.7! Lower means more points at the middle of the field line
      u3Y = (/ ( real ( i ) / ( nU3 / 2 + 1 ) + 1.0 / ( nU3 / 2 + 1 ), i = 0, nU3 / 2 ) /)
      A = -1.0 / log ( 1.0 / ( nU3 / 2 + 1 ) )
      u3_half = ( -A * log ( u3Y ) ) ** ( 1.0 / u3SpaceFac ) 
      u3(1:nU3 / 2 + 1) = -u3_half
      u3(nU3 / 2 + 2:)  = u3_half ( (/ ( i, i = nU3 / 2, 1, -1 ) /) )

      do i = 1, nU1
        do j = 1, nU2
          do k = 1, nU3

            r(i,j,k)  = newtonsRule ( rI, u1(i), u3(k), theta0(i), rI )
            theta(i,j,k)  = aCos ( -u3(k) * r(i,j,k) ** 2 * cos ( theta0(i) ) / rI ** 2 )
            phi(i,j,k)  = u2(j)

          end do
        end do
      end do
      
      x = r * cos ( phi ) * sin ( theta )
      y = r * sin ( theta ) * sin ( phi )
      z = r * cos ( theta )

			theta0_e2	= theta0 + dTheta0 / 2.0 * degToRad!	e2 theta0 coord
			u2_e1	= u2 + dU2 / 2.0!	e1 u2 coord
			u2_e1_wrap (2:nU2+1)	= u2_e1
			u2_e1_wrap (nU2+2)	= u2_e1(1) + 2.0 * pi
  		u2_e1_wrap (1)	= u2_e1(nU2) - 2.0 * pi
  

			do j = 1, nU2 
				fitTheta0_h3(:,j)	= theta0
				fitTheta0_e2(:,j)	= theta0_e2
				fitTheta0_e1(:,j)	= theta0
			end do
     
			do i = 1, nU1
				fitPhi_h3(i,:)	= u2
				fitPhi_e2(i,:)	= u2
				fitPhi_e1(i,:)	= u2_e1
			end do

    end subroutine make_grid

end module mhd_grid

module metric
  use constants
	use mhd_grid

	implicit none

  real(kind=DBL), dimension ( nU1, nU3 ) :: sqrtG, &
    g11cov, g12cov, g13cov, g21cov, g22cov, &
    g23cov, g31cov, g32cov, g33cov, g11con, & 
    g12con, g13con, g21con, g22con, g23con, &
    g31con, g32con, g33con
  real(kind=DBL) :: iono_sf_hTh ( nU1 ), &
	iono_sf_hR ( nU1 ), &
	iono_sf_hR_2D ( nU1, nU2 ), iono_sf_hTh_2D ( nU1, nU2 )

	contains
	subroutine make_metric
	
		integer :: i, k

		u1Loop: do i = 1, nU1 

      		u3Loop: do k = 1, nU3

			  sqrtG(i,k)  = ( r(i,1,k) ** 6 * sqrt ( ( 4.0 * r(i,1,k) - 2.0 * rI + &
			    2.0 * rI * cos ( 2.0 * theta(i,1,k) ) ) /  r(i,1,k) ) ) &
			    / ( rI ** 3 * ( 5.0 +  3 * cos ( 2.0 * theta(i,1,k) ) ) )

				!	Covariant metric tensor components
				!	See Lysak, JGR, v109, A07201, 2004
				!	or mathematica worksheet.
		
				g11cov(i,k)	= ( 4.0 * r(i,1,k) ** 4 * ( ( 4.0 * ( r(i,1,k) - rI ) ** 2 ) &
					/ ( 5.0 + 3.0 * cos ( 2.0 * theta(i,1,k) ) ) ** 2 + &
			           ( cos ( theta(i,1,k) ) ** 2 * ( 8* r(i,1,k) - 3.0 * rI + &
								 3.0 * rI * cos ( 2.0 * theta(i,1,k) ) ) ** 2 ) &
								 / ( 7.0 * sin ( theta(i,1,k) ) + 3.0 * sin ( 3.0 * theta(i,1,k) ) ) ** 2 ) ) / &
			       ( rI ** 2 * ( 2.0 * r(i,1,k) - rI + rI * cos ( 2.0 * theta(i,1,k) ) ) ** 2 )
				g12cov(i,k)	= 0.0
				g13cov(i,k)	= ( sqrt ( 2.0 ) * r(i,1,k) ** 4 * cos ( theta(i,1,k) ) ) / &
					( rI ** 2 * ( 5.0 + 3.0 * cos ( 2.0 * theta(i,1,k) ) ) * &
					sqrt ( ( 2.0 * r(i,1,k) - rI + rI * cos ( 2.0 * theta(i,1,k) ) ) / r(i,1,k) ) ) 
			
				g21cov(i,k)	= 0.0
				g22cov(i,k)	=  r(i,1,k) ** 2.0 * sin ( theta(i,1,k) ) ** 2
				g23cov(i,k)	= 0.0
			
				g31cov(i,k)	= g13cov(i,k)
				g32cov(i,k)	= 0.0
				g33cov(i,k)	= ( r(i,1,k) ** 5 * ( 2.0 * r(i,1,k) - rI + rI * cos ( 2.0 * theta(i,1,k) ) ) ) / &
					( rI ** 4 * ( 5.0 + 3.0 * cos ( 2.0 * theta(i,1,k) ) ) )
		
				!	Contravariant metric tensor components
		
				g11con(i,k)	= ( rI ** 2 * ( 5.0 + 3.0 * cos ( 2.0 * theta(i,1,k) ) ) &
					* sin( theta(i,1,k) ) ** 2 ) / ( 2.0 * r(i,1,k) ** 4 )
				g12con(i,k)	= 0.0
				g13con(i,k)	= ( rI ** 4 * cos ( theta(i,1,k) ) * sin ( theta(i,1,k) ) ** 2 * &
					( -4.0 + 3.0 * sin ( theta(i,1,k) ) ** 2 ) &
					* sqrt ( 1.0 - ( rI * sin ( theta(i,1,k) ) ** 2 ) / r(i,1,k) ) ) /&
			    ( 2d0 * r(i,1,k) ** 4 * ( r(i,1,k) - rI * sin ( theta(i,1,k) ) ** 2 ) ** 2 )
			
				g21con(i,k)	= 0.0
				g22con(i,k)	= 1.0 / ( sin ( theta(i,1,k) ) ** 2 *  r(i,1,k) ** 2 )
				g23con(i,k)	= 0.0
			
				g31con(i,k)	= g13con(i,k)
				g32con(i,k)	= 0.0
				g33con(i,k)	= ( rI ** 4 * ( 5.0 + 3.0 * cos ( 2.0 * theta(i,1,k) ) ) &
					* ( 32.0 * r(i,1,k) ** 2 - 32.0 * r(i,1,k) * rI + &
					13.0 * rI ** 2 +rI * ( -16.0 * ( -2.0 * r(i,1,k) + rI ) * cos ( 2.0 * theta(i,1,k) ) + &
					3.0 * rI * cos ( 4.0 * theta(i,1,k) ) ) ) ) &
					/ ( 64.0 * r(i,1,k) ** 5 * ( r(i,1,k) - rI * sin ( theta(i,1,k) ) ** 2 ) ** 3 )
		
		  end do u3Loop
		end do u1Loop

end subroutine make_metric

end module metric

