Program mhd3d
    use constants
    use mhd_grid
    use vA_profile
    use metric
    use spherHarmFns
    use imsl_libraries
    use dlg_timer
    use dislin
    use lin_sol_svd_int
    use linear_operators
    use netcdf
    use dlg

  implicit none

    !   Internal variables

    integer, parameter :: SGL_  = selected_real_kind ( p = 6, r = 37 )
    integer, parameter :: DBL_  = selected_real_kind ( p = 13, r = 200 )

    integer :: nBFns, status

    real(kind=DBL_), allocatable, dimension ( :,: ) :: YBFnArr_h3, &
        hrBFnArr_h3, hThBFnArr_h3, hPhBFnArr_h3, &
        hrBFnArr_e2, hThBFnArr_e2, hPhBFnArr_e2, &
        hrBFnArr_e1, hThBFnArr_e1, hPhBFnArr_e1, &
        hrBFnArr_h3_gnd, hThBFnArr_h3_gnd, hPhBFnArr_h3_gnd, &
        YBFnArr_h3_gnd, YBFnArr_e1, hrBFnArrT, YBFnArr_e2! Basis function arrays

    real(kind=DBL_) :: hrRealVector ( 1, nU1 * nU2 ), hrRealVector_ ( nU1 * nU2, 1 )

    real(kind=DBL_), dimension ( nU1, nU2 ) :: sigmaP_N, sigmaH_N, sigmaD_N, &
        sigmaP_S, sigmaH_S, sigmaD_S, alpha_N, sigmaZZ_N, &
        sigma11_N, sigma21_N, sigma12_N, sigma22_N, alpha_S, sigmaZZ_S, &
        sigma11_S, sigma21_S, sigma12_S, sigma22_S! Ionosphere conductance

    real(kind=DBL_), dimension ( nU1, nU2 ) :: sigmaP_N_e2, sigmaH_N_e2, sigmaD_N_e2, &
        sigmaP_S_e2, sigmaH_S_e2, sigmaD_S_e2, alpha_N_e2, sigmaZZ_N_e2, &
        sigma11_N_e2, sigma21_N_e2, sigma12_N_e2, sigma22_N_e2, alpha_S_e2, sigmaZZ_S_e2, &
        sigma11_S_e2, sigma21_S_e2, sigma12_S_e2, sigma22_S_e2! Ionosphere conductance

    real(kind=DBL_) :: dTArray ( nU1, nU2, nU3 ), &
        derivU1 ( nU1 ), derivU2 ( nU2 ), derivU3 ( nU3 ), &
        innerSum1, innerSum2, innerSum3, totalSum, dT ! dT
    
    integer :: i, j, k, ii, jj, kk
    integer, dimension ( 0:nU2+1 ) :: jjj ! Azimuthally wrapped index array

    real(kind=DBL_) :: u3_ ( nU3 - 1 ), &
        x_ ( nU1, nU2, nU3 - 1 ), &
        y_ ( nU1, nU2, nU3 - 1 ), &
        z_ ( nU1, nU2, nU3 - 1 ), &
        sqrtG_ ( nU1, nU3 - 1 ), &
        g11cov_ ( nU1, nU3 - 1 ), &
        g12cov_ ( nU1, nU3 - 1 ), &
        g13cov_ ( nU1, nU3 - 1 ), &
        g22cov_ ( nU1, nU3 - 1 ), &
        g21cov_ ( nU1, nU3 - 1 ), &
        g23cov_ ( nU1, nU3 - 1 ), &
        g31cov_ ( nU1, nU3 - 1 ), &
        g32cov_ ( nU1, nU3 - 1 ), &
        g33cov_ ( nU1, nU3 - 1 ), &
        g11con_ ( nU1, nU3 - 1 ), &
        g22con_ ( nU1, nU3 - 1 ), &
        g33con_ ( nU1, nU3 - 1 )    ! Interpolated metric components

    !   Metric components on the ionosphere

    real(kind=DBL_) :: g33covIono ( nU1, nU2 ), &
                    g11conIono ( nU1, nU2 ), &
                    g22conIono ( nU1, nU2 ), &
                    g33conIono ( nU1, nU2 )

    real(kind=DBL_) :: g33covIonoS ( nU1, nU2 ), &
                    g11conIonoS ( nU1, nU2 ), &
                    g22conIonoS ( nU1, nU2 ), &
                    g33conIonoS ( nU1, nU2 )

    real(kind=DBL_), allocatable :: h3conReal (:,:)

    real(kind=DBL_), allocatable :: alpha (:,:), beta_(:,:), &
        coeffs (:,:), coeffsOut(:,:), beta_S(:,:), coeffsS (:), &
        coeffsOutS(:,:), coeffsT (:,:), &
        beta_SVSOL (:,:), coeffs_ (:,:)

    real(kind=DBL_) :: u3MinRemoved ( nU3 ), u3MinRemoved_ ( nU3 - 1 )

    type ( timer ) :: startTime, startTimeOLD

    real(kind=DBL_) :: t, driverFreq, driverAmp, u3Variation, u2Variation, &
        u2DiffF ( nU2 ), u2DiffB ( nU2 ), timeTaken, eta

    real(kind=DBL_) :: h1con ( nU1, nU2, nU3 - 1 ), &
        h2con ( nU1, nU2, nU3 - 1 ), &
        h3con ( nU1, nU2, nU3 ), &
        h1cov ( nU1, nU2, nU3 - 1 ), &
        h2cov ( nU1, nU2, nU3 - 1 ), &
        h3cov ( nU1, nU2, nU3 ), &
        e1con ( nU1, nU2, nU3 ), &
        e2con ( nU1, nU2, nU3 ), &
        e3con ( nU1, nU2, nU3 - 1 ), &
        e1cov ( nU1, nU2, nU3 ), &
        e2cov ( nU1, nU2, nU3 ), &
        e3cov ( nU1, nU2, nU3 - 1 ), &
        j2cov ( nU1, nU2, nU3 )

    !   Fitted variables above and within the ionosphere

    integer :: nSingular

    real(kind(1d0)), parameter :: zero = 1.0d0
    type ( d_options ) :: iOpt_ ( 1 ) = d_options ( 0, zero )
    real(kind=DBL_), allocatable, dimension ( : ) :: S_
    real(kind(1d0)), parameter :: small = 1.0d-34

    real(kind=DBL_), allocatable :: A_ (:,:)

    real(kind=DBL_), dimension ( nU1, nU2 ) :: hrFit_h3, &
        hThFit_e1, hThFit_e2, hPhFit_e1, hPhFit_e2, &
        hThFit_h3, hPhFit_h3, YFit_h3, hThFit_h3_manual, &
        hPhFit_h3_manual, hThFit_e2_manual, hPhFit_e2_manual, &
        hThFit_e1_manual, hPhFit_e1_manual, YFit_e1, YFit_e2, &
        hrFit_h3_manual, hrFit_e1, hrFit_e2

    real(kind=DBL_) :: YFit_h3_wrap ( nU1, nU2+2 ), &
        YFit_e2_wrap ( nU1, nU2+2 ), YFit_e1_wrap ( nU1, nU2+2 )

    real(kind=DBL_), dimension ( nU1, nU2 ) :: hrFit_h3S, &
        hThFit_e1S, hThFit_e2S, hPhFit_e1S, hPhFit_e2S, &
        hThFit_h3S, hPhFit_h3S, YFit_h3S

    !   Variables for interpolation at the ionosphere 

    real(kind=DBL_) :: h1cov_h2 ( nU1, nU2 ), h2cov_h1 ( nU1, nU2 ), &
        t1, t2, t3, h1B1, h1B2, h2B1, h2B2, &
        dhPh_e1 ( nU1, nU2 ), dhTh_e1 ( nU1, nU2 ), &
        dhPh_e2 ( nU1, nU2 ), dhTh_e2 ( nU1, nU2 ), &
        JTh_e1 ( nU1, nU2 ), JPh_e1 ( nU1, nU2 ), &
        JTh_e2 ( nU1, nU2 ), JPh_e2 ( nU1, nU2 ), &
        E1covIono ( nU1, nU2 ), E2covIono ( nU1, nU2 ), &
        E1covIono_dU2 ( nU1, nU2 ), E2covIono_dU1 ( nU1, nU2 )

    real(kind=DBL_) :: h1cov_h2S ( nU1, nU2 ), h2cov_h1S ( nU1, nU2 ), &
        t1S, t2S, t3S, h1B1S, h1B2S, h2B1S, h2B2S, &
        dhPh_e1S ( nU1, nU2 ), dhTh_e1S ( nU1, nU2 ), &
        dhPh_e2S ( nU1, nU2 ), dhTh_e2S ( nU1, nU2 ), &
        JTh_e1S ( nU1, nU2 ), JPh_e1S ( nU1, nU2 ), &
        JTh_e2S ( nU1, nU2 ), JPh_e2S ( nU1, nU2 ), &
        E1covIonoS ( nU1, nU2 ), E2covIonoS ( nU1, nU2 )

    !   Plotting variables
    
    integer, parameter ::  nLevsECov=41, nhLevsAuto=21, neLevsAuto = 21
    integer :: nTri, nTri_, nTriN, nTriN_, plotSlice, testFn = 10
    integer :: i1Ray (3*nU1*nU3+1), i2Ray (3*nU1*nU3+1), &
        i3Ray (3*nU1*nU3+1), i1Ray_ (3*nU1*(nU3-1)+1), i2Ray_ (3*nU1*(nU3-1)+1), &
        i3Ray_ (3*nU1*(nU3-1)+1)
    integer :: i1RayN (3*nU1*nU3+1), i2RayN (3*nU1*nU3+1), &
        i3RayN (3*nU1*nU3+1), i1RayN_ (3*nU1*(nU3-1)+1), i2RayN_ (3*nU1*(nU3-1)+1), &
        i3RayN_ (3*nU1*(nU3-1)+1)
    real(kind=DBL_) :: plotFreq, progFreq, eCovLevSpacing, hLevsAuto ( nhLevsAuto ), &
        eLevsAuto ( neLevsAuto )
    real(kind=DBL_) :: xPlot (nU1*nU3+3), zPlot (nU1*nU3+3), &
        xPlot_ (nU1*(nU3-1)+3), zPlot_ (nU1*(nU3-1)+3), &
        eCovLevs (nLevsECov), eCovColors (nLevsECov)

    !   netcdf variables
    
    character(len=100) :: ncFileName
    integer :: nc_id, nU1_id, nU2_id, Y_id, theta0_id, phi_id, &
        hR_id, hTh_id, hPh_id, hRFit_id, hThFit_id, hPhFit_id, &
        hThFitMan_id, hPhFitMan_id

    write (*,*) DBL_

    call make_grid 
    call make_vA_profile
    call make_metric
    
    !   Find the number of basis functions and allocate the arrays
    !   Could avoid this if i found an analytical expression for the
    !   number of basis functions, just too lazy

    nBFns   = numberBFns ()

    write (*,*) nBFns
    !read (*,*)

    allocate ( YBFnArr_h3 ( nBFns, nU1 * nU2 ), &
                            hrBFnArr_h3 ( nBFns, nU1 * nU2 ), &
                            hThBFnArr_h3 ( nBFns, nU1 * nU2 ), &
                            hPhBFnArr_h3 ( nBFns, nU1 * nU2 ), &
                            hrBFnArr_e2 ( nBFns, nU1 * nU2 ), &
                            hThBFnArr_e2 ( nBFns, nU1 * nU2 ), &
                            hPhBFnArr_e2 ( nBFns, nU1 * nU2 ), &
                            hrBFnArr_e1 ( nBFns, nU1 * nU2 ), &
                            hThBFnArr_e1 ( nBFns, nU1 * nU2 ), &
                            hPhBFnArr_e1 ( nBFns, nU1 * nU2 ), &
                            YBFnArr_e1 ( nBFns, nU1 * nU2 ), &
                            YBFnArr_e2 ( nBFns, nU1 * nU2 ), &
                            hrBFnArr_h3_gnd ( nBFns, nU1 * nU2 ), &
                            hThBFnArr_h3_gnd ( nBFns, nU1 * nU2 ), &
                            hPhBFnArr_h3_gnd ( nBFns, nU1 * nU2 ), &
                            YBFnArr_h3_gnd ( nBFns, nU1 * nU2 ), &
                            S_ ( nBFns ), &
                            A_ ( nU1 * nU2, 0:nBFns-1 ), &
                            stat = status )
    
    write (*,*) 'Allocated arrays'
    
    !   Generate the SH fit basis function arrays

    call setupSHFns ( nU1 * nU2, nBFns,  rI, &
        reshape ( fitTheta0_h3, (/ nU1 * nU2 /) ), reshape ( fitPhi_h3, (/ nU1 * nU2 /) ), &
        hrBFnArr_h3, hThBFnArr_h3, hPhBFnArr_h3, YBFnArr_h3, &
        minCoLat, maxCoLat + dTheta0 / 2.0, dTheta0, dU2 )
    
    !read (*,*)
!   call metaFl ( 'XWIN' )
!   call disIni () 
!   call frame ( 1 ) ! Set plot frame thickness to 0
!   call setGrf ('NAME', 'NAME', 'TICKS', 'TICKS')
!   call noClip ()
!   call noChek ()
!   call incMrk (1)
!   !call labels ( 'FEXP', 'XYZ' )
!   call graf ( 10.0, 80.0, 0.0, 20.0, &
!       real ( -maxVal ( abs ( hPhBFnArr_h3(2,:) ) ) ), &
!       real ( maxVal ( abs ( hPhBFnArr_h3(2,:) ) )), &
!       real ( -maxVal ( abs ( hPhBFnArr_h3(2,:) ) )), &
!       real ( maxVal ( abs ( hPhBFnArr_h3(2,:) ) )) /  4.0 )
!
!   call curve ( real ( fitTheta0_h3 * radToDeg ), &
!           real ( hPhBFnArr_h3(2,:) ), nU1 * nU2 )
!   call endGrf ()  
!
!   write (*,*)
!   write ( *,* ) hrBFnArr_h3 ( 1,:)
!   write (*,*)
!   write ( *,* ) hThBFnArr_h3 ( 1,:)
!   write (*,*)
!   write (*,*) fitTheta0_h3 * radToDeg
!   write (*,*)
!   write (*,*) fitPhi_h3 * radToDeg

!   write ( *,* ) hPhBFnArr_h3 ( 1,: )
!   read ( *,* )
    
    call setupSHFns ( nU1 * nU2, nBFns,  rI, &
        reshape ( fitTheta0_e2, (/ nU1 * nU2 /) ), reshape ( fitPhi_e2, (/ nU1 * nU2 /) ), &
        hrBFnArr_e2, hThBFnArr_e2, hPhBFnArr_e2, YBFnArr_e2, &
        minCoLat, maxCoLat + dTheta0 / 2d0, dTheta0, dU2 )

    call setupSHFns ( nU1 * nU2, nBFns,  rI, &
        reshape ( fitTheta0_e1, (/ nU1 * nU2 /) ), reshape ( fitPhi_e1, (/ nU1 * nU2 /) ), &
        hrBFnArr_e1, hThBFnArr_e1, hPhBFnArr_e1, YBFnArr_e1, &
        minCoLat, maxCoLat + dTheta0 / 2d0, dTheta0, dU2 )

    call setupSHFns ( nU1 * nU2, nBFns,  1.d0, &
        reshape ( fitTheta0_h3, (/ nU1 * nU2 /) ), reshape ( fitPhi_h3, (/ nU1 * nU2 /) ), &
        hrBFnArr_h3_gnd, hThBFnArr_h3_gnd, hPhBFnArr_h3_gnd, YBFnArr_h3_gnd, &
        minCoLat, maxCoLat + dTheta0 / 2d0, dTheta0, dU2 )

    write (*,*) 'Basis functions created'

    sigmaP_N    = 1.0d0 * rE ** 2!
    sigmaH_N    = 0.0d0 * rE ** 2!
    sigmaD_N    = 1d6 * rE ** 2!

    sigmaP_S    = 1.0d0 * rE ** 2!
    sigmaH_S    = 2.5d0 * rE ** 2!
    sigmaD_S    = 1d6 * rE ** 2!

    sigmaP_N_e2 = 1.0d0 * rE ** 2!
    sigmaH_N_e2 = 0.0d0 * rE ** 2!
    sigmaD_N_e2 = 1d6 * rE ** 2!

    sigmaP_S_e2 = 1.0d0 * rE ** 2!
    sigmaH_S_e2 = 2.5d0 * rE ** 2!
    sigmaD_S_e2 = 1d6 * rE ** 2!

    alpha_N =   aCos ( -2.0 * cos ( fitTheta0_h3 ) / sqrt ( 1.0 + 3.0 * cos ( fitTheta0_h3 ) ** 2 ) )! 
    sigmaZZ_N   = sigmaD_N * cos ( alpha_N ) ** 2 + sigmaP_N * sin ( alpha_N ) ** 2!
    
    sigma11_N   = sigmaD_N * sigmaP_N / sigmaZZ_N!
    sigma21_N   = -sigmaD_N * sigmaH_N * cos ( alpha_N ) / sigmaZZ_N!
    sigma12_N   = sigmaD_N * sigmaH_N * cos ( alpha_N ) / sigmaZZ_N!
    sigma22_N   = sigmaP_N + sigmaH_N ** 2 * sin ( alpha_N ) ** 2 / sigmaZZ_N! 

    alpha_S =   aCos ( -2.0 * cos ( pi - fitTheta0_h3 ) / sqrt ( 1.0 + 3.0 * cos ( pi - fitTheta0_h3 ) ** 2 ) )! 
    sigmaZZ_S   = sigmaD_S * cos ( alpha_S ) ** 2 + sigmaP_S * sin ( alpha_S ) ** 2!
    
    sigma11_S   = sigmaD_S * sigmaP_S / sigmaZZ_S!
    sigma21_S   = -sigmaD_S * sigmaH_S * cos ( alpha_S ) / sigmaZZ_S!
    sigma12_S   = sigmaD_S * sigmaH_S * cos ( alpha_S ) / sigmaZZ_S!
    sigma22_S   = sigmaP_S + sigmaH_S ** 2 * sin ( alpha_S ) ** 2 / sigmaZZ_S! 

    alpha_N_e2  =   aCos ( -2.0 * cos ( fitTheta0_e2 ) / sqrt ( 1.0 + 3.0 * cos ( fitTheta0_e2 ) ** 2 ) )! 
    sigmaZZ_N_e2    = sigmaD_N * cos ( alpha_N_e2 ) ** 2 + sigmaP_N * sin ( alpha_N_e2 ) ** 2!
    
    sigma11_N_e2    = sigmaD_N_e2 * sigmaP_N_e2 / sigmaZZ_N_e2!
    sigma21_N_e2    = -sigmaD_N_e2 * sigmaH_N_e2 * cos ( alpha_N_e2 ) / sigmaZZ_N_e2!
    sigma12_N_e2    = sigmaD_N_e2 * sigmaH_N_e2 * cos ( alpha_N_e2 ) / sigmaZZ_N_e2!
    sigma22_N_e2    = sigmaP_N_e2 + sigmaH_N_e2 ** 2 * sin ( alpha_N_e2 ) ** 2 / sigmaZZ_N_e2! 

    alpha_S_e2  =   aCos ( -2.0 * cos ( pi - fitTheta0_e2 ) / sqrt ( 1.0 + 3.0 * cos ( pi - fitTheta0_e2 ) ** 2 ) )! 
    sigmaZZ_S_e2    = sigmaD_S_e2 * cos ( alpha_S_e2 ) ** 2 + sigmaP_S_e2 * sin ( alpha_S_e2 ) ** 2!
    
    sigma11_S_e2    = sigmaD_S_e2 * sigmaP_S_e2 / sigmaZZ_S_e2!
    sigma21_S_e2    = -sigmaD_S_e2 * sigmaH_S_e2 * cos ( alpha_S_e2 ) / sigmaZZ_S_e2!
    sigma12_S_e2    = sigmaD_S_e2 * sigmaH_S_e2 * cos ( alpha_S_e2 ) / sigmaZZ_S_e2!
    sigma22_S_e2    = sigmaP_S_e2 + sigmaH_S_e2 ** 2 * sin ( alpha_S_e2 ) ** 2 / sigmaZZ_S_e2! 



    write (*,*) 'Sigma arrays created'

    !   Calculate the maximum stable dT for this grid

    write (*,*) 'Calculating derivatives of coordinate variables ...'   

    do i = 1, nU1 
        derivU1(i)  = qdder ( 1, real ( i, kind=DBL_ ), (/ ( real ( j, kind=DBL_ ), j = 1, nU1 ) /), u1 )
    end do
    do j = 1, nU2 
        derivU2(j)  = qdder ( 1, real ( j, kind=DBL_ ), (/ ( real ( i, kind=DBL_ ), i = 1, nU2 ) /), u2 )
    end do
    do k = 1, nU3 
        derivU3(k)  = qdder ( 1, real ( k, kind=DBL_ ), (/ ( real ( j, kind=DBL_ ), j = 1, nU3 ) /), u3 )
    end do

    write (*,*) 'Calculating minimum dT value ... '

    do i = 1, nU1 
        do j = 1, nU2 
            do k = 1, nU3 
                
                innerSum1 = g11con(i,k) / ( derivU1(i) * derivU1(i) ) &
                    + g12con(i,k) / ( derivU1(i) * derivU2(j) ) &
                    + g13con(i,k)   / ( derivU1(i) * derivU3(k) )!
                innerSum2 = g21con(i,k) / ( derivU2(j) * derivU1(i) ) &
                    + g22con(i,k) / ( derivU2(j) * derivU2(j) ) &
                    + g23con(i,k)   / ( derivU2(j) * derivU3(k) )!
                innerSum3 = g31con(i,k) / ( derivU3(k) * derivU1(i) ) &
                    + g32con(i,k) / ( derivU3(k) * derivU2(j) ) &
                    + g33con(i,k)   / ( derivU3(k) * derivU3(k) )!
                
                totalSum    = innerSum1 + innerSum2 + innerSum3!
                dTArray(i,j,k)  = 1.0d0 / ( vA(i,k) * sqrt ( totalSum ) )!
    
            end do
        end do
    end do

    dT  = minVal ( dTArray ) * 0.9

    write (*,*) 'dT: ', dT

    !   Create metric components at interpolated grid points

    write (*,*) 'Calculating metric components at interpolated grid points ... '

    do j = 1, nU3 - 1 
        u3_(j)  = qdval ( j + 0.5d0, (/ ( real ( i, kind=DBL_ ), i = 1, nU3 ) /), u3 )
    end do
    
    do i = 1, nU1 
        do j = 1, nU2
            do k = 1, nU3 - 1

                x_(i,j,k)   = qd3vl ( real(i, kind=DBL_), real(j, kind=DBL_), &
                    k + 0.5d0, (/ (real(ii, kind=DBL_),ii=1,nU1) /), &
					(/ (real(jj, kind=DBL_),jj=1,nU2) /), (/ (real(kk, kind=DBL_),kk=1,nU3) /), x )
                y_(i,j,k)   = qd3vl ( real(i, kind=DBL_), real(j, kind=DBL_), &
                    k + 0.5d0, (/ (real(ii, kind=DBL_),ii=1,nU1) /), &
					(/ (real(jj, kind=DBL_),jj=1,nU2) /), (/ (real(kk, kind=DBL_),kk=1,nU3) /), y )
                z_(i,j,k)   = qd3vl ( real(i, kind=DBL_), real(j, kind=DBL_), &
                    k + 0.5d0, (/ (real(ii, kind=DBL_),ii=1,nU1) /), &
					(/ (real(jj, kind=DBL_),jj=1,nU2) /), (/ (real(kk, kind=DBL_),kk=1,nU3) /), z )

            end do
        end do
    end do

    do i = 1, nU1
        do k = 1, nU3 - 1

            sqrtG_(i,k) = qd2vl ( real(i, kind=DBL_), k+0.5d0, &
              (/ (real(ii, kind=DBL_),ii=1,nU1) /), (/ (real(kk, kind=DBL_),kk=1,nU3) /), sqrtG )
            g11cov_(i,k)    = qd2vl ( real(i, kind=DBL_), k+0.5d0, &
              (/ (real(ii, kind=DBL_),ii=1,nU1) /), (/ (real(kk, kind=DBL_),kk=1,nU3) /), g11cov )
            g12cov_(i,k)    = qd2vl ( real(i, kind=DBL_), k+0.5d0, &
              (/ (real(ii, kind=DBL_),ii=1,nU1) /), (/ (real(kk, kind=DBL_),kk=1,nU3) /), g12cov )
            g13cov_(i,k)    = qd2vl ( real(i, kind=DBL_), k+0.5d0, &
              (/ (real(ii, kind=DBL_),ii=1,nU1) /), (/ (real(kk, kind=DBL_),kk=1,nU3) /), g13cov )
            g22cov_(i,k)    = qd2vl ( real(i, kind=DBL_), k+0.5d0, &
              (/ (real(ii, kind=DBL_),ii=1,nU1) /), (/ (real(kk, kind=DBL_),kk=1,nU3) /), g22cov )
            g21cov_(i,k)    = qd2vl ( real(i, kind=DBL_), k+0.5d0, &
              (/ (real(ii, kind=DBL_),ii=1,nU1) /), (/ (real(kk, kind=DBL_),kk=1,nU3) /), g21cov )
            g23cov_(i,k)    = qd2vl ( real(i, kind=DBL_), k+0.5d0, &
              (/ (real(ii, kind=DBL_),ii=1,nU1) /), (/ (real(kk, kind=DBL_),kk=1,nU3) /), g23cov )
            g31cov_(i,k)    = qd2vl ( real(i, kind=DBL_), k+0.5d0, &
              (/ (real(ii, kind=DBL_),ii=1,nU1) /), (/ (real(kk, kind=DBL_),kk=1,nU3) /), g31cov )
            g32cov_(i,k)    = qd2vl ( real(i, kind=DBL_), k+0.5d0, &
              (/ (real(ii, kind=DBL_),ii=1,nU1) /), (/ (real(kk, kind=DBL_),kk=1,nU3) /), g32cov )
            g33cov_(i,k)    = qd2vl ( real(i, kind=DBL_), k+0.5d0, &
              (/ (real(ii, kind=DBL_),ii=1,nU1) /), (/ (real(kk, kind=DBL_),kk=1,nU3) /), g33cov )
            
            g11con_(i,k)    = qd2vl ( real(i, kind=DBL_), k+0.5d0, &
              (/ (real(ii, kind=DBL_),ii=1,nU1) /), (/ (real(kk, kind=DBL_),kk=1,nU3) /), g11con )
            g22con_(i,k)    = qd2vl ( real(i, kind=DBL_), k+0.5d0, &
              (/ (real(ii, kind=DBL_),ii=1,nU1) /), (/ (real(kk, kind=DBL_),kk=1,nU3) /), g22con )
            g33con_(i,k)    = qd2vl ( real(i, kind=DBL_), k+0.5d0, &
              (/ (real(ii, kind=DBL_),ii=1,nU1) /), (/ (real(kk, kind=DBL_),kk=1,nU3) /), g33con )
        
        end do
    end do

    !   Create ionosphere copies of the metric components

    write (*,*) 'Creating ionosphere copies of the metric components ...'
    
    do j = 1, nU2 
            
        g33covIono(:,j) = g33cov(:,1)
        g11conIono(:,j) = g11con(:,1)
        g22conIono(:,j) = g22con(:,1)
        g33conIono(:,j) = g33con(:,1)

        g33covIonoS(:,j)    = g33cov(:,nU3)
        g11conIonoS(:,j)    = g11con(:,nU3)
        g22conIonoS(:,j)    = g22con(:,nU3)
        g33conIonoS(:,j)    = g33con(:,nU3)

    end do

	!	Create the ionosphere scale factors

	iono_sf_hR	= - 2.0 * rI * cos ( theta0 )**2 / ( 1.0 + 3.0 * cos ( theta0 )**2 )
	iono_sf_hTh	= -rI / ( 2.0 * sin ( theta0 ) * cos ( theta0 ) )

	do j=1,nU2
		iono_sf_hR_2D(:,j) = iono_sf_hR
		iono_sf_hTh_2D(:,j) = iono_sf_hTh
	end do


    write (*,*) 'Allocating ionosphere fit variables ... ' 

    allocate ( hrBFnArrT ( nU1 * nU2, nBFns ), &
        alpha ( nBFns, nBFns ), &
        coeffs ( 1, nBFns ), &
        coeffsOut ( nBFns, 1 ), &
        coeffs_ ( 0:nBFns-1, 1 ), &
        coeffsOutS ( nBFns, 1 ), &
        coeffsT ( 1, nBFns ), &
        beta_ ( 1, nBFns ), &
        beta_S ( nBFns, 1 ), &
        h3conReal ( nU1, nU2 ), &
        beta_SVSOL ( nBFns, 1 ) )

    write (*,*) 'Creating svd alpha array ...'

    hrBFnArrT   = transpose ( hrBFnArr_h3 )
    alpha   = matMul ( hrBFnArr_h3, transpose ( hrBFnArr_h3 ) )

    !write ( *,* ) hrBFnArr_h3

    u3MinRemoved    = u3 - minVal ( u3 )
    u3MinRemoved_   = u3_ - minVal ( u3_ )

    !   Create the azimuthally wrapping index array

    write (*,*) 'Creating azimuthally wrapped index array ... ' 

    jjj(1:nU2)  = (/ (j,j=1,nU2) /)
    jjj(0)  = nU2
    jjj(nU2+1)  = 1

    ! Create the azimuthal coord differenece array, 
    !   i.e., u2Diff

    do j = 1, nU2

        u2DiffF(j)  = ( u2(j) - u2( jjj(j+1) ) )
        if ( j == nU2 ) u2DiffF(j) = u2DiffF(j) - 2.0 * pi

        u2DiffB(j)  = ( u2( jjj(j-1) ) - u2( j ) )
        if ( j == 1 ) u2DiffB(j) = u2DiffB(j) - 2.0 * pi

    end do

    !write (*,*) u2DiffB, u2DiffF, jjj

    !   BEGIN TIME LOOP

    call metaFl ( 'CONS' )
    call x11Mod ( 'STORE' )
    call winSiz ( 1200, 600 )
    call page ( 4800, 2400 )
    call disIni () 
    call winMod ( 'NONE' )
    call setVlt ( 'RAIN' )

    plotFreq    = 1d0
    progFreq    = 10d0

    write (*,*) 'STARTING TIME LOOP ... '

    timeLoop: &
    do t = 0.0d0, 1000.0d0, dT
        
        if ( mod ( t, progFreq  ) < dT ) then
        
            startTimeOLD    = startTime
            call start_timer ( startTime )
        
        end if
    
        !write (*,*) dT
    
        !   Apply h boundary conditions
    
        !   Driver
    
        driverFreq  = 50d-3!1./ (200.0 * dT )!
        driverAmp   = 10d-9 / u0!   
    
        if ( t <= ( 0.5 / driverFreq ) ) then
    
        !   print, 'DRIVING...'!
            
            do k = nU3/2-1, nU3/2+1
                do j = 1, nU2 
            
                    u3Variation = exp ( -u3(k) ** 2 / 0.000001 )!
                    u2Variation = cos ( 2.0 * u2 ( j ) ) !exp   ( -( u2(j) / ( 2 * pi ) - 0.5 ) ** 2 / 0.005 )! 
                    !h3cov(1,j,k)   = driverAmp * sin ( 2.0 * pi * driverFreq ) &
                    !   * u3Variation * u2Variation / sqrt ( g33con(1,k) )! 
                    h3cov(1,j,k)    = driverAmp * &
                        exp ( - ( t - 10.0 ) ** 2 / ( 2d0 * 10.0 ) ) &
                        * u3Variation * u2Variation / sqrt ( g33con(1,k) )! 

                end do
            end do
    
        else
        !   print, 'NOT DRIVING...'!
            h3cov(1,:,:)    = 0.0
        end if

        h2cov(1,:,:)    = 0.0!
        !h1cov(nU1,:,:) = ( 4.0 * h1cov(nU1-1,:,:) - h1cov(nU1-2,:,:) ) / 3.0! 

        !e2cov(1,:,1)   = ( 4.0 * e2cov(2,:,1) - e2cov(3,:,1) ) / 3.0!
        e2cov(nU1,:,2:nU3-1)    = 0.0!

    !   Update contravariant e field

    !   e1con except ionospheres
    
        do k = 2, nU3-1
            do j = 1, nU2
                do i = 1, nU1 
                
                        e1con(i,j,k)    = e1con(i,j,k) + &
                            dT / ( epsilon_(i,k) * sqrtG(i,k) ) * &
                            ( ( h3cov(i,j,k) - h3cov(i,jjj(j+1),k) ) / u2DiffF(j) &
                            - ( h2cov(i,j,k-1) - h2cov(i,j,k) ) / ( u3MinRemoved_(k-1) - u3MinRemoved_(k) ) )!
                    
                end do
            end do
        end do

        !   e2con except inner boundary and ionospheres

        do k = 2, nU3-1
            do j = 1, nU2
                do i = 1, nU1-1 

                    e2con(i,j,k)    = e2con(i,j,k) + &
                        dT / ( epsilon_(i,k) * sqrtG(i,k) ) * &
                        ( ( h1cov(i,j,k-1) - h1cov(i,j,k) ) / ( u3MinRemoved_(k-1) - u3MinRemoved_(k) ) &
                        - ( h3cov(i,j,k) - h3cov(i+1,j,k) ) / ( u1(i) - u1(i+1) ) )!
    
                end do
            end do
        end do

    !   Calculate covariant e from contravariant e

    !   Do not update the points at the ionosphere as they
    !   have already been filled by the BC
    
    !   e1cov except ionospheres
    
        do k = 2, nU3-1
            do j = 1, nU2
                do i = 1, nU1 

                    e1cov(i,j,k)    = 1.0 / g11con(i,k) * e1con(i,j,k)! &
                    
                    !   Do not update the points at the ionosphere as they
                    !   have already been filled by the BC
            
                end do
            end do
        end do
    
        !   e2cov except inner boundary and ionospheres

        do k = 2, nU3-1
            do j = 1, nU2
                do i = 1, nU1-1 

                        e2cov(i,j,k)    = g22cov(i,k) * e2con(i,j,k)!
    
                end do
            end do
        end do
        
        !   Update contravariant h field
    
        !   h1con except inner boundary

        do k = 1, nU3-1
            do j = 1, nU2
                do i = 1, nU1-1

                        h1con(i,j,k)    = h1con(i,j,k) - &
                            dT / ( u0 * sqrtG_(i,k) ) * &
                            (   -( e2cov(i,j,k) - e2cov(i,j,k+1) ) / ( u3MinRemoved(k) - u3MinRemoved(k+1) ) )!
    
                end do
            end do
        end do

        !   h2con all locations

        do k = 1, nU3-1
            do j = 1, nU2
                do i = 1, nU1
        
                    h2con(i,j,k)    = h2con(i,j,k) - &
                        dT / ( u0 * sqrtG_(i,k) ) * &
                        ( ( e1cov(i,j,k) - e1cov(i,j,k+1) ) / ( u3MinRemoved(k) - u3MinRemoved(k+1) ) )!&
                
                    end do
            end do
        end do
    
        !   h3con except outer boundary

        do k = 1, nU3
            do j = 1, nU2
                do i = 2, nU1
    
                        h3con(i,j,k)    = h3con(i,j,k) - &
                            dT / ( u0 * sqrtG(i,k) ) * &
                            ( ( e2cov(i-1,j,k) - e2cov(i,j,k) ) / ( u1(i-1) - u1(i) ) &
                            - ( e1cov(i,jjj(j-1),k) - e1cov(i,j,k) ) / u2DiffB(j)   )!
    
                end do
            end do
        end do

    !write (*,*)
        !write (*,*) sqrtG(:,1)
        !read (*,*)

        !   Set the dh3/du1=0 @ outer boundary in the two ionospheres
    
        h3con(1,:,:)    = ( 4.0 * h3con(2,:,:) - h3con(3,:,:) ) / 3.0!
        h3con(1,:,:)    = ( 4.0 * h3con(2,:,:) - h3con(3,:,:) ) / 3.0! 

        !write (*,*) u2 * radToDeg
        !write (*,*)
        !write (*,*) h3con(20,:,1)
        !write (*,*)
    
        !h3con(nU1,:,1) = ( 4.0 * h3con(nU1-1,:,1) - h3con(nU1-2,:,1) ) / 3.0!

        do k = 1, nU3
            do j = 1, nU2
                do i = 1, nU1

                    !   Calculate covariant h from contravariant h
            
                    !   h1cov except inner boundary
    
                    if ( i < nU1 .AND. k < nU3 ) &
                    h1cov(i,j,k)    = g11cov_(i,k) * h1con(i,j,k) &
                        + 0.25 * g13cov(i,k) * &
                            ( h3con(i,j,k) + h3con(i,j,k+1) + h3con(i+1,j,k) + h3con(i+1,j,k+1) )!
                    if ( i == nU1 .AND. k < nU3 ) &
                        h1cov(i,j,k)    = 0.0!
    
                    !   h2cov everywhere
    
                    if ( k < nU3 ) &
                    h2cov(i,j,k)    = g22cov_(i,k) * h2con(i,j,k)! 
    
                    !   h3cov except outer boundary and ionospheres
    
                    if ( i > 1 .AND. k > 1 .AND. k < nU3 ) &
                    h3cov(i,j,k)    = g33cov(i,k) * h3con(i,j,k) &
                        + 0.25 * g31cov(i,k) * &
                            ( h1con(i-1,j,k) + h1con(i-1,j,k-1) + h1con(i,j,k) + h1con(i,j,k-1) )! 
                    
                    !   h3cov at ionospheres except outer boundary
                    
                    if ( i > 1 .AND. k == nU3 ) &
                    h3cov(i,j,k)    = g33cov(i,k) * h3con(i,j,k) &
                        + 0.5 * g31cov(i,k) * &
                            ( h1con(i-1,j,k-1) + h1con(i,j,k-1) )! &
                    if ( i > 1 .AND. k == 1 ) &
                    h3cov(i,j,k)    = g33cov(i,k) * h3con(i,j,k) &
                        + 0.5 * g31cov(i,k) * &
                            ( h1con(i-1,j,k) + h1con(i,j,k) )! &
                
                end do
            end do
        end do

        h1cov(nU1,:,1)  = ( 4.0 * h1cov(nU1-1,:,1) - h1cov(nU1-2,:,1) ) / 3.0! 

    !   NORTHERN IONOSPHERE
    !   Reconstruct the horizontal b field component below 
    !   the ionosphere using a Spherical Harmonic fit.

    !hrRealVector    = reshape ( h3con(:,:,1) * sqrt ( g33covIono ), (/ 1, nU1 * nU2 /) )
    !hrRealVector_   = reshape ( h3con(:,:,1) * sqrt ( g33covIono ), (/ nU1 * nU2, 1 /) )
    hrRealVector    = reshape ( h3con(:,:,1) * iono_sf_hR_2D, (/ 1, nU1 * nU2 /) )
    hrRealVector_   = reshape ( h3con(:,:,1) * iono_sf_hR_2D, (/ nU1 * nU2, 1 /) )

    ! TEST CODE 
    !hrRealVector   = reshape ( hrBFnArr_h3(testFn,:), (/ 1, nU1 * nU2 /) )
    !hrRealVector_  = reshape ( hrBFnArr_h3(testFn,:), (/ nU1 * nU2, 1 /) )
    !hrRealVector    = reshape ( 2.0 * 8.0d22 * cos ( fitTheta0_h3 ) / ( rE ** 3 ) * 1d9, (/ 1, nU1 * nU2 /) )
    !hrRealVector_   = reshape ( 2.0 * 8.0d22 * cos ( fitTheta0_h3 ) / ( rE ** 3 ) * 1d9, (/ nU1 * nU2, 1 /) )


    beta_   = matMul ( hrRealVector, transpose ( hrBFnArr_h3 ) )

    beta_SVSOL = transpose ( beta_ )
    iOpt_(1)    = d_options ( d_lin_sol_svd_set_small, small )
    call lin_sol_svd ( alpha, beta_SVSOL, coeffsOut, rank = nSingular, S = S_, &
        iOpt = iOpt_)

    coeffs(1,:) = coeffsOut(:,1)

    !write (*,*) coeffs(1,:)
    !write (*,*)

    A_  = transpose ( hrBFnArr_h3 )
    !write (*,*) hrBFnArr_h3

    call lin_sol_lsq ( A_, hrRealVector_, coeffs_, nRHS = 1 )

    !write (*,*)
    !write (*,*) coeffs_
    !coeffs (1,:)   = coeffs_(:,1)

    !write (*,*) 'Non singular values: ', nSingular, nBFns
    !write (*,*) nBFns, nU1, nU2, nU1*nU2
    !read (*,*)
    !write (*,*)
    !write (*,*) h3con(:,2,1) * sqrt ( g33covIono(:,2) )
    !write (*,*)
    !write (*,*) beta_
    !write (*,*)
    !write (*,*) coeffs

    hrFit_h3    = reshape ( matMul ( coeffs, hrBfnArr_h3 ), (/ nU1, nU2 /) )
    hrFit_e1    = reshape ( matMul ( coeffs, hrBfnArr_e1 ), (/ nU1, nU2 /) )
    hrFit_e2    = reshape ( matMul ( coeffs, hrBfnArr_e2 ), (/ nU1, nU2 /) )

    hThFit_e1   = reshape ( matMul ( coeffs, hThBFnArr_e1 ), (/ nU1, nU2 /) ) 
    hThFit_e2   = reshape ( matMul ( coeffs, hThBFnArr_e2 ), (/ nU1, nU2 /) )
    hThFit_h3   = reshape ( matMul ( coeffs, hThBFnArr_h3 ), (/ nU1, nU2 /) )
    
    hPhFit_e1   = reshape ( matMul ( coeffs, hPhBFnArr_e1 ), (/ nU1, nU2 /) )
    hPhFit_e2   = reshape ( matMul ( coeffs, hPhBFnArr_e2 ), (/ nU1, nU2 /) )
    hPhFit_h3   = reshape ( matMul ( coeffs, hPhBFnArr_h3 ), (/ nU1, nU2 /) )
    
    YFit_h3 = reshape ( matMul ( coeffs, YBFnArr_h3 ), (/ nU1, nU2 /) ) 
    YFit_e1 = reshape ( matMul ( coeffs, YBFnArr_e1 ), (/ nU1, nU2 /) ) 
    YFit_e2 = reshape ( matMul ( coeffs, YBFnArr_e2 ), (/ nU1, nU2 /) ) 

    YFit_h3_wrap(:,2:nU2+1) = YFit_h3
    YFit_h3_wrap(:,1)   = YFit_h3(:,nU2)
    YFit_h3_wrap(:,nU2+2)   = YFit_h3(:,1)

    YFit_e2_wrap(:,2:nU2+1) = YFit_e2
    YFit_e2_wrap(:,1)   = YFit_e2(:,nU2)
    YFit_e2_wrap(:,nU2+2)   = YFit_e2(:,1)

    YFit_e1_wrap(:,2:nU2+1) = YFit_e1
    YFit_e1_wrap(:,1)   = YFit_e1(:,nU2)
    YFit_e1_wrap(:,nU2+2)   = YFit_e1(:,1)

    !   Try reconstructing the potential and then differentiate numberically.

    do j = 1, nU2
        do i = 1, nU1
            
                hThFit_h3_manual(i,j)   = &
                    1d0 / rI * qdDer ( 1, fitTheta0_h3 ( i,j ), fitTheta0_h3 (:,j), YFit_h3(:,j) )
                hPhFit_h3_manual(i,j)   = &
                    1d0 / ( rI * sin ( fitTheta0_h3 ( i,j ) ) ) * &
                    qdDer ( 1, u2 ( j ), u2_wrap, YFit_h3_wrap(i,:) )

                hThFit_e1_manual(i,j)   = &
                    1d0 / rI * qdDer ( 1, fitTheta0_e1 ( i,j ), fitTheta0_e1 (:,j), YFit_e1(:,j) )
                hPhFit_e1_manual(i,j)   = &
                    1d0 / ( rI * sin ( fitTheta0_e1 ( i,j ) ) ) * &
                    qdDer ( 1, u2_e1 ( j ), u2_e1_wrap, YFit_e1_wrap(i,:) )

                hThFit_e2_manual(i,j)   = &
                    1d0 / rI * qdDer ( 1, fitTheta0_e2 ( i,j ), fitTheta0_e2(:,j), YFit_e2(:,j) )
                hPhFit_e2_manual(i,j)   = &
                    1d0 / ( rI * sin ( fitTheta0_e2 ( i,j ) ) ) * &
                    qdDer ( 1, u2 ( j ), u2_wrap, YFit_e2_wrap(i,:) )

        end do
    end do

!   WRITE NETCDF FILE FOR DEBUGGING SH FITTING IN IDL

!        write(*,*) 'Writing sh_debug.nc ...'
!        ncFileName = 'sh_debug.nc'
!
!        call dlg_check ( nf90_create ( ncFileName, nf90_clobber, nc_id ) )
!        call dlg_check ( nf90_def_dim ( nc_id, "nU1", nU1, nU1_id ) )
!        call dlg_check ( nf90_def_dim ( nc_id, "nU2", nU2, nU2_id ) )
!
!        call dlg_check ( nf90_def_var ( nc_id, "Y", NF90_REAL, &
!            (/ nU1_id, nU2_id /), Y_id ) )
!        call dlg_check ( nf90_def_var ( nc_id, "theta0", NF90_REAL, &
!            (/ nU1_id, nU2_id /), theta0_id ) )
!        call dlg_check ( nf90_def_var ( nc_id, "phi", NF90_REAL, &
!            (/ nU1_id, nU2_id /), phi_id ) )
!        call dlg_check ( nf90_def_var ( nc_id, "hR", NF90_REAL, &
!            (/ nU1_id, nU2_id /), hR_id ) )
!        call dlg_check ( nf90_def_var ( nc_id, "hTh", NF90_REAL, &
!            (/ nU1_id, nU2_id /), hTh_id ) )
!        call dlg_check ( nf90_def_var ( nc_id, "hPh", NF90_REAL, &
!            (/ nU1_id, nU2_id /), hPh_id ) )
!        call dlg_check ( nf90_def_var ( nc_id, "hRFit", NF90_REAL, &
!            (/ nU1_id, nU2_id /), hRFit_id ) )
!        call dlg_check ( nf90_def_var ( nc_id, "hThFit", NF90_REAL, &
!            (/ nU1_id, nU2_id /), hThFit_id ) )
!        call dlg_check ( nf90_def_var ( nc_id, "hPhFit", NF90_REAL, &
!            (/ nU1_id, nU2_id /), hPhFit_id ) )
!        call dlg_check ( nf90_def_var ( nc_id, "hThFitMan", NF90_REAL, &
!            (/ nU1_id, nU2_id /), hThFitMan_id ) )
!        call dlg_check ( nf90_def_var ( nc_id, "hPhFitMan", NF90_REAL, &
!            (/ nU1_id, nU2_id /), hPhFitMan_id ) )
! 
!        call dlg_check ( nf90_enddef ( nc_id ) )
!
!        call dlg_check ( nf90_put_var ( nc_id, Y_id, &
!            reshape ( YBfnArr_h3(testFn,:), (/ nU1, nU2 /) )  ) )
!        call dlg_check ( nf90_put_var ( nc_id, hR_id, &
!            reshape ( hRBfnArr_h3(testFn,:), (/ nU1, nU2 /) )  ) )
!        call dlg_check ( nf90_put_var ( nc_id, hTh_id, &
!            reshape ( hThBfnArr_h3(testFn,:), (/ nU1, nU2 /) )  ) )
!        call dlg_check ( nf90_put_var ( nc_id, hPh_id, &
!            reshape ( hPhBfnArr_h3(testFn,:), (/ nU1, nU2 /) )  ) )
!        call dlg_check ( nf90_put_var ( nc_id, hRFit_id, &
!            hRFit_h3 ) )
!        call dlg_check ( nf90_put_var ( nc_id, hThFit_id, &
!            hThFit_h3 ) )
!        call dlg_check ( nf90_put_var ( nc_id, hPhFit_id, &
!            hPhFit_h3 ) )
!        call dlg_check ( nf90_put_var ( nc_id, hThFitMan_id, &
!            hThFit_h3_manual ) )
!         call dlg_check ( nf90_put_var ( nc_id, hPhFitMan_id, &
!            hPhFit_h3_manual ) )
!       
!        call dlg_check ( nf90_put_var ( nc_id, theta0_id, &
!            fitTheta0_h3 ) )
!        call dlg_check ( nf90_put_var ( nc_id, phi_id, &
!            fitPhi_h3 ) )
!
!        call dlg_check ( nf90_close ( nc_id ) )
!        write(*,*) 'DONE'
!        !read(*,*) 
!
!   hThFit_h3   = hThFit_h3_manual
!   hPhFit_h3   = hPhFit_h3_manual
!   hThFit_e1   = hThFit_e1_manual
!   hPhFit_e1   = hPhFit_e1_manual
!   hThFit_e2   = hThFit_e2_manual
!   hPhFit_e2   = hPhFit_e2_manual

!   write (*,*) fitTheta0_e1(:,1) * radToDeg
!   write (*,*) 
!   write (*,*) fitTheta0_h3(:,1) * radToDeg
    !read (*,*)
    !hThFit_e1  = hThFit_e1_manual
    !hThFit_e2  = hThFit_e2_manual
    !hPhFit_e1  = hPhFit_e1_manual
    !hPhFit_e2  = hPhFit_e2_manual
    !hPhFit_h3  = hPhFit_h3_manual

    !write (*,*)
    !write (*,*) hPhFit_h3_manual(:,2) * u0 * 1e9
    !write (*,*)    
    !write (*,*) hPhFit_h3(:,2) * u0 * 1e9

    !hThFit_h3  = -reshape ( hThBFnArr_h3 .x. transpose ( coeffs ), (/ nU1, nU2 /) );
    !hPhFit_h3  = reshape ( hPhBFnArr_h3 .x. transpose ( coeffs ), (/ nU1, nU2 /) ); 

    !   Calculate the eField boundary values from reconstructed b 
    !   and assumed sigma values. First need reconstructed h2 and h1
    !   at e1 and e1 locations. Here use a simple averaging ... but 
    !   could use the SH reconstruction to get them at those points too.

    !   JUST ABOVE THE IONOSPHERE
    !   Get h1 @ h2 locations


    do j = 1, nU2
        do i = 1, nU1

            jj  = j + 1
            if ( jj > nU2 ) jj = jj - nU2

            if ( i > 1 ) then 

                h1cov_h2(i,j)   = ( h1cov(i,j,1) &
                    + h1cov(i,jj,1) &
                    + h1cov(i-1,j,1) &
                    + h1cov(i-1,jj,1) ) / 4d0!
            
            else 
            
            !   @ i = 0 need to extrapolate h1 to h2 locations.
            !   Use Taylor expansion to get h1 above and below
            !   h2 location, then linear interpolate :)

                t1  = h1cov(1,j,1)!
                t2  = dTheta0 / 2d0 * ( -3d0 * h1cov(1,j,1) + 4d0 * h1cov(2,j,1) &
                    - h1cov(3,j,1) ) / ( 2d0 * dTheta0 )!
                t3  = ( dTheta0 / 2d0 ) ** 2 / 2d0 * ( 2d0 * h1cov(1,j,1) &
                    - 5d0 * h1cov(2,j,1) + 4d0 * h1cov(3,j,1) - h1cov(4,j,1) ) / dTheta0 ** 2!

                h1B1    = t1 + t2 + t3!
    
                t1  = h1cov(1,jj,1)!
                t2  = dTheta0 / 2d0 * ( -3d0 * h1cov(1,jj,1) &
                    + 4d0 * h1cov(2,jj,1) &
                    - h1cov(3,jj,1) ) / ( 2d0 * dTheta0 )!
                t3  = ( dTheta0 / 2d0 ) ** 2 / 2d0 * ( 2d0 * h1cov(1,jj,1) &
                    - 5d0 * h1cov(2,jj,1) &
                    + 4d0 * h1cov(3,jj,1) &
                    - h1cov(4,jj,1) ) / dTheta0 ** 2!

                h1B2    = t1 + t2 + t3!
            
                h1cov_h2(1,j)   = ( h1B1 + h1B2 ) / 2d0!
            
            end if

        end do
    end do
    
    !   Get h2 @ h1 locations above the ionosphere

    do j = 1, nU2
        do i = 1, nU1

            jj  = j - 1
            if ( jj < 1 ) jj = jj + nU2

            if ( i < nU1 ) then

                h2cov_h1(i,j)   = ( h2cov(i,j,1) &
                    + h2cov(i,jj,1) &
                    + h2cov(i+1,j,1) &
                    + h2cov(i+1,jj,1) ) / 4d0!
            
            else
            
            !   @ i = nU1 - 1 need to extrapolate h2 to h1 locations.
            !   Use Taylor expansion to get h2 above and below
            !   h1 location, then linear interpolate :)

                t1  = h2cov(i,j,1)!
                t2  = dTheta0 / 2d0 * ( 3d0 * h2cov(i,j,1) - 4d0 * h2cov(i-1,j,1) &
                    + h2cov(i-2,j,1) ) / ( 2d0 * dTheta0 )!
                t3  = ( dTheta0 / 2d0 ) ** 2 / 2d0 * ( 2d0 * h2cov(i,j,1) &
                    - 5d0 * h2cov(i-1,j,1) + 4d0 * h2cov(i-2,j,1) - h2cov(i-3,j,1) ) / dTheta0 ** 2!

                h2B1    = t1 + t2 + t3!
    
                t1  = h2cov(i,jj,1)!
                t2  = dTheta0 / 2d0 * ( -3d0 * h2cov(i,jj,1) &
                    + 4d0 * h2cov(i-1,jj,1) &
                    - h2cov(i-2,jj,1) ) / ( 2d0 * dTheta0 )!
                t3  = ( dTheta0 / 2d0 ) ** 2 / 2d0 * ( 2d0 * h2cov(i,jj,1) &
                    - 5d0 * h2cov(i-1,jj,1) &
                    + 4d0 * h2cov(i-2,jj,1) &
                    - h2cov(i-3,jj,1) ) / dTheta0 ** 2!

                h2B2    = t1 + t2 + t3!
            
                h2cov_h1(i,j)   = ( h2B1 + h2B2 ) / 2d0!
            
            end if

        end do
    end do
    
    dhPh_e1 = (h2cov(:,:,1) * sqrt ( g22conIono ) - hPhFit_e1)!
    !dhTh_e1 = (h1cov_h2 * sqrt ( g11conIono ) - hThFit_e1)!
	dhTh_e1 = (h1cov_h2 / iono_sf_hTh_2D - hThFit_e1)!

    dhPh_e2 = (h2cov_h1 * sqrt ( g22conIono ) - hPhFit_e2)!
    !dhTh_e2 = (h1cov(:,:,1) * sqrt ( g11conIono ) - hThFit_e2)!
	dhTh_e2 = (h1cov(:,:,1) / iono_sf_hTh_2D - hThFit_e2)!

!   JTh_e1  = -dhPh_e1!
!   JPh_e1  = dhTh_e1!
!
!   JTh_e2  = -dhPh_e2!
!   JPh_e2  = dhTh_e2!

    E1covIono   = ( dhTh_e1 * sigma21_N + dhPh_e1 * sigma22_N ) / &
                ( sigma12_N * sigma21_N - sigma11_N * sigma22_N )!

    E2covIono   = ( dhTh_e2 * sigma11_N_e2 + dhPh_e2 * sigma12_N_e2 ) / &
                ( -sigma12_N_e2 * sigma21_N_e2 + sigma11_N_e2 * sigma22_N_e2 )!

!   E1covIono   = ( JPh_e1 * sigma21_N - JTh_e1 * sigma22_N ) / &
!               ( sigma12_N * sigma21_N - sigma11_N * sigma22_N )!
!
!   E2covIono   = ( JPh_e2 * sigma11_N - JTh_e2 * sigma12_N ) / &
!               ( -sigma12_N * sigma21_N + sigma11_N * sigma22_N )!

!   NEED TO SET SIGMAS SLIGHTLY DIFFERENT do GRID OFFSET!

    !write (*,*)
    !write (*,*) ( E1covIono(:,2) / sqrt ( g11conIono(:,2) ) ) / e1cov(:,2,1)
    
    !e1cov(:,:,1)    = E1covIono / sqrt ( g11conIono )! 
    !e2cov(:,:,1)    = E2covIono / sqrt ( g22conIono )!
    e1cov(:,:,1)    = E1covIono * iono_sf_hTh_2D 
    e2cov(:,:,1)    = E2covIono / sqrt ( g22conIono )!
    
    !if ( t < 200 ) then

    !   e1cov(:,:,1)    = 0d0 
    !   e2cov(:,:,1)    = 0d0

    !else

    !   plotFreq    = dT

    !end if 

    !e2cov(nU1,:,:) = 0.0!

    !   These prints should be real values in V/m and A/m

    !write (*,*) 
    !write (*,*) sqrt ( ( E2covIono(:,2) * rE ) ** 2 + ( E1covIono(:,2) * rE ) ** 2 ) * 1e3
    !write (*,*) sqrt ( ( JTh_e1(:,2) / rE ) ** 2 + ( JPh_e1(:,2) / rE ) ** 2 )
    !read (*,*)

!   !   SOUTHERN IONOSPHERE
!   !   Reconstruct the horizontal b field component below 
!   !   the ionosphere using a Spherical Harmonic fit.
!
!   beta_S  = matMul ( transpose ( hrBFnArr_h3 ), reshape ( h3con(:,:,nU3) * sqrt ( g33covIonoS ), (/ nU1 * nU2, 1 /) ) )
!   
!   call lin_sol_svd ( alpha, beta_S, coeffsOutS )
!
!   coeffsS = coeffsOutS(:,1)
!   
!   hrFit_h3S   = reshape ( matMul ( coeffsS, hrBfnArr_h3 ), (/ nU1, nU2 /) )
!   hThFit_e1S  = reshape ( matMul ( coeffsS, hThBFnArr_e1 ), (/ nU1, nU2 /) ) 
!   hThFit_e2S  = reshape ( matMul ( coeffsS, hThBFnArr_e2 ), (/ nU1, nU2 /) )
!   hPhFit_e1S  = reshape ( matMul ( coeffsS, hPhBFnArr_e1 ), (/ nU1, nU2 /) ) 
!   hPhFit_e2S  = reshape ( matMul ( coeffsS, hPhBFnArr_e2 ), (/ nU1, nU2 /) ) 
!
!   !hThFit_h3  = -reshape ( hThBFnArr_h3 .x. transpose ( coeffs ), (/ nU1, nU2 /) );
!   !hPhFit_h3  = reshape ( hPhBFnArr_h3 .x. transpose ( coeffs ), (/ nU1, nU2 /) ); 
!   !YFit_h3    = reshape ( YBFnArr_h3 .x. transpose ( coeffs ), (/ nU1, nU2 /) ); 
!
!   !   Calculate the eField boundary values from reconstructed b 
!   !   and assumed sigma values. First need reconstructed h2 and h1
!   !   at e1 and e1 locations. Here use a simple averaging ... but 
!   !   could use the SH reconstruction to get them at those points too.
!
!   !   JUST ABOVE THE IONOSPHERE
!   !   Get h1 @ h2 locations
!
!
!   do j = 1, nU2
!       do i = 1, nU1
!
!           jj  = j + 1
!           if ( jj > nU2 ) jj = jj - nU2
!
!           if ( i > 1 ) then 
!
!               h1cov_h2S(i,j)  = ( h1cov(i,j,nU3-1) &
!                   + h1cov(i,jj,nU3-1) &
!                   + h1cov(i-1,j,nU3-1) &
!                   + h1cov(i-1,jj,nU3-1) ) / 4.0!
!           
!           else 
!           
!           !   @ i = 0 need to extrapolate h1 to h2 locations.
!           !   Use Taylor expansion to get h1 above and below
!           !   h2 location, then linear interpolate :)
!
!               t1S = h1cov(1,j,nU3-1)!
!               t2S = dTheta0 / 2.0 * ( -3.0 * h1cov(1,j,nU3-1) + 4.0 * h1cov(2,j,nU3-1) &
!                   - h1cov(3,j,nU3-1) ) / ( 2.0 * dTheta0 )!
!               t3S = ( dTheta0 / 2.0 ) ** 2 / 2.0 * ( 2.0 * h1cov(1,j,nU3-1) &
!                   - 5.0 * h1cov(2,j,nU3-1) + 4.0 * h1cov(3,j,nU3-1) - h1cov(4,j,nU3-1) ) / dTheta0 ** 2!
!
!               h1B1S   = t1S + t2S + t3S!
!   
!               t1S = h1cov(1,jj,nU3-1)!
!               t2S = dTheta0 / 2.0 * ( -3.0 * h1cov(1,jj,nU3-1) &
!                   + 4.0 * h1cov(2,jj,nU3-1) &
!                   - h1cov(3,jj,nU3-1) ) / ( 2.0 * dTheta0 )!
!               t3S = ( dTheta0 / 2.0 ) ** 2 / 2.0 * ( 2.0 * h1cov(1,jj,nU3-1) &
!                   - 5.0 * h1cov(2,jj,nU3-1) &
!                   + 4.0 * h1cov(3,jj,nU3-1) &
!                   - h1cov(4,jj,nU3-1) ) / dTheta0 ** 2!
!
!               h1B2S   = t1S + t2S + t3S!
!           
!               h1cov_h2S(1,j)  = ( h1B1S + h1B2S ) / 2.0!
!           
!           end if
!
!       end do
!   end do
!   
!   !   Get h2 @ h1 locations above the ionosphere
!
!   do j = 1, nU2
!       do i = 1, nU1
!
!           jj  = j - 1
!           if ( jj < 1 ) jj = jj + nU2
!
!           if ( i < nU1 ) then
!
!               h2cov_h1S(i,j)  = ( h2cov(i,j,nU3-1) &
!                   + h2cov(i,jj,nU3-1) &
!                   + h2cov(i+1,j,nU3-1) &
!                   + h2cov(i+1,jj,nU3-1) ) / 4.0!
!           
!           else
!           
!           !   @ i = nU1 - 1 need to extrapolate h2 to h1 locations.
!           !   Use Taylor expansion to get h2 above and below
!           !   h1 location, then linear interpolate :)
!
!               t1S = h2cov(i,j,nU3-1)!
!               t2S = dTheta0 / 2.0 * ( 3.0 * h2cov(i,j,nU3-1) - 4.0 * h2cov(i-1,j,nU3-1) &
!                   + h2cov(i-2,j,nU3-1) ) / ( 2.0 * dTheta0 )!
!               t3S = ( dTheta0 / 2.0 ) ** 2 / 2.0 * ( 2.0 * h2cov(i,j,nU3-1) &
!                   - 5.0 * h2cov(i-1,j,nU3-1) + 4.0 * h2cov(i-2,j,nU3-1) - h2cov(i-3,j,nU3-1) ) / dTheta0 ** 2!
!
!               h2B1S   = t1S + t2S + t3S!
!   
!               t1S = h2cov(i,jj,nU3-1)!
!               t2S = dTheta0 / 2.0 * ( -3.0 * h2cov(i,jj,nU3-1) &
!                   + 4.0 * h2cov(i-1,jj,nU3-1) &
!                   - h2cov(i-2,jj,nU3-1) ) / ( 2.0 * dTheta0 )!
!               t3S = ( dTheta0 / 2.0 ) ** 2 / 2.0 * ( 2.0 * h2cov(i,jj,nU3-1) &
!                   - 5.0 * h2cov(i-1,jj,nU3-1) &
!                   + 4.0 * h2cov(i-2,jj,nU3-1) &
!                   - h2cov(i-3,jj,nU3-1) ) / dTheta0 ** 2!
!
!               h2B2S   = t1S + t2S + t3S!
!           
!               h2cov_h1S(i,j)  = ( h2B1S + h2B2S ) / 2.0!
!           
!           end if
!
!       end do
!   end do
!   
!   dhPh_e1S    = (h2cov(:,:,nU3-1) * sqrt ( g22conIonoS ) - hPhFit_e1S)!
!   dhTh_e1S    = (h1cov_h2S * sqrt ( g11conIonoS ) - hThFit_e1S)!
!
!   dhPh_e2S    = (h2cov_h1S * sqrt ( g22conIonoS ) - hPhFit_e2S)!
!   dhTh_e2S    = (h1cov(:,:,nU3-1) * sqrt ( g11conIonoS ) - hThFit_e2S)!
!
!   JTh_e1S = -dhPh_e1S!
!   JPh_e1S = dhTh_e1S!
!
!   JTh_e2S = -dhPh_e2S!
!   JPh_e2S = dhTh_e2S!
!
!   E1covIonoS  = ( dhTh_e1S * sigma21_S + dhPh_e1S * sigma22_S ) / &
!               ( sigma12_S * sigma21_S - sigma11_S * sigma22_S )!
!
!   E2covIonoS  = ( dhTh_e2S * sigma11_S + dhPh_e2S * sigma12_S ) / &
!               ( -sigma12_S * sigma21_S + sigma11_S * sigma22_S )!
!
!   !write ( *,* ) E1covIono
!   !write ( *,* ) E2covIono
!
!   !e1cov(:,:,nU3) = E1covIonoS / sqrt ( g11conIonoS )! 
!   !e2cov(:,:,nU3) = E2covIonoS / sqrt ( g22conIonoS )!
        
        if ( mod ( t, plotFreq ) < dT ) then  

            write (*,*) t*dT
        
            plotSlice = nU2 / 2 + 1 
            
            xPlot(1:nU1*nU3)    = reshape ( sqrt ( x(:,plotSlice,:)**2 + y(:,plotSlice,:)**2 ), (/ nU1 * nU3 /) )
            zPlot(1:nU1*nU3)    = reshape ( z ( :, plotSlice, : ), (/ nU1 * nU3 /) )
        
            call triang ( real ( xPlot ), real ( zPlot ), &
                    nU1 * nU3, i1Ray, i2Ray, i3Ray, 3 * nU1 * nU3 + 1, nTri )
            call triang ( real ( -xPlot ), real ( zPlot ), &
                    nU1 * nU3, i1RayN, i2RayN, i3RayN, 3 * nU1 * nU3 + 1, nTriN )
            
            xPlot_(1:nU1*(nU3-1))   = reshape ( sqrt ( x_(:,plotSlice,:)**2 + y_(:,plotSlice,:)**2 ), (/ nU1 * (nU3-1) /) )
            zPlot_(1:nU1*(nU3-1))   = reshape ( z_( :, plotSlice, : ), (/ nU1 * (nU3-1) /) )
    
            call triang ( real ( xPlot_ ), real ( zPlot_ ), &
                    nU1 * (nU3-1), i1Ray_, i2Ray_, i3Ray_, 3 * nU1 * (nU3-1) + 1, nTri_ )
            call triang ( real ( -xPlot_ ), real ( zPlot_ ), &
                    nU1 * (nU3-1), i1RayN_, i2RayN_, i3RayN_, 3 * nU1 * (nU3-1) + 1, nTriN_ )
                    
            call erase ()
            call height ( 18 )
            call pageRa ()
            call pagFll ( 255 ) ! Set page background color
            call frame ( 0 ) ! Set plot frame thickness to 0
            call setGrf ( 'NONE', 'NONE', 'NONE', 'NONE' )

            eCovLevSpacing  = 0.005 ! [mVm^-1]
            eCovLevs    = ( (/ ( real(i), i=0,nLevsECov-1 ) /) - real ( int ( nLevsECov / 2 ) ) ) &
                * eCovLevSpacing
            eCovColors  = ( eCovLevs + abs ( minVal ( eCovLevs ) ) ) &
                / maxVal ( eCovLevs + abs ( minVal ( eCovLevs ) ) ) &
                * 253 + 1.0
        
            call axsPos ( 50, 800 )
            call axsLen ( 1500, 800 )
            call graf ( -9.0, 9.0, -9.0, 1.0, -4.0, 4.0, -4.0, 1.0 )
            call conFll ( real ( xPlot(1:nU1*nU3) ), real ( zPlot(1:nU1*nU3) ), &
                real ( reshape ( e1cov (:,plotSlice,:) * sqrt ( g11con ) * rE * 1e3, (/ nU1 * nU3 /) ) ), &
                nU1 * nU3, i1Ray, i2Ray, i3Ray, nTri, real ( eCovLevs ), nLevsECov )
            call conFll ( real ( -xPlot(1:nU1*nU3) ), real ( zPlot(1:nU1*nU3) ), &
                real ( reshape ( e1cov (:,1,:) * sqrt ( g11con ) * rE * 1e3, (/ nU1 * nU3 /) ) ), &
                nU1 * nU3, i1RayN, i2RayN, i3RayN, nTriN, real ( eCovLevs ), nLevsECov )
            call endGrf ()

            call axsPos ( 1650, 800 )
            call axsLen ( 1500, 800 )
            call graf ( -9.0, 9.0, -9.0, 1.0, -4.0, 4.0, -4.0, 1.0 )
            call conFll ( real ( xPlot(1:nU1*nU3) ), real ( zPlot(1:nU1*nU3) ), &
                real ( reshape ( e2cov (:,plotSlice,:) * sqrt ( g22con ) * rE * 1e3, (/ nU1 * nU3 /) ) ), &
                nU1 * nU3, i1Ray, i2Ray, i3Ray, nTri, real ( eCovLevs ), nLevsECov )
            call conFll ( real ( -xPlot(1:nU1*nU3) ), real ( zPlot(1:nU1*nU3) ), &
                real ( reshape ( e2cov (:,1,:) * sqrt ( g22con ) * rE * 1e3, (/ nU1 * nU3 /) ) ), &
                nU1 * nU3, i1RayN, i2RayN, i3RayN, nTriN, real ( eCovLevs ), nLevsECov )
            call endGrf ()

            call axsPos ( 3250, 800 )
            call axsLen ( 1500, 800 )
        !   call graf ( -9.0, 9.0, -9.0, 1.0, -4.0, 4.0, -4.0, 1.0 )
        !   call conFll ( real ( xPlot_(1:nU1*(nU3-1)) ), real ( zPlot_(1:nU1*(nU3-1)) ), &
        !       real ( reshape ( e3cov (:,plotSlice,:) * sqrt ( g33con_ ) * rE * 1e3, (/ nU1 * (nU3-1) /) ) ), &
        !       nU1 * (nU3-1), i1Ray_, i2Ray_, i3Ray_, nTri_, real ( eCovLevs ), nLevsECov )
        !   call conFll ( real ( -xPlot_(1:nU1*(nU3-1)) ), real ( zPlot_(1:nU1*(nU3-1)) ), &
        !       real ( reshape ( e3cov (:,1,:) * sqrt ( g33con_ ) * rE * 1e3, (/ nU1 * (nU3-1) /) ) ), &
        !       nU1 * (nU3-1), i1RayN_, i2RayN_, i3RayN_, nTriN_, real ( eCovLevs ), nLevsECov )
        !   call graf ( -80.0, 80.0, -80.0, 10.0, -80.0, 80.0, -80.0, 10.0 )
            do i = 1, neLevsAuto
                eLevsAuto ( i ) = -maxVal ( abs ( E1covIono * rE * 1e3 ) ) / 2.0 + &
                    i * maxVal ( abs ( E1covIono * rE * 1e3 ) ) / &
                    neLevsAuto
            end do
        
            do i = 1, nhLevsAuto
                hLevsAuto ( i ) = -maxVal ( abs ( h3con(:,:,1) * sqrt ( g33covIono ) * u0 * 1d9 ) ) / 2.0 + &
                    i * maxVal ( abs ( h3con(:,:,1) * sqrt ( g33covIono ) * u0 * 1d9 ) ) / &
                    nhLevsAuto
                !hLevsAuto ( i ) = -0.1 + i * 0.2 / nhLevsAuto
            end do
            !hLevsAuto  = hLevsAuto * 10000.0
!           do i = 1, nhLevsAuto
!               hLevsAuto ( i ) = -maxVal ( abs ( hrBFnArr_h3(testFn,:) * u0 * 1d9 ) ) + &
!                   i * 2.0 * maxVal ( abs ( hrBFnArr_h3(testFn,:) * u0 * 1d9 ) ) / &
!                   nhLevsAuto
!           end do
            do i = 1, nhLevsAuto
                hLevsAuto ( i ) = -maxVal ( abs ( 2.0 * 8.0d22 * cos ( fitTheta0_h3) / ( rE ** 3 ) * 1d9 ) ) + &
                    i * 2.0 * maxVal ( abs ( 2.0 * 8.0d22 * cos ( fitTheta0_h3 ) / ( rE ** 3 ) * 1d9 ) ) / &
                    nhLevsAuto
            end do


            !write (*,*) hLevsAuto
            !write (*,*)
            !write (*,*) h3con(:,:,1) * g33covIono * u0 * 1d9 
            !read (*,*)

            !call polar ( 90., 0., 10., 0., 60. )   
            call projct ( 'LAMBERT' )
            call mapPol ( 0.0, 90.0 )
            call grafMp ( -180.0, 180.0, -180.0, 30.0, -90.0, 90.0, -90.0, 30.0 )
!           call conShd ( real ( u2 * radToDeg ), nU2, &
!               real ( 90.0 - theta(:,1,1) * radToDeg ), nU1, &
!               real ( transpose ( h3con(:,:,1) * sqrt ( g33covIono ) * u0 * 1d9 ) ), &
!               real ( hLevsAuto ), nhLevsAuto )
!           call conShd ( real ( phi(1,:,1) * radToDeg ), nU2, &
!               real ( 90.0 - fitTheta0_h3(:,1) * radToDeg ), nU1, &
!               real ( transpose ( reshape ( hrBFnArr_h3(testFn,:), (/ nU1, nU2 /) ) * u0 * 1d9 ) ), &
!               real ( hLevsAuto ), nhLevsAuto )
            call conShd ( real ( phi(1,:,1) * radToDeg ), nU2, &
                real ( 90.0 - fitTheta0_h3(:,1) * radToDeg ), nU1, &
                real ( transpose ( reshape ( 2.0 * 8.0d22 * cos ( fitTheta0_h3 ) / ( rE ** 3 ) * 1d9, (/ nU1, nU2 /) ) ) ), &
                real ( hLevsAuto ), nhLevsAuto )
    
            call color ( 'BLACK' )
            call gridMp ( 1, 1)
            call endGrf ()

            call axsPos ( 3250, 1600 )
            call axsLen ( 1500, 800 )
            call projct ( 'LAMBERT' )
            call mapPol ( 0.0, 90.0 )
            call grafMp ( -180.0, 180.0, -180.0, 30.0, -90.0, 90.0, -90.0, 30.0 )
            call conShd ( real ( phi(1,:,1) * radToDeg ), nU2, &
                real ( 90.0 - fitTheta0_h3(:,1) * radToDeg ), nU1, &
                !real ( transpose ( hrFit_h3 * u0 * 1d9 ) ), &
                real ( transpose ( hrFit_h3 ) ), &
                real ( hLevsAuto ), nhLevsAuto )
            call color ( 'BLACK' )
            call gridMp ( 1, 1)
            call endGrf ()
            
            call axsPos ( 3250, 2399 )
            call axsLen ( 1500, 800 )
            call projct ( 'LAMBERT' )
            call mapPol ( 0.0, 90.0 )
            call grafMp ( -180.0, 180.0, -180.0, 30.0, -90.0, 90.0, -90.0, 30.0 )
            call conShd ( real ( phi(1,:,1) * radToDeg ), nU2, &
                real ( 90.0 - fitTheta0_e2(:,1) * radToDeg ), nU1, &
                !real ( transpose ( hrFit_e2 * u0 * 1d9 ) ), &
                real ( transpose ( hrFit_e2 ) ), &
                real ( hLevsAuto ), nhLevsAuto )
            call color ( 'BLACK' )
            call gridMp ( 1, 1)
            call endGrf ()
    
            call axsPos ( 2550, 1600 )
            call axsLen ( 1500, 800 )
            call projct ( 'LAMBERT' )
            call mapPol ( 0.0, 90.0 )
            call grafMp ( -180.0, 180.0, -180.0, 30.0, -90.0, 90.0, -90.0, 30.0 )
            call conShd ( real ( phi(1,:,1) * radToDeg ), nU2, &
                real ( 90.0 - theta(:,1,1) * radToDeg ), nU1, &
                !real ( transpose ( hThFit_h3 * u0 * 1d9 ) ), &
                real ( transpose ( hThFit_h3 ) ), &
                real ( hLevsAuto ), nhLevsAuto )
            call color ( 'BLACK' )
            call gridMp ( 1, 1)
            call endGrf ()

            call axsPos ( 2550, 2399 )
            call axsLen ( 1500, 800 )
            call projct ( 'LAMBERT' )
            call mapPol ( 0.0, 90.0 )
            call grafMp ( -180.0, 180.0, -180.0, 30.0, -90.0, 90.0, -90.0, 30.0 )
            call conShd ( real ( phi(1,:,1) * radToDeg ), nU2, &
                real ( 90.0 - fitTheta0_h3(:,1) * radToDeg ), nU1, &
                !real ( transpose ( hThFit_h3_manual * u0 * 1d9 ) ), &
                real ( transpose ( 8.0d22 * sin ( fitTheta0_h3 ) / ( rE ** 3 ) * 1d9 ) ), &
                real ( hLevsAuto ), nhLevsAuto )
            call color ( 'BLACK' )
            call gridMp ( 1, 1)
            call endGrf ()
        
            call axsPos ( 1850, 1600 )
            call axsLen ( 1500, 800 )
            call projct ( 'LAMBERT' )
            call mapPol ( 1.0, 90.0 )
            call grafMp ( -180.0, 180.0, -180.0, 30.0, -90.0, 90.0, -90.0, 30.0 )
            call conShd ( real ( phi(1,:,1) * radToDeg ), nU2, &
                real ( 90.0 - fitTheta0_h3(:,1) * radToDeg ), nU1, &
                !real ( transpose ( hPhFit_h3 * u0 * 1d9 ) ), &
                real ( transpose ( hPhFit_h3 ) ), &
                real ( hLevsAuto ), nhLevsAuto )
            call color ( 'BLACK' )
            call gridMp ( 1, 1)
            call endGrf ()

            call axsPos ( 1850, 2399 )
            call axsLen ( 1500, 800 )
            call projct ( 'LAMBERT' )
            call mapPol ( 0.0, 90.0 )
            call grafMp ( -180.0, 180.0, -180.0, 30.0, -90.0, 90.0, -90.0, 30.0 )
            call conShd ( real ( phi(1,:,1) * radToDeg ), nU2, &
                real ( 90.0 - fitTheta0_h3(:,1) * radToDeg ), nU1, &
                real ( transpose ( hPhFit_h3_manual * u0 * 1d9 ) ), &
                real ( hLevsAuto ), nhLevsAuto )
            call color ( 'BLACK' )
            call gridMp ( 1, 1)
            call endGrf ()

            !write (*,*)
            !write (*,*) hPhFit_h3(20,:) * u0 *1d9
            !read (*,*)

            do i = 1, nU1
                do j = 1, nU2
                    
                    E1covIono_dU2(i,j)  = 1d0 / ( rI * sin ( fitTheta0_e1(i,j) ) ) * &
                        qdDer ( 1, u2 (j), u2, E1covIono (i,:) )

                    E2covIono_dU1(i,j)  = 1d0 / rI * &
                        qdDer ( 1, fitTheta0_e2 (i,j), fitTheta0_e2(:,j), E2covIono(:,j) )

                end do
            end do
    
            call axsPos ( 1050, 2399 )
            call axsLen ( 1500, 800 )
            call projct ( 'LAMBERT' )
            call mapPol ( 0.0, 90.0 )
            call grafMp ( -180.0, 180.0, -180.0, 30.0, -90.0, 90.0, -90.0, 30.0 )
            call conShd ( real ( u2 * radToDeg ), nU2, &
                real ( 90.0 - fitTheta0_e2(:,1) * radToDeg ), nU1, &
                real ( transpose ( E2covIono_dU1 * rE * 1e3 ) ), &
                real ( eLevsAuto ), neLevsAuto )
            call color ( 'BLACK' )
            call gridMp ( 1, 1)
            call endGrf ()

        !write (*,*)
        !   write (*,*) e1covIono_du2(:,plotSlice)
        !   write (*,*)
        !   write (*,*) e1covIono(:,plotSlice)
        !   write (*,*)
            
            call axsPos ( 1050, 600 )
            call axsLen ( 1500, 800 )
            call projct ( 'LAMBERT' )
            call mapPol ( 0.0, 90.0 )
            call grafMp ( -180.0, 180.0, -180.0, 30.0, -90.0, 90.0, -90.0, 30.0 )
            call conShd ( real ( u2_e1 * radToDeg ), nU2, &
                real ( 90.0 - fitTheta0_e1(:,1) * radToDeg ), nU1, &
                real ( transpose ( E1covIono_dU2 * rE * 1e3 ) ), &
                real ( eLevsAuto ), neLevsAuto )
            call color ( 'BLACK' )
            call gridMp ( 1, 1)
            call endGrf ()

            call frame ( 1 ) ! Set plot frame thickness to 0
            call setGrf ('NAME', 'NAME', 'TICKS', 'TICKS')
            call color ( 'BLACK' )
            call axsPos ( 200, 1500 )
            call axsLen ( 4000, 700 )
            call noClip ()
            call noChek ()
            call incMrk (1)
            call labels ( 'FEXP', 'XYZ' )
            call graf ( 10.0, 170.0, 0.0, 20.0, &
                real ( -maxVal ( abs ( hrFit_h3(:,:) ) * u0 * 1e9 ) ), &
                real ( maxVal ( abs ( hrFit_h3(:,:) ) * u0 * 1e9 ) ), &
                real ( -maxVal ( abs ( hrFit_h3(:,:) * u0 * 1e9 ) ) ), &
                real ( maxVal ( abs ( hrFit_h3(:,:) * u0 * 1e9 ) ) ) /  4.0 )
            !call graf ( 10.0, 170.0, 0.0, 20.0, -0.5, 0.5, -0.5, 0.1 ) 

!           call graf ( 10.0, 170.0, 0.0, 20.0, &
!               real ( -maxVal ( abs ( hrBFnArr_h3(16,nU1+1:2*nU1) ) ) ), &
!               real ( maxVal ( abs ( hrBFnArr_h3(16,nU1+1:2*nU1) ) ) ), &
!               real ( -maxVal ( abs ( hrBFnArr_h3(16,nU1+1:2*nU1) ) ) ), &
!               real ( maxVal ( abs ( hrBFnArr_h3(16,nU1+1:2*nU1) ) ) ) /   4.0 )
!           write (*,*) hrBFnArr_h3(16,nU1+1:2*nU1)
            !write (*,*)  theta (:,2,1) * radToDeg

            !   h3con
            !call curve ( real ( theta (:,2,1) * radToDeg ), &
            !   real ( hrBFnArr_h3(16,nU1+1:2*nU1) ), nU1 )
        !   call color ( 'GREEN' )
        !   call curve ( real ( theta (:,plotSlice,1) * radToDeg ), &
        !       real ( h3con(:,plotSlice,1) * u0 * 1e9 ), nU1 )
            call color ( 'BLACK' )  
            call curve ( real ( theta (:,plotSlice,1) * radToDeg ), &
                real ( hrRealVector ( 1, ( plotSlice - 1 ) * nU1 + 1 : plotSlice * nU1 ) * u0 * 1e9 ), nU1 )

        !   call curve ( real ( fitTheta0_e1(:,plotSlice) * radToDeg ), real ( h3con(:,plotSlice,2) * sqrt ( g33covIono(:,plotSlice) ) * u0 * 1e9 ), nU1 )
        !   call curve ( real ( theta (:,plotSlice,1) * radToDeg ), &
        !       real ( h3con(:,plotSlice,2) * sqrt ( g33covIono(:,plotSlice) ) * u0 * 1e9 ), nU1 )
    
        !   !call curve ( real ( theta (:,4,1) * radToDeg ), &
        !   !   real ( h3con(:,4,1) * sqrt ( g33covIono(:,4) ) * u0 * 1e9 ), nU1 )


        !   !call curve ( theta (:,2,nU3) * radToDeg, h3con(:,2,nU3) * sqrt ( g33covIonoS(:,2) ) * u0 * 1e9, nU1 )

        !   !   h3con fit (hrFit_h3)
        !   call color ( 'BLACK' )
            call curve ( real ( fitTheta0_h3(:,plotSlice) * radToDeg ), real ( hrFit_h3(:,plotSlice) * u0 * 1e9 ), nU1 )
            call curve ( real ( fitTheta0_h3(:,plotSlice+1) * radToDeg ), real ( hrFit_h3(:,plotSlice+1) * u0 * 1e9 ), nU1 )
        !   call curve ( real ( fitTheta0_e2(:,plotSlice) * radToDeg ), real ( hrFit_e2(:,plotSlice) * u0 * 1e9 ), nU1 )
            
            call color ( 'MAGENTA' )    
            call curve ( real ( fitTheta0_h3(:,plotSlice) * radToDeg ), real ( hrFit_h3(:,plotSlice) * u0 * 1e9 ), nU1 )
            call curve ( real ( fitTheta0_e2(:,plotSlice+1) * radToDeg ), real ( hrFit_e2(:,plotSlice+1) * u0 * 1e9 ), nU1 )

            !call curve ( real ( theta (:,plotSlice,1) * radToDeg ), real ( hrFit_h3(:,plotSlice) ), nU1 )
            !call curve ( real ( fitTheta0_h3 * radToDeg ), real ( hrFit_h3(:,4) * u0 * 1e9 ), nU1 )
            !!call curve ( theta (:,plotSlice,nU3) * radToDeg, hrFit_h3S(:,plotSlice) * u0 * 1e9, nU1 )

            !   hThFit
            call color ( 'BLUE' )
        !   call curve ( real ( fitTheta0_e2(:,plotSlice) * radToDeg ), real ( -dhPh_e2(:,plotSlice) * u0 * 1e9 ), nU1 )
            call curve ( real ( fitTheta0_e2(:,plotSlice+1) * radToDeg ), real ( hThFit_e2(:,plotSlice+1) * u0 * 1e9 ), nU1 )
            call curve ( real ( fitTheta0_e2(:,plotSlice) * radToDeg ), real ( hThFit_e2(:,plotSlice) * u0 * 1e9 ), nU1 )
        !   call curve ( real ( fitTheta0_e2(:,plotSlice) * radToDeg ), real ( h1cov(:,plotSlice,1) * sqrt ( g11conIono(:,plotSlice) ) * u0 * 1e9  ), nU1 )
        !   call curve ( real ( fitTheta0_e1(:,plotSlice) * radToDeg ), &
        !       real ( h1cov(:,plotSlice,1) * sqrt ( g11conIono(:,plotSlice) ) * u0 * 1e9 ), nU1 )

    !       call curve ( real ( fitTheta0_e2(:,plotSlice) * radToDeg ), real ( hThFit_e2_manual(:,plotSlice) * u0 * 1e9  ), nU1 )
    !       !call color ( 'BLUE' )
    !       call curve ( real ( fitTheta0_h3(:,plotSlice) * radToDeg ), real ( hThFit_h3(:,plotSlice) * u0 * 1e9 ), nU1 )
    !       call curve ( real ( fitTheta0_h3(:,plotSlice) * radToDeg ), real ( hThFit_h3_manual(:,plotSlice) * u0 * 1e9  ), nU1 )
    !       !call color ( 'ORANGE' )
    !       call curve ( real ( fitTheta0_e1(:,plotSlice) * radToDeg ), real ( hThFit_e1(:,plotSlice) * u0 * 1e9 ), nU1 )
    !       call curve ( real ( fitTheta0_e1(:,plotSlice) * radToDeg ), real ( hThFit_e1_manual(:,plotSlice) * u0 * 1e9  ), nU1 )



    !       !call curve ( ( pi - fitTheta0_e1 ) * radToDeg, hThFit_e1S(:,plotSlice) * u0 * 1e9, nU1 )
    !       !
            !   hPhFit
            call color ( 'RED' )
            !call curve ( real ( fitTheta0_e2(:,plotSlice) * radToDeg ), real ( dhTh_e2(:,plotSlice) * u0 * 1e9 ), nU1 )
            call curve ( real ( fitTheta0_e2(:,plotSlice) * radToDeg ), real ( hPhFit_e2(:,plotSlice) * u0 * 1e9 ), nU1 )
            call curve ( real ( fitTheta0_e2(:,plotSlice+1) * radToDeg ), real ( hPhFit_e2(:,plotSlice+1) * u0 * 1e9 ), nU1 )
        !   call curve ( real ( fitTheta0_e1(:,plotSlice) * radToDeg ), real ( h2cov(:,plotSlice,1) * sqrt ( g22conIono(:,plotSlice) ) * u0 * 1e9  ), nU1 )

            !call curve ( real ( fitTheta0_e1(:,plotSlice) * radToDeg ), &
            !   real ( h2cov(:,plotSlice,1) * sqrt ( g22conIono(:,plotSlice) ) * u0 * 1e9 ), nU1 )

        !   call curve ( real ( fitTheta0_e2(:,plotSlice) * radToDeg ), real ( h2cov_h1(:,plotSlice) * u0 * 1e9 ), nU1 )

    !       call curve ( real ( fitTheta0_e2(:,plotSlice) * radToDeg ), real ( hPhFit_e2_manual(:,plotSlice) * u0 * 1e9  ), nU1 )
    !       
    !       call curve ( real ( fitTheta0_e1(:,plotSlice) * radToDeg ), real ( hPhFit_e1(:,plotSlice) * u0 * 1e9 ), nU1 )
    !       call curve ( real ( fitTheta0_e1(:,plotSlice) * radToDeg ), real ( hPhFit_e1_manual(:,plotSlice) * u0 * 1e9  ), nU1 )
    !       
    !       call curve ( real ( fitTheta0_h3(:,plotSlice) * radToDeg ), real ( hPhFit_h3(:,plotSlice) * u0 * 1e9 ), nU1 )
    !       call curve ( real ( fitTheta0_h3(:,plotSlice) * radToDeg ), real ( hPhFit_h3_manual(:,plotSlice) * u0 * 1e9  ), nU1 )
            
            !!call curve ( real ( fitTheta0_e2 * radToDeg ), real ( h2cov_h1(:,plotSlice) * sqrt ( g22conIono(:,plotSlice) ) * u0 * 1e9  ), nU1 )

            !call curve ( ( pi - fitTheta0_e1 ) * radToDeg, hPhFit_e1S(:,plotSlice) * u0 * 1e9, nU1 )

            call endGrf ()
        
            !   IONOSPHERE EFIELD PLOTS
            call color ( 'BLACK')
            call axsPos ( 200, 2300 )
            call axsLen ( 4000, 700 )
        !   call graf ( 10.0, 170.0, 0.0, 20.0, &
        !       real ( -maxVal ( abs ( sqrt ( (E1covIono(:,plotSlice) * rE * 1e3)**2 + (E2covIono(:,plotSlice) * rE * 1e3)**2 ) ) ) )/1.0, &
        !       real ( maxVal ( abs ( sqrt ( (E1covIono(:,plotSlice) * rE * 1e3)**2 + (E2covIono(:,plotSlice) * rE * 1e3)**2 ) ) ) )/1.0, &
        !       real ( -maxVal ( abs ( sqrt ( (E1covIono(:,plotSlice) * rE * 1e3)**2 + (E2covIono(:,plotSlice) * rE * 1e3)**2 ) ) ) )/1.0, &
        !       real ( maxVal ( abs ( sqrt ( (E1covIono(:,plotSlice) * rE * 1e3)**2 + (E2covIono(:,plotSlice) * rE * 1e3)**2 ) ) ) )/ 1.0 )
            call graf ( 10.0, 170.0, 0.0, 20.0, -0.1, 0.1, -0.1, 0.05 )


        !   call graf ( 10.0, 170.0, 0.0, 20.0, &
        !       real ( -maxVal ( abs ( sqrt ( (sigma11_N(:,plotSlice) )**2 + (E2covIono(:,plotSlice) )**2 ) ) ) ), &
        !       real ( maxVal ( abs ( sqrt ( (sigma11_N(:,plotSlice) )**2 + (E2covIono(:,plotSlice) )**2 ) ) ) ), &
        !       real ( -maxVal ( abs ( sqrt ( (sigma11_N(:,plotSlice) )**2 + (E2covIono(:,plotSlice) )**2 ) ) ) ), &
        !       real ( maxVal ( abs ( sqrt ( (sigma11_N(:,plotSlice) )**2 + (E2covIono(:,plotSlice) )**2 ) ) ) )/ 4.0 )
        !   call curve ( real ( fitTheta0_e1(:,plotSlice) * radToDeg ), real (sigma11_N(:,plotSlice) ), nU1 )
        !   call curve ( real ( fitTheta0_e1(:,plotSlice) * radToDeg ), real (sigma12_N(:,plotSlice) ), nU1 )
        !   call curve ( real ( fitTheta0_e1(:,plotSlice) * radToDeg ), real (sigma21_N(:,plotSlice) ), nU1 )
        !   call curve ( real ( fitTheta0_e1(:,plotSlice) * radToDeg ), real (sigma22_N(:,plotSlice) ), nU1 )

        !   call graf ( 10.0, 170.0, 0.0, 20.0, &
        !       real ( -maxVal ( abs ( sqrt ( (dhTh_e1(:,plotSlice) )**2 + (dhPh_e1(:,plotSlice) )**2 ) ) ) ), &
        !       real ( maxVal ( abs ( sqrt ( (dhTh_e1(:,plotSlice) )**2 + (dhPh_e1(:,plotSlice) )**2 ) ) ) ), &
        !       real ( -maxVal ( abs ( sqrt ( (dhTh_e1(:,plotSlice) )**2 + (dhPh_e1(:,plotSlice) )**2 ) ) ) ), &
        !       real ( maxVal ( abs ( sqrt ( (dhTh_e1(:,plotSlice) )**2 + (dhPh_e1(:,plotSlice) )**2 ) ) ) )/ 4.0 )
        !   call curve ( real ( fitTheta0_e1(:,plotSlice) * radToDeg ), real (dhTh_e1(:,plotSlice) ), nU1 )
        !   call curve ( real ( fitTheta0_e1(:,plotSlice) * radToDeg ), real (dhPh_e1(:,plotSlice) ), nU1 )
        !   call curve ( real ( fitTheta0_e2(:,plotSlice) * radToDeg ), real (dhTh_e2(:,plotSlice) ), nU1 )
        !   call curve ( real ( fitTheta0_e2(:,plotSlice) * radToDeg ), real (dhPh_e2(:,plotSlice) ), nU1 )

            !   eTh
            call color ( 'BLUE' )
        !   call curve ( real ( fitTheta0_e1(:,plotSlice) * radToDeg ), real ( E1covIono(:,plotSlice) * rE * 1e3 ), nU1 )
            call color ( 'GREEN' )
        !   call curve ( real ( fitTheta0_e1(:,plotSlice) * radToDeg ), real ( E1covIono_dU2(:,plotSlice) * rE * 1e3 / 30.0), nU1 )
    !       call curve ( real ( fitTheta0_e1(:,plotSlice) * radToDeg ), real ( e1cov(:,plotSlice,1) * sqrt ( g11conIono(:,plotSlice ) ) * rE * 1e3 ), nU1 )
    !       call curve ( real ( fitTheta0_e1(:,plotSlice) * radToDeg ), real ( e1cov(:,plotSlice,2) * sqrt ( g11conIono(:,plotSlice ) ) * rE * 1e3 ), nU1 )
            !call curve ( ( pi - fitTheta0_e1 ) * radToDeg, E1covIonoS(:,plotSlice) * rE * 1e3, nU1 )

            !   ePh 
            call color ( 'RED' )
        !   call curve ( real ( fitTheta0_e2(:,plotSlice) * radToDeg ), real ( E2covIono(:,plotSlice) * rE * 1e3 ), nU1 )
        !   call curve ( real ( fitTheta0_e2(:,plotSlice) * radToDeg ), real ( e2cov(:,plotSlice,1) * sqrt ( g22conIono(:,plotSlice ) ) * rE * 1e3 ), nU1 )
            call color ( 'MAGENTA' )
        !   call curve ( real ( fitTheta0_e2(:,plotSlice) * radToDeg ), real ( e2cov(:,plotSlice,2) * sqrt ( g22conIono(:,plotSlice ) ) * rE * 1e3 ), nU1 )
        !   call curve ( real ( fitTheta0_e2(:,plotSlice) * radToDeg ), real ( E2covIono_dU1(:,plotSlice) * rE * 1e3 / 30.0), nU1 )

            !call curve ( ( pi - fitTheta0_e2 ) * radToDeg, E2covIonoS(:,plotSlice) * rE * 1e3, nU1 )

            call endGrf ()
            !read (*,*)
        !   write (*,*) theta (:,2,1) * radToDeg
        !   write (*,*) h3con(:,2,1) * sqrt ( g33covIono(:,2) ) * u0 * 1e9
        !   write (*,*) hrFit_h3(:,2) * u0 * 1e9
        !   write (*,*) coeffsOut


        !if ( t > 200 ) read (*,*)

                !write (*,*) E2covIono(:,2) * rE * 1e3

        end if
    
        if ( mod ( t, progFreq ) < dT ) then  

            timeTaken   = end_timer ( startTimeOLD )
            eta = timeTaken / progFreq * ( 1000.0 - t ) / 3600.0
            
            if ( eta < 1 ) then 
                write (*,'(a12,f5.1,a6,f4.1,a8)') 'Model time: ', t, &
                    ' ETA: ', eta * 60.0, &
                    ' [MINS]'
            else
                write (*,'(a12,f5.1,a6,f4.1,a8)') 'Model time: ', t, &
                    ' ETA: ', eta, &
                    ' [HOURS]'
            end if

        end if

        !read (*,*)

    end do timeLoop

    !call disFin () 

  !read (*,*)

end program mhd3d
