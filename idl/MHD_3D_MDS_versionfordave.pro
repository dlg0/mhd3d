
;	Model 3D Magnetosphere with Realistic Ionospheric Boundaries
;
;	Version 0.2		19 Nov 2008
;
;   Murray D Sciffer
; 	Centre for Space Physics
; 	University of Newcastle
;
; 	Mods:   search for "<Dave>" for comments
;
; ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** *****

pro MHD_3D_MDS_versionfordave

; Directories

;	<DAVE>	To run create a directory called Blah with a subdirectory called "PNG" (you get the idea!)

 Plot_Dir	='D:\Data\MHD_3D\Blah\'						;	Directory for data files

; 	Parameters for the Model
 N_sec		= 1;00;0001									;	Number of Seconds of model run

;	Number of Grid Points
 Num_u1		= 50                                       	; 	Number of field lines, u1 coords.
 Num_u2		= 36                                       	; 	Number of Aziamuthal sliaces, u2 coords .
 Num_u3		= 51                                       	; 	Number of points on field lines, u3 coords.

; 	<DAVE> not required for this version to run
;	Ionospheric Data directory + filename
 Fn_Iono_N 	= 'D:\MurraysDocuments\Pen Drive 8G\Data\2D_iono\COND_Yuki_16UT_N.sig'		; 	File containing Northern Ionospheres Sigma
 Fn_Iono_S 	= 'D:\MurraysDocuments\Pen Drive 8G\Data\2D_iono\COND_Yuki_16UT_S.sig'		; 	File containing Southern Ionospheres Sigma

; 	<DAVE> not sure yet as to how many basis function to use ona particular grid size. There must be a number which optimises comp time and
;	coverage but haven't come to a final conclusion... any ideas?

;	Number of basis function
 Num_L		= 15										; 	Number of Legendre's used in Spherical harmonic expansion
 Num_M		= 10										;	Number of Azimuthal basis sets (including Zero, actually use -(m-1),...,0,1,..,m-1

;	Time Driver Parameters
 Amp		= 1.0e-8									;	Max amplitude of Driver (in T)

;	Plot Settings
 Plot_Switch= 0											;	Plot to screen 1 = Plots to window
 Plot_PNG 	= 1											;	Switch to create PNG files of each plot,  1 = Create Files File
 XSz		= 750										;	X size of Plot Window Size
 YSz		= 750										;	Y size of Plot Window Size
 win_num	= 1											;	Plot window number

;	Distances (for scaling)
 Re			= double(6370.e3)							 ;	One Earth Radius (in m)
 Ri			= 1.02                                       ; 	Ionos. radius in Re
 Ree		= 1.00                                       ; 	Earth. radius in Re

;	Grid space paramters
 LMin_u1	= 1.2                                        ; 	Min L Value for u1 Lines
 LMax_u1	= 10.0                                       ; 	Max L Value
 LMin_u3	= Ri
 LMax_u3	= 55.0

;	Save file names for Data
 Iono_file  = 'Mhd3d_Iono'
 Grnd_file  = 'Mhd3d_Grnd'
 Parm_file 	= 'Mhd3d_Parm'
 Plot_file 	= 'Mhd3d_Tme_'

;	Physical constants (unit for length in model is Re)

 c			= 2.997e8									;	Speed of light in a Vacuum (in SI units)
 u0_u		= 4.0*!Dpi*1.0e-7							;	Magnetic Permeability (in SI units)
 e0 		= 1.0/(u0_u*c^2)							;	Dielectric Constant for a vacuum (in SI units)
 im			= dcomplex(0.0,1.0)							;	Imaginary Number

 c			= c/Re										;	Scaled speed of light (in Re per sec)
 u0			= u0_u/Re									;	Scaled Permablity
 e0			= e0*Re^3									;	Scaled Permativity

 Num_u1		= 2*Fix(Num_u1/2) 							;	Must be even
 Num_u2		= 2*Fix(Num_u2/2) 							;	Must be even
 Num_u3		= 2*Fix(Num_u3/2)+1 						;	Must be odd

 Starttime 	= SYSTIME(/SECONDS)							;	Mark this computations start time

;	Call to generate Grid for 3D Model (using Lysak 2004 - Tilt Dipole)

gen_grid,	Num_u1,Num_u1_2,Num_u3,Num_u3_2,Num_u2,LMin_u1,LMax_u1, LMin_u3, LMax_u3,$
            Ri,del_Th_0_u1,del_Th,$
            u1Arr,u2Arr,u3Arr,RArr,ThArr,PhiArr,LVArr,XArr,YArr,ZArr,$
            g11_con,g13_con,g22_con,g33_con,$
            g11_cov,g13_cov,g22_cov,g33_cov,$
            Jacb,Colat_N,Colat_S,h_nu,h_ph,h_mu,$
           	h_ra_n,h_th_n,h_ph_n,h_ra_s,h_th_s,h_ph_s

;	Call to generate Va profile (include Oxygen and plume is needed)

gen_Va,		Lvarr,Num_u1,Num_u2,Num_u3,Rarr,ThArr,PhiArr,Xarr,rho_a,Va

Courant,	Xarr,Zarr,Va,Dt,Nt,Num_u1,Num_u2,Num_u3,Delta_L,Trans_T,N_sec,Eperp,V2,Sig_a

Iono_Cond,	Num_u1,Num_u2,Num_u3,Colat_N,Colat_S,Sd_N,Sp_N,Sh_N,INVSS_N,Sd_S,Sp_S,Sh_S,$
			INVSS_S,Fn_Iono_N,Fn_Iono_S,ThArr,PhiArr,Iono_file,Plot_Dir

Sph_Harm,	Num_u1,Num_u2,Num_L,Num_M,LMin_u1,LMax_u1,Ri,Ree,SHF,SHF_dr,Alpha,Alpha_dr,Zeros,PhiArr,del_Th


IO_Parm,	Xarr,Yarr,Zarr,Num_u1,Num_u2,Num_u3,Delta_L,Trans_T,rho_a,Va,Plot_Dir,Parm_file,Iono_file,$
 			g11_con,g13_con,g22_con,g33_con, g11_cov,g13_cov,g22_cov,g33_cov,N_sec,Eperp,V2,Sig_a,$
            Jacb,Colat_N,Colat_S,h_nu,h_ph,h_mu,LMin_u1,LMax_u1, LMin_u3, LMax_u3,RArr,ThArr,PhiArr,$
            LVArr,Sd_S,Sp_S,Sh_S,Sd_N,Sp_N,Sh_N,SHF,SHF_dr,Alpha,Alpha_dr,Zeros,$
           	h_ra_n,h_th_n,h_ph_n,h_ra_s,h_th_s,h_ph_s

;	Call to time iterate ULF wave

Iterate,	Num_u1,Num_u2,Num_u3,g11_cov,g13_cov,g22_cov,g33_cov,g11_con,jacb,u1Arr,u2Arr,u3Arr,V2,Dt,Nt,$
			Xarr,Yarr,Zarr,Rarr,ThArr,PhiArr,Plot_Dir,Plot_File,h_nu,h_ph,h_mu,SHF,SHF_dr,Alpha,Alpha_dr,$
			del_Th,u0,Colat_N,Colat_S,INVSS_S,INVSS_N,Ri,h_ra_n,h_th_n,h_ph_n,h_ra_s,h_th_s,h_ph_s


PRINT, 'Memory required (Kb): ', MEMORY(/HIGHWATER)/1024; - start_mem


stop
END			; END PROGRAM (MHD_3D Pro)


; ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** *****

 Pro gen_grid,	Num_u1,Num_u1_2,Num_u3,Num_u3_2,Num_u2,LMin_u1,LMax_u1, LMin_u3, LMax_u3,$
            	Ri,del_Th_0_u1,del_Th,$
           	 	u1Arr,u2Arr,u3Arr,RArr,ThArr,PhiArr,LVArr,XArr,YArr,ZArr,$
            	g11_con,g13_con,g22_con,g33_con,$
            	g11_cov,g13_cov,g22_cov,g33_cov,$
            	Jacb,Colat_N,Colat_S,h_nu,h_ph,h_mu,$
            	h_ra_n,h_th_n,h_ph_n,h_ra_s,h_th_s,h_ph_s

; Set up spatial grid and scale factors for the Non-Orthogonal Coords (using Lysak 2004)

; Set_up Parameters for u1 lines
 Num_u1_2	=fix(Num_u1/2.0)  									; Number of u1 points per cell (see discretisation scheme)
 Num_u1		=fix(2*Num_u1_2)
; Set_up Parameters for u2 lines
 Num_u2_2	=fix(Num_u2/2.0)  									; Number of u1 points per cell (see discretisation scheme)
 Num_u2		=fix(2*Num_u2_2)
; Set_up Parameters for u3 lines
 Num_u3_2	=fix(Num_u3/2.0)                                  	; Number of u3 points per hemisphere
 Num_u3		=fix(2*Num_u3_2)+1									; Plus 1/2 cell for boundaries (in schematic)
; Determine del_Theta for u1 calcs...
 CLat_Min_u1	= asin(sqrt(RI/LMax_u1))*180.0/!pi                       ; Min Co_Lat in degrees for Nth Hemisphere (Max Latitude), specifies outer boundary
 CLat_Max_u1	= asin(sqrt(RI/LMin_u1))*180.0/!pi                       ; Max Co_Lat in degrees for Nth Hemisphere (Min Latitude), specifies inner boundary
 del_Th_0_u1	= double(CLat_Max_u1-CLat_Min_u1)/double(Num_u1-1)       ; del_Theta for u1 points
 del_Th			= del_Th_0_u1*!dpi/180.0								;  del_Thet in radians

 Print,'Non-Orthogonal Grid Generation:'
 Print,'L_Value_Min, L_Value_Max, CLat_Min_u1, CLat_Max_u1 [Deg] : ',LMin_u1,LMax_u1,CLat_Min_u1,CLat_Max_u1
 Print,'Number of u1 Lines, Del_Theta_u1 : ',Num_u1,del_Th_0_u1

 CLat_Min_u3	= asin(sqrt(RI/LMax_u3))*180.0/!pi
 CLat_Max_u3	= asin(sqrt(RI/LMin_u3))*180.0/!pi
 del_Th_0_u3	= double(CLat_Max_u3-CLat_Min_u3)/double(Num_u3_2-1)

 Print,'Number of u3 Lines, Del_Theta_u3 : ',Num_u3,del_Th_0_u3

 u1Arr=dblarr(Num_u1,Num_u2,Num_u3) & u3Arr=dblarr(Num_u1,Num_u2,Num_u3)

 RArr =dblarr(Num_u1,Num_u2,Num_u3) & ThArr=dblarr(Num_u1,Num_u2,Num_u3) & LVArr=dblarr(Num_u1)

 xArr =dblarr(Num_u1,Num_u2,Num_u3) & yArr =dblarr(Num_u1,Num_u2,Num_u3) & zArr = dblarr(Num_u1,Num_u2,Num_u3)

 g11_con=dblarr(Num_u1,Num_u2,Num_u3) & g13_con=dblarr(Num_u1,Num_u2,Num_u3)
 g22_con=dblarr(Num_u1,Num_u2,Num_u3) & g33_con=dblarr(Num_u1,Num_u2,Num_u3)

 g11_cov=dblarr(Num_u1,Num_u2,Num_u3) & g13_cov=dblarr(Num_u1,Num_u2,Num_u3)
 g22_cov=dblarr(Num_u1,Num_u2,Num_u3) & g33_cov=dblarr(Num_u1,Num_u2,Num_u3)

 Jacb=dblarr(Num_u1,Num_u2,Num_u3)

 Colat_N=dblarr(Num_u1)									; Dip Angle at Along Northern Ionosphere
 Colat_S=dblarr(Num_u1)									; Dip Angle at Along Southern Ionosphere

 h_nu=DblArr(Num_u1,Num_u2,Num_u3) & h_ph=DblArr(Num_u1,Num_u2,Num_u3) & h_mu=DblArr(Num_u1,Num_u2,Num_u3)

;  u2 (L shell) Loop
 For jj=0,Num_u2-1 do $									; Starts at inner-most field line (i.e. smallest L)
 Begin

;  u1 (L shell) Loop
 For ii=0,Num_u1-1 do $									; Starts at inner-most field line (i.e. smallest L)
 Begin

; Northern Hemisphere
  CLat_0_u1				= CLat_Max_u1-float(ii)*del_Th_0_u1           ; Co_Lat at Nth hemisphere Ionosphere (in degrees)
  CLat_0r_u1			= CLat_0_u1*!pi/180.0                        ; Convert CLat_0 to radians
  u1					=-sin(CLat_0r_u1)*sin(CLat_0r_u1)                   ; Determine u1 (which is constant along an L shell) for this L value
  Colat_N(ii) 			= CLat_0r_u1
  Colat_S(ii) 			= !Dpi-CLat_0r_u1
  LVArr(ii)				= 1.0/(cos(!pi/2.0-CLat_0r_u1)^2)

  For kk=0,Num_u3_2-1 do $
  Begin

   ; CLat_0_u3=CLat_Max_u3-float(ii)*del_Th_0_u3          ; Co_Lat at Nth hemisphere Ionosphere (in degrees)
   ;	@@@@ Changed here for better grid in u3
   CLat_0_u3		= CLat_Max_u3-(CLat_Max_u3- CLat_Min_u3)* (float(kk)/(Num_u3_2-1))^0.5      ; Co_Lat at Nth hemisphere Ionosphere (in degrees)
   CLat_0r_u3		= CLat_0_u3*!pi/180.0                       ; Convert CLat_0 to radians
   u3				= sin(CLat_0r_u3)^4                                 ; Determine u3 for Nth (determined at CLat=0)
   u1Arr(ii,jj,kk)	= u1
   u3Arr(ii,jj,kk)	= u3
   rg				= RI
   New_r,rg,u1,u3,CLat_0r_u1,RI,ans
   RArr(ii,jj,kk)		= ans                                      ; r value for the intersection
   cos_Th				= u3*cos(CLat_0r_u1)*ans*ans/(RI*RI)            ; see eqn (11) Lysak, 2004
   one_p_3cos_sq		= 1.0+3.0*cos_Th*cos_Th
   Theta				= acos(cos_Th)
   sin_Th				= sin(Theta)
   ThArr(ii,jj,kk)		= Theta
   Lat					= !pi/2.0-Theta
   g11_con(ii,jj,kk)	= (RI^2/ans^4)*sin_Th^2*one_p_3cos_sq
   g13_con(ii,jj,kk)	=-(RI^3/ans^5)*sin(CLat_0r_u1)^2*cos_Th*one_p_3cos_sq/(2.0*cos(CLat_0r_u1)^3)
   g22_con(ii,jj,kk)	= 1.0/(ans^2*sin_Th^2)
   g33_con(ii,jj,kk)	= RI^4/(ans^6*cos(CLat_0r_u1)^6)*(0.25*cos_Th^2*(1.0+3.0*cos(CLat_0r_u1)^2)^2 + sin_Th^2*(1.0-RI/ans)^2)
   g11_cov(ii,jj,kk)	= ans^4/(RI^2*cos(CLat_0r_u1)^4*one_p_3cos_sq^2)*((1.0-RI/ans)^2 + 0.25*(cos_Th/sin_Th)^2*(1.0+3.0*cos(CLat_0r_u1)^2)^2)
   g13_cov(ii,jj,kk)	= ans^4*cos_Th/(2.0*RI^2*cos(CLat_0r_u1)*one_p_3cos_sq)
   g22_cov(ii,jj,kk)	= ans*ans*sin_Th*sin_Th
   g33_cov(ii,jj,kk)	= ans^6*cos(CLat_0r_u1)^2/(RI^4*one_p_3cos_sq)
   Jacb(ii,jj,kk)		= ans^6*cos(CLat_0r_u1)/(RI^3*one_p_3cos_sq)
; assumes Re=1
   h_nu(ii,jj,kk)	= ans^2/(sin_Th*sqrt(one_p_3cos_sq))   	; At the Equator this is in the R direction
   h_ph(ii,jj,kk)	= ans*sin_Th                          	; E_W, Azimuthal
   h_mu(ii,jj,kk)	= ans^3/(sqrt(one_p_3cos_sq))         	; At the Equator this is in the Theta direction (FA)

; Sothern Hemisphere
; Mirror hemisphere by changing Theta
   u1Arr(ii,jj,Num_u3-1-kk)		= u1                               ; Keep same u3 value
   u3Arr(ii,jj,Num_u3-1-kk)		=-u3                               ; Keep same u3 value
   RArr(ii,jj,Num_u3-1-kk)		= ans                             ; Keep same r value
   Theta						= !pi-ThArr(ii,jj,kk)                               ; Change Theta here for Sth hemisphere
   ThArr(ii,jj,Num_u3-1-kk)		= Theta
   Lat							= (!pi/2.0-Theta)
   cos_Th						= cos(Theta)
   one_p_3cos_sq				= 1.0+3.0*cos_Th*cos_Th
   sin_Th						= sin(Theta)
   g11_con(ii,jj,Num_u3-1-kk)	= (RI^2/ans^4)*sin_Th^2*one_p_3cos_sq
   g13_con(ii,jj,Num_u3-1-kk)	= -(RI^3/ans^5)*sin(CLat_0r_u1)^2*cos_Th*one_p_3cos_sq/(2.0*cos(CLat_0r_u1)^3)
   g22_con(ii,jj,Num_u3-1-kk)	= 1.0/(ans^2*sin_Th^2)
   g33_con(ii,jj,Num_u3-1-kk)	= RI^4/(ans^6*cos(CLat_0r_u1)^6)*(0.25*cos_Th^2*(1.0+3.0*cos(CLat_0r_u1)^2)^2 + sin_Th^2*(1.0-RI/ans)^2)
   g11_cov(ii,jj,Num_u3-1-kk)	= ans^4/(RI^2*cos(CLat_0r_u1)^4*one_p_3cos_sq^2)*((1.0-RI/ans)^2 + 0.25*(cos_Th/sin_Th)^2*(1.0+3.0*cos(CLat_0r_u1)^2)^2)
   g13_cov(ii,jj,Num_u3-1-kk)	= ans^4*cos_Th/(2.0*RI^2*cos(CLat_0r_u1)*one_p_3cos_sq)
   g22_cov(ii,jj,Num_u3-1-kk)	= ans*ans*sin_Th*sin_Th
   g33_cov(ii,jj,Num_u3-1-kk)	= ans^6*cos(CLat_0r_u1)^2/(RI^4*one_p_3cos_sq)
   Jacb(ii,jj,Num_u3-1-kk)		= ans^6*cos(CLat_0r_u1)/(RI^3*one_p_3cos_sq)
; assumes Re=1
   h_nu(ii,jj,Num_u3-1-kk)		= Ans^2/(sin_Th*sqrt(one_p_3cos_sq))   ; At the Equator this is in the R direction
   h_ph(ii,jj,Num_u3-1-kk)		= Ans*sin_Th                          ; E_W, Azimuthal
   h_mu(ii,jj,Num_u3-1-kk)		= Ans^3/(sqrt(one_p_3cos_sq))         ; At the Equator this is in the Theta direction (FA)
  end
 end

;  Do Equator, located at Num_u3_2 position
 For ii=0,Num_u1-1 do $
 Begin
  CLat_0_u1				= CLat_Max_u1-float(ii)*del_Th_0_u1            ; Co_Lat at Nth hemisphere Ionosphere (in degrees)
  CLat_0r_u1			= CLat_0_u1*!pi/180.0                         ; Convert CLat_0 to radians
  ThArr(ii,jj,Num_u3_2)	= !pi/2.0
  RArr(ii,jj,Num_u3_2)	= RI/sin(CLat_0r_u1)^2
  u1Arr(ii,jj,Num_u3_2)	= -sin(CLat_0r_u1)^2  ; CHECK
; u3 should all be zero along the equator
  Theta=!Dpi/2                                			; Change Theta here for Equator
  ThArr(ii,jj,Num_u3_2)	= Theta
  Lat					= !pi/2.0-Theta
  cos_Th				= cos(Theta)
  one_p_3cos_sq			= 1.0+3.0*cos_Th*cos_Th
  sin_Th				= sin(Theta)
  ans 					= Rarr(ii,jj,Num_u3_2)
  g11_con(ii,jj,Num_u3_2)	= (RI^2/ans^4)*sin_Th^2*one_p_3cos_sq
  g13_con(ii,jj,Num_u3_2)	= -(RI^3/ans^5)*sin(CLat_0r_u1)^2*cos_Th*one_p_3cos_sq/(2.0*cos(CLat_0r_u1)^3)
  g22_con(ii,jj,Num_u3_2)	= 1.0/(ans^2*sin_Th^2)
  g33_con(ii,jj,Num_u3_2)	= RI^4/(ans^6*cos(CLat_0r_u1)^6)*(0.25*cos_Th^2*(1.0+3.0*cos(CLat_0r_u1)^2)^2 + sin_Th^2*(1.0-RI/ans)^2)
  g11_cov(ii,jj,Num_u3_2)	= ans^4/(RI^2*cos(CLat_0r_u1)^4*one_p_3cos_sq^2)*((1.0-RI/ans)^2 + 0.25*(cos_Th/sin_Th)^2*(1.0+3.0*cos(CLat_0r_u1)^2)^2)
  g13_cov(ii,jj,Num_u3_2)	= ans^4*cos_Th/(2.0*RI^2*cos(CLat_0r_u1)*one_p_3cos_sq)
  g22_cov(ii,jj,Num_u3_2)	= ans*ans*sin_Th*sin_Th
  g33_cov(ii,jj,Num_u3_2)	= ans^6*cos(CLat_0r_u1)^2/(RI^4*one_p_3cos_sq)
  Jacb(ii,jj,Num_u3_2)		= ans^6*cos(CLat_0r_u1)/(RI^3*one_p_3cos_sq)
; assumes Re=1
   h_nu(ii,jj,Num_u3_2)		= ans^2/(sin_Th*sqrt(one_p_3cos_sq))   	; At the Equator this is in the R direction
   h_ph(ii,jj,Num_u3_2)		= ans*sin_Th                          	; E_W, Azimuthal
   h_mu(ii,jj,Num_u3_2)		= ans^3/(sqrt(one_p_3cos_sq))         	; At the Equator this is in the Theta direction (FA)
 end
 end

 PhiArr = 2.0*!Dpi/(Num_u2)*Dindgen(num_u2)					; Phi = 0 is at Midday
 u2Arr	= PhiArr


;	Spherical Coordinate for Grid Points

 For ii=0,Num_u1-1 do $
 Begin
 For jj=0,Num_u2-1 do $
 Begin
 For kk=0,Num_u3-1 do $
 Begin

 	xArr(ii,jj,kk)	= Rarr(ii,jj,kk)*sin(ThArr(ii,jj,kk))*cos(PhiArr(jj))
 	yArr(ii,jj,kk)	= Rarr(ii,jj,kk)*sin(ThArr(ii,jj,kk))*sin(PhiArr(jj))
	zArr(ii,jj,kk)	= Rarr(ii,jj,kk)*cos(ThArr(ii,jj,kk))
 end
 end
 end

;	<DAVE> Needed to general the spherical coordinate scale functions (see Lysak 2004 appendix) so added this section to my orginal code
;	Compare the expression for h_mu ( = sqrt(g_33) ) at theta_0 to h_r at theta_0 (equation A6 in Lysak 04)
;	I also just noticed my expression for h_mu assumes Ri = 1 so my scaling into dipole coordinates is off near the ionospheres... I'll correct that!

;	Generate spherical scale function at Ionosphere for use in ionospheric BC
h_ra_n = Dblarr(Num_u1,Num_u2)				;	Note these are different h_nu, h_ph, h_th
h_ph_n = Dblarr(Num_u1,Num_u2)
h_th_n = Dblarr(Num_u1,Num_u2)

h_ra_s = Dblarr(Num_u1,Num_u2)
h_ph_s = Dblarr(Num_u1,Num_u2)
h_th_s = Dblarr(Num_u1,Num_u2)

For jj = 0,Num_u2-1 do $
begin
 For ii = 0,Num_u1-1 do $
 begin
  ; Northern Hemipshere											;	These where coded in my version of the IDL 2D code so just copied them
  Theta_0	 	=  Colat_N(ii)
  h_ph_n(ii,jj)	=  Ri*sin(Theta_0)
  h_th_n(ii,jj)	= -Ri/(2.d0*sin(Theta_0)*cos(Theta_0))
  h_ra_n(ii,jj)	= -(2.d0*Ri*(cos(Theta_0))^2)/(1.d0+3.d0*(cos(theta_0))^2)

  ; Southern Hemipshere
  h_ph_s(ii,jj)	=  Ri*sin(Theta_0)
  h_th_s(ii,jj)	=  Ri/(2.d0*sin(Theta_0)*cos(Theta_0))
  h_ra_s(ii,jj)	=  (2.d0*Ri*(cos(Theta_0))^2)/(1.d0+3.d0*(cos(theta_0))^2)
  end
end





end

; ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** *****

; ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** *****
;
Pro New_r,r0,u1,u3,Th_0r,RI,ans
;	Newtons Rule Procedure
; 	Solve for R using Newton's method
 MaxN=100
 Tol=1.e-3
 err=1.0
 ii=0
 While (ii lt MaxN) AND (err gt Tol) do $
 Begin
  fr=u3*u3*r0^4*cos(Th_0r)*cos(Th_0r)/RI^3 - u1*r0 - RI
  df_dr=4.0*u3*u3*r0^3*cos(Th_0r)*cos(Th_0r)/RI^3 - u1
  ans=r0-fr/df_dr
  err=abs(ans-r0)
  r0=ans
 end
end
;
; ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** *****
;
; ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** *****
;
Pro gen_Va,		Lvarr,Num_u1,Num_u2,Num_u3,Rarr,ThArr,PhiArr,Xarr,rho_a,Va

Va		= DblArr(Num_u1,Num_u2,Num_u3)
Rho_a	= DblArr(Num_u1,Num_u2,Num_u3)


 For jj=0,Num_u2-1 do $        										; chooses field line
 Begin
 For ii=0,Num_u1-1 do $        										; chooses field line
 Begin

 Lshell = LVarr(ii)

; Calc of B_0

 Re		= 6378000.0           											; Earth Radii in m
 Lat	= 0.0 				 											; Set LAtitude at equator
 K_0	= 8.0e15              											; Earth magnetic moment
 Beq_0	= sqrt(4.0-3.0*cos(Lat)^2)/(cos(Lat)^6) * K_0/(Lshell*Re)^3


;	Select a Equatorial Va profile (Differing Plasmapause profiles)

; tanhm=4.0*(Lshell-5.6)										;	Plasmapause at L = 5.6
; tterm=(157./Lshell^2.6)*(tanh(tanhm)+1.0)/1.6


 tanhm=2.0*(Lshell-4.5)											;	Plasmapause at L = 4.5 (Used for 1st 2D Model Paper)
 tterm=(55.0/Lshell^2.6)*(tanh(tanhm)+1.0)/2.0

; tanhm=3.0*(Lshell-3.5)										;	Plasmapause at L = 3.5
; tterm=(27.0/Lshell^2.6)*(tanh(tanhm)+1.0)/1.6


  va_equator	= ((2.09/Lshell^1.1+tterm)*990.0-50.0)*1.0e3

  PowL		= 4.0												; Power Law index on Va Profile

  nH		= Beq_0^2/((4.0*!Dpi*1.0e-7)*va_equator^2)      	; Calc Hydrogen plasma number density at the equator
  nO		= 4.0e7												; Oxygen density (O2 Kg/Re^3) at ionosphere
  nO_Ht		= 250e3/Re											; Scale Height of Oxygen

  EnH		= 0.0*1.6e-21										; Density Enhancement for Plume (SI units)
  Posn		= 6.5												; Location of plume in equatorial plane (in Re)
  Radn		= 1.0												; Half width of plume

  For kk=0,Num_u3-1 do $       									; Points along field line
  Begin
   	Lat=!pi/2.0-ThArr(ii,jj,kk)

	bmg = sqrt(4.0-3.0*cos(Lat)^2)/(cos(Lat)^6) * K_0/(Lshell*Re)^3

   	rho_a(ii,jj,kk)=nH*(Lshell/RArr(ii,jj,kk))^PowL + nO/Re^3*Exp(-(Rarr(ii,jj,kk)-1.0)/nO_Ht) + $
   				(Lshell/RArr(ii,jj,kk))^PowL*EnH* EXP(-((Rarr(ii,jj,Fix(Num_u3/2))-Posn)^2)^2/Radn)^2

   	Va(ii,jj,kk) = bmg/sqrt((4.0*!Dpi*1.0e-7)*rho_a(ii,jj,kk))/Re

  end
end
end


End					; End of Va profile

; ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** *****
;
; ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** *****

;	Calculate Courant condition for Model
Pro Courant,	Xarr,Zarr,Va,Dt,Nt,Num_u1,Num_u2,Num_u3,Delta_L,Trans_T,N_sec,Eperp,V2,Sig_a


 Delta_L 	= DblArr(num_u1,Num_u3-1)
 Delta_T  	= DblArr(num_u1,Num_u3-1)
 Trans_T	= DblArr(num_u1)

 Nt			= long(0)

 For ii = 0, Num_u1-1 do $
 begin
  For kk = 0, Num_u3-2 do $
  begin
   Delta_x 			= (xArr(ii,0,kk+1) - xArr(ii,0,kk))			; 	Change in X coordinate along field line
   Delta_y 			= (zArr(ii,0,kk+1) - zArr(ii,0,kk))			; 	Change in Y coordinate along field line
   Delta_L(ii,kk)	= sqrt(Delta_x^2 + Delta_y^2)					; 	Line segment Length (in Re)

   Av_Va			= max([Va(ii,0,kk+1),Va(ii,0,kk)])

   Delta_T(ii,kk)	= Delta_L(ii,kk)/Av_Va							;	Transit time across Delta_L segment
   Trans_T(ii)		= Trans_T(ii) + Delta_T(ii,kk)					;	Transit time along field line

  end
 end

	Dt 				= 0.99*min(Delta_T)								; Maximum allowable Time step (Courant Condition)
	Nt				= Long((N_sec/(2.0*Dt))+2)								; Number of iterations for this run

 print,'Courant Condition =',min(Delta_T)
 print,'Dt =',Dt
 print,'# of interations per second = ',1.0/(2.0*Dt)
 print,'Total Number of interations = ',Nt

 Re			= 6370.0e3
  c			= 3.0e8/Re										;	Scaled speed of light (in Re per sec)
 u0			= 4.0*!Dpi*1.0e-7/Re							;	Scaled Permablity
 e0			= 8.85e-12*Re^3									;	Scaled Permativity


 Eperp		= DblArr(Num_u1,Num_u2,Num_u3)
 V2			= DblArr(Num_u1,Num_u2,Num_u3)
 Sig_a		= DblArr(Num_u1,Num_u2,Num_u3)

 For ii = 0,Num_u1-1 do $
 Begin
  For jj = 0,Num_u2-1 do $
  Begin
  For kk = 0,Num_u3-1 do $
  Begin
   Eperp(ii,jj,kk)		=	e0*(1+c^2/Va(ii,jj,kk)^2)
   V2(ii,jj,kk)			=	1.0/(u0*Eperp(ii,jj,kk))
   Sig_a(ii,jj,kk)		= 	1.0/(u0*Va(ii,jj,kk))
  end
  end
 end



end			; End of Courant calculation

; ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** *****
;
; ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** *****


Pro IO_Parm,	Xarr,Yarr,Zarr,Num_u1,Num_u2,Num_u3,Delta_L,Trans_T,rho_a,Va,Plot_Dir,Parm_file,Iono_file,$
 				g11_con,g13_con,g22_con,g33_con, g11_cov,g13_cov,g22_cov,g33_cov,N_sec,Eperp,V2,Sig_a,$
            	Jacb,Colat_N,Colat_S,h_nu,h_ph,h_mu,LMin_u1,LMax_u1, LMin_u3, LMax_u3,RArr,ThArr,PhiArr,$
            	LVArr,Sd_S,Sp_S,Sh_S,Sd_N,Sp_N,Sh_N,SHF,SHF_dr,Alpha,Alpha_dr,Zeros,$
   	           	h_ra_n,h_th_n,h_ph_n,h_ra_s,h_th_s,h_ph_s


Xord 	= Reform(Xarr(*,0,*))
Zord 	= Reform(Zarr(*,0,*))

Xord_Eq = DblArr(Num_u1,Num_u2+1)
Yord_Eq = DblArr(Num_u1,Num_u2+1)
Va_Eq 	= DblArr(Num_u1,Num_u2+1)
Rho_Eq 	= DblArr(Num_u1,Num_u2+1)

Xord_Eq(*,0:Num_u2-1) 	= Reform(Xarr(*,*,Fix(Num_u3/2)))
Yord_Eq(*,0:Num_u2-1) 	= Reform(Yarr(*,*,Fix(Num_u3/2)))

Xord_Eq(*,Num_u2)		= Xarr(*,0,Fix(Num_u3/2))
Yord_Eq(*,Num_u2)		= Yarr(*,0,Fix(Num_u3/2))

Va_Eq(*,0:Num_u2-1)		= Va(*,*,Num_u3/2)
Va_Eq(*,Num_u2)			= Va(*,0,Num_u3/2)

Rho_eq(*,0:Num_u2-1)	= Rho_a(*,*,Num_u3/2)
Rho_eq(*,Num_u2)		= Rho_a(*,0,Num_u3/2)


Device,decompose=0
 LoadCT,13
 TVLCT,r,g,b,/get
 r(0)	= 255	& g(0)	=255 	& b(0)	=255
 r(255)	=0		& g(255)=0 		& b(255)=0
 TVLCT,r,g,b

Window,11,Title='Parameters 3D Model',Xsize=650,Ysize=800
!P.charsize=2.0
!P.multi=[0,2,4]


contour,Reform(Va(*,0,*)),Xord,Zord,nlevels=255,/fill,/isotropic,$
		Title='V!dA!n',Xtitle='X ( R!dE!n)',Ytitle='Z ( R!dE!n)'

contour,Reform(Alog10(rho_a(*,0,*))),Xord,Zord,nlevels=255,/fill,/isotropic,$
		Title='Log!d10!n(!7q!3)',Xtitle='X ( R!dE!n)',Ytitle='Z ( R!dE!n)'

contour,Va_Eq,Xord_Eq,Yord_Eq,nlevels=255,/fill,/isotropic,$
		Title='V!dA!n',Xtitle='X ( R!dE!n)',Ytitle='Y ( R!dE!n)'

contour,Alog10(Rho_eq),Xord_Eq,Yord_Eq,nlevels=255,/fill,/isotropic,$
		Title='Log!d10!n(!7q!3)',Xtitle='X ( R!dE!n)',Ytitle='Y ( R!dE!n)'

Plot,Xord,Zord,Psym=3,symsize=7,/isotropic,$
		Title='Grid (X Z plane Y = 0)',Xtitle='X ( R!dE!n)',Ytitle='Z ( R!dE!n)'

Plot,Xord_Eq,Yord_Eq,Psym=3,symsize=7,/isotropic,$
		Title='Grid (X Y plane Z = 0)',Xtitle='X ( R!dE!n)',Ytitle='Y ( R!dE!n)'


Plot,Xarr(*,0,Fix(Num_u3/2)),Trans_T,/Ylog,$
		Title='Transit Time along B!d0!n',Xtitle='L shell ( R!dE!n)',Ytitle='T (sec)'

Plot,Xarr(*,0,Fix(Num_u3/2)),(2.0*!Dpi*Xarr(*,0,Fix(Num_u3/2)))/Va(*,Fix(Num_u3/2)),/Ylog,$
		Title='Travel Time Round Earth (Equator)',Xtitle='L shell ( R!dE!n)',Ytitle='T (sec)'



 Filename = Plot_dir+Parm_file+'.PNG'
 image = TVRD(0,0,!d.X_size,!d.Y_size,true=1)
 Write_PNG,filename,image,r,g,b
 Print,'File written to ', Filename


 Window,12,title='Conductivity',Xsize=450,Ysize=900
 !P.multi=[0,1,3]
 !P.charsize=2.0

  Re = 6370e3

  Nlvls=254

  Sdmax	= max([max(Sd_S),max(Sd_N)])/Re^2
  LvlSd = 1.01*Sdmax*Findgen(Nlvls)/(Nlvls)
  contour,transpose(Sd_S/Re^2),PhiArr*180/!Pi,(Reform(ThArr(*,0,Num_u3-1)))*180/!Pi,/fill,levels=LvlSd,$
   Title='!7R!3!dD!n Max = '+String(Sdmax,Format='(F9.1)')+' S',$
   Xtitle='Longitiude (deg)',Ytitle='Co-Latitude (deg)',$
   Xstyle=1,Ystyle=1,Xrange=[min(PhiArr*180/!Pi),max(PhiArr*180/!Pi)],Yrange=[180.0,0.0]
  contour,transpose(Sd_N/Re^2),PhiArr*180/!Pi,(Reform(ThArr(*,0,0)))*180/!Pi,/fill,levels=LvlSd,$
   Xstyle=1,Ystyle=1,Xrange=[min(PhiArr*180/!Pi),max(PhiArr*180/!Pi)],Yrange=[180.0,0.0],/overplot

  Smax	= max([max(Sp_S),max(Sp_N),max(Sh_S),max(Sh_N)])/Re^2
  Spmax	= max([max(Sp_S),max(Sp_N)])/Re^2
  Shmax	= max([max(Sh_S),max(Sh_N)])/Re^2

  LvlS = 1.01*Smax*Findgen(Nlvls)/(Nlvls)
  contour,transpose(Sp_S/Re^2),PhiArr*180/!Pi,(Reform(ThArr(*,0,Num_u3-1)))*180/!Pi,/fill,levels=LvlS,$
   Title='!7R!3!dP!n  (Max = '+String(Spmax,Format='(F6.2)')+' S)',$
   Xtitle='Longitiude (deg)',Ytitle='Co-Latitude (deg)',$
   Xstyle=1,Ystyle=1,Xrange=[min(PhiArr*180/!Pi),max(PhiArr*180/!Pi)],Yrange=[180.0,0.0]
  contour,transpose(Sp_N/Re^2),PhiArr*180/!Pi,(Reform(ThArr(*,0,0)))*180/!Pi,/fill,levels=LvlS,$
   Xstyle=1,Ystyle=1,Xrange=[min(PhiArr*180/!Pi),max(PhiArr*180/!Pi)],Yrange=[180.0,0.0],/overplot

  contour,transpose(Sh_S/Re^2),PhiArr*180/!Pi,(Reform(ThArr(*,0,Num_u3-1)))*180/!Pi,/fill,levels=lvlS,$
   Title='!7R!3!dH!n (Max = '+String(Shmax,Format='(F6.2)')+' S)',$
   Xtitle='Longitiude (deg)',Ytitle='Co-Latitude (deg)',$
   Xstyle=1,Ystyle=1,Xrange=[min(PhiArr*180/!Pi),max(PhiArr*180/!Pi)],Yrange=[180.0,0.0]
  contour,transpose(Sh_N/Re^2),PhiArr*180/!Pi,(Reform(ThArr(*,0,0)))*180/!Pi,/fill,levels=lvlS,$
   Xstyle=1,Ystyle=1,Xrange=[min(PhiArr*180/!Pi),max(PhiArr*180/!Pi)],Yrange=[180.0,0.0],/overplot



 Filename = Plot_dir+Iono_file+'.PNG'
 image = TVRD(0,0,!d.X_size,!d.Y_size,true=1)
 Write_PNG,filename,image,r,g,b
 Print,'File written to ', Filename



 Filename1 = Plot_dir+Parm_file+'.sav'
 Save,Filename=Filename1,	Xarr,Yarr,Zarr,Num_u1,Num_u2,Num_u3,Delta_L,Trans_T,rho_a,Va,Plot_Dir,Parm_file,$
 							g11_con,g13_con,g22_con,g33_con, g11_cov,g13_cov,g22_cov,g33_cov,N_sec,Eperp,V2,Sig_a,$
            				Jacb,Colat_N,Colat_S,h_nu,h_ph,h_mu,LMin_u1,LMax_u1, LMin_u3, LMax_u3,RArr,ThArr,PhiArr,$
           					LVArr,Sd_S,Sp_S,Sh_S,Sd_N,Sp_N,Sh_N,SHF,SHF_dr,Alpha,Alpha_dr,Zeros,$
           					h_ra_n,h_th_n,h_ph_n,h_ra_s,h_th_s,h_ph_s



End					;	 End of Output of Paramters

; ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** *****
;
; ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** *****

pro Iterate,	Num_u1,Num_u2,Num_u3,g11_cov,g13_cov,g22_cov,g33_cov,g11_con,jacb,u1Arr,u2Arr,u3Arr,V2,Dt,Nt,$
				Xarr,Yarr,Zarr,Rarr,ThArr,PhiArr,Plot_Dir,Plot_File,h_nu,h_ph,h_mu,SHF,SHF_dr,Alpha,Alpha_dr,$
				del_Th,u0,Colat_N,Colat_S,INVSS_S,INVSS_N,Ri,h_ra_n,h_th_n,h_ph_n,h_ra_s,h_th_s,h_ph_s


Ct	= 0

;	Spatial Extent of Driver
Dspatial	= DblArr(Num_u2,Num_u3)
Driver_S	= DblArr(Num_u2,Num_u3)

For jj = 0,Num_u2-2 do $
begin
For kk = 1,Num_u3-2 do $
begin

D_th	= sqrt((Xarr(Num_u1-1,jj,kk) - Xarr(Num_u1-1,0,Num_u3/2))^2+(Zarr(Num_u1-1,jj,kk) - Zarr(Num_u1-1,0,Num_u3/2))^2)
;D_ph	= sqrt((ABS(Yarr(Num_u1-1,jj,kk)) - Yarr(Num_u1-1,Num_u2/4,Num_u3/2))^2 )

D_th	= (Rarr(Num_u1-1,jj,kk)-max(Rarr))

Hwidth	= 0.5
;D0		= sqrt((D_th^2+D_ph^2)/Hwidth^2)
;
D0		= sqrt(D_th^2)/Hwidth
Ph0		= sqrt((!Dpi-Phiarr(jj))^2)/Hwidth

;Dspatial(jj,kk) = Exp(-D0^2)*sin(2.0*PhiArr(jj))
Dspatial(jj,kk) = Exp(-D0^2)*Exp(-Ph0^2)

end
end


For kk = 0,Num_u3/2 do $
begin
For jj = 0,Num_u2/2-1 do $
begin
 Driver_S(2*jj,2*kk) = Dspatial(2*jj,2*kk)/max(Dspatial)	;	place driver at B3 nodes on outer Boundary + Normalise
end
end


D_u1	= Shift(u1Arr,-1,0,0) - Shift(u1Arr,1,0,0)		; Spacing Difference in u1
D_u2	= 2.0*!Dpi/(Num_u2/2) ;Shift(u2Arr_3D,0,0,-1) - Shift(u2Arr_3D,0,0,1)		; Spacing Difference in u2
D_u3	= Shift(u3Arr,0,0,-1) - Shift(u3Arr,0,0,1)		; Spacing Difference in u3



 E1_loc	= Dblarr(Num_u1,Num_u2)
 E2_loc	= Dblarr(Num_u1,Num_u2)

 For jj = 0,Num_u2/2-1 do $
 begin
 For ii = 0,Num_u1/2-1 do $
 begin
 E1_loc(2*ii+1,2*jj+1)	= 1.0			;	Used in Ionopsheric Boundary condition for E's location on ionosphere
 E2_loc(2*ii,2*jj)		= 1.0
 end
 end

 B1_con		= Dblarr(Num_u1,Num_u2,Num_u3)
 B2_con		= Dblarr(Num_u1,Num_u2,Num_u3)
 B3_con		= Dblarr(Num_u1,Num_u2,Num_u3)
 E1_con		= Dblarr(Num_u1,Num_u2,Num_u3)
 E2_con		= Dblarr(Num_u1,Num_u2,Num_u3)

 B1_cov		= Dblarr(Num_u1,Num_u2,Num_u3)
 B2_cov		= Dblarr(Num_u1,Num_u2,Num_u3)
 B3_cov		= Dblarr(Num_u1,Num_u2,Num_u3)
 E1_cov		= Dblarr(Num_u1,Num_u2,Num_u3)
 E2_cov		= Dblarr(Num_u1,Num_u2,Num_u3)

 B3_N_C 	= DblArr(Num_u1*Num_u2)
 B3_S_C 	= DblArr(Num_u1*Num_u2)

 E1_N		= Dblarr(Num_u1,Num_u2)
 E2_N		= Dblarr(Num_u1,Num_u2)
 E1_S		= Dblarr(Num_u1,Num_u2)
 E2_S		= Dblarr(Num_u1,Num_u2)

 DPsi_dr_N	= Dblarr(Num_u1,Num_u2)
 Psi_N		= Dblarr(Num_u1,Num_u2)
 DPsi_dr_S	= Dblarr(Num_u1,Num_u2)
 Psi_S		= Dblarr(Num_u1,Num_u2)



	Norm 		= 10.0							; Normalisation for Time driver
	Amp			= 1.0							; Amplitude of Driver

;	If continuing from restore file
;Number = 20										; Store New number of times (for total time)
;Restore,filename=Plot_Dir+'ReRun_'+StrTrim(Fix(Number),2)+'.sav'			; Restoring sav file
;tt1 = tt										; Storing sav files number of iterations

Starttime = SYSTIME(/SECONDS)

For tt = Long(0),Long(Nt)-1 do $
begin

;If tt eq 0 then tt = tt1						; Restoring iteration number for this run first time threw (continuing from sav file)

Time_Iter = Float(tt)*Dt*2.0

;	Time Dependence of driver
  	T_dep		= 1000*Time_Iter/Norm*(1.0-Time_Iter/Norm)*Exp(-5.0*Time_Iter/Norm)

;	Constructing Driver on outer boundary (Space and Time variations)
	OuterL 					=  Driver_S*T_dep*Amp
;	B3_cov(Num_u1-1,*,*) 	=  OuterL

;	Start Iteration
;	Start E field iteration

  E1_con		= 2.d0*Dt*V2/Jacb*$
  				  [(Shift(B3_cov,0,-1,0) - Shift(B3_cov,0,1,0))/D_u2 - (Shift(B2_cov,0,0,-1) - Shift(B2_cov,0,0,1))/D_u3] +  E1_con

  E2_con		= 2.d0*Dt*V2/Jacb*$
  				  [(Shift(B1_cov,0,0,-1) - Shift(B1_cov,0,0,1))/D_u3 - (Shift(B3_cov,-1,0,0) - Shift(B3_cov,1,0,0))/D_u1] +  E2_con


  E1_con(0,*,*) 		= 0.d0				; Cleanup after shifts at inner L shell
  E1_con(*,*,0) 		= 0.d0				; Cleanup after shifts at Northern Ionopshere
  E1_con(*,*,Num_u3-1)	= 0.d0				; Cleanup after shifts at Southern Ionopshere

  E2_con(0,*,*) 		= 0.d0				; Cleanup after shifts at inner L shell
  E2_con(*,*,0) 		= 0.d0				; Cleanup after shifts at Northern Ionopshere
  E2_con(*,*,Num_u3-1)	= 0.d0				; Cleanup after shifts at Southern Ionopshere


;	Evolve E_con into E_cov

  E1_cov		= E1_con/g11_con

  E2_cov		= E2_con*g22_cov


;	Ionospheric Boundary Conditons

;	Northern Ionosphere
;	~~~~~~~~~~~~~~~~~~~~

  E1_cov(*,*,0) 		= E1_N						; Northern ionospheric Boundary
  E2_cov(*,*,0) 		= E2_N

;	Southern Ionosphere
;	~~~~~~~~~~~~~~~~~~~~

  E1_cov(*,*,Num_u3-1)	= E1_S						; Southern ionospheric Boundary
  E2_cov(*,*,Num_u3-1)	= E2_S

;	Inner L shell Boundary
;	~~~~~~~~~~~~~~~~~~~~~~
  E1_cov(0,*,*) 		= 0.0				; Perfectly Reflecting
  E2_cov(0,*,*) 		= 0.0				; Perfectly Reflecting


;	Start B field iteration

    B1_con	=   2.d0*Dt/Jacb*[(Shift(E2_cov,0,0,-1) - Shift(E2_cov,0,0,1))/D_u3 ] + B1_con

    B2_con	=  -2.d0*Dt/Jacb*[(Shift(E1_cov,0,0,-1) - Shift(E1_cov,0,0,1))/D_u3 ] + B2_con

    B3_con	=   2.d0*Dt/Jacb*[(Shift(E1_cov,0,-1,0) - Shift(E1_cov,0,1,0))/D_u2 - (Shift(E2_cov,-1,0,0) - Shift(E2_cov,1,0,0))/D_u1] + B3_con

;	Tidy up after shifts B3_con at corner between Ionosphere and outer L ?????????????????

	B3_con(Num_u1-1,*,0)		= 4.d0/3.d0*B3_con(Num_u1-3,*,0)        - 1.d0/3.d0*B3_con(Num_u1-5,*,0)		; ?
	B3_con(Num_u1-1,*,Num_u3-1)	= 4.d0/3.d0*B3_con(Num_u1-3,*,Num_u3-1) - 1.d0/3.d0*B3_con(Num_u1-5,*,Num_u3-1)		; ?


	B3_con(0,*,*) = 0.0					;	Inner L shell Boundary


;	Evolve B_con into B_cov

;	Evolve B1

   	Av_B3 	= ( Shift(B3_con, 1,0,1) + Shift(B3_con, 1,0,-1)+ $
   				Shift(B3_con,-1,0,1) + Shift(B3_con,-1,0,-1)  )/4.0

   	Av_B3(*,*,0)	= 0.0	&  Av_B3(*,*,Num_u3-1)	= 0.0	& Av_B3(0,*,*) = 0.0

   	B1_cov 	= B1_con*g11_cov + Av_B3*g13_cov

   	B1_cov(*,*,0) 			= 0.0   ;	Clean up after shift at the ends of Arrays  (Northern ionosphere)
   	B1_cov(*,*,num_u3-1) 	= 0.0   ;	Clean up after shift at the ends of Arrays  (Southern ionosphere)
   	B1_cov(0,*,*) 			= 0.0   ;	Clean up after shift at the ends of Arrays  (inner L shell)

;	Evolve B2

   	B2_cov 	= B2_con*g22_cov


;	Evolve B3_cov
   	Av_B1 	= (Shift(B1_con,1,0,1) + Shift(B1_con,-1,0,1) + Shift(B1_con,1,0,-1) + Shift(B1_con,-1,0,-1))/ 4.d0
   	Av_B1(*,*,0)	= 0.0 		& 	Av_B1(*,*,Num_u3-1)	= 0.0			; Clean up after shifts due to odd number of points
   	B3_cov  = Av_B1*g13_cov + B3_con*g33_cov

	Av_B1	= (Shift(B1_con,1,0,-1)+Shift(B1_con,-1,0,-1))/2.0
   	B3_cov(*,*,0)  = Av_B1(*,*,0)*g13_cov(*,*,0) + B3_con(*,*,0)*g33_cov(*,*,0)		;	Northern Ionosphere

	Av_B1	= (Shift(B1_con,1,0,1)+Shift(B1_con,-1,0,1))/2.0
   	B3_cov(*,*,Num_u3-1)  = Av_B1(*,*,Num_u3-1)*g13_cov(*,*,Num_u3-1) + B3_con(*,*,Num_u3-1)*g33_cov(*,*,Num_u3-1)		;	Southern Ionosphere

;	Boundary Conditons @@@
;	Along Outer L shell Boundary
	B3_cov(Num_u1-1,*,*) 	= OuterL						; This is the driver function

;	Along Inner L shell Boundary
;    Av_B1	= (Shift(B1_con,-1,0,-1)+Shift(B1_con,-1,0,1))/2.0
	B3_cov(0,*,*) 			= 0.0

;	Ionospheric Boundary Condition
;	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

;	FOR NORTHERN HEMISPHERE
  	B3_N		= B3_con(*,*,0)*h_ra_n									;	<DAVE> I orginally had h_mu(*,*,0) here which is different to h_ra_n
 	B3_N		= (Shift(B3_N,1,0)+Shift(B3_N,-1,0) )/2.0 + B3_N		;	Linear Iterpolate in U1 direction
	B3_N		= (Shift(B3_N,0,1)+Shift(B3_N, 0,-1) )/2.0 + B3_N
; 	B3_N(0,*)	= 0.0
 	B3_N(0,*)	= 4.0/3.0*B3_N(1,*) - 1.0/3.0*B3_N(2,*)												; d1 B3 = 0

;	Converting B3_N into a column vector for the LU decomposition
 	For jj	= 0,Num_u2-1 do $
 	begin
 	For ii	= 0,Num_u1-1 do $
 	begin
	B3_N_C(ii+jj*Num_u1)	=  B3_N(ii,jj)
 	end
 	end

;	<DAVE> Although coded to use the LU_complex routine the arrays and answer are actually just real (I've checked),
;	so probable should change to a routine which does the SVD only with real arrays

  	SphdotB_N		= Conj(SHF_dr) ## B3_N_C
  	SphdotB_N		= transpose(SphdotB_N)
 	Coeffs_B3_N		= LU_Complex(Alpha_dr,SphdotB_N,/double)
 	DPsi_dr_N_C 	= Transpose(SHF_dr)##Coeffs_B3_N
	Psi_N_C			= Transpose(SHF   )##Coeffs_B3_N

;	Converting Psi and dPsi/dr into 2D arrays for the LU decomposition
 	For jj	= 0,Num_u2-1 do $
 	begin
 	For ii	= 0,Num_u1-1 do $
 	begin
	Psi_N(ii,jj)	=  Psi_N_C(ii+jj*Num_u1)
	DPsi_dr_N(ii,jj)= DPsi_dr_N_C(ii+jj*Num_u1)
 	end
 	end

;	Calculate the Components of B below the ionosphere
;	B_theta
	DPsi_dTh_N				= (Shift(Psi_N,-1,0) - Shift(Psi_N,1,0))/(-2.d0*del_th)		; B_theta in Atmosphere
;	Tide up end point after shift
    DPsi_dTh_N(0,*)   		= (-3.0*Psi_N(0,*)        + 4.d0*Psi_N(1,*)        - Psi_N(2,*))/(-2.d0*del_th)
  	DPsi_dTh_N(Num_u1-1,*) 	= ( 3.0*Psi_N(Num_u1-1,*) - 4.d0*Psi_N(Num_u1-2,*) + Psi_N(Num_u1-3,*))/(-2.d0*del_th)
	DPsi_dTh_N 				= 1.0/Ri*DPsi_dTh_N					; Scaling for Spherical derivative in theta direction = 1/Ri

;	B_phi
	DPsi_dPh_N		= 1.0/h_ph_n*((Shift(Psi_N,0,-1) - Shift(Psi_N,0,1))/(D_u2))	; B_phi = 1/rsin(th) dPsi/dphi in Atmosphere

;;	Tide up end point after shift
;    DPsi_dPh_N(0,*)   		= (-3.0*Psi_N(0,*)        + 4.d0*Psi_N(1,*)        - Psi_N(2,*))/(-2.d0*del_th)
;  	DPsi_dPh_N(Num_u1-1,*) 	= ( 3.0*Psi_N(Num_u1-1,*) - 4.d0*Psi_N(Num_u1-2,*) + Psi_N(Num_u1-3,*))/(-2.d0*del_th)

; 	Interpolate B1 and B2 to ALL points just above the ionosphere from model solution to calculate currents in Ionosphere
;	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
;	<DAVE> I orginally had h_th(*,*,0) here which is different to h_th_n

	B1_N			= B1_cov(*,*,1)/h_th_n
	B1_N			= (Shift(B1_N,1,0) + Shift(B1_N,-1, 0))/2.0 + B1_N
	B1_N			= (Shift(B1_N,0,1) + Shift(B1_N, 0,-1))/2.0 + B1_N

;	interpolation using Second order Taylor expansion from end points at outer boundary
  	h				= CoLat_N(Num_u1-2)-CoLat_N(Num_u1-4) 																; step size along ionosphere is uniform on B1_N grid
  	FirD 			= ( B1_N(Num_u1-6,*) - 4.0*B1_N(Num_u1-4,*) + 3.0*B1_N(Num_u1-2,*))/(2.0*h)							; First Derivative Backwards difference O(h^2)
  	SecD 			= (-B1_N(Num_u1-8,*) + 4.0*B1_N(Num_u1-6,*) - 5.0*B1_N(Num_u1-4,*) + 2.0*B1_N(Num_u1-2,*))/(h^2) 	; Second Derivative Backwards difference O(h^2)
  	h1				= CoLat_N(Num_u1-2)-CoLat_N(Num_u1-1)																; Step size on B1_n_interp grid
  	B1_N(Num_u1-1,*)= B1_N(Num_u1-2,*) + h1*FirD + ((h1^2)/2.0)*SecD


	B2_N			= B2_cov(*,*,1)/h_ph_n								;	<DAVE> I orginally had h_ph(*,*,0) here which isn't much different to h_ph_n but I changed it anyway!
	B2_N			= (Shift(B2_N,1,0) + Shift(B2_N,-1, 0))/2.0 + B2_N  ;	Obviously did these changes to the southern ionosphere as well using the h_*_s arrays
	B2_N			= (Shift(B2_N,0,1) + Shift(B2_N, 0,-1))/2.0 + B2_N

;	interpolation using Second order Taylor expansion from end points at outer boundary
  	h				= CoLat_N(3)-CoLat_N(1) 																			; step size along ionosphere is uniform
  	FirD 			= (-B2_N(5,*) + 4.0*B2_N(3,*) - 3.0*B2_N(1,*))/(2.0*h)													; First Derivative Forwards difference O(h^2)
  	SecD 			= (-B2_N(7,*) + 4.0*B2_N(5,*) - 5.0*B2_N(3,*) + 2.0*B2_N(1,*))/(h^2) 										; Second Derivative Forwards difference O(h^2)
 	h1				= CoLat_N(0)-CoLat_N(1)
  	B2_N(0,*)		= B2_N(1,*) + h1*FirD + ((h1^2)/2.0)*SecD


;	Currents in the Northern Ionosphere
  	J_Th_N			= -(B2_N - DPsi_dPh_N)/u0
  	J_Ph_N			=  (B1_N - DPsi_dTh_N)/u0

;	Calculate electric fields in Ionosphere from discontinuity in B's for next time step...
 	E_Th_N			= (INVSS_N(*,*,0,0)*J_Th_N + INVSS_N(*,*,1,0)*J_Ph_N)
  	E_Ph_N			= (INVSS_N(*,*,0,1)*J_Th_N + INVSS_N(*,*,1,1)*J_Ph_N)

 	E1_N			= E_Th_N*h_th_n*E1_Loc			;	need to put on corect gridpoints
  	E2_N			= E_Ph_N*h_ph_n*E2_Loc


;	FOR SOUTHERN HEMISPHERE
  	B3_S		= B3_con(*,*,Num_u3-1)*h_ra_s
 	B3_S		= (Shift(B3_S,1,0)+Shift(B3_S,-1,0) )/2.0 + B3_S		;	Linear Iterpolate in U1 direction
	B3_S		= (Shift(B3_S,0,1)+Shift(B3_S, 0,-1) )/2.0 + B3_S
; 	B3_S(0,*)	= 0.0
 	B3_S(0,*)	= 4.0/3.0*B3_S(1,*) - 1.0/3.0*B3_S(2,*)													; d1 B3 = 0

;	Converting B3_S into a column vector for the LU decomposition
 	For jj	= 0,Num_u2-1 do $
 	begin
 	For ii	= 0,Num_u1-1 do $
 	begin
	B3_S_C(ii+jj*Num_u1)	=  B3_S(ii,jj)
 	end
 	end

  	SphdotB_S		= Conj(SHF_dr) ## B3_S_C
  	SphdotB_S		= transpose(SphdotB_S)
 	Coeffs_B3_S		= LU_Complex(Alpha_dr,SphdotB_S,/double)
 	DPsi_dr_S_C 	= Transpose(SHF_dr)##Coeffs_B3_S
	Psi_S_C			= Transpose(SHF   )##Coeffs_B3_S

;	Converting Psi and dPsi/dr into 2D arrays for the LU decomposition
 	For jj	= 0,Num_u2-1 do $
 	begin
 	For ii	= 0,Num_u1-1 do $
 	begin
	Psi_S(ii,jj)	=  Psi_S_C(ii+jj*Num_u1)
	DPsi_dr_S(ii,jj)= DPsi_dr_S_C(ii+jj*Num_u1)
 	end
 	end

;	Calculate the Components of B below the ionosphere
;	B_theta
	DPsi_dTh_S				= (Shift(Psi_S,-1,0) - Shift(Psi_S,1,0))/(2.d0*del_th)		; B_theta in Atmosphere
;	Tide up end point after shift
    DPsi_dTh_S(0,*)   		= (-3.0*Psi_S(0,*)        + 4.d0*Psi_S(1,*)        - Psi_S(2,*))/(2.d0*del_th)
  	DPsi_dTh_S(Num_u1-1,*) 	= ( 3.0*Psi_S(Num_u1-1,*) - 4.d0*Psi_S(Num_u1-2,*) + Psi_S(Num_u1-3,*))/(2.d0*del_th)
	DPsi_dTh_S				= 1.0/Ri*DPsi_dTh_S

;	B_phi
	DPsi_dPh_S				= 1.0/h_ph_s*((Shift(Psi_S,0,-1) - Shift(Psi_S,0,1))/(D_u2))	; B_phi = 1/rsin(th) dPsi/dphi in Atmosphere

;;	Tide up end point after shift
;    DPsi_dPh_S(0,*)   		= (-3.0*Psi_S(0,*)        + 4.d0*Psi_S(1,*)        - Psi_S(2,*))/(-2.d0*del_th)
;  	DPsi_dPh_S(Num_u1-1,*) 	= ( 3.0*Psi_S(Num_u1-1,*) - 4.d0*Psi_S(Num_u1-2,*) + Psi_S(Num_u1-3,*))/(-2.d0*del_th)

; 	Interpolate B1 and B2 to ALL points just above the ionosphere from model solution to calculate currents in Ionosphere
;	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	B1_S			= B1_cov(*,*,Num_u3-2)/h_th_s
	B1_S			= (Shift(B1_S,1,0) + Shift(B1_S,-1, 0))/2.0 + B1_S
	B1_S			= (Shift(B1_S,0,1) + Shift(B1_S, 0,-1))/2.0 + B1_S

;	interpolation using Second order Taylor expansion from end points at outer boundary
  	h				= CoLat_S(Num_u1-2)-CoLat_S(Num_u1-4) 																; step size along ionosphere is uniform on B1_S grid
  	FirD 			= ( B1_S(Num_u1-6,*) - 4.0*B1_S(Num_u1-4,*) + 3.0*B1_S(Num_u1-2,*))/(2.0*h)							; First Derivative Backwards difference O(h^2)
  	SecD 			= (-B1_S(Num_u1-8,*) + 4.0*B1_S(Num_u1-6,*) - 5.0*B1_S(Num_u1-4,*) + 2.0*B1_S(Num_u1-2,*))/(h^2) 	; Second Derivative Backwards difference O(h^2)
  	h1				= CoLat_S(Num_u1-2)-CoLat_S(Num_u1-1)																; Step size on B1_S_interp grid
  	B1_S(Num_u1-1,*)= B1_S(Num_u1-2,*) + h1*FirD + ((h1^2)/2.0)*SecD


	B2_S			= B2_cov(*,*,Num_u3-2)/h_ph_s
	B2_S			= (Shift(B2_S,1,0) + Shift(B2_S,-1, 0))/2.0 + B2_S
	B2_S			= (Shift(B2_S,0,1) + Shift(B2_S, 0,-1))/2.0 + B2_S

;	interpolation using Second order Taylor expansion from end points at outer boundary
  	h				= CoLat_S(3)-CoLat_S(1) 																			; step size along ionosphere is uniform
  	FirD 			= (-B2_S(5,*) + 4.0*B2_S(3,*) - 3.0*B2_S(1,*))/(2.0*h)													; First Derivative Forwards difference O(h^2)
  	SecD 			= (-B2_S(7,*) + 4.0*B2_S(5,*) - 5.0*B2_S(3,*) + 2.0*B2_S(1,*))/(h^2) 										; Second Derivative Forwards difference O(h^2)
 	h1				= CoLat_S(0)-CoLat_S(1)
  	B2_S(0,*)		= B2_S(1,*) + h1*FirD + ((h1^2)/2.0)*SecD

;	Currents in the Ionosphere
  	J_Th_S			= -(B2_S - DPsi_dPh_S)/u0
  	J_Ph_S			=  (B1_S - DPsi_dTh_S)/u0

;	Calculate electric fields in Ionosphere from discontinuity in B's for next time step...
  	E_Th_S			= (INVSS_S(*,*,0,0)*J_Th_S + INVSS_S(*,*,1,0)*J_Ph_S)
  	E_Ph_S			= (INVSS_S(*,*,0,1)*J_Th_S + INVSS_S(*,*,1,1)*J_Ph_S)

;	need to put on corect gridpoints
  	E1_S			= E_Th_S*h_th_s*E1_Loc
  	E2_S			= E_Ph_S*h_ph_s*E2_Loc


;	End of Ionospheric Boundary Calculation


 If (tt Mod Fix(1.0/(2.0*Dt))) EQ 0 then $
 begin

 OUTPUT,	Xarr, Yarr, Zarr, Num_u1, Num_u2, Num_u3, B1_con, B2_con, B3_cov, E1_con, E2_con,$
			tt,Time_Iter,Plot_dir,Plot_file,Ct,h_nu,h_ph,h_mu,Rarr,ThArr,PhiArr,$
			DPsi_dr_S,DPsi_dPh_S,DPsi_dTh_S,Psi_S,B1_S,B2_S,B3_S,E_Th_S,E_Ph_S,Coeffs_B3_S,$
			DPsi_dr_N,DPsi_dPh_N,DPsi_dTh_N,Psi_N,B1_N,B2_N,B3_N,E_Th_N,E_Ph_N,Coeffs_B3_N,Colat_N,Colat_S


 end


 end							; end of one time step (2 Dt) in iteration



 Finishtime = SYSTIME(/SECONDS)
 Print,'Average Computational Time per time iteration = ',(Finishtime - starttime)/Nt,' seconds'
 Print,'Total Iterations = ',Nt
 print,'Simulation Time  = ',Nt*(2.0*Dt),' seconds'
 Print,'Total Computational time = ',Finishtime - Starttime,' seconds'
 Print,'Memory required (Kb): ', MEMORY(/HIGHWATER)/1024; - start_mem
 print,'Grid Dimensions =',Num_u1,Num_u2,Num_u3


;	Save file for running the code after stop
  save,filename=Plot_Dir+'ReRun_'+StrTrim(Fix(time_iter),2)+'.sav',E1_con,E1_cov,E2_con,E2_cov,ct,$
  															B1_con,B1_cov,B2_con,B2_cov,B3_con,B3_cov,Time_Iter,tt,Ct


;  stop


end			; End time iteration Pro

; ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** *****
;
; ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** *****

Pro Output,	Xarr, Yarr, Zarr, Num_u1, Num_u2, Num_u3, B1_con, B2_con, B3_cov, E1_con, E2_con,$
			tt,Time_Iter,Plot_dir,Plot_file,Ct,h_nu,h_ph,h_mu,Rarr,ThArr,PhiArr,$
			DPsi_dr_S,DPsi_dPh_S,DPsi_dTh_S,Psi_S,B1_S,B2_S,B3_S,E_Th_S,E_Ph_S,Coeffs_B3_S,$
			DPsi_dr_N,DPsi_dPh_N,DPsi_dTh_N,Psi_N,B1_N,B2_N,B3_N,E_Th_N,E_Ph_N,Coeffs_B3_N,Colat_N,Colat_S

 Re_M = 6370e3

 b_nu_a = DblArr(Num_u1/2,Num_u2/2,Num_u3/2)
 b_ph_a = DblArr(Num_u1/2,Num_u2/2,Num_u3/2)
 b_mu_a = DblArr(Num_u1/2,Num_u2/2,Num_u3/2+1)
 e_nu_a = DblArr(Num_u1/2,Num_u2/2,Num_u3/2+1)
 e_ph_a = DblArr(Num_u1/2,Num_u2/2,Num_u3/2+1)

 Xarr1  = DblArr(Num_u1/2,Num_u2/2,Num_u3/2)
 Xarr2  = DblArr(Num_u1/2,Num_u2/2,Num_u3/2+1)
 Yarr1  = DblArr(Num_u1/2,Num_u2/2,Num_u3/2)
 Yarr2  = DblArr(Num_u1/2,Num_u2/2,Num_u3/2+1)
 Zarr1  = DblArr(Num_u1/2,Num_u2/2,Num_u3/2)
 Zarr2  = DblArr(Num_u1/2,Num_u2/2,Num_u3/2+1)

 Xarr1_eq = DblArr(Num_u1/2,Num_u2/2+1)
 Yarr1_eq = DblArr(Num_u1/2,Num_u2/2+1)
 Xarr2_eq = DblArr(Num_u1/2,Num_u2/2+1)
 Yarr2_eq = DblArr(Num_u1/2,Num_u2/2+1)


 b_nu_eq  = DblArr(Num_u1/2,Num_u2/2+1)
 b_ph_eq  = DblArr(Num_u1/2,Num_u2/2+1)
 b_mu_eq  = DblArr(Num_u1/2,Num_u2/2+1)
 e_nu_eq  = DblArr(Num_u1/2,Num_u2/2+1)
 e_ph_eq  = DblArr(Num_u1/2,Num_u2/2+1)


For ii = 0,Num_u1/2-1 do $
 begin
For jj = 0,Num_u2/2-1 do $
 begin
For kk = 0,Num_u3/2-1 do $
 begin
     b_nu_a(ii,jj,kk) = 1.0e9*    h_nu(2*ii  ,2*jj  ,2*kk+1)*B1_con(2*ii,  2*jj  ,2*kk+1)/Re_m        ; dB in nT
     b_ph_a(ii,jj,kk) = 1.0e9*    h_ph(2*ii+1,2*jj+1,2*kk+1)*B2_con(2*ii+1,2*jj+1,2*kk+1)/Re_m
     b_mu_a(ii,jj,kk) = 1.0e9*1.0/h_mu(2*ii+1,2*jj  ,2*kk)  *B3_cov(2*ii+1,2*jj  ,2*kk)/Re_m
     e_nu_a(ii,jj,kk) = 1.0e3*    h_nu(2*ii+1,2*jj+1,2*kk)  *E1_con(2*ii+1,2*jj+1,2*kk)             ;dE in mV/m
     e_ph_a(ii,jj,kk) = 1.0e3*    h_ph(2*ii  ,2*jj  ,2*kk)  *E2_con(2*ii,  2*jj  ,2*kk)
end
end
end

For ii = 0,Num_u1/2-1 do $	;	uses all cells in U1 direction
Begin
For jj = 0,Num_u2/2-1 do $	;	uses all cells in U2 direction
 Begin
 	 b_mu_a(ii,jj,Num_u3/2) = 1.0e9*1.0/h_mu(2*ii+1,2*jj  ,Num_u3-1)*B3_cov(2*ii+1,2*jj  ,Num_u3-1)/Re_m
     e_nu_a(ii,jj,Num_u3/2) = 1.0e3*    h_nu(2*ii+1,2*jj+1,Num_u3-1)*E1_con(2*ii+1,2*jj+1,Num_u3-1)
     e_ph_a(ii,jj,Num_u3/2) = 1.0e3*    h_ph(2*ii  ,2*jj  ,Num_u3-1)*E2_con(2*ii  ,2*jj  ,Num_u3-1)
 end
end

For ii = 0,Num_u1/2-1 do $
begin
For jj = 0,Num_u2/2-1 do $
 begin
For kk = 0,Num_u3/2-1 do $
 begin

	Xarr1(ii,jj,kk) = Xarr(2*ii+1,2*jj,2*kk+1)
	Xarr2(ii,jj,kk) = Xarr(2*ii+1,2*jj,2*kk)
	Yarr1(ii,jj,kk) = Yarr(2*ii+1,2*jj,2*kk+1)
	Yarr2(ii,jj,kk) = Yarr(2*ii+1,2*jj,2*kk)
	Zarr1(ii,jj,kk) = Zarr(2*ii+1,2*jj,2*kk+1)
	Zarr2(ii,jj,kk) = Zarr(2*ii+1,2*jj,2*kk)
;    Rarr1(ii,jj)	= Rarr(2*ii,Num_u3/2+1) 				; For equatorial Polar plots
;    Parr1(ii,jj)	= Phiarr(2*jj)
 end
 end
end

For ii = 0,Num_u1/2-1 do $	;	uses all cells in U1 direction
Begin
For jj = 0,Num_u2/2-1 do $	;	uses all cells in U2 direction
 Begin
 Xarr2(ii,jj,Num_u3/2) = Xarr(2*ii,2*jj,Num_u3-1)
 Yarr2(ii,jj,Num_u3/2) = Yarr(2*ii,2*jj,Num_u3-1)
 Zarr2(ii,jj,Num_u3/2) = Zarr(2*ii,2*jj,Num_u3-1)
 end
end

Xarr1_Eq(*,0:Num_u2/2-1) = Xarr1(*,*,Num_u3/4)
Yarr1_Eq(*,0:Num_u2/2-1) = Yarr1(*,*,Num_u3/4)
Xarr2_Eq(*,0:Num_u2/2-1) = Xarr2(*,*,Num_u3/4)
Yarr2_Eq(*,0:Num_u2/2-1) = Yarr2(*,*,Num_u3/4)

Xarr1_Eq(*,Num_u2/2) = Xarr1(*,0,Num_u3/4)
Yarr1_Eq(*,Num_u2/2) = Yarr1(*,0,Num_u3/4)
Xarr2_Eq(*,Num_u2/2) = Xarr2(*,0,Num_u3/4)
Yarr2_Eq(*,Num_u2/2) = Yarr2(*,0,Num_u3/4)

b_nu_eq(*,0:Num_u2/2-1) = b_nu_a(*,*,Num_u3/4)
b_ph_eq(*,0:Num_u2/2-1)	= b_ph_a(*,*,Num_u3/4)
b_mu_eq(*,0:Num_u2/2-1) = b_mu_a(*,*,Num_u3/4)
e_nu_eq(*,0:Num_u2/2-1) = e_nu_a(*,*,Num_u3/4)
e_ph_eq(*,0:Num_u2/2-1) = e_ph_a(*,*,Num_u3/4)

b_nu_eq(*,Num_u2/2)  	= b_nu_a(*,0,Num_u3/4)
b_ph_eq(*,Num_u2/2)  	= b_ph_a(*,0,Num_u3/4)
b_mu_eq(*,Num_u2/2)  	= b_mu_a(*,0,Num_u3/4)
e_nu_eq(*,Num_u2/2)  	= e_nu_a(*,0,Num_u3/4)
e_ph_eq(*,Num_u2/2)  	= e_ph_a(*,0,Num_u3/4)




LoadCT,4,/silent
;Device,decomposed=0
;LoadCT,23;13
TVLCT,r,g,b,/get
r(0)=255
g(0)=255
b(0)=255

r(127)=r(0)
g(127)=r(0)
b(127)=r(0)

r(255)=0
g(255)=0
b(255)=0

TVLCT,r,g,b

 Window,2,title='Plots Tme = '+StrTrim(Time_Iter,2)+' tt ='+StrTrim(tt,2),Xsize=900,Ysize=900
 !P.multi=[0,3,5]
 !P.charsize=1.5

 Nlvls 	= 253
 Pscale = 1.0

 Nlvls	= Fix(2*(Nlvls/2))
 Mid 	= Double((Nlvls)/2)

 Bmax 	= max([max(ABS(b_nu_a)),max(ABS(b_ph_a)),max(ABS(b_mu_a))])
 Emax 	= max([max(ABS(E_nu_a)),max(ABS(E_ph_a))])
 Blvls 	= Pscale*Bmax*(Findgen(Nlvls)-Mid)/Mid
 Elvls 	= Pscale*Emax*(Findgen(Nlvls)-Mid)/Mid

 B1Max  = Max(ABS(b_nu_a))
 B1lvls = Pscale*B1max*(Findgen(Nlvls)-Mid+0.5)/Mid
 B2Max  = Max(ABS(b_ph_a))
 B2lvls = Pscale*B2max*(Findgen(Nlvls)-Mid+0.5)/Mid
 B3Max  = Max(ABS(b_mu_a))
 B3lvls = Pscale*B3max*(Findgen(Nlvls)-Mid+0.5)/Mid

 E1Max  = Max(ABS(e_nu_a))
 E1lvls = Pscale*E1max*(Findgen(Nlvls)-Mid+0.5)/Mid
 E2Max  = Max(ABS(e_ph_a))
 E2lvls = Pscale*E2max*(Findgen(Nlvls)-Mid+0.5)/Mid

 If tt EQ 0 then B1lvls=0.0
 If tt EQ 0 then B2lvls=0.0
 If tt EQ 0 then B3lvls=0.0
 If tt EQ 0 then E1lvls=0.0
 If tt EQ 0 then E2lvls=0.0

 contour,Reform(b_nu_a(*,0,*)),Reform(Xarr1(*,0,*)),Reform(Zarr1(*,0,*)),/fill,levels=B1lvls,/isotropic,$
 title='B!d!7t!3!n T = '+String(Time_Iter,Format='(F8.3)'),Xtitle='X (R!dE!n)',Ytitle='Z (R!dE!n)',$
 Xrange=[-max(Xarr1(*,0,*)),max(Xarr1(*,0,*))],YRange=[-max(Zarr1(*,0,*)),max(Zarr1(*,0,*))]

 contour,Reform(b_nu_a(*,Num_u2/4,*)),Reform(Xarr1(*,Num_u2/4,*)),Reform(Zarr1(*,Num_u2/4,*)),/fill,$
 levels=B1lvls,/isotropic,/overplot

 contour,Reform(b_nu_a(*,Num_u2/8,*)),Reform(Yarr1(*,Num_u2/8,*)),Reform(Zarr1(*,Num_u2/8,*)),/fill,levels=B1lvls,/isotropic,$
 title='B!d!7t!3!n Max = '+String(B1max,Format='(F8.3)'),Xtitle='Y (R!dE!n)',Ytitle='Z (R!dE!n)',$
 Xrange=[-max(Yarr1(*,Num_u2/8,*)),max(Yarr1(*,Num_u2/8,*))],YRange=[-max(Zarr1(*,Num_u2/8,*)),max(Zarr1(*,Num_u2/8,*))]

 contour,Reform(b_nu_a(*,3*Num_u2/8,*)),Reform(Yarr1(*,3*Num_u2/8,*)),Reform(Zarr1(*,3*Num_u2/8,*)),/fill,$
 levels=B1lvls,/isotropic,/overplot

 contour,b_nu_eq,Xarr1_eq,Yarr1_eq,/fill,levels=B1lvls,/isotropic,$
 title='B!d!7t!3!n T = '+String(Time_Iter,Format='(F8.3)'),Xtitle='X (R!dE!n)',Ytitle='Y (R!dE!n)'


; ~~~~~~~~~~~~~~~~~

 contour,Reform(b_ph_a(*,0,*)),Reform(Xarr1(*,0,*)),Reform(Zarr1(*,0,*)),/fill,levels=B2lvls,/isotropic,$
 title='B!d!7u!3!n T = '+String(Time_Iter,Format='(F8.3)'),Xtitle='X (R!dE!n)',Ytitle='Z (R!dE!n)',$
 Xrange=[-max(Xarr1(*,0,*)),max(Xarr1(*,0,*))],YRange=[-max(Zarr1(*,0,*)),max(Zarr1(*,0,*))]

 contour,Reform(b_ph_a(*,Num_u2/4,*)),Reform(Xarr1(*,Num_u2/4,*)),Reform(Zarr1(*,Num_u2/4,*)),/fill,$
 levels=B2lvls,/isotropic,/overplot

 contour,Reform(b_ph_a(*,Num_u2/8,*)),Reform(Yarr1(*,Num_u2/8,*)),Reform(Zarr1(*,Num_u2/8,*)),/fill,levels=B2lvls,/isotropic,$
 title='B!d!7u!3!n Max = '+String(B2max,Format='(F8.3)'),Xtitle='Y (R!dE!n)',Ytitle='Z (R!dE!n)',$
 Xrange=[-max(Yarr1(*,Num_u2/8,*)),max(Yarr1(*,Num_u2/8,*))],YRange=[-max(Zarr1(*,Num_u2/8,*)),max(Zarr1(*,Num_u2/8,*))]

 contour,Reform(b_ph_a(*,3*Num_u2/8,*)),Reform(Yarr1(*,3*Num_u2/8,*)),Reform(Zarr1(*,3*Num_u2/8,*)),/fill,$
 levels=B2lvls,/isotropic,/overplot

 contour,b_ph_eq,Xarr1_eq,Yarr1_eq,/fill,levels=B2lvls,/isotropic,$
 title='B!d!7u!3!n T = '+String(Time_Iter,Format='(F8.3)'),Xtitle='X (R!dE!n)',Ytitle='Y (R!dE!n)'

; ~~~~~~~~~~~~~~~~~~~~~~~`

 contour,Reform(b_mu_a(*,0,*)),Reform(Xarr2(*,0,*)),Reform(Zarr2(*,0,*)),/fill,levels=B3lvls,/isotropic,$
 title='B!d!7l!3!n T = '+String(Time_Iter,Format='(F8.3)'),Xtitle='X (R!dE!n)',Ytitle='Z (R!dE!n)',$
 Xrange=[-max(Xarr2(*,0,*)),max(Xarr2(*,0,*))],YRange=[-max(Zarr2(*,0,*)),max(Zarr2(*,0,*))]

 contour,Reform(b_mu_a(*,Num_u2/4,*)),Reform(Xarr2(*,Num_u2/4,*)),Reform(Zarr2(*,Num_u2/4,*)),/fill,$
 levels=B3lvls,/isotropic,/overplot

 contour,Reform(b_mu_a(*,Num_u2/8,*)),Reform(Yarr2(*,Num_u2/8,*)),Reform(Zarr2(*,Num_u2/8,*)),/fill,levels=B3lvls,/isotropic,$
 title='B!d!7l!3!n Max = '+String(B3max,Format='(F8.3)'),Xtitle='Y (R!dE!n)',Ytitle='Z (R!dE!n)',$
 Xrange=[-max(Yarr2(*,Num_u2/8,*)),max(Yarr2(*,Num_u2/8,*))],YRange=[-max(Zarr2(*,Num_u2/8,*)),max(Zarr2(*,Num_u2/8,*))]

 contour,Reform(b_mu_a(*,3*Num_u2/8,*)),Reform(Yarr2(*,3*Num_u2/8,*)),Reform(Zarr2(*,3*Num_u2/8,*)),/fill,$
 levels=B3lvls,/isotropic,/overplot

 contour,b_mu_eq,Xarr1_eq,Yarr1_eq,/fill,levels=B3lvls,/isotropic,$
 title='B!d!7l!3!n T = '+String(Time_Iter,Format='(F8.3)'),Xtitle='X (R!dE!n)',Ytitle='Y (R!dE!n)'

; ~~~~~~~~~~~~~~~~~

 contour,Reform(e_nu_a(*,0,*)),Reform(Xarr2(*,0,*)),Reform(Zarr2(*,0,*)),/fill,levels=E1lvls,/isotropic,$
 title='E!d!7t!3!n T = '+String(Time_Iter,Format='(F8.3)'),Xtitle='X (R!dE!n)',Ytitle='Z (R!dE!n)',$
 Xrange=[-max(Xarr2(*,0,*)),max(Xarr2(*,0,*))],YRange=[-max(Zarr2(*,0,*)),max(Zarr2(*,0,*))]

 contour,Reform(e_nu_a(*,Num_u2/4,*)),Reform(Xarr2(*,Num_u2/4,*)),Reform(Zarr2(*,Num_u2/4,*)),/fill,$
 levels=E1lvls,/isotropic,/overplot

 contour,Reform(e_nu_a(*,Num_u2/8,*)),Reform(Yarr2(*,Num_u2/8,*)),Reform(Zarr2(*,Num_u2/8,*)),/fill,levels=E1lvls,/isotropic,$
 title='E!d!7t!3!n Max = '+String(E1max,Format='(F8.3)'),Xtitle='Y (R!dE!n)',Ytitle='Z (R!dE!n)',$
 Xrange=[-max(Yarr2(*,Num_u2/8,*)),max(Yarr2(*,Num_u2/8,*))],YRange=[-max(Zarr2(*,Num_u2/8,*)),max(Zarr2(*,Num_u2/8,*))]

 contour,Reform(e_nu_a(*,3*Num_u2/8,*)),Reform(Yarr2(*,3*Num_u2/8,*)),Reform(Zarr2(*,3*Num_u2/8,*)),/fill,$
 levels=E1lvls,/isotropic,/overplot

 contour,e_nu_eq,Xarr1_eq,Yarr1_eq,/fill,levels=E1lvls,/isotropic,$
 title='E!d!7t!3!n T = '+String(Time_Iter,Format='(F8.3)'),Xtitle='X (R!dE!n)',Ytitle='Y (R!dE!n)'


; ~~~~~~~~~~~~~~~~~

 contour,Reform(e_ph_a(*,0,*)),Reform(Xarr2(*,0,*)),Reform(Zarr2(*,0,*)),/fill,levels=E2lvls,/isotropic,$
 title='E!d!7u!3!n T = '+String(Time_Iter,Format='(F8.3)'),Xtitle='X (R!dE!n)',Ytitle='Z (R!dE!n)',$
 Xrange=[-max(Xarr2(*,0,*)),max(Xarr2(*,0,*))],YRange=[-max(Zarr2(*,0,*)),max(Zarr2(*,0,*))]

 contour,Reform(e_ph_a(*,Num_u2/4,*)),Reform(Xarr2(*,Num_u2/4,*)),Reform(Zarr2(*,Num_u2/4,*)),/fill,$
 levels=E2lvls,/isotropic,/overplot

 contour,Reform(e_ph_a(*,Num_u2/8,*)),Reform(Yarr2(*,Num_u2/8,*)),Reform(Zarr2(*,Num_u2/8,*)),/fill,levels=E2lvls,/isotropic,$
 title='E!d!7u!3!n Max = '+String(E2max,Format='(F8.3)'),Xtitle='Y (R!dE!n)',Ytitle='Z (R!dE!n)',$
 Xrange=[-max(Yarr2(*,Num_u2/8,*)),max(Yarr2(*,Num_u2/8,*))],YRange=[-max(Zarr2(*,Num_u2/8,*)),max(Zarr2(*,Num_u2/8,*))]

 contour,Reform(e_ph_a(*,3*Num_u2/8,*)),Reform(Yarr2(*,3*Num_u2/8,*)),Reform(Zarr2(*,3*Num_u2/8,*)),/fill,$
 levels=E2lvls,/isotropic,/overplot

 contour,e_ph_eq,Xarr1_eq,Yarr1_eq,/fill,levels=E2lvls,/isotropic,$
 title='E!d!7u!3!n T = '+String(Time_Iter,Format='(F8.3)'),Xtitle='X (R!dE!n)',Ytitle='Y (R!dE!n)'


 Filename = Plot_dir+'PNG\'+Plot_file+StrTrim(Ct,2)+'.PNG'
 image = TVRD(0,0,!d.X_size,!d.Y_size,true=1)
 Write_PNG,filename,image,r,g,b
 Print,'File written to ', Filename


 Window,3,title='Ionosphere Time = '+StrTrim(Time_Iter,2)+'sec tt ='+StrTrim(tt,2),Xsize=900,Ysize=900
 !P.multi=[0,4,4]
 !P.charsize=1.5


 contour,Transpose(B3_N),PhiArr*180.0/!Dpi,(Colat_N*180.0/!Dpi),/fill,nlevels=255,Xstyle=1,Ystyle=1,$
  Title='B!dr!n North T = '+String(Time_Iter,Format='(F8.3)'),Xtitle='!7u!3 (deg) Max ='+StrTrim(max(ABS(B3_N)),2),Ytitle='!7h!3 (deg)'

 contour,Transpose(Dpsi_dr_N),PhiArr*180.0/!Dpi,(Colat_N*180.0/!Dpi),/fill,nlevels=255,Xstyle=1,Ystyle=1,$
  Title='d!7W!3/dr North',Xtitle='!7u!3 (deg) Max = '+StrTrim(max(ABS(Dpsi_dr_N)),2),Ytitle='!7h!3 (deg)'

 contour,Transpose(B3_S),PhiArr*180.0/!Dpi,(Colat_S*180.0/!Dpi),/fill,nlevels=255,Xstyle=1,Ystyle=1,$
  Title='B!dr!n South',Xtitle='!7u!3 (deg) Max ='+StrTrim(max(ABS(B3_S)),2),Ytitle='!7h!3 (deg)'

 contour,Transpose(Dpsi_dr_S),PhiArr*180.0/!Dpi,(Colat_S*180.0/!Dpi),/fill,nlevels=255,Xstyle=1,Ystyle=1,$
  Title='d!7W!3/dr South',Xtitle='!7u!3 (deg) Max ='+StrTrim(max(ABS(Dpsi_dr_S)),2),Ytitle='!7h!3 (deg)'


 contour,Transpose(B1_N),PhiArr*180.0/!Dpi,(Colat_N*180.0/!Dpi),/fill,nlevels=255,Xstyle=1,Ystyle=1,$
  Title='B!d!7h!3!n North',Xtitle='!7u!3 (deg) Max = '+StrTrim(max(ABS(B1_N)),2),Ytitle='!7h!3 (deg)'

 contour,Transpose(Dpsi_dth_N),PhiArr*180.0/!Dpi,(Colat_N*180.0/!Dpi),/fill,nlevels=255,Xstyle=1,Ystyle=1,$
  Title='d!7W!3/d!7h!3 North',Xtitle='!7u!3 (deg) Max = '+StrTrim(max(ABS(Dpsi_dth_N)),2),Ytitle='!7h!3 (deg)'

 contour,Transpose(B1_S),PhiArr*180.0/!Dpi,(Colat_S*180.0/!Dpi),/fill,nlevels=255,Xstyle=1,Ystyle=1,$
  Title='B!d!7h!3!n South',Xtitle='!7u!3 (deg) Max = '+StrTrim(max(ABS(B1_S)),2),Ytitle='!7h!3 (deg)'

 contour,Transpose(Dpsi_dth_S),PhiArr*180.0/!Dpi,(Colat_S*180.0/!Dpi),/fill,nlevels=255,Xstyle=1,Ystyle=1,$
  Title='d!7W!3/d!7h!3 South',Xtitle='!7u!3 (deg) Max ='+StrTrim(max(ABS(Dpsi_dth_S)),2),Ytitle='!7h!3 (deg)'


 contour,Transpose(B2_N),PhiArr*180.0/!Dpi,(Colat_N*180.0/!Dpi),/fill,nlevels=255,Xstyle=1,Ystyle=1,$
  Title='B!d!7h!3!n North',Xtitle='!7u!3 (deg) Max = '+StrTrim(max(ABS(B2_N)),2),Ytitle='!7h!3 (deg)'

 contour,Transpose(Dpsi_dPh_N),PhiArr*180.0/!Dpi,(Colat_N*180.0/!Dpi),/fill,nlevels=255,Xstyle=1,Ystyle=1,$
  Title='d!7W!3/d!7u!3 North',Xtitle='!7u!3 (deg) Max = '+StrTrim(max(ABS(Dpsi_dPh_N)),2),Ytitle='!7h!3 (deg)'

 contour,Transpose(B2_S),PhiArr*180.0/!Dpi,(Colat_S*180.0/!Dpi),/fill,nlevels=255,Xstyle=1,Ystyle=1,$
  Title='B!d!7h!3!n South',Xtitle='!7u!3 (deg) Max = '+StrTrim(max(ABS(B2_S)),2),Ytitle='!7h!3 (deg)'

 contour,Transpose(Dpsi_dPh_S),PhiArr*180.0/!Dpi,(Colat_S*180.0/!Dpi),/fill,nlevels=255,Xstyle=1,Ystyle=1,$
  Title='d!7W!3/d!7u!3 South',Xtitle='!7u!3 (deg) Max = '+StrTrim(max(ABS(Dpsi_dPh_S)),2),Ytitle='!7h!3 (deg)'


contour,Transpose(Psi_N),PhiArr*180.0/!Dpi,(Colat_N*180.0/!Dpi),/fill,nlevels=255,Xstyle=1,Ystyle=1,$
  Title='!7W!3 North',Xtitle='!7u!3 (deg) Max = '+StrTrim(max(ABS(Psi_N)),2),Ytitle='!7h!3 (deg)'

Plot,Coeffs_B3_N,title='Coeff North'

contour,Transpose(Psi_S),PhiArr*180.0/!Dpi,(Colat_S*180.0/!Dpi),/fill,nlevels=255,Xstyle=1,Ystyle=1,$
  Title='!7W!3 South',Xtitle='!7u!3 (deg) Max = '+StrTrim(max(ABS(Psi_S)),2),Ytitle='!7h!3 (deg)'

Plot,Coeffs_B3_S,title='Coeff South'

 Filename = Plot_dir+'PNG\Iono_'+Plot_file+StrTrim(Ct,2)+'.PNG'
 image = TVRD(0,0,!d.X_size,!d.Y_size,true=1)
 Write_PNG,filename,image,r,g,b
 Print,'File written to ', Filename


 Filename1 = Plot_dir+Plot_file+StrTrim(Ct,2)+'.sav'
 Save,Filename=Filename1,b_nu_a,b_ph_a,b_mu_a,e_nu_a,e_ph_a,Ct,time_iter


;	Saving Iono
 Filename2 = Plot_dir+Plot_file+'Iono_'+StrTrim(Ct,2)+'.sav'
 Save,Filename=Filename2,DPsi_dr_S,DPsi_dPh_S,DPsi_dTh_S,Psi_S,B1_S,B2_S,B3_S,E_Th_S,E_Ph_S,Coeffs_B3_S,$
		                 DPsi_dr_N,DPsi_dPh_N,DPsi_dTh_N,Psi_N,B1_N,B2_N,B3_N,E_Th_N,E_Ph_N,Coeffs_B3_N,Colat_N,Colat_S



 Ct = Ct+1


;  If Ct EQ 60 then stop

 End					; End of output Pro

; ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** *****
;
; ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** *****

Pro Iono_Cond,		Num_u1,Num_u2,Num_u3,Colat_N,Colat_S,Sd_N,Sp_N,Sh_N,INVSS_N,Sd_S,Sp_S,Sh_S,$
					INVSS_S,Fn_Iono_N,Fn_Iono_S,ThArr,PhiArr,Iono_file,Plot_Dir

 Re			= Double(6370.e3)						; One Earth Radius (in m)
 SS			= DcomplexArr(2,2)						; 2D form of conductivity Tensor
 InvSS_N	= DcomplexArr(Num_u1,Num_u2,2,2) 		; Inverse of the Conductivity Tensor for Northern Ionosphere
 InvSS_S	= DcomplexArr(Num_u1,Num_u2,2,2) 		; Inverse of the Conductivity Tensor for Southern Ionosphere
 Dip_N		= DblArr(Num_u1,Num_u2)
 Dip_S		= DblArr(Num_u1,Num_u2)

 Sd_N		= DblArr(Num_u1,Num_u2)
 Sp_N		= DblArr(Num_u1,Num_u2)
 Sh_N		= DblArr(Num_u1,Num_u2)

 Sd_S		= DblArr(Num_u1,Num_u2)
 Sp_S		= DblArr(Num_u1,Num_u2)
 Sh_S		= DblArr(Num_u1,Num_u2)

;
 Sig0		= dcomplex(0.0,0.0)
 Sig1		= dcomplex(0.0,0.0)
 Sig2		= dcomplex(0.0,0.0)
 y1			= 0.0
 y2			= 0.0


 For jj = 0,Num_u2-1 do $
 begin
 ;	Northern Ionosphere
 ; ~~~~~~~~~~~~~~~~~~~~~
; OpenR,u7,Fn_Iono_N,/Get_Lun						; Open Northern Hemisphere Conductivity File  ; need to change when we get proper 2D surface
; ReadF,u7,y1										; Reading Info line of file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 For ii = 0,Num_u1-1 do $
 begin

  Theta	 = Colat_N(ii)								; Colatitude
  I = !Dpi/2.d0-Acos(-2.d0*cos(Theta)/(sqrt(1.d0+3.d0*cos(theta)^2)))		; Dip Angle (as per sciffer and waters 2002)

  Dip_N(ii,jj) = I

;  ReadF,u7,y1,y1,y1,y2,y1,y1,y1
;  ReadF,u7,Sig0
;  ReadF,u7,Sig1
;  ReadF,u7,Sig2

  Sd_N(ii,jj) =  1.0e6*Re^2										; Height Integrated Direct Conductivity in Northern Hemipshere
  Sp_N(ii,jj) =  10.0e0*Re^2									; Height Integrated Pederson Conductivity in Northern Hemipshere
  Sh_N(ii,jj) =  12.5e0*Re^2									; Height Integrated Hall Conductivity in Noprthern Hemipshere

;  Sd_N(ii,jj) =  Double(Sig0)*Re^2								; Height Integrated Direct Conductivity in Northern Hemipshere
;  Sp_N(ii,jj) =  Double(Sig1)*Re^2								; Height Integrated Pederson Conductivity in Northern Hemipshere
;  Sh_N(ii,jj) =  Double(Sig2)*Re^2								; Height Integrated Hall Conductivity in Noprthern Hemipshere
;
  Alpha 	= Acos(-2.d0*cos(Theta)/(sqrt(1.d0+3.d0*cos(Theta)^2)))
  Sz_N 		= 	Sd_N(ii,jj)*cos(Alpha)*cos(Alpha) + Sp_N(ii,jj)*sin(Alpha)*sin(Alpha)

  SS(0,0)	=  Sd_N(ii,jj)*Sp_N(ii,jj)/Sz_N       				   	; Lysak 2004
  SS(1,0)	= -Sd_N(ii,jj)*Sh_N(ii,jj)*cos(Alpha)/Sz_N
  SS(0,1)	=  Sd_N(ii,jj)*Sh_N(ii,jj)*cos(Alpha)/Sz_N
  SS(1,1)	=  Sp_N(ii,jj)+(Sh_N(ii,jj)^2*sin(Alpha)^2)/Sz_N

  Inverse_SS = LA_Invert(SS,/DOUBLE)								; Inverse of Conductivity Tensor used to compute E's in Ionosphere

  InvSS_N(ii,jj,0,0) = Inverse_SS(0,0)
  InvSS_N(ii,jj,1,0) = Inverse_SS(1,0)
  InvSS_N(ii,jj,0,1) = Inverse_SS(0,1)
  InvSS_N(ii,jj,1,1) = Inverse_SS(1,1)

  end
;Free_Lun,u7
  end

 For jj = 0,Num_u2-1 do $
 begin
 ;	Southern Ionosphere
 ; ~~~~~~~~~~~~~~~~~~~~~
; OpenR,u7,Fn_Iono_S,/Get_Lun						; Open Southern Hemisphere Conductivity File  ; need to change when we get proper 2D surface
; ReadF,u7,y1										; Reading Info line of file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 For ii = 0,Num_u1-1 do $
 begin

  Theta	 = Colat_S(ii)								; Colatitude
  I = !Dpi/2.d0-Acos(-2.d0*cos(Theta)/(sqrt(1.d0+3.d0*cos(theta)^2)))		; Dip Angle (as per sciffer and waters 2002)

  Dip_S(ii,jj) = I

;  ReadF,u7,y1,y1,y1,y2,y1,y1,y1
;  ReadF,u7,Sig0
;  ReadF,u7,Sig1
;  ReadF,u7,Sig2
;
  Sd_S(ii,jj) =  1.0e6*Re^2										; Height Integrated Direct Conductivity in Northern Hemipshere
  Sp_S(ii,jj) =  10.0e0*Re^2									; Height Integrated Pederson Conductivity in Northern Hemipshere
  Sh_S(ii,jj) =  12.5e0*Re^2									; Height Integrated Hall Conductivity in Noprthern Hemipshere

;  Sd_S(ii,jj) =  Double(Sig0)*Re^2								; Height Integrated Direct Conductivity in Northern Hemipshere
;  Sp_S(ii,jj) =  Double(Sig1)*Re^2								; Height Integrated Pederson Conductivity in Northern Hemipshere
;  Sh_S(ii,jj) =  Double(Sig2)*Re^2								; Height Integrated Hall Conductivity in Noprthern Hemipshere
;
  Alpha 	= Acos(-2.d0*cos(Theta)/(sqrt(1.d0+3.d0*cos(Theta)^2)))
  Sz_S 		= 	Sd_S(ii,jj)*cos(Alpha)*cos(Alpha) + Sp_S(ii,jj)*sin(Alpha)*sin(Alpha)

  SS(0,0)	=  Sd_S(ii,jj)*Sp_S(ii,jj)/Sz_S       				   	; Lysak 2004
  SS(1,0)	= -Sd_S(ii,jj)*Sh_S(ii,jj)*cos(Alpha)/Sz_S
  SS(0,1)	=  Sd_S(ii,jj)*Sh_S(ii,jj)*cos(Alpha)/Sz_S
  SS(1,1)	=  Sp_S(ii,jj)+(Sh_S(ii,jj)^2*sin(Alpha)^2)/Sz_S

  Inverse_SS = LA_Invert(SS,/DOUBLE)								; Inverse of Conductivity Tensor used to compute E's in Ionosphere

  InvSS_S(ii,jj,0,0) = Inverse_SS(0,0)
  InvSS_S(ii,jj,1,0) = Inverse_SS(1,0)
  InvSS_S(ii,jj,0,1) = Inverse_SS(0,1)
  InvSS_S(ii,jj,1,1) = Inverse_SS(1,1)

  end
;  Free_Lun,u7
  end


End						; End of Iono Pro

; ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** *****
;
; ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** *****


Pro Sph_Harm,	Num_u1,Num_u2,Num_L,Num_M,LMin_u1,LMax_u1,Ri,Ree,SHF,SHF_dr,Alpha,Alpha_dr,Zeros,PhiArr,del_Th

 Num_Bfnc	= 2*Num_L * (Num_M-1)	+ Num_L ;	Number of SH basis functions (+,- for each M Ne 0)
; Num_Bfnc	= Num_L * (Num_M-1)	+ Num_L ;	Number of SH basis functions (each M Ne 0)

 im = Dcomplex(0.0,1.0)

 SHF 		= DcomplexArr(Num_u1*Num_u2,Num_Bfnc)
 SHF_dr 	= DcomplexArr(Num_u1*Num_u2,Num_Bfnc)

 NPts		= Num_u1
 IPts		= Npts-2

 x			= DblArr(Npts)                              ; Array for cos(theta)
 Theta		= DblArr(Npts)
 Array		= DcomplexArr(Ipts,Ipts)
 Lp			= DblArr(Ipts)
 Lm			= DblArr(Ipts)

 Plm		= DcomplexArr(Num_L,Num_M,Npts)
 Plm_S		= DcomplexArr(Num_L,Num_M,Npts)

 RTerm		= DComplexArr(Num_L,Num_M)
 RTerm_dr	= DComplexArr(Num_L,Num_M)

 zeros		= DblArr(Num_L,Num_M)
 ortho_n	= DblArr(Num_L,Num_L,Num_M)

For mm = 0,Num_M-1 do $
begin

M_nm		= Double(mm)
Re_eimphi 	= Dcomplex(cos(M_nm*Phiarr))
Im_eimphi 	= Dcomplex(sin(M_nm*Phiarr))

; Construct basis set of pV/pr

 CLat_Min	= asin(sqrt(RI/LMax_u1))*180.0/!dpi                   ; Min Co_Lat in degrees, specifies outer boundary
 CLat_Max	= asin(sqrt(RI/LMin_u1))*180.0/!dpi                   ; Max Co_Lat in degrees for Nth Hemisphere (Min Latitude), specifies inner boundary
 del_Th	 	= double(CLat_Max-CLat_Min)/double(NPts-1.0)          ; del_Theta in degrees

; Low here is low theta or poleward end
; High here is high theta or equatorward point
; If Thi=0 then B2(Azimuthal) is zero else B1 (theta) is zero
 LowBC_switch   = 0				; '1' -> Thi = 0	'0' -> Dthi/Dtheta = 0, PoleWard condition : If Thi=0 then this also sets B2=0
 HighBC_switch	= 0				; '1' -> Thi = 0	'0' -> Dthi/Dtheta = 0, Equatorward condition has B1=0 in the inner boundary condition

; Get cos(theta) points for Legendres
 For ii	= 0,Npts-1 do $
 Begin
  Theta(ii)	= (CLat_min+(ii*del_th))*!Dpi/180.0									; Northern Hemisphere
  X(ii)	= cos(Theta(ii))
 end
; For ii	= 0,Npts-1 do x(ii)	= cos(Theta(ii))
 del_th = del_th*!Dpi/180.0
 h 		= del_th

; Numerical differentiation approx of Laplace Eqn in Spherical coords
; Interior Points
 For ii	= 1,Npts-4 do $
 Begin
  SinT	= Sin(Theta(ii+1))
  CosT	= Cos(Theta(ii+1))
  Array(ii,ii-1)= ( 1.d0/h^2) + CosT/SinT*(-1.d0/(2.d0*h))
  Array(ii,ii) 	= (-2.d0/h^2) - (m_nm^2)/(SinT^2)
  Array(ii,ii+1)= ( 1.d0/h^2) + CosT/SinT*( 1.d0/(2.d0*h))
 end

; Low Latitude End - Forward differences ( dp/dT = 0  -> P(0) =  4/3 P(1) - 1/3 P(2) )
 SinT= Sin(Theta(1))
 CosT= Cos(Theta(1))
 If LowBC_switch eq 0 then $                         ; For deriv=0 BC
 Begin
  P0 			= ( 1.0/h^2) + CosT/SinT*(-1.0/(2.0*h))
  Array(0,0) 	= (-2.0/h^2) - (m_nm^2)/(SinT^2) + 4.0/3.0*P0
  Array(0,1) 	= ( 1.0/h^2) + CosT/SinT*( 1.0/(2.0*h)) - 1.0/3.0*P0
 endif else $
 Begin												; For Thi=0 BC
  P0 = 0.0
  Array(0,0) 	= (-2.0/h^2)                			- (m_nm^2)/(SinT^2) 	+ 4.0/3.0*P0
  Array(0,1) 	= ( 1.0/h^2) + CosT/SinT*( 1.0/(2.0*h))						- 1.0/3.0*P0
 endelse

; High Latitude End - Backwards differences ( dp/dT = 0  -> P(n) =  4/3 P(n-1) - 1/3 P(n-2) )
 SinT			= Sin(Theta(NPts-2))
 CosT			= Cos(Theta(NPts-2))

 If HighBC_switch eq 0 then $
 Begin
  Pn = ( 1.0/h^2) + CosT/SinT*( 1.0/(2.0*h))													; Derivative = 0
  Array(Npts-3,Npts-4)= ( 1.0/h^2) + CosT/SinT*(-1.0/(2.0*h)) - 1.0/3.0*Pn
  Array(Npts-3,Npts-3) 	= (-2.0/h^2) - (m_nm^2)/(SinT^2) + 4.0/3.0*Pn
 endif else $
 Begin
  Pn = 0.0																	; Thi = 0
  Array(Npts-3,Npts-4) = ( 1.0/h^2) + CosT/SinT*(-1.0/(2.0*h)) - 1.0/3.0*Pn
  Array(Npts-3,Npts-3) = (-2.0/h^2) - (m_nm^2)/(SinT^2) + 4.0/3.0*Pn
 endelse

 Array1 = Transpose(Array)
 Eval = LA_EIGENPROBLEM(Array1, EIGENVECTORS=Evec, /Double, scale_result=Escale)            ; Solve Ax=lmbda x Eqn for lmbda and x, lbmda=-l(l+1)

 For ii	= 0,Npts-3 do $
 Begin
  Lp(ii) 	= (-1+sqrt(1-4.0*(EVAL(ii))))/2.0                  ; Solve quadratic lmbda=-l^2-l
  Lm(ii) 	= (-1-sqrt(1-4.0*(EVAL(ii))))/2.0
 end
 Index		= (Sort(Lp))                                        ; Use the +ve solutions

 If Lp(Index(0)) le 0.0 then Index = SHIFT(Index, -1)		; for m_nm = 0 remove the constant solution from basis vector set
 Zeros(0:Num_L-1,mm) = Lp(Index(0:Num_L-1))                 ; Store sorted solns to eigenvalue quadratic

; Now get the sorted eigenvectors - these form the potential basis function set
 For ll 	= 0,Num_L-1 do $
 Begin
  For ii 	= 1,Npts-2 do $
  Begin
   Plm(ll,mm,ii) = Evec(ii-1,Index(ll))
  end
 end

 If LowBC_switch eq 0 then $
 begin
  For ll = 0,Num_L-1 do $												; Derivative = 0
  begin
   Plm(ll,mm,0) = 4.0/3.0*Plm(ll,mm,1)  - 1.0/3.0*Plm(ll,mm,2)			; Low  Lat Boundary
  end
 endif else $
 begin
  Plm(0:Num_L-1,mm,0) 		= 0.0										; Low  Lat Boundary, Thi=0
 endelse

 If HighBC_switch eq 0 then $
 begin
  For ll = 0,Num_L-1 do $															; Derivative = 0
  begin
   Plm(ll,mm,Npts-1) = 4.0/3.0*Plm(ll,mm,Npts-2) - 1.0/3.0*Plm(ll,mm,Npts-3)	; High Lat Boundary
  end
 endif else $
 begin
  Plm(0:Num_L-1,mm,Npts-1) = 0.0												; High Lat Boundary
 endelse

;	Normalise Basis Functions Numerically (using trapezoidal rule)
 For ll = 0 ,Num_L-1 do $
 begin
  For ll1 = 0 ,Num_L-1 do $
  begin
   End1 =(Plm(ll,mm,0)*sqrt(sin(theta(0))))*(Plm(ll1,mm,0)*sqrt(sin(theta(0))))
   End2 =(Plm(ll,mm,Npts-1)*sqrt(sin(theta(Npts-1))))*(Plm(ll1,mm,Npts-1)*sqrt(sin(theta(Npts-1))))
   res=(2.0*(reform(Plm(ll,mm,*)*sqrt(sin(theta))))##Transpose(Reform((Plm(ll1,mm,*)*sqrt(sin(theta))))) - End1 - End2)*(del_th/2.0)
   ortho_n(ll,ll1,mm) = res
  end
 end


;	<DAVE> for the m > 0 -> P cos(m theta) Normalised in phi direction is = Pi (itegrate cos^2 between 0 and 2!pi) hence the Pi term
 For ll	= 0,Num_L-1 do $
 begin
;  Plm(ll,mm,*) = Plm(ll,mm,*)/sqrt(2.0*!Dpi*(ortho_n(ll,ll,mm)))	; Basis functions are now an orthonormal basis set
  Plm(ll,mm,*) = Plm(ll,mm,*)/sqrt(1.0*!Dpi*(ortho_n(ll,ll,mm)))	; Basis functions are now an orthonormal basis set
;  Plm(ll,mm,*) = Plm(ll,mm,*)/sqrt((ortho_n(ll,ll,mm)))				; Basis functions are now an orthonormal basis set

  Plm(ll,mm,*) = Reverse(Plm(ll,mm,*))								; Need to reverse the order in Theta direction due to grid direction (high Colatitude to Low Colatitude)
 end

 For ll = 0,Num_L-1 do $
 Begin
  l	= Zeros(ll,mm)
  RTerm_dr(ll,mm)	= l*Ri^(l-1.0)*(1.0-(Ree/Ri)^(2.0*l+1.0))					; Radial Derivative
  RTerm(ll,mm)		= Ri^(l) + (l/(l+1.0))*Ree^(2.0*l+1.0)*Ri^(-(l+1.0))			; Radial Term
  Plm_S(ll,mm,*) 	= Plm(ll,mm,*)*RTerm_dr(ll,mm)
 end


 For ll	= 0,Num_L-1 do $
 begin
 For jj	= 0,Num_u2-1 do $
 begin
 For ii	= 0,Num_u1-1 do $
 begin

  If mm Eq 0 then SHF(ii+jj*Num_u1,ll) $
     = Plm(ll,mm,ii); * (cos(M_nm*Phiarr(jj)))

  If mm Ne 0 then SHF(ii+jj*Num_u1,ll+Fix(M_nm-1)*Num_L+Num_L) $
     = Plm(ll,mm,ii) * (cos(M_nm*Phiarr(jj)) ); + im*sin(M_nm*Phiarr(jj)))
;     = Plm(ll,mm,ii) * (cos(M_nm*Phiarr(jj))  + im*sin(M_nm*Phiarr(jj)))

  If mm Ne 0 then SHF(ii+jj*Num_u1,ll+Fix(M_nm-1)*Num_L+ Fix(Num_M-1)*Num_L+Num_L) $
     = Plm(ll,mm,ii) * (sin(M_nm*Phiarr(jj)) ); - im*sin(M_nm*Phiarr(jj)))
;     = Plm(ll,mm,ii) * (cos(M_nm*Phiarr(jj))  - im*sin(M_nm*Phiarr(jj)))

  If mm Eq 0 then SHF_dr(ii+jj*Num_u1,ll) $
     = Plm_S(ll,mm,ii); * (cos(M_nm*Phiarr(jj)))

  If mm Ne 0 then SHF_dr(ii+jj*Num_u1,ll+Fix(M_nm-1)*Num_L+Num_L) $
     = Plm_S(ll,mm,ii) * (cos(M_nm*Phiarr(jj)) ); + im*sin(M_nm*Phiarr(jj)))
;     = Plm_S(ll,mm,ii) * (cos(M_nm*Phiarr(jj))  + im*sin(M_nm*Phiarr(jj)))

  If mm Ne 0 then SHF_dr(ii+jj*Num_u1,ll+Fix(M_nm-1)*Num_L+ Fix(Num_M-1)*Num_L+Num_L) $
     = Plm_S(ll,mm,ii) * (sin(M_nm*Phiarr(jj)) ); - im*sin(M_nm*Phiarr(jj)))
;     = Plm_S(ll,mm,ii) * (cos(M_nm*Phiarr(jj))  - im*sin(M_nm*Phiarr(jj)))

 end							; end Colat loop
 end							; end azi loop
 end							; end of L loop

End				; End of M loop for P calculation


 Alpha			= conj(SHF)   ##transpose(SHF)					;	Metric used in SVD
 Alpha_dr		= conj(SHF_dr)##transpose(SHF_dr)				;	Metric used in SVD for Bradial



End 				; End of Spherical harmonic generation


; ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** *****
;
; ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** *****




