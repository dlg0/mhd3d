pro sh_debug

	fileName	= '/home/superdarn/mhd3d/output/sh_debug.nc'
	cdfId	= ncdf_open ( fileName, /noWrite )
	glob	= ncdf_inquire ( cdfId )
	ncdf_varGet, cdfId, 'Y', Y
	ncdf_varGet, cdfId, 'theta0', theta0
	ncdf_varGet, cdfId, 'phi', phi
	ncdf_varGet, cdfId, 'hR', hR
	ncdf_varGet, cdfId, 'hTh', hTh 
	ncdf_varGet, cdfId, 'hPh', hPh 
	ncdf_varGet, cdfId, 'hRFit', hRFit
	ncdf_varGet, cdfId, 'hThFit', hThFit 
	ncdf_varGet, cdfId, 'hPhFit', hPhFit
	ncdf_varGet, cdfId, 'hThFitMan', hThFitMan 
	ncdf_varGet, cdfId, 'hPhFitMan', hPhFitMan


	ncdf_close, cdfId
	
	@daves_paths
	device, decomposed = 0
	loadct, 13, file = colortb_path
	!p.multi = [0,4,3]
	contourFac, transpose(Y), transpose(theta0)*!radeg, transpose(phi)*!radeg
	contourFac, transpose(hR), transpose(theta0)*!radeg, transpose(phi)*!radeg
	contourFac, transpose(hTh), transpose(theta0)*!radeg, transpose(phi)*!radeg
	contourFac, transpose(hPh), transpose(theta0)*!radeg, transpose(phi)*!radeg
	contourFac, transpose(hRFit), transpose(theta0)*!radeg, transpose(phi)*!radeg
	contourFac, transpose(hThFit), transpose(theta0)*!radeg, transpose(phi)*!radeg
	contourFac, transpose(hPhFit), transpose(theta0)*!radeg, transpose(phi)*!radeg
	contourFac, transpose(hThFitMan), transpose(theta0)*!radeg, transpose(phi)*!radeg
	contourFac, transpose(hPhFitMan), transpose(theta0)*!radeg, transpose(phi)*!radeg



	!p.multi = 0
	

stop
end
