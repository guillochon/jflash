pro make_flash_1dprof,basename,num,var,xrange,yrange,zrange,nbins=nbins,sample=sample,diff=diff,log=log,frac=frac
	load_flash_var, data, basename+'_'+strn(num, length=4, padtype=1, padchar='0'), var, xrange, yrange, zrange, time=t, sample=sample
	if n_elements(diff) ne 0 then begin
		load_flash_var, data2, basename+'_'+strn(diff, length=4, padtype=1, padchar='0'), var, xrange, yrange, zrange, time=t, sample=sample
		data3 = data2 - data
		if keyword_set(frac) then diffdata = data3 / data else diffdata = data3
	endif else diffdata = data
	dims = size(data, /dimensions)
	if n_elements(dims) eq 1 then dims = [dims[0], 1, 1]
	if n_elements(dims) eq 2 then dims = [dims[0], dims[1], 1]
	if n_elements(nbins) eq 0 then nbins = floor(dims[0]/2.)
	sqsize = (xrange[1]-xrange[0])/dims[0]
	vol = sqsize^3.
	midx = floor(dims[0]/2.)
	midy = floor(dims[1]/2.)
	midz = floor(dims[2]/2.)
	prof = fltarr(nbins)
	binvol = intarr(nbins)
	ratio = nbins/midx
	mass_interior = fltarr(nbins)
	for i=0-midx,dims[0]-1-midx do begin
		for j=0-midy,dims[1]-1-midy do begin
			for k=0-midz,dims[2]-1-midz do begin
				ind = floor(sqrt((i+0.5)^2. + (j+0.5)^2. + (k+0.5)^2.)*ratio)
				if ind lt nbins then begin
					prof[ind] = prof[ind] + diffdata[i+midx,j+midy,k+midz]
					binvol[ind] = binvol[ind] + 1
					mass_interior[ind:nbins-1] = mass_interior[ind:nbins-1] + data[i+midx,j+midy,k+midz]*sqsize^3.
				endif
			endfor
		endfor
	endfor

	;window, xsize=800, ysize=800
	;set_plot, 'x'
	jps_start, filename="onedprof.ps"
	rmax = (xrange[1]-xrange[0])/2.
	if keyword_set(log) then y = alog10(abs(prof/binvol)) else y = prof/binvol
	;plot, rmax*findgen(nbins)/nbins, y ;radial coordinate
	mass_interior = mass_interior/total(data)/sqsize^3.
	plot, mass_interior, y, xtitle="M(<r)/M", ytitle="Log !MD!Mr!3"  ;mass coordinate
	jps_end

end
