pro flash_histogram,filename,var,xrange,yrange,zrange,sample=sample,simsize=simsize,stride=stride,ps=ps,$
	cellsize=cellsize,fprefix=fprefix
	
	histmin = 1e5
	histmax = 1e8
	histstep = 0.1
    !p.color = 255
    !p.background = 0
	;load_flash_var, slice, filename, 'dens', xrange, yrange, zrange, time=time, sample=sample, simsize=simsize, dens=dens
	load_flash_var, slice, filename, var, xrange, yrange, zrange, time=time, sample=sample, simsize=simsize, dens=dens
	dims = size(slice, /dim)
	if n_elements(cellsize) eq 1 then begin
		vol = cellsize^3.
	endif else begin
		vol = double((xrange[1]-xrange[0])/dims[0])^3.
	endelse
	totmass = vol*total(dens)
	histdata = dblarr(floor((alog10(histmax) - alog10(histmin))/histstep))
	histlabels = dindgen(floor((alog10(histmax) - alog10(histmin))/histstep))*histstep + alog10(histmin)
	for i=0L,n_elements(slice)-1 do begin
		if slice[i] lt histmin or slice[i] gt histmax then continue
		histind = where(alog10(slice[i]) ge histlabels and alog10(slice[i]) lt histlabels+histstep)
		histdata[histind] = histdata[histind] + slice[i]*vol/totmass
	endfor
	print, histlabels
	print, histdata
end
