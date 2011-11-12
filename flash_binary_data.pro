pro flash_binary_data,filename,var,xrange,yrange,zrange,ures,$
    simsize=simsize,rangemin=rangemin,subdiv=subdiv,$
    rangemax=rangemax,$
    thrvar=thrvar,thrval=thrval,thrtype=thrtype,sample=sample,$
    fprefix=fprefix,extdata=extdata,datatime=datatime,$
    log=log,colorbarcolor=colorbarcolor,charsize=charsize,hidetime=hidetime,$
    hidecolorbar=hidecolorbar,ambval=ambval,exactmult=exactmult,symrange=symrange,lmin=lmin,xticks=xticks,yticks=yticks,$
	output=output,special=special,negative=negative,subtractavg=subtractavg,$
	excision=excision,product=product

	compile_opt idl2

	if n_elements(subdiv) eq 0 then subdiv = 1
	intdata = dblarr(ures, ures, ures)

	for v=0, n_elements(var)-1 do begin
	for i=0, subdiv-1 do begin
		for j=0, subdiv-1 do begin
			for k=0, subdiv-1 do begin
				xr = [xrange[0] + i*(xrange[1] - xrange[0])/subdiv, xrange[0] + (i+1)*(xrange[1] - xrange[0])/subdiv]
				yr = [yrange[0] + j*(yrange[1] - yrange[0])/subdiv, yrange[0] + (j+1)*(yrange[1] - yrange[0])/subdiv]
				zr = [zrange[0] + k*(zrange[1] - zrange[0])/subdiv, zrange[0] + (k+1)*(zrange[1] - zrange[0])/subdiv]
				load_flash_var, data, filename, var[v], xr, yr, zr, simsize=simsize, subtractavg=subtractavg, sample=sample, time=time
				s = size(data, /dimensions)
				imin = i*ures/subdiv
				imax = (i+1)*ures/subdiv-1
				jmin = j*ures/subdiv
				jmax = (j+1)*ures/subdiv-1
				kmin = k*ures/subdiv
				kmax = (k+1)*ures/subdiv-1
				for ii=imin, imax do begin
					for jj=jmin, jmax do begin
						for kk=kmin, kmax do begin
							intdata[ii,jj,kk] = interpolate(data, float(ii)/ures*s[0], float(jj)/ures*s[1], float(kk)/ures*s[2])
						endfor
					endfor
				endfor
			endfor
		endfor
	endfor

	openw, lun, filename + '_' + var[v] + '_' + strtrim(ures, 2) + '.bin', /get_lun
	writeu, lun, intdata
	free_lun, lun
	close, lun
	if var[v] eq 'gpot' then negative = 1 else negative = 0
	make_flash_slice,'test',var[v],3,xrange,yrange,zrange,'z',extdata=reform(intdata[*,*,ures/2], ures, ures), $
		datatime=time, /log, /exactsize, negative=negative
	fsc_undefine, data
	endfor

end
