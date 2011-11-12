pro flash_var_max,basename,start,finish,vars,xrange,yrange,zrange,ps=ps,sample=sample,indexlength=indexlength
	if n_elements(ps) eq 0 then begin
		set_plot, 'x'
		window,1,xsize=1000,ysize=1000
	endif else set_plot, 'ps'
	if n_elements(indexlength) eq 0 then begin
		indexlength = 4
	endif
	y = dblarr(finish-start+1, n_elements(vars))
	x = dblarr(finish-start+1)
	for i=start,finish do begin
		num = strn(i, length=indexlength, padtype=1, padchar='0')
		filename = basename + '_' + num
		var_string = ''
		for j=0,n_elements(vars)-1 do begin
			load_flash_var, slice, filename, vars[j], xrange, yrange, zrange, time=time, $
				dens=dens, temp=temp, velx=velx, vely=vely, velz=velz, gpot=gpot, sample=sample
			x[i-start:finish-start] = time
			y[i-start:finish-start, j] = max(slice)
			print, y
			var_string = var_string + vars[j] + '_'
		endfor
		if keyword_set(ps) then device, filename='max_'+var_string+basename+'.eps', xsize=10.5, ysize=7.5, $
			xoffset=10.75, yoffset=0.25, /inches
		if n_elements(vars) eq 1 then plot, x, y[*,0], xrange=[x[0], max(x)+1e-6], /xstyle $
			else plot, x, (y[*,0] - min(y[*,0]))/(max(y[*,0]) - min(y[*,0])), yrange=[0.0,1.0], /xstyle
		for j=1,n_elements(vars)-1 do oplot, x, (y[*,j] - min(y[*,j]))/(max(y[*,j]) - min(y[*,j])), linestyle=j
		if keyword_set(ps) then device, /close
		fsc_undefine, slice, dens, temp, velx, vely, velz, gpot
	endfor
end
