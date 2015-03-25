pro flash_multipole_expansion, basename, start,finish, xrange, yrange, zrange, stride=stride, indexlength=indexlength, trackfile=trackfile, sample=sample, maxl=maxl, denscut=denscut
	compile_opt idl2
	if n_elements(indexlength) eq 0 then indexlength = 5
	if n_elements(maxl) eq 0 then maxl = 2
	if n_elements(denscut) eq 0 then denscut = 0.0
	if n_elements(stride) eq 0 then begin
		if (start gt finish) then begin
			stride = -1
		endif else begin
			stride = 1
		endelse
	endif

	poly_gam = 2.0

	if n_elements(trackfile) ne 0 then begin
		track = read_ascii(trackfile)
		trackt = track.field1[0,*]
		trackx = track.field1[1,*]
		tracky = track.field1[2,*]
		trackz = track.field1[3,*]
	endif

	numformat = '(I' + string(indexlength) + '.' + string(indexlength) + ')'

	times = dblarr(floor((finish - start) / stride) + 1)
	mode_kine = dblarr(floor((finish - start) / stride) + 1, maxl+1, 2*maxl+1)
	log_mode_kine = dblarr(floor((finish - start) / stride) + 1, maxl+1, 2*maxl+1)
	norm_mode_kine = dblarr(floor((finish - start) / stride) + 1, maxl+1, 2*maxl+1)
	mode_tote = dblarr(floor((finish - start) / stride) + 1, maxl+1, 2*maxl+1)
	log_mode_tote = dblarr(floor((finish - start) / stride) + 1, maxl+1, 2*maxl+1)
	mode_pres = dblarr(floor((finish - start) / stride) + 1, maxl+1, 2*maxl+1)
	mode_eint = dblarr(floor((finish - start) / stride) + 1, maxl+1, 2*maxl+1)
	mode_gpot = dblarr(floor((finish - start) / stride) + 1, maxl+1, 2*maxl+1)
	cnt = 0
	for i=start,finish,stride do begin
		filename = basename + '_' + string(i, format=numformat)
		print, i, format='("Processing ", I5)'
		if file_test(filename) eq 0 then begin
			print, "Can't find " + filename + ", skipping."
			continue
		endif
		if n_elements(trackfile) eq 0 then begin
			xr = xrange
			yr = yrange
			zr = zrange
			tx = total(xrange) / 2.
			ty = total(yrange) / 2.
			tz = total(zrange) / 2.
		endif else begin
			read_parameters, filename, parameters=params
			time = params.time
			tx = interpol(trackx, trackt, time)
			ty = interpol(tracky, trackt, time)
			tz = interpol(trackz, trackt, time)
			xr = xrange + tx
			yr = yrange + ty
			zr = zrange + tz
		endelse
		dens = jloaddata(filename,'dens',xrange=xr,yrange=yr,zrange=zr,sample=sample,time=time,xcoords=xcoords,ycoords=ycoords,zcoords=zcoords)
		pres = jloaddata(filename,'pres',xrange=xr,yrange=yr,zrange=zr,sample=sample,time=time,xcoords=xcoords,ycoords=ycoords,zcoords=zcoords)
		gpot = jloaddata(filename,'gpot',xrange=xr,yrange=yr,zrange=zr,sample=sample,time=time,xcoords=xcoords,ycoords=ycoords,zcoords=zcoords)/2.0
		eint = jloaddata(filename,'eint',xrange=xr,yrange=yr,zrange=zr,sample=sample,time=time,xcoords=xcoords,ycoords=ycoords,zcoords=zcoords)
		velx = jloaddata(filename,'velx',xrange=xr,yrange=yr,zrange=zr,sample=sample,time=time,xcoords=xcoords,ycoords=ycoords,zcoords=zcoords)
		vely = jloaddata(filename,'vely',xrange=xr,yrange=yr,zrange=zr,sample=sample,time=time,xcoords=xcoords,ycoords=ycoords,zcoords=zcoords)
		velz = jloaddata(filename,'velz',xrange=xr,yrange=yr,zrange=zr,sample=sample,time=time,xcoords=xcoords,ycoords=ycoords,zcoords=zcoords)

		vol = (xcoords[2] - xcoords[1])^3

		cut_index = where(dens lt denscut, count)
		if count ne 0 then dens[cut_index] = 0.0

		dims = size(dens, /dimensions)

		xcoords = xcoords - tx
		ycoords = ycoords - ty
		zcoords = zcoords - tz

		xc = cmreplicate(xcoords, [dims[1], dims[2]])
		yc = transpose(cmreplicate(ycoords, [dims[0], dims[2]]), [1, 0, 2])
		zc = transpose(cmreplicate(zcoords, [dims[1], dims[0]]), [2, 1, 0])
		rad = sqrt(xc^2 + yc^2 + zc^2)

		thc = acos(zc/sqrt(xc^2 + yc^2 + zc^2))
		phc = atan(yc,xc)

		mass = total(dens)
		velx = velx - total(dens*velx)/mass
		vely = vely - total(dens*vely)/mass
		velz = velz - total(dens*velz)/mass
		;vrad = (xc*velx + yc*vely + zc*velz)/rad

		kine = 0.5*dens*(velx^2+vely^2+velz^2)
		;kine = 0.5*dens*vrad^2
		tote = (kine + dens*gpot + dens*eint)
		for l = 0, maxl do begin
			for m = -l, l do begin
				print, l, m, format='("Calculating l = ", I3, ", m = ", I3)' 
				if m gt 0 then begin
					harm = (spher_harm(thc, phc, l, m, /double) + (-1)^m*spher_harm(thc, phc, l, -m, /double))/sqrt(2)
				endif
				if m eq 0 then begin
					harm = spher_harm(thc, phc, l, 0, /double)
				endif
				if m lt 0 then begin
					harm = (spher_harm(thc, phc, l, -m, /double) - (-1)^m*spher_harm(thc, phc, l, m, /double))/sqrt(2)/complex(0.0,1.0)
				endif
				mode_kine[cnt,l,m+l] = vol*total(kine*harm)
				mode_pres[cnt,l,m+l] = vol*total(pres*harm)
				mode_eint[cnt,l,m+l] = vol*total(eint*harm)
				mode_gpot[cnt,l,m+l] = vol*total(gpot*harm)
				mode_tote[cnt,l,m+l] = vol*total(tote*harm)
			endfor
		endfor
		times[cnt] = time
		;norm_mode_kine[cnt,*,*] = mode_kine[cnt,*,*]/total(mode_kine[cnt,*,*])
		;for l = 0, maxl do begin
		;	for m = -l, l do begin
		;		log_mode_kine[cnt,l,m+l] = absalog10(mode_kine[cnt,l,m+l])
		;		log_mode_tote[cnt,l,m+l] = alog10(mode_tote[cnt,l,m+l])
		;		norm_mode_kine[cnt,l,m+l] = alog10(norm_mode_kine[cnt,l,m+l])
		;	endfor
		;endfor
		if cnt ne 0 then begin
			;set_plot, 'ps'
			;device, file='mode_kine.ps'
			;min_mode_kine = !VALUES.F_INFINITY
			;max_mode_kine = -!VALUES.F_INFINITY
			;for l = 0, maxl do begin
			;	for m = -l, l do begin
			;		lmax_mode_kine = max(log_mode_kine[0:cnt,l,m+l])
			;		lmin_mode_kine = min(log_mode_kine[0:cnt,l,m+l])
			;		if lmax_mode_kine gt max_mode_kine then max_mode_kine = lmax_mode_kine
			;		if lmin_mode_kine lt min_mode_kine then min_mode_kine = lmin_mode_kine
			;	endfor
			;endfor
			;plot, times[0:cnt], log_mode_kine[0:cnt, 0, 0], xrange=[min(times[0:cnt]), max(times[0:cnt])], yrange=[min_mode_kine, max_mode_kine]
			;for l = maxl, 1, -1 do begin
			;	if l eq 1 then color = fsc_color('red')
			;	if l eq 2 then color = fsc_color('blue')
			;	if l eq 3 then color = fsc_color('green')
			;	if l eq 4 then color = fsc_color('yellow')
			;	for m = -l, l do begin
			;		oplot, times[0:cnt], log_mode_kine[0:cnt, l, m+l], thick=m+l+1, color=color
			;	endfor
			;endfor
			;device, /close

			openw, 1, 'mode_kine.dat'
			str = strjoin(string(times[0:cnt]))
			printf, 1, str, format='(A'+string(strlen(str))+')'
			for c = 0, cnt do begin
				for l = 0, maxl do begin
					arr = mode_kine[c, l, 0:2*l]
					str = strjoin(string(reform(arr, 2*l+1)))
					printf, 1, str, format='(A'+string(strlen(str))+')'
				endfor
			endfor
			close, 1

			;set_plot, 'ps'
			;device, file='norm_mode_kine.ps'
			;min_norm_mode_kine = !VALUES.F_INFINITY
			;max_norm_mode_kine = -!VALUES.F_INFINITY
			;for l = 0, maxl do begin
			;	for m = -l, l do begin
			;		lmax_norm_mode_kine = max(norm_mode_kine[0:cnt,l,m+l])
			;		lmin_norm_mode_kine = min(norm_mode_kine[0:cnt,l,m+l])
			;		if lmax_norm_mode_kine gt max_norm_mode_kine then max_norm_mode_kine = lmax_norm_mode_kine
			;		if lmin_norm_mode_kine lt min_norm_mode_kine then min_norm_mode_kine = lmin_norm_mode_kine
			;	endfor
			;endfor
			;plot, times[0:cnt], norm_mode_kine[0:cnt, 0, 0], xrange=[min(times[0:cnt]), max(times[0:cnt])], yrange=[min_norm_mode_kine, max_norm_mode_kine]
			;for l = maxl, 1, -1 do begin
			;	if l eq 1 then color = fsc_color('red')
			;	if l eq 2 then color = fsc_color('blue')
			;	if l eq 3 then color = fsc_color('green')
			;	if l eq 4 then color = fsc_color('yellow')
			;	for m = -l, l do begin
			;		oplot, times[0:cnt], norm_mode_kine[0:cnt, l, m+l], thick=m+l+1, color=color
			;	endfor
			;endfor
			;device, /close

			openw, 1, 'norm_mode_kine.dat'
			str = strjoin(string(times[0:cnt]))
			printf, 1, str, format='(A'+string(strlen(str))+')'
			for c = 0, cnt do begin
				for l = 0, maxl do begin
					arr = norm_mode_kine[c, l, 0:2*l]
					str = strjoin(string(reform(arr, 2*l+1)))
					printf, 1, str, format='(A'+string(strlen(str))+')'
				endfor
			endfor
			close, 1

			;set_plot, 'ps'
			;device, file='mode_tote.ps'
			;min_mode_tote = !VALUES.F_INFINITY
			;max_mode_tote = -!VALUES.F_INFINITY
			;for l = 0, maxl do begin
			;	for m = -l, l do begin
			;		lmax_mode_tote = max(log_mode_tote[0:cnt,l,m+l])
			;		lmin_mode_tote = min(log_mode_tote[0:cnt,l,m+l])
			;		if lmax_mode_tote gt max_mode_tote then max_mode_tote = lmax_mode_tote
			;		if lmin_mode_tote lt min_mode_tote then min_mode_tote = lmin_mode_tote
			;	endfor
			;endfor
			;plot, times[0:cnt], log_mode_tote[0:cnt, 0, 0], xrange=[min(times[0:cnt]), max(times[0:cnt])], yrange=[min_mode_tote, max_mode_tote]
			;for l = maxl, 1, -1 do begin
			;	if l eq 1 then color = fsc_color('red')
			;	if l eq 2 then color = fsc_color('blue')
			;	if l eq 3 then color = fsc_color('green')
			;	if l eq 4 then color = fsc_color('yellow')
			;	for m = -l, l do begin
			;		oplot, times[0:cnt], log_mode_tote[0:cnt, l, m+l], thick=m+l+1, color=color
			;	endfor
			;endfor
			;device, /close

			openw, 1, 'mode_tote.dat'
			str = strjoin(string(times[0:cnt]))
			printf, 1, str, format='(A'+string(strlen(str))+')'
			for c = 0, cnt do begin
				for l = 0, maxl do begin
					arr = mode_tote[c, l, 0:2*l]
					str = strjoin(string(reform(arr, 2*l+1)))
					printf, 1, str, format='(A'+string(strlen(str))+')'
				endfor
			endfor
			close, 1
			openw, 1, 'mode_pres.dat'
			str = strjoin(string(times[0:cnt]))
			printf, 1, str, format='(A'+string(strlen(str))+')'
			for c = 0, cnt do begin
				for l = 0, maxl do begin
					arr = mode_pres[c, l, 0:2*l]
					str = strjoin(string(reform(arr, 2*l+1)))
					printf, 1, str, format='(A'+string(strlen(str))+')'
				endfor
			endfor
			close, 1
			openw, 1, 'mode_eint.dat'
			str = strjoin(string(times[0:cnt]))
			printf, 1, str, format='(A'+string(strlen(str))+')'
			for c = 0, cnt do begin
				for l = 0, maxl do begin
					arr = mode_eint[c, l, 0:2*l]
					str = strjoin(string(reform(arr, 2*l+1)))
					printf, 1, str, format='(A'+string(strlen(str))+')'
				endfor
			endfor
			close, 1
			openw, 1, 'mode_gpot.dat'
			str = strjoin(string(times[0:cnt]))
			printf, 1, str, format='(A'+string(strlen(str))+')'
			for c = 0, cnt do begin
				for l = 0, maxl do begin
					arr = mode_gpot[c, l, 0:2*l]
					str = strjoin(string(reform(arr, 2*l+1)))
					printf, 1, str, format='(A'+string(strlen(str))+')'
				endfor
			endfor
			close, 1
		endif

		fsc_undefine, filename, dens, velx, vely, velz, xc, yc, zc, thc
		fsc_undefine, gpot, kine, tote, phc

		cnt = cnt + 1
	endfor
end
