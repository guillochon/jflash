pro make_flash_revolve,filename,xrange,yrange,zrange,texp,npts=npts,sample=sample,norcut=norcut,force_homology=force_homology
	compile_opt idl2

	if (keyword_set(sample)) then begin
		if (sample gt 3) then begin
			oversample = sample - 3
			sample = 3
		endif
	endif

	load_flash_var, velx, filename, 'velx', xrange, yrange, zrange, sample=sample, xcoord=xcoord
	load_flash_var, vely, filename, 'vely', xrange, yrange, zrange, sample=sample
	load_flash_var, dens, filename, 'dens', xrange, yrange, zrange, sample=sample

	npts = n_elements(xcoord)

	if (n_elements(oversample) ne 0) then begin
		npts = round(npts/2.0^oversample)
		dens = congrid(dens,npts,newnpts,newnpts)
		velx = congrid(velx,npts,npts,npts)
		vely = congrid(vely,npts,npts,npts)
		xcoord = congrid(xcoord,npts,/interp)
	endif

	dx0 = (xcoord[1] - xcoord[0])

	nr = round(npts/2)
	print, nr, npts

	r = dindgen(nr)
	r[0] = 1.e-3
	theta = dindgen(npts)/(npts+1)*2*!pi
	dens2d = dblarr(nr)
	flatvelxy2d = dblarr(nr)
	mid = npts/2.0 - 0.5

	for t=0,n_elements(theta)-1 do begin
		xpos = r*cos(theta[t]) + mid
		ypos = r*sin(theta[t]) + mid
		zpos = make_array(nr, value=mid)

		intdens = interpolate(dens,xpos,ypos,zpos)
		dens2d = dens2d + intdens
		flatvelxy2d = flatvelxy2d + intdens*(interpolate(velx,xpos,ypos,zpos)^2.0 + interpolate(vely,xpos,ypos,zpos)^2.0)
	endfor

	flatvelxy2d = sqrt(flatvelxy2d/dens2d)

	print, flatvelxy2d
	fit = linfit(r, flatvelxy2d, measure_errors=(dens2d/max(dens2d))^(-1.0), /double)
	dflatvelxy2d = dblarr(nr-1)
	for i=0,nr-2 do begin
		dflatvelxy2d[i] = flatvelxy2d[i+1] - flatvelxy2d[i]
	endfor

	if keyword_set(norcut) then begin
		rcut = 0
	endif else begin
		for i=nr-2,1,-1 do begin
			if (2.0*abs(dflatvelxy2d[i] - dflatvelxy2d[i-1])/abs(dflatvelxy2d[i] + dflatvelxy2d[i-1]) - 1.0 gt 0.5) then begin
				rcut = i
				break
			endif
		endfor
	endelse
	print, 'rcut', r[rcut]

	;vouter = flatvelxy2d[nr-1]
	vouter = fit[1]*nr
	print, 'vouter', vouter
	
	load_flash_var, temp, filename, 'temp', xrange, yrange, zrange, sample=sample
	load_flash_var, pres, filename, 'pres', xrange, yrange, zrange, sample=sample
	load_flash_var, velz, filename, 'velz', xrange, yrange, zrange, sample=sample

	load_flash_var, he4,  filename, 'he4',  xrange, yrange, zrange, sample=sample
	load_flash_var, c12,  filename, 'c12',  xrange, yrange, zrange, sample=sample
	load_flash_var, o16,  filename, 'o16',  xrange, yrange, zrange, sample=sample
	load_flash_var, ne20, filename, 'ne20', xrange, yrange, zrange, sample=sample
	load_flash_var, mg24, filename, 'mg24', xrange, yrange, zrange, sample=sample
	load_flash_var, si28, filename, 'si28', xrange, yrange, zrange, sample=sample
	load_flash_var, s32,  filename, 's32',  xrange, yrange, zrange, sample=sample
	load_flash_var, ar36, filename, 'ar36', xrange, yrange, zrange, sample=sample
	load_flash_var, ca40, filename, 'ca40', xrange, yrange, zrange, sample=sample
	load_flash_var, ti44, filename, 'ti44', xrange, yrange, zrange, sample=sample
	load_flash_var, cr48, filename, 'cr48', xrange, yrange, zrange, sample=sample
	load_flash_var, fe52, filename, 'fe52', xrange, yrange, zrange, sample=sample
	load_flash_var, ni56, filename, 'ni56', xrange, yrange, zrange, sample=sample

	if (n_elements(oversample) ne 0) then begin
		npts = round(npts/2.0^oversample)

		temp = congrid(temp,npts,newnpts,newnpts)
		pres = congrid(pres,npts,newnpts,newnpts)
		velz = congrid(velz,npts,newnpts,newnpts)

		he4  = congrid(he4,npts,newnpts,newnpts)
		c12  = congrid(c12,npts,newnpts,newnpts)
		o16  = congrid(o16,npts,newnpts,newnpts)
		ne20 = congrid(ne20,npts,newnpts,newnpts)
		mg24 = congrid(mg24,npts,newnpts,newnpts)
		si28 = congrid(si28,npts,newnpts,newnpts)
		s32  = congrid(s32,npts,newnpts,newnpts)
		ar36 = congrid(ar36,npts,newnpts,newnpts)
		ca40 = congrid(ca40,npts,newnpts,newnpts)
		ti44 = congrid(ti44,npts,newnpts,newnpts)
		cr48 = congrid(cr48,npts,newnpts,newnpts)
		fe52 = congrid(fe52,npts,newnpts,newnpts)
	endif

	dx = vouter*texp/nr
	arad = 7.5657e-15

	if (n_elements(npts) eq 0) then npts = n_elements(xcoord)

	dens = reform(dens)
	zpos = cmreplicate(dindgen(npts), nr) + 0.5

	dens2d = dblarr(npts,nr)
	rade2d = dens2d
	pres2d = dens2d
	velxy2d = dens2d
	velz2d = dens2d

	he4_2d  = dens2d
	c12_2d  = dens2d
	o16_2d  = dens2d
	ne20_2d = dens2d
	mg24_2d = dens2d
	si28_2d = dens2d
	s32_2d  = dens2d
	ar36_2d = dens2d
	ca40_2d = dens2d
	ti44_2d = dens2d
	cr48_2d = dens2d
	fe52_2d = dens2d
	ni56_2d = dens2d
	dumm_2d = dens2d

	zeroind = where(sqrt(transpose(cmreplicate(r, npts))^2 + (zpos - mid)^2) lt r[rcut], count)

	rpos = transpose(cmreplicate(r,npts))
	for t=0,n_elements(theta)-1 do begin
		xpos = transpose(cmreplicate(r*cos(theta[t]),npts)) + mid
		ypos = transpose(cmreplicate(r*sin(theta[t]),npts)) + mid
		
		intdens = interpolate(dens,xpos,ypos,zpos)
		dens2d = dens2d + intdens
		rade2d = rade2d + intdens*interpolate(temp,xpos,ypos,zpos)
		pres2d = pres2d + intdens*interpolate(pres,xpos,ypos,zpos)
		if (~keyword_set(force_homology)) then begin
			velxy2d = velx2d + intdens*sqrt(interpolate(velx,xpos,ypos,zpos)^2.0 + interpolate(vely,xpos,ypos,zpos)^2.0)
			velz2d = velz2d + intdens*interpolate(velz,xpos,ypos,zpos)
		endif else begin
			locr = dx*sqrt((xpos - mid)^2.0 + (ypos - mid)^2.0 + (zpos - mid)^2.0)
			locv = locr/texp
			velxy2d = velxy2d + intdens*locv*dx*rpos/locr
			velz2d = velz2d + intdens*locv*dx*(zpos - mid)/locr
		endelse

		he4_2d  = he4_2d  + intdens*interpolate(he4 ,xpos,ypos,zpos)
		c12_2d  = c12_2d  + intdens*interpolate(c12 ,xpos,ypos,zpos)
		o16_2d  = o16_2d  + intdens*interpolate(o16 ,xpos,ypos,zpos)
		ne20_2d = ne20_2d + intdens*interpolate(ne20,xpos,ypos,zpos)
		mg24_2d = mg24_2d + intdens*interpolate(mg24,xpos,ypos,zpos)
		si28_2d = si28_2d + intdens*interpolate(si28,xpos,ypos,zpos)
		s32_2d  = s32_2d  + intdens*interpolate(s32 ,xpos,ypos,zpos)
		ar36_2d = ar36_2d + intdens*interpolate(ar36,xpos,ypos,zpos)
		ca40_2d = ca40_2d + intdens*interpolate(ca40,xpos,ypos,zpos)
		ti44_2d = ti44_2d + intdens*interpolate(ti44,xpos,ypos,zpos)
		cr48_2d = cr48_2d + intdens*interpolate(cr48,xpos,ypos,zpos)
		fe52_2d = fe52_2d + intdens*interpolate(fe52,xpos,ypos,zpos)
		ni56_2d = ni56_2d + intdens*interpolate(ni56,xpos,ypos,zpos)
	endfor

	if (count ne 0) then begin
		rade2d[zeroind]  = 0.0
		pres2d[zeroind]  = 0.0
		velxy2d[zeroind]  = 0.0
		velz2d[zeroind]  = 0.0

		he4_2d[zeroind]  = 0.0
		c12_2d[zeroind]  = 0.0
		o16_2d[zeroind]  = 0.0
		ne20_2d[zeroind] = 0.0
		mg24_2d[zeroind] = 0.0
		si28_2d[zeroind] = 0.0
		s32_2d[zeroind]  = 0.0
		ar36_2d[zeroind] = 0.0
		ca40_2d[zeroind] = 0.0
		ti44_2d[zeroind] = 0.0
		cr48_2d[zeroind] = 0.0
		fe52_2d[zeroind] = 0.0
		ni56_2d[zeroind] = 0.0
	endif

	dens2d = transpose(dens2d)
	rade2d = transpose(rade2d)/dens2d
	pres2d = transpose(pres2d)/dens2d
	velxy2d = transpose(velxy2d)/dens2d
	velz2d = transpose(velz2d)/dens2d

	he4_2d  = transpose(he4_2d)/dens2d
	c12_2d  = transpose(c12_2d)/dens2d
	o16_2d  = transpose(o16_2d)/dens2d
	ne20_2d = transpose(ne20_2d)/dens2d
	mg24_2d = transpose(mg24_2d)/dens2d
	si28_2d = transpose(si28_2d)/dens2d
	s32_2d  = transpose(s32_2d)/dens2d
	ar36_2d = transpose(ar36_2d)/dens2d
	ca40_2d = transpose(ca40_2d)/dens2d
	ti44_2d = transpose(ti44_2d)/dens2d
	cr48_2d = transpose(cr48_2d)/dens2d
	fe52_2d = transpose(fe52_2d)/dens2d
	ni56_2d = transpose(ni56_2d)/dens2d
	dumm_2d = transpose(dumm_2d)

	dens2d = dens2d/npts

	if (count ne 0) then begin
		dens2d = transpose(dens2d)
		he4_2d = transpose(he4_2d)
		dens2d[zeroind] = 1.0e-6*min(dens2d)
		he4_2d[zeroind] = 1.0
		dens2d = transpose(dens2d)
		he4_2d = transpose(he4_2d)
	endif

	dens2d = dens2d*(dx0/dx)^3.0
	rade2d = arad*rade2d^4.0*(dx0/dx)^4.0
	print, 'vel', vouter, max(sqrt(velxy2d^2.0 + velz2d^2.0)), dx*nr/texp

	totmass = 0.e0
	tothe4  = 0.e0
	totc12  = 0.e0
	toto16  = 0.e0
	totne20 = 0.e0
	totmg24 = 0.e0
	totsi28 = 0.e0
	tots32  = 0.e0
	totar36 = 0.e0
	totca40 = 0.e0
	totti44 = 0.e0
	totcr48 = 0.e0
	totfe52 = 0.e0
	totni56 = 0.e0
	totkine = 0.e0

	for i=0,nr-1 do begin
		for j=0,npts-1 do begin
			lmass = 2.e0*!pi*(i+0.5)*dx^3*dens2d[i,j]
			totmass = totmass + lmass
			tothe4  = tothe4  + lmass*he4_2d[i,j]
			totc12  = totc12  + lmass*c12_2d[i,j]
			toto16  = toto16  + lmass*o16_2d[i,j]
			totne20 = totne20 + lmass*ne20_2d[i,j]
			totmg24 = totmg24 + lmass*mg24_2d[i,j]
			totsi28 = totsi28 + lmass*si28_2d[i,j]
			tots32  = tots32  + lmass*s32_2d[i,j]
			totar36 = totar36 + lmass*ar36_2d[i,j]
			totca40 = totca40 + lmass*ca40_2d[i,j]
			totti44 = totti44 + lmass*ti44_2d[i,j]
			totcr48 = totcr48 + lmass*cr48_2d[i,j]
			totfe52 = totfe52 + lmass*fe52_2d[i,j]
			totni56 = totni56 + lmass*ni56_2d[i,j]
			totkine = totkine + lmass*0.5*(velxy2d[i,j]^2.0 + velz2d[i,j]^2.0)
		endfor
	endfor

	print, 'totmass', totmass
	print, 'tothe4 ', tothe4 
	print, 'totc12 ', totc12 
	print, 'toto16 ', toto16 
	print, 'totne20', totne20
	print, 'totmg24', totmg24
	print, 'totsi28', totsi28
	print, 'tots32 ', tots32 
	print, 'totar36', totar36
	print, 'totca40', totca40
	print, 'totti44', totti44
	print, 'totcr48', totcr48
	print, 'totfe52', totfe52
	print, 'totni56', totni56
	print, 'totkine', totkine

	openw, 1, '2d_revolve.dat'

	;Header
	printf, 1, 2, nr, dx, texp, 2, 16, 3, format='(2I4, 2E14.5, 3I4)'
	;Atomic numbers
	printf, 1, 2, 6, 8, 10, 12, 14, 16, 18, 20, 22, 23, 24, 25, 26, 27, 28, format='(16I4)'
	printf, 1, 24, 26, 28, format='(3I4)'
	
	for i=0,nr-1 do begin
		for j=0,npts-1 do begin
			printf, 1, dens2d[i,j], pres2d[i,j], velxy2d[i,j], velz2d[i,j], rade2d[i,j], $
				he4_2d[i,j], c12_2d[i,j], o16_2d[i,j], ne20_2d[i,j], mg24_2d[i,j], $
				si28_2d[i,j], s32_2d[i,j], ar36_2d[i,j], ca40_2d[i,j], ti44_2d[i,j], $
				dumm_2d[i,j], cr48_2d[i,j], dumm_2d[i,j], fe52_2d[i,j], dumm_2d[i,j], ni56_2d[i,j], $
				cr48_2d[i,j], fe52_2d[i,j], ni56_2d[i,j], $
				format='($,24(E14.5,3X))'
			printf, 1, ''
		endfor
	endfor
	close, 1
end
