pro make_flash_cube,filename,xrange,yrange,zrange,texp,npts=npts,sample=sample,norcut=norcut,force_homology=force_homology
	compile_opt idl2

	load_flash_var, velx, filename, 'velx', xrange, yrange, zrange, sample=sample, xcoord=xcoord, ycoord=ycoord, zcoord=zcoord
	load_flash_var, vely, filename, 'vely', xrange, yrange, zrange, sample=sample
	load_flash_var, velz, filename, 'velz', xrange, yrange, zrange, sample=sample
	load_flash_var, dens, filename, 'dens', xrange, yrange, zrange, sample=sample

	npts = n_elements(xcoord)

	if (n_elements(oversample) ne 0) then begin
		npts = round(npts/2.0^oversample)
		dens = congrid(dens,npts,npts,npts)
		velx = congrid(velx,npts,npts,npts)
		vely = congrid(vely,npts,npts,npts)
		velz = congrid(velz,npts,npts,npts)
		xcoord = congrid(xcoord,npts,/interp)
		ycoord = congrid(ycoord,npts,/interp)
		zcoord = congrid(zcoord,npts,/interp)
	endif

	dx0 = (xcoord[1] - xcoord[0])

	print, total(0.5*dens*(velx^2.0 + vely^2.0 + velz^2.0))*dx0^3.0
	return

	nr = round(npts/2)
	print, nr, npts

	r = dindgen(nr)
	r[0] = 1.e-3
	phi = dindgen(npts)/(npts+1)*2*!dpi
	theta = dindgen(npts)/(npts+1)*!dpi
	dens3d = dblarr(nr)
	flatvelxyz3d = dblarr(nr)
	mid = npts/2.0 - 0.5

	for t=0,n_elements(theta)-1 do begin
		for p=0,n_elements(phi)-1 do begin
			xpos = r*sin(theta[t])*cos(phi[p]) + mid
			ypos = r*sin(theta[t])*sin(phi[p]) + mid
			zpos = r*cos(theta[t]) + mid 

			intdens = interpolate(dens,xpos,ypos,zpos)
			dens3d = dens3d + intdens
			flatvelxyz3d = flatvelxyz3d + intdens*(interpolate(velx,xpos,ypos,zpos)^2.0 + interpolate(vely,xpos,ypos,zpos)^2.0 + interpolate(velz,xpos,ypos,zpos)^2.0)
		endfor
	endfor

	flatvelxyz3d = sqrt(flatvelxyz3d/dens3d)

	print, flatvelxyz3d
	fit = linfit(r, flatvelxyz3d, measure_errors=(dens3d/max(dens3d))^(-1.0), /double)
	dflatvelxyz3d = dblarr(nr-1)
	for i=0,nr-2 do begin
		dflatvelxyz3d[i] = flatvelxyz3d[i+1] - flatvelxyz3d[i]
	endfor

	if keyword_set(norcut) then begin
		rcut = 0
	endif else begin
		for i=nr-2,1,-1 do begin
			if (2.0*abs(dflatvelxyz3d[i] - dflatvelxyz3d[i-1])/abs(dflatvelxyz3d[i] + dflatvelxyz3d[i-1]) - 1.0 gt 0.5) then begin
				rcut = i
				break
			endif
		endfor
	endelse
	print, 'rcut', r[rcut]

	vouter = fit[1]*nr
	print, 'vouter', vouter
	
	load_flash_var, temp, filename, 'temp', xrange, yrange, zrange, sample=sample
	load_flash_var, pres, filename, 'pres', xrange, yrange, zrange, sample=sample

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

		temp = congrid(temp,npts,npts,npts)
		pres = congrid(pres,npts,npts,npts)

		he4  = congrid(he4,npts,npts,npts)
		c12  = congrid(c12,npts,npts,npts)
		o16  = congrid(o16,npts,npts,npts)
		ne20 = congrid(ne20,npts,npts,npts)
		mg24 = congrid(mg24,npts,npts,npts)
		si28 = congrid(si28,npts,npts,npts)
		s32  = congrid(s32,npts,npts,npts)
		ar36 = congrid(ar36,npts,npts,npts)
		ca40 = congrid(ca40,npts,npts,npts)
		ti44 = congrid(ti44,npts,npts,npts)
		cr48 = congrid(cr48,npts,npts,npts)
		fe52 = congrid(fe52,npts,npts,npts)
	endif

	dx = vouter*texp/nr
	arad = 7.5657e-15

	if (keyword_set(force_homology)) then begin
		for i=0,npts-1 do begin
		for j=0,npts-1 do begin
		for k=0,npts-1 do begin
			xpos = double(i) - mid
			ypos = double(j) - mid
			zpos = double(k) - mid
			
			rpos = sqrt(xpos^2.0+ypos^2.0+zpos^2.0)
			locv = dx*rpos/texp
			velx[i,j,k] = locv*xpos/rpos
			vely[i,j,k] = locv*ypos/rpos
			velz[i,j,k] = locv*zpos/rpos
		endfor
		endfor
		endfor
	endif

	dens = dens*(dx0/dx)^3.0
	rade = arad*temp^4.0*(dx0/dx)^4.0
	print, 'vel', vouter, max(sqrt(velx^2.0 + vely^2.0 + velz^2.0)), dx*nr/texp

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

	lmass = dens*dx^3.0
	totmass = total(lmass)
	tothe4  = total(lmass*he4)
	totc12  = total(lmass*c12)
	toto16  = total(lmass*o16)
	totne20 = total(lmass*ne20)
	totmg24 = total(lmass*mg24)
	totsi28 = total(lmass*si28)
	tots32  = total(lmass*s32)
	totar36 = total(lmass*ar36)
	totca40 = total(lmass*ca40)
	totti44 = total(lmass*ti44)
	totcr48 = total(lmass*cr48)
	totfe52 = total(lmass*fe52)
	totni56 = total(lmass*ni56)
	totkine = total(lmass*0.5*(velx^2.0 + vely^2.0 + velz^2.0))

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

	openw, 1, 'cube.dat'

	;Header
	printf, 1, 2, npts, dx, texp, 2, 16, 3, format='(2I4, 2E14.5, 3I4)'
	;Atomic numbers
	printf, 1, 2, 6, 8, 10, 12, 14, 16, 18, 20, 22, 23, 24, 25, 26, 27, 28, format='(16I4)'
	printf, 1, 24, 26, 28, format='(3I4)'
	
	;Dummy array for unused isotopes
	dumm = make_array(npts,npts,npts,value=1.d-30)

	for i=0,npts-1 do begin
	for j=0,npts-1 do begin
	for k=0,npts-1 do begin
		printf, 1, dens[i,j,k], pres[i,j,k], velx[i,j,k], vely[i,j,k], velz[i,j,k], rade[i,j,k], $
			he4[i,j,k], c12[i,j,k], o16[i,j,k], ne20[i,j,k], mg24[i,j,k], $
			si28[i,j,k], s32[i,j,k], ar36[i,j,k], ca40[i,j,k], ti44[i,j,k], $
			dumm[i,j,k], cr48[i,j,k], dumm[i,j,k], fe52[i,j,k], dumm[i,j,k], ni56[i,j,k], $
			cr48[i,j,k], fe52[i,j,k], ni56[i,j,k], $
			format='($,24(E14.5,3X))'
		printf, 1, ''
	endfor
	endfor
	endfor
	close, 1
end
