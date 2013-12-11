; make_flash_profile.pro
; Assumes xrange,yrange,and zrange are centered about origin.

function intgpot, pos
	common flash_prof, gpot, xsize, ysize, zsize, s, xcoords, ycoords, zcoords, maxgpot
	return, interpolate(gpot, pos[0], pos[1], pos[2])/maxgpot
end

pro make_flash_profile,basename,start,finish,xrange,yrange,zrange,sample=sample,simsize=simsize,stride=stride,ps=ps,$
	fprefix=fprefix,usegpot=usegpot,gpotfile=gpotfile,outfile=outfile,calcgpot=calcgpot,indexlength=indexlength,vars=vars,speclist=speclist

	compile_opt idl2
	common flash_prof, gpot, xsize, ysize, zsize, s, xcoords, ycoords, zcoords, maxgpot

	G = 6.673e-8
	amu = 1.66053886d-24
	kb = 1.3806503d-16
	ab = 7.56576738d-15
	if ~keyword_set(simsize) then begin
		print, 'Need to set simsize keyword!'
		return
	endif
	if ~keyword_set(stride) then stride = 1
	if ~keyword_set(indexlength) then indexlength = 5
	if ~keyword_set(outfile) then outfile = 'profile.dat'

	!p.color = fsc_color('white',/nodisplay)
	rsamp = 2000
    atmpts = 200
	arrlen = 500
	angres = 16.
    nspec = n_elements(speclist)
	muspec = dblarr(nspec)
    for i=0,nspec-1 do begin
        case speclist[i] of
            'h1' : muspec[i] = 1./2.
            'he4': muspec[i] = 4./3.
            'c12': muspec[i] = 12./7.
            'n14': muspec[i] = 14./8.
            'o16': muspec[i] = 16./9.
            else: print, 'Element not found!'
        endcase
    endfor
	for f=start,finish,stride do begin
		numformat = '(I' + string(indexlength) + '.' + string(indexlength) + ')'
		filename = basename + '_' + string(f, format=numformat)

	if keyword_set(usegpot) then begin
		if keyword_set(calcgpot) then begin
			dens = loaddata(filename,'dens',xrange=xrange,yrange=yrange,zrange=zrange,sample=sample,time=time,xcoords=xcoords,ycoords=ycoords,zcoords=zcoords)
			delt = xcoords[1] - xcoords[0]
			com = make_array(3, value=simsize/2.)
			orig = [xcoords[0], ycoords[0], zcoords[0]]
			gpot = G*delt^3.*multipole_potential(dens,com,orig,delt)
		endif else begin
			if ~keyword_set(gpotfile) then gpotfile = filename
			gpot = abs(loaddata(gpotfile,'gpot',xrange=xrange,yrange=yrange,zrange=zrange,sample=sample,time=time,xcoords=xcoords,ycoords=ycoords,zcoords=zcoords))/2.
		endelse
		s = size(gpot, /dimensions)
	endif

	xsize = xcoords[s[0]-1]-xcoords[0]
	ysize = ycoords[s[1]-1]-ycoords[0]
	zsize = zcoords[s[2]-1]-zcoords[0]
	csize = xcoords[1]-xcoords[0]
	;parinfo = replicate({fixed:0, limited:[1,1], limits:[0.,0.], step:0.01}, 3)
	;parinfo[0].limits = [0,s[0]-1]
	;parinfo[1].limits = [0,s[1]-1]
	;parinfo[2].limits = [0,s[2]-1]
	;maxgpot = max(gpot)
	;print, maxgpot
	;parms = tnmin('intgpot',s/2.,parinfo=parinfo,errmsg=errmsg,/maximize,/autoderivative,status=status)
	;print, status
	;print, errmsg
	;print, parms

	rmax = max([xsize,ysize,zsize])/2.
	r_prof = dblarr(arrlen)
	dens_prof = dblarr(arrlen)
	mass_prof = dblarr(arrlen)
	gpot_prof = dblarr(arrlen)
	eint_prof = dblarr(arrlen)
	pres_prof = dblarr(arrlen)
	temp_prof = dblarr(arrlen)
	max_temp_prof = make_array(arrlen, value=0.0, /double)
	vtot_prof = dblarr(arrlen)
	vinf_prof = dblarr(arrlen)
	angmom_prof = dblarr(arrlen,3)
	angvel_prof = dblarr(arrlen)
	radvel_prof = dblarr(arrlen)
    X_prof = dblarr(arrlen,nspec)
	
	rmin = csize/5.
	rstep = (rmax - rmin)/rsamp
	thstep = !pi/angres
	phstep = thstep

	if ~keyword_set(usegpot) then begin
		npts_prof = make_array(arrlen, /long, value=floor(2.*!pi^2./thstep/phstep))
	endif else begin
		npts_prof = lonarr(arrlen)
	endelse

	maxgpot = max(gpot, gpotmax1dpos)
	gpotmaxpos = double(array_indices(gpot, gpotmax1dpos))
	maxx = gpotmaxpos[0]/s[0]*xsize
	maxy = gpotmaxpos[1]/s[1]*ysize
	maxz = gpotmaxpos[2]/s[2]*zsize
	print, maxx, maxy, maxz
	if keyword_set(usegpot) then begin
		ming = !values.f_infinity
		maxg = 0
		gpotmap = dblarr(rsamp,floor(2*!pi/thstep),floor(!pi/phstep))

		for r_i=0,rsamp-1 do begin
			r = r_i*rstep + rmin
			for th_i=0,floor(2*!pi/thstep)-1 do begin
				th = th_i*thstep
				costh = cos(th)
				sinth = sin(th)
				for ph_i=0,floor(!pi/phstep)-1 do begin
					ph = ph_i*phstep
					rsinph = r*sin(ph)
					xc = maxx + costh*rsinph
					yc = maxy + sinth*rsinph
					zc = maxz + r*cos(ph)
					x = (xc/xsize)*s[0]
					y = (yc/ysize)*s[1]
					z = (zc/zsize)*s[2]
					if z lt 0 or z ge s[2] then continue

					locg = interpolate(gpot,x,y,z)
					gpotmap[r_i,th_i,ph_i] = locg
					if locg gt maxg then maxg = locg
					if locg lt ming then ming = locg	
				endfor
			endfor
		endfor

		gstep = double(arrlen)/(maxg - ming)
		print, maxg
		print, ming

		for r_i=0,rsamp-1 do begin
			r = r_i*rstep + rmin
			for th_i=0,floor(2*!pi/thstep)-1 do begin
				th = th_i*thstep
				costh = cos(th)
				sinth = sin(th)
				for ph_i=0,floor(!pi/phstep)-1 do begin
					ph = ph_i*phstep
					rsinph = r*sin(ph)
					xc = maxx + costh*rsinph
					yc = maxy + sinth*rsinph
					zc = maxz + r*cos(ph)
					x = (xc/xsize)*s[0]
					y = (yc/ysize)*s[1]
					z = (zc/zsize)*s[2]
					if z lt 0 or z ge s[2] then continue

					locg = gpotmap[r_i,th_i,ph_i]
					rg_i = min([max([arrlen - ceil((locg-ming)*gstep),0]),arrlen-1])
					npts_prof[rg_i] = npts_prof[rg_i] + 1
					gpot_prof[rg_i] = gpot_prof[rg_i] + locg
				endfor
			endfor
		endfor

		fsc_undefine, gpot
	endif

	dens = loaddata(filename,'dens',xrange=xrange,yrange=yrange,zrange=zrange,sample=sample,time=time)
	velx = loaddata(filename,'velx',xrange=xrange,yrange=yrange,zrange=zrange,sample=sample)
	vely = loaddata(filename,'vely',xrange=xrange,yrange=yrange,zrange=zrange,sample=sample)
	velz = loaddata(filename,'velz',xrange=xrange,yrange=yrange,zrange=zrange,sample=sample)

    slice_dims = size(dens, /dimensions)

	if ~keyword_set(usegpot) then s = size(dens, /dimensions)

	for r_i=0,rsamp-1 do begin
		r = r_i*rstep + rmin
		for th_i=0,floor(2*!pi/thstep)-1 do begin
			th = th_i*thstep
			costh = cos(th)
			sinth = sin(th)
			for ph_i=0,floor(!pi/phstep)-1 do begin
				ph = ph_i*phstep
				rsinph = r*sin(ph)
				xc = maxx + costh*rsinph
				yc = maxy + sinth*rsinph
				zc = maxz + r*cos(ph)
				x = (xc/xsize)*s[0]
				y = (yc/ysize)*s[1]
				z = (zc/zsize)*s[2]
				if z lt 0 or z ge s[2] then continue

				if keyword_set(usegpot) then begin
					locg = gpotmap[r_i,th_i,ph_i]
					rg_i = min([max([arrlen - ceil((locg-ming)*gstep),0]),arrlen-1])
				endif else rg_i = r_i

				rho = interpolate(dens,x,y,z)
				vx = interpolate(velx,x,y,z)
				vy = interpolate(vely,x,y,z)
				vz = interpolate(velz,x,y,z)

				r_prof[rg_i] = r_prof[rg_i] + r
				mass_prof[rg_i] = mass_prof[rg_i] + r^2.*rho
				dens_prof[rg_i] = dens_prof[rg_i] + rho
				vtot_prof[rg_i] = vtot_prof[rg_i] + rho*sqrt(vx^2. + vy^2. + vz^2.)
				vinf_prof[rg_i] = vinf_prof[rg_i] + rho*(vx*xc + vy*yc + vz*zc)/r

				radvel_prof[rg_i] = rho*(xc*vx + yc*vy + zc*vz)/sqrt(xc^2. + yc^2. + zc^2.)
				angmom_prof[rg_i,0] = angmom_prof[rg_i,0] + rho*(yc*vz - zc*vy)
				angmom_prof[rg_i,1] = angmom_prof[rg_i,1] + rho*(zc*vx - xc*vz)
				angmom_prof[rg_i,2] = angmom_prof[rg_i,2] + rho*(xc*vy - yc*vx)
			endfor
		endfor
	endfor
	

	fsc_undefine, velx, vely, velz

	eint = loaddata(filename,'eint',xrange=xrange,yrange=yrange,zrange=zrange,sample=sample)
	pres = loaddata(filename,'pres',xrange=xrange,yrange=yrange,zrange=zrange,sample=sample)
	temp = loaddata(filename,'temp',xrange=xrange,yrange=yrange,zrange=zrange,sample=sample)
	
	for r_i=0,rsamp-1 do begin
		r = r_i*rstep + rmin
		for th_i=0,floor(2*!pi/thstep)-1 do begin
			th = th_i*thstep
			costh = cos(th)
			sinth = sin(th)
			for ph_i=0,floor(!pi/phstep)-1 do begin
				ph = ph_i*phstep
				rsinph = r*sin(ph)
				xc = maxx + costh*rsinph
				yc = maxy + sinth*rsinph
				zc = maxz + r*cos(ph)
				x = (xc/xsize)*s[0]
				y = (yc/ysize)*s[1]
				z = (zc/zsize)*s[2]
				if z lt 0 or z ge s[2] then continue

				if keyword_set(usegpot) then begin
					locg = gpotmap[r_i,th_i,ph_i]
					rg_i = min([max([arrlen - ceil((locg-ming)*gstep),0]),arrlen-1])
				endif else rg_i = r_i

				rho = interpolate(dens,x,y,z)
				eint_prof[rg_i] = eint_prof[rg_i] + rho*interpolate(eint,x,y,z)
				pres_prof[rg_i] = pres_prof[rg_i] + rho*interpolate(pres,x,y,z)
				int_temp = interpolate(temp,x,y,z)
				temp_prof[rg_i] = temp_prof[rg_i] + rho*int_temp
				if (int_temp gt max_temp_prof[rg_i]) then max_temp_prof[rg_i] = int_temp
			endfor
		endfor
	endfor

	fsc_undefine, eint, temp
	
    species = dblarr(slice_dims[0], slice_dims[1], slice_dims[2], nspec)
    for i=0,nspec-1 do begin
        species[*,*,*,i] = loaddata(filename, speclist[i], xrange=xrange,yrange=yrange,zrange=zrange,sample=sample)
    endfor

	for sp_i=0, nspec-1 do begin
		curspec = species[*,*,*,sp_i]
		for r_i=0,rsamp-1 do begin
			r = r_i*rstep + rmin
			for th_i=0,floor(2*!pi/thstep)-1 do begin
				th = th_i*thstep
				costh = cos(th)
				sinth = sin(th)
				for ph_i=0,floor(!pi/phstep)-1 do begin
					ph = ph_i*phstep
					rsinph = r*sin(ph)
					xc = maxx + costh*rsinph
					yc = maxy + sinth*rsinph
					zc = maxz + r*cos(ph)
					x = (xc/xsize)*s[0]
					y = (yc/ysize)*s[1]
					z = (zc/zsize)*s[2]
					if z lt 0 or z ge s[2] then continue

					if keyword_set(usegpot) then begin
						locg = gpotmap[r_i,th_i,ph_i]
						rg_i = min([max([arrlen - ceil((locg-ming)*gstep),0]),arrlen-1])
					endif else rg_i = r_i

					rho = interpolate(dens,x,y,z)
					X_prof[rg_i,sp_i] = X_prof[rg_i,sp_i] + rho*interpolate(curspec,x,y,z)
				endfor
			endfor
		endfor
	endfor
	fsc_undefine, species

	ind = where(npts_prof gt 0)
	npts_prof = npts_prof[ind]
	r_prof = r_prof[ind] / npts_prof
	sort_r = sort(r_prof)
	npts_trim = n_elements(r_prof)

	gpot_prof = gpot_prof[ind] / npts_prof
	angmom_prof = [[angmom_prof[ind,0]],[angmom_prof[ind,1]],[angmom_prof[ind,2]]]
	dens_prof = dens_prof[ind]
	eint_prof = eint_prof[ind] / dens_prof
	pres_prof = pres_prof[ind] / dens_prof
	temp_prof = temp_prof[ind] / dens_prof
	tmparr = dblarr(npts_trim,nspec)
	for i=0, nspec-1 do begin
		tmparr[*,i] = X_prof[ind,i] / dens_prof
	endfor
	X_prof = tmparr
	radvel_prof = radvel_prof[ind] / dens_prof
	angvel_prof = sqrt(angmom_prof[*,0]^2. + angmom_prof[*,1]^2. + angmom_prof[*,2]^2.)
	angvel_prof = angvel_prof / dens_prof / r_prof^2.
	vtot_prof = vtot_prof[ind] / dens_prof
	vinf_prof = vinf_prof[ind] / dens_prof
	dens_prof = dens_prof / npts_prof

	r_prof = r_prof[sort_r]
	dens_prof = dens_prof[sort_r]
	gpot_prof = gpot_prof[sort_r]
	eint_prof = eint_prof[sort_r]
	pres_prof = pres_prof[sort_r]
	temp_prof = temp_prof[sort_r]
	angvel_prof = angvel_prof[sort_r]
	radvel_prof = radvel_prof[sort_r]
	vtot_prof = vtot_prof[sort_r]
	vinf_prof = vinf_prof[sort_r]
	tmparr = dblarr(npts_trim,nspec)
	for i=0, nspec-1 do begin
		tmparr[*,i] = X_prof[sort_r,i]
	endfor
	X_prof = tmparr

	dr_prof = dblarr(npts_trim)
	dr_prof[0] = 0.5*(r_prof[0] + r_prof[1])
	for i=1,npts_trim-2 do begin
		dr_prof[i] = 0.5*(r_prof[i+1] - r_prof[i-1])
	endfor
	dr_prof[npts_trim-1] = r_prof[npts_trim-1] - r_prof[npts_trim-2]
	mass_prof = mass_prof[ind]
	mass_prof = 4.*!pi*mass_prof[sort_r]/npts_prof*dr_prof
	mass_int_prof = total(mass_prof, /cumulative)
	kine_prof = dens_prof*vtot_prof^2.

	print, npts_prof
	print, dr_prof
	print, r_prof

    ;Now, tack on an isothermal atmosphere.
    mu_bar = 1./total(X_prof[npts_trim-1]/muspec)
    atm_t = temp_prof[npts_trim-1]
	pres_atm = dblarr(atmpts)
	temp_atm = dblarr(atmpts)
	dens_atm = dblarr(atmpts)
	mass_int_atm = dblarr(atmpts)
	mass_atm = dblarr(atmpts)
	eint_atm = dblarr(atmpts)
	angvel_atm = dblarr(atmpts)
	radvel_atm = dblarr(atmpts)
	r_atm = dblarr(atmpts)

	dr = rmax / 100.
    temp_atm[*] = atm_t
	dens_atm[0] = dens_prof[npts_trim-1]
    pres_atm[0] = dens_atm[0]/mu_bar/amu*kb*atm_t + 1./3.*ab*atm_t^4.
	eint_atm[0] = 1.5/mu_bar/amu*kb*atm_t + ab*atm_t^4./dens_atm[0]
	mass_int_atm[0] = mass_int_prof[npts_trim-1]
	r_atm[0] = r_prof[npts_trim-1]
	;angvel_atm[0] = sqrt(g*mass_int_prof[npts_trim-1]/r_prof[npts_trim-1])/r_prof[npts_trim-1]
	angvel_atm[0] = angvel_prof[npts_trim-1]

    lasti = npts_trim - 1
    for i=1,atmpts-1 do begin
        mass_int_atm[i] = mass_int_atm[i-1] + 4.*!pi*dens_atm[i-1]*r_atm[i-1]*dr
        r_atm[i] = r_atm[i-1] + dr
        pres_atm[i] = pres_atm[i-1] - G*dens_atm[i-1]*dr*mass_int_atm[i]/r_atm[i]^2.
		;angvel_atm[i] = sqrt(g*mass_int_atm[i]/r_atm[i])/r_atm[i]
		angvel_atm[i] = angvel_atm[0] * (r_atm[0]/r_atm[i])^1.5
        if pres_atm[i] lt 0.0 or pres_atm[i] gt 0.99*pres_atm[i-1] then begin
			lasti = i - 1
            break
        endif
		lasti = i
        dens_atm[i] = mu_bar*amu/kb*(pres_atm[i]/temp_atm[i] - ab*temp_atm[i]^3./3.)
		eint_atm[i] = 1.5/mu_bar/amu*kb*temp_atm[i] + ab*temp_atm[i]^4./dens_atm[i]
		mass_atm[i] = dens_atm[i]*dr*4.*!pi*r_atm[i]^2.
    endfor

	print, r_atm
	print, mass_int_atm
	print, pres_atm
	print, dens_atm

	plot, [r_prof, r_atm[1:lasti]], alog10([dens_prof, dens_atm[1:lasti]]), yrange=[-15,15]
	oplot, [r_prof, r_atm[1:lasti]], alog10([angvel_prof, angvel_atm[1:lasti]]), linestyle=1
	oplot, [r_prof, r_atm[1:lasti]], alog10([eint_prof, eint_atm[1:lasti]]), linestyle=2
	oplot, [r_prof, r_atm[1:lasti]], alog10([pres_prof, pres_atm[1:lasti]]), linestyle=3
	;; For KEPLER
    format_str1 = '(' + string(nspec + 7) + 'A-16)'
    format_str2 = '(' + string(nspec + 7) + 'E-16.6)'
	if f eq start then begin
		openw, lun, outfile, /get_lun
		printf, lun, $
			'mass', $
			'radius', $
			'radial vel.', $
            speclist,$
			'density', $
			'temperature', $
			'int. ener.', $
			'ang. vel.', $
			format=format_str1
		free_lun, lun
	endif
	openw, 1, outfile, /append
	for i=0,n_elements(r_prof)-1 do begin
		printf, 1, $
			mass_prof[i], $
			r_prof[i], $
			radvel_prof[i], $
            X_prof[i,0:nspec-1], $
			dens_prof[i], $
			temp_prof[i], $
			eint_prof[i], $
			angvel_prof[i], $
			format=format_str2
	endfor
	;for i=1,lasti do begin
	;	printf, 1, $
	;		mass_atm[i], $
	;		r_atm[i], $
	;		radvel_atm[i], $
    ;        X_prof[npts_trim-1,*], $
	;		dens_atm[i], $
	;		temp_atm[i], $
	;		eint_atm[i], $
	;		angvel_atm[i], $
	;		format=format_str2
	;endfor
	close, 1

	; For Peter's Code
    ;format_str = '(6E-12.4)'
	;endpt = (n_elements(r_prof)-1)
	;lum = 4.0d0*!pi*5.67e-5*temp_prof[endpt]^4.*r_prof[endpt]^2.
	;print, lum, mass_int_prof[endpt]
	;lum_prof = total(mass_int_prof*temp_prof, /cumulative)
	;lum_prof = lum_prof / lum_prof[endpt]*lum/3.85e33
	;openw, lun, outfile, /get_lun
	;for i=0,endpt do begin
	;	printf, lun, $
	;		mass_int_prof[i]/mass_int_prof[endpt], $
	;		r_prof[i]/r_prof[endpt], $
	;		temp_prof[i], $
	;		lum_prof[i], $
	;		dens_prof[i], $
	;		angvel_prof[i], $
	;		format=format_str
	;endfor
	;printf, lun, mass_int_prof[endpt], format='(1E-12.4)'
	;free_lun, lun
	;close, lun
	endfor
end
