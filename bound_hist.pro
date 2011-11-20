;pruned_orbit.dat needs to be generated from make_track_file in mode 3.
pro bound_hist, basename, nb=nb, resamp=resamp, rhocut=rhocut, min_hist=min_hist, max_hist=max_hist, $
	oversamp=oversamp, ofile=ofile, linear=linear, calc_extras=calc_extras, include_selfbound=include_selfbound, $
	include_unbound=include_unbound, range=range, fprefix=fprefix, showplot=showplot

g = 6.67428d-8
if (~keyword_set(nb)) then nb = 100
if (~keyword_set(resamp)) then resamp = 0
if (~keyword_set(rhocut)) then rhocut = 1.d-10
if (~keyword_set(oversamp)) then oversamp = 0
if (~keyword_set(ofile)) then ofile = 'pruned_orbit.dat'
if (~keyword_set(efile)) then efile = 'extras.dat'
if (keyword_set(include_unbound)) then linear=1
if (~keyword_set(indexlength)) then indexlength = 5
if (~keyword_set(fprefix)) then begin
	fprefix = ''
endif else begin
	fprefix = fprefix + '_'
endelse
if (~keyword_set(range)) then begin
	filename=basename
	rng = [0,0,1]
endif else begin
	if (n_elements(range) eq 2) then begin
		rng = [range[0], range[1], 1]
	endif else begin
		rng = range
	endelse
endelse

nrows = file_lines(ofile)
openr, lun, 'pruned_orbit.dat', /get_lun
line = ""
readf, lun, line
ncols = n_elements(strsplit(line, /regex, /extract))
odata = dblarr(ncols,nrows)
point_lun, lun, 0
readf, lun, odata
free_lun, lun
numformat = '(I' + string(indexlength) + '.' + string(indexlength) + ')'

for f=rng[0],rng[1],rng[2] do begin
	if (keyword_set(range)) then filename = basename + '_' + string(f, format=numformat)
	unk_names = get_var_list(filename)
	dvar = where(strmatch(unk_names, 'dens') eq 1)
	vxvar = where(strmatch(unk_names, 'velx') eq 1)
	vyvar = where(strmatch(unk_names, 'vely') eq 1)
	vzvar = where(strmatch(unk_names, 'velz') eq 1)
	gvar = where(strmatch(unk_names, 'gpot') eq 1)
	jread_amr, filename, VAR_NAME='all', TREE=tree, DATA=unk, PARAMETERS=params

	time = params.time

	nrows = file_lines(efile)
	openr, lun, 'extras.dat', /get_lun
	line = ""
	readf, lun, line
	ncols = n_elements(strsplit(line, /regex, /extract))
	edata = dblarr(ncols,nrows)
	point_lun, lun, 0
	readf, lun, edata
	free_lun, lun
	m = edata[6]

	t = odata[0,*]
	dattime = min(abs(t - time), tindex)

	ptvec = odata[1:6,tindex]
	obvec = odata[7:12,tindex]
	bndvec = odata[13:18,tindex]
	peakvec = odata[25:30,tindex]
	totvec = odata[19:24,tindex]

	ptvec = totvec - obvec + ptvec

	numblocks = params.totblocks

	;First, do the calculation with no sampling to find min and max.
	if (~keyword_set(min_hist) or ~keyword_set(max_hist)) then begin
		mihist = !values.d_infinity
		mahist = -!values.d_infinity
		for n=0L, (numblocks-1) do begin
			if (tree[n].nodetype NE 1) then continue
				
			if (max(unk[dvar,n,*,*,*]) LT rhocut) then continue
		
			nc = 8
			
			xcoords = ((dindgen(nc) + 1./2. - nc/2.)*(tree[n].size[0]/nc)) + tree[n].coord[0]
			ycoords = ((dindgen(nc) + 1./2. - nc/2.)*(tree[n].size[1]/nc)) + tree[n].coord[1]
			zcoords = ((dindgen(nc) + 1./2. - nc/2.)*(tree[n].size[2]/nc)) + tree[n].coord[2]
		
			ptpos = dblarr(3,nc,nc,nc)
			ptx = ptvec[0] - xcoords
			ptpos[0,*,*,*] = cmreplicate(ptx, [nc, nc])
			pty = ptvec[1] - ycoords
			ptpos[1,*,*,*] = transpose(cmreplicate(pty, [nc, nc]), [1, 0, 2])
			ptz = ptvec[2] - zcoords
			ptpos[2,*,*,*] = transpose(cmreplicate(ptz, [nc, nc]), [2, 1, 0])
			ptr = sqrt(total(ptpos^2, 1))

			gptpot = g*m/ptr

			peakpos = dblarr(3,nc,nc,nc)
			peakx = peakvec[0] - xcoords
			peakpos[0,*,*,*] = cmreplicate(peakx, [nc, nc])
			peaky = peakvec[1] - ycoords
			peakpos[1,*,*,*] = transpose(cmreplicate(peaky, [nc, nc]), [1, 0, 2])
			peakz = peakvec[2] - zcoords
			peakpos[2,*,*,*] = transpose(cmreplicate(peakz, [nc, nc]), [2, 1, 0])
			peakr = sqrt(total(peakpos^2, 1))

			if (~keyword_set(include_selfbound)) then begin
				gpot = reform(unk[gvar,n,*,*,*])
		
				avelx = reform(unk[vxvar,n,*,*,*]) - peakvec[3]
				avely = reform(unk[vyvar,n,*,*,*]) - peakvec[4]
				avelz = reform(unk[vzvar,n,*,*,*]) - peakvec[5]

				;ptfrac = peakpos[0,*,*,*]*ptpos[0,*,*,*] + peakpos[1,*,*,*]*ptpos[1,*,*,*] + peakpos[2,*,*,*]*ptpos[2,*,*,*] 
				;ptfrac = transpose(peakpos) # ptpos
				bindbal = gpot/peakr + abs(gptpot/ptr)
				selfbound = 0.5*(avelx^2 + avely^2 + avelz^2) + gpot	
				
				;if (max(selfbound) LT 0) then continue
			endif
		
			dens = reform(unk[dvar,n,*,*,*])
			hvelx = reform(unk[vxvar,n,*,*,*]) - ptvec[3]
			hvely = reform(unk[vyvar,n,*,*,*]) - ptvec[4]
			hvelz = reform(unk[vzvar,n,*,*,*]) - ptvec[5]
			
			temp = -0.5*(hvelx^2 + hvely^2 + hvelz^2) + gptpot
		
			;if (max(temp) LT 0) then continue
		
			if (~keyword_set(include_selfbound)) then begin
				cells = where(selfbound LT 0); AND bindbal LT 0)
				if (cells[0] ne -1) then temp[cells] = 0.0
			endif
			cells = where(dens LT rhocut)
			if (cells[0] ne -1) then temp[cells] = 0.0
			if (~keyword_set(include_unbound)) then begin
				cells = where(temp LT 0)
				if (cells[0] ne -1) then temp[cells] = 0.0
			endif
		
			temp_ind = where(abs(temp) ne 0.0)
			if (temp_ind[0] ne -1) then begin
				min_temp = min(temp[temp_ind])
				if (min_temp lt mihist) then mihist = min_temp
				max_temp = max(temp[temp_ind])
				if (max_temp gt mahist) then mahist = max_temp
			endif
		endfor
	endif else begin

		mihist = min_hist
		mahist = max_hist
	endelse

	if (~keyword_set(linear)) then begin
		mihist = alog10(mihist)
		mahist = alog10(mahist)
	endif ;else begin
	;	mihist = mihist/10.0 ;Add some buffer in case higher sampling yields different min/max
	;	mahist = mahist*10.0
	;endelse

	print, mihist, mahist

	max_ref = max(tree[*].lrefine)

	;Now resample
	bs = (mahist-mihist)/nb
	hist_bins = mihist + (dindgen(nb) + 0.5)*bs
	histo = dblarr(nb)
	mbound = double(0.)
	if (keyword_set(calc_extras)) then begin
		bndx = double(0.)
		bndy = double(0.)
		bndz = double(0.)
	endif

	for n=0L, (numblocks-1) do begin
		if (tree[n].nodetype NE 1) then continue
			
		;if (max(unk[dvar,n,*,*,*]) LT rhocut) then continue

		lresamp = 2^(min([max_ref - tree[n].lrefine, resamp]) + oversamp)
		nc = 8*lresamp
		
		xcoords = ((dindgen(nc) + 1./2. - nc/2.)*(tree[n].size[0]/nc)) + tree[n].coord[0]
		ycoords = ((dindgen(nc) + 1./2. - nc/2.)*(tree[n].size[1]/nc)) + tree[n].coord[1]
		zcoords = ((dindgen(nc) + 1./2. - nc/2.)*(tree[n].size[2]/nc)) + tree[n].coord[2]

		ptpos = dblarr(3,nc,nc,nc)
		ptx = ptvec[0] - xcoords
		ptpos[0,*,*,*] = cmreplicate(ptx, [nc, nc])
		pty = ptvec[1] - ycoords
		ptpos[1,*,*,*] = transpose(cmreplicate(pty, [nc, nc]), [1, 0, 2])
		ptz = ptvec[2] - zcoords
		ptpos[2,*,*,*] = transpose(cmreplicate(ptz, [nc, nc]), [2, 1, 0])
		ptr = sqrt(total(ptpos^2, 1))

		gptpot = g*m/ptr

		peakpos = dblarr(3,nc,nc,nc)
		peakx = peakvec[0] - xcoords
		peakpos[0,*,*,*] = cmreplicate(peakx, [nc, nc])
		peaky = peakvec[1] - ycoords
		peakpos[1,*,*,*] = transpose(cmreplicate(peaky, [nc, nc]), [1, 0, 2])
		peakz = peakvec[2] - zcoords
		peakpos[2,*,*,*] = transpose(cmreplicate(peakz, [nc, nc]), [2, 1, 0])
		peakr = sqrt(total(peakpos^2, 1))

		if (~keyword_set(include_selfbound)) then begin
			gpot = unk[gvar,n,*,*,*]
			gpot = congrid(reform(gpot), nc, nc, nc, /center, /minus_one)

			avelx = unk[vxvar,n,*,*,*] - peakvec[3]
			avely = unk[vyvar,n,*,*,*] - peakvec[4]
			avelz = unk[vzvar,n,*,*,*] - peakvec[5]
			avelx = congrid(reform(avelx), nc, nc, nc, /center, /minus_one)
			avely = congrid(reform(avely), nc, nc, nc, /center, /minus_one)
			avelz = congrid(reform(avelz), nc, nc, nc, /center, /minus_one)

			;ptfrac = peakpos[0,*,*,*]*ptpos[0,*,*,*] + peakpos[1,*,*,*]*ptpos[1,*,*,*] + peakpos[2,*,*,*]*ptpos[2,*,*,*] 
			;ptfrac = transpose(peakpos) # ptpos
			bindbal = gpot/peakr + abs(gptpot/ptr)
			selfbound = 0.5*(avelx^2 + avely^2 + avelz^2) + gpot	
		endif

		dens = congrid(reform(unk[dvar,n,*,*,*]), nc, nc, nc, /center, /minus_one)
		mass = dens*(tree[n].size[0]/nc)^3.

		if (~keyword_set(include_selfbound)) then begin
			bnd_cut = where(dens GE rhocut AND selfbound LT 0); AND bindbal LT 0)
		endif else begin
			bnd_cut = where(dens GE rhocut)
		endelse
		if (bnd_cut[0] ne -1) then begin ;Count the mass bound to the object, might be different from .dat file because of resampling.
			mbound = mbound + total(mass[bnd_cut])
		endif

		if (keyword_set(calc_extras)) then begin
			if (bnd_cut[0] ne -1) then begin ;Count the mass bound to the object, might be different from .dat file because of resampling.
				xc_arr = cmreplicate(xcoords, [nc, nc])
				yc_arr = transpose(cmreplicate(ycoords, [nc, nc]), [1, 0, 2])
				zc_arr = transpose(cmreplicate(zcoords, [nc, nc]), [2, 1, 0])
				bndx = bndx + total(mass[bnd_cut] * xc_arr[bnd_cut])
				bndy = bndy + total(mass[bnd_cut] * yc_arr[bnd_cut])
				bndz = bndz + total(mass[bnd_cut] * zc_arr[bnd_cut])
			endif
		endif
		
		;if (max(selfbound) LT 0) then continue

		hvelx = congrid(reform(unk[vxvar,n,*,*,*]), nc, nc, nc, /center, /minus_one)
		hvely = congrid(reform(unk[vyvar,n,*,*,*]), nc, nc, nc, /center, /minus_one)
		hvelz = congrid(reform(unk[vzvar,n,*,*,*]), nc, nc, nc, /center, /minus_one)
		hvelx = hvelx - ptvec[3]
		hvely = hvely - ptvec[4]	
		hvelz = hvelz - ptvec[5]	
		temp = -0.5*(hvelx^2 + hvely^2 + hvelz^2) + gptpot

		;if (max(temp) LT 0) then continue

		if (~keyword_set(include_selfbound)) then begin
			cells = where(selfbound LT 0); AND bindbal LT 0)
			if (cells[0] ne -1) then temp[cells] = 0.0
		endif
		cells = where(dens LT rhocut)
		if (cells[0] ne -1) then temp[cells] = 0.0
		if (~keyword_set(include_unbound)) then begin
			cells = where(temp LT 0)
			if (cells[0] ne -1) then temp[cells] = 0.0
		endif

		temp_ind = where(abs(temp) ne 0.0)
		if (temp_ind[0] ne -1) then begin
			mass = mass[temp_ind]
			mass = reform(mass, n_elements(mass), /overwrite)
			temp = temp[temp_ind]
			if (~keyword_set(linear)) then temp = alog10(temp)
			temp = reform(temp, n_elements(temp), /overwrite)
			hdata = histogram(temp, min=mihist, max=mahist, nbins=nb, reverse_indices=rev_ind)
			for i = 0, nb-1 do begin
				if (rev_ind[i] ge rev_ind[i+1]) then continue
				ind = rev_ind[rev_ind[i] : rev_ind[i+1]-1]
				histo(i) = histo(i) + total(mass[ind])
			endfor
		endif
	endfor

	if (keyword_set(calc_extras)) then begin
		bndx = bndx / mbound
		bndy = bndy / mbound
		bndz = bndz / mbound
	endif

	openw,lun,fprefix+'dmde_'+filename+'.dat', /get_lun
	printf, lun, nb, format='(I0,X)'
	format_str = '(' + string(nb) + '(E0,X))'
	printf, lun, total(histo), format='(E0,X)'
	printf, lun, mbound, format='(E0,X)'
	if (keyword_set(calc_extras)) then printf, lun, bndx, bndy, bndz, format='(3(E0,X))'
	printf, lun, time, format='(E0,X)'
	printf, lun, -hist_bins, format=format_str 
	printf, lun, histo, format=format_str
	free_lun,lun

	if (keyword_set(showplot) then begin
		set_plot, 'x'
		plot, histo
	endif
endfor

end
