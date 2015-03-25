pro load_flash_var, slice, filename, var, xrange, yrange, $
	zrange, sliceplane=sliceplane, slicetype=slicetype, $
	dens=dens, temp=temp, velx=velx, vely=vely, velz=velz, pres=pres, gpot=gpot,$
	simsize=simsize, time=time, log=log, sample=sample,lwant=lwant, xcoords=xcoords, $
	ycoords=ycoords, zcoords=zcoords, memefficient=memefficient, $
	subtractavg=subtractavg,orbinfo=orbinfo,refcoor=refcoor, special=special, base_state=base_state,$
	trackfile=trackfile, particles=particles,dt=dt, num_particles=num_particles

	compile_opt idl2
	xr = xrange
	yr = yrange
	zr = zrange
	mp = 1.67262158d-24
	kb = 1.3806503d-16
	h = 6.626068d-27
	;FLASH 3 g value
	;g = 6.673d-8
	;FLASH 4 g value
	g = 6.67428d-8

	if n_elements(sliceplane) eq 0 then sliceplane = 'z'
	if n_elements(slicetype) eq 0 then slicetype = 'plane'
	;orbinfo order
	;1: ecc
	;2: mh
	;3: q
	;4: tp

	;First load base variables if necessary
	if var eq 'temp' or var eq 'entropy' or var eq 'machz' or var eq 'ipres' or $
   		var eq 'gpresz' or var eq 'gpresyzmag' or var eq 'pp' or var eq 'cno' or var eq 'he' or var eq 'nuc' or $
		var eq 'gpresz_tidez' or var eq 'imagbeta' then begin
		if n_elements(temp) eq 0 then begin
			temp = (jloaddata(filename,'temp',xrange=xrange,yrange=yrange,zrange=zrange,sample=sample,lwant=lwant,time=time,xcoords=xcoords,ycoords=ycoords,zcoords=zcoords,particles=particles,dt=dt))
		endif
		dims = size(temp, /dimensions)
		if var eq 'temp' then begin
			if keyword_set(memefficient) then slice = temporary(temp) else slice = temp
		endif
	endif
	if total(var eq ['dens','entropy','angmom','ipres','tide','gpresz_tidez','gprez','gpresymag','pp','cno','kine','kineskew','momentum',$
					 'he','nuc','gdensz','tidex','tidey','tidez','h1dens','he4dens','c12dens','n14dens','o16dens','ne20dens','mg24dens',$
					 'si28dens','s32dens','ar36dens','ca40dens','ti44dens','cr48dens','fe52dens','fe54dens','ni56dens','neutdens',$
					 'protdens','gpotener','eintener','torque','torqmom','mach','csnd','acom','acomx','acomy','acomz','mass',$
					 'selfbound','kinpresratio','kinpresdiff','enuctot','a20a32dens','a36a56dens','cfl','mominertia','jeans',$
					 'dens2', 'ecoodens', 'valf', 'imagbeta']) ge 1 then begin
		if n_elements(dens) eq 0 then begin
			dens = jloaddata(filename,'dens',xrange=xrange,yrange=yrange,zrange=zrange,sample=sample,lwant=lwant,time=time,xcoords=xcoords,ycoords=ycoords,zcoords=zcoords,particles=particles,dt=dt)
		endif
		dims = size(dens, /dimensions)
		if var eq 'dens' then begin
			if keyword_set(memefficient) then slice = temporary(dens) else slice = dens
		endif
	endif
	if var eq 'velx' or var eq 'velxy' or var eq 'angvel' or var eq 'angmom' or var eq 'angmomdens' or var eq 'torqmom' or var eq 'cfl' or $
		var eq 'angvelx' or var eq 'angvely' or var eq 'mach' or var eq 'vtot' or var eq 'vtotxy' or var eq 'vtot2' or var eq 'shock' or $
		var eq 'bhbound' or var eq 'kine' or var eq 'kineskew' or var eq 'momentum' or var eq 'selfbound' or var eq 'kinpresratio' or var eq 'kinpresdiff' then begin
		if n_elements(velx) eq 0 then begin
			velx = double(jloaddata(filename,'velx',xrange=xrange,yrange=yrange,zrange=zrange,sample=sample,lwant=lwant,time=time,$
				xcoords=xcoords,ycoords=ycoords,zcoords=zcoords,particles=particles,dt=dt))
		endif
		if var eq 'velx' then slice = velx
		dims = size(velx, /dimensions)
	endif
	if var eq 'vely' or var eq 'velxy' or var eq 'angvel' or var eq 'angmom' or var eq 'angmomdens' or var eq 'torqmom' or var eq 'cfl' or $
		var eq 'angvelx' or var eq 'angvely' or var eq 'mach' or var eq 'vtot' or var eq 'vtotxy' or var eq 'vtot2' or var eq 'shock' or $
		var eq 'bhbound' or var eq 'kine' or var eq 'kineskew' or var eq 'momentum' or var eq 'selfbound' or var eq 'kinpresratio' or var eq 'kinpresdiff' then begin
		if n_elements(vely) eq 0 then begin
			vely = double(jloaddata(filename,'vely',xrange=xrange,yrange=yrange,zrange=zrange,sample=sample,lwant=lwant,time=time,$
				xcoords=xcoords,ycoords=ycoords,zcoords=zcoords,particles=particles,dt=dt))
		endif
		if var eq 'vely' then slice = vely
		dims = size(vely, /dimensions)
	endif
	if var eq 'velz' or var eq 'absvelz' or var eq 'mach' or var eq 'machz' or var eq 'vtot' or var eq 'vtot2' or var eq 'shock' or var eq 'cfl' or $
		var eq 'gvelzz' or var eq 'bhbound' or var eq 'kine' or var eq 'momentum' or var eq 'selfbound' or $
		var eq 'kinpresratio' or var eq 'kinpresdiff' then begin
		if n_elements(velz) eq 0 then begin
			velz = double(jloaddata(filename,'velz',xrange=xrange,yrange=yrange,zrange=zrange,sample=sample,lwant=lwant,time=time,$
				xcoords=xcoords,ycoords=ycoords,zcoords=zcoords,particles=particles,dt=dt))
		endif
		if var eq 'velz' then slice = velz
		if var eq 'absvelz' then slice = abs(velz)
		dims = size(velz, /dimensions)
	endif
	if var eq 'gpot' or var eq 'selfbound' or var eq 'gpotener' then begin
		if n_elements(gpot) eq 0 then begin
			gpot = double(jloaddata(filename,'gpot',xrange=xrange,yrange=yrange,zrange=zrange,sample=sample,lwant=lwant,time=time,xcoords=xcoords,ycoords=ycoords,zcoords=zcoords,particles=particles,dt=dt))
			gpot = gpot/2. ;FLASH doubles this for some reason...
		endif
		if var eq 'gpot' then slice = gpot
		if var eq 'gpot' then begin
			if keyword_set(log) then slice = -slice
		endif
		dims = size(gpot, /dimensions)
	endif
	if total(strcmp(var, ['pres','kinpresratio','kinpresdiff','mach','csnd','cfl','jeans']), /pre) eq 1 then begin
		if n_elements(pres) eq 0 then begin
			pres = (jloaddata(filename,'pres',xrange=xrange,yrange=yrange,zrange=zrange,sample=sample,lwant=lwant,time=time,xcoords=xcoords,ycoords=ycoords,zcoords=zcoords,particles=particles,dt=dt))
		endif
		if var eq 'pres' then slice = pres
		dims = size(pres, /dimensions)
	endif
	if var eq 'enuc' or var eq 'enuctot' or var eq 'enucratio' then begin
		if n_elements(enuc) eq 0 then begin
			enuc = (jloaddata(filename,'enuc',xrange=xrange,yrange=yrange,zrange=zrange,sample=sample,lwant=lwant,time=time,xcoords=xcoords,ycoords=ycoords,zcoords=zcoords,particles=particles,dt=dt))
		endif
		if var eq 'enuc' then slice = enuc
		dims = size(enuc, /dimensions)
	endif
	if var eq 'eint' or var eq 'eintener' or var eq 'enucratio' or var eq 'ecooratio' then begin
		if n_elements(eint) eq 0 then begin
			eint = (jloaddata(filename,'eint',xrange=xrange,yrange=yrange,zrange=zrange,sample=sample,lwant=lwant,time=time,xcoords=xcoords,ycoords=ycoords,zcoords=zcoords,particles=particles,dt=dt))
		endif
		if var eq 'eint' then slice = eint
		dims = size(eint, /dimensions)
	endif
	if var eq 'ecoo' or var eq 'ecooratio' or var eq 'ecoodens' then begin
		if n_elements(ecoo) eq 0 then begin
			ecoo = (jloaddata(filename,'ecoo',xrange=xrange,yrange=yrange,zrange=zrange,sample=sample,lwant=lwant,time=time,xcoords=xcoords,ycoords=ycoords,zcoords=zcoords,particles=particles,dt=dt))
		endif
		if var eq 'ecoo' then slice = ecoo
		dims = size(ecoo, /dimensions)
	endif

	;Special cases
	if var eq 'dens2' then begin
		slice = dens^2.
		dims = size(slice, /dimensions)
	endif
	if var eq 'UT_ratio' then begin
		slice = double(jloaddata(filename,'eint',xrange=xrange,yrange=yrange,zrange=zrange,sample=sample,lwant=lwant,time=time,xcoords=xcoords,particles=particles,dt=dt))
		kine = double(jloaddata(filename,'velx',xrange=xrange,yrange=yrange,zrange=zrange,sample=sample,lwant=lwant,time=time,xcoords=xcoords,particles=particles,dt=dt))^2.
		kine = kine + double(jloaddata(filename,'vely',xrange=xrange,yrange=yrange,zrange=zrange,sample=sample,lwant=lwant,time=time,xcoords=xcoords,particles=particles,dt=dt))^2.
		kine = kine + double(jloaddata(filename,'velz',xrange=xrange,yrange=yrange,zrange=zrange,sample=sample,lwant=lwant,time=time,xcoords=xcoords,particles=particles,dt=dt))^2.
		slice = 2.*slice/kine
		dims = size(slice, /dimensions)
	endif
	if var eq 'enuctot' then begin
		slice = enuc*dens
		dims = size(slice, /dimensions)
	endif
	if var eq 'mass' then begin
		vol = (xcoords[1]-xcoords[0])^3.
		slice = vol*dens
		dims = size(slice, /dimensions)
	endif
	if var eq 'enucratio' then begin
		slice = enuc/eint*dt
		dims = size(slice, /dimensions)
	endif
	if var eq 'ecooenucratio' then begin
		slice = double(jloaddata(filename,'ecoo',xrange=xrange,yrange=yrange,zrange=zrange,sample=sample,lwant=lwant,time=time,particles=particles,dt=dt))
		slice = slice/double(jloaddata(filename,'enuc',xrange=xrange,yrange=yrange,zrange=zrange,sample=sample,lwant=lwant,time=time,particles=particles,dt=dt))
		slice[where(~finite(slice))] = 0.0
		dims = size(slice, /dimensions)
	endif
	if var eq 'ecooratio' then begin
		slice = ecoo/eint
		dims = size(slice, /dimensions)
	endif
	if var eq 'ecoodens' then begin
		slice = ecoo*dens
		dims = size(slice, /dimensions)
	endif
	if var eq 'coonucdiff' then begin
		slice = double(jloaddata(filename,'enuc',xrange=xrange,yrange=yrange,zrange=zrange,sample=sample,lwant=lwant,time=time,particles=particles,dt=dt))
		slice = slice-double(jloaddata(filename,'ecoo',xrange=xrange,yrange=yrange,zrange=zrange,sample=sample,lwant=lwant,time=time,particles=particles,dt=dt))
		dims = size(slice, /dimensions)
	endif
	if var eq 'a20a32dens' then begin
		slice = double(jloaddata(filename,'ne20',xrange=xrange,yrange=yrange,zrange=zrange,sample=sample,lwant=lwant,time=time,particles=particles,dt=dt))
		slice = slice + double(jloaddata(filename,'mg24',xrange=xrange,yrange=yrange,zrange=zrange,sample=sample,lwant=lwant,time=time,particles=particles,dt=dt))
		slice = slice + double(jloaddata(filename,'si28',xrange=xrange,yrange=yrange,zrange=zrange,sample=sample,lwant=lwant,time=time,particles=particles,dt=dt))
		slice = slice + double(jloaddata(filename,'s32',xrange=xrange,yrange=yrange,zrange=zrange,sample=sample,lwant=lwant,time=time,particles=particles,dt=dt))
		slice = slice*dens
		dims = size(slice, /dimensions)
	endif
	if var eq 'a36a56dens' then begin
		slice = double(jloaddata(filename,'ar36',xrange=xrange,yrange=yrange,zrange=zrange,sample=sample,lwant=lwant,time=time,particles=particles,dt=dt))
		slice = slice + double(jloaddata(filename,'ca40',xrange=xrange,yrange=yrange,zrange=zrange,sample=sample,lwant=lwant,time=time,particles=particles,dt=dt))
		slice = slice + double(jloaddata(filename,'ti44',xrange=xrange,yrange=yrange,zrange=zrange,sample=sample,lwant=lwant,time=time,particles=particles,dt=dt))
		slice = slice + double(jloaddata(filename,'cr48',xrange=xrange,yrange=yrange,zrange=zrange,sample=sample,lwant=lwant,time=time,particles=particles,dt=dt))
		slice = slice + double(jloaddata(filename,'fe52',xrange=xrange,yrange=yrange,zrange=zrange,sample=sample,lwant=lwant,time=time,particles=particles,dt=dt))
		slice = slice + double(jloaddata(filename,'ni56',xrange=xrange,yrange=yrange,zrange=zrange,sample=sample,lwant=lwant,time=time,particles=particles,dt=dt))
		slice = slice*dens
		dims = size(slice, /dimensions)
	endif
	if var eq 'composition' then begin
		slice = double(jloaddata(filename,'he4',xrange=xrange,yrange=yrange,zrange=zrange,sample=sample,lwant=lwant,time=time,particles=particles,dt=dt))          * 4.0
		slice = slice + double(jloaddata(filename,'c12',xrange=xrange,yrange=yrange,zrange=zrange,sample=sample,lwant=lwant,time=time,particles=particles,dt=dt))  * 12.0
		slice = slice + double(jloaddata(filename,'o16',xrange=xrange,yrange=yrange,zrange=zrange,sample=sample,lwant=lwant,time=time,particles=particles,dt=dt))  * 16.0
		slice = slice + double(jloaddata(filename,'ne20',xrange=xrange,yrange=yrange,zrange=zrange,sample=sample,lwant=lwant,time=time,particles=particles,dt=dt)) * 20.0
		slice = slice + double(jloaddata(filename,'mg24',xrange=xrange,yrange=yrange,zrange=zrange,sample=sample,lwant=lwant,time=time,particles=particles,dt=dt)) * 24.0
		slice = slice + double(jloaddata(filename,'si28',xrange=xrange,yrange=yrange,zrange=zrange,sample=sample,lwant=lwant,time=time,particles=particles,dt=dt)) * 28.0
		slice = slice + double(jloaddata(filename,'s32',xrange=xrange,yrange=yrange,zrange=zrange,sample=sample,lwant=lwant,time=time,particles=particles,dt=dt))  * 32.0
		slice = slice + double(jloaddata(filename,'ar36',xrange=xrange,yrange=yrange,zrange=zrange,sample=sample,lwant=lwant,time=time,particles=particles,dt=dt)) * 36.0
		slice = slice + double(jloaddata(filename,'ca40',xrange=xrange,yrange=yrange,zrange=zrange,sample=sample,lwant=lwant,time=time,particles=particles,dt=dt)) * 40.0
		slice = slice + double(jloaddata(filename,'ti44',xrange=xrange,yrange=yrange,zrange=zrange,sample=sample,lwant=lwant,time=time,particles=particles,dt=dt)) * 44.0
		slice = slice + double(jloaddata(filename,'cr48',xrange=xrange,yrange=yrange,zrange=zrange,sample=sample,lwant=lwant,time=time,particles=particles,dt=dt)) * 48.0
		slice = slice + double(jloaddata(filename,'fe52',xrange=xrange,yrange=yrange,zrange=zrange,sample=sample,lwant=lwant,time=time,particles=particles,dt=dt)) * 52.0
		slice = slice + double(jloaddata(filename,'ni56',xrange=xrange,yrange=yrange,zrange=zrange,sample=sample,lwant=lwant,time=time,particles=particles,dt=dt)) * 56.0
		dims = size(slice, /dimensions)
	endif
	if var eq 'ash' then begin
		slice = double(jloaddata(filename,'ne20',xrange=xrange,yrange=yrange,zrange=zrange,sample=sample,lwant=lwant,time=time,particles=particles,dt=dt))
		slice = slice + double(jloaddata(filename,'mg24',xrange=xrange,yrange=yrange,zrange=zrange,sample=sample,lwant=lwant,time=time,particles=particles,dt=dt))
		slice = slice + double(jloaddata(filename,'si28',xrange=xrange,yrange=yrange,zrange=zrange,sample=sample,lwant=lwant,time=time,particles=particles,dt=dt))
		slice = slice + double(jloaddata(filename,'s32',xrange=xrange,yrange=yrange,zrange=zrange,sample=sample,lwant=lwant,time=time,particles=particles,dt=dt))
		slice = slice + double(jloaddata(filename,'ar36',xrange=xrange,yrange=yrange,zrange=zrange,sample=sample,lwant=lwant,time=time,particles=particles,dt=dt))
		slice = slice + double(jloaddata(filename,'ca40',xrange=xrange,yrange=yrange,zrange=zrange,sample=sample,lwant=lwant,time=time,particles=particles,dt=dt))
		slice = slice + double(jloaddata(filename,'ti44',xrange=xrange,yrange=yrange,zrange=zrange,sample=sample,lwant=lwant,time=time,particles=particles,dt=dt))
		slice = slice + double(jloaddata(filename,'cr48',xrange=xrange,yrange=yrange,zrange=zrange,sample=sample,lwant=lwant,time=time,particles=particles,dt=dt))
		slice = slice + double(jloaddata(filename,'fe52',xrange=xrange,yrange=yrange,zrange=zrange,sample=sample,lwant=lwant,time=time,particles=particles,dt=dt))
		slice = slice + double(jloaddata(filename,'ni56',xrange=xrange,yrange=yrange,zrange=zrange,sample=sample,lwant=lwant,time=time,particles=particles,dt=dt))
		dims = size(slice, /dimensions)
	endif
	if var eq 'forsdiff' then begin
		fors = double(jloaddata(filename,'fors',xrange=xrange,yrange=yrange,zrange=zrange,sample=sample,lwant=lwant,time=time,particles=particles,dt=dt))
		for2 = double(jloaddata(filename,'for2',xrange=xrange,yrange=yrange,zrange=zrange,sample=sample,lwant=lwant,time=time,particles=particles,dt=dt))
		dens = double(jloaddata(filename,'dens',xrange=xrange,yrange=yrange,zrange=zrange,sample=sample,lwant=lwant,time=time,particles=particles,dt=dt))
		slice = dens*(fors - for2)
		dims = size(slice, /dimensions)
	endif
	if var eq 'c12o16ratio' then begin
		c12 = double(jloaddata(filename,'c12',xrange=xrange,yrange=yrange,zrange=zrange,sample=sample,lwant=lwant,time=time,particles=particles,dt=dt))
		o16 = double(jloaddata(filename,'o16',xrange=xrange,yrange=yrange,zrange=zrange,sample=sample,lwant=lwant,time=time,particles=particles,dt=dt))
		slice = c12/o16
		dims = size(slice, /dimensions)
	endif
	if var eq 'forsdens' then begin
		fors = double(jloaddata(filename,'fors',xrange=xrange,yrange=yrange,zrange=zrange,sample=sample,lwant=lwant,time=time,xcoords=xcoords,particles=particles,dt=dt))
		dens = double(jloaddata(filename,'dens',xrange=xrange,yrange=yrange,zrange=zrange,sample=sample,lwant=lwant,time=time,particles=particles,dt=dt))
		vol = (xcoords[1]-xcoords[0])^3.
		slice = fors*dens*vol
		dims = size(slice, /dimensions)
	endif
	if var eq 'for2dens' then begin
		for2 = double(jloaddata(filename,'for2',xrange=xrange,yrange=yrange,zrange=zrange,sample=sample,lwant=lwant,time=time,xcoords=xcoords,particles=particles,dt=dt))
		dens = double(jloaddata(filename,'dens',xrange=xrange,yrange=yrange,zrange=zrange,sample=sample,lwant=lwant,time=time,particles=particles,dt=dt))
		vol = (xcoords[1]-xcoords[0])^3.
		slice = for2*dens*vol
		dims = size(slice, /dimensions)
	endif
	if var eq 'forsratio' then begin
		fors = double(jloaddata(filename,'fors',xrange=xrange,yrange=yrange,zrange=zrange,sample=sample,lwant=lwant,time=time,particles=particles,dt=dt))
		for2 = double(jloaddata(filename,'for2',xrange=xrange,yrange=yrange,zrange=zrange,sample=sample,lwant=lwant,time=time,particles=particles,dt=dt))
		slice = (fors - for2)/fors
		dims = size(slice, /dimensions)
	endif
	if var eq 'h1prot' then begin
		slice = jloaddata(filename,'h1',xrange=xrange,yrange=yrange,zrange=zrange,sample=sample,lwant=lwant,time=time,particles=particles,dt=dt)
		slice = slice + jloaddata(filename,'prot',xrange=xrange,yrange=yrange,zrange=zrange,sample=sample,lwant=lwant,time=time,particles=particles,dt=dt)
		dims = size(slice, /dimensions)
	endif
	if var eq 'h1dens' then begin
		h1 = jloaddata(filename,'h1',xrange=xrange,yrange=yrange,zrange=zrange,sample=sample,lwant=lwant,time=time,particles=particles,dt=dt)
		slice = dens*h1
	endif
	if var eq 'he4dens' then begin
		he4 = jloaddata(filename,'he4',xrange=xrange,yrange=yrange,zrange=zrange,sample=sample,lwant=lwant,time=time,particles=particles,dt=dt)
		slice = dens*he4
	endif
	if var eq 'c12dens' then begin
		c12 = jloaddata(filename,'c12',xrange=xrange,yrange=yrange,zrange=zrange,sample=sample,lwant=lwant,time=time,particles=particles,dt=dt)
		slice = dens*c12
	endif
	if var eq 'o16dens' then begin
		o16 = jloaddata(filename,'o16',xrange=xrange,yrange=yrange,zrange=zrange,sample=sample,lwant=lwant,time=time,particles=particles,dt=dt)
		slice = dens*o16
	endif
	if var eq 'ne20dens' then begin
		ne20 = jloaddata(filename,'ne20dens',xrange=xrange,yrange=yrange,zrange=zrange,sample=sample,lwant=lwant,time=time,particles=particles,dt=dt)
		slice = dens*ne20
	endif
	if var eq 'mg24dens' then begin
		mg24 = jloaddata(filename,'mg24dens',xrange=xrange,yrange=yrange,zrange=zrange,sample=sample,lwant=lwant,time=time,particles=particles,dt=dt)
		slice = dens*mg24
	endif
	if var eq 'si28dens' then begin
		si28 = jloaddata(filename,'si28',xrange=xrange,yrange=yrange,zrange=zrange,sample=sample,lwant=lwant,time=time,particles=particles,dt=dt)
		slice = dens*si28
	endif
	if var eq 's32dens' then begin
		s32 = jloaddata(filename,'s32',xrange=xrange,yrange=yrange,zrange=zrange,sample=sample,lwant=lwant,time=time,particles=particles,dt=dt)
		slice = dens*s32
	endif
	if var eq 'ar36dens' then begin
		ar36 = jloaddata(filename,'ar36',xrange=xrange,yrange=yrange,zrange=zrange,sample=sample,lwant=lwant,time=time,particles=particles,dt=dt)
		slice = dens*ar36
	endif
	if var eq 'ca40dens' then begin
		ca40 = jloaddata(filename,'ca40',xrange=xrange,yrange=yrange,zrange=zrange,sample=sample,lwant=lwant,time=time,particles=particles,dt=dt)
		slice = dens*ca40
	endif
	if var eq 'ti44dens' then begin
		ti44 = jloaddata(filename,'ti44',xrange=xrange,yrange=yrange,zrange=zrange,sample=sample,lwant=lwant,time=time,particles=particles,dt=dt)
		slice = dens*ti44
	endif
	if var eq 'cr48dens' then begin
		cr48 = jloaddata(filename,'cr48',xrange=xrange,yrange=yrange,zrange=zrange,sample=sample,lwant=lwant,time=time,particles=particles,dt=dt)
		slice = dens*cr48
	endif
	if var eq 'fe52dens' then begin
		fe52 = jloaddata(filename,'fe52',xrange=xrange,yrange=yrange,zrange=zrange,sample=sample,lwant=lwant,time=time,particles=particles,dt=dt)
		slice = dens*fe52
	endif
	if var eq 'fe54dens' then begin
		fe54 = jloaddata(filename,'fe54',xrange=xrange,yrange=yrange,zrange=zrange,sample=sample,lwant=lwant,time=time,particles=particles,dt=dt)
		slice = dens*fe54
	endif
	if var eq 'ni56dens' then begin
		ni56 = jloaddata(filename,'ni56',xrange=xrange,yrange=yrange,zrange=zrange,sample=sample,lwant=lwant,time=time,particles=particles,dt=dt)
		slice = dens*ni56
	endif
	if var eq 'neutdens' then begin
		neut = jloaddata(filename,'neut',xrange=xrange,yrange=yrange,zrange=zrange,sample=sample,lwant=lwant,time=time,particles=particles,dt=dt)
		slice = dens*neut
	endif
	if var eq 'protdens' then begin
		prot = jloaddata(filename,'prot',xrange=xrange,yrange=yrange,zrange=zrange,sample=sample,lwant=lwant,time=time,particles=particles,dt=dt)
		slice = dens*prot
	endif
	if total(var eq ['magx','valf','imagp','imagbeta','divbratio']) ge 1 then begin
		if n_elements(magx) eq 0 then begin
			magx = jloaddata(filename,'magx',xrange=xrange,yrange=yrange,zrange=zrange,sample=sample,lwant=lwant,time=time,xcoords=xcoords,ycoords=ycoords,zcoords=zcoords,particles=particles,dt=dt)
		endif
		dims = size(magx, /dimensions)
		if var eq 'magx' then begin
			if keyword_set(memefficient) then slice = temporary(magx) else slice = magx
		endif
	endif
	if total(var eq ['magy','valf','imagp','imagbeta','divbratio']) ge 1 then begin
		if n_elements(magy) eq 0 then begin
			magy = jloaddata(filename,'magy',xrange=xrange,yrange=yrange,zrange=zrange,sample=sample,lwant=lwant,time=time,xcoords=xcoords,ycoords=ycoords,zcoords=zcoords,particles=particles,dt=dt)
		endif
		dims = size(magy, /dimensions)
		if var eq 'magy' then begin
			if keyword_set(memefficient) then slice = temporary(magy) else slice = magy
		endif
	endif
	if total(var eq ['magz','valf','imagp','imagbeta','divbratio']) ge 1 then begin
		if n_elements(magz) eq 0 then begin
			magz = jloaddata(filename,'magz',xrange=xrange,yrange=yrange,zrange=zrange,sample=sample,lwant=lwant,time=time,xcoords=xcoords,ycoords=ycoords,zcoords=zcoords,particles=particles,dt=dt)
		endif
		dims = size(magz, /dimensions)
		if var eq 'magz' then begin
			if keyword_set(memefficient) then slice = temporary(magz) else slice = magz
		endif
	endif
	if total(var eq ['divb','divbratio']) ge 1 then begin
		if n_elements(divb) eq 0 then begin
			divb = jloaddata(filename,'divb',xrange=xrange,yrange=yrange,zrange=zrange,sample=sample,lwant=lwant,time=time,xcoords=xcoords,ycoords=ycoords,zcoords=zcoords,particles=particles,dt=dt)
		endif
		dims = size(divb, /dimensions)
		if var eq 'divb' then begin
			if keyword_set(memefficient) then slice = temporary(divb) else slice = divb
		endif
	endif

	if n_elements(dims) gt 0 then begin
		if sliceplane eq 'x' then begin
			if n_elements(dims) eq 2 then dims = [1, dims[0], dims[1]]
		endif
		if sliceplane eq 'y' then begin
			if n_elements(dims) eq 2 then dims = [dims[0], 1, dims[1]]
		endif
		if sliceplane eq 'z' then begin
			if n_elements(dims) eq 2 then dims = [dims[0], dims[1], 1]
		endif
	endif

	;Calculate trajectory information for variables that need it.
	if var eq 'tide' or var eq 'tidex' or var eq 'tidey' or var eq 'tidez' or var eq 'bhbound' or $
		var eq 'gpresz_tidez' or var eq 'kineskew' or var eq 'selfbound' then begin
		jread_amr, filename, VAR_NAME='none', TREE=tree, DATA=unk, PARAMETERS=params

		if ~keyword_set(orbinfo) then begin
			print, 'Error: Set "orbinfo" equal to the mass of the point mass in grams.'
			return
		endif
		m = orbinfo
		time = params.time

		data = read_ascii('pruned_orbit.dat')
		data = data.field01
		t = data[0,*]
		dattime = min(abs(t - time), tindex)
		hxcm = data[1,tindex]
		hycm = data[2,tindex]
		hzcm = data[3,tindex]

		hvx = data[4,tindex]
		hvy = data[5,tindex]
		hvz = data[6,tindex]

		oxcm = data[7,tindex]
		oycm = data[8,tindex]
		ozcm = data[9,tindex]
				   
		ovx = data[10,tindex]
		ovy = data[11,tindex]

		mxcm = data[25,tindex]
		mycm = data[26,tindex]
		mzcm = data[27,tindex]

		mvx = data[28,tindex]
		mvy = data[29,tindex]
		mvz = data[30,tindex]

		dx = -mxcm + oxcm - hxcm
		dy = -mycm + oycm - hycm
		dz = -mzcm + ozcm - hzcm

		vx = -mvx + ovx - hvx
		vy = -mvy + ovy - hvy
	endif

	;Now create derived variables if necessary
	if var eq 'entropy' then begin ;this is the entropy per unit mass
		slice = kb/mp*(alog(mp/double(dens)*(2*mp*!pi*kb*double(temp)/h^2)^(3.0/2.0)) + (5.0/2.0))
	endif
	if var eq 'velxy' then begin 
		slice = sqrt(velx^2 + vely^2)
	endif
	if var eq 'rigid_deviation' then begin 
		if n_elements(velx) eq 0 then $
			velx = jloaddata(filename,'velx',xrange=xrange,yrange=yrange,zrange=zrange,sample=sample,lwant=lwant,time=time,xcoords=xcoords,ycoords=ycoords,particles=particles,dt=dt)
		if n_elements(vely) eq 0 then $
			vely = jloaddata(filename,'vely',xrange=xrange,yrange=yrange,zrange=zrange,sample=sample,lwant=lwant,time=time,particles=particles,dt=dt)
		xcoords = xcoords - simsize/2.
		ycoords = ycoords - simsize/2.
		dims = size(velx, /dimensions)
		xcoordarr = dblarr(dims)
		for j=0,dims[0]-1 do xcoordarr[*,j] = xcoords
		ycoordarr = dblarr(dims)
		for j=0,dims[1]-1 do ycoordarr[j,*] = reverse(ycoords)
		distmat = sqrt(xcoordarr^2.+ycoordarr^2.)
		newvelx = -velx*ycoordarr/distmat
		newvely = vely*xcoordarr/distmat
		slice = sqrt(newvelx^2. + newvely^2.)
		slice = slice/distmat
		print, dims, size(slice, /dimensions)
	endif
	if var eq 'divbratio' then begin
		slice = divb / sqrt(magx^2 + magy^2 + magz^2)
	endif
	if var eq 'imagp' then begin
		slice = 0.5/(4.0*!pi)*(magx^2 + magy^2 + magz^2)
	endif
	if var eq 'valf' then begin
		; not sure on units
		magcon = 0.125663706
		slice = sqrt(magx^2 + magy^2 + magz^2)/sqrt(magcon*dens)
	endif
	if var eq 'gpotener' then begin
		slice = dens*gpot
	endif
	if var eq 'vtot' then begin
		if keyword_set(subtractavg) then begin
			velx = velx - total(velx)/n_elements(velx)
			vely = vely - total(vely)/n_elements(vely)
			velz = velz - total(velz)/n_elements(velz)
		endif
		slice = sqrt(velx^2 + vely^2 + velz^2)
	endif
	if var eq 'vtotxy' then begin
		if keyword_set(subtractavg) then begin
			velx = velx - total(velx)/n_elements(velx)
			vely = vely - total(vely)/n_elements(vely)
		endif
		slice = sqrt(velx^2 + vely^2)
	endif
	if var eq 'vtot2' then begin
		if keyword_set(subtractavg) then begin
			velx = velx - total(velx)/n_elements(velx)
			vely = vely - total(vely)/n_elements(vely)
			velz = velz - total(velz)/n_elements(velz)
		endif
		slice = velx^2 + vely^2 + velz^2
	endif
	if var eq 'momentum' then begin
		if keyword_set(subtractavg) then begin
			velx = velx - total(velx)/n_elements(velx)
			vely = vely - total(vely)/n_elements(vely)
			velz = velz - total(velz)/n_elements(velz)
		endif
		slice = dens*sqrt(velx^2. + vely^2. + velz^2.)
	endif
	if var eq 'kine' then begin
		if keyword_set(subtractavg) then begin
			velx = velx - total(velx)/n_elements(velx)
			vely = vely - total(vely)/n_elements(vely)
			velz = velz - total(velz)/n_elements(velz)
		endif
		slice = 0.5*dens*(velx^2. + vely^2. + velz^2.)
	endif
	if var eq 'kineskew' then begin
		if keyword_set(subtractavg) then begin
			velx = velx - total(velx)/n_elements(velx)
			vely = vely - total(vely)/n_elements(vely)
		endif
		slice = 0.5*dens*(velx*vx + vely*vy)
	endif
	if var eq 'eintener' then begin
		slice = dens*eint
	endif
	if var eq 'kinpresratio' then begin
		slice = (0.5*dens*(velx^2. + vely^2. + velz^2.))/pres
	endif
	if var eq 'kinpresdiff' then begin
		slice = 0.5*dens*(velx^2. + vely^2. + velz^2.)-pres
	endif
	if var eq 'ipres' then begin
		slice = dens * temp * 82544092.3
	endif
	if var eq 'imagbeta' then begin
		slice = dens * temp * 82544092.3 / (0.5/(4.0*!pi)*(magx^2 + magy^2 + magz^2))
	endif
	if var eq 'csnd' then begin
		gamc = reform(jloaddata(filename,'gamc',xrange=xrange,yrange=yrange,zrange=zrange,sample=sample,lwant=lwant,time=time,xcoords=xcoords,ycoords=ycoords,zcoords=zcoords,particles=particles,dt=dt))
		;gamc = 2.0
		slice = sqrt(gamc*pres/dens)
	endif
	if var eq 'cfl' then begin
		gamc = reform(jloaddata(filename,'gamc',xrange=xrange,yrange=yrange,zrange=zrange,sample=sample,lwant=lwant,time=time,xcoords=xcoords,ycoords=ycoords,zcoords=zcoords,particles=particles,dt=dt))
		csnd = sqrt(gamc*pres/dens)
		if n_elements(xcoords) gt 1 then begin
			dx = xcoords[1] - xcoords[0]
		endif else if n_elements(ycoords) gt 1 then begin
			dx = ycoords[1] - ycoords[0]
		endif else begin
			dx = zcoords[1] - zcoords[0]
		endelse

		slice = dx/(csnd + max(abs([[[velx]], [[vely]], [[velz]]]), dimension=3))
	endif
	if var eq 'jeans' then begin
		gamc = reform(jloaddata(filename,'gamc',xrange=xrange,yrange=yrange,zrange=zrange,sample=sample,lwant=lwant,time=time,xcoords=xcoords,ycoords=ycoords,zcoords=zcoords,particles=particles,dt=dt))
		slice = sqrt(!pi*gamc*pres/dens^2./g) ; lambda_jeans
		slice = (xcoords[1] - xcoords[0])/slice
	endif
	if var eq 'mach' then begin
		;slice = sqrt(velx^2 + vely^2 + velz^2) / sqrt(138574524.*temp)
		gam = 4./3.
		slice = sqrt(velx^2 + vely^2 + velz^2) / sqrt(gam*2.52192246e-15*pres/dens)
	endif
	if var eq 'machz' then begin
		slice = abs(velz) / sqrt(138574524.*temp)
	endif
	if var eq 'acom' or var eq 'acomx' or var eq 'acomy' or var eq 'acomz' then begin
		csize = xcoords[1]-xcoords[0]
		coeff = g*csize^3.
		dx = xcoords - simsize/2.
		dx = cmreplicate(dx, [dims[1], dims[2]])
		dy = ycoords - simsize/2.
		dy = transpose(cmreplicate(dy, [dims[0], dims[2]]), [1, 0, 2])
		dz = zcoords - simsize/2.
		dz = transpose(cmreplicate(dz, [dims[1], dims[0]]), [2, 1, 0])
		dr = sqrt(dx^2. + dy^2. + dz^2.)
		if var eq 'acom' then slice = coeff*dens/dr^2.
		if var eq 'acomx' then slice = coeff*dens*dx/dr^3.
		if var eq 'acomy' then slice = coeff*dens*dy/dr^3.
		if var eq 'acomz' then slice = coeff*dens*dz/dr^3.
	endif	
	if var eq 'bhbound' then begin
		coeff = g*m
		dr = (dx+xcoords)^2.
		dr = cmreplicate(dr, [dims[1], dims[2]])
		dy2 = (dy+ycoords)^2.
		temp1 = cmreplicate(temporary(dy2), [dims[0], dims[2]])
		if n_elements(size(temp1, /dimensions)) eq 2 then begin
			dr = temporary(dr) + transpose(temp1, [1, 0])
		endif else begin
			dr = temporary(dr) + transpose(temp1, [1, 0, 2])
		endelse
		dz2 = (dz+zcoords)^2.
		temp1 = cmreplicate(temporary(dz2), [dims[1], dims[0]])
		if n_elements(size(temp1, /dimensions)) eq 2 then begin
			dr = temporary(dr) + transpose(temp1, [1, 0])
		endif else begin
			dr = temporary(dr) + transpose(temp1, [2, 1, 0])
		endelse
		;if keyword_set(subtractavg) then begin
		;	avelx = velx - total(velx)/n_elements(velx)
		;	avely = vely - total(vely)/n_elements(vely)
		;	avelz = velz - total(velz)/n_elements(velz)
		;endif
		slice = 0.5*((vx+velx)^2 + (vy+vely)^2 + velz^2) - coeff/sqrt(temporary(dr))
		;slice = dens*(temporary(slice) - gpot)
	endif	
	if var eq 'selfbound' then begin
		;Changed this to ignore the black hole, the commented out parts are for the BH
		;coeff = g*m
		;dx2 = (dx+(findgen(dims[0])+0.5)*sqsize)^2.
		;dy2 = (dy+(findgen(dims[1])+0.5)*sqsize)^2.
		;dz2 = (dz+(findgen(dims[2])+0.5)*sqsize)^2.
		if keyword_set(subtractavg) then begin
			totmass = total(dens)
			avelx = velx - total(dens*velx)/totmass
			avely = vely - total(dens*vely)/totmass
			avelz = velz - total(dens*velz)/totmass
		endif else begin
			avelx = double(velx) - mvx
			avely = double(vely) - mvy
			avelz = double(velz) - mvz
		endelse
		slice = 0.5*(avelx^2 + avely^2 + avelz^2) + gpot
	endif	
	if var eq 'tide' or var eq 'tidex' or var eq 'tidey' or var eq 'tidez' or var eq 'gpresz_tidez' then begin
		d3 = double(norm([x,y]))^3.
		xd3 = x/d3
		yd3 = y/d3

		slice = double(g*m*dens)
		dx2 = (dx+(findgen(dims[0])+0.5)*sqsize)
		dy2 = (dy+(findgen(dims[1])+0.5)*sqsize)
		dz2 = (dz+(findgen(dims[2])+0.5)*sqsize)
		if var eq 'tide' then begin
			for j=0,dims[0]-1 do begin
				for k=0,dims[1]-1 do begin
					for l=0,dims[2]-1 do begin
						celld3 = (dx2[j]^2.+dy2[k]^2.+dz2[l]^2.)^(1.5)
						slice[j,k,l] = slice[j,k,l]*sqrt((dx2[j]/celld3-xd3)^2.+(dy2[k]/celld3-yd3)^2.+(dz2[l]/celld3)^2.)
					endfor
				endfor
			endfor
		endif
		if var eq 'tidez' or var eq 'gpresz_tidez' then begin
			for j=0,dims[0]-1 do begin
				for k=0,dims[1]-1 do begin
					for l=0,dims[2]-1 do begin
						slice[j,k,l] = abs(slice[j,k,l]*dz2[l]/(dx2[j]^2.+dy2[k]^2.+dz2[l]^2.)^1.5)
					endfor
				endfor
			endfor
		endif
	endif
	if var eq 'angvel' or var eq 'angvelx' or var eq 'angvely' then begin
		sqsize = (xrange[1] - xrange[0]) / dims[0]
		midx = (xrange[0] - xrange[1])/2.
		midy = (yrange[0] - yrange[1])/2.
		slice = fltarr(dims[0], dims[1])
		softening = 5*sqsize
		xpos = midx + (dindgen(dims[0])+0.5)*sqsize
		ypos = midy + (dindgen(dims[1])+0.5)*sqsize
		if var eq 'angvelx' then begin
			for i=0,dims[0]-1 do begin
				for j=0,dims[1]-1 do begin
					d2 = xpos[i]^2 + ypos[j]^2
					slice[i,j] = velx[i,j]/sqrt(softening^2.+d2)
				endfor
			endfor
		endif
		if var eq 'angvely' then begin
			for i=0,dims[0]-1 do begin
				for j=0,dims[1]-1 do begin
					d2 = xpos[i]^2 + ypos[j]^2
					slice[i,j] = vely[i,j]/sqrt(softening^2.+d2)
				endfor
			endfor
		endif
		if var eq 'angvel' then begin
			for i=0,dims[0]-1 do begin
				for j=0,dims[1]-1 do begin
					d2 = xpos[i]^2 + ypos[j]^2
					slice[i,j] = sqrt(velx[i,j]^2. + vely[i,j]^2.)/sqrt(softening^2.+d2)
				endfor
			endfor
		endif
	endif
	if var eq 'torque' then begin
		sqsize = (xrange[1] - xrange[0]) / dims[0]
		midx = 1.5e10
		midy = 1.5e10
		leverx = 1.2e10
		levery = 1.5e10
		fulcrumx = 1.4e10
		fulcrumy = 1.5e10
		slice = fltarr(dims[0], dims[1])

		for i=0,dims[0]-1 do begin
			for j=0,dims[1]-1 do begin
				xpos = double(xrange[0]+i*sqsize-midx)
				ypos = -double(yrange[0]+j*sqsize-midy)
				xpos2 = double(xrange[0]+i*sqsize-leverx)
				ypos2 = -double(yrange[0]+j*sqsize-levery)
				xpos3 = double(xrange[0]+i*sqsize-fulcrumx)
				ypos3 = -double(yrange[0]+j*sqsize-fulcrumy)
				r = sqrt(xpos^2. + ypos^2.)
				r2 = sqrt(xpos2^2. + ypos2^2.)
				r3 = sqrt(xpos3^2. + ypos3^2.)
				slice[i,j] = dens[i,j]*g*((xpos3*ypos2 - ypos3*xpos2)/r2^3. + (xpos3*ypos - ypos3*xpos)/r^3.)
			endfor
		endfor
		slice = abs(slice)
	endif
	if var eq 'torqmom' then begin
		sqsize = (xrange[1] - xrange[0]) / dims[0]
		midx = 1.5e10
		midy = 1.5e10
		leverx = 1.165e10
		levery = 1.5e10
		fulcrumx = 1.39e10
		fulcrumy = 1.5e10
		slice = fltarr(dims[0], dims[1])

		for i=0,dims[0]-1 do begin
			for j=0,dims[1]-1 do begin
				xpos = double(xrange[0]+i*sqsize-midx)
				ypos = double(yrange[0]+j*sqsize-midy)
				xpos2 = double(xrange[0]+i*sqsize-leverx)
				ypos2 = double(yrange[0]+j*sqsize-levery)
				xpos3 = double(xrange[0]+i*sqsize-fulcrumx)
				ypos3 = double(yrange[0]+j*sqsize-fulcrumy)
				r = sqrt(xpos^2. + ypos^2.)
				r2 = sqrt(xpos2^2. + ypos2^2.)
				r3 = sqrt(xpos3^2. + ypos3^2.)
				slice[i,j] = dens[i,j]*g*((xpos3*ypos2 - ypos3*xpos2)/r2^3.)*(xpos*vely[i,j] - ypos*velx[i,j])
			endfor
		endfor
		newslice = fltarr(dims[0], dims[1]/2)
		newslice = slice[0:dims[0]-1, 0:dims[1]/2-1] + reverse(slice[0:dims[0]-1, dims[1]/2:dims[1]-1],2)
		slice = abs(newslice)
		dims = size(slice, /dimensions)
		return
	endif
	if var eq 'angmom' then begin
		sqsize = (xrange[1] - xrange[0]) / dims[0]
		midx = (xrange[0] - xrange[1])/2.
		midy = (yrange[0] - yrange[1])/2.
		xpos = midx + (dindgen(dims[0])+0.5)*sqsize
		ypos = midy + (dindgen(dims[1])+0.5)*sqsize
		if n_elements(dims) eq 2 then begin
			slice = dblarr(dims[0], dims[1])
		endif else begin
			slice = dblarr(dims[0], dims[1], dims[2])
		endelse

		for i=0,dims[0]-1 do begin
			for j=0,dims[1]-1 do begin
				if n_elements(dims) eq 2 then begin
					slice[i,j] = dens[i,j]*(xpos[i]*vely[i,j] - ypos[j]*velx[i,j])
				endif else begin
					for k=0,dims[2]-1 do begin
						slice[i,j,k] = dens[i,j,k]*(xpos[i]*vely[i,j,k] - ypos[j]*velx[i,j,k])
					endfor
				endelse
			endfor
		endfor
	endif
	if var eq 'mominertia' then begin
		sqsize = (xrange[1] - xrange[0]) / dims[0]
		midx = (xrange[0] - xrange[1])/2.
		midy = (yrange[0] - yrange[1])/2.
		xpos = midx + (dindgen(dims[0])+0.5)*sqsize
		ypos = midy + (dindgen(dims[1])+0.5)*sqsize
		if n_elements(dims) eq 2 then begin
			slice = dblarr(dims[0], dims[1])
		endif else begin
			slice = dblarr(dims[0], dims[1], dims[2])
		endelse

		for i=0,dims[0]-1 do begin
			for j=0,dims[1]-1 do begin
				if n_elements(dims) eq 2 then begin
					slice[i,j] = dens[i,j]*(xpos[i]^2. + ypos[j]^2.)
				endif else begin
					for k=0,dims[2]-1 do begin
						slice[i,j,k] = dens[i,j,k]*(xpos[i]^2. + ypos[j]^2.)
					endfor
				endelse
			endfor
		endfor
	endif
	if var eq 'angmomdens' then begin
		sqsize = (xrange[1] - xrange[0]) / dims[0]
		midx = simsize/2.
		midy = simsize/2.
		slice = fltarr(dims[0], dims[1])
		for i=0,dims[0]-1 do begin
			for j=0,dims[0]-1 do begin
				xpos = double(xrange[0]+j*sqsize-midx)
				ypos = double(yrange[0]+i*sqsize-midy)
				slice[i,j] = (xpos*vely[i,j] - ypos*velx[i,j])
			endfor
		endfor
	endif
	if var eq 'shock' then begin
		csnd = slice
		csnd = reform(csnd)
		slice = dblarr(dims[0], dims[1]) 
		for i=0,dims[0]-1 do begin
			for j=0,dims[1]-1 do begin
				vx = velx[i,j]
				vy = vely[i,j]
				vz = velz[i,j]
				cs = csnd[i,j]
				ss = 0.3
				if  (i-1 ge 0 && sqrt((velx[i-1,j]-vx)^2 + (vely[i-1,j]-vy)^2 + (velz[i-1,j]-vz)^2)/cs ge ss) or $
					(j-1 ge 0 && sqrt((velx[i,j-1]-vx)^2 + (vely[i,j-1]-vy)^2 + (velz[i,j-1]-vz)^2)/cs ge ss) or $
					(i-1 ge 0 && j-1 ge 0 && $
						sqrt((velx[i-1,j-1]-vx)^2 + (vely[i-1,j-1]-vy)^2 + (velz[i-1,j-1]-vz)^2)/cs ge ss) or $
					(i+1 lt dims[0] && sqrt((velx[i+1,j]-vx)^2 + (vely[i+1,j]-vy)^2 + (velz[i+1,j]-vz)^2)/cs ge ss) or $
					(j+1 lt dims[1] && sqrt((velx[i,j+1]-vx)^2 + (vely[i,j+1]-vy)^2 + (velz[i,j+1]-vz)^2)/cs ge ss) or $
					(i+1 lt dims[0] && j+1 lt dims[1] && $
						sqrt((velx[i+1,j+1]-vx)^2 + (vely[i+1,j+1]-vy)^2 + (velz[i+1,j+1]-vz)^2)/cs ge ss) or $
					(i-1 ge 0 && j+1 lt dims[1] && $
						sqrt((velx[i-1,j+1]-vx)^2 + (vely[i-1,j+1]-vy)^2 + (velz[i-1,j+1]-vz)^2)/cs ge ss) or $
					(j-1 ge 0 && i+1 lt dims[0] && $
						sqrt((velx[i+1,j-1]-vx)^2 + (vely[i+1,j-1]-vy)^2 + (velz[i+1,j-1]-vz)^2)/cs ge ss) $
					then begin
						slice[i,j]=1.0
				endif
			endfor
		endfor
	endif
	if var eq 'pp' or var eq 'cno' or var eq 'he' or var eq 'nuc' then begin
		t6 = temp / 1e6
		if var eq 'pp' or var eq 'nuc' then begin
			pp = 2.38e6 * (1 + 0.0123*t6^(1./3.) + 0.0109*t6^(2./3) + 0.0009*t6) * 0.7^2.$
				* dens * t6^(-2./3) * exp(-33.8*t6^(-1./3)) * dens
		endif
		if var eq 'cno' or var eq 'nuc' then begin
			cno = 8.67e27 * (1 + 0.0027*t6^(1./3.) - 0.0078*t6^(2./3) - 0.000149*t6) * (0.01) * (0.7)$
				* dens * t6^(-2./3) * exp(-152.28*t6^(-1./3)) * dens
			cno = reform(cno)
			for i=0,dims[0]-1 do begin
				for j=0,dims[1]-1 do begin
					if cno[i,j] le 0.0 then cno[i,j] = 1e-20 
				endfor
			endfor
		endif
		if var eq 'he' or var eq 'nuc' then begin
			t8 = t6 / 1e2
			he = 5.09e11 * (0.3)^3. * dens*dens / t8^3. * exp(-44.027/t8) * dens
			he = reform(he)
			for i=0,dims[0]-1 do begin
				for j=0,dims[1]-1 do begin
					if he[i,j] le 0.0 then he[i,j] = 1e-20 
				endfor
			endfor
		endif
		case 1 of
			var eq 'pp': slice = pp
			var eq 'cno': slice = cno
			var eq 'he': slice = he
			else: slice = pp + cno + he
		endcase
	endif
	if var eq 'gdensz' then begin
		print, size(dens, /dimensions)
		slice = fltarr(dims[1], dims[2])
		for i=0,dims[1]-1 do begin
			for j=0,dims[2]-1 do begin
				if j ne 0 and j ne dims[2]-1 then begin
					slice[i,j] = (dens[0,i,j+1] - dens[0,i,j-1]) / dens[0,i,j]
				endif
			endfor
		endfor
	endif
	if var eq 'gvelzz' then begin
		slice = fltarr(dims[0],dims[1])
		for i=0,dims[0]-1 do begin
			for j=0,dims[1]-1 do begin
				if j ne 0 and j ne dims[1]-1 then begin
					slice[i,j] = (velz[i,j+1] - velz[i,j-1])
				endif
			endfor
		endfor
	endif
	if var eq 'gmachz' then begin
		csnd = sqrt(138574524.*temp)
		mach = sqrt(velx^2 + vely^2 + velz^2) / csnd
		slice = fltarr(dims[0],dims[1])
		for i=0,dims[0]-1 do begin
			for j=0,dims[1]-1 do begin
				if j ne 0 and j ne dims[1]-1 then begin
					slice[i,j] = (mach[i,j+1] - mach[i,j-1])
				endif
			endfor
		endfor
	endif
	if var eq 'gpresz' or var eq 'gpresz_tidez' then begin
		if var eq 'gpresz_tidez' then begin
			tidezslice = slice
			tidezslice[*,*,dims[2]-1] = tidezslice[*,*,dims[2]-2]
			tidezslice[*,*,0] = tidezslice[*,*,1]
		endif
		pres = dens * temp * 82544092.3
		slice = dblarr(dims)
		for i=0,dims[0]-1 do begin
			for j=0,dims[1]-1 do begin
				for k=1,dims[2]-2 do begin
					slice[i,j,k] = abs(pres[i,j,k+1] - pres[i,j,k]) / sqsize
				endfor
			endfor
		endfor
	endif
	if var eq 'gpresyzmag' then begin
		pres = dens * temp * 82544092.3
		slice = fltarr(dims[0], dims[1])
		for i=1,dims[0]-2 do begin
			for j=1,dims[1]-2 do begin
				if ~((pres[i,j+1] - pres[i,j-1]) eq 0 or (pres[i+1,j] - pres[i-1,j]) eq 0) then begin
					slice[i,j] = sqrt((pres[i,j+1] - pres[i,j-1])^2 + (pres[i+1,j] - pres[i-1,j])^2) / pres[i,j]
				endif
			endfor
		endfor
	endif

	if var eq 'gpresyzmag' then begin
		for i=1,dims[0]-2 do begin
			for j=1,dims[1]-2 do begin
				if slice[i,j] eq 0 then begin
					for k=i-1,1,-1 do begin
						if slice[k,j] ne 0 then begin
							slice[i,j] = slice[k,j]
							break
						endif
					endfor
				endif
			endfor
		endfor
		for i=1,dims[0]-2 do begin
			for j=1,dims[1]-2 do begin
				if slice[i,j] eq 0 then begin
					for k=j-1,1,-1 do begin
						if slice[i,k] ne 0 then begin
							slice[i,j] = slice[i,k]
							break
						endif
					endfor
				endif
			endfor
		endfor
					;if ~found then begin
					;	for k=i-1,1,-1 do begin
					;		if slice[k,j] ne 0 then begin
					;			slice[i,j] = slice[k,j]
					;			found = 1
					;			break
					;		endif
					;	endfor
					;	if ~found then begin
					;	for k=i
					;endif
		slice[0,1:dims[1]-2] = slice[1,1:dims[1]-2]
		slice[dims[0]-1,1:dims[1]-2] = slice[dims[0]-2,1:dims[1]-2]
		slice[*,0] = slice[*,1]
		slice[*,dims[1]-1] = slice[*,dims[1]-2]
	endif
	if var eq 'gpresz' or var eq 'gpresz_tidez' then begin
		for i=0,dims[0]-1 do begin
			for j=0,dims[1]-1 do begin
				for k=1,dims[2]-2 do begin
					if slice[i,j,k] eq 0 then begin
						for l=k+1,dims[2]-2 do begin
							if slice[i,j,l] ne 0 then begin
								slice[i,j,k:l] = slice[i,j,l]/double(l-k+1)
								break
							endif
						endfor
						if slice[i,j,k] eq 0 then slice[i,j,k] = slice[i,j,k-1]
					endif
					;if slice[i,j] ne 0 then begin
					;	for k=j-1,1,-1 do begin
					;		if slice[i,k] ne 0 then begin
					;			d = j-k
					;			if d gt 1 then slice[i,k:j] = slice[i,j] / (d + 1.)
					;		endif
					;	endfor
					;endif
					;if slice[i,j] eq 0 then begin
					;	for k=j-1,1,-1 do begin
					;		if slice[i,k] ne 0 then begin
					;			slice[i,j] = slice[i,k]
					;			break
					;		endif
					;	endfor
					;endif
				endfor
			endfor
		endfor
		slice[*,*,dims[2]-1] = slice[*,*,dims[2]-2]
		slice[*,*,0] = slice[*,*,1]
	endif
	;print, slice
	;print, tidezslice
	if var eq 'gpresz_tidez' then slice = double(slice) / double(tidezslice)

	if n_elements(slice) eq 0 then begin
		varname = var ;jloaddata will alter var to be 4 characters, this is not desired
		slice = reform(jloaddata(filename,varname,xrange=xrange,yrange=yrange,zrange=zrange,sample=sample,lwant=lwant,time=time,xcoords=xcoords,ycoords=ycoords,zcoords=zcoords,particles=particles,dt=dt))
		dims = size(slice, /dimensions)
	endif

	if n_elements(temp) ne 0 then temp = reform(temp, dims)
	if n_elements(dens) ne 0 then dens = reform(dens, dims)
	if n_elements(pres) ne 0 then pres = reform(pres, dims)
	if n_elements(velx) ne 0 then velx = reform(velx, dims)
	if n_elements(vely) ne 0 then vely = reform(vely, dims)
	if n_elements(velz) ne 0 then velz = reform(velz, dims)
	if n_elements(gpot) ne 0 then gpot = reform(gpot, dims)

	if keyword_set(special) then begin
		if special eq 'diff' then begin
			mc = double(xcoords[n_elements(xcoords)-1] - xcoords[0])/n_elements(xcoords)
			xcoords = (dindgen(n_elements(xcoords)) + 0.5) * mc
			ycoords = (dindgen(n_elements(xcoords)) + 0.5) * mc
			x0 = (xr[0] - refcoor[0,0]) mod mc
			y0 = (yr[0] - refcoor[0,1]) mod mc
			nxc = ((indgen(n_elements(xcoords))) * mc + x0)/mc
			nyc = ((indgen(n_elements(xcoords))) * mc + y0)/mc

			slice = interpolate(slice, nxc, nyc, /grid, cubic=-0.5)
			;slice = slice/base_state - 1.0
			slice = slice-base_state
			if keyword_set(log) then slice = -slice
		endif
	endif

	slice = reform(slice, dims)
end
