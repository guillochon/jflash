function wave_kernel, x, y, z, xc, yc, zc, a
	;if sqrt((x-xc)^2.+(y-yc)^2.+(z-zc)^2.) gt a then return, 0
	;return, 1
	return, (3.0 - ((x-xc)^2. + (y-yc)^2. + (z-zc)^2.)/a^2.)*exp(-((x-xc)^2.+(y-yc)^2.+(z-zc)^2.)/(2.*a^2.))
end

pro flash_wavelet,filename,var,xrange,yrange,zrange,$
    simsize=simsize,sample=sample,fprefix=fprefix,extdata=extdata,datatime=datatime,$
    log=log,annotatepos=annotatepos,output=output,special=special,subtractavg=subtractavg

    compile_opt idl2
    fname = filename
	if n_elements(simsize) eq 0 then simsize = xrange[1] - xrange[0]
	if n_elements(fprefix) eq 0 then fprefix = ''

	load_flash_var, slice, filename, var, xrange, yrange, zrange, dens=dens, temp=temp, $
		velx=velx, vely=vely, velz=velz, gpot=gpot, log=log, sample=sample, time=time, simsize=simsize, subtractavg=subtractavg, $
		xcoords=xcoords, ycoords=ycoords, zcoords=zcoords

	slice = reform(slice)
	gc = xcoords[1]-xcoords[0]

	slice_dims = size(slice, /dimensions)
	print, slice_dims
	if(size(slice_dims,/n_elements) eq 2) then begin
		slice_dims = [slice_dims, 1]
	endif

	smat = dblarr(floor(alog(slice_dims[0])/alog(2.0)))
	kwave = dblarr(floor(alog(slice_dims[0])/alog(2.0)))
	
	for s=0,n_elements(smat)-1 do begin
		ss = 2^s
		ss1 = 2^s-1
		tot = 0
		cnt = 0
		for i=0,ss1 do begin
		for j=0,ss1 do begin
		for k=0,ss1 do begin
			iw = ceil(double(i)/ss*slice_dims[0])
			iw1 = ceil(double(i+1)/ss*slice_dims[0])-1
			jw = ceil(double(j)/ss*slice_dims[1])
			jw1 = ceil(double(j+1)/ss*slice_dims[1])-1
			kw = ceil(double(k)/ss*slice_dims[2])
			kw1 = ceil(double(k+1)/ss*slice_dims[2])-1
			subvelx = velx[iw:iw1,jw:jw1,kw:kw1] - total(velx[iw:iw1,jw:jw1,kw:kw1])/n_elements(velx[iw:iw1,jw:jw1,kw:kw1])
			subvely = vely[iw:iw1,jw:jw1,kw:kw1] - total(vely[iw:iw1,jw:jw1,kw:kw1])/n_elements(vely[iw:iw1,jw:jw1,kw:kw1])
			subvelz = velz[iw:iw1,jw:jw1,kw:kw1] - total(velz[iw:iw1,jw:jw1,kw:kw1])/n_elements(velz[iw:iw1,jw:jw1,kw:kw1])
			;tot = tot + total(subvelx^2. + subvelz^2.)
			tot = tot + total(subvelx^2. + subvely^2. + subvelz^2.)
			cnt = cnt + 1
		endfor
		endfor
		endfor
		smat[s] = tot*slice_dims[0]*gc/ss/(slice_dims[0]^3.)
		kwave[s] = 2*!pi/(slice_dims[0]*gc/ss)
	endfor
	
	;smat = dblarr(floor(min(slice_dims)/2.))
	;kwave = dblarr(floor(min(slice_dims)/2.))

	;for bc=1,n_elements(smat)-1 do begin
	;	sbc = 2.*bc+1
	;	kern = make_array((2*bc+1)*[1,1,1],value=1.0,/double)
	;	;kern = dblarr((2*bc+1)*[1,1,1])
	;	;for ii=-bc,bc do begin 
	;	;	for jj=-bc,bc do begin 
	;	;		for kk=-bc,bc do begin 
	;	;			kern[ii+bc,jj+bc,kk+bc] = wave_kernel(ii,jj,kk,0,0,0,(bc+1)/2.)
	;	;		endfor
	;	;	endfor
	;	;endfor

	;	print, 100.*bc/n_elements(smat), '%'
	;	velxmat = velx - convol(velx, kern, sbc^3.)
	;	velymat = vely - convol(vely, kern, sbc^3.)
	;	velzmat = velz - convol(velz, kern, sbc^3.)
	;	vtot2mat = velxmat^2. + velymat^2. + velzmat^2.
	;	cmat = convol(vtot2mat, kern, sbc^3.)
	;	smat[bc] = total(cmat)/sbc
	;	smat[bc] = total(cmat)/total(where(cmat))/sbc
	;	kwave[bc] = 2*!pi/(sbc*gc)
	;endfor
	;print, kwave
	;print, smat

	kwave = reverse(kwave)
	smat = reverse(smat)
	pfit = poly_fit(alog10(kwave), alog10(smat), 1)
	set_plot, 'x'
	plot, alog10(kwave), alog10(smat), xsty=1+2, ysty=1+2
	print, pfit[1]
	;oplot, kwave, pfit[0] + pfit[1]*kwave + pfit[2]*kwave^2.0, linestyle=2
end

