; make_flash_frames.pro
; Generates a series of frames from one or more FLASH data files.
;
; Written by James Guillochon, jfg@ucolick.org
; (LAST UPDATED 5/04/2012)
;
; PARAMETERS
; **Required**
; basename      (str)             - Base file name, e.g. 'simulation_hdf5_plt_cnt'.
; start         (int)             - Which file number to start with.
; finish        (int)             - Which file number to finish with.
; var           (str)             - Variable to generate frames from, e.g. 'dens'.
; xrange        ([dbl, dbl])      - 2 values indicating range of x values to extract.
; yrange        ([dbl, dbl])      - 2 values indicating range of y values to extract.
; zrange        ([dbl, dbl])      - 2 values indicating range of z values to extract.
; sample        (int)             - Number of times to downsample data from highest refinment, must be <= 3.
; simsize       (dbl)             - Width of the box in cgs units (assumes box is cubical, used for certain variables).
; stride        (int)             - Stride of frame iteration, default 1.
; ps            (bool)            - Output plot to .ps file.
; ytitle        (str)             - y-axis label.
; thrvar        ([str,str,...])   - Variables used as a "threshold," data in regions where thresholds aren't
;                                   satisfied are set to rangemin.
; thrval        ([dbl,dbl,...])   - Threshold values, see above.
; thrtype       ([str,str,...])   - Whether the threshold is a minimum or a maximum, valid values 'min' or 'max'.
; volnorm       (bool)            - Multiply sum by the volume.
; fprefix       (str)             - Prefix to append to output image files, e.g. 'myprefix_restoffilename'
; log           (bool)            - Scale data logarithmically.
; subtractavg   (bool)            - Subtract the average value from all grid cells when loading variable.
; usetime       (bool)            - Multiply variable by dt.
; sum           (bool)            - Sum y values and print.
; indexlength   (int)             - Number of digits in filename index, e.g. _0000 is 4, _00000 is 5, etc.
; memefficient  (bool)            - Deallocate variables as soon as they are not needed. Incompatible with revolve plots.
; subdivide     (int)             - Divide volume into subchunks. Useful for large volumes that cannot be loaded all at once.
; trackfile     (str)             - File with x-y coordinates used to keep view centered.

pro flash_sum_var,basename,start,finish,var,xrange,yrange,zrange,sample=sample,simsize=simsize,stride=stride,ps=ps,$
	ytitle=ytitle,thrvar=thrvar,thrval=thrval,volnorm=volnorm,cellsize=cellsize,fprefix=fprefix,log=log,thrtype=thrtype,$
	subtractavg=subtractavg,usetime=usetime,sum=sum,indexlength=indexlength,memefficient=memefficient,subdivide=subdivide,$
	retx=retx,rety=rety,orbinfo=orbinfo,trackfile=trackfile
	
    xtitle = 'Time (s)'
    if ~n_elements(ytitle) then ytitle = var[0]
	if ~n_elements(fprefix) then fprefix = var[0]
	if ~keyword_set(indexlength) then indexlength = 4
	if ~keyword_set(thrtype) then thrtype = 'min'
	if ~keyword_set(subdivide) then subdivide = 1

    !p.color = 255
    !p.background = 0
    if n_elements(stride) eq 0 then stride = 1
    if ~keyword_set(ps) then begin
        set_plot, 'x'	
        window,1,xsize=800,ysize=800
    endif
	if n_elements(trackfile) ne 0 then begin
		track = read_ascii(trackfile)
		trackt = track.field1[0,*]
		trackx = track.field1[1,*]
		tracky = track.field1[2,*]
		trackz = track.field1[3,*]
	endif
	for i=start,finish,stride do begin
		num = strn(i, length=indexlength, padtype=1, padchar='0')
		filename = basename + '_' + num
        if i gt start then begin
			told = time
            yold = y
            xold = x
            l = floor(double(i-start)/stride+1)
            y = dblarr(l,n_elements(var))
            x = dblarr(l)
            y[0:(l-2),*] = yold
            x[0:(l-2)] = xold
        endif else begin
            y = dblarr(1,n_elements(var))
            x = dblarr(1)
			if keyword_set(usetime) then begin
				load_flash_var, slice, filename, 'dens', xrange, yrange, zrange, time=time, sample=3
				continue
			endif
        endelse
		if keyword_set(ps) then begin
			jps_start, filename=fprefix+'_sum_'+basename+'.eps',xsize=10.5,ysize=8.5,xoffset=10.75,yoffset=0.25,/inches
		endif

		read_amr, filename, var_name='dens', parameters=params
		if n_elements(trackfile) ne 0 then begin
			bxr = xrange + interpol(trackx, trackt, params.time)
			byr = yrange + interpol(tracky, trackt, params.time)
			bzr = zrange + interpol(trackz, trackt, params.time)
		endif else begin
			bxr = xrange
			byr = yrange
			bzr = zrange
		endelse

		xrsize = (xrange[1] - xrange[0])/subdivide
		yrsize = (yrange[1] - yrange[0])/subdivide
		zrsize = (zrange[1] - zrange[0])/subdivide

		y[(i-start)/stride,*] = 0
		for sdx=0,subdivide-1 do begin
		for sdy=0,subdivide-1 do begin
		for sdz=0,subdivide-1 do begin
			for v=0,n_elements(var)-1 do begin
				xr = bxr[0]+[0, xrsize] + xrsize*sdx
				yr = byr[0]+[0, yrsize] + yrsize*sdy
				zr = bzr[0]+[0, zrsize] + zrsize*sdz

				load_flash_var, slice, filename, var[v], xr, yr, zr, time=time, sample=sample, simsize=simsize, $
					dens=dens, temp=temp, velx=velx, vely=vely, velz=velz, pres=pres, gpot=gpot, xcoords=xcoords, ycoords=ycoords, zcoords=zcoords,$
					subtractavg=subtractavg, memefficient=memefficient, orbinfo=orbinfo
				if keyword_set(usetime) then slice = slice*(time-told)
				slice_dims = size(slice, /dimensions)
				if (n_elements(slice_dims) eq 2) then begin
					slice = reform(slice, slice_dims[0], slice_dims[1], 1)
					slice_dims = size(slice, /dimensions)
				endif
				thrslice = fltarr(slice_dims[0],slice_dims[1],slice_dims[2],n_elements(uniq(thrvar)))
				for t=0,n_elements(thrvar)-1 do begin
					if var[v] eq thrvar[t] then begin
						thrslice[*,*,*,t] = slice
					endif else begin
						thrloaded = 0
						for j=0,t-1 do begin
							if thrvar[t] eq thrvar[j] then begin
								thrslice[*,*,*,t] = thrslice[*,*,*,j]
								thrloaded = 1
							endif
						endfor
						if thrloaded eq 0 then begin
							load_flash_var, newthrslice, filename, thrvar[t], xr, yr, zr, dens=dens, temp=temp, pres=pres, $
								velx=velx, vely=vely, velz=velz, gpot=gpot, log=log, sample=sample, simsize=simsize, subtractavg=subtractavg
							thrslice[*,*,*,t] = temporary(newthrslice)
						endif	
					endelse
				endfor
				dims = size(slice, /dim)
				if n_elements(volnorm) eq 1 then begin
					if n_elements(cellsize) eq 1 then begin
						vol = cellsize^3.
					endif else begin
						vol = double(xcoords[1]-xcoords[0])^3.0
					endelse
				endif else begin
					vol = 1.0
				endelse
				if v eq 0 then x[(i-start)/stride] = time
				if n_elements(thrvar) ne 0 then begin
					thrind = make_array(n_elements(slice), val=1, /integer)
					for t=0,n_elements(thrvar)-1 do begin
						if thrtype[t] eq 'min' then thrind = logical_and(thrind, thrslice[*,*,*,t] ge thrval[t]) $
							else thrind = logical_and(thrind, thrslice[*,*,*,t] le thrval[t])
					endfor
					y[(i-start)/stride,v] = y[(i-start)/stride,v] + vol*total(slice[where(thrind, /l64)],/double)
				endif else y[(i-start)/stride,v] = y[(i-start)/stride,v] + vol*total(slice,/double)
				print, y[(i-start)/stride,v]
				if keyword_set(sum) then print, total(y, /double)
				if v eq 0 then begin
					if ~keyword_set(log) then begin
						plot, x, y[*,v], charsize=1.5, xtitle=xtitle, ytitle=ytitle, /xsty, ysty=1+2, yrange=[0,max(y)]
					endif else begin
						plot, x, alog10(y[*,v]), charsize=1.5, xtitle=xtitle, ytitle=ytitle, /xsty, ysty=1+2, yrange=[min(alog10(y)),max(alog10(y))]
					endelse
				endif else begin
					oplot, x, y[*,v], linestyle=v
				endelse
			endfor
			if (subdivide gt 1 or start ne finish) then fsc_undefine, slice, dens, temp, velx, vely, velz, gpot, pres
		endfor
		endfor
		endfor
		retx = x
		rety = y
        if keyword_set(ps) then jps_end
	endfor
end
