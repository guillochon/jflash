; make_flash_slice.pro
; Generates a plot of a slice from a FLASH data file.
;
; Written by James Guillochon, jfg@ucolick.org
; (LAST UPDATED 5/14/2009)
;
; PARAMETERS
; **Required**
; filename      (str)           - File name, e.g. 'simulation_hdf5_plt_cnt_0000'.
; var           (str)           - Variable to generate frames from, e.g. 'dens'.
; my_ct         (int)           - IDL color table to use, 1-45 are valid selections.
; xrange        ([dbl, dbl])    - 2 values indicating range of x values to extract.
; yrange        ([dbl, dbl])    - 2 values indicating range of y values to extract.
; zrange        ([dbl, dbl])    - 2 values indicating range of z values to extract.
; NOTE: One of either xrange, yrange, or zrange must have min = max, or must be a single value.
;
; **Optional**
; simsize       (dbl)               - Width of the box in cgs units (assumes box is cubical, used for certain variables).
; slicetype     (str)               - UNUSED.
; rangemin      (dbl)               - Minimum value to plot.
; rangemax      (dbl)               - Maximum value to plot.
; imgsizex      (dbl)               - Size of x dimension of output image, in pixels (applies to PNG only).
; imgsizey      (dbl)               - Size of y dimension of output image, in pixels (applies to PNG only).
; fprefix       (str)               - Prefix to append to output image files, e.g. 'myprefix_restoffilename'
; exactsize     (bool)              - Smallest grid cells will span exactly 1 pixel in output image.
; exactmult     (dbl)               - Exactsize must be true; smallest grid cells will span 1 pixel * exactmult.
; xticks        (int)               - Number of ticks to draw on x axis
; yticks        (int)               - Number of ticks to draw on y axis.
; zticks        (int)               - Number of ticks to draw on z axis (3D only).
; log           (bool)              - Scale data logarithmically.
; colorbarcolor (str)               - Color of colorbar annotation. 
; hidecolorbar  (bool)              - Hides the colorbar.
; contours      ([str,dbl,dbl,int]) - Parameters of contours to overplot, 'varname', min_cont_val, max_cont_val,
;	                                  nconts
; thrvar        ([str,str,...])     - Variables used as a "threshold," data in regions where thresholds aren't
;                                     satisfied are set to rangemin.
; thrval        ([dbl,dbl,...])     - Threshold values, see above.
; thrtype       ([str,str,...])     - Whether the threshold is a minimum or a maximum, valid values 'min' or 'max'.
; sample        (int)               - Number of times to downsample data from highest refinment, must be <= 3.
; stride        (int)               - Stride of frame iteration, default 1.
; charsize      (dbl)               - Multiplier of default char size, default 1.0.
; symrange      (bool)              - Makes the data range symmetric depending on min and max values of data.
; lmin          (dbl)               - Used with symrange only, specifies min log value for
;                                     symmetric ranges (Since log can't have negative argument).
; annotatepos   (str)               - Annotation position (colorbar, timestamp). Currently can only set to 'flip'.
; output        (str)               - Specifies output file type, can be 'ps' or 'png' (default).
; special       (str)               - Flag used to indicate "special" plots to generate.
; cellsize      (dbl)               - Size of the smallest grid cells in CGS units, used for some variables.
; memefficient  (bool)              - Deallocate variables as soon as they are not needed.
; hideaxes      (bool)              - Do not show axes.
; showblocks    (bool)              - Show block boundaries
; regrid
; gausswidth

pro make_flash_slice,filename,var,my_ct,xr,yr,zr,$
    simsize=simsize,slicetype=slicetype,rangemin=rangemin,$
    rangemax=rangemax,pos=pos,contours=contours,fieldvarx=fieldvarx,fieldvary=fieldvary,fieldmax=fieldmax,$
    thrvar=thrvar,thrval=thrval,thrtype=thrtype,sample=sample,lwant=lwant,imgsizex=imgsizex,imgsizey=imgsizey,$
    exactsize=exactsize,fprefix=fprefix,extdata=extdata,datatime=datatime,$
    log=log,colorbarcolor=colorbarcolor,charsize=charsize,hidetime=hidetime,$
    hidecolorbar=hidecolorbar,ambval=ambval,exactmult=exactmult,symrange=symrange,lmin=lmin,xticks=xticks,yticks=yticks,$
	annotatepos=annotatepos,output=output,special=special,hideaxes=hideaxes,negative=negative,subtractavg=subtractavg,$
	ctswitch=ctswitch,excision=excision,product=product,refcoor=refcoor,absval=absval,showblocks=showblocks,relaxes=relaxes,$
	base_state=base_state,orbinfo=orbinfo,trackfile=trackfile,memefficient=memefficient,ptpos=ptpos,ptradius=ptradius,$
	timeunit=timeunit,hideimage=hideimage,useextrema=useextrema,scaleval=scaleval,regrid=regrid,gausswidth=gausswidth,$
	fsuffix=fsuffix,minpartmass=minpartmass

    compile_opt idl2
    fname = filename
	if n_elements(pos) eq 0 then begin
		if keyword_set(hideaxes) then pos = [0.0, 0.0, 1.0, 1.0] else pos = [0.08, 0.07, 0.95, 0.94]
	endif
	if n_elements(simsize) eq 0 then begin
		if n_elements(xr) eq 2 then begin
			simsize = xr[1] - xr[0]
		endif
		if n_elements(yr) eq 2 then begin
			simsize = yr[1] - yr[0]
		endif
		if n_elements(zr) eq 2 then begin
			simsize = zr[1] - zr[0]
		endif
	endif
	if n_elements(slicetype) eq 0 then slicetype = 'plane'
	if n_elements(colorbarcolor) eq 0 then colorbarcolor = 'white'
	if n_elements(imgsizex) eq 0 then imgsizex = 1000
	if n_elements(fprefix) eq 0 then fprefix = ''
	if n_elements(fsuffix) eq 0 then fsuffix = ''
	if n_elements(thrtype) eq 0 and n_elements(thrvar) gt 0 then thrtype=make_array(n_elements(thrvar),/string,value='min')
   	if n_elements(exactmult) eq 0 then exactmult = 1.0
	if n_elements(datatime) ne 0 then time = datatime
	if n_elements(charsize) eq 0 then charsize = 1.0
	if n_elements(annotatepos) eq 0 then annotatepos = 'ur'
	if n_elements(output) eq 0 then output = 'png'
	if n_elements(special) eq 0 then special = ''
	if n_elements(timeunit) eq 0 then timeunit = 's'
	if n_elements(contours) ne 0 then begin
		if n_tags(contours) lt 6 then begin
			print, 'Error: Wrong contours specification.'
			return
		endif
	endif
	if n_elements(minpartmass) eq 0 then minpartmass = 0.0
   	if output ne 'x' then nowindow = 1
	if (n_elements(xr) eq 1 or xr[0] eq xr[1]) then begin
		sliceplane = 'x'
		xrange = [xr[0],xr[0]]
		yrange = yr
		zrange = zr
	endif else if (n_elements(yr) eq 1 or yr[0] eq yr[1]) then begin
		sliceplane = 'y'
		xrange = xr
		yrange = [yr[0],yr[0]]
		zrange = zr
	endif else if (n_elements(zr) eq 1 or zr[0] eq zr[1]) then begin
		sliceplane = 'z'
		xrange = xr
		yrange = yr
		zrange = [zr[0],zr[0]]
	endif else if special eq 'revolve_z' then begin
		sliceplane = 'z'
		xrange = xr
		yrange = yr
		zrange = zr
	endif else if special eq 'column_x' then begin
		sliceplane = 'x'
		xrange = xr
		yrange = yr
		zrange = zr
	endif else if special eq 'column_y' then begin
		sliceplane = 'y'
		xrange = xr
		yrange = yr
		zrange = zr
	endif else if special eq 'column_z' then begin
		sliceplane = 'z'
		xrange = xr
		yrange = yr
		zrange = zr
	endif else begin
		print, 'Either xrange, yrange, or zrange must have a single value!'
		return
	endelse

	;if slicetype eq 'major' or slicetype eq 'minor' then begin
	;    idens = loaddata(fname,'dens',xrange=[0,simsize],yrange=[0,simsize],zrange=[zrange,zrange],sample=2)
	;    idens = reform(idens)
	;    flatslice = flatten_flash_slice(idens,sliceplane)
	;    angle = mom_of_inertia(flatslice,slicetype)
	;    print, angle*180/!pi
    ;     
    ;    boxwidth = min([xr[1]-xr[0],yr[1]-yr[0]])
    ;     
	;    xrange = fltarr(2)
	;    yrange = fltarr(2)
	;    xrange[1] = xr[1] - (1. - abs(sin(angle)))*boxwidth/2.
   	;    xrange[0] = xr[0] + (1. - abs(sin(angle)))*boxwidth/2.
	;    if xrange[0] gt xrange[1] then begin
	;        tmp = xrange[0]
   	;        xrange[0] = xrange[1]
   	;        xrange[1] = tmp
  	;    endif
  	;    yrange[1] = yr[1] - (1. - abs(cos(angle)))*boxwidth/2.
    ;    ;Begin older code
	;	yrange[0] = yr[0] + (1 - abs(cos(angle)))*boxwidth/2
	;	if yrange[0] gt yrange[1] then begin
	;		tmp = yrange[0]
	;		yrange[0] = yrange[1]
	;		yrange[1] = tmp
	;	endif
	;endif else begin
	;	xrange = xr
	;	yrange = yr
	;endelse
	;zrange = zr

	if (n_elements(extdata) eq 0) then begin
		if (n_elements(var) gt 1) then begin
			for i=0,n_elements(var)-1 do begin
				load_flash_var, tmpslice, filename, var[i], xrange, yrange, zrange, dens=dens, temp=temp, $
					velx=velx, vely=vely, velz=velz, gpot=gpot, sample=sample, lwant=lwant, time=time, simsize=simsize, subtractavg=subtractavgi, $
					xcoords=xcoords, ycoords=ycoords, zcoords=zcoords, refcoor=refcoor, special=special, base_state=base_state, orbinfo=orbinfo, $
					trackfile=trackfile, memefficient=memefficient,particles=particles
				if i eq 0 then slice = tmpslice else slice = slice*tmpslice
				if (keyword_set(scaleval)) then slice = slice*scaleval
			endfor
		endif else begin
			load_flash_var, slice, filename, var, xrange, yrange, zrange, dens=dens, temp=temp, $
				velx=velx, vely=vely, velz=velz, gpot=gpot, sample=sample, lwant=lwant, time=time, simsize=simsize, subtractavg=subtractavgi, $
				xcoords=xcoords, ycoords=ycoords, zcoords=zcoords, refcoor=refcoor, special=special, base_state=base_state, orbinfo=orbinfo, $
				trackfile=trackfile, memefficient=memefficient,particles=particles
			if (keyword_set(scaleval)) then slice = slice*scaleval
		endelse
	endif else begin
		slice = extdata
	endelse

	if (n_elements(rangemax) ne 0) then begin
		if (keyword_set(useextrema)) then begin
			maxslice = max(slice)
			rngmax = maxslice
			rngmin = rangemax*maxslice
			print, rngmin, rngmax
		endif else rngmax = rangemax
	endif
	if (n_elements(rangemin) ne 0) then begin
		if (keyword_set(useextrema)) then begin
			minslice = min(slice)
			rngmin = minslice
			rngmax = rangemin*minslice
		endif else rngmin = rangemin
	endif

	if n_elements(my_ct) gt 1 then begin
		custom_ct = intarr(256,3)
		ct_part_len = make_array(n_elements(my_ct)+1, value=0)
		for i=0,n_elements(my_ct)-1 do begin
			reverse_ct = 0
			if my_ct[i] lt 0 then begin
				reverse_ct = 1
			endif
			loadct, abs(my_ct[i]), rgb_table=colors
			;if i eq 0 then colors = reverse(colors,1)
			if keyword_set(ctswitch) then begin
				if i eq 0 then begin	
					if keyword_set(log) then begin
						ct_part_len[i+1] = floor(256.*(alog10(ctswitch)-alog10(rngmin))/(alog10(rngmax) - alog10(rngmin)))	
					endif else begin
						ct_part_len[i+1] = floor(256.*(ctswitch - rngmin)/(rngmax - rngmin))	
					endelse
				endif else begin
					if keyword_set(log) then begin
						ct_part_len[i+1] = floor(256.*(alog10(rngmax)-alog10(ctswitch))/(alog10(rngmax) - alog10(rngmin)))	
					endif else begin
						ct_part_len[i+1] = floor(256.*(rngmax - ctswitch)/(rngmax - rngmin))	
					endelse
				endelse
			endif else begin
				ct_part_len[i+1] = floor(256./n_elements(my_ct))
			endelse
			if reverse_ct eq 1 then begin
				custom_ct[ct_part_len[i]:ct_part_len[i]+ct_part_len[i+1]-1,*] = reverse(colors[floor(dindgen(ct_part_len[i+1])*256./ct_part_len[i+1]),*])
			endif else begin
				custom_ct[ct_part_len[i]:ct_part_len[i]+ct_part_len[i+1]-1,*] = colors[floor(dindgen(ct_part_len[i+1])*256./ct_part_len[i+1]),*]
			endelse
		endfor
	endif

	if (special eq 'revolve_z') then begin
		tmpslice = reform(slice)
		npts = n_elements(xcoords)
		r = dindgen(floor(npts/2))
		zpos = cmreplicate(dindgen(npts), floor(npts/2))
		theta = dindgen(npts)/(npts+1)*2*!pi
		slice = dblarr(npts,floor(npts/2),1)
		for t=0,n_elements(theta)-1 do begin
			xpos = transpose(cmreplicate(r*cos(theta[t]),npts)) + npts/2.0
			ypos = transpose(cmreplicate(r*sin(theta[t]),npts)) + npts/2.0
			slice = slice + interpolate(tmpslice,xpos,ypos,zpos)
		endfor
		slice = transpose(slice)/npts
	endif

	if keyword_set(fieldvarx) then begin
		load_flash_var, fieldslicex, filename, fieldvarx, xrange, yrange, zrange, dens=dens, temp=temp, $
			velx=velx, vely=vely, velz=velz, gpot=gpot, sample=sample, lwant=lwant, time=time, simsize=simsize, subtractavg=subtractavgi, refcoor=refcoor,particles=particles
		load_flash_var, fieldslicey, filename, fieldvary, xrange, yrange, zrange, dens=dens, temp=temp, $
			velx=velx, vely=vely, velz=velz, gpot=gpot, sample=sample, lwant=lwant, time=time, simsize=simsize, subtractavg=subtractavgi, refcoor=refcoor,particles=particles
	endif

	slice = reform(slice)
	if keyword_set(negative) then slice = (-1.0)*temporary(slice)
	if keyword_set(absval) then slice = abs(temporary(slice))

	slice_dims = size(slice, /dimensions)
	if(size(slice_dims,/n_elements) eq 2) then begin
		slice_dims = [slice_dims, 1]
	endif

	if keyword_set(thrvar) then begin
		thrslice = dblarr(slice_dims[0],slice_dims[1],slice_dims[2],n_elements(thrvar))
	endif

	;below may not work for diagonal slices
	for i=0,n_elements(thrvar)-1 do begin
		load_flash_var, newthrslice, filename, thrvar[i], xrange, yrange, zrange, dens=dens, temp=temp, $
			velx=velx, vely=vely, velz=velz, gpot=gpot, sample=sample, lwant=lwant, simsize=simsize, $
			refcoor=refcoor, orbinfo=orbinfo, xcoords=xcoords, ycoords=ycoords, zcoords=zcoords, $
			particles=particles
		if (n_elements(scaleval) ne 0) then newthrslice = newthrslice*scaleval
		thrslice[*,*,*,i] = newthrslice
	endfor

	if keyword_set(contours) then begin
		;if contours.var eq var then contourslice = slice else begin
		;	for i=0,n_elements(thrvar)-1 do if thrvar[i] eq contours.var then contourslice = thrslice[*,*,*,i]
		;endelse
		;if n_elements(contourslice) eq 0 then begin
			load_flash_var, contourslice, filename, contours.var, xrange, yrange, zrange, dens=dens, temp=temp, $
				velx=velx, vely=vely, velz=velz, gpot=gpot, sample=sample, lwant=lwant, simsize=simsize, refcoor=refcoor, $
				orbinfo=orbinfo,particles=particles
			if (n_elements(scaleval) ne 0) then contourslice = contourslice*scaleval
			contourslice = reform(contourslice)
		;endif
	endif
		
	if keyword_set(excision) then begin
		exx = excision[0]
		exy = excision[1]
		exr = excision[2]
		for i=0,slice_dims[0]-1 do begin
		for j=0,slice_dims[1]-1 do begin
		for k=0,slice_dims[2]-1 do begin
			d = sqrt((xcoords[i]-exx)^2. + (ycoords[j]-exy)^2.)
			if (d le exr) then slice[i,j,k] = rngmin
		endfor
		endfor
		endfor
	endif

	if (~keyword_set(rangemin) or ~keyword_set(rangemax)) and keyword_set(thrvar) then begin
		indices = indgen(n_elements(slice), /long)
		for i=0,n_elements(thrvar)-1 do begin
			if thrtype[i] eq 'max' then begin
				newindices = where(reform(thrslice[*,*,*,i]) le thrval[i], /null, /l64)
			endif else begin
				newindices = where(reform(thrslice[*,*,*,i]) ge thrval[i], /null, /l64)
			endelse
			indices = cmset_op(indices, 'AND', newindices)
		endfor
	endif

	if keyword_set(rangemin) or (keyword_set(rangemax) and keyword_set(useextrema)) then begin
		min_val = rngmin
	endif else begin
		if keyword_set(thrvar) then begin
			min_val = min(slice[indices])
		endif else begin
			min_val = min(slice)
		endelse
	endelse

	if keyword_set(rangemax) then begin
		max_val = rngmax
	endif else begin
		if keyword_set(thrvar) then begin
			max_val = max(slice[indices])
		endif else begin
			max_val = max(slice)
		endelse
	endelse

	if n_elements(contourslice) ne 0 then begin
		min_cont_val = contours.min
		max_cont_val = contours.max
	endif
	
	if n_elements(symrange) eq 0 then begin
		if n_elements(rangemin) ne 0 then min_val = rngmin
		if n_elements(rangemax) ne 0 then max_val = rngmax
	endif else begin
		if n_elements(rangemin) ne 0 then begin
			min_val = rngmin
			max_val = -rngmin
		endif else if n_elements(rangemax) ne 0 then begin
			min_val = -rngmax
			max_val = rngmax
		endif
		if n_elements(rangemin) eq 0 and n_elements(rangemax) eq 0 then begin
			min_val = -max([abs(min_val),abs(max_val)])
			max_val = -min_val
		endif
	endelse

	if n_elements(thrvar) gt 0 then begin
		;for i=0,slice_dims[0]-1 do begin
		;	for j=0,slice_dims[1]-1 do begin
		;		for k=0,slice_dims[2]-1 do begin
					for l=0,n_elements(thrvar)-1 do begin
						if thrtype[l] eq 'max' then begin
							slice[where(reform(thrslice[*,*,*,l]) gt thrval[l], /null, /l64)] = 0.d0
							;if thrslice[i,j,k,l] gt thrval[l] then begin
							;	if keyword_set(log) then begin
							;		slice[i,j,k] = 1.d-100
							;	endif else begin
							;		slice[i,j,k] = 0.d0
							;	endelse
							;	;if n_elements(rangemin) eq 0 then begin
							;	;	slice[i,j,k] = min_val
							;	;endif else begin
							;	;	slice[i,j,k] = rngmin
							;	;endelse
							;	if keyword_set(contours) ne 0 then contourslice[i,j,k] = min_cont_val
							;	if n_elements(fieldvarx) ne 0 then begin
							;		fieldslicex[i,j,k] = min_fieldx_val
							;		fieldslicey[i,j,k] = min_fieldy_val
							;	endif
							;endif
						endif else begin
							slice[where(reform(thrslice[*,*,*,l]) lt thrval[l], /null, /l64)] = 0.d0
							;if thrslice[i,j,k,l] lt thrval[l] then begin
							;	if keyword_set(log) then begin
							;		slice[i,j,k] = 1.d-100
							;	endif else begin
							;		slice[i,j,k] = 0.d0
							;	endelse
							;	;if n_elements(rangemin) eq 0 then begin
							;	;	slice[i,j,k] = min_val
							;	;endif else begin
							;	;	slice[i,j,k] = rngmin
							;	;endelse
							;	if keyword_set(contours) ne 0 then contourslice[i,j,k] = min_cont_val
							;endif
						endelse
					endfor
		;		endfor
		;	endfor
		;endfor
	endif

	if slicetype eq 'major' or slicetype eq 'minor' then begin
		newslice = fltarr(slice_dims[2],slice_dims[2])
		for i=0,(slice_dims[2]-1) do begin
			if cos(angle) ge 0 then begin
				x = min([round(cos(angle)*double(i)),slice_dims[0]-1])
			endif else begin
				x = max([round(double(slice_dims[0]-1) + cos(angle)*double(i)),0])
			endelse
			if sin(angle) ge 0 then begin
				y = min([round(sin(angle)*double(i)),slice_dims[1]-1])
			endif else begin
				y = max([round(double(slice_dims[1]-1) + sin(angle)*double(i)),0])
			endelse
			newslice[i,*] = slice[x,y,*]
		endfor
		slice = reform(newslice)
	endif

	;if (special eq 'column_x' or special eq 'column_y' or special eq 'column_z') then begin
	;	if n_elements(rangemin) ne 0 then begin
	;		indices = where(slice eq rngmin, count, /l64)
	;		if count ne 0 then slice[indices] = 0.e0
	;	endif
	;endif

	if (special eq 'column_x') then begin
		slice = total(slice, 1)*(xcoords[1]-xcoords[0])
		if n_elements(rangemin) ne 0 then begin
			indices = where(slice eq 0.e0, count, /l64)
			if count ne 0 then slice[indices] = rngmin
		endif
		if keyword_set(contours) then begin
			indices = where(contourslice eq min_cont_val, count, /l64)
			if count ne 0 then contourslice[indices] = 0.e0
			contourslice = total(contourslice, 1)*(xcoords[1]-xcoords[0])
			indices = where(contourslice eq 0.e0, count, /l64)
			if count ne 0 then contourslice[indices] = min_cont_val
		endif
	endif
	if (special eq 'column_y') then begin
		slice = total(slice, 2)*(ycoords[1]-ycoords[0])
		if n_elements(rangemin) ne 0 then begin
			indices = where(slice eq 0.e0, count, /l64)
			if count ne 0 then slice[indices] = rngmin
		endif
		if keyword_set(contours) then begin
			indices = where(contourslice eq min_cont_val, count, /l64)
			if count ne 0 then contourslice[indices] = 0.e0
			contourslice = total(contourslice, 2)*(ycoords[1]-ycoords[0])
			indices = where(contourslice eq 0.e0, count, /l64)
			if count ne 0 then contourslice[indices] = min_cont_val
		endif
	endif
	if (special eq 'column_z') then begin
		slice = total(slice, 3)*(zcoords[1]-zcoords[0])
		slice_dims = size(slice, /dimensions)
		slice_dims = [slice_dims, 1]
		;if n_elements(rangemin) ne 0 then begin
		;	indices = where(slice eq 0.e0, count, /l64)
		;	if count ne 0 then slice[indices] = rngmin
		;endif
		if keyword_set(contours) then begin
			indices = where(contourslice eq min_cont_val, count, /l64)
			if count ne 0 then contourslice[indices] = 0.e0
			contourslice = total(contourslice, 3)*(zcoords[1]-zcoords[0])
			indices = where(contourslice eq 0.e0, count, /l64)
			if count ne 0 then contourslice[indices] = min_cont_val
		endif
	endif

	if keyword_set(gausswidth) then begin
		kerw = gausswidth / (xrange[1] - xrange[0]) * double(slice_dims[0])
		slice = filter_image(slice,fwhm_gaussian=kerw,/all)
	endif

	if keyword_set(regrid) then begin
		if (sliceplane eq 'x') then begin
			ncx = 1
			ncy = round((yrange[1] - yrange[0])/regrid)
			ncz = round((zrange[1] - zrange[0])/regrid)
		endif
		if (sliceplane eq 'y') then begin
			ncx = round((xrange[1] - xrange[0])/regrid)
			ncy = 1
			ncz = round((zrange[1] - zrange[0])/regrid)
		endif
		if (sliceplane eq 'z') then begin
			ncx = round((xrange[1] - xrange[0])/regrid)
			ncy = round((yrange[1] - yrange[0])/regrid)
			ncz = 1
		endif
		print, slice_dims
		slice = congrid(reform(slice), ncx, ncy, ncz, /center, cubic=-0.5, /interp)		
		slice_dims = size(slice, /dimensions)
		print, slice_dims
	endif

	if n_elements(rangemin) ne 0 then begin
		minindices = where(slice lt rngmin, count)
		if count ne 0 then slice[minindices] = rngmin
	endif
	if n_elements(rangemax) ne 0 then begin
		if n_elements(lmin) ne 0 then begin
			maxindices = where(slice gt rngmax, count)
			if count ne 0 then slice[maxindices] = rngmax
			minindices = where(-slice gt rngmax, count)
			if count ne 0 then slice[minindices] = rngmax
		endif else begin
			maxindices = where(slice gt rngmax, count)
			if count ne 0 then slice[maxindices] = rngmax
		endelse
	endif

	if keyword_set(log) then begin
		if n_elements(lmin) eq 0 then begin
			min_pos = min(slice[where(slice gt 0.0, /L64)])
			tmp_indx = where(slice le 0.0)
			if tmp_indx[0] ne -1 then slice[tmp_indx] = min_pos
			slice = alog10(slice)
		endif else begin
			tmp_slice = dblarr(slice_dims)
			tmp_ind = where(slice gt 0.0 and abs(slice) ge 10.0^lmin)
			tmp_slice[tmp_ind] = alog10(slice[tmp_ind])
			if min(tmp_slice[tmp_ind]) lt 0 then tmp_slice[tmp_ind] = tmp_slice[tmp_ind] - min(tmp_slice[tmp_ind])
			tmp_ind = where(slice lt 0.0 and abs(slice) ge 10.0^lmin)
			tmp_slice[tmp_ind] = -alog10(abs(slice[tmp_ind]))
			if max(tmp_slice[tmp_ind]) gt 0 then tmp_slice[tmp_ind] = tmp_slice[tmp_ind] - max(tmp_slice[tmp_ind])
			slice = tmp_slice
		endelse
	endif

	if n_elements(rangemin) ne 0 then begin
		plot_min = rngmin
		if keyword_set(log) then plot_min = alog10(plot_min)	
	endif else plot_min = min(slice)
	if n_elements(rangemax) ne 0 then begin
		plot_max = rngmax
		if keyword_set(log) then plot_max = alog10(plot_max)	
	endif else plot_max = max(slice)

	thisDevice = !D.NAME
	set_plot, 'z'
	plotcolor = fsc_color('white',/nodisplay)
	!P.background = fsc_color('black',/nodisplay)

	if (special eq 'column_x') then begin
		if keyword_set(exactsize) then imgsizex = floor(double(slice_dims[1])/(pos[2]-pos[0]))*exactmult
		if ~keyword_set(imgsizey) then imgsizey = floor(double(slice_dims[2])/double(slice_dims[1])*imgsizex)
	endif else if (special eq 'column_y') then begin
		if keyword_set(exactsize) then imgsizex = floor(double(slice_dims[0])/(pos[2]-pos[0]))*exactmult
		if ~keyword_set(imgsizey) then imgsizey = floor(double(slice_dims[2])/double(slice_dims[0])*imgsizex)
	endif else begin
		if keyword_set(exactsize) then imgsizex = floor(double(slice_dims[0])/(pos[2]-pos[0]))*exactmult
		if ~keyword_set(imgsizey) then imgsizey = floor(double(slice_dims[1])/double(slice_dims[0])*imgsizex)
	endelse

	if output ne 'x' then begin
		logf = ''
		if keyword_set(log) then logf = '_log'
		varname = ''
		for i=0,n_elements(var)-1 do begin
			varname = varname + var[i]	
		endfor
		output_name = fprefix+varname+'_'+slicetype+'_'+sliceplane+logf+fsuffix+'_'+filename+'.ps'
		jps_start, filename=output_name,xsize=imgsizex/100.,ysize=imgsizey/100.
	endif else begin
		set_plot, thisDevice
	endelse
	if n_elements(custom_ct) ne 0 then tvlct, custom_ct else begin
		loadct, abs(my_ct)
		if my_ct lt 0 then begin
			tvlct, r, g, b, /get
			tvlct, reverse(r), reverse(g), reverse(b)
		endif
	endelse
   	tvlct, red, green, blue, /get
	if n_elements(lmin) gt 0 then begin
		slice_dims = size(slice, /dimensions)
		thisImage = bytarr(slice_dims)
		thisImage[where(slice gt 0)] = 128 + floor(bytscl(slice[where(slice gt 0)],min=lmin,max=max(slice),top=!D.N_Colors-1)/2)
		thisImage[where(slice lt 0)] = bytscl(slice[where(slice lt 0)],min=min(slice),max=-lmin,top=!D.N_Colors/2)
		thisImage[where(slice eq 0)] = 128
	endif else begin
		print, min(slice), max(slice), plot_min, plot_max
		thisImage = bytscl(slice,max=plot_max,min=plot_min,top=!D.N_Colors-1)
	endelse

	s = SIZE(thisImage)
	image3d = BYTARR(3, s[1], s[2])
	image3d[0, *, *] = red[thisImage]
	image3d[1, *, *] = green[thisImage]
	image3d[2, *, *] = blue[thisImage]
	;loadct, 1
   	;tvlct, red, green, blue, /get
	;image3d2 = BYTARR(3, s[1], s[2])
	;image3d2[0, *, *] = red[thisImage2]
	;image3d2[1, *, *] = green[thisImage2]
	;image3d2[2, *, *] = blue[thisImage2]
	if ~keyword_set(nowindow) then begin
		window,1,xsize=1000,ysize=1000
	endif
	erase, color=fsc_color('black',/nodisplay)
	polyfill, [1.0,1.0,0.0,0.0,1.0], [1.0,0.0,0.0,1.0,1.0], /normal, color=fsc_color('black',/nodisplay)

	case sliceplane of
		'x': begin
			;pos[0] = pos[0] + (pos[2] - pos[0])*(ycoords[0] - yrange[0])/(ycoords[n_elements(ycoords)-1] - ycoords[0])
			;pos[1] = pos[1] + (pos[3] - pos[1])*(zcoords[0] - zrange[0])/(zcoords[n_elements(zcoords)-1] - zcoords[0])
			pos[2] = pos[0] + (pos[2] - pos[0])*(ycoords[n_elements(ycoords)-1] - yrange[0])/(ycoords[n_elements(ycoords)-1] - ycoords[0])
			pos[3] = pos[1] + (pos[3] - pos[1])*(zcoords[n_elements(zcoords)-1] - zrange[0])/(zcoords[n_elements(zcoords)-1] - zcoords[0])
		end
		'y': begin
			;pos[0] = pos[0] + (pos[2] - pos[0])*(xcoords[0] - xrange[0])/(xcoords[n_elements(xcoords)-1] - xcoords[0])
			;pos[1] = pos[1] + (pos[3] - pos[1])*(zcoords[0] - zrange[0])/(zcoords[n_elements(zcoords)-1] - zcoords[0])
			pos[2] = pos[0] + (pos[2] - pos[0])*(xcoords[n_elements(xcoords)-1] - xrange[0])/(xcoords[n_elements(xcoords)-1] - xcoords[0])
			pos[3] = pos[1] + (pos[3] - pos[1])*(zcoords[n_elements(zcoords)-1] - zrange[0])/(zcoords[n_elements(zcoords)-1] - zcoords[0])
		end
		'z': begin
			;pos[0] = pos[0] + (pos[2] - pos[0])*(xcoords[0] - xrange[0])/(xcoords[n_elements(xcoords)-1] - xcoords[0])
			;pos[1] = pos[1] + (pos[3] - pos[1])*(ycoords[0] - yrange[0])/(ycoords[n_elements(ycoords)-1] - ycoords[0])
			pos[2] = pos[0] + (pos[2] - pos[0])*(xcoords[n_elements(xcoords)-1] - xrange[0])/(xcoords[n_elements(xcoords)-1] - xcoords[0])
			pos[3] = pos[1] + (pos[3] - pos[1])*(ycoords[n_elements(ycoords)-1] - yrange[0])/(ycoords[n_elements(ycoords)-1] - ycoords[0])
		end
	endcase

	!p.position = pos
	if ~keyword_set(hideimage) then begin
		if slicetype eq 'minor' or slicetype eq 'major' then begin
			tvimage, image3d, position=pos,/nointerpolation
		endif else begin
			tvimage, image3d, image3d2, position=pos;,/nointerpolation
			;blendimage, image3d, image3d2, alpha=0.0, /keep_aspect, position=pos,/nointerpolation
		endelse
	endif

	!p.color = fsc_color("white",/nodisplay)
	if ~keyword_set(hideaxes) then begin
		case sliceplane of
			'x': begin
				;range = max([yrange[1]-yrange[0],zrange[1]-zrange[0]])
				;midy = yrange[0]+(yrange[1]-yrange[0])/2
				;midz = zrange[0]+(zrange[1]-zrange[0])/2
				;contour, tri_surf(slice), levels=(1-(reverse(dindgen(20)+1)/21.)^2.)*(max_cont_val-min_cont_val)+min_val, /follow, /noerase, xstyle=1+4, ystyle=1+4, color=0
				;mylevels = (dindgen(10)+1)/10.
				;mythick = make_array(10, value=1.)
				;mythick[9] = 3.
				;contour, contourslice, levels=mylevels, c_thick=mythick, /follow, /noerase, xstyle=1+4, ystyle=1+4, color=0, position=pos
				;yvec = reverse(dindgen(slice_dims[0]))/double(slice_dims[0])*range + midy-range/2
				;yvals = dblarr(slice_dims[0],slice_dims[1])
				;END OLDER CODE
				;for i=0,slice_dims[0]-1 do zvals[i,*] = zvec
																																												   
				;contour, slice, yvals, zvals, levels=(1-(reverse(dindgen(20)+1)/21.)^2.)*(max_val-min_val)+min_val,$

				;   xrange=[midy-range/2, midy+range/2], yrange=[midz-range/2, midz+range/2], /noerase, xstyle=1+4, ystyle=1+4, /follow, color=0, /irregular

				plot, yrange, zrange, /noerase, /nodata, /xsty, /ysty, charsize=charsize*imgsizex/1000., $
					color=plotcolor, position=pos, xticks=xticks, yticks=yticks

				if special eq 'velyzrot' then begin
					load_flash_var, dens, filename, 'dens', xrange, yrange, zrange, dens=dens, temp=temp, $
						velx=velx, vely=vely, velz=velz, gpot=gpot, log=log, sample=sample, lwant=lwant, time=time, simsize=simsize, refcoor=refcoor,particles=particles
					load_flash_var, vely, filename, 'vely', xrange, yrange, zrange, dens=dens, temp=temp, $
						velx=velx, vely=vely, velz=velz, gpot=gpot, log=log, sample=sample, lwant=lwant, time=time, simsize=simsize, $
						zcoords=zcoords, ycoords=ycoords, refcoor=refcoor,particles=particles
					load_flash_var, velz, filename, 'velz', xrange, yrange, zrange, dens=dens, temp=temp, $
						velx=velx, vely=vely, velz=velz, gpot=gpot, log=log, sample=sample, lwant=lwant, time=time, simsize=simsize, refcoor=refcoor,particles=particles
					vely = reform(vely)
					vely = vely - avg(vely[where(temp lt 1e5)])
					velz = reform(velz)
					velz = velz - avg(velz[where(temp lt 1e5)])
					dims = size(vely, /dimensions)
					ycoordarr = dblarr(dims)
					for j=0,dims[0]-1 do ycoordarr[*,j] = ycoords
					zcoordarr = dblarr(dims)
					for j=0,dims[1]-1 do zcoordarr[j,*] = reverse(zcoords)
					;ycom = total(dens*ycoordarr)/total(dens)
					;zcom = total(dens*zcoordarr)/total(dens)
					;ycom = 1.455e10
					;zcom = 1.5e10
					;print, ycom, zcom
					csize = ycoords[1]-ycoords[0]
					vely[where(temp ge 1e6)] = 0.0
					velz[where(temp ge 1e6)] = 0.0
					;ycoordarr = ycoordarr-ycom
					;zcoordarr = zcoordarr-zcom
					;distmat = sqrt(ycoordarr^2. + zcoordarr^2)
					;
					;velmag = (vely*zcoordarr - velz*ycoordarr)/distmat
					;vely = velmag*zcoordarr/distmat
					;velz = velmag*ycoordarr/distmat
					velovect, vely, velz, ycoords, zcoords, /noerase, position=pos, xstyle=5, ystyle=5, $
						xrange=[ycoords[0]-0.5*csize,ycoords[dims[0]-1]+0.5*csize], yrange=[zcoords[0]-0.5*csize,zcoords[dims[1]-1]+0.5*csize]
				endif
			end
			'y': begin
				plot, xrange, zrange, /noerase, /nodata, /xsty, /ysty, charsize=charsize*imgsizex/1000., $
					color=plotcolor, position=pos, xticks=xticks, yticks=yticks

			end
			'z': begin
				if slicetype eq 'major' or slicetype eq 'minor' then begin
					if xr[1]-xr[0] lt yr[1]-yr[0] then xyrange = xr else xyrange = yr
					plot, xyrange, zrange, /noerase, /nodata, /xsty, /ysty, color=plotcolor, $
						position=pos, charsize=charsize*imgsizex/1000., xticks=xticks, yticks=yticks

				;if slicetype eq 'major' or slicetype eq 'minor' then begin
				;   xyrange = [0,sqrt((xrange[1]-xrange[0])^2+(yrange[1]-yrange[0])^2)]
				;   range = max([xyrange[1]-xyrange[0],zrange[1]-zrange[0]])
				;   midxy = xyrange[0]+(xyrange[1]-xyrange[0])/2
				;   midz = zrange[0]+(zrange[1]-zrange[0])/2
				;   plot, [midxy-range/2, midxy+range/2], [midz-range/2, midz+range/2], /noerase, /nodata, /xsty, /ysty, charsize=1.0

				endif else begin
					if keyword_set(fieldvarx) then begin
						csize = xcoords[1]-xcoords[0]
						xcoordmat = dblarr(slice_dims[0],slice_dims[1])
						ycoordmat = dblarr(slice_dims[0],slice_dims[1])
						for i=0,slice_dims[1]-1 do begin
							xcoordmat[*,i] = xcoords
							ycoordmat[i,*] = ycoords
						endfor
						partvelvec, reform(fieldslicex), reform(fieldslicey), xcoordmat, ycoordmat, position=pos, xstyle=5, ystyle=5, $
							xrange=[xcoords[0]-0.5*csize,xcoords[slice_dims[0]-1]+0.5*csize], yrange=[ycoords[0]-0.5*csize,ycoords[slice_dims[1]-1]+0.5*csize], /noerase, $
							fraction=0.001, maxmag=fieldmax, seed=10, length=0.02
					endif

					if (keyword_set(relaxes)) then begin
						plot, [relaxes[0,0],relaxes[1,0]], [relaxes[0,1],relaxes[1,1]], /noerase, /nodata, /xsty, /ysty, charsize=charsize*imgsizex/1000., $
							color=plotcolor, position=pos, xticks=xticks, yticks=yticks
					endif else begin
						plot, xrange, yrange, /noerase, /nodata, /xsty, /ysty, charsize=charsize*imgsizex/1000., $
							color=plotcolor, position=pos, xticks=xticks, yticks=yticks
					endelse
				end
			end
		endcase
	endif

	if (keyword_set(contours)) then begin
		if contours.log eq 1 then begin
			contourlevels = 1.e1^((1-(reverse(dindgen(contours.num)+1)/(double(contours.num)+1.)))*(alog10(max_cont_val)-alog10(min_cont_val))+alog10(min_cont_val))
		endif else begin
			contourlevels = (1-(reverse(dindgen(contours.num)+1)/(double(contours.num)+1.)))*(max_cont_val-min_cont_val)+min_cont_val
		endelse

		if n_elements(contourslice) gt 0 then begin
			loadct, contours.ct, ncolors=contours.num, bottom=3, rgb_table = colortable
			tvlct, colortable[*,0], colortable[*,1], colortable[*,2]
			colorindices = indgen(contours.num)
			if contours.fill eq 1 then begin
				contour, contourslice, levels=contourlevels, $
					/noerase, xstyle=1+4, ystyle=1+4, position=pos, c_colors=colorindices, /fill, /closed
				contour, contourslice, levels=contourlevels, $
					/noerase, xstyle=1+4, ystyle=1+4, position=pos, color=cgColor('black'), /closed
			endif else begin
				contour, contourslice, levels=contourlevels, $
					/noerase, xstyle=1+4, ystyle=1+4, position=pos, c_colors=colorindices, /closed
			endelse
		endif
	end

	case sliceplane of
		'x': begin
			slice_dir = 2
		end
		'y': begin
			slice_dir = 1
		end
		'z': begin
			slice_dir = 0
		end
	endcase

	if keyword_set(showblocks) then begin
		jread_amr, fname, var_name=var, tree=tree, parameters=params
		draw3d_blocks, fsc_color('white'), orientation=0, slice_dir=slice_dir, $
			parameters=params, tree=tree, xrange=xrange, yrange=yrange, zrange=zrange
	endif

	if special eq 'markcenter' then begin
		mysym = findgen(49) * (!pi*2/48.)  
		usersym, cos(mysym), sin(mysym), thick=3.0
		plots, 0.0d0, 0.0d0, psym=8, symsize=min([xrange[1]-xrange[0],yrange[1]-yrange[0]])*0.05, color=fsc_color('blue'), /data
	endif

	num_particles = 0
	if keyword_set(ptpos) then begin
		num_particles = 1
	endif

	if keyword_set(particles) then begin
		if tag_exist(particles, 'posx') then begin
			num_particles = n_elements(particles)
		endif
	endif

	for i=0,num_particles-1 do begin
		if (particles[i].mass lt minpartmass) then continue
		if keyword_set(ptpos) then begin
			point = ptpos
		endif else if num_particles gt 0 then begin
			point = [particles[i].posx, particles[i].posy, particles[i].posz]
		endif
		mysym = findgen(49) * (!pi*2/48.)  
		usersym, cos(mysym), sin(mysym), thick=3.0, /fill
		if (sliceplane eq 'x') then begin
			sympos = [pos[0] + (pos[2] - pos[0])*(point[1] - yrange[0])/(yrange[1] - yrange[0]), $
					  pos[1] + (pos[3] - pos[1])*(point[2] - zrange[0])/(zrange[1] - zrange[0])]
		endif else if (sliceplane eq 'y') then begin
			sympos = [pos[0] + (pos[2] - pos[0])*(point[0] - xrange[0])/(xrange[1] - xrange[0]), $
					  pos[1] + (pos[3] - pos[1])*(point[2] - zrange[0])/(zrange[1] - zrange[0])]
		endif else if (sliceplane eq 'z') then begin
			sympos = [pos[0] + (pos[2] - pos[0])*(point[0] - xrange[0])/(xrange[1] - xrange[0]), $
					  pos[1] + (pos[3] - pos[1])*(point[1] - yrange[0])/(yrange[1] - yrange[0])]
		endif
		th = 2.0*!dpi/100.0*dindgen(101) 
		if n_elements(ptradius) EQ 0 then begin
			r = 0.02 
		endif else begin
			r = ptradius/(xrange[1]-xrange[0])
		endelse
		ptx = r*cos(th)
		pty = r*sin(th)*(xrange[1]-xrange[0])/(yrange[1]-yrange[0])

		polyfill, sympos[0]+ptx, sympos[1]+pty, color=fsc_color('blue'), /normal, _EXTRA=extra
	endfor

    if ~keyword_set(hidetime) then begin
		if annotatepos eq 'flip' then begin
			xyouts, 0.8, 0.95, 'Time: ' + strcompress(string(time, format='(G-15.7)'), /remove_all) + 's', $
				charsize=2.0*charsize*imgsizey/1000., color=plotcolor, /normal
		endif else begin
			if (timeunit eq 's') then begin
				xyouts, 0.08, 0.95, 'Time: ' + strcompress(string(time, format='(G-15.7)'), /remove_all) + ' seconds', $
					charsize=2.0*charsize*imgsizey/1000., color=plotcolor, /normal
			endif else if (timeunit eq 'h') then begin
				xyouts, 0.08, 0.95, 'Time: ' + strcompress(string(time/3600., format='(G-15.7)'), /remove_all) + ' hours', $
					charsize=2.0*charsize*imgsizey/1000., color=plotcolor, /normal
			endif else if (timeunit eq 'd') then begin
				xyouts, 0.08, 0.95, 'Time: ' + strcompress(string(time/86400., format='(G-15.7)'), /remove_all) + ' days', $
					charsize=2.0*charsize*imgsizey/1000., color=plotcolor, /normal
			endif
		endelse
	endif
		
	if n_elements(custom_ct) ne 0 then tvlct, custom_ct else tvlct, red, green, blue
	if min_val ne max_val and ~keyword_set(hidecolorbar) then begin
		if annotatepos eq 'flip' then begin
			colorbar_pos=[0.18, 0.75, 0.21, 0.91]
		endif else begin
			colorbar_pos=[0.89, 0.75, 0.92, 0.91]
		endelse
		print, colorbar_pos
		fsc_colorbar, /vertical, minrange=plot_min, $
			maxrange=plot_max, position=colorbar_pos, /nodisplay,$
			annotatecolor=colorbarcolor, format='(G10.3)', charsize=charsize*min([imgsizex,imgsizey])/1000.
	end

    if output eq 'png' then begin
		jps_end, /png 
		file_delete, output_name
	endif else if output eq 'ps' then jps_end
    fsc_undefine, slice, dens, temp, velx, vely, velz, gpot, time
end
