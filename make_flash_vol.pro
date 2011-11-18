; make_flash_vol.pro
; Generates a volume plot from a FLASH data file.
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
;
; **Optional**
; simsize       (dbl)           - Width of the box in cgs units (assumes box is cubical, used for certain variables).
; slicetype     (str)           - UNUSED.
; rangemin      (dbl)           - Minimum value to plot.
; rangemax      (dbl)           - Maximum value to plot.
; imgsizex      (dbl)           - Size of x dimension of output image, in pixels (applies to PNG only).
; imgsizey      (dbl)           - Size of y dimension of output image, in pixels (applies to PNG only).
; fprefix       (str)           - Prefix to append to output image files, e.g. 'myprefix_restoffilename'
; exactsize     (bool)          - Smallest grid cells will span exactly 1 pixel in output image.
; exactmult     (dbl)           - Exactsize must be true; smallest grid cells will span 1 pixel * exactmult.
; xticks        (int)           - Number of ticks to draw on x axis
; yticks        (int)           - Number of ticks to draw on y axis.
; zticks        (int)           - Number of ticks to draw on z axis (3D only).
; log           (bool)          - Scale data logarithmically.
; colorbarcolor (str)           - Color of colorbar annotation. 
; contourvar    (str)           - Variable used to draw contours on plot.
; thrvar        ([str,str,...]) - Variables used as a "threshold," data in regions where thresholds aren't
;                                 satisfied are set to rangemin.
; thrval        ([dbl,dbl,...]) - Threshold values, see above.
; thrtype       ([str,str,...]) - Whether the threshold is a minimum or a maximum, valid values 'min' or 'max'.
; sample        (int)           - Number of times to downsample data from highest refinment, must be <= 3.
; stride        (int)           - Stride of frame iteration, default 1.
; charsize      (dbl)           - Multiplier of default char size, default 1.0.
; symrange      (bool)          - Makes the data range symmetric depending on min and max values of data.
; lmin          (dbl)           - Used with symrange only, specifies min log value for
;                                 symmetric ranges (Since log can't have negative argument).
; annotatepos   (str)           - Annotation position (colorbar, timestamp). Currently can only set to 'flip'.
; output        (str)           - Specifies output file type, can be 'ps' or 'png' (default).
; ax            (dbl)           - Rotation about the x-axis (Volume plot only).
; az            (dbl)           - Rotation about the z-axis (Volume plot only).
; revolvestep   (dbl)           - Generate a series of frames for each data file, revolving this many degrees along
;								  the axis specified by revolvetype.
; revolvetype   (str)           - The axis to revolve around.
; special       (str)           - Flag used to indicate "special" plots to generate.
; cellsize      (dbl)           - Size of the smallest grid cells in CGS units, used for some variables.
; boxaxes       ([dbl,dbl])     - Draw data in a bounding box of this size (CGS units), volume plot only.
; memefficient  (bool)          - Deallocate variables as soon as they are not needed. Incompatible with revolve plots.
; hideaxes      (bool)          - Do not show axes.
; zoom          (dbl)           - Zoom plot by this factor, volume plot only.
; boxscale      (dbl)           - Scale each axis by this factor, volume plot only.
; indexlength   (int)           - Number of digits in filename index, e.g. _0000 is 4, _00000 is 5, etc.

pro make_flash_vol,filename,var,my_ct,xrange,yrange,zrange,ax=ax,az=az,hidetime=hidetime,charsize=charsize, $
	rangemin=rangemin,rangemax=rangemax,log=log,colorbarcolor=colorbarcolor,thrvar=thrvar,thrval=thrval,$
	sample=sample,output=output,fprefix=fprefix,fsuffix=fsuffix,revolvestep=revolvestep,revolvetype=revolvetype,$
	imgsizex=imgsizex,imgsizey=imgsizey,thrtype=thrtype,special=special,simsize=simsize,cellsize=cellsize,$
	boxaxes=boxaxes,memefficient=memefficient,hideaxes=hideaxes,zoom=zoom,boxscale=boxscale,$
	xticks=xticks,yticks=yticks,zticks=zticks,custom_rot=custom_rot,hidecolorbar=hidecolorbar,mirror=mirror,$
	ctswitch=ctswitch,hideticklabels=hideticklabels,extdata=extdata,relaxes=relaxes,oversample=oversample

    compile_opt idl2
	if n_elements(ax) eq 0 then ax = 30
	if n_elements(az) eq 0 then az = 30
    if n_elements(charsize) eq 0 then charsize = 1.0
	if n_elements(output) eq 0 then output = 'png'
	if n_elements(imgsizex) eq 0 then imgsizex = 1000
	if n_elements(imgsizey) eq 0 then imgsizey = imgsizex
	if n_elements(fprefix) eq 0 then fprefix = ''
	if n_elements(fsuffix) eq 0 then fsuffix = ''
	if n_elements(revolvetype) eq 0 then revolvetype = 'z'
	if n_elements(thrtype) eq 0 and n_elements(thrvar) gt 0 then thrtype=make_array(n_elements(thrvar),/string,value='min')
	if n_elements(special) eq 0 then special = ''
	if n_elements(zoom) eq 0 then zoom = 1.0
	if n_elements(boxscale) eq 0 then boxscale=[1.0,1.0,1.0]
	if n_elements(colorbarcolor) eq 0 then colorbarcolor='white'
	if n_elements(oversample) eq 0 then oversample = 1
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
						ct_part_len[i+1] = floor(256.*(alog10(ctswitch)-alog10(rangemin))/(alog10(rangemax) - alog10(rangemin)))	
					endif else begin
						ct_part_len[i+1] = floor(256.*ctswitch/(rangemax - rangemin))	
					endelse
				endif else begin
					if keyword_set(log) then begin
						ct_part_len[i+1] = floor(256.*(alog10(rangemax)-alog10(ctswitch))/(alog10(rangemax) - alog10(rangemin)))	
					endif else begin
						ct_part_len[i+1] = floor(256.*(rangemax - ctswitch)/(rangemax - rangemin))	
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
	
	;set_plot, 'z'
	;erase
	;device, set_resolution=[1000,1000]
	;device, /z_buffer

	if (n_elements(extdata) eq 0) then begin
		if (n_elements(var) gt 1) then begin
			for i=0,n_elements(var)-1 do begin
				load_flash_var, tmpslice, filename, var[i], xrange, yrange, zrange, dens=dens, temp=temp, $
					velx=velx, vely=vely, velz=velz, gpot=gpot, log=log, sample=sample, time=time, simsize=simsize, subtractavg=subtractavgi, $
					xcoord=xcoord, ycoord=ycoord, zcoord=zcoord, memefficient=memefficient
				if i eq 0 then slice = tmpslice else slice = slice*tmpslice
			endfor
		endif else begin
			load_flash_var, slice, filename, var, xrange, yrange, zrange, dens=dens, temp=temp, $
				velx=velx, vely=vely, velz=velz, gpot=gpot, log=log, sample=sample, time=time, simsize=simsize, subtractavg=subtractavgi, $
				xcoord=xcoord, ycoord=ycoord, zcoord=zcoord, memefficient=memefficient
		endelse
	endif else begin
		slice = extdata
	endelse
	if keyword_set(mirror) then begin
		mirrorslice = [slice,reverse(slice, 1)]
		slice = temporary(mirrorslice)
		oldrange = xrange
		xrange[1] = 2.*xrange[1] - xrange[0]
	endif
	s = size(slice,/dim)  
	
	if n_elements(thrvar) gt 0 then thrslice = fltarr(s[0],s[1],s[2],n_elements(thrvar))

	;below may not work for diagonal slices
	for i=0,n_elements(thrvar)-1 do begin
		if var eq thrvar[i] then begin
			thrslice[*,*,*,i] = slice
		endif else begin
			thrloaded = 0
			for j=0,i-1 do begin
				if thrvar[i] eq thrvar[j] then begin
					thrslice[*,*,*,i] = thrslice[*,*,*,j]
					thrloaded = 1
				endif
			endfor
			if thrloaded eq 0 then begin
				load_flash_var, newthrslice, filename, thrvar[i], xrange, yrange, zrange, dens=dens, temp=temp, $
					velx=velx, vely=vely, velz=velz, gpot=gpot, sample=sample
				thrslice[*,*,*,i] = newthrslice
			endif	
		endelse
	endfor

	if n_elements(thrvar) eq 0 then begin
		min_val = min(slice)
		max_val = max(slice)
	endif else begin
		for i=0,n_elements(thrvar)-1 do begin
			if thrtype[i] eq 'min' then begin
				indices = where(reform(thrslice[*,*,*,i]) ge thrval[i], count)
			endif else begin
				indices = where(reform(thrslice[*,*,*,i]) le thrval[i], count)
			endelse
			if count ne 0 then begin
				if i eq 0 then combindices = indices else begin
					combindices = cmset_op(combindices,'and',indices)
				endelse
			endif
		endfor
		min_val = min(slice[combindices])
		max_val = max(slice[combindices])
	endelse
	
	if n_elements(thrvar) gt 0 then begin
		for i=0,n_elements(thrvar)-1 do begin
			if thrtype[i] eq 'min' then begin
				indices = where(reform(thrslice[*,*,*,i]) lt thrval[i], count)
			endif else begin
				indices = where(reform(thrslice[*,*,*,i]) gt thrval[i], count)
			endelse
			if count ne 0 then begin
				if i eq 0 then thrindices = indices else begin
					thrindices = cmset_op(thrindices,'and',indices)
				endelse
			endif
		endfor
		if n_elements(rangemin) eq 0 then begin
			slice[thrindices] = min_val
		endif else begin
			slice[thrindices] = min([min_val, rangemin])
		endelse
	endif

	if n_elements(rangemin) ne 0 then begin
		indices = where(slice lt rangemin, count)
		if count ne 0 then slice[indices] = rangemin
	end
	if n_elements(rangemax) ne 0 then begin
		indices = where(slice gt rangemax, count)
		if count ne 0 then slice[indices] = rangemax
	end

	if keyword_set(log) then slice = alog10(slice)
	
	;slicemin = min(slice2)
	;slicemax = max(slice2)
	;shadesmin = min(myshades)
	;shadesmax = max(myshades)
	;minr = slicemin + shadesmin/255.*(slicemax-slicemin)
	;maxr = slicemin + shadesmax/255.*(slicemax-slicemin)
	;myshades = floor(double(myshades - shadesmin)/double(shadesmax-shadesmin)*255.)
	;image = polyshade(v, p, shades=myshades, /T3D)  
	if n_elements(rangemin) eq 0 then begin
		if n_elements(log) eq 0 then plotmin = min_val else plotmin = alog10(min_val)
	endif else begin
		if n_elements(log) eq 0 then plotmin = rangemin else plotmin = alog10(rangemin)
	endelse
	if n_elements(rangemax) eq 0 then begin
		if n_elements(log) eq 0 then plotmax = max_val else plotmax = alog10(max_val)
	endif else begin
		if n_elements(log) eq 0 then plotmax = rangemax else plotmax = alog10(rangemax)
	endelse
	vol = bytscl(slice, max=plotmax, min=plotmin)
	
	fsc_undefine, dens, temp, velx, vely, velz, gpot
	if n_elements(revolvestep) eq 0 then begin
   		revolvestep = 360
	endif
	for ang=0,360-revolvestep,revolvestep do begin
		if revolvestep lt 360 then begin
			fsuffix='_'+string(ang, format='(I3.3)')
			if revolvetype eq 'z' then begin
				aaz = az+ang
				aax = ax
			endif else begin
				aax = ax+ang
				aaz = az
			endelse
		endif else begin
			aax = ax
			aaz = az
		endelse
		if n_elements(custom_ct) ne 0 then tvlct, custom_ct else begin
			loadct, my_ct
		endelse
		tvlct, red, green, blue, /get
		thisDevice = !d.name
		set_plot, 'z'
		plotcolor = fsc_color('white',/nodisplay)
		!P.background = fsc_color('black',/nodisplay)
		erase, color=fsc_color('black',/nodisplay)
		if output ne 'x' then begin
			logf = ''
			if keyword_set(log) then logf = '_log'
			varname = ''
			for i=0,n_elements(var)-1 do begin
				varname = varname + var[i]	
			endfor
			jps_start, filename=fprefix+varname+'_vol'+logf+'_'+filename+fsuffix+'.ps',xsize=imgsizex/100.,ysize=imgsizey/100., /inches
		endif else begin
		    set_plot, 'x'
			window,1,xsize=1000,ysize=1000
		endelse

		if n_elements(boxaxes) ne 0 then begin
			scl = ((boxaxes[1]-boxaxes[0])/cellsize-s[0])/2.
			scale3, xrange=[-scl,s[0]+scl], yrange=[-scl,s[1]+scl], zrange=[-scl,s[2]+scl], ax=aax, az=aaz
		endif else begin
			scl = 0.0
			normdim = 1./sqrt(s[0]^2.+s[1]^2.+s[2]^2.)
			t3d, /reset, translate=-[0.5,0.5,0.5], $
				scale=[double(s[0])*normdim*boxscale[0],double(s[1])*normdim*boxscale[1],double(s[2])*normdim*boxscale[2]]*zoom
			;t3d, scale=replicate(1./sqrt(3),3)
			if n_elements(custom_rot) then begin
				t3d, rotate=custom_rot
			endif else begin
				t3d, rotate=[-90,aaz,0]
				t3d, rotate=[aax,0,0]
			endelse
			;t3d, translate=[double(s[0])/s[2]/2,double(s[1])/s[2]/2,0.0]
			t3d, translate=[0.5,0.5,0.5]
		endelse

		surface, dist(s[1]), /nodata, $
			xrange=[-scl,s[0]+scl], yrange=[-scl,s[1]+scl], zrange=[-scl,s[2]+scl], skirt=0.0, $
			xstyle=5, ystyle=5, zstyle=5, xmarg=[0,0], ymarg=[0,0], /t3d
		;surface, dist(s[1]), /nodata, $
		;	xrange=[0,s[0]], yrange=[0,s[1]], zrange=[0,s[2]], $
		;	/save, xstyle=5, ystyle=5, zstyle=5, charsize=charsize, skirt=0.0, $
		;	ax=ax, az=az, xmarg=[0,0], ymarg=[0,0]
		;t3d, /reset
		;t3d, translate=-[0.5,0.5,0.5]
		;t3d, rotate=[-ax,0,0]
		;t3d, rotate=[0,ang,0]
		;t3d, rotate=[ax,0,0]
		;t3d, translate=[0.5,0.5,0.5]
		;t3d, scale=[0.8,0.8,0.8]
		;t3d, translate=[0.1,0.1,0.1]

		;t3d, scale=[0.9, 0.9, 0.9], translate=[0.05,0.1,0.0]
		;scale3, xrange=[0,s[0]], yrange=[0,s[1]], zrange=[0,s[2]]

		rgbo = bytarr(256,2)
		rgbo[*,0] = 0.0
		rgbo[1:255,0] = ceil(dindgen(255)/255.)
		;rgbo[1:254,0] = ceil(5.*(make_array(254,value=1000.0,/double)))
		rgbo[*,1] = rgbo[*,0]
		;zpix = tvrd()
		;zbuffer = tvrd(/words,/chan)

		device,xsize=imgsizex/1000.,ysize=imgsizey/1000.,/inches ;voxel_proj assumes 1000 pixels per inch, for some reason
		;vox = voxel_proj(vol, background=[0,0,0], /interpolate, /maximum_intensity)
		slice = (slice - min(slice))/(max(slice) - min(slice))
		vox = jproject_vol(slice, s[0]*oversample, s[1]*oversample, s[2]*oversample)
		vox = round((vox - min(vox))/(max(vox) - min(vox))*255.0)
		
		device,xsize=imgsizex/100.,ysize=imgsizey/100.,/inches

		loadct, my_ct
		if imgsizex gt imgsizey then begin
			tvimage, vox, position=[0.0,(imgsizey-imgsizex)/imgsizey/2.,$
									1.0,1.0+(imgsizex-imgsizey)/imgsizey/2.]
		endif else begin
			tvimage, vox, position=[(imgsizex-imgsizey)/imgsizex/2.,0.0,$
									1.0+(imgsizey-imgsizex)/imgsizex/2.,1.0]
		endelse
			
		
		if special eq 'bh' then begin
			g = 6.673d-8
			c = 2.997d10

			; beta = 1 parameters (multitidal01)
			;m = 2.d39
			;q = 5.d12
			;tp = 25000.d0
			; beta = 7 parameters (tidal22/perisim)
			m = 2.d39
			q = 9.93571429d11
			tp = 10200.d0

			t = double(time)
			a = q^3.*sqrt(8.*q^3. + 9.*g*m*(t - tp)^2.) + 3.*sqrt(g*m)*q^3.*(t - tp)
			u = 2.*atan((-2.*q^3. + a^(2./3.))/(sqrt(2.)*q^1.5*a^(1./3.))) 
			r = 2.*q/(1+cos(u))
			x = -double(r*sin(u))
			y = double(r*cos(u))
			v = sqrt(2.*g*m/r)
			vx = cos(atan((y-q)/x))*v
			vy = sin(atan((y-q)/x))*v
			dx = x + simsize/2. + xrange[0] ;determine distance to bh from lower-left corner of frame
			dy = y + simsize/2. + yrange[0]
			dz = simsize/2. + zrange[0]
			if keyword_set(boxaxes) then begin
				dx = dx+simsize-xrange[1]
				dy = dy+simsize-yrange[1]
				dz = dz+simsize-zrange[1]
			endif

			xpos = dx/cellsize
			ypos = dy/cellsize
			zpos = dz/cellsize

			plots, [s[0]/2.,s[0]/2.], [s[1]/2.,s[1]/2.], [-scl,s[2]/2.], /t3d, $
				color=fsc_color('blue'), /data, linestyle=2
			if keyword_set(boxaxes) then begin
				xpos = xpos - s[0]/2.
				ypos = ypos - s[1]/2.
				zpos = zpos - s[2]/2.
				fac = simsize/(boxaxes[1]-boxaxes[0])
			endif else begin
				fac = 1.0
			endelse
			
			plots, [xpos,xpos], [ypos,ypos], [-scl,zpos], /t3d, $
				color=fsc_color('blue'), /data, linestyle=2
			plots, [xpos,xpos], [ypos,s[1]/2.], [-scl,-scl], /t3d, $
				color=fsc_color('blue'), /data, linestyle=2
			plots, [xpos,s[0]/2.], [s[1]/2.,s[1]/2.], [-scl,-scl], /t3d, $
				color=fsc_color('blue'), /data, linestyle=2
				
			mysym = findgen(49) * (!pi*2/48.)  
			usersym, cos(mysym), sin(mysym), thick=3.0
			if n_elements(cellsize) eq 0 then cellsize = simsize/s[0]
			rs = 2*g*m/c^2.
			plots, xpos, ypos, zpos, psym=8, symsize=fac*rs/cellsize*0.42, /t3d, color=fsc_color('red'), /data
			usersym, cos(mysym), sin(mysym), /fill
			plots, xpos, ypos, zpos, psym=8, symsize=fac*rs/cellsize/10., /t3d, color=fsc_color('blue'), /data
		endif

		if n_elements(boxaxes) ne 0 then begin
			xr = boxaxes
			yr = boxaxes
			zr = boxaxes
		endif else begin
			xr = xrange
			yr = yrange
			zr = zrange
		endelse

		if keyword_set(relaxes) then begin
			xr = relaxes[*,0]
			yr = relaxes[*,1]
			zr = relaxes[*,2]
		endif

		if keyword_set(hideticklabels) then begin
			xtickname = replicate(' ', 20)
			ytickname = replicate(' ', 20)
			ztickname = replicate(' ', 20)
		endif
		if ~keyword_set(hideaxes) then begin
			axis, -scl, s[1]+scl, -scl, zaxis=2, /t3d, zrange=zr, charsize=2.0*charsize, zsty=1, $
				color=fsc_color('white',/nodisplay), zticks=zticks, xtickname=xtickname, ytickname=ytickname, ztickname=ztickname
			axis, -scl, -scl, -scl, yaxis=0, /t3d, yrange=yr, charsize=2.0*charsize, ysty=1, $
				color=fsc_color('white',/nodisplay), yticks=yticks, xtickname=xtickname, ytickname=ytickname, ztickname=ztickname
			axis, -scl, -scl, -scl, xaxis=0, /t3d, xrange=xr, charsize=2.0*charsize, xsty=1, $
				color=fsc_color('white',/nodisplay), xticks=xticks, xtickname=xtickname, ytickname=ytickname, ztickname=ztickname
		endif

		tvlct, red, green, blue
		if ~keyword_set(hidecolorbar) then begin
			fsc_colorbar, /vertical, minrange=plotmin, $
				maxrange=plotmax, position=[0.89, 0.75, 0.92, 0.91], $
				annotatecolor=colorbarcolor, color=255,format='(G10.3)'
		endif

		if ~keyword_set(hidetime) then xyouts, 0.03, 0.96, 'Time' + string(time) + 's', charsize=1.5*charsize, $
			/normal, color=fsc_color('white',/nodisplay)
		if output eq 'png' then jps_end, /png else if output eq 'ps' then jps_end
	endfor
	fsc_undefine, vol

	if keyword_set(mirror) then xrange = oldrange
	;set_plot, thisDevice
	;tv, snapshot
end
