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
; my_ct         (int)             - IDL color table to use, 1-45 are valid selections.
; xrange        ([dbl, dbl])      - 2 values indicating range of x values to extract.
; yrange        ([dbl, dbl])      - 2 values indicating range of y values to extract.
; zrange        ([dbl, dbl])      - 2 values indicating range of z values to extract.
;
; **Optional**
; simsize       (dbl)             - Width of the box in cgs units (assumes box is cubical, used for certain variables).
; slicetype     (str)             - UNUSED.
; rangemin      (dbl)             - Minimum value to plot.
; rangeminstep  (dbl)             - Step the rangemin this much each frame.
; rangeminlstep (dbl)             - Step the log of rangemin this much each frame.
; rangemax      (dbl)             - Maximum value to plot.
; imgsizex      (dbl)             - Size of x dimension of output image, in pixels (applies to PNG only).
; imgsizey      (dbl)             - Size of y dimension of output image, in pixels (applies to PNG only).
; fprefix       (str)             - Prefix to append to output image files, e.g. 'myprefix_restoffilename'
; exactsize     (bool)            - Smallest grid cells will span exactly 1 pixel in output image.
; exactmult     (dbl)             - Exactsize must be true; smallest grid cells will span 1 pixel * exactmult.
; xticks        (int)             - Number of ticks to draw on x axis
; yticks        (int)             - Number of ticks to draw on y axis.
; zticks        (int)             - Number of ticks to draw on z axis (3D only).
; log           (bool)            - Scale data logarithmically.
; colorbarcolor (str)             - Color of colorbar annotation. 
; contours      (str)             - Variable used to draw contours on plot.
; thrvar        ([str,str,...])   - Variables used as a "threshold," data in regions where thresholds aren't
;                                   satisfied are set to rangemin.
; thrval        ([dbl,dbl,...])   - Threshold values, see above.
; thrtype       ([str,str,...])   - Whether the threshold is a minimum or a maximum, valid values 'min' or 'max'.
; sample        (int)             - Number of times to downsample data from highest refinment, must be <= 3.
; stride        (int)             - Stride of frame iteration, default 1.
; charsize      (dbl)             - Multiplier of default char size, default 1.0.
; symrange      (bool)            - Makes the data range symmetric depending on min and max values of data.
; lmin          (dbl)             - Used with symrange only, specifies min log value for
;                                   symmetric ranges (Since log can't have negative argument).
; annotatepos   (str)             - Annotation position (colorbar, timestamp). Currently can only set to 'flip'.
; output        (str)             - Specifies output file type, can be 'ps' or 'png' (default).
; ax            (dbl)             - Rotation about the x-axis (Volume plot only).
; az            (dbl)             - Rotation about the z-axis (Volume plot only).
; revolvestep   (dbl)             - Generate a series of frames for each data file, revolving this many degrees along
;							        the axis specified by revolvetype (Volume plot only).
; revolvetype   (str)             - The axis to revolve around (Volume plot only).
; special       (str)             - Flag used to indicate "special" plots to generate.
; cellsize      (dbl)             - Size of the smallest grid cells in CGS units, used for some variables.
; boxaxes       ([dbl,dbl])       - Draw data in a bounding box of this size (CGS units, volume plot only).
; memefficient  (bool)            - Deallocate variables as soon as they are not needed. Incompatible with revolve plots.
; hideaxes      (bool)            - Do not show axes.
; zoom          (dbl)             - Zoom plot by this factor, (Volume plot only).
; boxscale      (dbl)             - Scale each axis by this factor, (Volume plot only).
; indexlength   (int)             - Number of digits in filename index, e.g. _0000 is 4, _00000 is 5, etc.
; flythrough    ([[dbl,dbl,
;                  dbl,int],...]) - Table of coordinates used for generating flythrough movies. First two doubles are rotations
;                                   about the x and z axes, third double is zoom, integer is frame to interpolate the motion
;                                   to.
; trackfile     (str)             - File with x-y coordinates used to keep view centered.
; negative      (bool)            - Multiply the requested variable by -1.
; subtractavg   (bool)            - Subtract the average value from all grid cells when loading variable.
; showblocks    (bool)            - Show block boundaries
; oversample    (int)             - Interpolate the data to generate an oversampled image.
;
; EXAMPLES
; FLASH output files have the general syntax "simulationname_hdf5_plt_cnt_0000". Basename should be set to everything
; before the "_0000" in the filename, e.g. "simulationname_hdf5_plt_cnt". Here is a sample call that will generate
; frames of the logarithm of the density variable from every 5th data file numbered 0000 to 0100:
;
; make_flash_frames, 'mysimulation_hdf5_plt_cnt', 0, 100, 'dens', 3, [0,100], [0,100], 50, stride=5, /log
;
; Now lets do the same thing, but this time generate a volume plot:
;
; make_flash_frames, 'mysimulation_hdf5_plt_cnt', 0, 100, 'dens', 3, [0,100], [0,100], [0,100], /vol, stride=5, /log
;

pro make_flash_frames,basename,start,finish,var,my_ct,xrange,yrange,zrange,simsize=simsize,$
	slicetype=slicetype,rangemin=rangemin,rangemax=rangemax,imgsizex=imgsizex,imgsizey=imgsizey,$
	fprefix=fprefix,exactsize=exactsize,exactmult=exactmult,xticks=xticks,yticks=yticks,zticks=zticks,$
	log=log,colorbarcolor=colorbarcolor,contours=contours,thrvar=thrvar,thrval=thrval,sample=sample,lwant=lwant,stride=stride,$
	charsize=charsize,symrange=symrange,lmin=lmin,annotatepos=annotatepos,output=output,hidetime=hidetime,$
	ax=ax,az=az,revolvestep=revolvestep,revolvetype=revolvetype,thrtype=thrtype,special=special,cellsize=cellsize,$
	boxaxes=boxaxes,memefficient=memefficient,hideaxes=hideaxes,zoom=zoom,boxscale=boxscale,indexlength=indexlength,$
	flythrough=flythrough,negative=negative,subtractavg=subtractavg,custom_rot=custom_rot,hidecolorbar=hidecolorbar,$
	ctswitch=ctswitch,mirror=mirror,hideticklabels=hideticklabels,rminstep=rminstep,rminlstep=rminlstep,excision=excision,$
	fieldvarx=fieldvarx,fieldvary=fieldvary,fieldmax=fieldmax,refcoor=refcoor,absval=absval,trackfile=trackfile,showblocks=showblocks,$
	showrelaxes=showrelaxes,oversample=oversample,orbinfo=orbinfo,showpt=showpt,ptradius=ptradius,timeunit=timeunit

	compile_opt idl2
	if n_elements(indexlength) eq 0 then indexlength = 4
	if n_elements(special) eq 0 then special = ''
	if n_elements(timeunit) eq 0 then timeunit = 's'

	if n_elements(stride) eq 0 then begin
		if (start gt finish) then begin
			stride = -1
		endif else begin
			stride = 1
		endelse
	endif

	if n_elements(imgsize) eq 0 then imgsize = 1000

	if ((n_elements(xrange) eq 2 and $
		n_elements(yrange) eq 2 and $
		n_elements(zrange) eq 2) and $
		special ne 'revolve_z' and special ne 'column_x' and $
		special ne 'column_y' and special ne 'column_z') then vol = 1

	if (n_elements(xrange) eq 1) then xrange = [xrange, xrange]
	if (n_elements(yrange) eq 1) then yrange = [yrange, yrange]
	if (n_elements(zrange) eq 1) then zrange = [zrange, zrange]

	if keyword_set(vol) then begin
		if n_elements(ax) eq 0 then ax = 30
		if n_elements(az) eq 0 then az = 30
		if n_elements(zoom) eq 0 then zoom = 1.0
		curfly = [ax,az,zoom,0]
	endif

	if n_elements(trackfile) ne 0 then begin
		track = read_ascii(trackfile)
		track = track.(0)
		trackt = track[0,*]
		trackx = track[1,*]
		tracky = track[2,*]
		trackz = track[3,*]
		if keyword_set(showpt) then begin
			pttrackx = track[4,*]
			pttracky = track[5,*]
			pttrackz = track[6,*]
		endif
	endif

	numformat = '(I' + string(indexlength) + '.' + string(indexlength) + ')'

	if keyword_set(special) then begin
		if special eq 'diff' then begin
			base_files = indgen(10) + 200
			base_file = basename + '_' + string(base_files[0], format=numformat)
			read_amr, base_file, var_name=var, parameters=params
			time = params.time
			bxr = xrange + interpol(trackx, trackt, time)
			byr = yrange + interpol(tracky, trackt, time)
			bzr = zrange + interpol(trackz, trackt, time)
			bptxr = xrange + interpol(pttrackx, trackt, time)
			bptyr = yrange + interpol(pttracky, trackt, time)
			bptzr = zrange + interpol(pttrackz, trackt, time)
			print, 'Loading base state'
			base_state = double(loaddata(base_file,var,xrange=bxr,yrange=byr,zrange=bzr,sample=sample,lwant=lwant,xcoords=xcoords,ycoords=ycoords,zcoords=zcoords,/double))
			refcoor = [[bxr[0], bxr[1]], [byr[0], byr[1]], [bzr[0], bzr[1]]]
			for i=1, n_elements(base_files)-1 do begin
				base_file = basename + '_' + string(base_files[i], format=numformat)
				read_amr, base_file, var_name=var, parameters=params
				time = params.time
				bxr = xrange + interpol(trackx, trackt, time)
				byr = yrange + interpol(tracky, trackt, time)
				bzr = zrange + interpol(trackz, trackt, time)
				bptxr = xrange + interpol(pttrackx, trackt, time)
				bptyr = yrange + interpol(pttracky, trackt, time)
				bptzr = zrange + interpol(pttrackz, trackt, time)
				tmp_state = double(loaddata(base_file,var,xrange=bxr,yrange=byr,zrange=bzr,sample=sample,lwant=lwant,xcoords=xcoords,ycoords=ycoords,zcoords=zcoords,/double))
				mc = double(xcoords[n_elements(xcoords)-1] - xcoords[0])/n_elements(xcoords)
				x0 = (bxr[0] - refcoor[0,0]) mod mc
				y0 = (byr[0] - refcoor[0,1]) mod mc
				z0 = (bzr[0] - refcoor[0,2]) mod mc
				nxc = ((indgen(n_elements(xcoords))) * mc + x0)/mc
				nyc = ((indgen(n_elements(ycoords))) * mc + y0)/mc
				nzc = ((indgen(n_elements(zcoords))) * mc + z0)/mc
				base_state = base_state + interpolate(tmp_state, nxc, nyc, /grid, cubic=-0.5)
			endfor
			base_state = base_state/n_elements(base_files)
		endif
	endif

	for i=start,finish,stride do begin
		filename = basename + '_' + string(i, format=numformat)
		if file_test(filename) eq 0 then begin
			print, "Can't find " + filename + ", skipping."
			continue
		endif
		if n_elements(trackfile) eq 0 then begin
			xr = xrange
			yr = yrange
			zr = zrange
		endif else begin
			jread_amr, filename, var_name='none', parameters=params
			time = params.time
			xr = xrange + interpol(trackx, trackt, time, /quad)
			yr = yrange + interpol(tracky, trackt, time, /quad)
			zr = zrange + interpol(trackz, trackt, time, /quad)
			if keyword_set(showpt) then begin
				ptpos = [interpol(pttrackx, trackt, time), $
						 interpol(pttracky, trackt, time), $
						 interpol(pttrackz, trackt, time)]
			endif
			relaxes = [[xrange[0], xrange[1]], [yrange[0], yrange[1]], [zrange[0], zrange[1]]]
		endelse
		if keyword_set(vol) then begin
			make_flash_vol,filename,var,my_ct,xr,yr,zr,ax=curfly[0],az=curfly[1],log=log,rangemin=rangemin,rangemax=rangemax,$
				output=output,sample=sample,lwant=lwant,thrvar=thrvar,thrval=thrval,charsize=charsize,revolvestep=revolvestep,hideticklabels=hideticklabels,$
				revolvetype=revolvetype,imgsizex=imgsizex,imgsizey=imgsizey,fprefix=fprefix,thrtype=thrtype,special=special,$
				simsize=simsize,cellsize=cellsize,boxaxes=boxaxes,hideaxes=hideaxes,zoom=curfly[2],boxscale=boxscale,ctswitch=ctswitch,$
				xticks=xticks,yticks=yticks,zticks=zticks,custom_rot=custom_rot,hidecolorbar=hidecolorbar,hidetime=hidetime,mirror=mirror,$
				relaxes=relaxes,oversample=oversample
			
			if keyword_set(flythrough) then begin
				if i eq start then flythrough = [[0.0,0.0,0.0,start],[flythrough]]
				for j=1,n_elements(flythrough)/3-1 do begin
					if flythrough[3,j] gt i then begin
						curfly = curfly + flythrough[0:2,j]/(flythrough[3,j]-flythrough[3,j-1])
						break
					endif
				endfor
			endif
		endif else begin
			make_flash_slice,filename,var,my_ct,xr,yr,zr,simsize=simsize,slicetype=slicetype,$
				rangemin=rangemin,rangemax=rangemax,contours=contours,thrvar=thrvar,thrval=thrval,sample=sample,lwant=lwant,$
				log=log,colorbarcolor=colorbarcolor,imgsizex=imgsizex,imgsizey=imgsizey,hidetime=hidetime,$
				fprefix=fprefix,exactsize=exactsize,exactmult=exactmult,charsize=charsize,xticks=xticks,yticks=yticks,$
				ambval=ambval,symrange=symrange,lmin=lmin,annotatepos=annotatepos,output=output,hidecolorbar=hidecolorbar,$
				special=special,negative=negative,subtractavg=subtractavg,thrtype=thrtype,hideaxes=hideaxes,ctswitch=ctswitch,$
				excision=excision,fieldvarx=fieldvarx,fieldvary=fieldvary,fieldmax=fieldmax,absval=absval,refcoor=refcoor,$
				showblocks=showblocks,relaxes=relaxes,base_state=base_state,orbinfo=orbinfo,trackfile=trackfile,$
				memefficient=memefficient,ptpos=ptpos,ptradius=ptradius,timeunit=timeunit
		endelse
		fsc_undefine, filename
		if keyword_set(rminstep) then rangemin = rangemin + rminstep
		if keyword_set(rminlstep) then rangemin = 10^(alog10(rangemin) + rminlstep)
	end
end
