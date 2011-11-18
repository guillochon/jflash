pro make_flash_isosurf,filename,var,surfvar,my_ct,xrange,yrange,zrange,cut,$
	rangemin=rangemin,rangemax=rangemax,log=log,colorbarcolor=colorbarcolor,thrvar=thrvar,thrval=thrval,sample=sample,ps=ps

	window,1,xsize=1000,ysize=1000
	
	loadct, my_ct
	tvlct, red, green, blue, /get
	thisDevice = !d.name
	set_plot, 'z'
	erase
	device, set_resolution=[1000,1000]
	slice = loaddata(filename,var,XRANGE=xrange,YRANGE=yrange,ZRANGE=zrange,sample=sample,TIME=time)
	slice = reform(slice)
	
	slice2 = loaddata(filename,surfvar,XRANGE=xrange,YRANGE=yrange,ZRANGE=zrange,sample=sample,TIME=time)
	slice2 = reform(slice2)
	if keyword_set(rangemin) then begin
		indices = where(slice2 lt rangemin, count)
		if count ne 0 then slice2[indices] = rangemin
	end
	if keyword_set(rangemax) then begin
		indices = where(slice2 gt rangemax, count)
		if count ne 0 then slice2[indices] = rangemax
	end
	if keyword_set(log) then slice2 = alog10(slice2)
	;slice2 = alog10(slice)
	myshades = bytscl(slice2)
	
	s = SIZE(slice)  
	surface, dist(30), /NoData, XRange=xrange, YRange=yrange, $
		ZRange=zrange, /NoErase, /Save, XStyle=1, $
	    YStyle=1, ZStyle=1
	shade_volume, slice, cut, v, p, shades=myshades, /LOW  
	; Obtain the dimensions of the volume.  
	; Variables S[1], S[2], and S[3] now contain  
	; the number of columns, rows, and slices in the volume:  
	; Use SCALE3 to establish the three-dimensional  
	; transformation matrix. Rotate 45 degrees about the z-axis:  
	scale3, xrange=[0,S[1]], yrange=[0,S[2]], zrange=[0,S[3]]  

	slicemin = min(slice2)
	slicemax = max(slice2)
	shadesmin = min(myshades)
	shadesmax = max(myshades)
	minr = slicemin + shadesmin/255.*(slicemax-slicemin)
	maxr = slicemin + shadesmax/255.*(slicemax-slicemin)
	myshades = floor(double(myshades - shadesmin)/double(shadesmax-shadesmin)*255.)
	image = polyshade(v, p, shades=myshades, /T3D)  
	
	colorbar, /vertical, minrange=minr, $
		maxrange=maxr, position=[0.89, 0.75, 0.92, 0.91], $
		annotatecolor=colorbarcolor, format='(E0.3)', color=255

	snapshot = tvrd()

	set_plot, thisDevice
	tv, snapshot
end
