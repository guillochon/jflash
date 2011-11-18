function mom_of_inertia,data,slicetype
	inertia = fltarr(2,2)
	dims = size(data, /DIMENSIONS)
	for i=0,dims[0]-1 do begin	
		for j=0,dims[1]-1 do begin	
			inertia[0,0] = inertia[0,0] + data[i,j]*(j-dims[1]*0.5)^2
			inertia[1,1] = inertia[1,1] + data[i,j]*(i-dims[0]*0.5)^2
			inertia[0,1] = inertia[0,1] - data[i,j]*(i-dims[0]*0.5)*(j-dims[1]*0.5)
		endfor
	endfor
	inertia[1,0] = inertia[0,1]
	;eval = eigenql(inertia, eigenvectors=evec, /absolute, /ascending)
	trired, inertia, eval, e
	triql, eval, e, inertia
	print, eval
	print, inertia
	if slicetype eq 'major' then begin
		if abs(eval[0]) gt abs(eval[1]) then begin
			angle = atan(inertia[0,1],inertia[0,0])
		endif else begin
			angle = atan(inertia[1,1],inertia[1,0])
		endelse
		;return,inertia[0,*]
	endif else if slicetype eq 'minor' then begin
		if abs(eval[0]) lt abs(eval[1]) then begin
			angle = atan(inertia[0,1],inertia[0,0])
		endif else begin
			angle = atan(inertia[1,1],inertia[1,0])
		endelse
		;return,inertia[1,*]
	endif
	if angle lt -!pi/2. then angle = angle + !pi
	if angle gt !pi/2. then angle = angle - !pi
	return, angle
end
