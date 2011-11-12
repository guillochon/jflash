pro make_track_file, filename, offset=offset, outfile=outfile, mode=mode, stride=stride, nsmooth=nsmooth
	if ~keyword_set(offset) then offset = 0
	if ~keyword_set(mode) then mode = 1
	if ~keyword_set(outfile) then outfile = 'track.dat'
	if ~keyword_set(stride) then stride = 1
	spawn, 'sed -e "/^\s*\#/d" -e "s/^\s*//" < ' + filename + ' > ' + 'temporary_track.dat'
	nrows = file_lines('temporary_track.dat')
	openr, lun, 'temporary_track.dat', /get_lun
	line = ""
	readf, lun, line
	ncols = n_elements(strsplit(line, /regex, /extract))
	data = dblarr(ncols,nrows)
	point_lun, lun, 0
	readf, lun, data
	free_lun, lun
	time = data[0,*]
	length = n_elements(time)
	i = 1
	res_cnt = 0
	while i lt length do begin
		if time[i] lt time[i-1] then begin
			include = make_array(length, value = 1, /integer)
			si = -1
			for j = 0, i - 1 do begin
				if time[j] eq time[i] then begin
					si = j
				endif
			endfor
			if si eq -1 then begin
				print, 'Error, restart detected but matching time not found!'
				return
			endif
			include[si:i-1] = 0

			data = data[*,where(include)]
			time = data[0,*]
			length = n_elements(time)

			i = 1
			res_cnt = res_cnt + 1
		endif
		i = i + 1
	endwhile

	if res_cnt gt 0 then print, string(res_cnt) + ' restarts removed.'

	time = time[0:length-1:stride]
	if mode eq 1 then begin
		; Old orbit.dat format
		path = data[25+offset:27+offset,0:length-1:stride]
		; New orbit.dat format
		;path = data[31+offset:33+offset,0:length:stride]
	endif else if mode eq 2 then begin
		path = data[25+offset:27+offset,0:length-1:stride]+data[1+offset:3+offset,0:length-1:stride]-data[7+offset:9+offset,0:length-1:stride]
	endif else if mode eq 3 then begin
		path = data[*,0:length-1:stride]
	endif else begin
		time = data[0,0:length-1:stride]
		path1 = data[25+offset:27+offset,0:length-1:stride]
		path2 = path1+data[1+offset:3+offset,0:length-1:stride]-data[7+offset:9+offset,0:length-1:stride]
	endelse

	openw, lun, outfile, /get_lun
	if mode eq 1 or mode eq 2 then begin
		newtime = time[0] + (time[length-1] - time[0])*dindgen(length)/(length - 1)
		px = interpol(path[0,*], time, newtime, /lsq)
		py = interpol(path[1,*], time, newtime, /lsq)
		pz = interpol(path[2,*], time, newtime, /lsq)
		time = newtime
		if (keyword_set(nsmooth)) then begin
			px = ts_smooth(px, nsmooth)
			py = ts_smooth(py, nsmooth)
			pz = ts_smooth(pz, nsmooth)
		endif
		for i=0L, length-1 do printf, lun, time[i], px[i], py[i], pz[i], format='(4G25.16)'
	endif else if mode eq 3 then begin
		ds = size(data, /dimensions)
		fs = '('+strtrim(string(ds[0]))+'G25.16)'
		for i=0L, length-1 do printf, lun, data[*,i], format=fs
	endif else begin
		fs = '(7G25.16)'
		for i=0L, length-1 do printf, lun, time[i], path1[*,i], path2[*,i], format=fs
	endelse
	free_lun, lun
	spawn, 'rm temporary_track.dat'
end
