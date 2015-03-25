pro flash_shrink, basename, start, finish, vars, xrange, yrange, zrange, sample=sample, stride=stride
	if n_elements(stride) eq 0 then stride = 1

	for i=start,finish,stride do begin
		filename = basename + '_' + string(i, format='(I4.4)')
		for j=0,n_elements(vars)-1 do begin
			var_data = loaddata(filename, vars[j], sample=sample, TIME=time, xrange=xrange, yrange=yrange, zrange=zrange)
			if j eq 0 then data = fltarr([n_elements(vars), size(var_data, /dimensions)])
			data[j, *, *, *] = var_data
		endfor

		fid = H5F_CREATE('shrunk_' + filename + '.h5')

		for j=0,n_elements(vars)-1 do begin
			datatype_id = H5T_IDL_CREATE(reform(data[j, *, *, *]))
			dataspace_id = H5S_CREATE_SIMPLE(size(reform(data[j, *, *, *]),/dimensions))
			dataset_id = H5D_CREATE(fid,vars[j],datatype_id,dataspace_id)

			H5D_WRITE,dataset_id,reform(data[j, *, *, *])

			H5D_CLOSE,dataset_id
			H5S_CLOSE,dataspace_id
			H5T_CLOSE,datatype_id
		endfor

		H5F_CLOSE,fid
	endfor
end
