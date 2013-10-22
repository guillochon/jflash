function flatten_flash_slice,data,sliceplane
	if sliceplane eq 'x' then newdata = total(data,1)
	if sliceplane eq 'y' then newdata = total(data,2)
	if sliceplane eq 'z' then newdata = total(data,3)
	newdata = reform(newdata)
	return, newdata
end
