pro make_flash_sixplot,filename1,filename2,filename3,var,my_ct,thresh_var,thresh_val,xrange,yrange,zrange
	set_plot, 'ps'
	device, bits_per_pixel=8, color=1, filename='tidal_6plots.ps', xsize=10.5, ysize=(21./3.), $
		xoffset=0.25, yoffset=10.75, /inches, /landscape
	;Window,1,XSIZE=1200,YSIZE=800
	!p.multi = [0, 2, 3]
	slicerange = fltarr(2)
	slicerange[0] = (zrange[0] + zrange[1]) / 2.
	slicerange[1] = (zrange[0] + zrange[1]) / 2.

	make_flash_slice,filename1,var,my_ct,thresh_var,thresh_val,xrange,yrange,slicerange,'z',[0.04, 0.545, 0.33, 0.98]
	make_flash_slice,filename1,var,my_ct,thresh_var,thresh_val,xrange,slicerange,zrange,'y',[0.04, 0.06, 0.33, 0.495]
	make_flash_slice,filename2,var,my_ct,thresh_var,thresh_val,xrange,yrange,slicerange,'z',[0.37, 0.545, 0.66, 0.98]
	make_flash_slice,filename2,var,my_ct,thresh_var,thresh_val,xrange,slicerange,zrange,'y',[0.37, 0.06, 0.66, 0.495]
	make_flash_slice,filename3,var,my_ct,thresh_var,thresh_val,xrange,yrange,slicerange,'z',[0.70, 0.545, 0.99, 0.98]
	make_flash_slice,filename3,var,my_ct,thresh_var,thresh_val,xrange,slicerange,zrange,'y',[0.70, 0.06, 0.99, 0.495]
	device, /close
end
