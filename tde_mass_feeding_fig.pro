pro tde_mass_feeding_fig
	file = 'thr_data.dat'
	readcol, file, filenum, minval, maxval, format='I,D,D', /silent
	nrows = n_elements(filenum)

	steps = 256
    greenVector = replicate(0, steps)
	scaleFactor = findgen(steps) / (steps - 1)
	beginNum = 255
	endNum = 0
	blueVector = beginNum + (endNum - beginNum) * scaleFactor
	redVector = 255 - blueVector

	for i=0,nrows-1 do begin
		make_flash_frames,'multipoly_hdf5_plt_cnt',filenum[i],filenum[i],'dens',0,[-4.e12,5.e11],[-1.e12,1.e12],[0,0],indexlength=5,/exactsize,fprefix='contbg_',$
			/hideaxes,/hidecolorbar,/hidetime,trackfile='track.dat',contours={var:'dens', min:0.99999e-6, max:1e-6, num:1, colortable:[[redVector], [greenVector], [blueVector]], colorindex:[round(255.*double(i)/(nrows-1))]}, $
			rangemin=1e-6,orbinfo=1.9889225d39,/log,output='ps'
		;make_flash_frames,'multipoly_hdf5_plt_cnt',filenum[i],filenum[i],'dens',0,[-4.e12,5.e11],[-1.e12,1.e12],[0,0],indexlength=5,/exactsize,fprefix='cont2_',$
		;	/hideaxes,/hidecolorbar,/hidetime,trackfile='track.dat',contours={var:'dens', min:0.99999e-6, max:1e-6, num:1, colortable:[[redVector], [greenVector], [blueVector]], colorindex:[round(255.*double(i)/(nrows-1))]}, $
		;	thrvar=['bhbound','bhbound','selfbound'],thrval=[minval[i],maxval[i],0.0],thrtype=['min','max','min'],$
		;	rangemin=1e-6,orbinfo=1.9889225d39,output='ps'
	endfor
end
