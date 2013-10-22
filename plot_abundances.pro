pro plot_abundances
	data = read_ascii('sm.dat',data_start=1)
	data = data.field01
	!p.color = 255
	elems = ['n',  'H',  'D',  '3He',  '4He',  '6Li', '7Li', '7Be', '9Be', '8B', $
			 '10B', '11B', '11C', '12C', '13C', '12N', '14N', '15N', '16O', '17O', $
			 '18O', '20Ne', '21Ne', '22Ne', '23Na', '24Mg', '25Mg', '26Mg', '27Al', $
			 '28Si', '29Si', '30Si', '56Fe', '19F']
	for i=13,46 do begin
		if i eq 13 then plot, data(1,*), alog10(data(13,*)), yrange=[-6,0]
		if i ne 13 then oplot, data(1,*), alog10(data(i,*)), linestyle=(i-12) mod 5
		ind = floor(n_elements(data(1,*))/50.*(i-13.))
		print, ind
		xyouts, data(1,ind), alog10(data(i,ind)), elems(i-13)
	endfor
end
