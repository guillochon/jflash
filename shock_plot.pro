pro shock_plot, filename, zrange

unk_names = get_var_list(filename)
dvar = where(strmatch(unk_names, 'dens') eq 1)
vzvar = where(strmatch(unk_names, 'velz') eq 1)
read_amr, filename, VAR_NAME=var, TREE=tree, DATA=unk, PARAMETERS=params

time = params.time
numblocks = params.totblocks

max_ref = max(tree[*].lrefine)
nc = 8
arr_dim = 8*2^(max_ref - 1)
cell_size = (tree[0].bndbox[1,0] - tree[0].bndbox[0,0])/arr_dim
arr_z_dim = floor((zrange[1] - zrange[0]) / cell_size)
shock_zs = indgen(arr_z_dim) + round(zrange[0]/cell_size)

shock_profs = dblarr(arr_dim, arr_dim, arr_z_dim, 2)

for n=0L, (numblocks-1) do begin
	if (tree[n].nodetype NE 1) then continue
	if (tree[n].lrefine NE max_ref) then continue
	if (tree[n].bndbox[0,2] lt zrange[0]) then continue
	dens = reform(unk[dvar,n,*,*,*])
	velz = reform(unk[vzvar,n,*,*,*])
	;if (max(dens) LT 1.e-2) then continue

	xlocs = [round(tree[n].bndbox[0,0]/cell_size),round(tree[n].bndbox[1,0]/cell_size)-1]
	ylocs = [round(tree[n].bndbox[0,1]/cell_size),round(tree[n].bndbox[1,1]/cell_size)-1]
	zlocs = [round(tree[n].bndbox[0,2]/cell_size),round(tree[n].bndbox[1,2]/cell_size)-1]

	if (zlocs[0] gt shock_zs[arr_z_dim-1]) then continue
	if (zlocs[1] lt shock_zs[0]) then continue
	zlocso = zlocs 
	zlocs = [round(max([shock_zs[0],zlocs[0]])),round(min([shock_zs[arr_z_dim-1], zlocs[1]]))] - shock_zs[0]

	shock_profs[xlocs[0]:xlocs[1],ylocs[0]:ylocs[1],zlocs[0]:zlocs[1],0] = $
		dens[*,*,round(min([max([0,shock_zs[0]-zlocso[0]]),nc-1])):round(nc-1-max([0,zlocso[1]-shock_zs[arr_z_dim-1]]))]
	shock_profs[xlocs[0]:xlocs[1],ylocs[0]:ylocs[1],zlocs[0]:zlocs[1],1] = $
		velz[*,*,round(min([max([0,shock_zs[0]-zlocso[0]]),nc-1])):round(nc-1-max([0,zlocso[1]-shock_zs[arr_z_dim-1]]))]
endfor

shock_profs = reform(shock_profs, long(arr_dim)*arr_dim, arr_z_dim, 2)
pruned_ind = where(shock_profs[*,0,0] ge 1.0e4 and shock_profs[*,0,1] gt 0.0) 
set_plot, 'ps'
device, file='shock_profs.ps'
;sample_ind = round(randomu(1,5000)*(n_elements(pruned_ind)-0.5))
sample_ind = pruned_ind
for n=0L, n_elements(sample_ind)-1 do begin
	;nn = sample_ind[n]
	nn = n
	cut = where(shock_profs[pruned_ind[nn],*] ne 0.0)
	if cut[0] eq -1 then continue
	if n eq 0 then begin
		plot, shock_profs[pruned_ind[nn],cut,1], yrange=[-1e9,1e9], thick=alog10(shock_profs[pruned_ind[nn],0,0])-4
	endif else begin
		oplot, shock_profs[pruned_ind[nn],cut,1], thick=alog10(shock_profs[pruned_ind[nn],0,0])-4
	endelse
	;print, shock_profs[pruned_ind[n],*]
	;if n gt 1000 then break
endfor
device, /close

end
