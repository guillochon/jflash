;==============================================================================
;
; read a FLASH HDF5 file into IDL.  
;
; This routine reads in FLASH data in the HDF 5 format.  
;
; This version is for IDL >= 5.6, as it uses the built in HDF5
; support.
;
; Arguments:  filename -- the name of the file `to read in
;
;             var_name -- a 4-character string specifying the variable
;                         name to read in.  If this keyword is not
;                         present, all variables are read in.
;
;             /VERBOSE  -- output diagnostics while running
;
;==============================================================================

pro read_parameters, filename, $
              VERBOSE=verbose, $ 
              TREE=tree, $
              DATA=unk, $
              PARAMETERS=params, $
              STORED_VARS=unk_names, $
              NUM_PARTICLES=numParticle, $
              PARTICLES=particles, $
              INT_PROP_NAMES=IntPropNames, $
              REAL_PROP_NAMES=RealPropNames, $
              GEOMETRY=geometry

; clean up any pre-existing arrays
if n_elements(unk) GT 0 then undefine, unk
if n_elements(particles) gt 0 then undefine, particles
if (not Keyword_set(verbose)) then verbose = 0

; find out what software version we're running
flash_version = determine_flash_version(filename)

; There are no "logical" variables in IDL -- zero = false
if (flash_version EQ 2) then begin 
  FLASH3 = 0
  FLASH2 = 1
endif else begin 
  FLASH3 = 1
  FLASH2 = 0
endelse 

if Keyword_Set(verbose) then print,' FLASH Version is ', flash_version

;  determine the filetype of HDF5(1) or NetCDF(2)
file_type = determine_file_type(filename)

if (file_type EQ -1) then begin
  print, 'ERROR: file does not exist'
  return
endif

if Keyword_Set(verbose) then print, ' File Type is ',file_type

; DOUBLE is a flag to read data in double precision (checkpoint files)
; try to determine the precision from the file
file_precision = determine_file_precision(filename)
if (file_precision eq 1) then double = 0
if (file_precision eq 2) then double = 1

; check to see which variable to read; if var is not defined, read them all
if  n_elements(var_name) EQ 0 then begin
  var_name = 'none'
endif else begin 
  ; check that var name has four characters -- append with blanks if
  ; necessary.  format string pads with blanks on right if too short
  var_name = string(var_name, FORMAT='(A-4)')
endelse 

if Keyword_Set(verbose) then print, ' Variable var_name requested = ', var_name

;-----------------------------------------------------
;  begin reading
;---------------------------------------------------

;--------------------------------------------------------------------
; get the list of variables and make sure that our requested variable
; (if any) exists in the dataset
;-------------------------------------------------------------------
unk_names = get_var_list(filename)

;------------------------------------------------------------------
;  open the file
;------------------------------------------------------------------
if (file_type EQ 1) then file_identifier = H5F_OPEN(filename)
if (file_type EQ 2) then file_identifier = NCDF_OPEN(filename)

;------------------------------------------------------------------
; we can also now set the number of variables
;------------------------------------------------------------------
nvar = (size(unk_names))[1]
if (var_name NE 'none') then begin
  var = (where(unk_names EQ var_name))[0]
  if (var EQ -1) then begin
    print, 'ERROR: requested variable ',var_name,' not found in dataset'
    print, '       therefore reading in all variables'
  endif else begin 
    print, ' Reading in only the variable ', var_name
  endelse 
endif else begin
  var = -1
endelse

;-------------------------------------------------------------------
; read in the simulation parameters
;-------------------------------------------------------------------
if (file_type EQ 1) then begin 
  if(FLASH2) then begin

    dataset = H5D_OPEN(file_identifier, "simulation parameters")
    datatype = H5D_GET_TYPE(dataset)

    idl_type = H5T_IDLTYPE(datatype, STRUCTURE=sim_params)
    sim_params = H5D_READ(dataset)
    H5D_CLOSE, dataset
    H5T_CLOSE, datatype
  endif

  ; read in the integer scalars
  if (FLASH3) then begin
    dataset = H5D_OPEN(file_identifier, "integer scalars")
    datatype = H5D_GET_TYPE(dataset)
    idl_type = H5T_IDLTYPE(datatype, STRUCTURE=int_scalars_list)

    ; int_scalars_list is a linked list which we will change
    ; into a struct whose elements and values are the tag and
    ; and value names of the linked list
    int_scalars_list = H5D_READ(dataset)
    for i=1, (size(int_scalars_list))[3] do begin
      if ((size(int_scalars))[0] eq 0) then begin  ; if int_scalars doesn't exist yet, create it
        int_scalars = create_struct(strtrim(int_scalars_list[i-1].name), int_scalars_list[i-1].value)
      endif else begin  ; otherwise, append to it
        int_scalars = create_struct(int_scalars, strtrim(int_scalars_list[i-1].name), int_scalars_list[i-1].value)
      endelse
    endfor

    H5T_CLOSE, datatype
    H5D_CLOSE, dataset  ; "integer scalars"

    ; read in the real scalars
    dataset = H5D_OPEN(file_identifier, "real scalars")
    datatype = H5D_GET_TYPE(dataset)
    idl_type = H5T_IDLTYPE(datatype, STRUCTURE=real_scalars_list)

    ; real_scalars_list is a linked list which we will change
    ; into a struct whose elements and values are the tag and
    ; and value names of the linked list
    real_scalars_list = H5D_READ(dataset)
    for  i=1, (size(real_scalars_list))[3] do  begin 
      if ((size(real_scalars))[0] EQ  0) then  begin ; if real_scalars doesn't exist yet, create it
        real_scalars = create_struct(strtrim(real_scalars_list[i-1].name), real_scalars_list[i-1].value)
      endif else begin  ; otherwise, append to it`
        real_scalars = create_struct(real_scalars, strtrim(real_scalars_list[i-1].name), real_scalars_list[i-1].value)
      endelse 
    endfor 
    H5T_CLOSE, datatype
    H5D_CLOSE, dataset  ; "real scalars"
  endif                 ; end of readin for FLASH3
  ; DEV don't know where these are read in FLASH2

  ; get the dimensionality
  ; NOTE: below is the old way of getting ndim based on coordinate data
  ;       new way is based off of gid.
  ;dataset = H5D_OPEN(file_identifier, "coordinates")
  ;dataspace = H5D_GET_SPACE(dataset)
  ;dims = H5S_GET_SIMPLE_EXTENT_DIMS(dataspace)
  ;ndim = dims[0]
  dataset = H5D_OPEN(file_identifier, "gid")
  dataspace = H5D_GET_SPACE(dataset)
  dims = H5S_GET_SIMPLE_EXTENT_DIMS(dataspace)
  switch dims[0] of
    5: begin
       ndim = 1
       break
    end
    9: begin
       ndim = 2
       break
    end
    15: begin
        ndim  = 3
        break
    end
    else: message, 'gid size per block does not match 1, 2 or 3d values. Read will Fail.'
    endswitch


  H5S_CLOSE, dataspace
  H5D_CLOSE,dataset
endif else begin  ; end of hdf5, start of netCDF

  ; read in the simulation parameters
  if (FLASH3) then begin
    ncdf_attget, file_identifier, 1, "globalnumblocks",globalNumBlocks
    ncdf_attget, file_identifier, 1, "time", time
    ncdf_attget, file_identifier, 1, "dt", dt
    ncdf_attget, file_identifier, 1, "nstep", nsteps 
    ncdf_attget, file_identifier, 1, "nxb" , nxb
    ncdf_attget, file_identifier, 1, "nyb" , nyb
    ncdf_attget, file_identifier, 1, "nzb" , nzb

    int_scalars = create_struct('globalNumBlocks' , globalNumBlocks, $
                                'nsteps', nsteps, 'redshift', 1.0, $		
                                'nxb', nxb, 'nyb', nyb, 'nzb', nzb)
    real_scalars = create_struct('time',time,'dt',dt,'timestep',1.0)
  endif else begin  ; FLASH2 file
                    ;DEV: kda - this is the old format way
                    ;ncdf_attget, file_identifier, /GLOBAL, "timestep", timestep
                    ;ncdf_attget, file_identifier, /GLOBAL, "redshift", redshift 
    ncdf_attget, file_identifier, /GLOBAL, "total_blocks",total_blocks
    ncdf_attget, file_identifier, /GLOBAL, "time", time
    ncdf_attget, file_identifier, /GLOBAL, "timestep", timestep
    ncdf_attget, file_identifier, /GLOBAL, "nsteps", nsteps
    ncdf_attget, file_identifier, /GLOBAL, "redshift", redshift
    ncdf_attget, file_identifier, /GLOBAL, "nxb" , nxb
    ncdf_attget, file_identifier, /GLOBAL, "nyb" , nyb
    ncdf_attget, file_identifier, /GLOBAL, "nzb" , nzb

    sim_params = create_struct('total_blocks' , total_blocks, $
                               'nsteps', nsteps, 'redshift', redshift, $		
                               'time', time, 'timestep', timestep, $
                               'nxb', nxb, 'nyb', nyb, 'nzb', nzb)
  endelse  ; end of flash2/3 if block for netcdf 

  ; get the dimensionality
  dimid = ncdf_dimid(file_identifier, 'dim_NDIM')
  ncdf_diminq, file_identifier, dimid, nndim, ndim
endelse  ; end of netcdf

; Common -- get dimensions
nfaces = 2*ndim
nchild = 2^ndim

;---------------------------------------------------------
; figure out if we are dealing with corners
; DEV a hack in FLASH2, FLASH3 should write them
;---------------------------------------------------------
corners = 0
if (strpos(filename, 'crn_') GT 0) then begin
  corners = 1
endif

;------------------------------------------------------
;  get GEOMETRY
;-------------------------------------------------------

if (n_elements(geometry) eq 0) then begin
  geometry = determine_geometry(filename)
endif 

;------------------------------------------------------------------------------
;   setup the structures to pass the data
;----------------------------------------------------------------------------

; simulation parameters and int_scalars are different in flash2/3
if (FLASH3) then begin
  totBlocks = int_scalars.globalNumBlocks
  nxb = int_scalars.nxb
  nyb = int_scalars.nyb
  nzb = int_scalars.nzb
  time = real_scalars.time  ; double precision
  dt = real_scalars.dt      ; double precision
endif else begin            ; end flash3, begin FLASH2
  totBlocks = sim_params.total_blocks
  nxb = sim_params.nxb
  nyb = sim_params.nyb
  nzb = sim_params.nzb
  time = sim_params.time    ; double precision
  dt = sim_params.timestep  ; double precision
endelse 

; set up params structure to hold everything
if (not double) then begin  ; single precision, plot files
  params = {totBlocks:totBlocks, $
            corners:corners, $
            ndim:ndim, $
            nvar:nvar, $
            nxb:nxb, $
            nyb:nyb, $
            nzb:nzb, $
            ntopx:1, $
            ntopy:1, $
            ntopz:1, $
            time:float(time), $
            dt: float(dt), $
            redshift:1.0, $
            geometry:"unknown"}
endif else begin   ; double precision, checkpoint files
  params = {totBlocks:totBlocks, $
            corners:corners, $
            ndim:ndim, $
            nvar:nvar, $
            nxb:nxb, $
            nyb:nyb, $
            nzb:nzb, $
            ntopx:1, $
            ntopy:1, $
            ntopz:1, $
            time:time, $
            dt:dt, $
            redshift:1.d0, $
            geometry:"unknown"}
endelse   ; end of double/single blocks
params.geometry = strcompress(geometry, /REMOVE_ALL)

; HACK FLASH2
; for some reason, the redshift field is not stored in the plotfiles
if (FLASH2) then begin 
  if ( (where(tag_names(sim_params) EQ "REDSHIFT"))[0] NE -1) then begin
    if (not double) then begin
      params.redshift = float(sim_params.redshift)
    endif else begin
      params.redshift = sim_params.redshift
    endelse
  endif                       
endif  ; end of flash2 redshift hack

end  ; of read_parameters.pro    
