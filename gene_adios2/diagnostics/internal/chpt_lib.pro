;##########################################################################
;# this library contains almost all checkpoint related internal functions #
;##########################################################################

FUNCTION chpt_renorm, var, prec_single=prec_single, Lref=Lref

  COMMON global_vars

  IF NOT KEYWORD_SET(Lref) THEN Lref = series.Lref

  rho_ref = SQRT(series.mref*series.Tref/series.Qref) / series.Bref
  v_ref = SQRT(series.Tref*series.Qref/series.mref)

  return_value = 1.0D / (series.nref / v_ref^3 * rho_ref / Lref)
  IF KEYWORD_SET(prec_single) THEN return_value = FLOAT(return_value)

  RETURN, return_value

END

;#############################################################################
FUNCTION get_chpt_filename, run
  COMMON global_vars

  ;get checkpoint file path *including extension*
  file_chpt = get_file_path('checkpoint',run)
  nruns = N_ELEMENTS(file_chpt.path)

  ;if all found check for empty checkpoints
  ;and try s_checkpoint instead
  IF file_chpt.exist THEN BEGIN
     FOR irun = 0, nruns - 1 DO BEGIN
        IF ((FILE_INFO(file_chpt.path[irun])).SIZE LE 6) THEN BEGIN
           temp_struct = get_file_path('s_checkpoint',run)
           IF (temp_struct.exist) THEN BEGIN
              printerror, 'no valid checkpoint found; using s_checkpoint instead'
              file_chpt.path[irun] = temp_struct.path[0]
           ENDIF
        ENDIF
     ENDFOR
  ENDIF ELSE BEGIN
    ;if no checkpoints have been found, try without extension
    ;which is default for non-scan runs
     file_chpt.exist = 1
     FOR irun = 0, nruns - 1 DO BEGIN
        IF NOT FILE_TEST(file_chpt.path[irun]) THEN BEGIN
           file_chpt.path[irun] = gui.out.data_path + $
                  'checkpoint' + (par.chpt_h5 ? '.h5' : '')
           IF ((FILE_INFO(file_chpt.path[irun])).SIZE LE 6) THEN BEGIN
                                ;try s_checkpoint instead
              file_chpt.path[irun] = gui.out.data_path + $
                  's_checkpoint' + (par.chpt_h5 ? '.h5' : '')
              IF FILE_TEST(file_chpt.path[irun]) THEN $
                 printerror, 'no valid checkpoint found; using s_checkpoint instead'
           ENDIF
           file_chpt.exist = file_chpt.exist AND ((FILE_INFO(file_chpt.path[irun])).SIZE GT 6)
           IF par.chpt_h5 THEN file_chpt.exist = file_chpt.exist AND $
                 H5F_IS_HDF5(file_chpt.path[irun])
        ENDIF
     ENDFOR
  ENDELSE

  RETURN, file_chpt

END

;#############################################################################
FUNCTION get_chpt_time, run, get_n_steps=get_n_steps, first=first, last=last

  COMMON global_vars

  file_chpt = get_chpt_filename(run)

  IF NOT file_chpt.exist THEN BEGIN
     printerror, 'no valid checkpoint file found'
     RETURN, 0.0
  ENDIF

  IF par.chpt_h5 THEN BEGIN
     err = 1 - H5F_IS_HDF5(file_chpt.path)
     IF NOT err THEN chpt_lun = H5F_OPEN(file_chpt.path)
  ENDIF ELSE OPENR, chpt_lun, file_chpt.path, /GET_LUN, $
     ERROR=err, SWAP_ENDIAN=series.swap_endian

  IF err NE 0 THEN BEGIN
     printerror, 'no checkpoint file found'
     RETURN, 0.0
  ENDIF

  IF par.chpt_h5 THEN BEGIN
    h5_dir_str = 'dist/g_'
    time_dset = H5D_OPEN(chpt_lun,h5_dir_str)
    time_attrib = H5A_OPEN_NAME(time_dset,'time')
    time = H5A_READ(time_attrib)
    H5A_CLOSE, time_attrib
    H5D_CLOSE, time_dset

    H5F_CLOSE, chpt_lun

    RETURN, time * time_renorm(0)
  ENDIF ELSE BEGIN
    chpt_prec = 'abcdef'

    READU, chpt_lun, chpt_prec
    IF !QUIET NE 1 THEN PRINT, 'checkpoint precision: ' + chpt_prec

    time_arr = chpt_prec EQ 'SINGLE' ? [0.0,0.0] : [0.0D,0.0D]

    READU, chpt_lun, time_arr
    FREE_LUN, chpt_lun

    RETURN, time_arr[0] * time_renorm(0)
  ENDELSE

END

;##########################################################################

PRO jump_to_chpt_step, jump_step, vsp_lun

END

;##########################################################################

PRO read_chpt_step, vsp_lun

END

;##########################################################################

PRO chpt_loop, diag0

  COMMON global_vars

  IF PTR_VALID(chpt) THEN PTR_FREE, chpt

  run_labels = STRSPLIT(series.run_labels,',',/EXTRACT)

  file_chpt = get_chpt_filename(run)
  IF NOT file_chpt.exist THEN BEGIN
     printerror, 'Skipping chpt loop due to missing checkpoint file'
     RETURN
  ENDIF

  chpt_prec_str = 'abcdef'

  FOR run = 0, N_ELEMENTS(run_labels) - 1 DO BEGIN
    IF par.chpt_h5 THEN BEGIN
      err = 1 - H5F_IS_HDF5(file_chpt.path[run])
      IF NOT err THEN chpt_lun = H5F_OPEN(file_chpt.path[run])
    ENDIF ELSE OPENR, chpt_lun, file_chpt.path[run], /GET_LUN, $
      ERROR=err, SWAP_ENDIAN=series.swap_endian

    IF err NE 0 THEN BEGIN
      printerror, file_chpt.path, ' does not exist'
      RETURN
    ENDIF ELSE BEGIN
      IF !QUIET NE 1 THEN PRINT, 'reading ', file_chpt.path ELSE $
        PRINT, 'reading run ' + run_labels[run] + ': chpt' + $
        (par.chpt_h5 ? ' (HDF5)' : '')
    ENDELSE

    IF par.chpt_h5 THEN BEGIN
      h5_dir_str = 'dist/g_'
      chpt_dset = H5D_OPEN(chpt_lun,h5_dir_str)

      time_attrib = H5A_OPEN_NAME(chpt_dset,'time')
      chpt_time = H5A_READ(time_attrib)
      H5A_CLOSE, time_attrib

      prec_attrib = H5A_OPEN_NAME(chpt_dset,'precision')
      chpt_prec_str = H5A_READ(prec_attrib)
      H5A_CLOSE, prec_attrib

      chpt_prec_single = chpt_prec_str EQ 'SINGLE'
    ENDIF ELSE BEGIN
      READU, chpt_lun, chpt_prec_str

      IF !QUIET NE 1 THEN PRINT, 'checkpoint precision: ' + chpt_prec_str
      chpt_prec_single = chpt_prec_str EQ 'SINGLE'

      time_arr = chpt_prec_single ? [0.0,0.0] : [0.0D,0.0D]
      READU, chpt_lun, time_arr
      chpt_time = time_arr[0] * time_renorm(0)
    ENDELSE

    res = MACHAR()
    IF (chpt_time GE (gui.out.start_t * (1.0 - res.EPS))) AND $
      (chpt_time LE (gui.out.end_t * (1.0 + res.EPS))) THEN BEGIN

      PRINT, 'time: ', rm0es(chpt_time)

      IF par.chpt_h5 THEN BEGIN
        attrib = H5A_OPEN_NAME(chpt_dset,'nx0')
        ni0 = H5A_READ(attrib)
        H5A_CLOSE, attrib
        attrib = H5A_OPEN_NAME(chpt_dset,'nky0')
        nkj0 = H5A_READ(attrib)
        H5A_CLOSE, attrib
        attrib = H5A_OPEN_NAME(chpt_dset,'nz0')
        nz0 = H5A_READ(attrib)
        H5A_CLOSE, attrib
        attrib = H5A_OPEN_NAME(chpt_dset,'nv0')
        nv0 = H5A_READ(attrib)
        H5A_CLOSE, attrib
        attrib = H5A_OPEN_NAME(chpt_dset,'nw0')
        nw0 = H5A_READ(attrib)
        H5A_CLOSE, attrib
        attrib = H5A_OPEN_NAME(chpt_dset,'n_spec')
        n_spec = H5A_READ(attrib)
        H5A_CLOSE, attrib
      ENDIF ELSE BEGIN
        resolution = LONARR(6)
        READU, chpt_lun, resolution
        ni0    = resolution[0]
        nkj0   = resolution[1]
        nz0    = resolution[2]
        nv0    = resolution[3]
        nw0    = resolution[4]
        n_spec = resolution[5]
      ENDELSE

      IF (ni0 NE (par.y_local?par.nx0:par.ny0)) THEN printerror, $
         'WARNING: inconsistent 1st dim. size in checkpoint and parameters!'
      IF (nkj0 NE (par.y_local?par.nky0:par.nkx0/2)) THEN printerror, $
         'WARNING: inconsistent 2nd dim. size in checkpoint and parameters!'
      IF (nz0 NE par.nz0) THEN printerror, $
         'WARNING: inconsistent nz0 in checkpoint and parameters!'
      IF (nv0 NE par.nv0) THEN printerror, $
         'WARNING: inconsistent nv0 in checkpoint and parameters!'
      IF (nw0 NE par.nw0) THEN printerror, $
         'WARNING: inconsistent nw0 in checkpoint and parameters!'
      IF (n_spec NE par.n_spec) THEN printerror, $
         'WARNING: inconsistent n_spec in checkpoint and parameters!'

;      mem_req = (8LL * ni0 * nkj0 * nz0 * nv0 * nw0 * n_spec) / (1024.0)^3
;      IF NOT chpt_prec_single THEN mem_req *= 2
;      PRINT, 'reading '+ file_chpt.path[run] + ' will require ' + $
;        rm0es(mem_req) + ' GB'

      IF !QUIET NE 1 THEN t_pre_read = SYSTIME(1)

      IF gui.out.n_spec_sel LT n_spec THEN BEGIN
        g1_dist = chpt_prec_single ? $
          COMPLEXARR(ni0,nkj0,nz0,nv0,nw0,gui.out.n_spec_sel,/NOZERO) : $
          DCOMPLEXARR(ni0,nkj0,nz0,nv0,nw0,gui.out.n_spec_sel,/NOZERO)

        IF par.chpt_h5 THEN BEGIN
          step_loc = '/dist/g_/'
          chpt_dset = H5D_OPEN(chpt_lun,step_loc)

          content_reim = H5D_READ(chpt_dset)
          H5D_CLOSE, chpt_dset

          g1_dist[0,0,0,0,0,0] = $
            COMPLEX(content_reim[*,*,*,*,*,*gui.out.spec_select].real,$
            content_reim[*,*,*,*,*,*gui.out.spec_select].imaginary,$
            DOUBLE=(1-par.prec_single))
          content_reim = 0
        ENDIF ELSE BEGIN
          g1_dist_in = chpt_prec_single ? $
            COMPLEXARR(ni0,nkj0,nz0,nv0,nw0,/NOZERO) : $
            DCOMPLEXARR(ni0,nkj0,nz0,nv0,nw0,/NOZERO)

          content_byte_size = 8LL * (2 - chpt_prec_single) * $
            ni0 * nkj0 * nz0 * nv0 * nw0

          isp = 0
          FOR sp = 0, n_spec - 1 DO BEGIN
            IF WHERE(*gui.out.spec_select EQ sp) NE -1 THEN BEGIN
              READU, chpt_lun, g1_dist_in
              g1_dist[0,0,0,0,0,isp] = isp EQ gui.out.n_spec_sel - 1 ? $
                TEMPORARY(g1_dist_in) : g1_dist_in

              isp += 1
            ENDIF ELSE BEGIN
              POINT_LUN, - chpt_lun, position_before
              position_after = position_before + content_byte_size
              POINT_LUN, chpt_lun, position_after
            ENDELSE
          ENDFOR
        ENDELSE
      ENDIF ELSE BEGIN
        g1_dist = chpt_prec_single ? $
          COMPLEXARR(ni0,nkj0,nz0,nv0,nw0,n_spec,/NOZERO) : $
          DCOMPLEXARR(ni0,nkj0,nz0,nv0,nw0,n_spec,/NOZERO)

        IF par.chpt_h5 THEN BEGIN
          step_loc = '/dist/g_/'
          chpt_dset = H5D_OPEN(chpt_lun,step_loc)

          content_reim = H5D_READ(chpt_dset)
          H5D_CLOSE, chpt_dset

          g1_dist[0,0,0,0,0,0] = COMPLEX(content_reim.real,$
            content_reim.imaginary,DOUBLE=(1-par.prec_single))
          content_reim = 0
        ENDIF ELSE READU, chpt_lun, g1_dist
      ENDELSE

      IF !QUIET NE 1 THEN PRINT, 'time to read checkpoint data: ' + $
        rm0es(SYSTIME(1)-t_pre_read) + ' s'

;      IF par.chpt_h5 THEN H5F_CLOSE, chpt_lun ELSE FREE_LUN, chpt_lun

      g1_dist *= chpt_renorm(1,prec_single=chpt_prec_single)
      chpt = PTR_NEW(g1_dist,/NO_COPY)

      call_diags, diag0, 'loop'

      IF PTR_VALID(chpt) THEN PTR_FREE, chpt

      series.step_count += 1
    ENDIF

    IF par.chpt_h5 THEN H5F_CLOSE, chpt_lun ELSE FREE_LUN, chpt_lun
  ENDFOR ; --run

END