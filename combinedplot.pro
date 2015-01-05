;It would be nice if there was an extra flag to have everything save
;to Downloads or some other drive, in case of stupid disk issues

PRO txttosav, openfile
;Convert txt files to IDL binary 

  ;Check if .sav file already exists before proceeding
  if(~file_test(openfile+'.sav')) then begin
     print, 'Now creating .sav for ',openfile
     
     ;Files to open
     posfile = openfile+'_pos.txt'
     velfile = openfile+'_velocity.txt'
     densfile = openfile+'_density.txt'
     pressurefile = openfile+'_pressure.txt'
     scalarfile = openfile+'_specific_scalar_0_.txt'

     ;Reading in files
     READCOL, posfile, format='D D D', xvec, yvec, zvec
     READCOL, velfile, format='D D D', vxvec, vyvec, vzvec
     READCOL, densfile, format='D', densvec
     READCOL, pressurefile, format='D', pvec
     READCOL, scalarfile, format='D', dratvec

     ;Converting from 1D to 3D arrays
     n = n_elements(densvec)^(1/3.)
     x = reform(xvec, n, n, n)
     y = reform(yvec, n, n, n)
     z = reform(zvec, n, n, n)
     vx = reform(vxvec, n, n, n)
     vy = reform(vyvec, n, n, n)
     vz = reform(vzvec, n, n, n)
     d = reform(densvec, n, n, n)
     p = reform(pvec, n, n, n)
     drat = reform(dratvec, n, n, n)
     
     save, /variables, filename=openfile+'.sav'
  endif
END

PRO mdptsav, openfile, ydim=ydim
;Create file with only mdpt values
;Maybe can be absorbed by plotting function, but for now separate file created

  restore, openfile+'.sav'

  mdpt = floor(n/2.)

  vmag = sqrt(vx^2. + vy^2. + vz^2.)
  cs = sqrt(p/d)

  if(~keyword_set(ydim)) then begin
     if(~file_test(openfile+'_mdptvals.sav')) then begin
        print, 'Now creating .sav for ',openfile,'_mdptvals.sav'
        xmdpt =  x[*,mdpt,mdpt]
        vmagmdpt = vmag[*,mdpt, mdpt]
        vxmdpt = vx[*,mdpt, mdpt]
        ;; vymdpt = vy[*,mdpt, mdpt]
        ;; vzmdpt = vz[*,mdpt, mdpt]
        csmdpt = cs[*,mdpt, mdpt]
        dmdpt = d[*,mdpt, mdpt]
        pmdpt = p[*,mdpt, mdpt]
        dratmdpt = drat[*,mdpt, mdpt]  
        save, /variables, filename= openfile+'_mdptvals.sav'
     endif

;;   mdpt = floor(n/2.)
;;   offmdpt = floor(n/2.)+1

;;   vmag = sqrt(vx^2. + vy^2. + vz^2.)
;;   cs = sqrt(p/d)

;;   if(~keyword_set(ydim)) then begin
;; ;     if(~file_test(openfile+'_mdptvals.sav')) then begin
;;         print, 'Now creating .sav for ',openfile,'_offmdptvals.sav'
;;         xmdpt1 =  x[*,offmdpt,mdpt]
;;         vmagmdpt1 = vmag[*,offmdpt, mdpt]
;;         vxmdpt1 = vx[*,offmdpt, mdpt]
;;         vymdpt1 = vy[*,offmdpt, mdpt]
;;         vzmdpt1 = vz[*,offmdpt, mdpt]
;;         csmdpt1 = cs[*,offmdpt, mdpt]
;;         dmdpt1 = d[*,offmdpt, mdpt]
;;         pmdpt1 = p[*,offmdpt, mdpt]
;;         dratmdpt1 = drat[*,offmdpt, mdpt]  

;;         xmdpt2 =  x[*,mdpt,offmdpt]
;;         vmagmdpt2 = vmag[*,mdpt, offmdpt]
;;         vxmdpt2 = vx[*,mdpt, offmdpt]
;;         vymdpt2 = vy[*,mdpt, offmdpt]
;;         vzmdpt2 = vz[*,mdpt, offmdpt]
;;         csmdpt2 = cs[*,mdpt, offmdpt]
;;         dmdpt2 = d[*,mdpt, offmdpt]
;;         pmdpt2 = p[*,mdpt, offmdpt]
;;         dratmdpt2 = drat[*,mdpt, offmdpt]  

;;         save, /variables, filename= openfile+'_offmdptvals.sav'
  endif else begin
     if(~file_test(openfile+'y_mdptvals.sav')) then begin
        print, 'Now creating .sav for ',openfile,'y_mdptvals.sav'
        xmdpt =  reform(y[mdpt,*,mdpt],n) ;Changed to accommodate the different dimension
        vmagmdpt = reform(vmag[mdpt,*,mdpt],n)
        vxmdpt =reform( vy[mdpt,*,mdpt],n) ;Changed to accommodate the different dimension
        csmdpt = reform(cs[mdpt,*,mdpt],n)
        dmdpt = reform(d[mdpt,*,mdpt],n)
        pmdpt = reform(p[mdpt,*,mdpt],n)
        dratmdpt = reform(drat[mdpt,*,mdpt],n)
        save, /variables, filename= openfile+'y_mdptvals.sav'
     endif
  endelse

END

PRO combinedplot, fileno, levs=levs, fromsaved=fromsaved, PS=PS
;Your one stop shop for plotting everything you could ever want

  filename = 'ioniz_sphere.'+string(fileno, format='(I04)')

  ;Currently only doing base level
  txttosav, filename
  mdptsav, filename

  if (keyword_set(levs)) then begin
     maxlev = 4
     for ilev = 1, maxlev do begin
        filename = 'ioniz_sphere-lev'+strtrim(ilev,2)+'.'+string(fileno, format='(I04)')
        txttosav, filename
        mdptsav, filename
     endfor
  endif

  if (keyword_set(atmos)) then begin
     atmfile1 = "/Volumes/Four1A/atripathi/hires_atmos_r1/ioniz_sphere.0060"
     atmfile2 = "/Volumes/Four1A/atripathi/hyp_recomb_r5/ioniz_sphere.0411"

     if (keyword_set(ps)) then begin
        atmoscompare, atmfile1, atmfile2, /ps
     endif else begin
        atmoscompare, atmfile1, atmfile2
     endelse
     stop
  endif

  if (keyword_set(ps)) then begin
     velplot_allresolutions, fileno, /ps
  endif else begin
     velplot_allresolutions, fileno
  endelse

  print, 'Did I make it?'

;;    rp = 1.5d10


;; ;Times to consider for mass loss case
;;   timesteps = [0,20,47,60,80,112]


;;    ;Do all file I/O first
;;    for tt = 0, n_elements(timesteps) -1 do begin
;;       openfile = basename+'.'+strn(timesteps[tt], length=4, padtype=1, padch='0')


;;    endfor
  stop
END
