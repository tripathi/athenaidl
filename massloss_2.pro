PRO massloss_2, basename, fromsaved=fromsaved, ps=ps

;Plotting setup
  plotsym, 0, .9, /fill
  set_plot, 'x'
  if (keyword_set(ps)) then begin
     !p.font=0
     !p.charsize=1.5
     !p.charthick=5
     !p.thick=8
     !x.thick=5
     !y.thick=5
     set_plot, 'ps'
;     device, filen='/Users/anjalitripathi/revised_mu_hires_massloss_geom.eps', xsize=10.5, ysize=7, /inches, /encapsulated
  endif else begin
     window, 0
     window, 1
     wset, 1
  endelse
!p.multi=0
;  !p.multi=[0,2,1]
;  !p.multi=[0,3,2]

   timesteps = [0, 10, 20, 30, 40, 50, 100, 150, 200]

   ;Do all file I/O first
   for tt = 0, n_elements(timesteps) -1 do begin
      openfile = basename+'.'+strn(timesteps[tt], length=4, padtype=1, padch='0')
      if (~keyword_set(fromsaved)) then begin
         densfile = openfile+'_density.txt'
         posfile = openfile+'_pos.txt'
         velfile = openfile+'_velocity.txt'
         scalarfile = openfile+'_specific_scalar_0_.txt'
         
         READCOL, posfile, format='D D D', xvec, yvec, zvec
         READCOL, velfile, format='D D D', vxvec, vyvec, vzvec
         READCOL, densfile, format='D', densvec
         READCOL, scalarfile, format='D', dratvec

         n = n_elements(densvec)^(1/3.)
         x = reform(xvec, n, n, n)
         y = reform(yvec, n, n, n)
         z = reform(zvec, n, n, n)
         vx = reform(vxvec, n, n, n)
         vy = reform(vyvec, n, n, n)
         vz = reform(vzvec, n, n, n)
         d = reform(densvec, n, n, n)
         drat = reform(dratvec, n, n, n)

         save, /variables, filename=openfile+'.sav'
      endif
   endfor

   ;Parameter setup
   rp = 1.5d10
   rreset = 0.75*rp
   rroche = 2.7*rp

   rint = (findgen(3)+3.)*rp
   radii = [rroche, 5.*rp]
;   radii = [rreset, rp, rroche, rint] ;Be sure to update!!

   dt = 5d3                     ;Be sure to update!!
   times = timesteps * dt
   mshell = fltarr(n_elements(radii), n_elements(times))
   mshellave = fltarr(n_elements(radii), n_elements(times))
   mbox = fltarr(n_elements(radii), n_elements(times))

   tottimes = n_elements(times)

   for tt = 0, tottimes -1 do begin
      openfile = basename+'.'+strn(timesteps[tt], length=4, padtype=1, padch='0')
   endfor

   for ts = 0, tottimes -1 do begin
      openfile = basename+'.'+strn(timesteps[ts], length=4, padtype=1, padch='0')
      restore, openfile+'.sav'
      radii = [rroche, 5.*rp]
      timesteps = [0, 10, 20, 30, 40, 50, 100, 150, 200]

      rad2 = fltarr(n, n, n)
      rad2 = (x^2. + y^2. + z^2.)
      rad = sqrt(rad2)
      dx = abs(xvec[1]-xvec[0])
      dA = dx^2.

;      radii[0]=rreset+dx

      vdotr = 0
      vdotr = (vx*x + vy*y + vz*z)/sqrt(rad2)
      vdotx = (vx)*x/abs(x)
      vdoty = (vy)*y/abs(y)
      vdotz = (vz)*z/abs(z)

      mflux = 0
      mflux = d*vdotr           ;Instantaneous mass flux
      mfluxx = d*vdotx
      mfluxy = d*vdoty
      mfluxz = d*vdotz

      for r=0,n_elements(radii)-1 do begin
         shell = where (rad gt radii[r]-dx/2. AND rad le radii[r]+dx/2.)
         ;; xbox = where (abs(x - radii[r]) eq min(abs(x - radii[r])) AND (abs(y) lt radii[r]) AND (abs(z) lt radii[r]))
         ;; ybox = where (abs(y - radii[r]) eq min(abs(y - radii[r])) AND (abs(x) lt radii[r]) AND (abs(z) lt radii[r]))
         ;; zbox = where (abs(z - radii[r]) eq min(abs(z - radii[r])) AND (abs(y) lt radii[r]) AND (abs(x) lt radii[r]))
;         xboxp = where (abs(x - radii[r]) eq min(abs(x - radii[r])) AND (abs(y) lt radii[r]) AND (abs(z) lt radii[r]) AND (x gt 0))

         xboxp = where (abs(x - radii[r]) lt dx AND (abs(y) lt radii[r]) AND (abs(z) lt radii[r]) AND (x lt radii[r]))
         xboxn = where (abs(x + radii[r]) lt dx AND (abs(y) lt radii[r]) AND (abs(z) lt radii[r]) AND (x gt -radii[r]))
         yboxp = where (abs(y - radii[r]) lt dx AND (abs(x) lt radii[r]) AND (abs(z) lt radii[r]) AND (y lt radii[r]))
         yboxn = where (abs(y + radii[r]) lt dx AND (abs(x) lt radii[r]) AND (abs(z) lt radii[r]) AND (y gt -radii[r]))
         zboxp = where (abs(z - radii[r]) lt dx AND (abs(y) lt radii[r]) AND (abs(x) lt radii[r]) AND (z lt radii[r]))
         zboxn = where (abs(z + radii[r]) lt dx AND (abs(y) lt radii[r]) AND (abs(x) lt radii[r]) AND (z gt -radii[r]))

;         print, 'Radius ', radii[r], ' ', min(x[xboxn]), max(x[xboxn])

         ;; if (ts lt 1) then begin
         ;;    plotsym, 0, .3, /fill
         ;;    plot, x/rp, y/rp, /nodata, xtitle = 'x [Rp]', ytitle = 'y [Rp]', title='Radius = '+ string(radii[r]/rp, format=' (F5.3)') + ' Rp'
         ;;    oplot, x[xboxn]/rp, y[xboxn]/rp, color=fsc_color("blue"), psym = 8
         ;;    oplot, x[xboxp]/rp, y[xboxp]/rp, color=fsc_color("cyan"), psym = 8
         ;;    oplot, x[yboxn]/rp, y[yboxn]/rp, color=fsc_color("red")    , psym = 8
         ;;    oplot, x[yboxp]/rp, y[yboxp]/rp, color=fsc_color("magenta") , psym = 8
         ;;    oplot, x[shell]/rp, y[shell]/rp, color=fsc_color("slate gray"), psym =8
         ;; endif

;         box = where( (abs(x - radii[r]) eq min(abs(x-radii[r]))) AND (abs(y - radii[r]) eq min(abs(y-radii[r]))) AND (abs(z - radii[r]) eq min(abs(z-radii[r]))))

         print, n_elements(xboxp), n_elements(yboxp), n_elements(zboxp)

         mbox[r,ts] = total(mfluxx(xboxp)*dA) + total(mfluxy(yboxp)*dA) + total(mfluxz(zboxp)*dA) + total(mfluxx(xboxn)*dA) + total(mfluxy(yboxn)*dA) + total(mfluxz(zboxn)*dA)

         mshell[r, ts]= total(mflux[shell]*dA) ;If using costheta, dA needs indices
         mshellave[r, ts] = mean(mflux[shell]) * (4.*!dpi*(radii[r])^2.)
      endfor
   endfor
   

   ;; if(keyword_set(ps)) then begin
   ;;    device, /close
   ;; endif

   if (keyword_set(ps)) then begin
      set_plot, 'ps'
      device, filen='/Users/anjalitripathi/revised_mu_hires_massloss_rates_2panel.eps', xsize=10.5, ysize=7.5, /inches, /encapsulated
;   !p.POSITION=[.1,.1,.9,.9] 
   endif else begin
      wset, 0
   endelse

   radii = [rroche, 5.*rp]
   pflag = 0
   colors=[ "red", "orange", "goldenrod", "green", "blue", "purple"]

   restore, '/Users/anjalitripathi/Atmospheric-Athena/bin/wind_correctmu.sav' ; Variables rr, rd, rv, rT, rnf, rtau
;Columns:  R(10^10 cm)  rho(10^-15 g/cm^3)  v(10^6 cm/s)  T(10^4 K)
;fp  tau
   rruth = rr*1d10
   druth = rd*1d-15
   vruth = rv*1d6

   for rr=0,n_elements(radii)-1 do begin
;      if (pflag eq 0) then begin
      plotsym, 0, .9, /fill
      if (rr eq 0) then begin
         plot, times, mshell[rr,*], xstyle = 2, ystyle = 2, xtitle='Time [s]', ytitle='Instantaneous mass loss [g/s]', /nodata, yra=[0, 2.5e11];, title='Radius = '+ string(radii[rr]/rp, format=' (F5.3)') + ' Rp'
;         pflag = 1
;      endif
      endif
      oplot, times, mshell[rr,*], psym = 6, color=fsc_color(colors[rr])
      oplot, times, mshell[rr,*], linestyle=0, color=fsc_color(colors[rr])
      oplot, times, mshellave[rr,*], psym = 5, color=fsc_color(colors[rr]), linestyle=1
      oplot, times, mshellave[rr,*], linestyle=1, color=fsc_color(colors[rr])
      oplot, times, mbox[rr,*], psym = 2, color=fsc_color(colors[rr])
      oplot, times, mbox[rr,*], linestyle = 2, color=fsc_color(colors[rr])
      al_legend,['1D estimate','Sphere','Mean', 'Box'],linestyle = [5, 0, 1, 2], /right_legend, /bottom, box=0, charsize=1.2
      al_legend,['','', ''],psym=[6, 5, 2], /right_legend, /bottom, box=0, charsize=1.2
      xyouts, 7e5, 6e10-rr*1e10, 'Radius = '+ string(radii[rr]/rp, format=' (F3.1)') + ' Rp', color=fsc_color(colors[rr]), /data

      ruthnear = where(abs(radii[rr] - rruth) eq min(abs(radii[rr]-rruth)))
      mruth = 4. * !dpi * rruth[ruthnear]^2. * druth[ruthnear] * vruth[ruthnear]
      print, 'At', radii[rr], ' Ruth nearest pt is ', rruth[ruthnear], 'Massloss', mruth
      oplot, [min(times)-1, max(times)+1], .26*[mruth, mruth], linestyle=5, color=fsc_color("gray")
      oplot, [min(times)-1, max(times)+1], .31*[mruth, mruth], linestyle=5, color=fsc_color("gray")
      print, rr
   endfor
   
   
   ;; inreset = where (rad2 le rreset^2.)



   ;; face1= where (x eq min(x))
   ;; face2 = where(x eq max(x))
   ;; face3 = where (y eq min(y))
   ;; face4 = where (y eq max(y))
   ;; face5 = where (z eq min(z))
   ;; face6 = where (z eq max(z))

   ;; mf1 = total(mshellx[face1]*dA)
   ;; mf2 = total(mshellx[face2]*dA)
   ;; mf3 = total(mshelly[face3]*dA)
   ;; mf4 = total(mshelly[face4]*dA)
   ;; mf5 = total(mshellz[face5]*dA)
   ;; mf6 = total(mshellz[face6]*dA)

   ;; f4rpx = where (abs(x)-4.*rp eq min(abs(x)-4.*rp))
   ;; f4rpy = where (abs(y)-4.*rp eq min(abs(y)-4.*rp))
   ;; f4rpz = where (abs(z)-4.*rp eq min(abs(z)-4.*rp))

   ;; mf4rp = total([mshellx(f4rpx),mshelly(f4rpy),mshellz(f4rpz)]*dA) 


   ;; totmass = total(d *(dx^3.))
   ;; resetmass = total(d[inreset]*(dx^3.))
   ;; resetdxmass = total(d[shellreset]*(dx^3.))


;;    dratvec = dratvec-9.9d-5
;; print, min(dratvec)
;;    drat = reform(dratvec, n, n, n)
;;    print, min(drat)

   if(keyword_set(ps)) then begin
      device, /close
   endif

   save, /variables, filename='/Users/anjalitripathi/revised_mu_hiresmassloss_2_'+basename+'.sav'
  STOP
END

