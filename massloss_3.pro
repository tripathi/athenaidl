PRO massloss_3, basename, fromsaved=fromsaved, ps=ps, readonly=readonly, justplot=justplot

;Plotting setup
  plotsym, 0, .9, /fill
  set_plot, 'x'
  if (keyword_set(ps)) then begin
     !p.font=0
     !p.charsize=1.7
!y.charsize = 1.1
!x.charsize = 1.1
     !p.charthick=5
     !p.thick=9
     !x.thick=6
     !y.thick=6
     set_plot, 'ps'
  endif
!p.multi=0


if (~keyword_set(justplot)) then begin

   timesteps = [20, 111, 211, 311, 411, 611, 811, 1011, 1211, 1237, 1247, 1257, 1267, 1578, 1687]
   times = (1+ indgen(n_elements(timesteps)-2))*1e5
   times=[times, 1.5e5, 2.5e5]
   times=times[sort(times)]
   timesfine = indgen(200)*1e4
   ;Parameter setup
   rp = 1.5d10
   rreset = 0.75*rp
   rroche = 2.7*rp

   rint = (findgen(3)+3.)*rp
   radii = [rroche, 5.*rp]
;   radii = [rreset, rp, rroche, rint] ;Be sure to update!!

   mshell = fltarr(n_elements(radii), n_elements(timesteps))
   mshellave = fltarr(n_elements(radii), n_elements(timesteps))
   mbox = fltarr(n_elements(radii), n_elements(timesteps))

   tottimes = n_elements(timesteps)




   for ts = 0, tottimes -1 do begin

      openfile = basename+'.'+strn(timesteps[ts], length=4, padtype=1, padch='0')
      restore, openfile+'.sav'

      rad2 = fltarr(n, n, n)
      rad2 = (x^2. + y^2. + z^2.)
      rad = sqrt(rad2)
      dx = abs(xvec[1]-xvec[0])
      dA = dx^2.

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
         xboxp = where (abs(x - radii[r]) lt dx AND (abs(y) lt radii[r]) AND (abs(z) lt radii[r]) AND (x lt radii[r]))
         xboxn = where (abs(x + radii[r]) lt dx AND (abs(y) lt radii[r]) AND (abs(z) lt radii[r]) AND (x gt -radii[r]))
         yboxp = where (abs(y - radii[r]) lt dx AND (abs(x) lt radii[r]) AND (abs(z) lt radii[r]) AND (y lt radii[r]))
         yboxn = where (abs(y + radii[r]) lt dx AND (abs(x) lt radii[r]) AND (abs(z) lt radii[r]) AND (y gt -radii[r]))
         zboxp = where (abs(z - radii[r]) lt dx AND (abs(y) lt radii[r]) AND (abs(x) lt radii[r]) AND (z lt radii[r]))
         zboxn = where (abs(z + radii[r]) lt dx AND (abs(y) lt radii[r]) AND (abs(x) lt radii[r]) AND (z gt -radii[r]))


         print, n_elements(xboxp), n_elements(yboxp), n_elements(zboxp)

         mbox[r,ts] = total(mfluxx(xboxp)*dA) + total(mfluxy(yboxp)*dA) + total(mfluxz(zboxp)*dA) + total(mfluxx(xboxn)*dA) + total(mfluxy(yboxn)*dA) + total(mfluxz(zboxn)*dA)

         mshell[r, ts]= total(mflux[shell]*dA) ;If using costheta, dA needs indices
         mshellave[r, ts] = mean(mflux[shell]) * (4.*!dpi*(radii[r])^2.)
      endfor
   endfor
   

   ;; if(keyword_set(ps)) then begin
   ;;    device, /close
   ;; endif

endif else begin
restore, '/Users/anjalitripathi//Users/anjalitripathi/Downloads/recomb_'+basename+'.sav'
endelse

   if (keyword_set(ps)) then begin
      set_plot, 'ps'
;      !y.charsize = 0.8
      device, filen='/Users/anjalitripathi/Downloads/recomb_massloss_rates.ps', xsize=10.5, ysize=7.5,  /inches,  /landscape;, /encapsulated,
;xsize=10.7, ysize=7.5, /inches, 
;   !p.POSITION=[.1,.1,.9,.9] 
   endif else begin
      wset, 0
   endelse

if (~keyword_set(justplot)) then begin
timesteps=times

endif
print, times      

   pflag = 0

   colors=["blu5", "grn7", "red3"] 
;   colors=[ "red", "orange", "goldenrod", "green", "blue", "purple"]

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
         plot, times, mshell[rr,*], xstyle = 2, ystyle = 10, xtitle='Time [s]', ytitle='Instantaneous mass loss [g/s]', /nodata, xmargin=[10,10] ;margin=[0.05, 0.05, 0.1, 0.05], ;, title='Radius = '+ string(radii[rr]/rp, format=' (F5.3)') + ' Rp',, yra=[0, 2.5e11]
;         pflag = 1
;      endif
      endif
;      oplot, times, mshell[rr,*], psym = 6, color=fsc_color(colors[rr])
;      oplot, times, mshell[rr,*], linestyle=0,
;      color=fsc_color(colors[rr])
;plotsym, 0, /fill
      oplot, times, mshellave[rr,*], psym = 4, color=fsc_color(colors[rr], /brewer);, linestyle=1
      oplot, times, mshellave[rr,*], color=fsc_color(colors[rr], /brewer);linestyle=1, 
;      oplot, times, mbox[rr,*], psym = 2, color=fsc_color(colors[rr], /brewer)
;      oplot, times, mbox[rr,*], linestyle = 2, color=fsc_color(colors[rr], /brewer)
;      al_legend,['1D estimate','Sphere','Mean', 'Box'],linestyle = [5, 0, 1, 2], /right_legend, /bottom, box=0, charsize=1.2
;      al_legend,['','', ''],psym=[6, 5, 2], /right_legend, /bottom, box=0, charsize=1.2
      al_legend,['2.7 Rp', '5.0 Rp', '1D estimate', 'Flux'],linestyle = [0,0,5,0], color=[fsc_color(colors[0], /brewer), fsc_color(colors[1], /brewer), fsc_color("gray"), fsc_color(colors[2], /brewer)],/right_legend, /bottom, box=0;, charsize=1.2
      al_legend,['', '', '', ''],color=[fsc_color(colors[0], /brewer), fsc_color(colors[1], /brewer), fsc_color("gray"), fsc_color(colors[2], /brewer)],psym=[4,4,3,3], /right_legend, /bottom, box=0;, charsize=1.3
;      al_legend,['','', ''],psym=[ 3, 5, 3], /right_legend, /bottom, box=0, charsize=1.2

;      xyouts, 1.65e6, 4e10-rr*1e10, 'Radius = '+ string(radii[rr]/rp, format=' (F3.1)') + ' Rp', color=fsc_color(colors[rr], /brewer), /data, charsize=1.2

      ruthnear = where(abs(radii[rr] - rruth) eq min(abs(radii[rr]-rruth)))
      mruth = 4. * !dpi * rruth[ruthnear]^2. * druth[ruthnear] * vruth[ruthnear]
      print, 'At', radii[rr], ' Ruth nearest pt is ', rruth[ruthnear], 'Massloss', mruth
;      oplot, [min(times)-1, max(times)+1], .26*[mruth, mruth], linestyle=5, color=fsc_color("gray")
;      oplot, [min(times)-1, max(times)+1], .31*[mruth, mruth], linestyle=5, color=fsc_color("gray")
      print, rr
   endfor
   
   axis, yaxis = 1, yrange= [0.01, 13.6], /save, ytitle = 'Flux/F!d0!n', ystyle=1;charsize = 1.3,
   oplot, timesfine,  5*(Erf((timesfine - 1.2d5)/(8d4)) + 1) + 0.1, color=fsc_color(colors[2], /brewer);, linestyle=1

   if(keyword_set(ps)) then begin
      device, /close
   endif

if (~keyword_set(justplot)) then begin
;   times = itimes
endif
   save, /variables, filename='/Users/anjalitripathi/Downloads/recomb_'+basename+'.sav'
  STOP
END

