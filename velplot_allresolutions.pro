PRO plotdvt, infile, inlev, plotifrac, plotheating
  common PHYSVARS, kb, m_h, mu, alphac, rp, km, sigmaph, egamma
  common PLOTPOS, p1, p2, p3, x1, x2, x3, y1, y2, y3, p4, x4, y4, p5, x5, y5, p6, x6, y6, p7, x7, y7
  common LEVLIMS, bounds, filenum
  restore, infile
  rp = 1.5e10
  print, 'For ', infile
  print, 'Min ', min(xmdpt), ' Max ', max(xmdpt)

  dx = abs(xmdpt[2]-xmdpt[1])
  if (inlev lt 3) then begin
     levlimn = where(xmdpt lt -1.*(bounds[inlev+1]-dx))
     levlimp = where(xmdpt gt (bounds[inlev+1]-dx))
  endif else begin
     levlimn = indgen(n_elements(xmdpt))
     levlimp = indgen(n_elements(xmdpt))
  endelse
  
  dneutral = dratmdpt*dmdpt
  nh = dneutral/m_h                        ;Number density of neutral H (atomic)
  nhp = (dmdpt - dneutral)/m_h             ;Number density of ionized H
  nel = nhp + dmdpt * alphac / (14. * m_h) ;Electron number density
  ifrac = nel/(nh + nhp)                   ;Ionization fraction
  mu_temp1 = (ifrac * m_h/2) + (1.-ifrac)*mu
  temp1 = pmdpt/dmdpt/kb*mu_temp1


  ;Find sonic point on day side
  neg = where( xmdpt lt 0)
  sonic = where (abs(vmagmdpt[neg]-csmdpt[neg]) eq min(abs(vmagmdpt[neg] -csmdpt[neg])))
  
  ;Calculate optical depth
  tau = dblarr(n)
  dx = abs(xmdpt[1]-xmdpt[0])
  for i=0,n-1 do begin   
     tau[i] = total(nh[0:i] *dx *sigmaph) ; - total(d[*, 0, floor(n/2)-1]* dx)
  endfor
   
   ;Find tau=1 location
   tau1 = where(abs(tau - 1) eq min(abs(tau -1)))

if (plotheating lt 1) then begin
   if (plotifrac lt 1) then begin
                                ;OVERPLOT
                                ;Plot 1D density
      !p=p1 & !x = x1 & !y = y1
      oplot, xmdpt[levlimn]/rp, dmdpt[levlimn] ;, color=fsc_color("green")
      oplot, xmdpt[levlimp]/rp, dmdpt[levlimp] ;, color=fsc_color("green")
      ;; plots, xmdpt[tau1]/rp,dmdpt(tau1), psym =7

                                ;Plot velocity
      !p = p2  & !x= x2  & !y = y2
      oplot, xmdpt[levlimn]/rp, vmagmdpt[levlimn]/km ;, color=fsc_color("green")
      oplot, xmdpt[levlimp]/rp, vmagmdpt[levlimp]/km ;, color=fsc_color("green")
      ;; plots, xmdpt[neg[sonic]]/rp, vmagmdpt[neg[sonic]]/km, psym = 8

                                ;Plot log T
      !p = p3  & !x= x3  & !y = y3
      oplot, xmdpt[levlimn]/rp, temp1[levlimn] ;,color=fsc_color("green")
      oplot, xmdpt[levlimp]/rp, temp1[levlimp] ;,color=fsc_color("green")
      ;; plots, xmdpt[neg[sonic]]/rp, temp1[neg[sonic]], psym = 8

   endif else begin
      !p = p4  & !x= x4  & !y = y4
;     levlimn = where(xmdpt lt -1.*(bounds[inlev+1]-dx))
;     levlimp = where(xmdpt gt (bounds[inlev+1]-dx))

      oplot, xmdpt[levlimn]/rp, ifrac[levlimn];, color=fsc_color("green")
      oplot, xmdpt[levlimp]/rp, ifrac[levlimp];, color=fsc_color("green")
      print,'Im in here'
   endelse
endif else begin
      ;Photoionization
      flux0 = 10.1 * 2d13              ;Accounting for the 10.1 factor
      flux = flux0 * exp(-tau)
      rate_ioniz = flux * nh * sigmaph
      rate_ionizheat = rate_ioniz * egamma

      ;Recombination
      alpharecomb = 2.59d-13*(temp1/1.0d4)^(-0.7)
      rate_recomb = -alpharecomb * nel * nhp

      ;Recomb cooling
      rate_recombcool = -6.11d-10 * temp1^(-0.89) * kB * temp1 * nel * nhp

      ;Lya cooling
      rate_lyacool = -7.5d-19 * exp (-118348/temp1) * nel * nh

      ;New advection heating
      dx = xmdpt[2]-xmdpt[1]

      ;X DIRECTION TERMS
      ;Internal energy contribution
      epsx = pmdpt/dmdpt/(5./3. - 1.) ;Specific internal energy
      depsx = shift(epsx,-1) - epsx    ;Finite difference for derivative
      advheatlhsx = dmdpt * vxmdpt * depsx/dx

      ;PdV contribution - Following sign convention of Ruth -
      ;i.e. making these negative
      advheatrhsx = pmdpt/dmdpt * vxmdpt * (shift(dmdpt,-1) - dmdpt)/dx

      ;Advection effect on ionization
      advionizx = dmdpt/m_h * vxmdpt * (shift(ifrac,-1) - ifrac)/dx


      ;Y DIRECTION TERMS
;      restore, infile
      filen = 'ioniz_sphere-lev'+strtrim(inlev,2)+'.'+string(filenum, format='(I04)')+'_offmdptvals.sav'
      restore, filen
      rp = 1.5e10

      epsy = pmdpt1/dmdpt1/(5./3. - 1.) ;Specific internal energy
      depsy = epsy-epsx
      advheatlhsy = dmdpt * vymdpt * depsy/dx ;Internal energy contribution
      advheatrhsy = pmdpt/dmdpt * vymdpt * (dmdpt1-dmdpt)/dx ;PdV contribution - Following sign convention of Ruth

      ;Advection effect on ionization
      dneutral1 = dratmdpt1*dmdpt1
      nh1 = dneutral1/m_h                        ;Number density of neutral H (atomic)
      nhp1 = (dmdpt1 - dneutral1)/m_h            ;Number density of ionized H
      nel1 = nhp1 + dmdpt1 * alphac / (14. * m_h) ;Electron number density
      ifrac1 = nel1/(nh1 + nhp1)                  ;Ionization fraction
      advionizy = dmdpt/m_h * vymdpt * (ifrac1-ifrac)/dx

      ; Z DIRECTION TERMS

      epsz = pmdpt2/dmdpt2/(5./3. - 1.) ;Specific internal energy
      depsz = epsz -epsx
      advheatlhsz = dmdpt * vzmdpt * depsz/dx
      advheatrhsz = pmdpt/dmdpt * vzmdpt * (dmdpt2-dmdpt)/dx

      dneutral2 = dratmdpt2*dmdpt2
      nh2 = dneutral2/m_h                        ;Number density of neutral H (atomic)
      nhp2 = (dmdpt2 - dneutral2)/m_h            ;Number density of ionized H
      nel2 = nhp2 + dmdpt2 * alphac / (14. * m_h) ;Electron number density
      ifrac2 = nel2/(nh2 + nhp2)                  ;Ionization fraction
      advionizz = dmdpt/m_h * vzmdpt * (ifrac2-ifrac)/dx

      ;Remove last element of array, since there are edge effects
      remove, uint(n-1), advionizx
      remove, uint(n-1), advionizy
      remove, uint(n-1), advionizz
      remove, uint(n-1), advheatlhsx
      remove, uint(n-1), advheatlhsy
      remove, uint(n-1), advheatlhsz
      remove, uint(n-1), advheatrhsx
      remove, uint(n-1), advheatrhsy
      remove, uint(n-1), advheatrhsz


      ;Advection totals
      advioniz = advionizx + advionizy + advionizz
      advheatlhs = advheatlhsx+ advheatlhsy + advheatlhsz
      advheatrhs = advheatrhsx+ advheatrhsy + advheatrhsz

      ;Range for each level, excluding gray area
      ilevmax = 3
      if inlev lt ilevmax then begin
         alimamt = max([(bounds[inlev+1]-dx), rp]) ;;where(xmdpt/rp gt -20 AND xmdpt/rp lt 20)
      endif else begin
         alimamt = rp ;;where(xmdpt/rp gt -20 AND xmdpt/rp lt 20)
      endelse
      alevlimn = where(xmdpt lt -1.*alimamt)
      alevlimp = where(xmdpt gt alimamt)

      !p = p6  & !x= x6  & !y = y6
      oplot, xmdpt[alevlimn]/rp, rate_ionizheat[alevlimn], color=fsc_color("red")
      oplot, xmdpt[alevlimp]/rp, rate_ionizheat[alevlimp], color=fsc_color("red")
      oplot, xmdpt[alevlimn]/rp, (rate_recombcool[alevlimn]), color=fsc_color("blue")
      oplot, xmdpt[alevlimp]/rp, (rate_recombcool[alevlimp]), color=fsc_color("blue")
      oplot, xmdpt[alevlimn]/rp, (rate_lyacool[alevlimn]), color=fsc_color("teal")
      oplot, xmdpt[alevlimp]/rp, (rate_lyacool[alevlimp]), color=fsc_color("teal")
      ;; oplot, xmdpt/rp, abs(advheatx), color=fsc_color("brown")
      ;; oplot, xmdpt/rp, abs(advheaty), color=fsc_color("orange")
      ;; oplot, xmdpt/rp, abs(advheatz), color=fsc_color("yellow")
      ;; oplot, xmdpt/rp, abs(advheatx+advheaty+advheatz);, color=fsc_color("brown")
      ;; oplot, xmdpt/rp, abs(advtermx), color=fsc_color("purple")
      ;; oplot, xmdpt/rp, abs(advtermy), color=fsc_color("green")
      ;; oplot, xmdpt/rp, abs(advtermz), color=fsc_color("magenta")
      ;; oplot, xmdpt/rp, abs(advtermx+advtermy+advtermz), color=fsc_color("cyan")
      oplot, xmdpt[alevlimn]/rp, (advheatlhs[alevlimn]), color=fsc_color("orange")
      oplot, xmdpt[alevlimp]/rp, (advheatlhs[alevlimp]), color=fsc_color("orange")
      oplot, xmdpt[alevlimn]/rp, (advheatrhs[alevlimn]), color=fsc_color("green")
      oplot, xmdpt[alevlimp]/rp, (advheatrhs[alevlimp]), color=fsc_color("green")
       outsidea = alevlimn;where ((xmdpt ) ge rp)
       outsideb = alevlimp;where ((xmdpt ) le rp)

       oplot, xmdpt(outsidea)/rp, (rate_ionizheat(outsidea) + rate_recombcool(outsidea) + rate_lyacool(outsidea) + advheatrhs(outsidea))
       oplot, xmdpt(outsideb)/rp, (rate_ionizheat(outsideb) + rate_recombcool(outsideb) + rate_lyacool(outsideb) + advheatrhs(outsideb))

      print, min(advheatrhs), max(rate_ionizheat)


      !p = p7  & !x= x7  & !y = y7
      oplot, xmdpt[alevlimn]/rp, rate_ioniz[alevlimn], color=fsc_color("red")
      oplot, xmdpt[alevlimp]/rp, rate_ioniz[alevlimp], color=fsc_color("red")
      oplot, xmdpt[alevlimn]/rp, abs(rate_recomb[alevlimn]), color=fsc_color("blue")
      oplot, xmdpt[alevlimp]/rp, abs(rate_recomb[alevlimp]), color=fsc_color("blue")
      oplot, xmdpt[alevlimn]/rp, advioniz[alevlimn], color=fsc_color("orange")
      oplot, xmdpt[alevlimp]/rp, advioniz[alevlimp], color=fsc_color("orange")





endelse

END


PRO velplot_allresolutions, fileno, ps=ps, recombtime=recombtime, yplot=yplot;, fromsaved=fromsaved
common PHYSVARS, kb, m_h, mu, alphac, rp, km, sigmaph, egamma
common PLOTPOS, p1, p2, p3, x1, x2, x3, y1, y2, y3, p4, x4, y4, p5, x5, y5, p6, x6, y6, p7, x7, y7
common LEVLIMS, bounds, filenum

filenum = fileno

;;Currently hardcoded. Can be changed/checked later
bounds=[7.4062500e+10, 2.390625d10, 2.2265625e+10, 2.0507812e+10, 1.8691406e+10]

;Plotting setup
plotsym, 0, .9, /fill
set_plot, 'x'
if (keyword_set(ps)) then begin
    !p.font=0
    !p.charsize=2
    !p.charthick=5
    !p.thick=8
    !x.thick=5
    !y.thick=5
endif


if (~keyword_set(fromsaved)) then begin
   ;; for i=3,9 do begin
   ;;    print, i
   ;;    vel1d, 'ioniz_sphere.0'+strtrim(i,2)+'00'
   ;; endfor
   ;; for i=11,19 do begin
   ;;    vel1d, 'ioniz_sphere.'+strtrim(i,2)+'00'
   ;; endfor
endif


;Physical setup
kb = 1.38e-16
m_h = 1.67d-24
mu = 1.67d-24                   ;Neutral gas mean molec weight
alphac = 1d-10
rp = 1.5d10                     ;CHANGE ACCORDINGLY
km = 1d6
sigmaph = 6.3d-18               ;CHANGE ACCORDINGLY
egamma = 3.85e-12               ;2.4eV in ergs


;Convert Ruth's .dat to a .sav
;readcol, 'tideflux4.5E+03_papermu.dat', format='D,D,D,D,D,D,D,D,D', rr, rd, rv, rT, rnf, rtau , rrq, rrz, rrtau
;save, /variables, filename='tideflux_papermu.sav'  

;Read in Ruth's wind solution
restore, '/Users/anjalitripathi/Atmospheric-Athena/bin/tideflux_papermu.sav'
;restore, '/Users/anjalitripathi/Atmospheric-Athena/bin/wind_correctmu.sav' ; Variables rr, rd, rv, rT, rnf, rtau
                                                                           ;Columns:  R(10^10 cm)  rho(10^-15 g/cm^3)  v(10^6 cm/s)  T(10^4 K)  fp  tau
rx = -rr/(rp/1d10)
rtemp = rT * 1d4
rdens = rd*1e-15
rcs = sqrt(kb/mu * rtemp)       ;Is neutral mu appropriate here?

;Main loop over files
;Adjust ii and string, to get filename formating correct !!!


   !p.multi=[0,1,3]

   if (keyword_set(ps)) then begin
      set_plot, 'ps'
      !p.charsize=1.7
      !y.charsize = 1.7
      !x.charsize = 1.7
      !p.thick=9
      !x.thick=6
      !y.thick=6

      ii = fileno

      if (keyword_set(yplot)) then begin
         device, filen='~/Downloads/velplot_'+string(fileno, format='(I04)')+'y_4lev.eps', ysize=10.5, xsize=8.2, /inches, /encapsulated
      endif else begin
         device, filen='~/Downloads/velplot_'+string(fileno, format='(I04)')+'_4lev.eps', ysize=10.5, xsize=8.2, /inches, /encapsulated
      endelse

;   !p.POSITION=[.1,.1,.9,.9] 
  endif else begin
     window
  endelse

; RESTORE HIRES FILE
  if (keyword_set(yplot)) then begin
     restore, 'ioniz_sphere.'+string(fileno, format='(I04)')+'y_mdptvals.sav'
  endif else begin
     restore, 'ioniz_sphere.'+string(fileno, format='(I04)')+'_mdptvals.sav'
  endelse

  ;Choose plotting bounds for the base level
  dx = abs(xmdpt[2]-xmdpt[1])
  limamt = (bounds[1]-dx) ;;where(xmdpt/rp gt -20 AND xmdpt/rp lt 20)
  levlimn = where(xmdpt lt - limamt)
  levlimp = where(xmdpt gt limamt)

  ;Physical quantities 
  dneutral = dratmdpt*dmdpt
  nh = dneutral/m_h                        ;Number density of neutral H (atomic)
  nhp = (dmdpt - dneutral)/m_h             ;Number density of ionized H
  nel = nhp + dmdpt * alphac / (14. * m_h) ;Electron number density
  ifrac = nel/(nh + nhp)                   ;Ionization fraction

  mu_temp1 = (ifrac * m_h/2) + (1.-ifrac)*mu

  temp1 = pmdpt/dmdpt/kb*mu_temp1


  ;Find sonic point on day side
  neg = where( xmdpt lt 0)
  sonic = where (abs(vmagmdpt[neg]-csmdpt[neg]) eq min(abs(vmagmdpt[neg] -csmdpt[neg])))
  rsonic = where (abs(rv*km - rcs) eq min(abs(rv*km - rcs))) ;For Ruth's output
  
  ;Calculate optical depth
  tau = dblarr(n)
  dx = abs(xmdpt[1]-xmdpt[0])
  for i=0,n-1 do begin   
     tau[i] = total(nh[0:i] *dx *sigmaph) ; - total(d[*, 0, floor(n/2)-1]* dx)
  endfor
   
  ;Find tau=1 location
  tau1 = where(abs(tau - 1) eq min(abs(tau -1)))
  rtau1 = where(abs(rtau - 1) eq min(abs(rtau -1))) ;For Ruth's output

  
  if (keyword_set(recombtime)) then begin
;      window, 2
     !p.multi=0
     !p.charsize = 1.5
     plot, xmdpt/rp, 1./(nhp+1d-12)/3d-13, xsty = 1, ysty=2, /ylog, yra=[1, 1e9], ytitle='Time', xtitle='X/1.5e10 [cm]'
     oplot, xmdpt/rp,(2.7*1.5d10)/sqrt(pmdpt/dmdpt), color=fsc_color("red")
     oplot, xmdpt/rp,(1.5d10)/sqrt(pmdpt/dmdpt), color=fsc_color("blue")
     stop
  endif 

   ;PLOTS

   ;Plot 1D density
   plot, xmdpt/rp, dmdpt, xsty = 1, ysty = 2, xtitle ='x/R!dp!n', ytitle='Density [g cm!e-3!n]', /ylog, yra=[7e-19, 2e-13], /nodata, xmargin=[15,5];, ycharsize = 1.1, xmargin=[5,5];, title='t/1e5s=0'+strtrim(ii-1,2)
   minval = 10.^(!y.crange[0])+1e-20
   maxval = 10.^(!y.crange[1])-1e-13
   polyfill, [-1,1,1,-1], [minval,minval, maxval, maxval], color=fsc_color("light gray"), /data
   oplot, xmdpt[levlimn]/rp,  dmdpt[levlimn]
   oplot, xmdpt[levlimp]/rp,  dmdpt[levlimp]
   oplot, rx, rdens, linestyle = 1
   plots, xmdpt[tau1]/rp,dmdpt(tau1), psym =7
   plots, rx[rtau1], rdens(rtau1), psym =7
   plots, xmdpt[neg[sonic]]/rp, dmdpt[neg[sonic]], psym = 8
   plots, rx[rsonic], rdens[rsonic], psym =8 ;, color=fsc_color("gray")

   ;Save plot coordinates for returning later
   p1 = !P & x1 = !X & y1 = !Y


   ;Plot 1D velocity
   plot, xmdpt/rp, vmagmdpt/km, xsty = 1, ysty= 1, ytitle='Velocity magnitude [10 km s!e-1!n]', xtitle ='x/R!dp!n',  yra=[0, 2.6], /nodata, xmargin=[15,5];ycharsize = 1.1 
   minval = (!y.crange[0])+2e-2
   maxval = (!y.crange[1])-2e-2
   polyfill, [-1,1,1,-1], [minval,minval, maxval, maxval], color=fsc_color("light gray"), /data
   oplot, xmdpt[levlimn]/rp, vmagmdpt[levlimn]/km
   oplot, xmdpt[levlimp]/rp, vmagmdpt[levlimp]/km
   ;oplot, xmdpt/rp, vxmdpt/km, linestyle = 4, thick = 3
   oplot, rx, rv, linestyle=1
   plots, xmdpt[neg[sonic]]/rp, vmagmdpt[neg[sonic]]/km, psym = 8
   plots, rx[rsonic], rv[rsonic], psym =8 ;, color=fsc_color("gray")
   ;oplot, xmdpt/rp, csmdpt/km, color=fsc_color("blue")
   ;oplot, [-.75, -.75], [0, 10], linestyle = 2, color=fsc_color("gray");, thick = 5
   ;oplot, [.75, .75], [0, 10], linestyle = 2, color=fsc_color("gray");, thick = 5
   ;plot, xmdpt/rp, vmagmdpt/csmdpt, xsty = 1, ysty = 1, xtitle ='x/R!dp!n', ytitle='Velocity/Sound speed'

   ;Save plot coordinates for returning later
   p2 = !P & x2 = !X & y2 = !Y


   ;Plot log T
   plot, xmdpt/rp, temp1, xsty = 1, ysty = 1, xtitle ='x/R!dp!n', /ylog, yra=[9e2, 3e4], /nodata, xmargin=[15,5], ymargin=[8,2], ytickformat='(A1)';, ycharsize = 1.1 ;, linestyle = 4, yra=[9e2, 4e4]
   axis, yaxis=0, yrange=[9e2,3e4], /ylog, ysty=1, ytickv=[1e3, 1e4], ytickname=['10!u3!n','10!u4!n'], ytitle='Temperature [K]'
   minval = 10.^(!y.crange[0])+2e1
   maxval = 10.^(!y.crange[1])-1e3
   polyfill, [-1,1,1,-1], [minval,minval, maxval, maxval], color=fsc_color("light gray"), /data
   oplot, xmdpt[levlimn]/rp, temp1[levlimn]
   oplot, xmdpt[levlimp]/rp, temp1[levlimp]
   oplot, rx, rtemp, linestyle=1
   plots, xmdpt[neg[sonic]]/rp, temp1[neg[sonic]], psym = 8
   plots, rx[rsonic], rtemp[rsonic], psym =8 ;, color=fsc_color("gray")

   p3 = !P & x3 = !X & y3 = !Y

;plot, xmdpt/rp, dratmdpt, xsty =1, ysty =1, xtitle='x/R!dp!n', ytitle='Neutral fraction', yra=[0,1]
;oplot, xmdpt/rp, ifrac, linestyle =4

   for ilev = 1, 3 do begin
      if (keyword_set(yplot)) then begin
         filen = 'ioniz_sphere-lev'+strtrim(ilev,2)+'.'+string(fileno, format='(I04)')+'y_mdptvals.sav'
      endif else begin
         filen = 'ioniz_sphere-lev'+strtrim(ilev,2)+'.'+string(fileno, format='(I04)')+'_mdptvals.sav'
      endelse

      plotdvt, filen, ilev, 0, 0
   endfor

   ;RESTORE LOW-RES FILE
   ;; if ((ii-1) mod 3 lt 1) then begin
      ;; !p.multi=0                ;So that the coordinates saved can be applied
      ;; restore, '/Volumes/Extra/revised_mu_0827/ioniz_sphere.0'+strtrim(ii/3,2)+'00_mdptvals.sav'

      ;; dneutral = dratmdpt*dmdpt
      ;; nh = dneutral/m_h                    ;Number density of neutral H (atomic)
      ;; nhp = (dmdpt - dneutral)/m_h         ;Number density of ionized H
      ;; nel = nhp + dmdpt * alphac / (14. * m_h) ;Electron number density
      ;; ifrac = nel/(nh + nhp)                   ;Ionization fraction

      ;; mu_temp1 = (ifrac * m_h/2) + (1.-ifrac)*mu
      ;; temp1 = pmdpt/dmdpt/kb*mu_temp1

      ;; ;Find sonic point on day side
      ;; neg = where(xmdpt lt 0)
      ;; sonic = where (abs(vmagmdpt[neg]-csmdpt[neg]) eq min(abs(vmagmdpt[neg] -csmdpt[neg])))
      ;; rsonic = where (abs(rv*km - rcs) eq min(abs(rv*km - rcs)))
  
      ;; ;Calculate optical depth
      ;; tau = dblarr(n)
      ;; dx = abs(xmdpt[1]-xmdpt[0])
      ;; for i=0,n-1 do begin   
      ;;    tau[i] = total(nh[0:i] *dx *sigmaph) ; - total(d[*, 0, floor(n/2)-1]* dx)
      ;; endfor
      ;; tau1 = where(abs(tau - 1) eq min(abs(tau -1)))
      ;; rtau1 = where(abs(rtau - 1) eq min(abs(rtau -1)))


      ;; ;OVERPLOT
      ;; ;Plot 1D density
      ;; !p=p1 & !x = x1 & !y = y1
      ;; oplot, xmdpt/rp, dmdpt, color=fsc_color("green")
      ;; plots, xmdpt[tau1]/rp,dmdpt(tau1), psym =7

      ;; ;Plot velocity
      ;; !p = p2  & !x= x2  & !y = y2
      ;; oplot, xmdpt/rp, vmagmdpt/km, color=fsc_color("green")
      ;; plots, xmdpt[neg[sonic]]/rp, vmagmdpt[neg[sonic]]/km, psym = 8

      ;; ;Plot log T
      ;; !p = p3  & !x= x3  & !y = y3
      ;; oplot, xmdpt/rp, temp1,color=fsc_color("green")
      ;; plots, xmdpt[neg[sonic]]/rp, temp1[neg[sonic]], psym = 8
   ;; endif


if(keyword_set(ps)) then begin
      device, /close
endif
;stop

print, 'Moving onto ifrac'

if (keyword_set(ps)) then begin
   set_plot, 'ps'
   !p.charsize=1.7
   !y.charsize = 1.7
   !x.charsize = 1.7
   !p.thick=9
   !x.thick=6
   !y.thick=6
   if (keyword_set(yplot)) then begin
      device, filen='~/Downloads/joined_ifrac_'+string(fileno, format='(I04)')+'y_anim.eps', ysize=10.5, xsize=8, /inches, /encapsulated
   endif else begin
      device, filen='~/Downloads/joined_ifrac_'+string(fileno, format='(I04)')+'_anim.eps', ysize=10.5, xsize=8, /inches, /encapsulated
   endelse

;   !p.POSITION=[.1,.1,.9,.9] 
endif
!p.multi=[0,1,2]

;Plot n
;plot, xmdpt/rp, dmdpt/mu_temp1, xsty = 1, ysty = 1, xtitle ='x/R!dp!n', ytitle='Number density [cm!e-3!n]', /ylog;, linestyle = 4, yra=[9e2, 4e4]
;oplot, xmdpt/rp, nh, linestyle = 2
;oplot, xmdpt/rp, nhp, linestyle = 3

;oplot, rx, rdens/mu, color=fsc_color("gray")
;oplot, rx, (1-rnf)*rdens/mu, linestyle = 3, color=fsc_color("gray")
;oplot, rx, rnf*rdens/mu, linestyle = 5, color=fsc_color("gray")

if (keyword_set(yplot)) then begin
   restore, 'ioniz_sphere.'+string(fileno, format='(I04)')+'y_mdptvals.sav'
endif else begin
   restore, 'ioniz_sphere.'+string(fileno, format='(I04)')+'_mdptvals.sav'
endelse


;Plot ifrac
plot, xmdpt/rp, ifrac, xsty =1, ysty =1, xtitle='x/R!dp!n', ytitle='Ionization fraction', yra=[0,1], /nodata, xmargin=[15,5];, title='t/1e5s=0'+strtrim(ii-1,2)
oplot, xmdpt[levlimn]/rp, ifrac[levlimn]
oplot, xmdpt[levlimp]/rp, ifrac[levlimp]
oplot, rx, rnf, linestyle = 1
p4 = !P & x4 = !X & y4 = !Y

plots, xmdpt[tau1]/rp, ifrac[tau1], psym = 7
plots, rx[rtau1], rnf[rtau1], psym = 7;, color=fsc_color("gray")

;Plot tau
plot, xmdpt/rp, tau, xsty = 1, ysty = 1, xtitle ='x/R!dp!n', ytitle='Optical depth', /ylog, yra=[8e-5, 1.5e4], xmargin=[15,5], ymargin=[8,2];
polyfill,  [min(xmdpt/rp)+5d-2,max(xmdpt/rp)-5d-2, max(xmdpt/rp)-5d-2, min(xmdpt/rp)+5d-2], [1., 1., 10.^(!y.crange[1])-40, 10.^(!y.crange[1]) -40 ], color=fsc_color("light gray"), /data
oplot, xmdpt/rp, tau
oplot, rx, rtau, linestyle = 1
p5 = !P & x5 = !X & y5 = !Y

   for ilev = 1, 3 do begin
      if (keyword_set(yplot)) then begin
         filen = 'ioniz_sphere-lev'+strtrim(ilev,2)+'.'+string(fileno, format='(I04)')+'y_mdptvals.sav'
      endif else begin
         filen = 'ioniz_sphere-lev'+strtrim(ilev,2)+'.'+string(fileno, format='(I04)')+'_mdptvals.sav'
      endelse


      plotdvt, filen, ilev, 1, 0
   endfor


;; if ((ii-1) mod 3 lt 1) then begin
;;    !p.multi=0
;;    restore, '/Volumes/Extra/revised_mu_0827/ioniz_sphere.0'+strtrim(ii/3,2)+'00_mdptvals.sav'
;;    !p=p4 & !x = x4 & !y = y4

;;   dneutral = dratmdpt*dmdpt
;;   nh = dneutral/m_h                      ;Number density of neutral H (atomic)
;;   nhp = (dmdpt - dneutral)/m_h           ;Number density of ionized H
;;   nel = nhp + dmdpt * alphac / (14. * m_h) ;Electron number density
;;   ifrac = nel/(nh + nhp)                   ;Ionization fraction

;;   mu_temp1 = (ifrac * m_h/2) + (1.-ifrac)*mu
;; ;mu_temp2 = (ifrac * m_h/2) + (1.-ifrac)*m_h

;;   temp1 = pmdpt/dmdpt/kb*mu_temp1
;; ;temp2 = pmdpt/dmdpt/kb*mu_temp2

;;   rcs = sqrt(kb/mu * rtemp)
  
;;   neg = where( xmdpt lt 0)
;;   sonic = where (abs(vmagmdpt[neg]-csmdpt[neg]) eq min(abs(vmagmdpt[neg] -csmdpt[neg])))
;;   rsonic = where (abs(rv*km - rcs) eq min(abs(rv*km - rcs)))
  
;; ;  rruth = where (xmdpt/rp le max(rx))
;; tau = dblarr(n)
;; dx = abs(xmdpt[1]-xmdpt[0])
;; sigmaph = 6.3d-18
;;    for i=0,n-1 do begin   
;;       tau[i] = total(nh[0:i] *dx *sigmaph); - total(d[*, 0, floor(n/2)-1]* dx)
;;    endfor

;; tau1 = where(abs(tau - 1) eq min(abs(tau -1)))
;; rtau1 = where(abs(rtau - 1) eq min(abs(rtau -1)))


;; ;Plot ifrac
;; oplot, xmdpt/rp, ifrac, color=fsc_color("green")
;; plots, xmdpt[tau1]/rp, ifrac[tau1], psym = 7

;; ;Plot tau
;; !p = p5  & !x= x5  & !y = y5
;; oplot, xmdpt/rp, tau, color=fsc_color("green")

;; endif




;plots, xmdpt[tau1]/rp, tau[tau1], psym = 8
;plots, rx[rtau1], rtau[rtau1], psym = 8, color=fsc_color("gray")

if(keyword_set(ps)) then begin
   device, /close
endif

     restore, 'ioniz_sphere.'+string(fileno, format='(I04)')+'_mdptvals.sav'

  dneutral = dratmdpt*dmdpt
  nh = dneutral/m_h                        ;Number density of neutral H (atomic)
  nhp = (dmdpt - dneutral)/m_h             ;Number density of ionized H
  nel = nhp + dmdpt * alphac / (14. * m_h) ;Electron number density
  ifrac = nel/(nh + nhp)                   ;Ionization fraction

  mu_temp1 = (ifrac * m_h/2) + (1.-ifrac)*mu

  temp1 = pmdpt/dmdpt/kb*mu_temp1


      ;Photoionization
      flux0 = 10.1 * 2d13              ;Accounting for the 10.1 factor
      flux = flux0 * exp(-tau)
      rate_ioniz = flux * nh * sigmaph
      rate_ionizheat = rate_ioniz * egamma

      ;Recombination
      alpharecomb = 2.59d-13*(temp1/1.0d4)^(-0.7)
      rate_recomb = -alpharecomb * nel * nhp

      ;Recomb cooling
      rate_recombcool = -6.11d-10 * temp1^(-0.89) * kB * temp1 * nel * nhp

      ;Lya cooling
      rate_lyacool = -7.5d-19 * exp (-118348/temp1) * nel * nh

      ;New advection heating
      dx = xmdpt[2]-xmdpt[1]

      ;X DIRECTION TERMS
      ;Internal energy contribution
      epsx = pmdpt/dmdpt/(5./3. - 1.) ;Specific internal energy
      depsx = shift(epsx,-1) - epsx    ;Finite difference for derivative
      advheatlhsx = dmdpt * vxmdpt * depsx/dx

      ;PdV contribution - Following sign convention of Ruth -
      ;i.e. making these negative
      advheatrhsx = pmdpt/dmdpt * vxmdpt * (shift(dmdpt,-1) - dmdpt)/dx

      ;Advection effect on ionization
      advionizx = dmdpt/m_h * vxmdpt * (shift(ifrac,-1) - ifrac)/dx


      ;Y DIRECTION TERMS
      restore, 'ioniz_sphere.'+string(fileno, format='(I04)')+'_offmdptvals.sav'

      epsy = pmdpt1/dmdpt1/(5./3. - 1.) ;Specific internal energy
      depsy = epsy-epsx
      advheatlhsy = dmdpt * vymdpt * depsy/dx ;Internal energy contribution
      advheatrhsy = pmdpt/dmdpt * vymdpt * (dmdpt1-dmdpt)/dx ;PdV contribution - Following sign convention of Ruth

      ;Advection effect on ionization
      dneutral1 = dratmdpt1*dmdpt1
      nh1 = dneutral1/m_h                        ;Number density of neutral H (atomic)
      nhp1 = (dmdpt1 - dneutral1)/m_h            ;Number density of ionized H
      nel1 = nhp1 + dmdpt1 * alphac / (14. * m_h) ;Electron number density
      ifrac1 = nel1/(nh1 + nhp1)                  ;Ionization fraction
      advionizy = dmdpt/m_h * vymdpt * (ifrac1-ifrac)/dx

      ; Z DIRECTION TERMS

      epsz = pmdpt2/dmdpt2/(5./3. - 1.) ;Specific internal energy
      depsz = epsz -epsx
      advheatlhsz = dmdpt * vzmdpt * depsz/dx
      advheatrhsz = pmdpt/dmdpt * vzmdpt * (dmdpt2-dmdpt)/dx

      dneutral2 = dratmdpt2*dmdpt2
      nh2 = dneutral2/m_h                        ;Number density of neutral H (atomic)
      nhp2 = (dmdpt2 - dneutral2)/m_h            ;Number density of ionized H
      nel2 = nhp2 + dmdpt2 * alphac / (14. * m_h) ;Electron number density
      ifrac2 = nel2/(nh2 + nhp2)                  ;Ionization fraction
      advionizz = dmdpt/m_h * vzmdpt * (ifrac2-ifrac)/dx

      ;Advection totals
      advioniz = advionizx + advionizy + advionizz
      advheatlhs = advheatlhsx+ advheatlhsy + advheatlhsz
      advheatrhs = advheatrhsx+ advheatrhsy + advheatrhsz

      if(keyword_set(ps)) then begin
         device, filen='~/Downloads/advheat_'+string(fileno, format='(I04)')+'.eps', ysize=10.85, xsize=8, /inches, /encapsulated
         !p.charsize = 1.5
         !p.charthick=4
         !x.thick=4
         !y.thick=4
      endif

      !p.multi=[0,2,4]
      ;Plot internal energy and PdV components
      plot, xmdpt/rp, advheatrhs, /nodata, title='Total', ytitle='Heating-Cooling Rates [erg cm!e-3!n s!e-1!n]'
      oplot, xmdpt/rp, advheatlhs, color=fsc_color("red")
      oplot, xmdpt/rp, advheatrhs, color=fsc_color("blue"), linestyle=2

      plot, xmdpt/rp, advioniz, title='Total',  ytitle='Ionization-Recombination Rates [cm!e-3!n s!e-1!n]'

      plot, xmdpt/rp, advheatlhsx, xsty=1, ysty=1, title='X component';, yra=[-1.5d-6, 1.5d-6]
      oplot, xmdpt/rp, advheatrhsx, linestyle=2

      plot, xmdpt/rp, advionizx, title='X component'

      plot, xmdpt/rp, advheatlhsy, title='Y component'
      oplot, xmdpt/rp, advheatrhsy, linestyle=2

      plot, xmdpt/rp, advionizy, title='Y component'

      plot, xmdpt/rp, advheatlhsz, title='Z component'
      oplot, xmdpt/rp, advheatrhsz, linestyle=2

      plot, xmdpt/rp, advionizz, title='Z component'


      if(keyword_set(ps)) then begin
         device, /close
      endif

      ;Plot heating/cooling rates
;      !p = p5  & !x= x5  & !y = y5

      if(keyword_set(ps)) then begin
         device, filen='~/Downloads/equilib_'+string(fileno, format='(I04)')+'.eps', ysize=10.85, xsize=8, /inches, /encapsulated
;         device, filen='~/Downloads/equilib_abs'+string(fileno, format='(I04)')+'.eps', ysize=10.85, xsize=8, /inches, /encapsulated

         !p.charsize = 1.5
         !p.charthick=4
         !x.thick=4
         !y.thick=4

      endif
      
      !p.multi=[0,1,2]

      alimamt = max([(bounds[1]-dx), rp]) ;;where(xmdpt/rp gt -20 AND xmdpt/rp lt 20)
      alevlimn = where(xmdpt lt - alimamt)
      alevlimp = where(xmdpt gt alimamt)

      rp = 1.5e10

      plot,  xmdpt/rp, rate_ionizheat, xsty =1, ysty =1, xtitle='x/R!dp!n', ytitle='Heating-Cooling Rates [erg cm!e-3!n s!e-1!n]', /nodata, xmargin=[15,5], yra=[-3d-7, 4d-7], xra=[-3, 3];, /ylog
minval = !y.crange[0]+1e-9
maxval = !y.crange[1]-1e-9

;      plot,  xmdpt/rp, rate_ionizheat, xsty =1, ysty =1, xtitle='x/R!dp!n', ytitle='Heating-Cooling Rates [erg cm!e-3!n s!e-1!n]', /nodata, xmargin=[15,5], yra=[1d-12, 4d-7], xra=[-3, 3], /ylog
;   minval = 10.^(!y.crange[0])+1e-20
;   maxval = 10.^(!y.crange[1])-1e-13

   polyfill, [-1,1,1,-1], [minval,minval, maxval, maxval], color=fsc_color("light gray"), /data

      oplot, xmdpt[alevlimn]/rp, rate_ionizheat[alevlimn], color=fsc_color("red")
      oplot, xmdpt[alevlimp]/rp, rate_ionizheat[alevlimp], color=fsc_color("red")
      oplot, xmdpt[alevlimn]/rp, (rate_recombcool[alevlimn]), color=fsc_color("blue")
      oplot, xmdpt[alevlimp]/rp, (rate_recombcool[alevlimp]), color=fsc_color("blue")
      oplot, xmdpt[alevlimn]/rp, (rate_lyacool[alevlimn]), color=fsc_color("teal")
      oplot, xmdpt[alevlimp]/rp, (rate_lyacool[alevlimp]), color=fsc_color("teal")
      ;; oplot, xmdpt/rp, abs(advheatx), color=fsc_color("brown")
      ;; oplot, xmdpt/rp, abs(advheaty), color=fsc_color("orange")
      ;; oplot, xmdpt/rp, abs(advheatz), color=fsc_color("yellow")
      ;; oplot, xmdpt/rp, abs(advheatx+advheaty+advheatz);, color=fsc_color("brown")
      ;; oplot, xmdpt/rp, abs(advtermx), color=fsc_color("purple")
      ;; oplot, xmdpt/rp, abs(advtermy), color=fsc_color("green")
      ;; oplot, xmdpt/rp, abs(advtermz), color=fsc_color("magenta")
      ;; oplot, xmdpt/rp, abs(advtermx+advtermy+advtermz), color=fsc_color("cyan")
      oplot, xmdpt[alevlimn]/rp, (advheatlhs[alevlimn]), color=fsc_color("orange")
      oplot, xmdpt[alevlimp]/rp, (advheatlhs[alevlimp]), color=fsc_color("orange")
      oplot, xmdpt[alevlimn]/rp, (advheatrhs[alevlimn]), color=fsc_color("green")
      oplot, xmdpt[alevlimp]/rp, (advheatrhs[alevlimp]), color=fsc_color("green")

      outsidea = alevlimn;where ((xmdpt ) ge rp)
      outsideb = alevlimp;where ((xmdpt ) le rp)

      oplot, xmdpt(outsidea)/rp, (rate_ionizheat(outsidea) + rate_recombcool(outsidea) + rate_lyacool(outsidea) + advheatrhs(outsidea))
      oplot, xmdpt(outsideb)/rp, (rate_ionizheat(outsideb) + rate_recombcool(outsideb) + rate_lyacool(outsideb) + advheatrhs(outsideb))

;      al_legend,['Photoioniz','Recomb','Lya', 'Int energy', 'PdV'],linestyle=[0,0,0,0,0],color=[fsc_color("red"), fsc_color("blue"), fsc_color("teal"), fsc_color("orange"), fsc_color("green")]

      p6 = !P & x6 = !X & y6 = !Y

      plot, xmdpt/rp, rate_ioniz, xsty =1, ysty =1, xtitle='x/R!dp!n', ytitle='Ionization-Recombination Rates [cm!e-3!n s!e-1!n]', /nodata, xmargin=[15,5], yra=[1d-1, 1d7], ymargin=[8,2], /ylog, xra=[-3,3]
      minval = 10.^(!y.crange[0])+1e-20
      maxval = 10.^(!y.crange[1])-1e-13
      polyfill, [-1,1,1,-1], [minval,minval, maxval, maxval], color=fsc_color("light gray"), /data

      oplot, xmdpt[alevlimn]/rp, rate_ioniz[alevlimn], color=fsc_color("red")
      oplot, xmdpt[alevlimp]/rp, rate_ioniz[alevlimp], color=fsc_color("red")
      oplot, xmdpt[alevlimn]/rp, abs(rate_recomb[alevlimn]), color=fsc_color("blue")
      oplot, xmdpt[alevlimp]/rp, abs(rate_recomb[alevlimp]), color=fsc_color("blue")
      oplot, xmdpt[alevlimn]/rp, advioniz[alevlimn], color=fsc_color("orange")
      oplot, xmdpt[alevlimp]/rp, advioniz[alevlimp], color=fsc_color("orange")

      ;; oplot, xmdpt/rp, advionizx, color=fsc_color("brown")
      ;; oplot, xmdpt/rp, advionizy, color=fsc_color("orange");, linestyle=2
      ;; oplot, xmdpt/rp, advionizz, color=fsc_color("yellow");, linestyle=5

      al_legend,['Photoioniz','Recomb','Advection'],color=[fsc_color("red"), fsc_color("blue"), fsc_color("teal"), fsc_color("orange")]

      p7 = !P & x7 = !X & y7 = !Y
;;       wset, 10
;; !p.multi=0
;;     plot, xmdpt/rp, rate_ioniz, xsty =1, ysty =1, xtitle='x/R!dp!n', ytitle='Ionization-Recombination Rates [cm!e-3!n s!e-1!n]', /nodata, xmargin=[15,5], yra=[-100,100], ymargin=[8,2]
;;       oplot, xmdpt/rp, rate_ioniz, color=fsc_color("red")
;;       oplot, xmdpt/rp, -rate_recomb, color=fsc_color("blue")
;;       oplot, xmdpt/rp, advionizx, color=fsc_color("brown")
;;       oplot, xmdpt/rp, advionizy, color=fsc_color("orange");, linestyle=2
;;       oplot, xmdpt/rp, advionizz, color=fsc_color("yellow");, linestyle=5
;;       oplot, xmdpt/rp, (advionizx+advionizy+advionizz);, color=fsc_color("yellow")

;; wset, 11
;;       plot,  xmdpt/rp, rate_ionizheat, xsty =1, ysty =1, xtitle='x/R!dp!n', ytitle='Heating-Cooling Rates [erg cm!e-3!n s!e-1!n]', /nodata, xmargin=[15,5], yra=[-1d-8, 1d-8]
;;    minval = 10.^(!y.crange[0])+1e-20
;;    maxval = 10.^(!y.crange[1])-1e-13
;;    polyfill, [-1,1,1,-1], [minval,minval, maxval, maxval], color=fsc_color("light gray"), /data

;;       oplot, xmdpt/rp, rate_ionizheat, color=fsc_color("red")
;;       oplot, xmdpt/rp, -rate_recombcool, color=fsc_color("blue")
;;       oplot, xmdpt/rp, -rate_lyacool, color=fsc_color("teal")
;;       oplot, xmdpt/rp, advheatx, color=fsc_color("brown")
;;       oplot, xmdpt/rp, advheaty, color=fsc_color("orange")
;;       oplot, xmdpt/rp, advheatz, color=fsc_color("yellow")
;;       oplot, xmdpt/rp, advheatx+advheaty+advheatz;, color=fsc_color("brown")
;;     oplot, xmdpt/rp, abs(advtermx), color=fsc_color("purple")
;;       oplot, xmdpt/rp, abs(advtermy), color=fsc_color("green")
;;       oplot, xmdpt/rp, abs(advtermz), color=fsc_color("magenta")
;;       oplot, xmdpt/rp, abs(advtermx+advtermy+advtermz) , color=fsc_color("cyan")

   for ilev = 1, 3 do begin
      if (keyword_set(yplot)) then begin
         filen = 'ioniz_sphere-lev'+strtrim(ilev,2)+'.'+string(fileno, format='(I04)')+'y_mdptvals.sav'
      endif else begin
         filen = 'ioniz_sphere-lev'+strtrim(ilev,2)+'.'+string(fileno, format='(I04)')+'_mdptvals.sav'
      endelse

      plotdvt, filen, ilev, 0, 1
   endfor


if(keyword_set(ps)) then begin
   device, /close
endif

stop


END
