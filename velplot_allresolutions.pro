FUNCTION alphab, T
a = 2.59e-13*(T/1.0d4)^(-0.7)
RETURN,a
END

PRO plotdvt, infile, inlev, plotifrac, plotheating
  common PHYSVARS, kb, m_h, mu, alphac, rp, km, sigmaph, egamma
  common PLOTPOS, p1, p2, p3, x1, x2, x3, y1, y2, y3, p4, x4, y4, p5, x5, y5, p6, x6, y6, p7, x7, y7, p8, x8, y8, p9, x9, y9
  common LEVLIMS, bounds, filenum, taucum
  common OFFAXIS, offa
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

  ;Adiabat
  csmdpt = csmdpt * sqrt(5./3.)

  ;Find sonic point on day side
  neg = where( xmdpt lt 0)
  sonic = where (abs(abs(vxmdpt[neg])-csmdpt[neg]) eq min(abs(abs(vxmdpt[neg]) -csmdpt[neg])))
 print,'Lev ', strtrim(inlev,2),' thinks vx/acs =', abs(vxmdpt[sonic])/csmdpt[sonic], ' at ', xmdpt[sonic]/rp  

  ;Calculate optical depth
  tau = dblarr(n)
  dx = abs(xmdpt[1]-xmdpt[0])
  for i=0,n-1 do begin   
     tau[i] = taucum[inlev-1]+total(nh[0:i] *dx *sigmaph) ; - total(d[*, 0, floor(n/2)-1]* dx)
  endfor
  taucum[inlev] = taucum[inlev-1]+total(nh[levlimn] *dx *sigmaph)
print, inlev, taucum[inlev]
;*2  taucum[inlev*2+1] = total(tau[levlimp]*dx *sigmaph)
   
   ;Find tau=1 location
   tau1 = where(abs(tau - 1) eq min(abs(tau -1)))
; print,'Lev ', strtrim(inlev,2),' thinks tau1 =', tau[tau1], ' at ', xmdpt[tau1]/rp


if (plotheating lt 1) then begin
   if (plotifrac lt 1) then begin
                                ;OVERPLOT
                                ;Plot 1D density
      !p=p1 & !x = x1 & !y = y1
      oplot, xmdpt[levlimn]/rp, dmdpt[levlimn] ;, color=fsc_color("green")
      oplot, xmdpt[levlimp]/rp, dmdpt[levlimp] ;, color=fsc_color("green")

      ;Based on the finding that Lev3 came the closest to tau=1
      if (inlev eq 3) then begin
         plots, xmdpt[tau1]/rp,dmdpt(tau1), psym =7
      endif

                                ;Plot velocity
      !p = p2  & !x= x2  & !y = y2
;      oplot, xmdpt[levlimn]/rp, vmagmdpt[levlimn]/km ;, color=fsc_color("green")
;      oplot, xmdpt[levlimp]/rp, vmagmdpt[levlimp]/km ;, color=fsc_color("green")

      oplot, xmdpt[levlimn]/rp, vxmdpt[levlimn]/km; , color=fsc_color("green")
      oplot, xmdpt[levlimp]/rp, vxmdpt[levlimp]/km; , color=fsc_color("green")
;oplot, xmdpt[levlimn]/rp, csmdpt[levlimn]/km, color=fsc_color("blue") ;, color=fsc_color("green")
;oplot, xmdpt[levlimp]/rp, csmdpt[levlimp]/km, color=fsc_color("blue") ;, color=fsc_color("green")
;       plots, xmdpt[neg[sonic]]/rp, vxmdpt[neg[sonic]]/km, psym = 8, color=fsc_color("blue")

                                ;Plot log T
      !p = p3  & !x= x3  & !y = y3
      oplot, xmdpt[levlimn]/rp, temp1[levlimn] ;,color=fsc_color("green")
      oplot, xmdpt[levlimp]/rp, temp1[levlimp] ;,color=fsc_color("green")
      ;; plots, xmdpt[neg[sonic]]/rp, temp1[neg[sonic]], psym = 8

      !p = p4  & !x= x4  & !y = y4
;     levlimn = where(xmdpt lt -1.*(bounds[inlev+1]-dx))
;     levlimp = where(xmdpt gt (bounds[inlev+1]-dx))
      if (inlev eq 3) then begin
         lt1 = where (xmdpt lt 1.*rp)
         oplot, xmdpt[lt1]/rp, ifrac[lt1];, color=fsc_color("blue")
         plots, xmdpt[tau1]/rp,ifrac[tau1], psym =7
      endif else begin
         oplot, xmdpt[levlimn]/rp, ifrac[levlimn]; , color=fsc_color("green")
      endelse
      



   endif else begin
      !p = p4  & !x= x4  & !y = y4
;     levlimn = where(xmdpt lt -1.*(bounds[inlev+1]-dx))
;     levlimp = where(xmdpt gt (bounds[inlev+1]-dx))

      oplot, xmdpt[levlimn]/rp, ifrac[levlimn];, color=fsc_color("green")
      oplot, xmdpt[levlimp]/rp, ifrac[levlimp];, color=fsc_color("green")



      !p = p5  & !x= x5  & !y = y5      
      oplot, xmdpt/rp, tau, color=fsc_color("red")
      oplot, xmdpt[levlimn]/rp, tau[levlimn], color=fsc_color("blue"), thick=5


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
      if (offa > 0) then begin
         filen = 'ioniz_sphere-lev'+strtrim(inlev,2)+'.'+string(filenum, format='(I04)')+'_a_offmdptvals.sav'
      endif else begin
         filen = 'ioniz_sphere-lev'+strtrim(inlev,2)+'.'+string(filenum, format='(I04)')+'_offmdptvals.sav'
      endelse
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
         alimamt = max([(bounds[inlev+1]-dx)]) ;;where(xmdpt/rp gt -20 AND xmdpt/rp lt 20)
      endif else begin
         alimamt = 0;rp ;;where(xmdpt/rp gt -20 AND xmdpt/rp lt 20)
      endelse
      alevlimn = where(xmdpt lt -1.*alimamt)
      alevlimp = where(xmdpt gt alimamt)

      !p = p6  & !x= x6  & !y = y6




      oplot, xmdpt[alevlimn]/rp, (rate_ionizheat[alevlimn] + rate_recombcool[alevlimn] + rate_lyacool[alevlimn] + advheatrhs[alevlimn]), thick =14
      oplot, xmdpt[alevlimn]/rp, rate_ionizheat[alevlimn], color=fsc_color("red6", /brewer), thick=10
      oplot, xmdpt[alevlimn]/rp, (rate_lyacool[alevlimn]), color=fsc_color("grn8", /brewer), thick=10
      oplot, xmdpt[alevlimn]/rp, (rate_recombcool[alevlimn]), color=fsc_color("ygb4", /brewer), thick=10
      oplot, xmdpt[alevlimn]/rp, (advheatlhs[alevlimn]), color=fsc_color("org3", /brewer), linestyle=5, thick = 10
      oplot, xmdpt[alevlimn]/rp, (advheatrhs[alevlimn]), color=fsc_color("blu7", /brewer), linestyle=5, thick = 10
      ;; oplot, xmdpt[alevlimn]/rp, (advheatlhsx[alevlimn]), color=fsc_color("gray")
      ;; oplot, xmdpt[alevlimn]/rp, (advheatlhsy[alevlimn]), color=fsc_color("gray"), linestyle=1
      ;; oplot, xmdpt[alevlimn]/rp, (advheatlhsz[alevlimn]), color=fsc_color("gray"), linestyle=2



      ;; oplot, xmdpt[alevlimn]/rp, rate_ionizheat[alevlimn], color=fsc_color("red")
      ;; oplot, xmdpt[alevlimn]/rp, (rate_recombcool[alevlimn]), color=fsc_color("blue")
      ;; oplot, xmdpt[alevlimn]/rp, (rate_lyacool[alevlimn]), color=fsc_color("teal")
      ;; oplot, xmdpt[alevlimn]/rp, (advheatlhs[alevlimn]), color=fsc_color("orange")
      ;; oplot, xmdpt[alevlimn]/rp, (advheatrhs[alevlimn]), color=fsc_color("green")
      ;; oplot, xmdpt[alevlimn]/rp, (rate_ionizheat[alevlimn] + rate_recombcool[alevlimn] + rate_lyacool[alevlimn] + advheatrhs[alevlimn])
      ;; oplot, xmdpt[alevlimn]/rp, (advheatlhsx[alevlimn]), color=fsc_color("gray")
      ;; oplot, xmdpt[alevlimn]/rp, (advheatlhsy[alevlimn]), color=fsc_color("gray"), linestyle=1
      ;; oplot, xmdpt[alevlimn]/rp, (advheatlhsz[alevlimn]), color=fsc_color("gray"), linestyle=2
      ttemp = (rate_ionizheat[alevlimn] + rate_recombcool[alevlimn] + rate_lyacool[alevlimn] + advheatrhs[alevlimn])

      
      print, min(ttemp), max(ttemp)

      !p = p7  & !x= x7  & !y = y7
      oplot, xmdpt[alevlimp]/rp, (rate_ionizheat[alevlimp] + rate_recombcool[alevlimp] + rate_lyacool[alevlimp] + advheatrhs[alevlimp]), thick=14
      oplot, xmdpt[alevlimp]/rp, rate_ionizheat[alevlimp], color=fsc_color("red6"), thick=10
      oplot, xmdpt[alevlimp]/rp, (rate_recombcool[alevlimp]), color=fsc_color("ygb4"), thick=10
      oplot, xmdpt[alevlimp]/rp, (rate_lyacool[alevlimp]), color=fsc_color("grn8"), thick=10
      oplot, xmdpt[alevlimp]/rp, (advheatlhs[alevlimp]), color=fsc_color("org3"), linestyle=5, thick = 10
      oplot, xmdpt[alevlimp]/rp, (advheatrhs[alevlimp]), color=fsc_color("blu7"), linestyle=5, thick = 10
      ;; oplot, xmdpt[alevlimp]/rp, (advheatlhsx[alevlimp]), color=fsc_color("gray")
      ;; oplot, xmdpt[alevlimp]/rp, (advheatlhsy[alevlimp]), color=fsc_color("gray"), linestyle=1
      ;; oplot, xmdpt[alevlimp]/rp, (advheatlhsz[alevlimp]), color=fsc_color("gray"), linestyle=2

      ttemp = (rate_ionizheat[alevlimp] + rate_recombcool[alevlimp] + rate_lyacool[alevlimp] + advheatrhs[alevlimp])
      print, min(ttemp), max(ttemp)
      ;; oplot, xmdpt/rp, abs(advheatx), color=fsc_color("brown")
      ;; oplot, xmdpt/rp, abs(advheaty), color=fsc_color("orange")
      ;; oplot, xmdpt/rp, abs(advheatz), color=fsc_color("yellow")
      ;; oplot, xmdpt/rp, abs(advheatx+advheaty+advheatz);, color=fsc_color("brown")
      ;; oplot, xmdpt/rp, abs(advtermx), color=fsc_color("purple")
      ;; oplot, xmdpt/rp, abs(advtermy), color=fsc_color("green")
      ;; oplot, xmdpt/rp, abs(advtermz), color=fsc_color("magenta")
      ;; oplot, xmdpt/rp, abs(advtermx+advtermy+advtermz), color=fsc_color("cyan")

      print, min(advheatrhs), max(rate_ionizheat)

;;       print, 'Annoying stop - temporary. Accept after all levs have plotted heating'
;; stop

      !p = p8  & !x= x8  & !y = y8
      oplot, xmdpt[alevlimn]/rp, rate_ioniz[alevlimn], color=fsc_color("red6"), thick=10
      oplot, xmdpt[alevlimn]/rp, (rate_recomb[alevlimn]), color=fsc_color("ygb4"), thick=10
      oplot, xmdpt[alevlimn]/rp, (advioniz[alevlimn]), color=fsc_color("orange"), thick=10, linestyle=5
;      oplot, xmdpt[alevlimn]/rp, (advionizx[alevlimn]), color=fsc_color("gray"), linestyle=3
;      oplot, xmdpt[alevlimn]/rp, (advionizy[alevlimn]), color=fsc_color("gray"), linestyle=1
;      oplot, xmdpt[alevlimn]/rp, (advionizz[alevlimn]), color=fsc_color("gray"), linestyle=2
;      oplot, xmdpt[alevlimn]/rp,
;      (rate_ioniz[alevlimn]+rate_recomb[alevlimn]), linestyle=3

;; fmt = '(E10.3, E10.3, E10.3, E10.3, E10.3)'
;; forprint, F=fmt, xmdpt[0:n_elements(xmdpt)-2], rate_ionizheat[0:n_elements(xmdpt)-2], rate_recombcool[0:n_elements(xmdpt)-2], advheatlhs, advheatrhs, text='~/Downloads/someratesforsimeon_lev'+strtrim(inlev,2)+'.dat'


      !p = p9  & !x= x9  & !y = y9
      oplot, xmdpt[alevlimp]/rp, rate_ioniz[alevlimp], color=fsc_color("red6"), thick=10
      oplot, xmdpt[alevlimp]/rp, (advioniz[alevlimp]), color=fsc_color("orange"), thick=10, linestyle=5
;      oplot, xmdpt[alevlimp]/rp, (advionizx[alevlimp]), color=fsc_color("gray"), linestyle=3
;      oplot, xmdpt[alevlimp]/rp, (advionizy[alevlimp]), color=fsc_color("gray"), linestyle=1
;      oplot, xmdpt[alevlimp]/rp, (advionizz[alevlimp]), color=fsc_color("gray"), linestyle=2
;      oplot, xmdpt[alevlimp]/rp, (rate_ioniz[alevlimp]+rate_recomb[alevlimp]), linestyle=3
      oplot, xmdpt[alevlimp]/rp, (rate_recomb[alevlimp]), color=fsc_color("ygb4"), thick=10;, linestyle =2      

;;       oplot, xmdpt[alevlimn]/rp, rate_ioniz[alevlimn], color=fsc_color("red")
;;       oplot, xmdpt[alevlimp]/rp, rate_ioniz[alevlimp], color=fsc_color("red")
;;       oplot, xmdpt[alevlimn]/rp, abs(rate_recomb[alevlimn]), color=fsc_color("blue")
;;       oplot, xmdpt[alevlimp]/rp, abs(rate_recomb[alevlimp]), color=fsc_color("blue")
;;       oplot, xmdpt[alevlimn]/rp, abs(advioniz[alevlimn]), color=fsc_color("orange"), thick=3
;;       oplot, xmdpt[alevlimp]/rp, abs(advioniz[alevlimp]), color=fsc_color("orange")
;;       oplot, xmdpt[alevlimn]/rp, abs(advionizx[alevlimn]), color=fsc_color("cyan")
;;       oplot, xmdpt[alevlimn]/rp, abs(advionizy[alevlimn]), color=fsc_color("magenta");, linestyle=2
;;       oplot, xmdpt[alevlimn]/rp, abs(advionizz[alevlimn]), color=fsc_color("yellow");, linestyle=5
;;       oplot, xmdpt[alevlimp]/rp, abs(advionizx[alevlimp]), color=fsc_color("cyan")
;;       oplot, xmdpt[alevlimp]/rp, abs(advionizy[alevlimp]), color=fsc_color("magenta");, linestyle=2
;;       oplot, xmdpt[alevlimp]/rp, abs(advionizz[alevlimp]), color=fsc_color("yellow");, linestyle=5
;; print, 'Probe here'
; stop


;      oplot, xmdpt
;stop




endelse

END



PRO velplot_allresolutions, fileno, ps=ps, recombtime=recombtime, yplot=yplot, offmdpt=offmdpt;, fromsaved=fromsaved
common PHYSVARS, kb, m_h, mu, alphac, rp, km, sigmaph, egamma
common PLOTPOS, p1, p2, p3, x1, x2, x3, y1, y2, y3, p4, x4, y4, p5, x5, y5, p6, x6, y6, p7, x7, y7, p8, x8, y8, p9, x9, y9
common LEVLIMS, bounds, filenum, taucum
common OFFAXIS, offa

filenum = fileno

if (keyword_set(offmdpt)) then begin
   offa = 1
endif else begin
   offa = 0
endelse

;;Currently hardcoded. Can be changed/checked later
bounds=[7.4062500e+10, 2.390625d10, 2.2265625e+10, 2.0507812e+10, 1.8691406e+10]

;Plotting setup
plotsym, 0, 1., /fill
set_plot, 'x'
if (keyword_set(ps)) then begin
    !p.font=0
    !p.charsize=2
    !p.charthick=5
    !p.thick=8
    !x.thick=5
    !y.thick=5
endif

;Physical setup
kb = 1.38e-16
m_h = 1.67d-24
mu = 1.67d-24                   ;Neutral gas mean molec weight
alphac = 1d-10
rp = 1.5d10                     ;CHANGE ACCORDINGLY
km = 1d6; This is actualy 10km in cm
sigmaph = 6.3d-18               ;CHANGE ACCORDINGLY
egamma = 3.85e-12               ;2.4eV in ergs
ggrav = 6.67e-8
mp = 1d30

;Convert Ruth's .dat to a .sav
;readcol, 'tideflux4.5E+03_papermu.dat', format='D,D,D,D,D,D,D,D,D', rr, rd, rv, rT, rnf, rtau , rrq, rrz, rrtau
;save, /variables, filename='tideflux_papermu.sav'  

;Read in Ruth's wind solution
;restore, '/Users/anjalitripathi/Atmospheric-Athena/bin/tideflux_papermu.sav'
restore, '/Users/anjalitripathi/Atmospheric-Athena/bin/wind_011315.sav' ; Variables rr, rd, rv, rT, rnf, rtau
                                ;Columns:  R(10^10 cm)  rho(10^-15
                                ;g/cm^3)  v(10^6 cm/s)  T(10^4 K)  fp
                                ;tau
;For new files, R is in units of 1.5e10
rx = -rr/(rp/1d10)
rtemp = rT * 1d4
rdens = rd*1e-15
rmu = (rnf * m_h/2) + (1.-rnf)*mu
rcs = sqrt(kb/rmu * rtemp)*sqrt(5./3.)       ;Is neutral mu appropriate here?

; RESTORE Lowest resolution FILE
  if (keyword_set(yplot)) then begin
     restore, 'ioniz_sphere.'+string(fileno, format='(I04)')+'y_mdptvals.sav'
  endif else begin
     ;!!! Off mdpt
     if (offa > 0) then begin
        restore, 'ioniz_sphere.'+string(fileno, format='(I04)')+'_a_mdptvals.sav'
     endif else begin
        restore, 'ioniz_sphere.'+string(fileno, format='(I04)')+'_mdptvals.sav'
     endelse
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

  ;Adiabat
  csmdpt = csmdpt * sqrt(5./3.)

  ;Find sonic point on day side
  neg = where( xmdpt lt 0, complement = nneg)
  sonic = where (abs(abs(vxmdpt[neg])-csmdpt[neg]) eq min(abs(abs(vxmdpt[neg]) -csmdpt[neg])))
  sonic2 = where (abs(vxmdpt[nneg]-csmdpt[nneg]) eq min(abs(vxmdpt[nneg] -csmdpt[nneg])))
  rsonic = where (abs(rv*km - rcs) eq min(abs(rv*km - rcs))) ;For Ruth's output
 print,'Lev 0 thinks vx/acs =', abs(vxmdpt[sonic])/csmdpt[sonic], ' at ', xmdpt[sonic]/rp  
  
  ;Calculate optical depth
  tau = dblarr(n)
  dx = abs(xmdpt[1]-xmdpt[0])
  for i=0,n-1 do begin   
     tau[i] = total(nh[0:i] *dx *sigmaph) ; - total(d[*, 0, floor(n/2)-1]* dx)
  endfor
  taucum = dblarr(5);*2)
  taucum[0] = total(nh[levlimn]*dx *sigmaph)   
  ;Find tau=1 location
  tau1 = where(abs(tau - 1) eq min(abs(tau -1)))
  rtau1 = where(abs(rtau - 1) eq min(abs(rtau -1))) ;For Ruth's output

  print,'Lev 0 thinks tau =',  tau[tau1], ' at ', xmdpt[tau1]/rp
  
  if (keyword_set(recombtime)) then begin
;DEPRECATED
;      window, 2
     ;; !p.multi=0
     ;; !p.charsize = 1.5
     ;; plot, xmdpt/rp, 1./(nhp+1d-12)/3d-13, xsty = 1, ysty=2, /ylog, yra=[1, 1e9], ytitle='Time', xtitle='X/1.5e10 [cm]'
     ;; oplot, xmdpt/rp,(2.7*1.5d10)/sqrt(pmdpt/dmdpt), color=fsc_color("red")
     ;; oplot, xmdpt/rp,(1.5d10)/sqrt(pmdpt/dmdpt),
     ;; color=fsc_color("blue")

     print, 'Plotting recombination time. Check length scale'
     tiny = 0.
     trecomb = 1./alphab(temp1)/(nhp + tiny)
     glocal = ggrav * mp/xmdpt^2.
     h = csmdpt*csmdpt/glocal
     window, 21
     !p.multi=0
     plot, xmdpt/rp, 1e11/abs(vxmdpt), yra=[1e3, 1e6], /ylog
     oplot, xmdpt/rp, abs(xmdpt)/abs(csmdpt), color=fsc_color("red")
     oplot, xmdpt/rp, trecomb, linestyle=1
     stop
  endif 


   ;PLOTS

   if (keyword_set(ps)) then begin
      set_plot, 'ps'
      !p.charsize=1.7
      !y.charsize = 1.7
      !x.charsize = 1.7
      !p.thick=9
      !x.thick=6
      !y.thick=6


      if (keyword_set(yplot)) then begin
         device, filen='~/Downloads/velplot_'+string(fileno, format='(I04)')+'y_4lev.eps', ysize=10.5, xsize=8.2, /inches, /encapsulated
      endif else begin
         if (offa > 0) then begin
            device, filen='~/Downloads/velplot_a_'+string(fileno, format='(I04)')+'_4lev.eps', ysize=10.5, xsize=8.2, /inches, /encapsulated
         endif else begin
            device, filen='~/Downloads/velplot_4panel'+string(fileno, format='(I04)')+'_4lev.eps', ysize=10.5, xsize=8.2, /inches, /encapsulated
         endelse
      endelse
  endif else begin
     window
     !p.charsize = 2
      print, 'Adjust window size now'
      stop
  endelse
   !p.multi=[0,1,4]

   ;Plot 1D density
   plot, xmdpt/rp, dmdpt, xsty = 1, ysty = 2, xtitle ='x/R!dp!n', ytitle='Density [g cm!e-3!n]', /ylog, yra=[7e-19, 2e-13], xmargin=[15,5], /nodata;, ycharsize = 1.1, xmargin=[5,5];, title='t/1e5s=0'+strtrim(ii-1,2)
   minval = 10.^(!y.crange[0])+1e-20
   maxval = 10.^(!y.crange[1])-1e-13
   polyfill, [-1,1,1,-1], [minval,minval, maxval, maxval], color=fsc_color("light gray"), /data
   oplot, xmdpt[levlimn]/rp,  dmdpt[levlimn]
   oplot, xmdpt[levlimp]/rp,  dmdpt[levlimp]
   oplot, rx, rdens, linestyle = 1
;   plots, xmdpt[tau1]/rp,dmdpt(tau1), psym =7 ;Let Level "3" plot it
   plots, rx[rtau1], rdens(rtau1), psym =7
   plots, xmdpt[neg[sonic]]/rp, dmdpt[neg[sonic]], psym = 8
   plots, xmdpt[nneg[sonic2]]/rp, dmdpt[nneg[sonic2]], psym = 8
   plotsym, 0, .8, thick=6
   plots, rx[rsonic], rdens[rsonic], psym =8 ;, color=fsc_color("gray")

   ;Save plot coordinates for returning later
   p1 = !P & x1 = !X & y1 = !Y


   ;Plot 1D velocity
   plot, xmdpt/rp, vxmdpt/km, xsty = 1, ysty= 1, ytitle='Velocity [10 km s!e-1!n]', xtitle ='x/R!dp!n',  yra=[-2.6, 2.6], /nodata, xmargin=[15,5];ycharsize = 1.1 
   minval = (!y.crange[0])+2e-2
   maxval = (!y.crange[1])-2e-2
   polyfill, [-1,1,1,-1], [minval,minval, maxval, maxval], color=fsc_color("light gray"), /data
;   oplot, xmdpt[levlimn]/rp, vmagmdpt[levlimn]/km
;   oplot, xmdpt[levlimp]/rp, vmagmdpt[levlimp]/km
   oplot, xmdpt[levlimn]/rp, vxmdpt[levlimn]/km ;, color=fsc_color("green")
   oplot, xmdpt[levlimp]/rp, vxmdpt[levlimp]/km ;, color=fsc_color("green")

   ;oplot, xmdpt/rp, vxmdpt/km, linestyle = 4, thick = 3
   oplot, rx, -rv, linestyle=1
   plots, rx[rsonic], -rv[rsonic], psym =8 ;, color=fsc_color("gray")
   plotsym, 0, 1., /fill
   plots, xmdpt[neg[sonic]]/rp, vxmdpt[neg[sonic]]/km, psym = 8
   plots, xmdpt[nneg[sonic2]]/rp, vxmdpt[nneg[sonic2]]/km, psym = 8

;   oplot, xmdpt/rp, csmdpt/km, color=fsc_color("blue")
   ;oplot, [-.75, -.75], [0, 10], linestyle = 2, color=fsc_color("gray");, thick = 5
   ;oplot, [.75, .75], [0, 10], linestyle = 2, color=fsc_color("gray");, thick = 5
   ;plot, xmdpt/rp, vmagmdpt/csmdpt, xsty = 1, ysty = 1, xtitle ='x/R!dp!n', ytitle='Velocity/Sound speed'

   ;Save plot coordinates for returning later
   p2 = !P & x2 = !X & y2 = !Y


   ;; Plot timescale
;;     dvdx = (shift(vxmdpt,-1) - vxmdpt)/(shift(xmdpt,-1) - xmdpt) 
;;     plot, xmdpt/rp, (1./dvdx), /xsty, /ysty, title='dx/dv vs. x';, xra=[1,5], /ylog
;;     oplot, xmdpt/rp, 1.5e11/csmdpt, color=fsc_color("red")


   ;Plot log T
   plot, xmdpt/rp, temp1, xsty = 1, ysty = 1, xtitle ='x/R!dp!n', /ylog, yra=[9e2, 3e4], /nodata, xmargin=[15,5], ytickformat='(A1)';, ycharsize = 1.1 ;, linestyle = 4, yra=[9e2, 4e4]
;, ymargin=[8,2]
   axis, yaxis=0, yrange=[9e2,3e4], /ylog, ysty=1, ytickv=[1e3, 1e4], ytickname=['10!u3!n','10!u4!n'], ytitle='Temperature [K]'
   minval = 10.^(!y.crange[0])+2e1
   maxval = 10.^(!y.crange[1])-1e3
   polyfill, [-1,1,1,-1], [minval,minval, maxval, maxval], color=fsc_color("light gray"), /data
   oplot, xmdpt[levlimn]/rp, temp1[levlimn]
   oplot, xmdpt[levlimp]/rp, temp1[levlimp]
   oplot, rx, rtemp, linestyle=1
   plots, xmdpt[neg[sonic]]/rp, temp1[neg[sonic]], psym = 8
   plots, xmdpt[nneg[sonic2]]/rp, temp1[nneg[sonic2]], psym = 8
   plotsym, 0, .8, thick=6
   plots, rx[rsonic], rtemp[rsonic], psym =8 ;, color=fsc_color("gray")

   p3 = !P & x3 = !X & y3 = !Y

;Plot ifrac
   plot, xmdpt/rp, ifrac, xsty =1, ysty =1, xtitle='x/R!dp!n', ytitle='Ionization fraction', yra=[0,1], /nodata, xmargin=[15,5], ymargin=[8,2] ;, title='t/1e5s=0'+strtrim(ii-1,2)
   minval = (!y.crange[0])+5e-3
   maxval = (!y.crange[1])-5e-3
   polyfill, [-1,1,1,-1], [minval,minval, maxval, maxval], color=fsc_color("light gray"), /data
   gt0 = where (xmdpt gt 0)
   oplot, xmdpt[gt0]/rp, ifrac[gt0]
   oplot, xmdpt[levlimn]/rp, ifrac[levlimn]
 ;  oplot, xmdpt[levlimp]/rp, ifrac[levlimp]
   oplot, rx, rnf, linestyle = 1
;   plots, xmdpt[tau1]/rp, ifrac[tau1], psym = 7
   plots, rx[rtau1], rnf[rtau1], psym = 7 ;, color=fsc_color("gray")
   p4 = !P & x4 = !X & y4 = !Y



   for ilev = 1, 3 do begin
      if (keyword_set(yplot)) then begin
         filen = 'ioniz_sphere-lev'+strtrim(ilev,2)+'.'+string(fileno, format='(I04)')+'y_mdptvals.sav'
      endif else begin
         ;!!!
         if (offa > 0) then begin
            filen = 'ioniz_sphere-lev'+strtrim(ilev,2)+'.'+string(fileno, format='(I04)')+'_a_mdptvals.sav'
         endif else begin
            filen = 'ioniz_sphere-lev'+strtrim(ilev,2)+'.'+string(fileno, format='(I04)')+'_mdptvals.sav'
         endelse
      endelse

      plotdvt, filen, ilev, 0, 0
   endfor

   if(keyword_set(ps)) then begin
      device, /close
   endif
    print, 'Will move onto ifrac after break'
    stop



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
         if (offa > 0) then begin
            device, filen='~/Downloads/joined_ifrac_a_'+string(fileno, format='(I04)')+'_anim.eps', ysize=10.5, xsize=8, /inches, /encapsulated
         endif else begin
            device, filen='~/Downloads/joined_ifrac_'+string(fileno, format='(I04)')+'_anim.eps', ysize=10.5, xsize=8, /inches, /encapsulated
         endelse
      endelse
      
;   !p.POSITION=[.1,.1,.9,.9] 
   endif
   !p.multi=[0,1,2]
   
   if (keyword_set(yplot)) then begin
      restore, 'ioniz_sphere.'+string(fileno, format='(I04)')+'y_mdptvals.sav'
   endif else begin
      if (offa>0) then begin
         restore, 'ioniz_sphere.'+string(fileno, format='(I04)')+'_a_mdptvals.sav'
      endif else begin
         restore, 'ioniz_sphere.'+string(fileno, format='(I04)')+'_mdptvals.sav'
      endelse
   endelse
   
print, 'STOP, ifrac plotting has been moved into the 1st plot. Move it back to continue'
;Plot ifrac
;;    plot, xmdpt/rp, ifrac, xsty =1, ysty =1, xtitle='x/R!dp!n', ytitle='Ionization fraction', yra=[0,1], /nodata, xmargin=[15,5] ;, title='t/1e5s=0'+strtrim(ii-1,2)
;;    oplot, xmdpt[levlimn]/rp, ifrac[levlimn]
;;    oplot, xmdpt[levlimp]/rp, ifrac[levlimp]
;;    oplot, rx, rnf, linestyle = 1
;;    plots, xmdpt[tau1]/rp, ifrac[tau1], psym = 7
;;    plots, rx[rtau1], rnf[rtau1], psym = 7 ;, color=fsc_color("gray")
;;    p4 = !P & x4 = !X & y4 = !Y

;Plot tau
   plot, xmdpt/rp, tau, xsty = 1, ysty = 1, xtitle ='x/R!dp!n', ytitle='Optical depth', /ylog, yra=[8e-5, 1.5e4], xmargin=[15,5], ymargin=[8,2] ;
   polyfill,  [min(xmdpt/rp)+5d-2,max(xmdpt/rp)-5d-2, max(xmdpt/rp)-5d-2, min(xmdpt/rp)+5d-2], [1., 1., 10.^(!y.crange[1])-40, 10.^(!y.crange[1]) -40 ], color=fsc_color("light gray"), /data
   oplot, xmdpt/rp, tau
   oplot, rx, rtau, linestyle = 1
   p5 = !P & x5 = !X & y5 = !Y
   
   for ilev = 1, 3 do begin
      if (keyword_set(yplot)) then begin
         filen = 'ioniz_sphere-lev'+strtrim(ilev,2)+'.'+string(fileno, format='(I04)')+'y_mdptvals.sav'
      endif else begin
         if (offa > 0) then begin
            filen = 'ioniz_sphere-lev'+strtrim(ilev,2)+'.'+string(fileno, format='(I04)')+'_a_mdptvals.sav'
         endif else begin
            filen = 'ioniz_sphere-lev'+strtrim(ilev,2)+'.'+string(fileno, format='(I04)')+'_mdptvals.sav'
         endelse
      endelse
      plotdvt, filen, ilev, 1, 0
   endfor

   if(keyword_set(ps)) then begin
      device, /close
   endif   
   
   print, 'Pausing before equilibrium rates'


;Determine equilib quantitites

   if (offa > 0) then begin
      restore, 'ioniz_sphere.'+string(fileno, format='(I04)')+'_a_mdptvals.sav'
   endif else begin
      restore, 'ioniz_sphere.'+string(fileno, format='(I04)')+'_mdptvals.sav'
   endelse

   dneutral = dratmdpt*dmdpt
   nh = dneutral/m_h                       ;Number density of neutral H (atomic)
   nhp = (dmdpt - dneutral)/m_h            ;Number density of ionized H
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

      epsx = pmdpt/dmdpt/(5./3. - 1.) ;Specific internal energy
      depsx = shift(epsx,-1) - epsx    ;Finite difference for derivative
      advheatlhsx = dmdpt * vxmdpt * depsx/dx            ;Internal energy contribution
      advheatrhsx = pmdpt/dmdpt * vxmdpt * (shift(dmdpt,-1) - dmdpt)/dx;PdV contribution - Following sign convention of Ruth -
      advionizx = dmdpt/m_h * vxmdpt * (shift(ifrac,-1) - ifrac)/dx ;Advection effect on ionization


      ;Y DIRECTION TERMS
      if (offa>0) then begin
         restore, 'ioniz_sphere.'+string(fileno, format='(I04)')+'_a_offmdptvals.sav'
      endif else begin
         restore, 'ioniz_sphere.'+string(fileno, format='(I04)')+'_offmdptvals.sav'
      endelse

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


      ;Ionization fraction changes
      dfx = (shift(ifrac,-1) - ifrac)
      dfy = (ifrac1-ifrac)
      dfz = (ifrac2-ifrac)
;      print, 'I have df. On Lev0, x:', dfx, ' y:', dfy, ' z:', dfz

      ;Kill edge effects
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

      ;Commented out because it was only the lowest level of resolution
      ;; if(keyword_set(ps)) then begin
      ;;    device, filen='~/Downloads/advheat_'+string(fileno, format='(I04)')+'.eps', ysize=10.85, xsize=8, /inches, /encapsulated
      ;;    !p.charsize = 1.5
      ;;    !p.charthick=4
      ;;    !x.thick=4
      ;;    !y.thick=4
      ;; endif

      ;; !p.multi=[0,2,4]
      ;; ;Plot internal energy and PdV components
      ;; plot, xmdpt/rp, advheatrhs, /nodata, title='Total', ytitle='Heating-Cooling Rates [erg cm!e-3!n s!e-1!n]'
      ;; oplot, xmdpt/rp, advheatlhs, color=fsc_color("red")
      ;; oplot, xmdpt/rp, advheatrhs, color=fsc_color("blue"), linestyle=2

      ;; plot, xmdpt/rp, advioniz, title='Total',  ytitle='Ionization-Recombination Rates [cm!e-3!n s!e-1!n]'

      ;; plot, xmdpt/rp, advheatlhsx, xsty=1, ysty=1, title='X component';, yra=[-1.5d-6, 1.5d-6]
      ;; oplot, xmdpt/rp, advheatrhsx, linestyle=2

      ;; plot, xmdpt/rp, advionizx, title='X component'

      ;; plot, xmdpt/rp, advheatlhsy, title='Y component'
      ;; oplot, xmdpt/rp, advheatrhsy, linestyle=2

      ;; plot, xmdpt/rp, advionizy, title='Y component'

      ;; plot, xmdpt/rp, advheatlhsz, title='Z component'
      ;; oplot, xmdpt/rp, advheatrhsz, linestyle=2

      ;; plot, xmdpt/rp, advionizz, title='Z component'


      ;; if(keyword_set(ps)) then begin
      ;;    device, /close
      ;; endif

      ;Plot heating/cooling rates
;      !p = p5  & !x= x5  & !y = y5

;stop
      if(keyword_set(ps)) then begin
         if (offa>0) then begin
            device, filen='~/Downloads/heat_a_'+string(fileno, format='(I04)')+'.eps', xsize=10.85, ysize=8, /inches, /encapsulated
         endif else begin
            device, filen='~/Downloads/heat_'+string(fileno, format='(I04)')+'.eps', xsize=10.85, ysize=8, /inches, /encapsulated
         endelse
         !p.charsize = 1.1
         !p.charthick=4         
         !p.thick=8
         !x.thick=5
         !y.thick=5
      endif
      
      !p.multi=[0,2,2]

      ;Bounds setup                                   
      alimamt = max([(bounds[1]-dx), rp]) ;;where(xmdpt/rp gt -20 AND xmdpt/rp lt 20)
      alevlimn = where(xmdpt lt - alimamt)
      alevlimp = where(xmdpt gt alimamt)
      alt0 = where (xmdpt lt 0)
      agt0 = where (xmdpt gt 0)
      outsidea = alevlimn;where ((xmdpt ) ge rp)
      outsideb = alevlimp;where ((xmdpt ) le rp)


      ;Physical param
      rp = 1.5e10

      ;PLOTS
      ;Plot heating/cooling rates

;      plot,  xmdpt/rp, rate_ionizheat, xsty =1, ysty =1, xtitle='x/R!dp!n', ytitle='Heating Rates [erg cm!e-3!n s!e-1!n]', /nodata, xmargin=[15,5], yra=[1d-12, 4d-7], xra=[-3, 3], /ylog
;   minval = 10.^(!y.crange[0])+1e-20
;   maxval = 10.^(!y.crange[1])-1e-13

;      polyfill, [-1,1,1,-1], [minval,minval, maxval, maxval], color=fsc_color("light gray"), /data

      plot,  xmdpt/rp, rate_ionizheat, xsty =1, ysty =1, xtitle='x/R!dp!n', ytitle='dE/dt [erg cm!e-3!n s!e-1!n]', /nodata, xmargin=[18,5], yra=[-3.3d-7, 3.3d-7], xra=[-2, -1] ;, /ylog
;yra=[-2.5d-7, 2.5d-7];offaxis
      minval = !y.crange[0]+1e-9
      maxval = !y.crange[1]-1e-9

      oplot, xmdpt(outsidea)/rp, (rate_ionizheat(outsidea) + rate_recombcool(outsidea) + rate_lyacool(outsidea) + advheatrhs(outsidea)), thick =14
      oplot, xmdpt[alevlimn]/rp, rate_ionizheat[alevlimn], color=fsc_color("red6", /brewer), thick=10
      oplot, xmdpt[alevlimn]/rp, (rate_lyacool[alevlimn]), color=fsc_color("grn8",/brewer), thick=10
      oplot, xmdpt[alevlimn]/rp, (rate_recombcool[alevlimn]), color=fsc_color("ygb4",/brewer), thick=10
      oplot, xmdpt[alevlimn]/rp, (advheatlhs[alevlimn]), color=fsc_color("org3",/brewer), linestyle=5, thick = 10
      oplot, xmdpt[alevlimn]/rp, (advheatrhs[alevlimn]), color=fsc_color("blu7",/brewer), linestyle=5, thick = 10
      ;; oplot, xmdpt[alevlimn]/rp, (advheatlhsx[alevlimn]), color=fsc_color("gray")
      ;; oplot, xmdpt[alevlimn]/rp, (advheatlhsy[alevlimn]), color=fsc_color("gray"), linestyle=1
      ;; oplot, xmdpt[alevlimn]/rp, (advheatlhsz[alevlimn]), color=fsc_color("gray"), linestyle=2


      xyouts, !x.crange[0]+0.14, 2.8e-7, 'Photoionization', color=fsc_color("red6", /brewer), /data, charsize=1.0
      xyouts, !x.crange[0]+0.14, 2.4e-7, 'Lyman alpha', color=fsc_color("grn8", /brewer), /data, charsize=1.0
      xyouts, !x.crange[0]+0.14, 2.0e-7, 'Recombination', color=fsc_color("ygb4", /brewer), /data, charsize=1.0
      xyouts, !x.crange[0]+0.14, 1.6e-7, 'PdV', color=fsc_color("blu7", /brewer), /data, charsize=1.0
      xyouts, !x.crange[0]+0.14, 1.2e-7, 'Internal energy', color=fsc_color("org3", /brewer), /data, charsize=1.0
      xyouts, !x.crange[0]+0.14, 0.8e-7, 'Total-internal energy', /data, charsize=1.0

      oplot, [!x.crange[0]+0.05,!x.crange[0]+0.12], [2.9e-7, 2.9e-7], color=fsc_color("red6", /brewer), thick=7
      oplot, [!x.crange[0]+0.05,!x.crange[0]+0.12], [2.5e-7, 2.5e-7], color=fsc_color("grn8",/brewer), thick=7
      oplot, [!x.crange[0]+0.05,!x.crange[0]+0.12], [2.1e-7, 2.1e-7], color=fsc_color("ygb4",/brewer), thick=7
      oplot, [!x.crange[0]+0.05,!x.crange[0]+0.12], [1.7e-7, 1.7e-7], color=fsc_color("blu7",/brewer), linestyle=2, thick = 7
      oplot, [!x.crange[0]+0.05,!x.crange[0]+0.12], [1.3e-7, 1.3e-7], color=fsc_color("org3",/brewer), linestyle=2, thick = 7
      oplot, [!x.crange[0]+0.05,!x.crange[0]+0.12], [0.9e-7, 0.9e-7], thick=12


      p6 = !P & x6 = !X & y6 = !Y

      plot,  xmdpt/rp, rate_ionizheat, xsty =1, ysty =1, xtitle='x/R!dp!n', ytitle='dE/dt [erg cm!e-3!n s!e-1!n]', /nodata, xmargin=[15,5], yra=[-5d-9, 5d-9], xra=[1, 3] ;, /ylog
;yra=[-5d-9, 5d-9];offaxis
      minval = !y.crange[0]+1e-11
      maxval = !y.crange[1]-1e-11
      atmedge=1.5703/1.5 ;Real value is 1.562, but that's 5.3 cells in from 0.
      oplot, xmdpt(outsideb)/rp, (rate_ionizheat(outsideb) + rate_recombcool(outsideb) + rate_lyacool(outsideb) + advheatrhs(outsideb)), thick=14
      oplot, xmdpt[alevlimp]/rp, rate_ionizheat[alevlimp], color=fsc_color("red6"), thick=10
      oplot, xmdpt[alevlimp]/rp, (rate_recombcool[alevlimp]), color=fsc_color("ygb4"), thick=10
      oplot, xmdpt[alevlimp]/rp, (rate_lyacool[alevlimp]), color=fsc_color("grn8"), thick=10
      oplot, xmdpt[alevlimp]/rp, (advheatlhs[alevlimp]), color=fsc_color("org3"), linestyle=5, thick = 10
      oplot, xmdpt[alevlimp]/rp, (advheatrhs[alevlimp]), color=fsc_color("blu7"), linestyle=5, thick = 10

      ;; oplot, xmdpt[alevlimp]/rp, (advheatlhsx[alevlimp]), color=fsc_color("gray")
      ;; oplot, xmdpt[alevlimp]/rp, (advheatlhsy[alevlimp]), color=fsc_color("gray"), linestyle=1
      ;; oplot, xmdpt[alevlimp]/rp, (advheatlhsz[alevlimp]),
      ;; color=fsc_color("gray"), linestyle=2
      polyfill, [1,atmedge,atmedge,1], [minval,minval, maxval, maxval], color=fsc_color("light gray"), /data

      p7 = !P & x7 = !X & y7 = !Y

;      al_legend,['Photoioniz','Recomb','Lya', 'Int energy',
;      'PdV'],linestyle=[0,0,0,0,0],color=[fsc_color("red"),
;      fsc_color("blue"), fsc_color("teal"), fsc_color("orange"),
;      fsc_color("green")]
;; print, 'Here here' 
;; stop
;; fmt = '(E10.3, E10.3, E10.3, E10.3, E10.3)'
;; forprint, F=fmt, xmdpt[0:78], rate_ionizheat[0:78], rate_recombcool[0:78], advheatlhs, advheatrhs, text='~/Downloads/someratesforsimeon_lev0.dat'



      plot, xmdpt/rp, rate_ioniz, xsty =1, ysty =1, xtitle='x/R!dp!n', ytitle='Ioniz. Rates [cm!e-3!n s!e-1!n]', /nodata, xmargin=[18,5], yra=[-3d4, 1.1d5], ymargin=[8,2], xra=[-2,-1];, /ylog
;yra=[-1d4, 7d4]
      minval = !y.crange[0]+1e-9
      maxval = !y.crange[1]-1e-9

;      minval = 10.^(!y.crange[0])+1e-20
 ;     maxval = 10.^(!y.crange[1])-1e-13
;      polyfill, [-1,1,1,-1], [minval,minval, maxval, maxval], color=fsc_color("light gray"), /data

      oplot, xmdpt[alevlimn]/rp, rate_ioniz[alevlimn], color=fsc_color("red6"), thick=10
      oplot, xmdpt[alevlimn]/rp, (rate_recomb[alevlimn]), color=fsc_color("ygb4"), thick=10
      oplot, xmdpt[alevlimn]/rp, (advioniz[alevlimn]), color=fsc_color("orange"), thick=10, linestyle=5
      ;; oplot, xmdpt[alevlimn]/rp, (advionizx[alevlimn]), color=fsc_color("gray")
      ;; oplot, xmdpt[alevlimn]/rp, (advionizy[alevlimn]), color=fsc_color("gray"), linestyle=1
      ;; oplot, xmdpt[alevlimn]/rp, (advionizz[alevlimn]), color=fsc_color("gray"), linestyle=2
;      oplot, xmdpt[alevlimn]/rp,
;      (rate_ioniz[alevlimn]+rate_recomb[alevlimn]), linestyle=3

      xyouts, !x.crange[0]+0.14, 9.5e4, 'Photoionization', color=fsc_color("red6", /brewer), /data, charsize=1.0
      xyouts, !x.crange[0]+0.14, 8.5e4, 'Recombination', color=fsc_color("ygb4", /brewer), /data, charsize=1.0
      xyouts, !x.crange[0]+0.14, 7.5e4, 'Advection', color=fsc_color("orange"), /data, charsize=1.0


      oplot, [!x.crange[0]+0.05,!x.crange[0]+0.12], [1.e5, 1.e5], color=fsc_color("red6", /brewer), thick=7
      oplot, [!x.crange[0]+0.05,!x.crange[0]+0.12], [9e4, 9e4], color=fsc_color("ygb4",/brewer), thick=7
      oplot, [!x.crange[0]+0.05,!x.crange[0]+0.12], [8e4, 8e4], color=fsc_color("orange"), linestyle=2, thick = 7




      p8 = !P & x8 = !X & y8 = !Y


      plot, xmdpt/rp, rate_ioniz, xsty =1, ysty =1, xtitle='x/R!dp!n', ytitle='Ioniz. Rates [cm!e-3!n s!e-1!n]', /nodata, xmargin=[15,5], yra=[-3.5d2, 5d1], ymargin=[8,2], xra=[1,3]
      minval = !y.crange[0]+1e-11
      maxval = !y.crange[1]-1e-11
      oplot, xmdpt[alevlimp]/rp, rate_ioniz[alevlimp], color=fsc_color("red6"), thick=10
      oplot, xmdpt[alevlimp]/rp, (rate_recomb[alevlimp]), color=fsc_color("ygb4"), thick=10
      oplot, xmdpt[alevlimp]/rp, (advioniz[alevlimp]), color=fsc_color("orange"), thick=10, linestyle=5
      ;; oplot, xmdpt[alevlimp]/rp, (advionizx[alevlimp]), color=fsc_color("gray")
      ;; oplot, xmdpt[alevlimp]/rp, (advionizy[alevlimp]), color=fsc_color("gray"), linestyle=1
      ;; oplot, xmdpt[alevlimp]/rp, (advionizz[alevlimp]), color=fsc_color("gray"), linestyle=2
;      oplot, xmdpt[alevlimp]/rp,
;      (rate_ioniz[alevlimp]+rate_recomb[alevlimp]), linestyle=3
      polyfill, [1,atmedge,atmedge,1], [minval,minval, maxval, maxval], color=fsc_color("light gray"), /data

      p9 = !P & x9 = !X & y9 = !Y

      
;      al_legend,['Photoioniz','Recomb','Advection'],color=[fsc_color("red"), fsc_color("blue"), fsc_color("teal"), fsc_color("orange")]


      ;; oplot, xmdpt[alevlimn]/rp, abs(rate_recomb[alevlimn]), color=fsc_color("blue")
      ;; oplot, xmdpt[alevlimp]/rp, abs(rate_recomb[alevlimp]), color=fsc_color("blue")
      ;; oplot, xmdpt[alevlimn]/rp, abs(advioniz[alevlimn]), color=fsc_color("orange")
      ;; oplot, xmdpt[alevlimp]/rp, abs(advioniz[alevlimp]), color=fsc_color("orange")

      ;; oplot, xmdpt[alevlimn]/rp, abs(advionizx[alevlimn]), color=fsc_color("cyan")
      ;; oplot, xmdpt[alevlimn]/rp, abs(advionizy[alevlimn]), color=fsc_color("magenta");, linestyle=2
      ;; oplot, xmdpt[alevlimn]/rp, abs(advionizz[alevlimn]), color=fsc_color("yellow");, linestyle=5
      ;; oplot, xmdpt[alevlimp]/rp, abs(advionizx[alevlimp]), color=fsc_color("cyan")
      ;; oplot, xmdpt[alevlimp]/rp, abs(advionizy[alevlimp]), color=fsc_color("magenta");, linestyle=2
      ;; oplot, xmdpt[alevlimp]/rp, abs(advionizz[alevlimp]), color=fsc_color("yellow");, linestyle=5




   for ilev = 1, 3 do begin
      if (keyword_set(yplot)) then begin
         filen = 'ioniz_sphere-lev'+strtrim(ilev,2)+'.'+string(fileno, format='(I04)')+'y_mdptvals.sav'
      endif else begin
         if (offa > 0) then begin
            filen = 'ioniz_sphere-lev'+strtrim(ilev,2)+'.'+string(fileno, format='(I04)')+'_a_mdptvals.sav'
         endif else begin
            filen = 'ioniz_sphere-lev'+strtrim(ilev,2)+'.'+string(fileno, format='(I04)')+'_mdptvals.sav'
         endelse
      endelse

      plotdvt, filen, ilev, 0, 1
   endfor


if(keyword_set(ps)) then begin
   device, /close
endif

stop


END
