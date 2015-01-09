;; FUNCTION blackbody, freq, temp
;;   RETURN 2.*h*freq^3./(c^2.) / (exp(h*freq/k/temp)-1.)
;; END

PRO simvoigt_valt, basename, fromsaved=fromsaved, plotobsc=plotobsc, plotxsec=plotxsec, ps=ps, cmap=cmap, plotaveobs=plotaveobs, coldens=coldens
  plotsym, 0, .4, /fill
  if (keyword_set(ps)) then begin
     !p.font=0
     !p.charsize=1.4
     !p.charthick=7
     !x.charsize=1.3
     !y.charsize=1.3
     !p.thick=8
     !x.thick=5
     !y.thick=5
     set_plot, 'ps'
  endif else begin
     set_plot, 'x'
  endelse

if (~keyword_set(fromsaved)) then begin
  densfile = basename+'_density.txt'
  pressurefile = basename+'_pressure.txt'
  posfile = basename+'_pos.txt'
  velfile = basename+'_velocity.txt'
  scalarfile = basename+'_specific_scalar_0_.txt'

  READCOL, posfile, format='D D D', xvec, yvec, zvec
  READCOL, velfile, format='D D D', vxvec, vyvec, vzvec
  READCOL, densfile, format='D', densvec
  READCOL, pressurefile, format='D', pressurevec
  READCOL, scalarfile, format='D', dratvec

  n = n_elements(densvec)^(1/3.)
  x = reform(xvec, n, n, n)
  y = reform(yvec, n, n, n)
  z = reform(zvec, n, n, n)
  vx = reform(vxvec, n, n, n)
  vy = reform(vyvec, n, n, n)
  vz = reform(vzvec, n, n, n)
  d = reform(densvec, n, n, n)
  pressure = reform(pressurevec, n, n, n)
  dratvec = dratvec-1d-4      ;Account for the fact that min
  fixdrat = where (dratvec lt 0)
  dratvec(fixdrat) = 0.
;  neutral frac is 1e-4, which is the background level
  drat = reform(dratvec, n, n, n)
   
  save, /variables, filename=basename+'_sv.sav'

endif else begin
   restore, basename+'_sv.sav'
endelse

;x =x + max(x)

  ;; Constants
   c = 2.99792458d10            ; Speed of light in cm/s from NIST CODATA - 9/2014
   h = 6.62606957d-27           ; Planck's constant in erg*s from NIST CODATA - 9/2014
   k = 1.3806488d-16            ; Boltzmann constant in erg/K from NIST CODATA - 9/2014
   qe = 4.803205d-10            ; Elementary charge from Wolfram alpha - 9/2014
   me = 9.10938291e-28          ; Electon mass in g from NIST CODATA - 9/2014

   m_H = 1.67d-24
   mu = 1.67d-24
   alpha_C = 1.0d-10



   ;;Frequency range of spectrum
   lambda0 = 1215.67d-8         ; cm based on Morton's value - rounded
   nu0 = c/lambda0

   nu_h = 2.471d15              ;1213 angstroms in Hz
   nu_l = 2.461d15              ;1218 angstroms in Hz

   nfreq = 200
   nu = dindgen(nfreq)*(nu_h-nu_l)/(nfreq-1)*1d0 + nu_l
   lambda = c / nu

   ;;Unnecessary step 1: Calculate the unobscured stellar blackbody spectrum
    Tstar = 6075                ; Effective temperature of HD209458 from exoplanet.eu on 8/28/14
    Bnustar = 2.*h*nu^3./(c^2.) / (exp(h*nu/k/Tstar)-1.) ;Blackbody spectrum for star

    rsun = 6.958d10             ;Solar radius
    rsun2 = rsun*rsun

   ;;Compute temperature & ionization fraction, as done in the code (ionrad_3d.c)
   dneutral = drat*d            ;Convert back to neutral density
   n_H = dneutral / m_H
   n_Hplus = (d - dneutral) / m_H
   n_e = n_Hplus + d * alpha_C / (14.0 * m_H) 
   ifrac = n_e / (n_H + n_Hplus)
   mu_temp = ifrac*0.5*m_H + (1.0-ifrac)*mu
   T = pressure/d* mu_temp/k

   nucorr = 1+(vx/c)
   lambdacorr = 1-(vx/c)

   ;;Calculate line profile  - Voigt
   gamma = 6.265d8              ;Einstein A for Lya in Hz from Morton, D. , 2003 - Atomic Data for Resonance Lines III
   fosc = 0.4164                ;Oscillator strength from Morton (2003) for gl = 2, gu = 6   
   dx = abs(xvec[1]-xvec[0])

   delnud = 1./lambda0*sqrt(2.*k*T/mu_temp) ; AT: Is mu_temp correct?  And T?
   a = gamma/4./!DPI/delnud

;   u = dblarr(n,n,n)
   sigmanu = dblarr(n,n,n,nfreq)
   tau = dblarr(n,n, nfreq)
   etau = dblarr(n,n, nfreq)
   etauave = dblarr(nfreq)

   ;;Calculate quantities at each location, as you loop over frequencies
   for i = 0, nfreq-1 do begin
      u = (nu[i] - nu0*nucorr)/delnud ;There's an indexing problem here!!
;      u = (nu[i]*nucorr - nu0*nucorr)/delnud ;There's an indexing problem here!!
      phinu = voigt(a,u) / sqrt(!dpi) / delnud ;Voigt profile
;qphi = where (finite(phinu, /nan))
      sigmanu[*,*,*,i] = !DPI * qe^2./me/c* fosc * phinu ;Cross section
;qsigma = where (finite(sigmanu, /nan))

      etautot=0
      counter = 0

      ;;Calculate optical depth along each ray
      ;Here we cast rays along the x direction
      for k = 0, n-1 do begin
         for j = 0, n-1 do begin
            tau[j,k,i] = total(n_H[*,j,k]*dx * sigmanu[*,j,k,i])
            etau[j,k,i] = exp(-tau[j,k,i])

            r2 = y[0,j,0]^2 + z[0,0,k]^2
            if (r2 le rsun2) then begin
               etautot = etautot + etau[j,k,i]
               counter=counter+1
            endif 
         endfor
      endfor
      etauave[i] = etautot/counter
      
   endfor


   ;Create channel map
   if (keyword_set(cmap)) then begin
;writefits, 'mytau3.fits', tau[*,*,findex]
      lcenter = where(abs(lambda-lambda0) eq min(abs(lambda-lambda0)))
      vlines=[0, 20, 25, 30,45]
;100,200, 400]
      nuvlines = nu0*(1+vlines*1d5/c)
      findex = dblarr(n_elements(vlines))
      rp = 1.5e10

      cgLoadCT, 9, /BREWER, RGB_TABLE=palette
      !p.multi=[0,n_elements(vlines),1]

      if (keyword_set(ps)) then begin
         set_plot, 'ps'
         device, filen='/Users/anjalitripathi/Downloads/'+basename+'_cmap_spacing_rsun.eps', xsize=9.5, ysize=5, /inches, /encapsulated
         !p.thick=6
      endif else begin
         set_plot, 'x'
      endelse

      for ff = 0, n_elements(vlines)-1 do begin
         findex[ff] = where (abs(nu-nuvlines[ff]) eq min(abs(nu-nuvlines[ff])))
         print, ff, ' v:', vlines[ff], ' index:', findex[ff]
                                ;reform(y[0,*,0],n), yvector=reform(z[0,0,*],n),
         cgimage,1-etau[*,*,findex[ff]], xrange=[min(y), max(y)]/rp, yrange=[min(z),max(z)]/rp, PALETTE=palette, /keep_aspect, /axes, minvalue = 0, maxvalue=1, title=strtrim(vlines[ff],2)+'km/s', xtitle='Position [Rp]' ;, stretch='Log'
;         cgcontour, 1-etau[*,*,findex[ff]], levels=[1-exp(-1)],
;         label=0, /onimage, c_thick=3 ;Creat tau=1 contour
         tvcircle, rsun/rp, 0, 0, linestyle=2, /data
      endfor
;ysize = 5
cgcolorbar, PALETTE=palette, range=[0,1], title='Obscuration', tlocation='top', tcharsize =.9


      if (keyword_set(ps)) then begin
         device, /close
      endif
      stop
   endif

   ;Create channel map
   if (keyword_set(coldens)) then begin
      rp = 1.5e10
      coldensx = dblarr(n,n)
      coldensy = dblarr(n,n)
      coldensz = dblarr(n,n)
      for k = 0, n-1 do begin
         for j = 0, n-1 do begin
            coldensx[j,k] = total(n_H[*,j,k]*dx)
            coldensy[j,k] = total(n_H[k,*,j]*dx)
            coldensz[j,k] = total(n_H[j,k,*]*dx)
         endfor
      endfor

;; openw,lun2,'coldensx.txt', /get_lun
;; printf, lun2, coldensx,FORMAT='(G10.4)' 
;; ;, coldensy, coldensz
;; close, lun2
;; free_lun, lun2
;; stop

      cgLoadCT, 9, /BREWER, RGB_TABLE=palette
;      !p.multi=[0,3,1]
!p.multi=0
      if (keyword_set(ps)) then begin
         set_plot, 'ps'
         device, filen='/Users/anjalitripathi/Downloads/'+basename+'_coldens.eps', ysize=8.5, xsize=4, /inches, /encapsulated
;         device, filen='/Users/anjalitripathi/Downloads/'+basename+'_coldenscolorbar.eps', ysize=8.5, xsize=4, /inches, /encapsulated
         !p.thick=6
      endif else begin
         set_plot, 'x'
      endelse

      ;; plotsym, 0, 0.75, /fill
      ;; plot, y/rp, coldensx, psym=8, /ylog, ysty=2, ytitle='X Column density', xtitle='Position (Rp)'
      ;; plot, x/rp, coldensy, psym=8, /ylog, ysty=2, ytitle='Y Column density', xtitle='Position (Rp)'
      ;; plot, x/rp, coldensz, psym=8, /ylog, ysty=2, ytitle='Z Column density', xtitle='Position (Rp)'



      cgimage,bytscl(alog10(coldensx/max(coldensx))), xrange=[min(y), max(y)]/rp, yrange=[min(z),max(z)]/rp, PALETTE=palette, /keep_aspect, /axes,  title='Coldens integ along x', xtitle='Y [Rp]', ytitle='Z [Rp]' ;,  minvalue=-8, maxvalue=1;,stretch='Log'
      cgimage,bytscl(alog10(coldensy/max(coldensy))), xrange=[min(y), max(y)]/rp, yrange=[min(z),max(z)]/rp, PALETTE=palette, /keep_aspect, /axes,  title='Coldens integ along y', xtitle='Z [Rp]', ytitle='X [Rp]' ;,  minvalue=-8, maxvalue=1
      cgimage,bytscl(alog10(coldensz/max(coldensz))), xrange=[min(y), max(y)]/rp, yrange=[min(z),max(z)]/rp, PALETTE=palette, /keep_aspect, /axes,  title='Coldens integ along z', xtitle='X [Rp]', ytitle='Y [Rp]' ;,  minvalue=-8, maxvalue=1
;cgcolorbar, PALETTE=palette, range=[1e-8,1], title='Normalized column density', tcharsize =.9, /ylog, yticks=0, /vertical

;ysize = 5



      if (keyword_set(ps)) then begin
         device, /close
      endif
      stop
   endif



   ;Total (averaged) obscuration
   if (keyword_set(plotaveobs)) then begin
      if (keyword_set(ps)) then begin
         set_plot, 'ps'
         device, filen=basename+'_aveobs3.eps', xsize=10.5, ysize=7.5, /inches, /encapsulated
      endif

      !p.multi=0
      plot, lambda* 1d8, 1-etauave, xstyle = 9, ystyle =2, /ylog, xtitle='Wavelength [angstroms]', ytitle='Average obscured fraction', xrange=[1213.44, 1217.90], xmargin=[12,3], ymargin= [5,5]
      oplot, ([lambda0,lambda0])*1d8, [1e-20,10], color=fsc_color("gray"), linestyle = 2
      axis, xaxis = 1, xrange=-c*(lambda0*1d8/!x.crange -1.) / 1d5, xtitle = 'Velocity [km s!e-1!n]', xstyle=1 ;charsize = 1.3,
;The velocity axis has been flipped 


      if (keyword_set(ps)) then begin
         device, /close
      endif
      stop
   endif








!p.multi = [0,5,5]
if (keyword_set(plotxsec)) then begin
   xsel = 9
   if (keyword_set(ps)) then begin
      device, filen=basename+'_simxsec_'+strtrim(xsel,2)+'.eps', xsize=10.5, ysize=7, /inches, /encapsulated
   endif
;print, x[xsel] - 0.5*max(x)
endif
if (keyword_set(plotobsc)) then begin
   device, filen=basename+'_simobsc.eps', xsize=10.5, ysize=7, /inches, /encapsulated
endif

for ay = 0, 4 do begin
   for az = 0, 4 do begin
      yy = ay*20 -1 
      if (yy lt 0) then yy=0;

      zz = az*20 -1 
      if (zz lt 0) then zz=0;
;sigmanu[39, yy, zz, *]
;1-etau[yy, zz, *]
;sigmanu[xsel, yy, zz, *]
;v='+strtrim(vx[xsel],2)
;+'Tcent='+strtrim(round(T[39,yy,zz]),2)+'K
      
      if (keyword_set(plotobsc)) then begin
         plot, (lambda-lambda0)* 1d8, 1-etau[yy, zz, *], ystyle =1, /ylog, xtitle='Lambda-Lambda0 [angstrom]', ytitle='Obscuration', title = 'Tcent='+string(T[39,yy,zz],format='(E9.2)')+'K'+' ('+strtrim(yy,2)+','+strtrim(zz,2)+')', yra=[1d-7, 4], xstyle=1
;, /nodata
;         oplot, (lambda-lambda0) * 1d8, 1-etau[yy, zz, *], color=fsc_color("gray"), psym=8
         oplot, ([lambda0,lambda0]-lambda0)*1d8, [1e-20,10], color=fsc_color("rosy brown"), linestyle = 2
      endif

      if (keyword_set(plotxsec)) then begin
         plot, (lambda)* 1d8, sigmanu[xsel, yy, zz, *], xstyle=1,ystyle =1, /ylog, xtitle='Wavelength [angstrom]', ytitle='X-sec [cm^2]', title = 'v='+string(vx[xsel, yy, zz], format='(E9.2)')+'cm/s'+' ('+strtrim(xsel,2)+','+strtrim(yy,2)+','+strtrim(zz,2)+')', xra=[lambda0*1d8-1, lambda0*1d8+1],xticks=2;, /nodata
         oplot, lambda * 1d8, sigmanu[xsel, yy, zz, *], color=fsc_color("gray"), psym=8
         oplot, [lambda0,lambda0]*1d8, [1e-23,10], color=fsc_color("rosy brown"),linestyle = 2
         ;; plot, (lambda)*lambdacorr[xsel,yy,zz]* 1d8, sigmanu[xsel, yy, zz, *], xstyle=1,ystyle =1, /ylog, xtitle='Wavelength [angstrom]', ytitle='X-sec [cm^2]', title = 'v='+string(vx[xsel, yy, zz], format='(E9.2)')+'cm/s'+' ('+strtrim(xsel,2)+','+strtrim(yy,2)+','+strtrim(zz,2)+')', xra=[lambda0*1d8-1, lambda0*1d8+1],xticks=2;, /nodata
         ;; oplot, lambda *lambdacorr[xsel,yy,zz]* 1d8, sigmanu[xsel, yy, zz, *], color=fsc_color("gray"), psym=8
         ;; oplot, [lambda0,lambda0]*1d8, [1e-23,10], color=fsc_color("rosy brown"),linestyle = 2
      endif

   endfor
endfor

   if(keyword_set(ps)) then begin
      device, /close
   endif


   ;; ;Solid angle
   ;; levfactor = 1.               ;Assuming we're on the lowest resolution level
   ;; solidangle = 1 *levfactor               
   

   ;; Fnu[i] = total(bnu*solidangle)






  
;;    rp = 1.5e10
;;    for i=0,n-1 do begin   
;;       cdens[i] = total(d[*, i, floor(n/2)-1]*drat[*, i, floor(n/2)-1]* dx); - total(d[*, 0, floor(n/2)-1]* dx)
;;       if (i lt floor(n/2)-1) then begin
;;          b[i] = sqrt(y[0,i,0]^2 + z[0,0,n/2 - 1]^2)*(-1)/rp
;;       endif else begin
;;          b[i] = sqrt(y[0,i,0]^2 + z[0,0,n/2 - 1]^2)/rp
;;       endelse
;;    endfor

;;    cdens = cdens / 1.67d-24
;;    save, /variables, filename='coldens_'+basename+'.sav'
  STOP
END

