; GRACES pipeline, by Andre-Nicolas Chene

; Main script is dg.pro
; It can be found at the end, as the last pro, after all the functions.


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;  FUNCTION APPLYING THE OVERSCAN CORRECTION (nice, but not necessary)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function overs_corr,im,h
s=size(im)
len=s[2] ;numer of pixels in one row

;Two approches depending on if it was observed in 1 or 2amp mode
if strcmp(strcompress(sxpar(h,'AMPLIST'),/remove),'a') then begin
  ;Gets the biassec and ccdsec from the header
  bsec=sxpar(h,'BIASSEC')
  bsec=strmid(bsec,1,strpos(bsec,',')-1)
  bsec1=strmid(bsec,0,strpos(bsec,':'))-1
  bsec2=strmid(bsec,strpos(bsec,':')+1)-1
  csec=sxpar(h,'CCDSEC')
  csec=strmid(csec,1,strpos(csec,',')-1)
  csec1=strmid(csec,0,strpos(csec,':'))-1
  csec2=strmid(csec,strpos(csec,':')+1)-1
  
  ;Gets the overscan value for each row
  vec_overs=fltarr(len)
  for i=0,len-1 do vec_overs[i]=mean(im[bsec1+1:bsec2,i])
  
  ;Fits a slope to the overscan values as a function of row #
  res=ladfit(findgen(len),vec_overs)
  
  ;Creates the overscan corrected resulting image
  im_corr=fltarr(csec2+1,len)
  for i=0,len-1 do im_corr[0:csec2,i]=im[csec1+8:csec2+8,i]-(res[0]+res[1]*i) ;the first 8 pixels have no information
endif else begin
  ;Gets the biassec and ccdsec from the header
  bsecA=sxpar(h,'BSECA')
  bsecA=strmid(bsecA,1,strpos(bsecA,',')-1)
  bsecA1=strmid(bsecA,0,strpos(bsecA,':'))-1
  bsecA2=strmid(bsecA,strpos(bsecA,':')+1)-1
  csecA=sxpar(h,'CSECA')
  csecA=strmid(csecA,1,strpos(csecA,',')-1)
  csecA1=strmid(csecA,0,strpos(csecA,':'))-1
  csecA2=strmid(csecA,strpos(csecA,':')+1)-1
  
  bsecB=sxpar(h,'BSECB')
  bsecB=strmid(bsecB,1,strpos(bsecB,',')-1)
  bsecB1=strmid(bsecB,0,strpos(bsecB,':'))-1
  bsecB2=strmid(bsecB,strpos(bsecB,':')+1)-1
  csecB=sxpar(h,'CSECB')
  csecB=strmid(csecB,1,strpos(csecB,',')-1)
  csecB1=strmid(csecB,0,strpos(csecB,':'))-1
  csecB2=strmid(csecB,strpos(csecB,':')+1)-1
  
  ;Gets the overscan value for each row
  vec_oversA=fltarr(len)
  vec_oversB=fltarr(len)
  for i=0,len-1 do begin
    ; for the A part
    vec_oversA[i]=mean(im[bsecA1:bsecA2-1,i])
    ; for the B part
    vec_oversB[i]=mean(im[bsecB1+1:bsecB2,i])
  endfor
  
  ;Fits a slope to the overscan values as a function of row #
  resA=ladfit(findgen(len),vec_oversA)
  resB=ladfit(findgen(len),vec_oversB)
  
  ;Creates the overscan corrected resulting image
  im_corr=fltarr(csecB2-bsecA2,len)
  for i=0,len-1 do begin
    im_corr[0:csecA2-bsecA2-1,i]=im[csecA1:csecA2,i]-(resA[0]+resA[1]*i)
    im_corr[csecB1-bsecA2-1:csecB2-bsecA2-1,i]=im[csecB1:csecB2,i]-(resB[0]+resB[1]*i)
  endfor
endelse

return,im_corr
end



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;  FUNCTION FINDING THE TRACES
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function find_trace,flat,spmd,ts,wdt,vcen,nco

;;;;;;;;;;;;;;;;
;Useful values
s=size(flat)         ;s contains dimensions x and y of the flat frame
len=s[2]             ;len is the number of pix in y in the flat frame
cutmh=flat[*,len/2]  ;trace is a cut at mid-height of the flat frame
nord=n_elements(vcen);number of orders to extract

;;;;;;;;;;;;;;;;
;Model (4 gaussian curves) of the sliced image (ref)
a=findgen(wdt) ;x-vector the size of the aperture
if strcmp(spmd,'1f') then begin
  ref=.7*exp(-1*((a-3.5)^2/(2*2.^2)))+exp(-1*((a-10.5)^2/(2*2.^2)))+exp(-1*((a-17.5)^2/(2*2.^2)))+.5*exp(-1*((a-24.5)^2/(2*2.^2)))
endif else begin
  ref=.7*exp(-1*((a-4.0)^2/(2*2.^2)))+exp(-1*((a-10.0)^2/(2*2.^2)))+.9*exp(-1*((a-20.5)^2/(2*2.^2)))+exp(-1*((a-27.0)^2/(2*2.^2)))
  vcen=vcen-6 ;slight shift measured in the position of the order with respect to 1-fiber mode
endelse

;first crude re-alignement to find the centroid of the traces (as they may move over time on the detector)
sp_ref=fltarr(n_elements(cutmh))
for i=0,n_elements(vcen)-1 do sp_ref[vcen[i]-n_elements(ref)/2:vcen[i]+n_elements(ref)/2]=ref
lag=findgen(101)-50 ;x-vector for the cross-correlation
cc=c_correlate(cutmh,sp_ref,lag)          ;cross-correlation (CC)
junk=max(cc,pos_lag)                      ;Simply picks the max of the cc
vcen=vcen-lag[pos_lag]                    ;shifts vcen

;;;;;;;;;;;;;;;;
;Quick loop to improve the centroid of each order as measured in "trace" (cut of the flat at mid-height)
lag=findgen(9)-4 ;x-vector for the cross-correlation
for k=0,nord-1 do begin
  b=cutmh[vcen[k]-2*ts:vcen[k]+2*ts]       ;isolates one order at a time
  a=findgen(n_elements(b))
  cc=c_correlate(b,ref,lag)                ;cross-correlation (CC)
  res=gaussfit(lag,cc-min(cc),co,nterms=3) ;gaussian fit to the CC profile
  vcen[k]=vcen[k]-co[1]                    ;correction of the centroid
endfor

;;;;;;;;;;;;;;;;
;Starting from the center, this loop does a CC at each row (first going down, then going up) to find the trace for each order
vy=findgen(len)  ;vector of the y-pixel position (from 0 to length-Y)
matcoeff=fltarr(nco+1,nord)
for k=0,nord-1 do begin
  ;borders of an order, centered on the measured centered at mid-height
  b1=fix(vcen[k]-2*ts)
  b2=b1+wdt
  
  templ=cutmh[b1:b2] ;template for the CC
  lag=findgen(11)-5  ;x-vector for the C-C (shifting over 11pix seems enough)
  cent=vcen[k]       ;just to remember the initial value of the center of the order at mid-height (for when we go on the other direction)
  vtr=fltarr(len)    ;vector of the y-pixels
  ;C-C going downward (we fit a gaussian to the CC profile to identify the "true" centroid at each y-pixel)
  for j=0,len/2-1 do begin
    prof=flat[fix(cent-wdt/2):fix(cent-wdt/2)+wdt,len/2-j] ;profil of the trace at current pixel
    if median(prof) gt 3 then begin ;if there is enough flux...
      
      lagt=lag           ;just a temporary variable, so we can remove some pixels when the CC does not behave perfectly
      cc=c_correlate(prof,templ,lag)
      pos=where(cc lt 0) ;points to all negative values in the CC (they will be removed, otherwise the gaussian fit can crash)
      if n_elements(pos) lt n_elements(cc) then begin
        if max(pos) ne -1 then remove,pos,lagt,cc ;where the negative values are removed
        cc=cc-min(cc)
        if n_elements(cc) gt 3 then begin
          res=gaussfit(lagt,cc,co,nterms=3) ;the gaussian fit of the CC profile, if it has at least 3 data points
          if abs(co[1]) gt 5 then cent=cent+lagt[where(cc eq max(cc))] else cent=cent-co[1] ;The centroid has to be within the range of "lag", else just takes the position of the max of the CC profil as the centroid
        endif else cent=cent+lagt[where(cc eq max(cc))]
        vtr[len/2-j]=cent ;records the centroid for that given line
      endif
    endif
  endfor
  cent=vcen[k] ;refreshes the centroid at the center of the flat image
  ;C-C going upward
  for j=1,len/2-1 do begin
    prof=flat[fix(cent-wdt/2):fix(cent-wdt/2)+wdt,len/2+j] ;profil
    if median(prof) gt 3 then begin
      ;Same as above
      lagt=lag
      cc=c_correlate(prof,templ,lag)
      pos=where(cc lt 0)
      if n_elements(pos) lt n_elements(cc) then begin
        if max(pos) ne -1 then remove,pos,lagt,cc
        cc=cc-min(cc)
        if n_elements(cc) gt 3 then begin
          res=gaussfit(lagt,cc,co,nterms=3)
          if abs(co[1]) gt 5 then cent=cent+lagt[where(cc eq max(cc))] else cent=cent-co[1]
        endif else cent=cent+lagt[where(cc eq max(cc))]
        vtr[len/2+j]=cent
      endif
    endif
  endfor
  
  ;We remove all the points for which the C-C could not be performed
  vyt=vy ;temporary y vector
  pos=where(vtr eq 0) ;failed measurements have a value of 0
  if max(pos) ne-1 then remove,pos,vyt,vtr
  
  ;Fit of the trace of the order (polynomial)
  res=poly_fit(vyt[1000:3500],vtr[1000:3500],nco) ;first selecting 2500 pix around the center
  ;Use the first fit to reject diviant points (mostly useful when an order is faint on the edges)
  xfit=findgen(len)
  yfit=0
  for i=0,nco do yfit=yfit+res[i]*xfit^i
  test=vtr-yfit
  sig=stddev(test)
  pos=where(abs(test-median(test)) gt 3*sig)
  while max(pos) ne -1 do begin
    remove,pos,test,vyt,vtr
    sig=stddev(test)
    pos=where(abs(test-median(test)) gt 3*sig)
  endwhile
  ;Final fit of the trace, including all the remaining points
  res=poly_fit(vyt,vtr,nco)  
  ;The fit of the trace is recorded in matcoeff
  matcoeff[*,k]=res
endfor

;mat coeff has the coefficients of the polynomial fits
;  along all the detected orders. They will be used to
;  get the aperture during extraction.
return,matcoeff
end



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;  FUNCTION EITHER USED TO FIND THE SLIT TILT OR THE LINE LIST (pix position)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; used to find the ThAr lines or to measure the slit tilt
function find_lines,sp_1d,sp_2d,spmd,wdt,resel,nonlin,find_tilt=find_tilt
s=size(sp_1d)
;the variable step is used to double the size of the matrix when also the sky
;  spectrum needs to be extracted in the 2-fiber mode
case spmd of
  '1f': step=1
  '2f': step=2
endcase

;These are first guesses of the wavelength calibration solution
;They were obtained from a ThAr spectrum from commissining, extracted
;   using DRAGRACES, and calibrated using IRAF. It is expected that
;   all ThAr would have comparable solutions, with a lateral shift
;   withing 0.1 A (i.e. few 100s km/s)
coeff_orders=[$
  ;[10336.0,0.147116,-5.51402e-05,2.90528e-08,-5.11140e-12],$ ;order 21
  [10071.5,0.103340,-3.04415e-06,-4.96906e-12,-1.24169e-15],$ ;order 22
  [9633.55,0.0988646,-2.89381e-06,-1.30645e-11,2.81250e-17],$ ;order 23
  [9232.09,0.0947520,-2.75468e-06,-2.04805e-11,9.89946e-16],$ ;order 24
  [8862.74,0.0909965,-2.65929e-06,-1.36005e-11,1.52931e-16],$ ;order 25
  [8521.81,0.0874968,-2.54023e-06,-1.89059e-11,8.05414e-16],$ ;order 26
  [8206.12,0.0842816,-2.45478e-06,-1.36680e-11,1.95089e-17],$ ;order 27
  [7912.99,0.0812925,-2.35935e-06,-1.90286e-11,9.70765e-16],$ ;order 28
  [7640.06,0.0785443,-2.31235e-06,-6.85041e-12,-3.56170e-16],$ ;order 29
  [7385.33,0.0759578,-2.24494e-06,-4.20288e-12,-5.77798e-16],$ ;order 30
  [7147.05,0.0734940,-2.14783e-06,-1.18985e-11,2.55312e-16],$ ;order 31
  [6923.65,0.0712137,-2.07533e-06,-1.44779e-11,6.29751e-16],$ ;order 32
  [6713.79,0.0690757,-2.01587e-06,-1.31607e-11,5.16460e-16],$ ;order 33
  [6516.27,0.0670719,-1.96697e-06,-9.50634e-12,1.13857e-16],$ ;order 34
  [6330.05,0.0651599,-1.90213e-06,-1.17891e-11,3.51085e-16],$ ;order 35
  [6154.16,0.0633734,-1.85263e-06,-1.13164e-11,3.64765e-16],$ ;order 36
  [5987.78,0.0616764,-1.80404e-06,-1.03363e-11,2.40366e-16],$ ;order 37
  [5830.16,0.0600725,-1.75770e-06,-1.02853e-11,2.77548e-16],$ ;order 38
  [5680.62,0.0585413,-1.70834e-06,-1.17815e-11,4.84551e-16],$ ;order 39
  [5538.55,0.0571003,-1.66974e-06,-1.05941e-11,3.72918e-16],$ ;order 40
  [5403.41,0.0557408,-1.64509e-06,-5.22269e-12,-2.10593e-16],$ ;order 41
  [5274.71,0.0544249,-1.60222e-06,-6.80346e-12,2.45565e-18],$ ;order 42
  [5151.99,0.0531873,-1.57105e-06,-5.57193e-12,-7.86697e-17],$ ;order 43
  [5034.86,0.0519739,-1.52342e-06,-8.94262e-12,2.61081e-16],$ ;order 44
  [4922.93,0.0508314,-1.48609e-06,-1.05462e-11,5.01844e-16],$ ;order 45
  [4815.86,0.0497462,-1.46099e-06,-7.50083e-12,1.08748e-16],$ ;order 46
  [4713.35,0.0486981,-1.42537e-06,-9.07843e-12,2.83862e-16],$ ;order 47
  [4615.10,0.0477166,-1.40866e-06,-5.60920e-12,-3.68590e-17],$ ;order 48
  [4520.87,0.0467529,-1.37851e-06,-5.73774e-12,-4.05722e-17],$ ;order 49
  [4430.41,0.0458190,-1.34127e-06,-8.75478e-12,2.76433e-16],$ ;order 50
  [4343.49,0.0449399,-1.31924e-06,-7.65794e-12,1.94136e-16],$ ;order 51
  [4259.91,0.0440855,-1.28792e-06,-1.00712e-11,5.02010e-16],$ ;order 52
  [4179.50,0.0432611,-1.26202e-06,-9.98627e-12,4.42814e-16],$ ;order 53
  [4102.05,0.0424826,-1.24638e-06,-7.64900e-12,2.01544e-16],$ ;order 54
  [4027.44,0.0416900,-1.20111e-06,-1.38099e-11,7.84916e-16],$ ;order 55
  [3955.45,0.0410027,-1.20714e-06,-6.74539e-12,1.71361e-16]$ ;order 56
  ;[] ;order 57
  ]

;ThAr line list (in Angstrom) from NOAO (iraf.noao.edu/specatlas/thar/thar_list)
readcol,'ThAr_list.lst',format='(f)',ThAr_list,/silent

lineHW=fix(resel*2)+1 ;value of a line half width at half max

;For each order, it guesses where the lines are and looks for the real
;  position using a convolution optimisation technique
px=findgen(s[2])
slit_tilt=fltarr(2,s[1])
line_list=-1
for i=0,s[1]-1 do begin
  sp=sp_1d[i,*]           ;1d spectrum of that given order
  
  ;Finds a rough continuum to the ThAr spectrum
  x_bin=findgen(s[2]/50)*50 ;x-vector for 50pix bins 
  y_bin=x_bin               ;y-vector for 50pix bins
  ;takes the median value for each bin
  for k=0,n_elements(x_bin)-1 do y_bin[k]=median(sp[k*50:(k+1)*50-1])
  ;removes deviant points (3sig)
  pos=0
  while max(pos) ne -1 do begin
    pos=where(abs(y_bin-median(y_bin)) gt 3*stddev(y_bin))
    if max(pos) ne -1 then remove,pos,x_bin,y_bin
  endwhile
  ;Estimates a continuum
  cont=spline(x_bin,y_bin,px)
  ;spN is the 1d spectrum "corrected" by the continuum
  spN=sp-cont
  
  ;Applies a mask that hides saturated lines and their ghosts
  case i/step of
    ;0:mask_range=[-1] ;order 21
    0:mask_range=[150,350,1550,1700] ;order 22
    1:mask_range=[150,350,1550,1700] ;order 23
    2:mask_range=[150,350] ;order 24
    3:mask_range=[3050,3300] ; order 25
    4:mask_range=[0,100,1700,1850] ; order 26
    5:mask_range=[650,800,2500,2950,4200,4300] ; order 27
    6:mask_range=[350,550,1150,1400,2450,2850] ; order 28
    7:mask_range=[1050,1150] ; order 29
    8:mask_range=[1550,1850,3600,3800] ; order 30
    9:mask_range=[1750,1850,3550,3700] ; order 31
    10:mask_range=[550,650,2100,2250] ; order 32
    11:mask_range=[4100,4300] ; order 33
    12:mask_range=[-1] ; order 34
    13:mask_range=[-1] ; order 35
    14:mask_range=[-1] ; order 36
    15:mask_range=[-1] ; order 37
    16:mask_range=[-1] ; order 38
    17:mask_range=[-1] ; order 39
    18:mask_range=[-1] ; order 40
    19:mask_range=[-1] ; order 41
    20:mask_range=[-1] ; order 42
    21:mask_range=[-1] ; order 43
    22:mask_range=[-1] ; order 44
    23:mask_range=[-1] ; order 45
    24:mask_range=[-1] ; order 46
    25:mask_range=[-1] ; order 47
    26:mask_range=[-1] ; order 48
    27:mask_range=[-1] ; order 49
    28:mask_range=[-1] ; order 50
    29:mask_range=[-1] ; order 51
    30:mask_range=[-1] ; order 52
    31:mask_range=[-1] ; order 53
    32:mask_range=[-1] ; order 54
    33:mask_range=[-1] ; order 55
    34:mask_range=[-1] ; order 56
    ;36:mask_range=[-1] ; order 57
  endcase
  if max(mask_range) ne -1 then begin
    mask=-1
    for j=0,n_elements(mask_range)/2-1 do if max(mask) eq -1 then mask=findgen(mask_range[j*2+1]-mask_range[j*2])+mask_range[j*2] else mask=[mask,findgen(mask_range[j*2+1]-mask_range[j*2])+mask_range[j*2],mask] 
    spN[mask]=0
  endif

  ;the guess solution comes from coeff_orders
  coeff=coeff_orders[*,i/step]
  vec_pix=findgen(n_elements(spN))
  wave_vs_pix=coeff[0]+coeff[1]*(vec_pix+1)+coeff[2]*(vec_pix+1)^2+coeff[3]*(vec_pix+1)^3+coeff[4]*(vec_pix+1)^4
  ThAr_list_order=ThAr_list[where(ThAr_list ge min(wave_vs_pix) and ThAr_list le max(wave_vs_pix))]
  pos_line_pix=spline(wave_vs_pix,vec_pix,ThAr_list_order)
  
  ;does a "mean" profile of all the lines
  ;  (convolving the ThAr spectrum with a "comb" at a moving
  ;  position around the guessed position in pixel of the lines)
  nlag=20
  optim_profile=fltarr(nlag*2+1)
  for j=0,2*nlag do begin
    pos=pos_line_pix+j-nlag
    pos=pos>0<4300
    spNint=spline(findgen(n_elements(spN)),spN,pos)
    optim_profile[j]=total(spNint)
  endfor
  
  ;finds the offset between the real line possitions and the guessed ones
  junk=max(optim_profile,pos)
  b1=pos-2*lineHW
  b1=b1>0
  b2=pos+2*lineHW
  b2=b2<n_elements(optim_profile)-1
  res=gaussfit(findgen(b2-b1+1)+b1,optim_profile[b1:b2],coeff,nterms=4)
  
  ;applies the shift
  pos_line_pix=pos_line_pix+(coeff[1]-nlag)
  pos=where(pos_line_pix lt 0 or pos_line_pix gt 4300)
  if max(pos) ne -1 then remove,pos,pos_line_pix,ThAr_list_order
  
  ;removes small lines
  amp_lines=spline(findgen(n_elements(spN)),spN,pos_line_pix)
  pos=where(amp_lines lt 50) ;50 is an semi-educated, but mostly arbitrary value
  if max(pos) ne -1 then remove,pos,pos_line_pix,ThAr_list_order
  
  ;does different things if you want to get the tilt of the line list
  if keyword_set(find_tilt) then begin
    vec_tilt=-1 ;vector for slit tilt values
    for j=0,n_elements(pos_line_pix)-1 do begin
      ;samples the 2d ThAr spectrum around the identified line
      yb=round(pos_line_pix[j]+lineHW*[-1,1])
      ;makes sure yb does not go outside the image borders
      yb=yb>0
      yb=yb<s[2]
      sample=sp_2d[i/step*wdt:(i/step+1)*wdt-1,yb[0]:yb[1]]
      ;uses the extracted shape of the line as reference (on non-continuum corrected spectrum)
      ref=spN[yb[0]:yb[1]]
      ;;would skip if a pixel-value is over the non-linearity threshold (in ADU)
      ;psat=where(sample[*,lineHw/2:lineHw*1.5] gt nonlin)
      
      ;cross-correlates the shape of the line accross "sample"
      vcen=fltarr(wdt)
      for k=0,wdt-1 do begin
        lag=findgen(9)-4
        b=sample[k,*]
        a=findgen(n_elements(b))
        cc=c_correlate(b,ref,lag)
        res=gaussfit(lag,cc-min(cc),co,nterms=3)
        vcen[k]=-co[1]
      endfor
      ;fits the slit tilt
      res=ladfit(findgen(n_elements(vcen)),vcen)
      if max(vec_tilt) eq -1 then vec_tilt=res[1] else vec_tilt=[vec_tilt,res[1]]
    endfor
    ;removes deviant measurement of tilt (3sig)
    pos=where(abs(vec_tilt-median(vec_tilt)) gt 3*stddev(vec_tilt))
    while max(pos) ne -1 do begin
      remove,pos,pos_line_pix,vec_tilt
      pos=where(abs(vec_tilt-median(vec_tilt)) gt 3*stddev(vec_tilt))
    endwhile
    ;fits the slit tilt as a function of pixel (gentle slope)
    res=ladfit(pos_line_pix,vec_tilt)
    slit_tilt[*,i]=res
  endif else begin
    ;improves the fit to individual lines
    fwhm_line_pix=fltarr(n_elements(pos_line_pix)) ;vector to record FWHM
    for j=0,n_elements(pos_line_pix)-1 do begin
      ;selects the 1D ThAr spectrum around the identified line
      yb=round(pos_line_pix[j]+lineHW*[-1,1])
      ;makes sure yb does not go outside the image borders
      yb=yb>0
      yb=yb<s[2]
      
      ;re-finds the line
      temp=spN[yb[0]:yb[1]]
      res=gaussfit(findgen(n_elements(temp)),temp,coeff,nterms=4)
      
      if coeff[1] ge 0 and coeff[1] le n_elements(temp) then begin
        fwhm_line_pix[j]=coeff[2]*2*sqrt(2*alog(2))
        pos_line_pix[j]=yb[0]+coeff[1]
      endif else pos_line_pix[j]=-2
    endfor
    for j=0,n_elements(pos_line_pix)-1 do if max(line_list) eq -1 then line_list=[i,pos_line_pix[j],ThAr_list_order[j],fwhm_line_pix[j]] else line_list=[[line_list],[i,pos_line_pix[j],ThAr_list_order[j],fwhm_line_pix[j]]] ; sets line_list
  endelse
endfor

if keyword_set(find_tilt) then return,slit_tilt else return,line_list
end



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;  FUNCTION WHERE THE SPECTRUM IS REDUCED
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function reduce,im,matcoeff,wdt,slit_tilt,mask=mask,normal=normal,notilt=notilt,disp=disp,dir=dir
s=size(im)
len=s[2] ;number of pixels along the spectrum
s=size(matcoeff)
nord=s[2] ;number of orders to extract
;final variable with the reduced 2D spectrum
sp2d=fltarr(nord*wdt,len)

;saves a display of the original ThAr frame with an overlay of the traces with the slit-tilt.
if keyword_set(disp) then begin
  sim=size(im)
  set_plot,'ps'
  device,filename=dir+'Trace'+disp+'.eps',xsize=sim[1]/300.,ysize=sim[2]/300.,/inches,/color,/encapsulated
  loadct,0
  imdisp,-1*im,/usepos,range=[-1000,0]
  loadct,13
endif

;variable to hold a strip corresponding to the rectified trace of the current order
band=fltarr(wdt,len)
for i=0,nord-1 do begin
  ;grabs the coefficient of the polynomial fit of the current order
  res=matcoeff[*,i]
  for k=0,len-1 do begin
    ;recreates the fit of the trace
    pos_center=0.
    for j=0,n_elements(res)-1 do pos_center=pos_center+res[j]*double(k)^j
    ;contains the tilt to add to the slit (=0 if no slit tilt is used)
    slope_extr=0.
    if keyword_set(notilt) ne 1 then slope_extr=slope_extr+(slit_tilt[0,i]+k*slit_tilt[1,i])
    
    ;interpolation along the trace
    b1=-4
    b2=4
    if k+b1 lt 0 then b1=-1*k
    if k+b2 gt len-1 then b2=len-1-k
    sample=im[pos_center-wdt/2:pos_center+wdt/2,k+b1:k+b2]
    temp=fltarr(wdt)
    b=(-1*b1)*slope_extr
    rge=dindgen(wdt+1)-wdt/2
    decim=pos_center-fix(pos_center)
    for l=0,wdt-1 do temp[l]=bilinear(sample,l+decim,l*slope_extr+b-b1) ; the interpolation is using the bilinear IDL function
    
    ;for the display
    if keyword_set(disp) then begin
      if i eq 0 and k eq 1 then plot,rge+pos_center,rge*slope_extr+b+k,xrange=[0,sim[1]-1],/xst,yrange=[0,sim[2]-1],/yst,/nodata,/noerase,position=[0,0,1,1]
      if k mod 10 eq 0 then oplot,rge+pos_center,rge*slope_extr+b+k,color=255
    endif
    
    ;if the pixels above the non-linearity threshold are masked
    if keyword_set(mask) then begin
      pos=where(temp gt mask)
      if max(pos) ne -1 then if n_elements(pos) eq n_elements(temp) then temp=fltarr(wdt) else temp[pos]=0
    endif
    band[*,k]=temp
  endfor
  ;If the option 'normal' is chosen, it normalizes
  ;  the 2d spectrum, line per line (only used for flatfield)
  if keyword_set(normal) then begin
    for j=0,wdt-1 do begin
      ytemp=band[j,*]
      ;Due to the vignetting of the 300 higher pix, we exclude them in the fit
      flat_fit_order=6 ;seems to work fine
      res=poly_fit(findgen(4301),ytemp[0:4300],flat_fit_order)
      yfit=res[0]
      for k=1,flat_fit_order do yfit=yfit+res[k]*findgen(len)^k
      band[j,*]=band[j,*]/yfit
    endfor
  endif
  ;sticks the "band" to the 2D spectrum, and loops to the next order
  sp2d[i*wdt:(i+1)*wdt-1,*]=band
endfor
;closes the figure displaying the traces
if keyword_set(disp) then begin
  device,/close
  set_plot,'x'
endif

return,sp2d
end



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;  FUNCTION TO CORRECT FOR ILLUMIATION
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function illumcorr,im,sp2d,matcoeff,wdt,slit_tilt
  s=size(matcoeff)
  nord=s[2] ;number of orders
  s=size(im)
  illum=fltarr(s[1],s[2]) ;Final 2D illumination (background)
  bckg_x=fltarr(nord+1,s[2]) ;pixels where the background level is measeured
  bckg_y=fltarr(nord+1,s[2]) ;background level at each order
  vp0=fltarr(s[2]) ;position at the start of the first order at each row
  
  for i=0,s[2]-1 do begin
    ;Uses the cut of the image at each row to determine the scattered light (illumination)
    illum_i=im[*,i]
    for k=0,nord-1 do begin
      if k eq 0 then begin
        ;points x and y where the background is measure for the current row
        sample_x=fltarr(nord+1)
        sample_y=fltarr(nord+1)
        
        res=matcoeff[*,0]
        pos_center=0. ;center of the first order (22)
        for j=0,n_elements(res)-1 do pos_center=pos_center+res[j]*double(i)^j
        
        p1l=round(pos_center-wdt/2-2) ;lowest point in x of the left side of the order
        vp0[i]=p1l
      endif else begin
        p1l=round(pos_center+wdt/2+.5) ;lowest point in x of the left side of the order
        
        res=matcoeff[*,k]
        pos_center=0. ;center of the current order
        for j=0,n_elements(res)-1 do pos_center=pos_center+res[j]*double(i)^j
      endelse
      if k eq nord-1 then begin
        p1r=round(pos_center+wdt/2+.5) ;the lowest point in x of the right side of the order
        p2r=n_elements(illum_i)-1 ;the highest point in x of the right side of the order
      endif
      
      p2l=fix(pos_center-wdt/2) ;the highest point in x of the left side of the order
      ;flips the points if they are inverted (happens if the order are really close)
      if p1l gt p2l then begin
        junk=p1l
        p1l=p2l
        p2l=junk
      endif

      if p2l-p1l gt 3 then begin
        res=poly_fit(findgen(p2l-p1l+1),illum_i[p1l:p2l],3,yfit=yfit)
        sample_y[k]=min(yfit) ;illumination value
      endif else sample_y[k]=min(illum_i[p1l:p2l]) ;illumination value
      sample_x[k]=(p1l+p2l)/2. ;x position of the illumination value
      if k eq nord-1 then begin
        sample_y[k+1]=median(illum_i[p1r:p2r]) ;illumination value
        sample_x[k+1]=(p1r+p2r)/2. ;x position of the illumination value
        bckg_x[*,i]=sample_x
        bckg_y[*,i]=sample_y
      endif
    endfor
  endfor
  ;This loop goes over the image again, but averages the background measured for nbin (40) rows before creating the background frame
  for i=0,s[2]-1 do begin
    illum_i=im[*,i]
    sample_x=fltarr(nord+1)
    sample_y=fltarr(nord+1)
    for k=0,nord do begin
      ;smoothens in y by taking min of all values within 2xnbin pixels
      nbin=40
      smth_y_low=i-nbin/2
      smth_y_high=i+nbin/2
      if smth_y_low lt 0 then smth_y_low=0
      if smth_y_high ge s[2] then smth_y_high=s[2]-1
      
      sample_y[k]=median(bckg_y[k,smth_y_low:smth_y_high])
      sample_x[k]=bckg_x[k,i]
    endfor
    ;fits the measures for a smoother, and resampled background
    xint=findgen(s[1]-vp0[i])+vp0[i]
    res=poly_fit(sample_x,sample_y,6)
    yint=0
    for j=0,6 do yint=yint+res[j]*xint^j
    illum_i[vp0[i]:s[1]-1]=yint
    illum[*,i]=illum_i
  endfor
  
  ;Creates the 2D background (illumination) frame
  illum2d=reduce(illum,matcoeff,wdt,slit_tilt)
  ;Fits a 2D surface to each order
  for i=0,nord-1 do illum2d[i*wdt:i*wdt+wdt-1,*]=sfit(illum2d[i*wdt:i*wdt+wdt-1,*],5)
  
  ;Returns the corrected 2D frame
  return,sp2d-illum2d
end



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;  FUNCTION WHERE THE SPECTRUM IS EXTRACTED
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function extract,sp2d,wdt,spmd
s=size(sp2d)
nord=s[1]/wdt ;number of orders extracted
len=s[2]      ;number of pixels along the spectrum

;this is simply adding up the columns corresponding to the trace. Nothing sophisticated there...
case spmd of
  '1f': begin
    sp1d=fltarr(nord,s[2])
    for i=0,nord-1 do begin
      ;Add rejection method HERE
      sp1d[i,*]=(total(sp2d[i*wdt:(i+1)*wdt-1,*],1))
    endfor
  end
  '2f': begin
    sp1d=fltarr(nord*2,s[2])
    for i=0,nord-1 do begin
      ;Add rejection method HERE
      sp1d[i*2,*]=(total(sp2d[i*wdt:((i+1)-.5)*wdt-1,*],1))
      sp1d[i*2+1,*]=(total(sp2d[(i+.5)*wdt:(i+1)*wdt-1,*],1))
    endfor
  end
endcase

;Kills the last 300pix (affected by vignetting)
sp1d[*,4300:len-1]=0
return,sp1d

end



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;  FUNCTION WHERE THE WAVELENGTH SOLUTION IS CALCULATED
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; used to find the ThAr lines or to measure the slit tilt
function wavel_sol,line_list,spmd,reddir,vec_lines_used=vec_lines_used

;depends on which spectral mode we are in
case spmd of
  '1f': step=1
  '2f': step=2
endcase

nord=(max(line_list[0,*])+1) ;number of orders extracted
res=fltarr(5,nord)           ;variable with the final wavelength solution for all the orders
vec_lines_used=[-1,-1]       ;vector of the lines that were used for the final solution
for i=0,nord-1 do begin
  pos=where(line_list[0,*] eq i)
  pos_line_pix=line_list[1,pos]
  line_wave=line_list[2,pos]
  pos=where(pos_line_pix lt 0)
  if max(pos) ne -1 then remove,pos,pos_line_pix,line_wave
    
  ;Does the final fit of the pixel position vs wavelength for all the remaining lines
  ;the while loop is to remove deviant points
  coeff=robust_poly_fit(pos_line_pix,line_wave,4,yfit)
  res[*,i]=coeff 
  
  for k=0,n_elements(pos_line_pix)-1 do vec_lines_used=[[vec_lines_used],[pos_line_pix[k],line_wave[k]]]
  vec_lines_used=[[vec_lines_used],[-1,-1]]
endfor

;Creates plots of the resolution power (R) as a function of wavelength
set_plot,'ps'
device,filename=reddir+'Resolution'+spmd+'.eps',xsize=6,ysize=4,/inches,/color,/encapsulated
loadct,13
plotsym,0,/fill
dispersion=(res[1,line_list[0,*]]+res[2,line_list[0,*]]*4608+res[3,line_list[0,*]]*4608^2+res[4,line_list[0,*]]*4608^3)
resolution=(line_list[3,*]*dispersion/line_list[2,*])^(-1)/1000
plot,line_list[2,*],resolution,psym=8,symsize=.5,yrange=median(resolution)+20*[-1,1],/xst,/yst,xtitle='wavelength (A)',ytitle='resolution power (x1000)'
if max(where(finite(resolution)) ne 1) ne -1 then remove,where(finite(resolution) ne 1),resolution
oplot,[3500,1.1e4],median(resolution)*[1,1],color=255
oplot,[3500,1.1e4],median(resolution)*[1,1]-stddev(resolution),linestyle=1,color=255
oplot,[3500,1.1e4],median(resolution)*[1,1]+stddev(resolution),linestyle=1,color=255
device,/close
set_plot,'x'


return,res
end



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;  FUNCTION THE WAVELENGTH CALIBRATION IS APPLIED
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function wavel_cal,sp,wavel_sol,param=param
s=size(sp)
final=fltarr(s[1],s[2]) ;final calibrated spectrum
param=fltarr(2,s[1])    ;contains CRVAL1 and CD1_1
vec_pix=findgen(s[2])   ;vector of pixel positions from 1 to #pix
;resamples the spectrum to linearize the wavelength solution
for i=0,s[1]-1 do begin
  coeff=wavel_sol[*,i]
  x_lamb=coeff[0]+coeff[1]*vec_pix+coeff[2]*vec_pix^2+coeff[3]*vec_pix^3+coeff[4]*vec_pix^4 ;wavelength solution
  param[0,i]=min(x_lamb)
  param[1,i]=(max(x_lamb)-min(x_lamb))/s[2]
  x_interp=findgen(s[2])*param[1,i]+param[0,i]
  final[i,*]=spline(x_lamb,sp[i,*],x_interp) ;resampling
endfor


return,final
end


pro dg,dir=dir,utdate=utdate,lbias=lbias,lflat=lflat,lthar=lthar,skip_wavel=skip_wavel,ascii_file=ascii_file,new=new,logo=logo,help=help

  print,''
  print,''
  print,'________ __________    _____    __________________    _____  _________ ___________ _________'
  print,'\______ \\______   \  /  _  \  /  _____/\______   \  /  _  \ \_   ___ \\_   _____//   _____/'
  print,' |    |  \|       _/ /  /_\  \/   \  ___ |       _/ /  /_\  \/    \  \/ |    __)_ \_____  \ '
  print,' |    `   \    |   \/    |    \    \_\  \|    |   \/    |    \     \____|        \/        \'
  print,'/_______  /____|_  /\____|__  /\______  /|____|_  /\____|__  /\______  /_______  /_______  /'
  print,'        \/       \/         \/        \/        \/         \/        \/        \/        \/ '
  print,' Data Reduction and Analysis for GRACES (Gemini Remote Access to CFHT ESPaDOnS Spectrograph)'
  if keyword_set(logo) then begin
    print,'   ______  ______  ______  ______  ______  ______  ______  ______  ______  ______  ______       ,****,,,.........                                           '
    print,'  /_____/ /_____/ /_____/ /_____/ /_____/ /_____/ /_____/ /_____/ /_____/ /_____/ /_____/     .*                 .*                                         '
    print,'                                                                                             *.  ,.........,,,,..    .,*                                    '
    print,'                                                                                           ,*  .,........      ..       *                                   '
    print,'                                                                                     .,  *,   ,,,..........    ,    ,/.,*.                                  '
    print,'                                                                                   ,    *.   .,,,,,.......... ,     ,.,      .,                             '
    print,'                     .*********.                                                  ,    * .   *,,,,,,,..........     . *           ,                         '
    print,'                 *,,*.           .**                                             ,,**,.,*,  ***,,,,,,,,,.....,     . ,/,.          .                        '
    print,'              **.*,                 .*.                                         .,  */,**,  *****,,,,,,,,,..,.     ,/*.    .,,**.   ,                       '
    print,'            **,,*                      *.                                        .**,,.,*  ********,,,,,,,,,*       *,*   ****.  .,*..                      '
    print,'           /*,,*                        .*                                        .** ***  /*********,,,,,,,,       *,,   *,  *,,,,  *                      '
    print,'         ./***,                           *                                       *** **,. ...,*********,,,*       .**   .,* .,* .,*,                       '
    print,'         /****    *     *      *           *                                      **,,**.                    *.    ,**   ,** *** *,*                        '
    print,'        */**/    *     ,.     *     .,     ,.                                    .*/ **,.                    *.    **,   **, **. **.                        '
    print,'        //**,   .,     *      *     *       *                                    */*.// .                    *     //.   **..** ,**                         '
    print,'       .//*/    *.    .*,    ,.    *        *                                   ,**,*// ,                    *    ,//   ,// */* /*,                         '
    print,'       .//*,                       *        *                                 ,..../  , .                    ,   ****.  */* //,.//                          '
    print,'        *                                  .,                    /&@/%#.      .    *  . .,*****,,...        .,   , *****,,.    .,*                          '
    print,'        * ....,,,,.....   ,,               *                  *..,(./&&((.    ,    ,..,         .....,,,,****,.  , *          ..,,,*   ,%@&/#&,             '
    print,'        .,,,,,,,,,,,,,,,...?, ,..,,,.     ,.                ...../,#,,,,%%/,,,,**,..*/.                         .,.,               ,.,..,/,(@@#(,           '
    print,'         *...,,,,,,....      ,, ..,,.   .*,                 ,...,/**...*.                     ...,,*****,,...   .,*               ,...../,#,,..(&/          '
    print,'       ,.,,,,....,,,,,,,,..,    ..,,,. .*.                  ,,,,&//,*                                                ,,  ...,******,...,/**......@#         '
    print,'        ..                   ,, ..,,.   .,                  @@@@@(***                                                             ,,,,*@//,       @%        '
    print,'       ,((/.                        ,#''+#,                  &@@@@%(%*                          *,*,*.*                            .@@@@@(//.      .@.       '
    print,'  ''@@*#@\#%&%                  ,`@@./%##&#                 .@@@@@#@*,...                   ,*****,*.*                             %@@@@&(%.. .    %&       '
    print,'  @#@@/(**&&*,%%               ,:@@(/#,&&%#%%                ,@@@@@%@&.....,&@.....,*/////***/**////*,...,,                         @@@@@(@&...   .#@       '
    print,' (*@@@&&.+?`.&,%%              /@@@@&%,`.+#%%,,*********,,,..,@@@@&&@@#,,*/##*,.                            ?,,.....,,,,****//*****/@@@@@%@@......&@.       '
    print,'  ,@@@%&++...*%.&+............. @@@@&&/+...#(/&+                .####%%@@@@@@@(                                                       .@@@@@%@@&,,,%.       '
    print,'   (@@@%&+,(&&%.&%              ,@@@@%&%.#&&(/&%                ,@@@#&@@@@(.                                                         (@@@@#@@@@@@@/         '
    print,'    .(@@/%&/#&&%                 ,/@@@&&%#&&&*                       .,*/*,.                                                             *&@@@@@@*          '
    print,'      .****(%&&.                  .*///#%&&%*                                                                                                .,,,.          '
    print,''
  endif
  print,''
  print,'~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*'
  print,''
  print,'Author: Andre-Nicolas Chene'
  print,''
  print,'version 1.4'
  print,''
  print,'Please cite the AJ paper:'
  print,'https://ui.adsabs.harvard.edu/abs/2021AJ....161..109C/abstract'
  print,''

  if keyword_set(help) then begin
    print,''
    print,'NAME:'
    print,'    dg.pro'
    print,'EXPLANATION:'
    print,'    This program is designed to reduce and extract spectra obtained with the Gemini/CFHT high-'
    print,'    resolution spectrograph GRACES (Gemini Remote Access to the CFHT ESPaDOnS Spectrograph).'
    print,'    (www.gemini.edu/sciops/instruments/graces)'
    print,'    '
    print,'CALLING SEQUENCE:'
    print,'    dg[,dir=dir,utdate=utdate,lbias=lbias,lflat=lflat,lthar=lthar,/skip_wavel,/ascii_file,/new,/help]'
    print,'    '
    print,'INPUTS (all optional):'
    print,'        dir - Path to where the data can be found. A new directory named ""Reduction"" will'
    print,'              be created, if it does not already exist.'
    print,'     utdate - Date when the spectra where observed. The format is YYYYMMDD, and can be'
    print,'              given as a number (without ''''). If it not provided, it is expected that the'
    print,'              data in the directory pointed by the input ""dir"" are all from the same date. '
    print,'      lbias - Name of a file where the bias frames are listed.'
    print,'      lflat - Name of a file where the flats frames are listed.'
    print,'      lthar - Name of a file where the ThAr frames are listed.'
    print,' skip_wavel - When provided, the wavelength solution is not calculated.'
    print,' ascii_file - (NOT YET WORKING)'
    print,'        new - Forces DRAGRACES to recalculate all the calibrations instead of using those'
    print,'              obtained in a previous run of the pipeline'
    print,'       logo - Displays the DRAGRACES logo when you run the pipeline!'
    print,'       help - Displays this help summary.'
    print,'OUTPUTS:'
    print,'    The program saves the extracted spectra in the fits format, in the directory ""Reduction"".'
    print,'    If the wavelength solution was calculated, the filenames start with ""ext_"". If not, '
    print,'    the filenames start with ""red_"", and an additional ""red_Thar_..." spectrum is provided.'
    print,''
    print,'PREREQUISITES:'
    print,'    - IDL (only tested under IDL 8+)'
    print,'    - Astrolib (http://idlastro.gsfc.nasa.gov/ftp/)'
    print,'    - imdisp.pro (https://www.idlcoyote.com/programs/coyoteplus/imdisp.pro)'
    print,'         . if you do no want imdisp, simply comment line 955'
    print,'    note1: make sure that Astrolib and imdisp are copied into your IDL libraries'
    print,''
    print,'IMPORTANT NOTE:'
    print,''
    print,'    The pipeline expects the files to be ""well behaved"". It expects all the files in a given'
    print,'    folder to be relevant to the extraction. If calibrations data of different nigts are in the'
    print,'    same folder, they will be blended and potentially wrong. Unless you use the file list'
    print,'    option, but they have not been fully tested yet. Visit www.gemini.edu/node/12552 for more'
    print,'    information.'
    print,''
    print,'ENJOY!'
    print,''
    goto,fin
  endif



  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;  INPUTS
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ;;;;;;;;;;;;;;;;
  ;Directory were the data are:
  if keyword_set(dir) then begin
    datadir=dir
    if strcmp(strmid(dir,strlen(dir)-1),'/') ne 1 then dir=dir+'/'
  endif else datadir='./'
  if strcmp(strmid(datadir,strlen(datadir)-1),'/') ne 1 then datadir=datadir+'/' ;adds the / if not included in the path
  reddir=datadir+'Reduction/'
  ;creates the reduction directory if it does not exist
  if file_test(reddir) ne 1 then spawn,'mkdir '+reddir
  ;converts the entered utdate into string
  if keyword_set(utdate) then utdate=strtrim(string(utdate),2)

  ;;;;;;;;;;;;;;;;
  ;Parameters for detector
  nonlin=46400 ;non-linearity threshold in ADU
  satura=65000 ;saturation limit in ADU


  ;;;;;;;;;;;;;;;;
  ;Parameters for the traces (can be tunned)
  ;;for 1-fiber more:
  ;;width of a single section of the sliced image (e.g.          |-----|            )
  ;;                                                        ----- ----- ----- ----
  ts1f=7

  ;;for 2-fiber more:
  ;;width of a single section of the sliced image (e.g.          |-----|              )
  ;;                                                        ----- -----    ----- ----
  ts2f=8

  ;;width of the trace
  wdt1f=4*ts1f+1 ;1-fiber
  wdt2f=4*ts2f+1 ;2-fiber
  ;;order of the polynomial fitting the trace
  nco=4
  ;;first order to extract
  first_order=22
  ;hardcoded x-pixel values used as a first guess where to find the traces
  ;vcen=[30, 60,  91,  121,  152,  183,  215,  246,  279,  312,  346,  381, 416, 453, 491, 529, 569, 610, 652, 696, 740, 786, 834, 883, 933, 985, 1039, 1094,  1151,  1210,  1270,  1333,  1397,  1464,  1533, 1604, 1678] ;from order 21 to 57
  vcen=[60,  91,  121,  152,  183,  215,  246,  279,  312,  346,  381, 416, 453, 491, 529, 569, 610, 652, 696, 740, 786, 834, 883, 933, 985, 1039, 1094,  1151,  1210,  1270,  1333,  1397,  1464,  1533, 1604] ;from order 22 to 56
  ;number of orders to extract
  nord=n_elements(vcen)

  ;;;;;;;;;;;;;;;;
  ;Parameters for the line identifycation (can be tunned)
  ;optical resolution element in pixels
  resel1f=1.74
  resel2f=2.88



  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;  LISTS CREATION
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  if keyword_set(utdate) then begin
    ;Finds all the GRACES frames observed that utdate
    lstot=findfile(datadir+'N'+utdate+'G*.fits')
  endif else begin
    ;Finds all the GRACES frames in the data directory
    lstot=findfile(datadir+'N????????G*.fits')
  endelse
  ;Now the date is known, tweaks vcen based on date to account for shifts due to the fiber being taken out and back into the spectrograph between Sept and Dec 2018
  if float(strmid(lstot[0],strpos(lstot[0],'.fits')-13,8)) gt 20181101 then vcen=vcen+15
  ;checks if it any frame exist
  if min(strcmp(lstot,'')) eq 1 then begin
    print,''
    print,'ERROR: Could not find any GRACES observations.'
    goto,fin
  endif
  ;creates the calibrations and the science lists
  lbiasNr=[]   ;Biases in normal read mode
  lbiasSr=[]   ;Biases in slow read mode
  lflat1f=[]   ;Flats for 1-fiber mode
  lflat2f=[]   ;Flats for 2-fiber mode
  lthar1f=[]   ;ThAr for 1-fiber mode
  lthar2f=[]   ;ThAr for 2-fiber mode
  lsci1f=[]    ;Objects for 1-fiber mode
  lsci2f=[]    ;Objects for 2-fiber mode
  for i=0,n_elements(lstot)-1 do begin
    h=headfits(lstot[i])
    if strcmp(strtrim(sxpar(h,'OBSTYPE'),2),'BIAS') then begin
      if sxpar(h,'RDNOISEA') eq 4.2 then begin
        if lbiasNr eq !NULL then lbiasNr=lstot[i] else lbiasNr=[lbiasNr,lstot[i]]
      endif else if lbiasSr eq !NULL then lbiasSr=lstot[i] else lbiasSr=[lbiasSr,lstot[i]]
    endif else begin
      if strcmp(sxpar(h,'GSLIPOS'),'FOURSLICE') then begin
        case strtrim(sxpar(h,'OBSTYPE'),2) of
          'FLAT': if lflat1f eq !NULL then lflat1f=lstot[i] else lflat1f=[lflat1f,lstot[i]]
          'ARC': if lthar1f eq !NULL then lthar1f=lstot[i] else lthar1f=[lthar1f,lstot[i]]
          'OBJECT': if lsci1f eq !NULL then lsci1f=lstot[i] else lsci1f=[lsci1f,lstot[i]]
        endcase
      endif else begin
        case strtrim(sxpar(h,'OBSTYPE'),2) of
          'FLAT': if lflat2f eq !NULL then lflat2f=lstot[i] else lflat2f=[lflat2f,lstot[i]]
          'ARC': if lthar2f eq !NULL then lthar2f=lstot[i] else lthar2f=[lthar2f,lstot[i]]
          'OBJECT': if lsci2f eq !NULL then lsci2f=lstot[i] else lsci2f=[lsci2f,lstot[i]]
        endcase
      endelse
    endelse
  endfor



  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;  BIAS FRAMES CREATION
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ;Takes the user's list of bias, if provided
  if keyword_set(lbias) then begin
    readcol,lbias,format='(a)',lbias,/silent
    h=headfits(lbias[0])
    if sxpar(h,'RDNOISEA') eq 4.2 then lbiasNr=lbias else lbiasSr=lbias
  endif

  ;checks if there was any bias observed in Normal readout mode.
  if lbiasNr ne !NULL then begin
    ;gets the date from the header of the first bias
    hb=headfits(lbiasNr[0])
    date=sxpar(hb,'DATE')
    biasdate=strmid(date,0,4)+strmid(date,5,2)+strmid(date,8,2)

    ;checks first if the Master has not been created before
    if file_test(reddir+'BiasNormal'+biasdate+'.fits') and keyword_set(new) ne 1 then begin
      ;if a Master Bias has already been processed, it takes this one.
      biasNr=readfits(reddir+'BiasNormal'+biasdate+'.fits',/silent)
    endif else begin
      ;IMPORTANT NOTE: MinMax rejection method used!
      for i=0,n_elements(lbiasNr)-1 do begin
        im=readfits(lbiasNr[i],/silent)
        im_overs=overs_corr(im,hb) ;overscan corrected frame
        if i eq 0 then begin
          ;if it is the first frame, it creates the matrix that will hold all the bias frames
          dim_bias_frame=size(im_overs)
          mat_biases=fltarr(dim_bias_frame[1],dim_bias_frame[2],n_elements(lbiasNr))
        endif
        mat_biases[*,*,i]=im_overs ;adds the bias frame to mat_biases
      endfor
      biasNr=fltarr(dim_bias_frame[1],dim_bias_frame[2]) ;finale master bias frame
      for i=0,dim_bias_frame[1]-1 do begin
        for j=0,dim_bias_frame[2]-1 do begin
          temp=mat_biases[i,j,*]
          ;removes minimum value
          junk=min(temp,pos)
          remove,min(pos),temp
          ;removes maximum value
          junk=max(temp,pos)
          remove,max(pos),temp
          ;records the mean
          biasNr[i,j]=mean(temp)
        endfor
      endfor
      writefits,reddir+'BiasNormal'+biasdate+'.fits',biasNr;,hb
    endelse
  endif
  ;checks if there was any bias observed in Slow readout mode.
  if lbiasSr ne !NULL then begin
    ;gets the date from the header of the first bias
    hb=headfits(lbiasSr[0])
    date=sxpar(hb,'DATE')
    biasdate=strmid(date,0,4)+strmid(date,5,2)+strmid(date,8,2)

    ;checks first if the Master has not been created before
    if file_test(reddir+'BiasSlow'+biasdate+'.fits') and keyword_set(new) ne 1 then begin
      ;if a Master Bias has already been processed, it takes this one.
      biasSr=readfits(reddir+'BiasSlow'+biasdate+'.fits',/silent)
    endif else begin
      ;IMPORTANT NOTE: MinMax rejection method used!
      for i=0,n_elements(lbiasSr)-1 do begin
        im=readfits(lbiasSr[i],/silent)
        im_overs=overs_corr(im,hb) ;overscan corrected frame
        if i eq 0 then begin
          ;if it is the first frame, it creates the matrix that will hold all the bias frames
          dim_bias_frame=size(im_overs)
          mat_biases=fltarr(dim_bias_frame[1],dim_bias_frame[2],n_elements(lbiasSr))
        endif
        mat_biases[*,*,i]=im_overs ;adds the bias frame to mat_biases
      endfor
      biasSr=fltarr(dim_bias_frame[1],dim_bias_frame[2]) ;finale master bias frame
      for i=0,dim_bias_frame[1]-1 do begin
        for j=0,dim_bias_frame[2]-1 do begin
          temp=mat_biases[i,j,*]
          ;removes minimum value
          junk=min(temp,pos)
          remove,min(pos),temp
          ;removes maximum value
          junk=max(temp,pos)
          remove,max(pos),temp
          ;records the mean
          biasSr[i,j]=mean(temp)
        endfor
      endfor
      writefits,reddir+'BiasSlow'+biasdate+'.fits',biasSr;,hb
    endelse
  endif
  ;cuts the program short if no bias is available
  if lbiasNr eq !NULL and lbiasSr eq !NULL then begin
    print,''
    print,'ERROR: Could not find any Bias frame'
    goto,fin
  endif



  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;  ThAr FRAMES CREATION
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ;Takes the user's list of ThAr, if provided
  if keyword_set(lthar) then begin
    readcol,lthar,format='(a)',lthar,/silent
    h=headfits(lthar[0])
    if strcmp(sxpar(h,'GSLIPOS'),'FOURSLICE') then lthar1f=lthar else lthar2f=lthar
  endif

  ;checks if there was any ThAr observed in 1-fiber mode.
  if lthar1f ne !NULL then begin
    ;gets the date from the header of the first ThAr
    ht=headfits(lthar1f[0])
    date=sxpar(ht,'DATE')
    thar1fdate=strmid(date,0,4)+strmid(date,5,2)+strmid(date,8,2)

    ;checks first if the Master has not been created before
    if file_test(reddir+'ThAr1f'+thar1fdate+'.fits') and keyword_set(new) ne 1 then begin
      ;if a Master ThAr has already been processed, it takes this one.
      thar1f=readfits(reddir+'ThAr1f'+thar1fdate+'.fits',/silent)
    endif else begin
      ;IMPORTANT NOTE: No rejection method used for ThAr spectra!
      for i=0,n_elements(lthar1f)-1 do begin
        im=readfits(lthar1f[i],/silent)
        im_overs=overs_corr(im,ht)
        if i eq 0 then thar1f=im_overs/n_elements(lthar1f) else thar1f=thar1f+im_overs/n_elements(lthar1f)
      endfor
      thar1f=thar1f-biasNr
      writefits,reddir+'ThAr1f'+thar1fdate+'.fits',thar1f;,ht
    endelse
  endif
  ;checks if there was any ThAr observed in 2-fiber mode.
  if lthar2f ne !NULL then begin
    ;gets the date from the header of the first ThAr
    ht=headfits(lthar2f[0])
    date=sxpar(ht,'DATE')
    thar2fdate=strmid(date,0,4)+strmid(date,5,2)+strmid(date,8,2)

    ;checks first if the Master has not been created before
    if file_test(reddir+'ThAr2f'+thar2fdate+'.fits') and keyword_set(new) ne 1 then begin
      ;if a Master ThAr has already been processed, it takes this one.
      thar2f=readfits(reddir+'ThAr2f'+thar2fdate+'.fits',/silent)
    endif else begin
      ;IMPORTANT NOTE: No rejection method used for ThAr spectra!
      for i=0,n_elements(lthar2f)-1 do begin
        im=readfits(lthar2f[i],/silent)
        im_overs=overs_corr(im,ht)
        if i eq 0 then thar2f=im_overs/n_elements(lthar2f) else thar2f=thar2f+im_overs/n_elements(lthar2f)
      endfor
      thar2f=thar2f-biasNr
      writefits,reddir+'ThAr2f'+thar2fdate+'.fits',thar2f;,ht
    endelse
  endif
  ;cuts the program short if no bias is available
  if lthar1f eq !NULL and lthar2f eq !NULL then begin
    print,''
    print,'ERROR: Could not find any ThAr spectrum'
    goto,fin
  endif



  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;  FLAT FRAMES CREATION
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ;Takes the user's list of flats, if provided
  if keyword_set(lflat) then begin
    readcol,lflat,format='(a)',lflat,/silent
    h=headfits(lflat[0])
    if strcmp(sxpar(h,'GSLIPOS'),'FOURSLICE') then lflat1f=lflat else lflat2f=lflat
  endif

  ;checks if there was any Flat observed in 1-fiber mode.
  if lflat1f ne !NULL then begin
    ;gets the date from the header of the first flat
    hf=headfits(lflat1f[0])
    date=sxpar(hf,'DATE')
    flat1fdate=strmid(date,0,4)+strmid(date,5,2)+strmid(date,8,2)

    ;checks first if the Master has not been created before
    if file_test(reddir+'Flat1f'+flat1fdate+'.fits') and keyword_set(new) ne 1 then begin
      ;if a Master Flat has already been processed, it takes this one.
      flat1f=readfits(reddir+'Flat1f'+flat1fdate+'.fits',/silent)
    endif else begin
      ;IMPORTANT NOTE: MinMax rejection method used!
      for i=0,n_elements(lflat1f)-1 do begin
        im=readfits(lflat1f[i],/silent)
        im_overs=overs_corr(im,hf) ;overscan corrected frame
        if i eq 0 then begin
          ;if it is the first frame, it creates the matrix that will hold all the flat frames
          dim_flat_frame=size(im_overs)
          mat_flats=fltarr(dim_flat_frame[1],dim_flat_frame[2],n_elements(lflat1f))
        endif
        mat_flats[*,*,i]=im_overs ;adds the flat frame to mat_flats
      endfor
      flat1f=fltarr(dim_flat_frame[1],dim_flat_frame[2]) ;finale master bias frame
      for i=0,dim_flat_frame[1]-1 do begin
        for j=0,dim_flat_frame[2]-1 do begin
          temp=mat_flats[i,j,*]
          ;removes minimum value
          junk=min(temp,pos)
          remove,min(pos),temp
          ;removes maximum value
          junk=max(temp,pos)
          remove,max(pos),temp
          ;records the mean
          flat1f[i,j]=mean(temp)
        endfor
      endfor
      flat1f=flat1f-biasNr
      writefits,reddir+'Flat1f'+flat1fdate+'.fits',flat1f;,hf
    endelse
  endif
  ;checks if there was any Flat observed in 2-fiber mode.
  if lflat2f ne !NULL then begin
    ;gets the date from the header of the first flat
    hf=headfits(lflat2f[0])
    date=sxpar(ht,'DATE')
    flat2fdate=strmid(date,0,4)+strmid(date,5,2)+strmid(date,8,2)

    ;checks first if the Master has not been created before
    if file_test(reddir+'Flat2f'+flat2fdate+'.fits') and keyword_set(new) ne 1 then begin
      ;if a Master Flat has already been processed, it takes this one.
      flat2f=readfits(reddir+'Flat2f'+flat2fdate+'.fits',/silent)
    endif else begin
      ;IMPORTANT NOTE: MinMax rejection method used!
      for i=0,n_elements(lflat2f)-1 do begin
        im=readfits(lflat2f[i],/silent)
        im_overs=overs_corr(im,hf) ;overscan corrected frame
        if i eq 0 then begin
          ;if it is the first frame, it creates the matrix that will hold all the flat frames
          dim_flat_frame=size(im_overs)
          mat_flats=fltarr(dim_flat_frame[1],dim_flat_frame[2],n_elements(lflat2f))
        endif
        mat_flats[*,*,i]=im_overs ;adds the flat frame to mat_flats
      endfor
      flat2f=fltarr(dim_flat_frame[1],dim_flat_frame[2]) ;finale master bias frame
      for i=0,dim_flat_frame[1]-1 do begin
        for j=0,dim_flat_frame[2]-1 do begin
          temp=mat_flats[i,j,*]
          ;removes minimum value
          junk=min(temp,pos)
          remove,min(pos),temp
          ;removes maximum value
          junk=max(temp,pos)
          remove,max(pos),temp
          ;records the mean
          flat2f[i,j]=mean(temp)
        endfor
      endfor
      flat2f=flat2f-biasNr
      writefits,reddir+'Flat2f'+flat2fdate+'.fits',flat2f;,hf
    endelse
  endif
  ;cuts the program short if no bias is available
  if lflat1f eq !NULL and lflat2f eq !NULL then begin
    print,''
    print,'ERROR: Could not find any Flat spectrum'
    goto,fin
  endif



  ;Loops on the 2 possible spectral modes
  for idx=1,2 do begin
    case idx of
      ;sets-up parameters for the 1-fiber mode
      1: begin
        if lflat1f ne !NULL and lthar1f ne !NULL then begin
          spmd='1f'
          lflat=lflat1f
          lthar=lthar1f
          lsci=lsci1f
          flat=flat1f
          thar=thar1f
          flatdate=flat1fdate
          thardate=thar1fdate
          ts=ts1f
          wdt=wdt1f
          resel=resel1f
        endif else begin
          print,''
          if lsci1f ne !NULL then print,'ERROR: missing calibrations for the 1-fiber mode'
          goto,skip_loop
        endelse
      end
      ;sets-up parameters for the 2-fiber mode
      2: begin
        if lflat2f ne !NULL and lthar2f ne !NULL then begin
          spmd='2f'
          lflat=lflat2f
          lthar=lthar2f
          lsci=lsci2f
          flat=flat2f
          thar=thar2f
          flatdate=flat2fdate
          thardate=thar2fdate
          ts=ts2f
          wdt=wdt2f
          resel=resel2f
        endif else begin
          print,''
          if lsci2f ne !NULL then print,'ERROR: missing calibrations for the 2-fiber mode'
          goto,fin
        endelse
      end
    endcase
    if strcmp(flatdate,thardate) ne 1 then begin
      print,''
      print,'WARNING!!! The flatfield image and the ThAr spectrum are not from the same date.'
      print,'           The real geometry on the 2D spectra may not correspond to the correction applied.'
      print,''
    endif



    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ;  FINDING THE TRACES
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

    ;Relies on the function "find_trace"
    if lflat ne !NULL then begin
      if file_test(reddir+'TraceCoeff_'+spmd+flatdate+'.dat') ne 1 or keyword_set(new) eq 1 then begin
        ;calls find_trace
        matcoeff=find_trace(flat,spmd,ts,wdt,vcen,nco)
        ;records the result, so we do not need to recalculate it next time we reduce spectra from the same date
        openw,lun,reddir+'TraceCoeff_'+spmd+flatdate+'.dat',/get_lun  ;the traces (x-position) are recorded here
        for i=0,nord-1 do begin
          temp=''
          for j=0,nco do temp=temp+strtrim(string(double(matcoeff[j,i])),2)+' '
          printf,lun,format='(a)',temp
        endfor
        close,lun
        free_lun,lun
      endif else begin
        ;reads the values if it has already been calculated
        openr,lun,reddir+'TraceCoeff_'+spmd+flatdate+'.dat',/get_lun
        matcoeff=dblarr(nco+1,nord)
        readf,lun,matcoeff
        close,lun
        free_lun,lun
      endelse
    endif



    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ; FINDING THE SLIT TILT
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

    ;Relies on the function "find_lines"
    if file_test(reddir+'Tilt_'+spmd+thardate+'.dat') ne 1 or keyword_set(new) eq 1 then begin
      ;reduces the ThAr frame, without applying a tilt
      tilt_2d=reduce(thar,matcoeff,wdt,mask=nonlin,/notilt)
      ;extracts the ThAr frame, so we can find the lines
      tilt_1d=extract(tilt_2d,wdt,spmd)
      ;calls find_lines
      slit_tilt=find_lines(tilt_1d,tilt_2d,spmd,wdt,resel,nonlin,/find_tilt)

      ;records the result, so we do not need to recalculate it next time we reduce spectra from the same date
      openw,lun,reddir+'Tilt_'+spmd+thardate+'.dat',/get_lun  ;the tilts are recorded here
      for i=0,nord-1 do begin
        temp=[strtrim(string(double(slit_tilt[0,i])),2),strtrim(string(double(slit_tilt[1,i])),2)]
        printf,lun,format='(a,a,a)',temp[0],'  ',temp[1]
      endfor
      close,lun
      free_lun,lun
    endif else begin
      ;reads the values if it has already been calculated
      openr,lun,reddir+'Tilt_'+spmd+thardate+'.dat',/get_lun
      slit_tilt=fltarr(2,nord)
      readf,lun,slit_tilt
      close,lun
      free_lun,lun
    endelse



    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ; REDUCTION OF CALIBRATIONS
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    wdt=wdt-1 ;the width is reduced by one pixel at this point; easier for reduction

    ;runs reduce on the ThAr wih the aim to produce an .eps figure showing where the
    ;  traces are extracted, and with which slit tilt.
    junk=reduce(thar,matcoeff,wdt,slit_tilt,mask=nonlin,disp=spmd+thardate,dir=reddir)

    ;Gets the normalized flat field
    flatN=reduce(flat,matcoeff,wdt,slit_tilt,mask=nonlin,/normal)
    writefits,reddir+'FlatNorm'+spmd+flatdate+'.fits',flatN;,hf



    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ; LINE IDENTIFICATION
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ;Extracts the ThAr spectrum for final line identification
    wavel2d=reduce(thar,matcoeff,wdt,slit_tilt)
    wavel2d=wavel2d/flatN
    wavel_sp=extract(wavel2d,wdt,spmd)
    if keyword_set(skip_wavel) then begin
      ;Saves the extracted ThAr spectrum
      s=size(wavel_sp)
      fits_open,reddir+'red_ThAr_sp'+thardate+spmd+'.fits',un,/write
      for i=0,nord-1 do begin
        spectrum=total(wavel_sp[i,*],1)
        mkhdr,h,spectrum,/image
        sxaddpar,h,'NAXIS',1
        sxaddpar,h,'NAXIS1',s[2]
        sxdelpar,h,'NAXIS2'
        fits_write,un,spectrum,h,extname='order '+strtrim(string(i+first_order),2),extver=i+1
      endfor
      fits_close,un
    endif else begin
      ;Finds the pix position of strongest lines
      line_list=find_lines(wavel_sp,wavel2d,spmd,wdt,resel,nonlin)
      wavel_sol_per_order=wavel_sol(line_list,spmd,reddir,vec_lines_used=vec_lines_used)
      ;displays the identified lines in a multiple pages .ps file.
      ; note that most of the follow variables are tuned to give acceptable plots.
      set_plot,'ps'
      device,filename=reddir+'Identification'+spmd+thardate+'.ps',/color;,/landscape
      pos=where(vec_lines_used[0,*] eq -1)
      for i=0,nord-1 do begin
        if i mod 5 eq 0 then begin
          multiplot,/default
          multiplot,[1,5],ygap=.01
        endif
        spectrum=wavel_sp[i,*]
        pix_pos=vec_lines_used[0,pos[i]+1:pos[i+1]-1]
        wav_pos=vec_lines_used[1,pos[i]+1:pos[i+1]-1]
        ylable=spectrum[pix_pos]
        posy=where(ylable gt 8000)
        if max(posy) ne -1 then ylable[posy]=8000
        if (i+1) mod 5 eq 0 then xttl='pixels' else xttl=''
        plot,spectrum,yrange=[0,10000],xtitle=xttl,ytitle='ADU',title='Order '+strtrim(string(i+first_order),2),/xst,xrange=[0,4300]
        plotsym,1
        oplot,pix_pos,ylable,psym=8,color=255
        for j=0,pos[i+1]-pos[i]-2 do xyouts,pix_pos[j]+5,ylable[j]+20,strtrim(string(wav_pos[j]),2),color=255,orientation=90,charsize=0.7
        multiplot
        if (i+1) mod 5 eq 0 then erase
      endfor
      multiplot,[1,1]
      device,/close
      set_plot,'x'
    endelse

    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ; EXTRACTION
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

    ;LOOP ON SCIENCE FRAMES to get list of objects
    for j=0,n_elements(lsci)-1 do begin
      h=headfits(lsci[j])
      temp_name=strcompress(sxpar(h,'OBJECT'),/remove)
      if j eq 0 then obs_name=temp_name else obs_name=[obs_name,temp_name]
    endfor
    ;sorts lsci by object name (obs_name) to make sure they are grouped by object
    s_ind=sort(obs_name)
    obs_name=obs_name[s_ind]
    lsci=lsci[s_ind]
    ;number of different objects oberved
    nobj=obs_name[uniq(obs_name)]
    
    ;Loop on each object
    for j=0,n_elements(nobj)-1 do begin
      print,'Extracting ',nobj[j]
      pos_obj=where(nobj[j] eq obs_name)
      ;arrays to store keywords' values
      v_RA=fltarr(n_elements(pos_obj))
      v_DEC=fltarr(n_elements(pos_obj))
      v_EPOCH=fltarr(n_elements(pos_obj))
      v_EXPTIME=fltarr(n_elements(pos_obj))
      v_AIRMASS=fltarr(n_elements(pos_obj))
      v_OBSID=strarr(n_elements(pos_obj))
      v_DATE=fltarr(n_elements(pos_obj))
      v_MJD=fltarr(n_elements(pos_obj))
      for k=0,n_elements(pos_obj)-1 do begin
        ;reads the science frame
        im=readfits(lsci[pos_obj[k]],h,/silent)
        ;adds the bzero value
        im+=float(sxpar(h,'BZERO'))
        ;corrects for the bias level using the overscan
        im=overs_corr(im,h)
        ;picks the right bias, readout noise, gain
        if sxpar(h,'RDNOISEA') eq 4.2 then begin
          bias_i=biasNr
          gain=1.3 & readn=4.2
        endif else begin
          bias_i=biasSr
          gain=1.2 & readn=2.9
        endelse
        ;corrects for the bias structure (usually pretty flat) using the bias frame
        im=im-bias_i
        ;Stacking images for cosmic ray removal
        if k eq 0 then begin
          dim_im=size(im)
          mat_im=fltarr(dim_im[1],dim_im[2],n_elements(pos_obj))
          lsci_obj=strarr(n_elements(pos_obj))
        endif
        mat_im[*,*,k]=im
        lsci_obj[k]=lsci[pos_obj[k]]
        v_AIRMASS[k]=sxpar(h,'AIRMASS')
        v_RA[k]=sxpar(h,'RA')
        v_DEC[k]=sxpar(h,'DEC')
        v_EPOCH[k]=sxpar(h,'EPOCH')
        v_EXPTIME[k]=sxpar(h,'EXPTIME')
        v_OBSID[k]=sxpar(h,'OBSID')
        v_DATE[k]=sxpar(h,'DATE')
        v_MJD[k]=sxpar(h,'MJDATE')
      endfor
      if n_elements(pos_obj) ge 2e6 then begin
        ;Cosmic ray rejection
        cr_reject,mat_im,readn,0.,gain,1.0,comb_im,mask_cube=mask_cube,/noskyadjust,/bias,/init_med,/median_loop
        for k=0,n_elements(ls_temp)-1 do begin
          im_temp=mat_im[*,*,k]
          im_mask=mask_cube[*,*,k]
          pos_cr=where(im_mask lt 1)
          if max(pos_cr) ne -1 then im_temp[pos_cr]=comb_im[pos_cr] ;replaces detected CR by the pixel value in the combined image
          mat_im[*,*,k]=im_temp
        endfor
      endif
      for k=0,n_elements(pos_obj)-1 do begin
        ;reduces the 2d spectrum
        sci2d=reduce(mat_im[*,*,k],matcoeff,wdt,slit_tilt)
        ;corrects the illumination
        sci2d=illumcorr(mat_im[*,*,k],sci2d,matcoeff,wdt,slit_tilt)
        ;corrects for the flat field
        sci2d=sci2d/flatN
        ;extracts the spectrum
        sci=extract(sci2d,wdt,spmd)
        
        s=size(im)
        ;would save the extracted spectrum if the wavelength calibration is skipped
        if keyword_set(skip_wavel) then begin
          fits_open,reddir+'red_'+strmid(lsci_obj[k],strlen(datadir)),un,/write
          ;would calibrate for the wavelength and save the final spectrum
        endif else begin
          final=wavel_cal(sci,wavel_sol_per_order,param=param)
          fits_open,reddir+'ext_'+strmid(lsci_obj[k],strlen(datadir)),un,/write
        endelse
;        mkhdr,h,sci
;        sxaddpar,h,'AIRMASS',v_AIRMASS[k]
;        sxaddpar,h,'RA',v_RA[k]
;        sxaddpar,h,'DEC',v_DEC[k]
;        sxaddpar,h,'EPOCH',v_EPOCH[k]
;        sxaddpar,h,'EXPTIME',v_EXPTIME[k]
;        sxaddpar,h,'OBSID',v_OBSID[k]
;        sxaddpar,h,'DATE',v_DATE[k]
;        sxaddpar,h,'MJDATE',v_MJD[k]
;        fits_write,un,sci,h,/no_data,extver=0
        
        for i=0,nord*idx-1 do begin
          ;if the spectrum is NOT calibrated
          if keyword_set(skip_wavel) then begin
            spectrum=total(sci[i,*],1)
            mkhdr,h,spectrum,/image
            sxaddpar,h,'NAXIS',1
            sxaddpar,h,'NAXIS1',s[2]
            sxdelpar,h,'NAXIS2'
            ;if the spectrum is calibrated
          endif else begin
            spectrum=total(final[i,*],1)
            mkhdr,h,spectrum,/image
            sxaddpar,h,'NAXIS',1
            sxaddpar,h,'NAXIS1',s[2]
            sxdelpar,h,'NAXIS2'
            sxaddpar,h,'CTYPE1','WAVE-WAV-PLY'
            sxaddpar,h,'CUNIT1','Angstrom'
            sxaddpar,h,'DC-FLAG',0
            sxaddpar,h,'CRPIX1',1
            sxaddpar,h,'CRVAL1',param[0,i]
            sxaddpar,h,'CD1_1',param[1,i]
          endelse
          if spmd eq 1 then ordname=strtrim(string(i+first_order),2) else if i mod 2 eq 0 then ordname=strtrim(string(i/2+first_order),2)+' sky' else ordname=strtrim(string(i/2+first_order),2)+' target'
          fits_write,un,spectrum,h,extname='order '+ordname,extver=i+1
        endfor
        fits_close,un
      endfor
    endfor
    skip_loop:
  endfor
  print,'Please cite the AJ paper:'
  print,'https://ui.adsabs.harvard.edu/abs/2021AJ....161..109C/abstract'

  fin:
end




