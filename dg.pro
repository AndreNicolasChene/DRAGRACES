; GRACES pipeline, by Andre-Nicolas Chene

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
  print,'        .,,,,,,,,,,,,,,,...�, ,..,,,.     ,.                ...../,#,,,,%%/,,,,**,..*/.                         .,.,               ,.,..,/,(@@#(,           '
  print,'         *...,,,,,,....      ,, ..,,.   .*,                 ,...,/**...*.                     ...,,*****,,...   .,*               ,...../,#,,..(&/          '
  print,'       ,.,,,,....,,,,,,,,..,    ..,,,. .*.                  ,,,,&//,*                                                ,,  ...,******,...,/**......@#         '
  print,'        ..                   ,, ..,,.   .,                  @@@@@(***                                                             ,,,,*@//,       @%        '
  print,'       ,((/.                        ,#''+#,                  &@@@@%(%*                          *,*,*.*                            .@@@@@(//.      .@.       '
  print,'  ''@@*#@\#%&%                  ,`@@./%##&#                 .@@@@@#@*,...                   ,*****,*.*                             %@@@@&(%.. .    %&       '
  print,'  @#@@/(**&&*,%%               ,:@@(/#,&&%#%%                ,@@@@@%@&.....,&@.....,*/////***/**////*,...,,                         @@@@@(@&...   .#@       '
  print,' (*@@@&&.+�`.&,%%              /@@@@&%,`.+#%%,,*********,,,..,@@@@&&@@#,,*/##*,.                            ?,,.....,,,,****//*****/@@@@@%@@......&@.       '
  print,'  ,@@@%&++...*%.&+............. @@@@&&/+...#(/&+                .####%%@@@@@@@(                                                       .@@@@@%@@&,,,%.       '
  print,'   (@@@%&+,(&&%.&%              ,@@@@%&%.#&&(/&%                ,@@@#&@@@@(.                                                         (@@@@#@@@@@@@/         '
  print,'    .(@@/%&/#&&%                 ,/@@@&&%#&&&*                       .,*/*,.                                                             *&@@@@@@*          '
  print,'      .****(%&&.                  .*///#%&&%*                                                                                                .,,,.          '
  print,''
endif
print,''
print,'~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*'
print,''
print,'Version 1.0.1'
print,'Author: Andre-Nicolas Chene'
print,'Release date: 20 September 2016'
print,''
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
  print,'       dir - Path to where the data can be found. A new directory named ""Reduction"" will'
  print,'             be created, if it does not already exist.'
  print,'    utdate - Date when the spectra where observed. The format is YYYYMMDD, and can be'
  print,'             given as a number (without ''''). If it not provided, it is expected that the'
  print,'             data in the directory pointed by the input ""dir"" are all from the same date. '
  print,'     lbias - Name of a file where the bias frames are listed.'
  print,'     lflat - Name of a file where the flats frames are listed.'
  print,'     lthar - Name of a file where the ThAr frames are listed.'
  print,'skip_wavel - When provided, the wavelength solution is not calculated.'
  print,'ascii_file - (NOT YET WORKING)'
  print,'       new - Forces DRAGRACES to recalculate all the calibrations instead of using those'
  print,'             obtained in a previous run of the pipeline'
  print,'      logo - Displays the DRAGRACES logo when you run the pipeline!'
  print,'      help - Displays this help summary.'
  print,'OUTPUTS:'
  print,'    The program saves the extracted spectra in the fits format, in the directory ""Reduction"".'
  print,'    If the wavelength solution was calculated, the filenames start with ""ext_"". If not, '
  print,'    the filenames start with ""red_"", and an additional ""red_Thar_..." spectrum is provided.'
  print,''
  print,'IMPORTANT NOTE:'
  print,''
  print,'    The pipeline expets the files to be ""well behaved"". It expects all the files in a given'
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
if keyword_set(dir) then datadir=dir else datadir='./'
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
overs=31     ;size of the overscan in pixel


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
;vcen=[38, 68,  99,  129,  160,  191,  223,  254,  287,  320,  354,  389, 424, 461, 499, 537, 577, 618, 660, 704, 748, 794, 842, 891, 941, 993, 1047, 1102,  1159,  1218,  1278,  1341,  1405,  1472,  1541, 1612, 1686] ;from order 21 to 57
vcen=[68,  99,  129,  160,  191,  223,  254,  287,  320,  354,  389, 424, 461, 499, 537, 577, 618, 660, 704, 748, 794, 842, 891, 941, 993, 1047, 1102,  1159,  1218,  1278,  1341,  1405,  1472,  1541, 1612] ;from order 22 to 56
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
  readcol,lbias,format='(a)',lbias
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
    ;IMPORTANT NOTE: No rejection method used at this point!
    for i=0,n_elements(lbiasNr)-1 do begin
      im=readfits(lbiasNr[i],/silent) 
      im_overs=overs_corr(im,overs)
      if i eq 0 then biasNr=im_overs/n_elements(lbiasNr) else biasNr=biasNr+im_overs/n_elements(lbiasNr)
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
    ;IMPORTANT NOTE: No rejection method used at this point!
    for i=0,n_elements(lbiasSr)-1 do begin
      im=readfits(lbiasSr[i],/silent) 
      im_overs=overs_corr(im,overs)
      if i eq 0 then biasSr=im_overs/n_elements(lbiasSr) else biasSr=biasSr+im_overs/n_elements(lbiasSr)
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
  readcol,lthar,format='(a)',lthar
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
    ;IMPORTANT NOTE: No rejection method used at this point!
    for i=0,n_elements(lthar1f)-1 do begin
      im=readfits(lthar1f[i],/silent) 
      im_overs=overs_corr(im,overs)
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
    ;IMPORTANT NOTE: No rejection method used at this point!
    for i=0,n_elements(lthar2f)-1 do begin
      im=readfits(lthar2f[i],/silent) 
      im_overs=overs_corr(im,overs)
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
  readcol,lflat,format='(a)',lflat
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
    ;IMPORTANT NOTE: No rejection method used at this point!
    for i=0,n_elements(lflat1f)-1 do begin
      im=readfits(lflat1f[i],/silent) 
      im_overs=overs_corr(im,overs)
      if i eq 0 then flat1f=im_overs/n_elements(lflat1f) else flat1f=flat1f+im_overs/n_elements(lflat1f)
    endfor
    flat1f=flat1f-biasNr
    writefits,reddir+'Flat1f'+flat1fdate+'.fits',flat1f;,hf
  endelse
endif
;checks if there was any Flat observed in 2-fiber mode.
if lflat2f ne !NULL then begin
  ;gets the date from the header of the first flat
  ht=headfits(lflat2f[0])
  date=sxpar(hf,'DATE')
  flat2fdate=strmid(date,0,4)+strmid(date,5,2)+strmid(date,8,2)
  
  ;checks first if the Master has not been created before
  if file_test(reddir+'Flat2f'+flat2fdate+'.fits') and keyword_set(new) ne 1 then begin
    ;if a Master Flat has already been processed, it takes this one.
    flat2f=readfits(reddir+'Flat2f'+flat2fdate+'.fits',/silent)
  endif else begin
    ;IMPORTANT NOTE: No rejection method used at this point!
    for i=0,n_elements(lflat2f)-1 do begin
      im=readfits(lflat2f[i],/silent) 
      im_overs=overs_corr(im,overs)
      if i eq 0 then flat2f=im_overs/n_elements(lflat2f) else flat2f=flat2f+im_overs/n_elements(lflat2f)
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
    print,'           The real geometry on the 2D spectra may not correspond to the correcion applied.'
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
  ;  traces are extracted, ans with which slit tilt.
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
    wavel_sol_per_order=wavel_sol(line_list,spmd,vec_lines_used=vec_lines_used)
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
    multiplot,/reset
    device,/close
    set_plot,'x'
  endelse

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; EXTRACTION
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ;LOOP ON SCIENCE FRAMES
  for j=0,n_elements(lsci)-1 do begin
    ;reads the science frame
    im=readfits(lsci[j],h,/silent)
    ;picks the right bias
    if sxpar(h,'RDNOISEA') eq 4.2 then bias_i=biasNr else bias_i=biasSr
    ;corrects for the bias level using the overscan
    im=overs_corr(im,overs)
    ;corrects for the bias structure (usually pretty flat) using the bias frame
    im=im-bias_i
    ;reduces the 2d spectrum
    sci2d=reduce(im,matcoeff,wdt,slit_tilt)
    ;corrects for the flat field
    sci2d=sci2d/flatN
    ;extracts the spectrum
    sci=extract(sci2d,wdt,spmd)
  
    s=size(im)
    ;would save the extracted spectrum if the wavelength calibration is skipped
    if keyword_set(skip_wavel) then begin
      fits_open,reddir+'red_'+strmid(lsci[j],strlen(datadir)),un,/write
    ;would calibrate for the wavelength and save the final spectrum
    endif else begin
      final=wavel_cal(sci,wavel_sol_per_order,param=param)
      fits_open,reddir+'ext_'+strmid(lsci[j],strlen(datadir)),un,/write
    endelse
    for i=0,nord-1 do begin
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
      fits_write,un,spectrum,h,extname='order '+strtrim(string(i+first_order),2),extver=i+1
    endfor
    fits_close,un
  endfor
  
  skip_loop:
endfor

fin:
end



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;  FUNCTION APPLYING THE OVERSCAN CORRECTION (nice, but not necessary)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function overs_corr,im,overs
s=size(im)
len=s[2] ;numer of pixels in one row

;Gets the overscan value for each row
vec_overs=fltarr(len)
for i=0,len-1 do vec_overs[i]=mean(im[s[1]-overs:s[1]-1,i])

;Fits a slope to the overscan values as a function of row #
res=ladfit(findgen(len),vec_overs)

;Creates the overscan corrected resulting image
im_corr=fltarr(s[1]-overs,len)
for i=0,len-1 do im_corr[*,i]=im[0:s[1]-overs-1,i]-(res[0]+res[1]*i)

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
  ref=.7*exp(-1*((a-4.0)^2/(2*2.^2)))+exp(-1*((a-10.0)^2/(2*2.^2)))+.9*exp(-1*((a-20.5)^2/(2*2.^2))+exp(-1*((a-27.0)^2/(2*2.^2))))
  vcen=vcen-6 ;slight shift measured in the position of the order with respect to 1-fiber mode
endelse

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

lineHW=fix(resel*2)+1 ;value of a line half width at half max

;For each order, it finds line on the 1d spectrum, and verifies if they
;  correspond to real lines on the 2d spectrum
px=findgen(s[2])
slit_tilt=fltarr(2,s[1])
line_list=-1
for i=0,s[1]/step-1 do begin
  sp=sp_1d[i*step,*]           ;1d spectrum of that given order
  
  ;Finds a rough coninuum to the ThAr spectrum
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
  
  ;Finds the lines and measures tilt
  vec_line=-1 ;vector for lines pix-position
  vec_tilt=-1 ;vector for slit tilt values
  ;we loop until we cannot find lines higher than 500ADU
  while max(spN) ge 500 do begin ;500 is an semi-educated, but mostly arbitrary value
    ;pos is the position of the highest peak (line?)
    bid=max(spN,pos)
    ;would skip if the line is too close to the ends of the spectrum
    if pos gt lineHW and pos lt s[2]-lineHW then begin
      ;would skip if a line has already been identified within the resolution element
      if max(where(spN[pos-lineHW:pos+lineHW] eq 0)) eq -1 then begin
        ;fits the line
        res=gaussfit(px[pos-lineHW:pos+lineHW],spN[pos-lineHW:pos+lineHW],nterms=3,coeff)
        ;would skip if the center if found way off bounds
        if coeff[1] gt px[pos]-lineHW/2. and coeff[1] lt px[pos]+lineHw/2. then begin
          ;samples the 2d ThAr spectrum around the identified potential line
          yb=coeff[1]+lineHW*[-1,1]
          sample=sp_2d[i*wdt:(i+1)*wdt-1,yb[0]:yb[1]]
          ;uses the extracted shape od the line as reference (on non-continuum corrected spectrum)
          ref=sp[pos-lineHW:pos+lineHW]
          ;would skip if a pixel-value is over the non-linearity threshold (in ADU)
          psat=where(sample[*,lineHw/2:lineHw*1.5] gt nonlin)
          if max(psat eq -1) then begin
            ;compares the profile of the 2d line with an approximative model to 
            ;  verify if it is a real line or not (model is spmd dependent)
            ytest=total(sample,2)
            ytest=ytest-min(ytest)
            go=0
            ;NOTE: this method could be improved
            case spmd of
              '1f': begin
                pnt1=max(ytest[2:5])
                pnt2=min(ytest[4:8])
                pnt3=max(ytest[9:12])
                pnt4=min(ytest[12:15])
                pnt5=max(ytest[16:19])
                pnt6=min(ytest[20:22])
                pnt7=max(ytest[23:26])
                if (pnt1 gt pnt2*1.2) and (pnt2 lt pnt3*.75) and (pnt3 gt pnt4*1.25) and (pnt4 lt pnt5*.75) and (pnt5 gt pnt6*1.25) and (pnt6 lt pnt7*.9) then go=1
              end
              '2f': begin
                pnt1=max(ytest[2:6])
                pnt2=min(ytest[5:8])
                pnt3=max(ytest[8:12])
                pnt4=min(ytest[13:17])
                pnt5=max(ytest[18:22])
                pnt6=min(ytest[21:25])
                pnt7=max(ytest[26:29])
                if (pnt1 gt pnt2*1.1) and (pnt2 lt pnt3*.75) and (pnt3 gt pnt4*1.25) and (pnt4 lt pnt5*.65) and (pnt5 gt pnt6*1.1) and (pnt6 lt pnt7*.9) then go=1
              end
            endcase
            
            ;would skip if the line profile does not match the (poor) model
            if go eq 1 then begin
              ;records the line pixel-position
              if max(vec_line) eq -1 then vec_line=coeff[1] else vec_line=[vec_line,coeff[1]]
              
              ;would skip if we are not interested to the tilt
              if keyword_set(find_tilt) then begin
                ;cross-correlates the shape of the line accross "sample"
                vcen=fltarr(wdt)
                for j=0,wdt-1 do begin
                  lag=findgen(9)-4
                  b=sample[j,*]
                  a=findgen(n_elements(b))
                  cc=c_correlate(b,ref,lag)
                  res=gaussfit(lag,cc-min(cc),co,nterms=3)
                  vcen[j]=-co[1]
                endfor
                
                ;fits the slit tilt
                if step eq 1 then begin
                  res=ladfit(findgen(n_elements(vcen)),vcen)
                  if max(vec_tilt) eq -1 then vec_tilt=res[1] else vec_tilt=[vec_tilt,res[1]]
                endif else begin
                  res1=ladfit(findgen(wdt/2),vcen[0:wdt/2-1])
                  res2=ladfit(findgen(round(wdt/2.)),vcen[wdt/2:n_elements(vcen)-1])
                  if max(vec_tilt) eq -1 then vec_tilt=[(res1[1]+res2[1])/2] else vec_tilt=[vec_tilt,(res1[1]+res2[1])/2]
                endelse
              endif
            endif
          endif
        endif
      endif
      ;"kills" the identified (potential) line by forcing the corresponding part of the spectrum to 0
      spN[pos-lineHW:pos+lineHW]=0
    endif else if pos le lineHW then spN[0:pos+lineHW]=0 else spN[pos-lineHW:s[2]-1]=0
  endwhile
  
  ;would skip if we are not interested to the tilt
  if keyword_set(find_tilt) then begin
    ;removes deviant measurement of tilt (3sig)
    pos=where(abs(vec_tilt-median(vec_tilt)) gt 3*stddev(vec_tilt))
    while max(pos) ne -1 do begin
      remove,pos,vec_line,vec_tilt
      pos=where(abs(vec_tilt-median(vec_tilt)) gt 3*stddev(vec_tilt))
    endwhile
    ;fits the slit tilt as a function of pixel (gentle slope)
    res=ladfit(vec_line,vec_tilt)
    slit_tilt[*,i]=res
  endif
  
  if max(vec_line) ne -1 then for j=0,n_elements(vec_line)-1 do if max(line_list) eq -1 then line_list=[i,vec_line[j]] else line_list=[[line_list],[i,vec_line[j]]]
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
function wavel_sol,line_list,spmd,vec_lines_used=vec_lines_used

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
ThAr_list=$
[3006.9314,3012.7078,3015.7252,3017.1322,3018.4899,3019.4219,3030.4877,3033.5083,3038.5987,3046.9513,3047.8163,3049.0924,3060.4389,3067.7294,3068.9075,3072.1150,3072.8252,3075.0565,3078.8280,3080.2170,3090.0937,3093.4019,3097.2672,3100.9382,3102.6639,3105.7474,3107.0261,3108.2968,3110.0209,3116.2630,3117.6840,3119.5262,3120.8808,3122.9634,3124.3874,3134.4257,3136.2161,3136.8293,3140.2720,3141.8472,3142.8356,3145.6370,3146.0434,3147.4486,3150.4551,3151.6466,3154.3009,3155.8302,3156.4006,3161.3726,3167.5606,3169.6685,3171.2755,3172.1112,3173.4266,3174.2043,3175.7257,3176.1671,3178.2433,3179.0482,3180.1937,3181.0376,3181.6722,3204.3210,3210.3088,3211.1943,3212.0220,3213.5741,3215.7781,3220.3508,3229.0096,3230.8679,3235.8400,3238.1157,3238.9343,3243.6887,3244.4488,3245.7603,3249.8003,3251.9159,3252.9644,3253.8659,3254.8099,3256.2738,3257.3667,3258.1054,3259.0608,3259.6164,3260.9417,3262.6684,3265.5778,3271.1199,3271.8912,3272.0268,3273.9157,3275.0676,3276.1712,3280.3713,3281.0488,3281.2731,3281.7016,3282.6166,3282.9795,3285.7525,3286.5829,3287.7893,3290.5996,3291.7394,3292.5209,3293.2257,3293.9481,3294.2452,3296.6073,3297.3748,3297.8333,3298.0498,3299.6687,3301.6511,3304.2383,3305.3036,3307.2283,3308.4740,3309.3654,3310.1988,3313.6784,3314.8268,3319.9101,3320.3004,3320.4763,3321.4508,3321.5741,3322.0933,3324.7527,3325.1207,3326.4652,3327.1931,3329.7284,3330.4770,3333.1290,3334.6041,3335.0639,3336.1618,3337.8703,3338.3970,3340.7254,3343.1644,3343.6183,3345.1709,3345.8804,3346.5560,3346.9648,3347.9923,3348.7684,3350.3514,3350.9243,3351.2286,3353.9484,3354.1796,3354.6185,3355.1061,3355.2587,3355.5631,3358.6020,3359.7569,3360.3735,3360.9982,3361.1979,3361.6190,3361.7367,3362.6762,3364.6855,3365.1376,3365.3383,3366.5171,3367.5822,3367.8189,3371.7967,3372.8230,3373.4925,3374.5814,3374.9749,3376.4359,3378.5734,3380.8595,3383.1068,3385.5316,3386.5006,3387.9205,3388.5309,3389.4632,3389.6405,3390.3681,3391.8713,3392.0349,3393.4215,3393.9929,3394.7970,3396.7278,3397.5161,3398.5448,3401.7110,3402.0264,3402.6952,3404.6508,3405.5584,3406.2418,3407.8302,3408.7499,3409.2699,3410.0756,3413.0130,3415.8846,3417.4978,3418.7755,3419.1733,3421.2100,3422.6561,3423.1290,3423.9897,3425.4348,3425.9436,3427.0923,3428.6221,3428.7144,3428.9992,3430.4175,3431.8104,3433.9988,3434.7271,3435.9771,3436.7272,3437.3071,3438.9503,3439.3987,3439.7118,3441.0359,3441.3647,3441.5268,3442.5790,3445.2173,3445.7441,3446.5474,3448.9487,3449.2871,3449.6448,3451.7023,3452.6820,3454.0952,3455.6130,3457.0691,3461.0187,3461.2172,3462.8505,3463.7197,3464.1272,3465.0626,3465.7650,3466.5383,3466.6459,3468.2198,3469.3454,3469.9208,3470.5675,3471.0005,3471.2186,3471.9593,3476.3852,3476.7474,3477.7042,3478.2324,3478.4635,3479.1725,3479.6847,3480.0525,3480.5055,3482.5482,3482.7613,3484.0800,3485.2124,3486.5512,3487.8436,3488.8338,3489.1841,3489.5076,3490.4524,3490.8733,3491.2439,3491.5360,3491.9000,3493.2513,3493.5185,3495.6998,3496.0604,3496.8107,3497.2628,3498.0098,3498.6210,3499.6829,3499.8218,3499.9866,$
3500.1032,3501.4559,3501.8666,3502.9637,3503.6144,3503.7859,3505.4941,3506.1328,3506.6458,3507.5464,3509.0889,3509.7785,3511.1574,3511.5620,3512.7425,3514.3877,3514.9621,3516.3545,3516.8243,3518.4040,3518.8858,3519.9936,3521.0595,3521.2601,3523.5061,3523.7585,3524.1787,3524.7091,3526.6342,3527.0005,3527.3224,3528.4116,3528.8202,3528.9544,3529.3859,3530.5148,3531.4505,3533.1826,3535.3196,3536.0108,3537.1597,3539.3223,3539.5872,3539.8401,3541.6159,3542.4979,3543.1475,3544.0179,3545.2851,3545.5956,3545.8450,3547.3376,3547.9175,3548.5144,3549.5959,3550.7184,3551.4019,3553.1103,3554.3058,3555.0135,3555.7046,3556.3144,3556.9041,3557.4645,3559.5081,3561.0304,3561.7809,3562.1930,3565.0298,3565.6042,3567.0471,3567.2638,3567.6564,3569.6204,3569.8204,3570.3577,3570.5239,3571.5731,3572.3923,3573.2196,3575.1263,3575.3018,3575.3224,3576.6156,3577.5994,3578.9497,3579.4575,3580.2314,3581.6084,3581.7579,3582.0091,3582.3546,3583.1022,3584.1756,3585.0509,3585.7702,3588.4407,3589.1538,3589.3603,3589.7495,3589.9941,3590.5250,3590.9253,3591.4524,3592.7794,3594.1111,3594.7216,3594.9856,3595.6177,3597.4955,3598.1199,3598.5247,3599.7240,3600.4323,3601.0344,3601.5093,3601.7716,3603.2042,3604.0656,3604.6819,3605.1936,3605.8792,3606.0909,3606.2481,3606.5218,3607.3933,3608.3770,3609.2251,3609.4452,3610.3971,3611.1604,3612.4275,3612.8666,3613.7795,3614.0097,3614.2133,3614.3534,3615.1327,3615.8503,3616.1778,3617.1173,3617.6716,3618.3633,3619.2119,3620.8387,3621.1187,3622.1375,3622.3358,3622.7954,3623.7725,3624.4720,3624.8953,3625.0293,3625.1496,3625.6280,3625.8933,3626.9390,3629.8503,3632.8303,3634.5822,3634.8124,3635.2419,3635.4195,3635.9433,3636.1728,3636.5667,3636.8345,3637.0310,3637.5557,3638.3192,3638.6444,3639.4469,3639.8329,3641.9647,3642.2490,3642.5729,3643.5123,3643.8250,3644.4075,3644.7138,3645.3638,3645.7065,3646.2677,3647.6455,3648.1741,3648.4212,3649.2496,3649.7349,3650.8896,3651.5716,3652.1683,3652.5372,3654.4618,3655.1256,3655.2782,3656.0498,3656.2006,3656.6938,3657.0554,3657.6427,3658.0637,3658.2275,3658.8087,3659.6294,3660.4370,3661.6214,3661.9833,3662.2855,3662.7504,3663.2025,3663.5444,3663.7023,3665.1811,3665.4407,3666.9811,3667.6218,3668.1398,3668.3901,3669.6024,3669.9684,3670.6370,3671.0048,3671.5398,3672.3003,3672.5218,3673.2645,3673.7935,3674.0147,3674.8910,3675.1372,3675.5675,3675.7893,3675.9592,3676.1424,3676.6897,3677.9139,3678.0483,3678.2701,3678.4804,3679.1343,3679.5444,3679.7105,3680.0609,3680.4477,3680.6053,3681.8836,3682.1793,3682.4863,3683.4931,3684.9329,3685.5870,3686.9812,3687.1943,3687.4936,3687.6700,3687.9841,3688.7604,3690.1158,3690.6238,3691.4117,3691.6141,3691.8759,3692.0763,3692.5664,3693.2473,3693.5238,3693.9946,3694.1785,3694.3651,3694.7854,3695.2889,3695.9737,3696.6468,3697.0308,3697.7436,3698.1061,3698.7545,3699.1811,3699.8808,3700.3118,3700.7653,3700.9872,3702.8638,3703.2299,3703.7743,3704.3373,3704.7103,3704.8616,3706.0238,3706.7672,3707.0036,3707.4282,3708.7530,3709.8620,3711.3041,3711.6229,3711.8332,3712.5591,3714.0705,3714.3254,$
3715.5606,3715.8615,3716.5836,3717.1713,3717.8293,3718.2065,3719.4347,3719.8362,3719.9772,3720.3073,3720.4265,3721.2152,3721.8254,3722.1154,3723.2898,3723.6561,3723.9199,3724.5165,3725.3932,3726.7246,3727.6120,3727.9027,3729.3087,3729.8361,3729.9459,3730.3683,3730.7485,3732.9854,3733.6724,3734.5977,3735.5109,3736.9026,3737.5125,3737.8890,3738.8439,3740.8551,3741.1830,3742.2773,3742.9234,3743.5083,3743.8958,3744.7367,3744.9109,3745.1763,3745.6591,3745.9706,3746.9133,3747.5390,3748.2838,3749.0843,3749.6183,3750.4939,3750.6602,3751.0219,3752.5689,3753.2421,3753.5177,3754.0308,3754.2884,3754.5930,3755.2121,3755.9869,3756.2941,3757.6941,3758.4671,3758.7063,3759.4342,3760.2724,3761.1019,3761.3164,3761.7044,3762.4166,3763.5053,3765.2700,3766.1186,3766.4473,3767.0819,3767.9007,3768.4359,3769.5854,3770.0560,3770.3687,3770.5200,3771.3708,3771.6150,3772.2369,3772.6498,3773.0521,3773.7574,3774.2037,3774.7331,3775.4234,3775.9028,3776.2711,3776.6239,3777.4167,3777.7473,3779.5624,3780.8398,3780.9662,3781.3189,3783.0127,3783.2964,3783.4936,3784.5749,3784.7942,3785.2804,3785.6002,3786.3824,3786.883,3789.168,3790.356,3790.795,3792.374,3792.730,3794.151,3794.318,3795.3800,3798.103,3799.3820,3800.198,3801.443,3803.0750,3803.984,3805.314,3805.820,3807.874,3808.129,3808.5748,3809.4561,3809.835,3813.0678,3816.1666,3817.4767,3818.6855,3820.7926,3821.4308,3822.1485,3822.8619,3823.0676,3823.4577,3823.5852,3824.5533,3825.1331,3825.6729,3826.3688,3826.8072,3828.3846,3830.0606,3830.5100,3830.7736,3831.6398,3832.3035,3833.0860,3834.6787,3835.7111,3836.5851,3837.8752,3838.2922,3839.6953,3840.3809,3840.8004,3841.5187,3841.9601,3842.5502,3842.8968,3843.5087,3844.7311,3845.4055,3846.2491,3846.6267,3846.8876,3847.6199,3849.1832,3849.9114,3850.1340,3850.5813,3851.1583,3852.1353,3852.9590,3853.8277,3854.5108,3856.3544,3857.5089,3858.9038,3859.8396,3861.3524,3862.0522,3862.4218,3863.4059,3866.9092,3867.5786,3868.2547,3868.5284,3869.6633,3871.4563,3871.6174,3872.1371,3873.1479,3873.4740,3873.8224,3874.2438,3874.8619,3875.3731,3875.6462,3878.6620,3879.2688,3879.4454,3879.6441,3880.1948,3880.3332,3881.4980,3882.1430,3883.7683,3884.5245,3884.8225,3885.4038,3885.7681,3886.0049,3886.2564,3886.9159,3889.9059,3891.2045,3891.4017,3891.7260,3891.9792,3893.6519,3894.6600,3895.4192,3898.4374,3898.7967,3899.0310,3900.5767,3900.8785,3901.1522,3901.6621,3902.1708,3903.1024,3903.4812,3904.0828,3905.1865,3907.1603,3907.5432,3908.7491,3908.9797,3909.1390,3910.5210,3910.7734,3911.5760,3911.9091,3912.2810,3913.6450,3914.1623,3914.7675,3915.2953,3915.8489,3916.4176,3917.2693,3918.0701,3918.5115,3919.0234,3919.4602,3922.2191,3923.7995,3924.4034,3925.0934,3925.7188,3926.0456,3926.8630,3927.1762,3927.8062,3928.6233,3929.2904,3929.6693,3930.3303,3931.2566,3931.9941,3932.2264,3932.5466,3932.9113,3933.2375,3933.6611,3934.2742,3935.2037,3935.9388,3936.3532,3937.0406,3937.9237,3938.1523,3938.4227,3939.1593,3939.3188,3940.3686,3940.8387,3942.0732,3942.3580,3944.2534,3945.1352,3945.5072,3945.8201,3946.1455,$
3946.3921,3947.1364,3947.3308,3947.5046,3948.0305,3948.9789,3950.3951,3950.8050,3951.1085,3951.5149,3952.2355,3952.7608,3953.4448,3954.7267,3955.1699,3955.8903,3956.6908,3958.2372,3959.3000,3960.8796,3961.5208,3962.4196,3963.2199,3963.4688,3964.0302,3964.7257,3966.9647,3967.3921,3968.4673,3969.0026,3969.6641,3972.1545,3972.6395,3973.1961,3974.4766,3974.7590,3975.2218,3975.4674,3976.4149,3979.0433,3979.3559,3980.0896,3980.7550,3981.1065,3981.8281,3982.1020,3982.8956,3984.0955,3984.3757,3984.8795,3986.0739,3987.2061,3987.7062,3988.0273,3988.5997,3988.8456,3990.4922,3991.7309,3992.0535,3992.2726,3993.2989,3994.5494,3994.7918,3996.0617,3996.6687,3997.4687,3997.8654,3998.4038,3998.7334,3998.9526,4000.2811,4001.0581,4001.8935,4003.1057,4003.3086,4003.5740,4005.0928,4005.3628,4005.9615,4006.3806,4007.0188,4008.2102,4009.0112,4011.5916,4011.7398,4012.4952,4013.8566,4014.5143,4014.7164,4017.0627,4018.0990,4019.1289,4020.3541,4021.1498,4021.7506,4022.0674,4023.3378,4024.8025,4025.6556,4027.0091,4028.6991,4029.8256,4030.2925,4030.8424,4031.3581,4032.4531,4032.5951,4034.2461,4034.9218,4035.4600,4036.0479,4036.5652,4037.5614,4038.2287,4038.8043,4039.8648,4040.3931,4041.2036,4042.8937,4043.1302,4043.3948,4044.4179,4044.9269,4045.2268,4045.9654,4047.4815,4048.0318,4048.2876,4048.4324,4049.5361,4049.9447,4050.8872,4051.1090,4051.4991,4052.9208,4053.5277,4054.3017,4054.5258,4057.3312,4057.9412,4059.2529,4060.0616,4061.6254,4063.4071,4063.6920,4064.3315,4066.8219,4067.4507,4067.8759,4068.4703,4069.2014,4069.4612,4069.7608,4070.2383,4070.7835,4071.7513,4072.0047,4072.3849,4072.6284,4073.8563,4074.6420,4075.5030,4075.9070,4076.6284,4076.9432,4078.3015,4078.8751,4079.5738,4080.3580,4080.7058,4081.3678,4081.5914,4082.0817,4082.3872,4083.4688,4084.6192,4085.0421,4085.2564,4085.4341,4086.5205,4087.2848,4088.7264,4089.1379,4091.3476,4091.7368,4093.3921,4093.6713,4094.7470,4096.0759,4097.7478,4098.7316,4098.9327,4100.3414,4102.6177,4102.8199,4103.2627,4103.9121,4104.3821,4104.8385,4105.3301,4105.9150,4107.0506,4107.8614,4108.4198,4109.3234,4110.8260,4112.7545,4113.5575,4115.7589,4116.3753,4116.7137,4118.4893,4124.0403,4127.4120,4128.0499,4128.6400,4129.5116,4130.3247,4131.0021,4131.4265,4131.7235,4132.7533,4134.0681,4134.3235,4135.4800,4138.0406,4140.2350,4140.3855,4142.4741,4142.7010,4143.6491,4147.5396,4148.1816,4151.1453,4154.7205,4156.0860,4156.2382,4156.5166,4157.2804,4158.5905,4159.1537,4159.6599,4161.7389,4162.5090,4163.6469,4163.9479,4164.1795,4164.9735,4165.7661,4168.6341,4169.0851,4170.5333,4170.7856,4171.3410,4172.6257,4174.4396,4178.0597,4178.3658,4178.8481,4179.9601,4180.3272,4180.7210,4181.8836,4182.9342,4183.5645,4184.1376,4184.6045,4185.1447,4185.4693,4187.1410,4190.7129,4191.0294,4191.8269,4192.3620,4193.0164,4194.0826,4194.9361,4195.8326,4198.3170,4199.8891,4200.6745,4201.9715,4203.5762,4204.0413,4208.4108,4208.8907,4209.9727,4210.4554,4210.9232,4213.0673,4213.3779,4214.8285,4215.6272,4216.0698,4216.3722,4217.4308,4218.6649,4220.0651,4220.7323,4222.6373,$
4223.4403,4224.2413,4226.2994,4226.7267,4226.9876,4227.3872,4228.1580,4229.1477,4229.4546,4229.8644,4230.4267,4230.8242,4233.2873,4235.4636,4236.0471,4237.2198,4239.7873,4240.5947,4241.0948,4241.9725,4242.7197,4243.0195,4243.2611,4243.9241,4247.5988,4247.9886,4248.3908,4248.7400,4249.4761,4249.6791,4250.3146,4250.7749,4251.1846,4253.5385,4253.8674,4255.2373,4255.7972,4256.2537,4257.4963,4257.8982,4258.5204,4258.9955,4259.3619,4260.3330,4260.9845,4261.2754,4261.4920,4261.9665,4262.6124,4263.3561,4264.1064,4264.3399,4266.2864,4266.5271,4267.5191,4269.0639,4269.9426,4270.3290,4270.7349,4271.4396,4272.1689,4272.8745,4273.3574,4274.0246,4275.1596,4276.8073,4277.3139,4277.5282,4278.3232,4279.0627,4280.5680,4281.0678,4281.4145,4282.0413,4283.5184,4284.9749,4285.1823,4286.2285,4287.0808,4288.4694,4288.6688,4289.3632,4289.6551,4290.3937,4291.8098,4292.3057,4293.7698,4294.7193,4295.5849,4295.8155,4297.3066,4297.9392,4298.1224,4298.9860,4299.6352,4299.8393,4300.1008,4300.6495,4302.5261,4303.3867,4303.9892,4304.4179,4304.9565,4306.3668,4307.1762,4307.7406,4308.1220,4308.6001,4309.2392,4311.7994,4312.9975,4313.3084,4313.6025,4314.3195,4315.2543,4315.9484,4316.1094,4316.3726,4317.8400,4318.4157,4318.6491,4319.0973,4320.1264,4325.2741,4325.8599,4327.7145,4328.6879,4328.9154,4330.4120,4330.8438,4331.1995,4332.0297,4332.3395,4333.5612,4333.9369,4335.3379,4337.0708,4337.2774,4338.1078,4339.8716,4340.8954,4342.4442,4343.3817,4343.6035,4343.9515,4344.3265,4344.6092,4345.1680,4346.4367,4347.6385,4348.0640,4348.5982,4349.0722,4350.2717,4351.2719,4352.2049,4352.6122,4352.8201,4353.4487,4354.4824,4355.3205,4355.5259,4356.0445,4357.6132,4358.3200,4359.3719,4360.1672,4360.7178,4361.3071,4362.0662,4362.4721,4363.7945,4364.0401,4365.9301,4366.9657,4367.2862,4367.4180,4367.8316,4369.8756,4370.7532,4371.3290,4372.4900,4374.1239,4374.7851,4375.9542,4376.5309,4378.1768,4379.6668,4380.2863,4381.4018,4381.8601,4382.8934,4384.6560,4385.0566,4385.7535,4388.1014,4388.9680,4391.1105,4391.3893,4392.9740,4393.7590,4394.8945,4396.1394,4396.4783,4397.0097,4397.9149,4399.0939,4399.4060,4400.0968,4400.3863,4400.9863,4401.5812,4402.2452,4402.6001,4402.9271,4404.9022,4408.4822,4408.8828,4409.8991,4413.3745,4413.6324,4414.4863,4414.7670,4416.2370,4416.8447,4418.4157,4419.7692,4420.2571,4421.5439,4422.0479,4423.7202,4424.8373,4426.0011,4428.5536,4429.6348,4429.8227,4430.1890,4430.9963,4432.2523,4432.9628,4433.4907,4433.8380,4434.8266,4434.9564,4435.6787,4436.0484,4436.2844,4436.9961,4438.7464,4439.1237,4439.4614,4439.8793,4440.5737,4440.8659,4441.6084,4443.0863,4443.6657,4445.0338,4445.9012,4447.2296,4447.8346,4448.8792,4449.5206,4450.8040,4452.5655,4454.1392,4454.7738,4456.7080,4457.2362,4458.0015,4458.7384,4459.5477,4460.5574,4461.2412,4461.5278,4463.6660,4463.8346,4464.1386,4465.3406,4465.9751,4466.6245,4468.3208,4469.1891,4469.5253,4470.9906,4472.8462,4474.0730,4474.7594,4475.2214,4478.4032,4478.5958,4479.6376,4480.2661,4480.5453,4481.8107,4482.1693,4483.3468,4486.8973,4487.4955,4488.3123,$
4488.6808,4489.6645,4490.6715,4490.9816,4492.1007,4493.3337,4494.6907,4494.9616,4495.2378,4496.3155,4496.6381,4497.9144,4498.5384,4498.9401,4499.3796,4499.7037,4499.9832,4500.5769,4502.9268,4504.6306,4505.2167,4506.4730,4507.4185,4507.8339,4509.9587,4510.5259,4510.7332,4512.4832,4513.2231,4513.6799,4514.0573,4514.5503,4515.1182,4515.6274,4517.0352,4519.2592,4521.1939,4521.7407,4522.3230,4522.7839,4524.1289,4524.8379,4526.0286,4526.9177,4527.7086,4528.0208,4528.9740,4530.3191,4530.5523,4531.7134,4532.2577,4533.0770,4533.3042,4534.1198,4534.4025,4535.2546,4535.4903,4535.9798,4537.6426,4540.9990,4542.0513,4543.4153,4543.8692,4544.5142,4545.0519,4545.9156,4546.6757,4547.2498,4547.7589,4548.3261,4549.6262,4549.8341,4551.4736,4552.1536,4553.8376,4555.8127,4556.3872,4558.3458,4559.0571,4559.3120,4561.3477,4562.0803,4562.4756,4564.4054,4564.8354,4566.6497,4567.2403,4568.1426,4568.5235,4570.9722,4571.2195,4571.5312,4573.7038,4574.0279,4574.3058,4575.2338,4575.4259,4576.7405,4577.5822,4578.5512,4579.3495,4579.8272,4581.1732,4584.6518,4585.8642,4588.4261,4589.1222,4589.2668,4589.8978,4592.6660,4593.6437,4595.4206,4596.0967,4596.3078,4598.7627,4599.7046,4600.8717,4602.8862,4603.1446,4606.5021,4607.3780,4607.9347,4608.6201,4609.5673,4610.0809,4610.6272,4611.8597,4612.5435,4613.6044,4615.0240,4615.3340,4616.4523,4620.2408,4620.4476,4621.1629,4622.1209,4623.1169,4623.8911,4624.1365,4624.3140,4624.6776,4625.0538,4625.2747,4625.4938,4627.2980,4628.2015,4628.4409,4631.7617,4633.7657,4637.2328,4638.1259,4638.6849,4639.5188,4639.7036,4640.0462,4643.5566,4644.7072,4646.6861,4647.2509,4647.9564,4649.9777,4650.2343,4651.5552,4651.9895,4654.0057,4655.2126,4656.5518,4657.9012,4659.5698,4660.3941,4662.5936,4663.2026,4664.9704,4666.0050,4666.5158,4666.7985,4668.1716,4669.3561,4669.6658,4669.9842,4672.4372,4673.6609,4675.3760,4676.0555,4676.9721,4678.2352,4680.2377,4680.6460,4682.7322,4683.3517,4686.1946,4686.5715,4690.3343,4690.6219,4691.0516,4691.3452,4691.6354,4691.8842,4694.0914,4695.0381,4695.4542,4698.2248,4700.1402,4700.7711,4702.3161,4703.9898,4705.2986,4705.7606,4706.2511,4707.0450,4708.2940,4711.4184,4712.0056,4712.4814,4712.8408,4714.6715,4715.4308,4717.3905,4719.4424,4719.9795,4720.4586,4720.7811,4721.2765,4721.5910,4722.0886,4722.5006,4723.4382,4723.7840,4726.8683,4727.8491,4728.1330,4729.1282,4729.8795,4730.0918,4730.6593,4730.8815,4732.0532,4734.0458,4735.9058,4737.9172,4739.6764,4740.5292,4740.9585,4741.3039,4742.1174,4742.5671,4743.6870,4745.3375,4747.6176,4748.5911,4749.2002,4749.9713,4752.4141,4757.2196,4761.1101,4761.8614,4762.5242,4764.3463,4764.8646,4765.5950,4766.6006,4768.0578,4773.2410,4775.3130,4775.7940,4776.7792,4777.1915,4778.2940,4779.7286,4780.7505,4781.2902,4782.7612,4783.6389,4783.8617,4784.0396,4786.5310,4786.7250,4787.1480,4788.6781,4789.3868,4790.4373,4791.5157,4792.0819,4793.2446,4793.9052,4795.9131,4796.2429,4800.1724,4801.0517,4802.6714,4803.9559,4805.6063,4806.0205,4808.1337,4809.6140,4812.3755,4813.0069,4813.7204,4813.8963,4817.0206,$
4818.6477,4819.1930,4820.4649,4820.8847,4821.5878,4821.8646,4822.8548,4823.1823,4823.6058,4823.9967,4824.6484,4826.7004,4828.6596,4829.7973,4831.1213,4831.5975,4832.4233,4832.8025,4833.6759,4833.8277,4834.0970,4837.4121,4838.3546,4840.4744,4840.8492,4842.1676,4843.9413,4844.1653,4844.5597,4844.7553,4844.9251,4845.1626,4847.3263,4847.8095,4848.3625,4849.1399,4849.8617,4850.4397,4852.8685,4857.5385,4858.0944,4858.3327,4861.2167,4861.7173,4863.1724,4865.4775,4865.9105,4868.2742,4868.8814,4871.2890,4872.0309,4872.9169,4874.3645,4876.4949,4877.0015,4878.0094,4878.7330,4879.3497,4879.8635,4881.2046,4881.8531,4882.2432,4886.8661,4889.0422,4889.4903,4889.8553,4890.4582,4891.0379,4891.6606,4892.7592,4893.4453,4894.9551,4898.8044,4899.2401,4902.0545,4902.7943,4903.0583,4904.7516,4907.2093,4909.8431,4910.1576,4910.5488,4910.7929,4911.1358,4911.3787,4912.0428,4912.5293,4914.1210,4915.4155,4917.8220,4919.8157,4921.6134,4922.9442,4924.4223,4925.4254,4925.9501,4926.2529,4927.2988,4927.7803,4929.0860,4929.9850,4933.2091,4933.8521,4934.3306,4935.9283,4936.7746,4937.8295,4938.5047,4939.0605,4939.2707,4939.6422,4941.4160,4943.0642,4943.7076,4945.4587,4946.6637,4947.5752,4950.2513,4950.6264,4952.6908,4953.5320,4957.2963,4958.0988,4958.7240,4961.3814,4961.7260,4962.0621,4963.1881,4965.0795,4965.7315,4967.0628,4967.6529,4968.7556,4970.0785,4972.1597,4975.0917,4975.9485,4976.5938,4977.3906,4980.1859,4980.7561,4980.9513,4982.4875,4983.5324,4985.3725,4985.9486,4987.1470,4987.5582,4989.3086,4990.0324,4991.2038,4992.1257,4992.6372,4993.7488,4994.1058,4995.2981,4996.5333,4997.8189,4999.2386,4999.9398,5000.2463,5002.0972,5002.8933,5003.5981,5004.1279,5005.9752,5008.1897,5009.3344,5009.9367,5010.4174,5011.4774,5012.2754,5013.1647,5014.4480,5014.7539,5015.8893,5016.5356,5017.1628,5017.5084,5018.0535,5019.3241,5019.8062,5020.5458,5021.2531,5022.0051,5022.1740,5023.7084,5024.6122,5026.9549,5028.6556,5029.0119,5029.6295,5029.8916,5031.8147,5034.2944,5035.3375,5036.7284,5038.3018,5039.2303,5039.5265,5040.6805,5041.1224,5041.5998,5043.5145,5044.7195,5045.2481,5046.3527,5046.6372,5047.0434,5047.4283,5048.9366,5049.7960,5050.7842,5051.3427,5051.8887,5052.5853,5055.3473,5057.9866,5058.3611,5058.5600,5059.8611,5061.2225,5061.6562,5062.0371,5062.9325,5063.5157,5063.9988,5064.6020,5064.9454,5065.1926,5066.1355,5066.7773,5067.1379,5067.9737,5068.8293,5069.3384,5070.7780,5071.7534,5072.6288,5074.6465,5075.4659,5079.1374,5079.9079,5081.4462,5082.6224,5084.9935,5085.2955,5085.8392,5086.1774,5086.8351,5089.2192,5090.0513,5090.5455,5090.7510,5094.6897,5095.0639,5096.4848,5098.0432,5100.6211,5101.1299,5101.3452,5103.7653,5106.6689,5109.7331,5110.0473,5110.5528,5111.0611,5111.2781,5114.0665,5115.0448,5116.2960,5117.2923,5118.2023,5119.8676,5121.3246,5122.4995,5124.2460,5125.4895,5125.9502,5128.4897,5129.5520,5130.2338,5131.0720,5131.9924,5133.1051,5134.7460,5136.1210,5137.4733,5139.4487,5140.7736,5141.2650,5141.7827,5143.2673,5143.9165,5145.0383,5145.3083,5146.0557,5147.1346,5148.2115,5149.2069,$
5149.8080,5150.7000,5151.6120,5153.5173,5154.2430,5154.8985,5157.4257,5158.6041,5159.4534,5159.6202,5161.5396,5162.2846,5163.4584,5164.4014,5165.7728,5166.6553,5168.5866,5168.9225,5170.2227,5172.4786,5172.6890,5172.9175,5173.6715,5174.3699,5174.7988,5175.3248,5175.9115,5176.2292,5176.4037,5176.9610,5177.6228,5178.4800,5179.1440,5180.7209,5182.5269,5183.0897,5183.3408,5183.6044,5183.9896,5184.4534,5184.7315,5185.2042,5186.4132,5187.1677,5187.3374,5187.7462,5188.3656,5188.8432,5189.6772,5190.8720,5193.8256,5194.4576,5195.8136,5197.2360,5198.7999,5199.1637,5202.0087,5203.8479,5205.1522,5205.7787,5206.4956,5206.6625,5206.7995,5207.8015,5209.7246,5211.2305,5213.3492,5214.1613,5216.5966,5216.8139,5218.5271,5219.1099,5220.3085,5220.7068,5220.9263,5221.2710,5222.2629,5222.7661,5224.1226,5226.5254,5228.2246,5228.9951,5229.7945,5231.1597,5232.3415,5233.2254,5234.1069,5237.6280,5237.9033,5238.8137,5239.5520,5240.1968,5241.8544,5242.0892,5245.3781,5247.1968,5247.6547,5250.8727,5252.5681,5253.4441,5254.2589,5255.5736,5255.9188,5258.3602,5259.7760,5260.1041,5261.4720,5261.6995,5262.2369,5264.2332,5264.7824,5265.5514,5266.7103,5269.7927,5270.2660,5272.6439,5272.9270,5273.1317,5274.1188,5274.3897,5276.4089,5277.1460,5277.5002,5279.1968,5280.0855,5280.3441,5281.0687,5281.6285,5282.3948,5283.6959,5284.5421,5286.8870,5289.6226,5289.8974,5291.8163,5294.3971,5295.0886,5296.2787,5297.4075,5297.7431,5298.2824,5298.5601,5300.5234,5301.4043,5303.4833,5306.9864,5307.4655,5308.3105,5309.6093,5310.2665,5310.8662,5312.0018,5312.5288,5312.9045,5315.2278,5317.4944,5318.1104,5318.4268,5320.0985,5320.7699,5322.8988,5323.4002,5325.1438,5325.4317,5326.2774,5326.9756,5329.3745,5329.7559,5330.0804,5331.7477,5333.3846,5337.0154,5338.6405,5339.8952,5340.4982,5343.5812,5345.3117,5346.3778,5347.9714,5348.3092,5349.0054,5349.4608,5350.5361,5351.1265,5351.8365,5354.6016,5355.6365,5358.7059,5359.8269,5360.1501,5361.1556,5361.4101,5362.2506,5362.5751,5363.6073,5365.4705,5366.7066,5369.2819,5369.4470,5370.7096,5371.1238,5372.7027,5374.8219,5375.3526,5375.7688,5376.1304,5376.7809,5378.8356,5379.1105,5381.3709,5382.9276,5384.0356,5386.6107,5387.8131,5388.0507,5388.9061,5391.0760,5392.5726,5393.5995,5393.9719,5394.7608,5396.1127,5397.5160,5398.2047,5398.7015,5398.9219,5399.1745,5399.6218,5400.1455,5402.6048,5403.1993,5405.2097,5407.3439,5407.6535,5410.7687,5412.2212,5414.5246,5415.5226,5417.4858,5419.1285,5421.3517,5421.6628,5421.8359,5422.9151,5424.0079,5424.7296,5425.6783,5426.4075,5427.6207,5429.1048,5430.5977,5431.1120,5431.6053,5432.7650,5433.2919,5433.7005,5434.1514,5435.1229,5435.7236,5435.8925,5436.4980,5437.3877,5438.1374,5439.8518,5439.9891,5440.6015,5441.0440,5441.2145,5443.1187,5447.1536,5448.2719,5449.4788,5449.8974,5450.5104,5451.6520,5452.2187,5455.1524,5458.9676,5460.7504,5461.7358,5462.3372,5462.6129,5463.3214,5463.7697,5464.2055,5467.1704,5470.1460,5470.7591,5476.4856,5477.0540,5477.5522,5478.7446,5479.0751,5480.8907,5481.4813,5481.9629,5483.1445,5484.1467,5484.5596,5487.8414,$
5488.6287,5488.9254,5489.0885,5490.1194,5492.6435,5493.2040,5494.3306,5495.8738,5496.1370,5498.1841,5499.2554,5499.6478,5501.2814,5501.9439,5504.3018,5506.1128,5506.9196,5507.5385,5508.5585,5509.9938,5510.6832,5511.0608,5512.7002,5512.9793,5514.3760,5514.8731,5517.7952,5518.9893,5519.2353,5521.7532,5523.1675,5524.2188,5524.5824,5524.9570,5527.2953,5528.2272,5529.0972,5530.0763,5530.6977,5535.9710,5537.1307,5537.5563,5538.6087,5539.2618,5539.9106,5541.1466,5541.5815,5541.9357,5542.8900,5545.0495,5545.8058,5546.1204,5547.1339,5548.1758,5551.3719,5552.6229,5555.5312,5556.1487,5557.0454,5557.9214,5558.3422,5558.7020,5559.8912,5564.2011,5567.9986,5568.5397,5570.9268,5571.1913,5572.0951,5572.4649,5573.3535,5573.8708,5574.9153,5576.2045,5577.6845,5579.3583,5580.0771,5580.3924,5580.7547,5581.9668,5582.3666,5582.6789,5582.8665,5583.7621,5585.8472,5587.0263,5587.7351,5588.7506,5590.1127,5593.2683,5593.6129,5594.4604,5595.0635,5595.6097,5597.4756,5598.4794,5599.6548,5601.1216,5601.6032,5602.6261,5604.5154,5605.2957,5606.3861,5606.7330,5609.5735,5610.6809,5612.0679,5612.6174,5615.3195,5615.7286,5616.6908,5617.2897,5617.6417,5619.9755,5620.9213,5621.7884,5624.9129,5625.6782,5626.0599,5626.7369,5629.3098,5630.2969,5633.2950,5635.8792,5639.7461,5640.3623,5641.5575,5641.7342,5643.1086,5645.5249,5645.6688,5645.8947,5646.4514,5647.7069,5648.6863,5648.9910,5650.7043,5652.9020,5654.0234,5654.5185,5655.4897,5656.7273,5657.1804,5657.9255,5659.1272,5660.6594,5663.0420,5664.6210,5665.1799,5665.6283,5666.4185,5667.1281,5671.0711,5673.8358,5674.9864,5675.5016,5677.0529,5679.0050,5681.9001,5685.1921,5687.3489,5690.6934,5691.6612,5696.3904,5698.2937,5700.4584,5700.9176,5702.6511,5703.9411,5705.6371,5707.1033,5708.6794,5711.9957,5715.5251,5715.7244,5717.1711,5717.5246,5719.0971,5719.6227,5720.1828,5721.4243,5724.2374,5724.4638,5725.0114,5725.3885,5725.8959,5727.0131,5727.7103,5729.9823,5731.6022,5732.9750,5733.8851,5736.0297,5738.9669,5739.5196,5741.1705,5741.8290,5742.0838,5745.6755,5747.1931,5747.3849,5748.7412,5749.3883,5749.7862,5753.0265,5760.5508,5761.5272,5762.7941,5763.5290,5767.7785,5768.1812,5771.7601,5772.1143,5773.9463,5776.1678,5777.4006,5779.8352,5781.6588,5782.3139,5784.8625,5789.6451,5791.7036,5792.4304,5796.0683,5797.3194,5798.4780,5800.8297,5802.0798,5804.1412,5805.7017,5807.6814,5808.6569,5811.8445,5812.9725,5815.4219,5817.7348,5818.4495,5819.1272,5819.6027,5819.9000,5820.3799,5822.5217,5822.7933,5829.1097,5829.8625,5830.8272,5832.3705,5834.2633,5838.9502,5839.8514,5840.6404,5842.0500,5843.8069,5844.7909,5845.9188,5852.6806,5853.4745,5854.1207,5855.4986,5857.4491,5858.7543,5859.6687,5860.3103,5863.7184,5866.8119,5868.3745,5869.8507,5870.5526,5871.1828,5872.6028,5872.8842,5874.0383,5874.3513,5874.9865,5878.2764,5879.1265,5882.6242,5883.6247,5884.0331,5885.7016,5886.5314,5888.2622,5888.5841,5889.9531,5891.4510,5892.4328,5893.4747,5894.6979,5895.2818,5895.9281,5896.3631,5896.7885,5899.8441,5902.6021,5904.1592,5905.1763,5905.5707,5911.2296,5912.0853,5913.3638,$
5914.3856,5914.6709,5916.7282,5918.9442,5925.4036,5926.2321,5928.8130,5929.4802,5929.9343,5932.4274,5936.3861,5936.8740,5937.1619,5937.6634,5937.9252,5938.4570,5938.8252,5940.5653,5942.3585,5943.4653,5944.6472,5946.2347,5948.7991,5950.4265,5953.5825,5954.5852,5954.8839,5955.5616,5956.2588,5957.5885,5960.7806,5967.5099,5968.3199,5969.7373,5970.4506,5971.6008,5972.4351,5973.6649,5975.0648,5976.2005,5981.8846,5982.0960,5986.2665,5987.3016,5989.0447,5989.7430,5991.0071,5993.4941,5994.1287,5995.2198,5996.6297,5998.9987,6000.7622,6001.2033,6001.7339,6005.1650,6005.7242,6007.0722,6010.1606,6010.5994,6011.5337,6013.6777,6014.0563,6015.4218,6016.3589,6017.8618,6018.9934,6019.6385,6020.2715,6021.0357,6021.4113,6023.2242,6025.1500,6025.7779,6026.0448,6028.2713,6029.2257,6029.6500,6030.4450,6032.1274,6032.8726,6033.4131,6035.1923,6036.7622,6037.6975,6038.6805,6042.5898,6043.2233,6044.4327,6046.8977,6049.0510,6049.7721,6050.9818,6052.7229,6053.3808,6055.5937,6058.1813,6058.6123,6059.3725,6061.5361,6064.7508,6065.7799,6066.4973,6067.0803,6069.0207,6070.3412,6073.1036,6077.1057,6077.8728,6078.4212,6079.2227,6079.4724,6085.3748,6086.7022,6087.2622,6088.0305,6090.4482,6090.7848,6097.1945,6098.1205,6098.8031,6099.0831,6099.9886,6101.1615,6101.7251,6105.6351,6106.8431,6107.5337,6112.8375,6113.4657,6114.5392,6114.9234,6116.1661,6119.6998,6120.5564,6121.4075,6122.2144,6123.3619,6123.8331,6124.0687,6124.4805,6125.7396,6127.4160,6129.5458,6131.9636,6132.8632,6133.0571,6133.8143,6137.9260,6138.6465,6144.7581,6145.2833,6145.4411,6150.6830,6151.5133,6151.9929,6154.0682,6154.5162,6155.2385,6155.5810,6157.0878,6161.3534,6162.1700,6162.9699,6164.4796,6165.1232,6166.4311,6169.0326,6169.8221,6170.1740,6172.2778,6173.0964,6178.4315,6179.8488,6180.7050,6182.6217,6187.1350,6188.1251,6189.1449,6191.9053,6193.8563,6195.3210,6196.3965,6198.2227,6200.4325,6201.1002,6203.4925,6204.1281,6205.0180,6205.2561,6205.8599,6207.2201,6207.7504,6208.6869,6208.9861,6209.9496,6212.5031,6212.7191,6215.9383,6217.5785,6219.7779,6220.0112,6221.0599,6221.3192,6223.5202,6224.5272,6226.3697,6231.1227,6231.3808,6232.9745,6234.1962,6234.8554,6239.7121,6240.9536,6242.9407,6243.1201,6243.3603,6245.0409,6245.6182,6246.1736,6246.9297,6248.4055,6250.4858,6257.4237,6258.6068,6261.0642,6261.4181,6264.7139,6265.6036,6266.1737,6268.2001,6271.5444,6274.1166,6276.1646,6277.2385,6279.1666,6279.9776,6285.2780,6286.0453,6287.2554,6289.4882,6291.1915,6291.6286,6292.0601,6292.8909,6293.2424,6296.8722,6298.9027,6300.9165,6301.4128,6303.2507,6304.2423,6307.6570,6309.1600,6310.8101,6312.6222,6314.2693,6315.7750,6317.1824,6321.3004,6321.8199,6323.8421,6324.3978,6326.3669,6327.2778,6331.4137,6332.5118,6333.1459,6335.4060,6337.6201,6339.3061,6339.6684,6342.8595,6346.1209,6348.2321,6348.7375,6355.6300,6355.9108,6357.0229,6357.6779,6358.6211,6359.1345,6359.6744,6362.2518,6364.8937,6369.1394,6369.5748,6371.9436,6372.4602,6374.8255,6376.9305,6378.7604,6379.6733,6381.3598,6381.7588,6384.7169,6387.3957,6388.4179,6389.3900,6390.1289,$
6391.1774,6392.3685,6393.7972,6394.0498,6396.9499,6397.4536,6399.2065,6399.9498,6400.6961,6403.0128,6406.4462,6408.9035,6411.8991,6413.6145,6414.7040,6415.5374,6416.3071,6418.3703,6422.1075,6422.8969,6424.8127,6427.5044,6428.7733,6431.5550,6436.6699,6437.7613,6438.3095,6438.9160,6439.0715,6441.8994,6443.8598,6445.6407,6446.7712,6449.8081,6450.0060,6450.9557,6452.0591,6455.2647,6455.5964,6457.2824,6460.4750,6461.0596,6466.3572,6466.5526,6466.7171,6467.6754,6471.2131,6471.6577,6481.3017,6483.0825,6484.0694,6485.3759,6487.4811,6488.8832,6490.7372,6493.1975,6493.7777,6495.2573,6496.0292,6497.4919,6499.6456,6500.6575,6501.9919,6503.5111,6506.9868,6508.3605,6509.0503,6510.9969,6512.3639,6522.0432,6522.4994,6526.9045,6531.3418,6534.6056,6536.1441,6537.1775,6537.6139,6538.1120,6542.0498,6542.5131,6545.7186,6550.1877,6550.5475,6551.7055,6554.1603,6558.8756,6560.0574,6564.4445,6565.0700,6569.6322,6576.1221,6577.2146,6577.6552,6580.2299,6583.9060,6584.6130,6585.6945,6588.5396,6591.4845,6593.4628,6593.9391,6596.1141,6596.9623,6598.6778,6599.4824,6600.7311,6602.7620,6603.6209,6604.8534,6605.4165,6608.4430,6613.3740,6617.0580,6617.5154,6618.1668,6619.9458,6620.9665,6632.0837,6633.4579,6638.2207,6638.9119,6639.7403,6643.6976,6644.6635,6646.5407,6648.4954,6648.9586,6654.3675,6655.4892,6656.9386,6658.6774,6660.6761,6662.2686,6663.6966,6664.0510,6666.3588,6668.8163,6673.5797,6674.6969,6677.2817,6678.7067,6680.0811,6683.3673,6684.2929,6687.5207,6692.7262,6694.0068,6694.4967,6696.1400,6697.7123,6700.7499,6703.8966,6704.0510,6711.2521,6713.9700,6715.1885,6717.3851,6719.1995,6726.3091,6727.4583,6728.1183,6728.7595,6729.9328,6732.6491,6733.7488,6738.1800,6741.9212,6742.8845,6751.4267,6752.8335,6753.6598,6756.4528,6757.1092,6758.2035,6760.8923,6766.6117,6770.1069,6772.1745,6773.0973,6778.3123,6779.3242,6780.1252,6780.4131,6787.7364,6788.8403,6791.2351,6795.7982,6798.4868,6800.4645,6807.3190,6808.5311,6809.1000,6809.3173,6809.5094,6809.6773,6810.5494,6812.4915,6812.7752,6815.6123,6823.5084,6824.6774,6827.2488,6829.0355,6829.3472,6832.8927,6834.9246,6839.2949,6846.4790,6851.8842,6852.3545,6853.5234,6854.1093,6854.5128,6855.3150,6855.6891,6861.2688,6862.8685,6863.5350,6866.3666,6866.7634,6868.4508,6870.7261,6871.2891,6874.7530,6879.5824,6882.8114,6883.3147,6886.4082,6887.0881,6888.1742,6888.8347,6889.3032,6892.2502,6902.2105,6908.9905,6909.8491,6911.2264,6914.7127,6916.1289,6920.0394,6924.6611,6925.0094,6936.6528,6937.6642,6942.5380,6943.6105,6945.4902,6946.2134,6947.9170,6948.2052,6948.3942,6951.4776,6952.9666,6954.6560,6955.3150,6960.2500,6962.3117,6965.4307,6969.2976,6981.0831,6985.4718,6986.0298,6989.6553,6992.2126,6992.6966,6993.0371,6993.9873,6996.7573,6999.2706,6999.6238,7000.8036,7002.8829,7007.0961,7014.9689,7015.3172,7018.5675,7020.4841,7021.2832,7025.2251,7026.4616,7028.0228,7030.2514,7032.9310,7033.3588,7033.6523,7036.2831,7038.7202,7045.7970,7053.6196,7054.4177,7055.9011,7058.4895,7059.5253,7060.0417,7060.6538,7061.3942,7064.4515,7067.2181,7068.7358,7071.0942,$
7071.4800,7072.3939,7074.2563,7075.3336,7084.1690,7086.7044,7088.8232,7089.3395,7100.5144,7102.5910,7107.4778,7109.1819,7109.8602,7112.9189,7113.9879,7114.3984,7122.0441,7124.5607,7125.5214,7125.8200,7130.1846,7130.7236,7131.3586,7132.6100,7140.4617,7142.3310,7147.0416,7148.1451,7148.5589,7150.2844,7153.5877,7154.7622,7154.9538,7156.9416,7158.8387,7159.9470,7162.5569,7167.2028,7168.8952,7170.3607,7173.3725,7176.1860,7176.7213,7179.7147,7188.5320,7191.1328,7200.0454,7201.8091,7202.1910,7202.5166,7206.4830,7206.9804,7208.0063,7212.6896,7218.0542,7219.1515,7220.9809,7225.1099,7229.9386,7230.8623,7233.5365,7240.1848,7242.0919,7244.6965,7246.1277,7248.9668,7250.5895,7253.6760,7255.3541,7256.9860,7258.1770,7265.1724,7270.6636,7272.9359,7284.9033,7285.4437,7296.2659,7298.1434,7305.4043,7308.6418,7311.7159,7315.0661,7315.6112,7316.0050,7323.2106,7324.8073,7326.1491,7328.2850,7329.4916,7335.5772,7337.7812,7339.6035,7341.1515,7342.2995,7342.5769,7345.3996,7346.3428,7348.0530,7350.8140,7353.2930,7358.3449,7361.3472,7362.1942,7370.8255,7372.1184,7376.8772,7380.4263,7383.9805,7385.5006,7386.3445,7392.9801,7393.4377,7402.2521,7403.9901,7404.2595,7410.9683,7411.7363,7412.3368,7413.5958,7417.7908,7418.5499,7422.3118,7423.8044,7425.2942,7428.9405,7430.2533,7435.3683,7436.2970,7440.4933,7443.8753,7444.7489,7447.8488,7449.6049,7449.8091,7455.2080,7458.7541,7461.8746,7462.9910,7465.4400,7469.0566,7471.1641,7481.3545,7483.6256,7484.3267,7487.9738,7493.4276,7495.5641,7499.0024,7500.6556,7503.8691,7508.4787,7510.4082,7511.3498,7511.7902,7514.6518,7518.7824,7523.1347,7525.5079,7528.4892,7531.1437,7536.4365,7537.4287,7549.3138,7555.3254,7557.7495,7565.8515,7566.5296,7567.7417,7569.5115,7570.5590,7571.0326,7573.3437,7580.3460,7585.7922,7586.2132,7586.9943,7589.3151,7598.2054,7603.6243,7607.8230,7616.6827,7618.3443,7620.0773,7623.5683,7624.3123,7625.1195,7625.7053,7627.1749,7628.8818,7630.3106,7632.1305,7635.1060,7636.1746,7637.3852,7638.7805,7644.7403,7647.3794,7651.0006,7651.7451,7652.3202,7653.8284,7654.6998,7658.3246,7660.0234,7660.8903,7666.5679,7668.9608,7670.0575,7672.2552,7676.2195,7678.1267,7683.0189,7685.3075,7686.1538,7691.7706,7693.8016,7697.9239,7699.4034,7699.7794,7701.1081,7703.6850,7704.8169,7709.5523,7709.8915,7710.2692,7711.1248,7712.4050,7713.9378,7721.2023,7723.7611,7724.2072,7728.9510,7731.7385,7742.5628,7743.9552,7761.7133,7762.7320,7771.9468,7776.6727,7778.4922,7782.3165,7787.8002,7788.9342,7793.1326,7798.3579,7804.1414,7804.5300,7807.8759,7808.3304,7810.6245,7813.4761,7813.9729,7814.3230,7816.1534,7817.7669,7822.3877,7834.4578,7835.6209,7836.4597,7840.2981,7840.4576,7841.7911,7842.2656,7847.5394,7848.4547,7849.6254,7861.9100,7864.0218,7865.9698,7867.6883,7868.1946,7872.6316,7875.4626,7877.0756,7886.2830,7891.0750,7893.6191,7896.4528,7899.6229,7900.3200,7904.4344,7916.4420,7920.8570,7922.5653,7924.9914,7925.7478,7937.7335,7941.7259,7945.0699,7948.1764,7954.5923,7956.9740,7961.9832,7970.2537,7972.5960,7974.1592,7978.9731,7981.2262,7987.9731,7991.3655,$
7992.1588,7993.6808,7996.9691,8000.0449,8003.4551,8006.1567,8014.7857,8017.5275,8022.2014,8024.2530,8025.7271,8026.2135,8030.2004,8032.4313,8032.9061,8035.6377,8037.2183,8046.1169,8050.6480,8053.3085,8054.5355,8062.6304,8066.8267,8075.6518,8085.2190,8089.4852,8092.2373,8093.6238,8094.0559,8096.2600,8103.6931,8105.3947,8115.3110,8119.1811,8122.7234,8129.4051,8137.9363,8138.4753,8139.9032,8143.1380,8149.7021,8152.3818,8157.5434,8159.7277,8162.0587,8163.1210,8166.4477,8169.7865,8172.1234,8177.1788,8186.9113,8190.8849,8194.3943,8198.4420,8202.1469,8203.2009,8205.1110,8207.4785,8209.4461,8214.1472,8217.2264,8227.4722,8227.6828,8231.4069,8234.9381,8252.3936,8253.6156,8254.7422,8259.5110,8261.0143,8264.5225,8275.6266,8288.4143,8292.5272,8295.5497,8297.1765,8304.4244,8311.6306,8320.8554,8330.4494,8335.7067,8341.4770,8346.5427,8356.0693,8357.8357,8358.7260,8360.4915,8366.0741,8367.3936,8369.3404,8379.2261,8379.7679,8383.0884,8384.7240,8385.7280,8387.1053,8388.5363,8398.1780,8399.2578,8401.9890,8403.7965,8408.2096,8411.9172,8416.7269,8417.9982,8421.2254,8424.6475,8445.4870,8446.5116,8450.6754,8456.3468,8464.2367,8465.6706,8471.8260,8478.3580,8490.3065,8498.0161,8500.6797,8501.4398,8510.6240,8511.9091,8513.4063,8516.5542,8519.3316,8521.4422,8530.9108,8531.4506,8532.9137,8534.6797,8539.7930,8542.0871,8543.7225,8544.5955,8554.9440,8556.3250,8558.4464,8560.4348,8568.2088,8573.1205,8575.3305,8577.2768,8587.6342,8591.8387,8593.1089,8599.7553,8604.0163,8605.7762,8616.2219,8620.4602,8621.3225,8629.1415,8631.3565,8638.3623,8639.4416,8645.3087,8649.1490,8655.8760,8662.1372,8665.4855,8667.9442,8675.3983,8678.4083,8687.8480,8701.1209,8703.7025,8704.8601,8707.3592,8709.2341,8710.4142,8712.8528,8713.6547,8719.6290,8721.6595,8722.4577,8724.3761,8732.4240,8734.0234,8739.7816,8748.0309,8749.1697,8758.2434,8760.4496,8761.6862,8766.7450,8771.8602,8772.8052,8775.5733,8782.7156,8784.5621,8790.3761,8792.0572,8798.1720,8799.0875,8804.5894,8810.2540,8812.5111,8816.1728,8817.7431,8820.2312,8820.4210,8821.7586,8829.6938,8841.1833,8842.0744,8848.3052,8849.3151,8852.7916,8854.9079,8857.8731,8859.0181,8859.4123,8860.9757,8862.2915,8866.7145,8868.8334,8875.2324,8881.9003,8882.5115,8889.1939,8890.2805,8892.9865,8893.5397,8899.2971,8905.6581,8907.0331,8910.8566,8912.7739,8917.5099,8919.7742,8920.2012,8924.2459,8927.7293,8928.0925,8941.6608,8948.4211,8949.1227,8955.8467,8957.9860,8959.2846,8962.1468,8962.8942,8964.0530,8967.0031,8967.6403,8968.9461,8969.8667,8979.7026,8980.7395,8985.2808,8987.4079,8990.4779,8990.8935,8995.1893,8997.8763,8999.5264,9008.4636,9009.8832,9011.5152,9012.5263,9013.9748,9016.5903,9017.5912,9031.8194,9035.9195,9037.8938,9038.6911,9040.1229,9045.3532,9048.2501,9053.4858,9056.0813,9058.2586,9062.5631,9063.9600,9066.1115,9068.0230,9069.5817,9072.2785,9075.3945,9090.8186,9094.8289,9101.0826,9103.3420,9105.4689,9107.2269,9109.3144,9111.5314,9113.4542,9118.1378,9119.6365,9122.9674,9126.3292,9129.1822,9130.9682,9132.2739,9134.6919,9135.9375,9137.7316,9140.5559,9153.3664,$
9155.2968,9156.9236,9157.4218,9159.0346,9165.8950,9167.7950,9170.8220,9178.7795,9181.8967,9187.5654,9192.5937,9193.5928,9194.6385,9197.2557,9199.6841,9203.9617,9208.0252,9208.5811,9215.9197,9221.4324,9224.4992,9227.5119,9232.4960,9233.2739,9233.8574,9234.3987,9239.3261,9240.2158,9243.7603,9245.2557,9248.1186,9249.9067,9250.5782,9260.3254,9263.6830,9266.2070,9266.9198,9267.6895,9270.1501,9271.1802,9276.2732,9279.7085,9286.3805,9289.5624,9291.5313,9294.9740,9300.0131,9307.8956,9310.4440,9317.7296,9323.0400,9327.2462,9336.1624,9340.7053,9349.2523,9354.2198,9355.9911,9360.9879,9366.7971,9379.7288,9380.6439,9383.2722,9384.0997,9388.9308,9390.5852,9392.2645,9393.7757,9399.0891,9406.8923,9409.3488,9413.6761,9414.0886,9417.4570,9420.4985,9422.3137,9422.9331,9431.5996,9432.2827,9435.1230,9436.8127,9446.9891,9450.4608,9455.2023,9456.0216,9461.0279,9467.1954,9470.6819,9474.8793,9486.9256,9495.4979,9497.1891,9500.2999,9504.2159,9505.3930,9507.6525,9508.4513,9510.9465,9514.7986,9535.6566,9536.4070,9548.0303,9553.9847,9561.2452,9565.5588,9567.2805,9567.8292,9570.4022,9571.5019,9577.3475,9582.8133,9587.0275,9588.8066,9590.3432,9595.3906,9605.8059,9607.5347,9608.4866,9608.9352,9613.6911,9619.2179,9620.9953,9625.1984,9627.6706,9629.5693,9630.7442,9632.6439,9636.9019,9641.1951,9642.4776,9643.3195,9657.7863,9663.6433,9664.6983,9674.7914,9676.8364,9685.6771,9695.0311,9700.5631,9702.2720,9713.1103,9716.1420,9718.4922,9736.2125,9738.6224,9739.7704,9743.5614,9746.4629,9753.5920,9757.2197,9764.6062,9769.5348,9779.4532,9784.5028,9789.5111,9796.2002,9797.2456,9800.3629,9801.7074,9812.6976,9814.9578,9819.1781,9826.4497,9833.4232,9837.1669,9838.0046,9840.9172,9845.6835,9855.7421,9864.5972,9865.4489,9867.8902,9868.9206,9871.9947,9872.6300,9873.8151,9878.5208,9880.7853,9896.0483,9898.3532,9902.3569,9906.3907,9907.4718,9910.0750,9911.1129,9912.1996,9913.6234,9923.3088,9927.3252,9932.7748,9934.7173,9935.1992,9938.8355,9943.0647,9952.3704,9952.8033,9963.4911,9970.4625,9974.6893,9985.0500,9987.6352,9989.9371,9992.6485,9993.8630,9998.5092,9998.9599,10011.3953,10022.2836,10033.2218,10037.1188,10039.3641,10045.3125,10048.0365,10054.9620,10056.2091,10079.5414,10082.8793,10083.7854,10086.4035,10089.1351,10102.5746,10105.0769,10117.9901,10133.5662,10137.3865,10140.4309,10141.3970,10144.2649,10175.0079,10178.5192,10180.5933,10184.5392,10188.4156,10211.5330,10218.4327,10223.6585,10236.0280,10241.7747,10247.5568,10250.6792,10255.5804,10257.3690,10258.1668,10271.1881,10273.6965,10276.7976,10283.1162,10288.9919,10293.0491,10301.1668,10308.5442,10316.8910,10343.8210,10346.5367,10349.0483,10358.1678,10412.9915,10419.5787,10429.6706,10450.4421,10459.7208,10470.0535,10492.2578,10494.8400,10498.4920,10512.2045,10527.7862,10556.4500,10565.3047]

;depends on which spectra mode we are in
case spmd of
  '1f': step=1
  '2f': step=2
endcase

nord=(max(line_list[0,*])+1)*step ;number of orders extracted
res=fltarr(5,nord)                ;variable with the final wavelength solution for all the orders
vec_lines_used=[-1,-1]            ;vector of the lines that were used for the final solution
for i=0,nord-1 do begin
  for ll=1,2 do begin
    ;takes the information for the current order
    pos=where(line_list[0,*] eq i/step)
    ;vector with the pixel position of the lines found
    pix_pos=line_list[1,pos]
    ;if we are in the first loop, the guess solution comes from coeff_orders. if not, it comes from the result from the first loop
    if ll eq 1 then coeff=coeff_orders[*,i/step]
    
    ;removes the lines from the vignetted area (last 300 pixels)
    pos=where(pix_pos gt 4300)
    if max(pos) ne -1 then remove,pos,pix_pos
    
    ;applies the current solution to convert pixel position into wavelength
    line_wave=coeff[0]+coeff[1]*(pix_pos+1)+coeff[2]*(pix_pos+1)^2+coeff[3]*(pix_pos+1)^3+coeff[4]*(pix_pos+1)^4
    ;compares the calculated wavelengths with the ThAr line list
    pos=where(ThAr_list ge min(line_wave) and ThAr_list le max(line_wave))
    ThAr_list_temp=Thar_list[pos]
    ;measures the difference between the calculated and the expected wavelengths
    vec_diff=fltarr(n_elements(line_wave))
    for k=1,2 do begin
      for j=0,n_elements(line_wave)-1 do begin
        junk=min(abs(line_wave[j]-ThAr_list_temp),pos)
        vec_diff[j]=line_wave[j]-ThAr_list_temp[pos]
        if k eq 2 then line_wave[j]=ThAr_list_temp[pos] ;if we are in the final (2nd) loop, it associates the line with its wavelength
      endfor
      if k eq 1 then line_wave=line_wave-median(vec_diff) ;if we are in the first loop, if shifts the calculated wavelength to get closer to the real ones
    endfor
    
    ;looks for difference in line central wavelength greater than 5 km/s
    ;  (5 x the expected stability; Barrick et al. 2016, SPIE)
    pos=where(abs(vec_diff) ge .1) 
    if max(pos) ne -1 then if n_elements(pos) lt n_elements(vec_diff)-3 then remove,pos,pix_pos,line_wave else if ll eq 2 then print,'WARNING!!! wavelength solution might be wrong for order '+strtrim(string(i/step+22),2)
    
    ;Does the final fit of the pixel position vs wavelength for all the remaining lines
    ;the while loop is to remove deviant points
    pos=0
    while max(pos) ne -1 and n_elements(pix_pos) ge 3 do begin
      new_coeff=poly_fit(pix_pos,line_wave,4,yfit=yfit)
      dev=line_wave-yfit
      pos=where(abs(dev) ge 3*stddev(dev)) ;removes deviant points
      if max(pos) ne -1 then if n_elements(pos) le n_elements(line_wave)-3 then remove,pos,pix_pos,line_wave else pos=-1
    endwhile
    ;either feeds with a better initial solution (first loop) or send the final solution (second loop)
    if ll eq 2 then begin
      res[*,i]=new_coeff 
      for k=0,n_elements(pix_pos)-1 do vec_lines_used=[[vec_lines_used],[pix_pos[k],line_wave[k]]]
      vec_lines_used=[[vec_lines_used],[-1,-1]]
    endif else coeff=new_coeff
  endfor
endfor

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


