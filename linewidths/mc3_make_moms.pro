pro mc3_make_moms,topdir=topdir,outdir=outdir,data=data,line=line,cloud=cloud,method=method $
                  ,name=name,force_redo=force_redo,rms=rms,emchans=emchans

;& mask tuning parameters
  
;&;&;&;&:&;&;&;&;&;&:&;&
;& Defaults
;&;&;&;&:&;&;&;&;&;&:&;&

  use_data='match_3p5'
  use_tag='3p5as'
  use_line='12CO'
  use_method='DIL'          ; 'BOSMA' | 'DIR' | 'SMO'
  use_cloud='GMC1'
  do_force_redo=0 ; don't regenerate map if it exists
  use_topdir='../orig_data/'
  use_outdir='../moment_maps/'
  use_name=use_cloud+use_line+use_tag+'_dil'
  ms_flag=0                     ; flag for vaxis 1= m/s, 0= km/s (default)
  
;&;&;&;&:&;&;&;&;&;&:&;&
;& Process User Inputs
;&;&;&;&:&;&;&;&;&;&:&;&

  if keyword_set(data) then use_data=data
  if keyword_set(rms) then use_rms=rms
  if keyword_set(topdir) then use_topdir=topdir
  if keyword_set(line) then use_line=line
  if keyword_set(method) then use_method=strupcase(method)
  if keyword_set(cloud) then use_cloud=cloud
  use_name=use_cloud+'_'+use_line+'_'+use_tag+'_'+strlowcase(use_method)
  if keyword_set(name) then use_name=name
  if keyword_set(force_redo) then do_force_redo=1

;&;&;&;&:&;&;&;&;&;&:&;&
;& Check sanity of inputs
;&;&;&;&:&;&;&;&;&;&:&;&

  ; need to redefine this for 13CO/CS and 12m7mMagma
  cubefile=use_topdir+'/'+use_cloud+'/cubes/'+use_data+'/'+use_cloud+'_'+use_line+'_'+use_tag+'.pbcor.fits'
  status=file_test(cubefile,/read)
  
  if status eq 0 then $
     message,'Cannot read input file '+cubefile

;&;&;&;&:&;&;&;&;&;&:&;&
;& Read Data
;&;&;&;&:&;&;&;&;&;&:&;&

  data=readfits(cubefile,hdr)

  szdata=size(data)
  if szdata[0] ne 3 then message,'Data is not a cube?'
  naxis1=szdata[1] & naxis2=szdata[2] &   naxis3=szdata[3] 

  cdelt3 = sxpar(hdr, 'CDELT3',count=ct)
  if ct eq 0 then message,'Not enough information about spectral axis',/info
  if cdelt3 gt 100 then begin
     message,'It looks like the velocity units are m/s. Will convert to km/s',/info
     ms_flag=1
  end

  use_emchans=[0,naxis3]
  if keyword_set(emchans) then use_emchans=emchans

  outpkfile=use_topdir+'/'+use_cloud+'/cubes/'+use_data+'/'+use_name+'.pk.fits'
  outm0file=use_topdir+'/'+use_cloud+'/cubes/'+use_data+'/'+use_name+'.mom0.fits'
  outm1file=use_topdir+'/'+use_cloud+'/cubes/'+use_data+'/'+use_name+'.mom1.fits'
  outm2file=use_topdir+'/'+use_cloud+'/cubes/'+use_data+'/'+use_name+'.mom2.fits'

;&;&;&;&:&;&;&;&;&;&:&;&
;& Will redo unless all the data exists and force-redo not requested
;&;&;&;&:&;&;&;&;&;&:&;&

  if do_force_redo eq 0 then begin
     pkstatus=file_test(outpkfile,/read)
     m0status=file_test(outm0file,/read)
     m1status=file_test(outm1file,/read)
     m2status=file_test(outm2file,/read)
     if pkstatus eq 1 and m0status eq 1 and m1status eq 1 and m2status eq 1 and method eq 'DIL' then goto,dil_end
     if pkstatus eq 1 and m0status eq 1 and m1status eq 1 and m2status eq 1 and method eq 'SMO' then goto,smo_end
     if pkstatus eq 1 and m0status eq 1 and m1status eq 1 and m2status eq 1 and method eq 'BOSMA' then goto,bosma_end
     if pkstatus eq 1 and m0status eq 1 and m1status eq 1 and m2status eq 1 and method eq 'DIRECT' then goto,direct_end
  end

  ; do not overwrite -- this could be removed later if TW agrees?
  outcovfile=use_outdir+'/'+use_name+'.cov.fits'
  outpkfile=use_outdir+'/'+use_name+'.pk.fits'
  outvpkfile=use_outdir+'/'+use_name+'.vpk.fits'
  outm0file=use_outdir+'/'+use_name+'.mom0.fits'
  outm1file=use_outdir+'/'+use_name+'.mom1.fits'
  outm2file=use_outdir+'/'+use_name+'.mom2.fits'
  outem0file=use_outdir+'/'+use_name+'.emom0.fits'
  outem1file=use_outdir+'/'+use_name+'.emom1.fits'
  outem2file=use_outdir+'/'+use_name+'.emom2.fits'

;&;&;&;&:&;&;&;&;&;&:&;&
;& Parse RMS information
;&;&;&;&:&;&;&;&;&;&:&;&

 if keyword_set(use_rms) then begin $
    rms_size=size(use_rms)
    if rms_size[0] eq 0 and rms_size[1] eq 7 then use_rms=readfits(use_rms,rhdr)
    rms_size=size(use_rms)
    if rms_size[0] eq 0 and rms_size[1] eq 1 then use_rms=reform(data[*,*,0])*0.+use_rms
    ; still need to implement case of rms cube!
;    if rms_size[0] eq 2 then use_rms=reform(data[*,*,0]*0.+use_rms)
;    if rms_size[0] eq 3 then use_rms=median(
 end else begin
    use_rms=reform(data[*,*,0])*0.+mad(data)
    message,'RMS estimated from data: '+strtrim(string(median(use_rms)),2),/info
 end

;&;&;&;&:&;&;&;&;&;&:&;&
;& Generate 2D headers -- currently only used for BOSMA
;&;&;&;&:&;&;&;&;&;&:&;&

 hdrvpk=twod_head(hdr)
 hdrpk=twod_head(hdr)
 hdrcov=twod_head(hdr)
 hdrmom0=twod_head(hdr)
 hdrmom1=twod_head(hdr)
 hdrmom2=twod_head(hdr)
 hdremom0=twod_head(hdr)
 hdremom1=twod_head(hdr)
 hdremom2=twod_head(hdr)
 
 sxaddpar, hdrcov, 'BUNIT', 'coverage'
 sxaddpar, hdrpk, 'BUNIT', 'K'
 sxaddpar, hdrvpk, 'BUNIT', 'km/s'
 sxaddpar, hdrmom0, 'BUNIT', 'K.km/s'
 sxaddpar, hdrmom1, 'BUNIT', 'km/s'
 sxaddpar, hdrmom2, 'BUNIT', 'km/s'
 sxaddpar, hdremom0, 'BUNIT', 'K.km/s'
 sxaddpar, hdremom1, 'BUNIT', 'km/s'
 sxaddpar, hdremom2, 'BUNIT', 'km/s'
 
 sxaddpar, hdrcov, 'COMMENT', 'Coverage map'
 sxaddpar, hdrpk, 'COMMENT', 'Peak brightness map'
 sxaddpar, hdrvpk, 'COMMENT', 'Velocity of line profile peak'
 sxaddpar, hdrmom0, 'COMMENT', 'Moment-0 map'
 sxaddpar, hdrmom1, 'COMMENT', 'Moment-1 map'
 sxaddpar, hdrmom2, 'COMMENT', 'Moment-2 map'
 sxaddpar, hdremom0, 'COMMENT', 'Error in Moment-0 map'
 sxaddpar, hdremom1, 'COMMENT', 'Error in Moment-1 map'
 sxaddpar, hdremom2, 'COMMENT', 'Error in Moment-2 map'

 ;&;&;&;&:&;&;&;&;&;&:&;&
;& CALCULATE AN RMS AND COVERAGE MAP
;&;&;&;&:&;&;&;&;&;&:&;&

 if keyword_set(use_rms) then covmap=finite(use_rms)
 sxaddpar,hdrcov,'DATAMIN',0
 sxaddpar,hdrcov,'DATAMAX',1
 writefits,outcovfile,fix(covmap),hdrcov
 
 
;&;&;&;&:&;&;&;&;&;&:&;&
;& CALCULATE THE MOMENT MAPS ACCORDING TO REQUESTED METHOD
;&;&;&;&:&;&;&;&;&;&:&;&
 
  case use_method of
     'DIRECT': begin

       moment_maps,infile=cubefile $
                    ,namestr=use_name,outdir=use_outdir+'/',snrmsk=[0,0] $
                    ,pkmsk=[0,0],edgeblank=[0,0],/nostop

      direct_end:
     end
     'BOSMA': begin
        make_axes,hdr,/vonly,vaxis=vaxis
        Nsig=[5]
        deltaV=vaxis[1]-vaxis[0]

        for Nsigmaloop=0, N_ELEMENTS(Nsig)-1 do begin ; <--- the different thresholds that we want to test: 3sigma, 5sigma, etc.

           vpk_Nsig=FLTARR(naxis1,naxis2)
           pk_Nsig=FLTARR(naxis1,naxis2)
           mom0_Nsig=FLTARR(naxis1,naxis2)
           emom0_Nsig=FLTARR(naxis1,naxis2) ; emom0,1,2 are the corresponding uncertainties (error propagation)
           mom1_Nsig=FLTARR(naxis1,naxis2)
           emom1_Nsig=FLTARR(naxis1,naxis2)
           mom2_Nsig=FLTARR(naxis1,naxis2)
           emom2_Nsig=FLTARR(naxis1,naxis2)
           window_Nsig=FLTARR(naxis1,naxis2) ; here we will store the number of channels used for the window in each px
           
           sigmacut=Nsig[Nsigmaloop] ; threshold imposed to identify 'significant' pixels

           
           for x = 0, naxis1 -1 do begin
              for y = 0, naxis2 -1 do begin
                 
                 spectrumnow=reform(data[x,y,*]) ; extract spectrum for pixel x,y
                 peak=MAX(spectrumnow,chanVpeak,/nan) ; identify channel corresponding to intensity peak
                 
                 Vpeak=vaxis[chanVpeak] ; velocity correponding to the peak
                 
                 i=1            ; initialise counters for while loop
                 continuumnow=1000000000000. ; initialise with an arbitrary, huge number
                 contdiff=1.
                 looplimit=0.
                 chanVpeakorig=chanVpeak
                 Vpeakorig=Vpeak
                 
                 if FINITE(peak) NE 0. AND peak NE 0. AND peak GT sigmacut*use_rms[x,y] AND chanVpeak GT use_emchans[0] AND chanVpeak LT use_emchans[1] then begin

                    while (contdiff GT looplimit) AND (chanVpeak GT use_emchans[0] AND chanVpeak LT use_emchans[1]) do begin
                       
                       WINDOWmin=chanVpeak-i ; in each step, expand the window (initially one channel) by one channel on each side: [WINDOWmin,WINDOWmax]
                       WINDOWmax=chanVpeak+i
                       
                       if WINDOWmin LE 4 then WINDOWmin=(WINDOWmin+1>0) ; in the first steps, expand only by one channel in each step (instead of 2, for a finer adjustment)
                       if WINDOWmax GE naxis3-4 then WINDOWmax=(WINDOWmax-1 < naxis3-1); stop expanding if we dangerously approach the end of the spectrum
                       
                       free1=[0] & free2=[naxis3-1]
                       if (WINDOWmin-1) le naxis3-1 then free1 = spectrumnow[0:WINDOWmin-1]
                       if (WINDOWmax+1) lt naxis3-1 then free2 = spectrumnow[WINDOWmax+1:naxis3-1]
                       freeT = [free1,free2] ; range of channels outside the current window 
                       
                       continuumpre=continuumnow
                       continuumnow=MEAN(freeT,/NAN) ; average "continuum" flux in the channels outside current window
                       
                       lineT=spectrumnow[WINDOWmin:WINDOWmax] ; I define these to compute mom1 inside current window
                       vaxisline=vaxis[WINDOWmin:WINDOWmax]
                       
                       totalnow=TOTAL(lineT,/NAN)
                       mom1now=TOTAL(vaxisline*lineT,/NAN)/totalnow ; compute mom1 inside current window
                       

                                ; Now find what channel the new mom1now velocity corresponds to (for next iteration):

                       mom1channel=chanVpeakorig+ROUND((mom1now-Vpeakorig)/deltaV)
                       Vpeak=mom1now ; this will be used in the next iteration
                       chanVpeak=mom1channel

                                ; Update the variables that define the convergence criterion:

                       contdiff=continuumpre-continuumnow ; absolute difference in the mean continuum between current and previous step
                       ncontchannels=naxis3-1-2*i
                       looplimit=ABS(use_rms[x,y]/ncontchannels) ; ==> empirical convergence criterion from Bosma (1981)
                       i=i+1
                       
                    endwhile
                    
                    pknow=peak
                    vpknow=Vpeakorig
                    mom0now=totalnow*0.001*cdelt3
                    mom2now=SQRT(TOTAL(lineT*(vaxisline^2),/NAN)/TOTAL(lineT,/NAN)- (TOTAL(lineT*vaxisline,/NAN)/TOTAL(lineT,/NAN))^2.)
   
                                ; Calculate uncertainties (error propagation on the moment formulae):
                    
                    emom0now=use_rms[x,y]*0.001*cdelt3*SQRT(N_ELEMENTS(lineT))
                    
                    windownow=N_ELEMENTS(lineT)
                    
                    terms=FLTARR(N_ELEMENTS(vaxisline))
                    
                    for i=0,N_ELEMENTS(vaxisline)-1 do begin
                       terms[i]=(vaxisline[i]*totalnow-TOTAL(vaxisline*lineT,/NAN))/totalnow^2.*use_rms[x,y]
                    endfor
                    
                    emom1now=SQRT(TOTAL(terms^2.,/NAN))
                    
                    
                    terms=FLTARR(N_ELEMENTS(vaxisline))
                    
                    for i=0,N_ELEMENTS(vaxisline)-1 do begin
                       terms[i]=0.5/mom2now*(    ( vaxisline[i]^2.*totalnow - TOTAL(vaxisline^2.*lineT,/NAN) )/totalnow^2.      -     2.*( TOTAL(vaxisline*lineT,/NAN) * (vaxisline[i]*totalnow-TOTAL(vaxisline*lineT,/NAN)) )/totalnow^3.     )*use_rms[x,y]
                    endfor
                    
                    emom2now=SQRT(TOTAL(terms^2.,/NAN))

                 endif else begin

;                    message,'NaN sightlines',/info
                    pknow=!VALUES.F_NAN ; if we are not in a 'significant' pixel, output NaN
                    vpknow=!VALUES.F_NAN ; if we are not in a 'significant' pixel, output NaN
                    mom0now=!VALUES.F_NAN ; if we are not in a 'significant' pixel, output NaN
                    mom1now=!VALUES.F_NAN
                    mom2now=!VALUES.F_NAN
                    emom0now=!VALUES.F_NAN
                    emom1now=!VALUES.F_NAN
                    emom2now=!VALUES.F_NAN
                    windownow=!VALUES.F_NAN
                    
                 endelse
                 
   ;if (indices[0] NE -1) then mom1_Nsig[x,y]=TOTAL(vaxisaboveNsig*aboveNsig,/NAN)/totalnow else mom1_Nsig[x,y]=!VALUES.F_NAN
                 
                 pk_Nsig[x,y]=pknow
                 vpk_Nsig[x,y]=vpknow
                 mom0_Nsig[x,y]=mom0now
                 mom1_Nsig[x,y]=mom1now
                 mom2_Nsig[x,y]=mom2now
                 emom0_Nsig[x,y]=emom0now
                 emom1_Nsig[x,y]=emom1now
                 emom2_Nsig[x,y]=emom2now
                 window_Nsig[x,y]=windownow
                 
              endfor
           endfor

           if ms_flag eq 1 then begin
              vpk_Nsig=1.e-3*vpk_Nsig
              mom1_Nsig=1.e-3*mom1_Nsig
              mom2_Nsig=1.e-3*mom2_Nsig
              emom1_Nsig=1.e-3*emom1_Nsig
              emom2_Nsig=1.e-3*emom2_Nsig
           end
           
           sxaddpar,hdrpk,'DATAMIN',min(pk_Nsig,/nan)
           sxaddpar,hdrpk,'DATAMAX',max(pk_Nsig,/nan)
           sxaddpar,hdrvpk,'DATAMIN',min(vpk_Nsig,/nan)
           sxaddpar,hdrvpk,'DATAMAX',max(vpk_Nsig,/nan)

           sxaddpar,hdrmom0,'DATAMIN',min(mom0_Nsig,/nan)
           sxaddpar,hdrmom0,'DATAMAX',max(mom0_Nsig,/nan)
           sxaddpar,hdremom0,'DATAMIN',min(emom0_Nsig,/nan)
           sxaddpar,hdremom0,'DATAMAX',max(emom0_Nsig,/nan)

           sxaddpar,hdrmom1,'DATAMIN',min(mom1_Nsig,/nan)
           sxaddpar,hdrmom1,'DATAMAX',max(mom1_Nsig,/nan)
           sxaddpar,hdremom1,'DATAMIN',min(emom1_Nsig,/nan)
           sxaddpar,hdremom1,'DATAMAX',max(emom1_Nsig,/nan)

           sxaddpar,hdrmom2,'DATAMIN',min(mom2_Nsig,/nan)
           sxaddpar,hdrmom2,'DATAMAX',max(mom2_Nsig,/nan)
           sxaddpar,hdremom2,'DATAMIN',min(emom2_Nsig,/nan)
           sxaddpar,hdremom2,'DATAMAX',max(emom2_Nsig,/nan)


           writefits,outvpkfile,vpk_Nsig,hdrvpk
           writefits,outpkfile,pk_Nsig,hdrpk

           writefits,outm0file,mom0_Nsig,hdrmom0
           writefits,outm1file,mom1_Nsig,hdrmom1
           writefits,outm2file,mom2_Nsig,hdrmom2

           writefits,outem0file,emom0_Nsig,hdremom0
           writefits,outem1file,emom1_Nsig,hdremom1
           writefits,outem2file,emom2_Nsig,hdremom2
        
        endfor                  ; <--- over Nsigmaloop (repeat whole process for the different thresholds that we want to test: 3sigma, 5sigma, etc.)

        bosma_end:
     end
     'DIL': begin

        make_cprops_mask, indata = data $
                          , outmask = sigmask $
                          , hi_thresh = 5 $
                          , lo_thresh = 2 $
                          , hi_nchan=3 $
                          , lo_nchan=4

        final_mask=(sigmask)
        sxaddpar,hdr,'DATAMIN',0
        sxaddpar,hdr,'DATAMAX',1 
        writefits,use_outdir+'/'+use_cloud+'_'+use_line+'_'+use_tag+'_emission_mask.fits',fix(final_mask),hdr

        blankidx=where(finite(final_mask) eq 0)
        final_mask[blankidx]=0
        writefits,use_outdir+'/'+use_cloud+'_'+use_line+'_'+use_tag+'_emission_mask_zeroes.fits',fix(final_mask),hdr

        print,'Flux in this mask: ', total(data*final_mask,/nan)/total(data,/nan)
           
        moment_maps,infile=cubefile $
                    ,inmask=use_outdir+'/'+use_cloud+'_'+use_line+'_'+use_tag+'_emission_mask.fits' $
                    ,namestr=use_name,outdir=use_outdir+'/',snrmsk=[0,0] $
                    ,pkmsk=[0,0],/auto_rng,edgeblank=[0,0],/nostop

        dil_end:
        stop
     end

     'SMO': begin
        make_cprops_mask, indata = data $
                     , outmask = sigmask $
                     , hi_thresh = 5 $
                     , lo_thresh = 2 $
                     , hi_nchan=3 $
                     , lo_nchan=4

        conv_with_gauss $
           , data=data $
           , hdr=hdr $
           , out_data=cube_cvl $
           , out_hdr=hdr_cvl $
           ,target_beam=[5.,5.,0] 

        sxaddpar,hdr,'DATAMIN',min(cube_cvl,/nan)
        sxaddpar,hdr,'DATAMAX',max(cube_cvl,/nan)
        writefits,use_outdir+'/'+use_cloud+'_'+use_line+'_'+use_tag+'_smo5as_cube.fits',cube_cvl,hdr_cvl

        make_cprops_mask, indata = cube_cvl $
                          , outmask = sigmask_cvl $
                          , hi_thresh = 4 $
                          , lo_thresh = 2 $
                          , hi_nchan=3 $
                          , lo_nchan=4

  
        final_mask=(sigmask_cvl or sigmask)
        sxaddpar,hdr,'DATAMIN',0
        sxaddpar,hdr,'DATAMAX',1 
        writefits,use_outdir+'/'+use_cloud+'_'+use_line+'_'+use_tag+'_emission_mask.fits',fix(final_mask),hdr

        blankidx=where(finite(final_mask) eq 0)
        final_mask[blankidx]=0
        writefits,use_outdir+'/'+use_cloud+'_'+use_line+'_'+use_tag+'_emission_mask_zeroes.fits',fix(final_mask),hdr

        print,'Flux in this mask: ', total(data*final_mask,/nan)/total(data,/nan)

        pk_cvl=max(cube_cvl,vpkidx,dim=3,/nan)
        make_simple_noisemap, cube_in=cube_cvl $
                              , out_map=rms2d_cvl $
                              , channels=[2,2] $
                              , box=3 $
                              , meannoise=noise_estimate_2dcvl
        message,'Estimate of average noise: '+strtrim(string(noise_estimate_2dcvl),2),/info
        snrpk_cvl=pk_cvl/rms2d_cvl
  
        moment_maps,incube=cube_cvl, inhdr=hdr_cvl $
                    ,inmask=sigmask_cvl $
                    , snrmsk=[0,0],pkmsk=[0,0],edgeblank=[0,0] $
                    ,namestr=use_name+'_5p0cvl',outdir=use_outdir+'/',/nostop
              
        moment_maps,infile=cubefile $
                    ,inmask=use_outdir+'/'+use_cloud+'_'+use_line+'_'+use_tag+'_emission_mask.fits' $
                    ,namestr=use_name,outdir=use_outdir+'/',snrmsk=[0,0] $
                    ,pkmsk=[0,0],/auto_rng,edgeblank=[0,0],/nostop

        smo_end:

                    STOP
                    
                 end
     
     
     
     else: begin
        message,'Unknown method?'
     end
     
  endcase
  
  

  
the_end:
end
 

