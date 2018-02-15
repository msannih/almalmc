pro mc3_sample_regions,sample_as=sample_as $
                       ,ico=ico,pk=pk,mom1=mom1,mom2=mom2,cube=cube $
                       ,save=save,scube=scube,mask=mask $
                       ,append=append,show=show,rebaseline=rebaseline,nostop=nostop

;&;&;&;&:&;&;&;&;&;&:&;&
;& Initialise LA common block
;&;&;&;&:&;&;&;&;&;&:&;&

   pdp_define_la_common

;&;&;&;&:&;&;&;&;&;&:&;&
;& Defaults
;&;&;&;&:&;&;&;&;&;&:&;&

   use_win=0L
   use_ico='../moment_maps/'+['GMC1_12CO_3p5as_bosma.mom0.fits']
   use_pkco='../moment_maps/'+['GMC1_12CO_3p5as_bosma.pk.fits']
   use_mom1='../moment_maps/'+['GMC1_12CO_3p5as_bosma.mom1.fits']
   use_mom2='../moment_maps/'+['GMC1_12CO_3p5as_bosma.mom2.fits']
   use_mask='../moment_maps/'+['GMC1_12CO_3p5as_bosma.cov.fits']
   use_cube='../orig_data/GMC1/cubes/match_3p5/'+['GMC1_12CO_3p5as.pbcor.fits']
   use_scube='../moment_maps/'+['GMC1_12CO_3p5as_bosma.scube.fits']
   use_save='../sav_data/'+['GMC1_12CO_3p5as_bosma.idl']
   use_sample_as=8.             ; minimum separation of sampling points
   use_minpix=5                 ; only do averaging for apertures with more than N pixels
   use_epsilon=1.e-15           ; a very small number for downweighting bad pixels
   do_show=0
   do_append=0
   do_rebaseline=0
   
;&;&;&;&:&;&;&;&;&;&:&;&
;&Process user inputs
;&;&;&;&:&;&;&;&;&;&:&;&

   if keyword_set(show) then do_show=1
   if keyword_set(append) then do_append=1
   if keyword_set(rebaseline) then do_rebaseline=1
   if keyword_set(sample_as) then use_sample_as=sample_as
   if keyword_set(ico) then use_ico=ico
   if keyword_set(pkco) then use_pkco=pkco
   if keyword_set(mom1) then use_mom1=mom1
   if keyword_set(mom2) then use_mom2=mom2
   if keyword_set(mask) then use_mask=mask
   if keyword_set(cube) then use_cube=cube
   if keyword_set(scube) then use_scube=scube
   if keyword_set(save) then use_save=save

;&;&;&;&:&;&;&;&;&;&:&;&
;& Read Data
;&;&;&;&:&;&;&;&;&;&:&;&

  icomap=readfits(use_ico,icohdr)
  pkcomap=readfits(use_pkco,pkcohdr)
  mom1map=readfits(use_mom1,mom1hdr)
  mom2map=readfits(use_mom2,mom2hdr)
  cocube=readfits(use_cube,chdr)
  mask=readfits(use_mask,mhdr)

  deltav=sxpar(chdr,'CDELT3',count=dct)
  nchans=sxpar(chdr,'NAXIS3',count=nct)
  nx=sxpar(chdr,'NAXIS1',count=xct)
  ny=sxpar(chdr,'NAXIS2',count=yct)
  cd1=sxpar(chdr,'CDELT1',count=c1ct)
  cd2=sxpar(chdr,'CDELT2',count=c2ct)

  if dct eq 0 or nct eq 0 then message,'Not enough information about velocity axis in cube hdr?'
  shuffle_vaxis = (findgen(nchans)-(nchans/2))*deltav

  use_field_deg=(abs(cd1)*nx) > (abs(cd2)*ny)
  use_sample_pix=round(use_sample_as/(abs(cd1)*3600))

  use_specplot_max=max(pkcomap,/nan)
  use_icoplot_max=max(icomap,/nan)

;&;&;&;&:&;&;&;&;&;&:&;&
;& Rebaseline if requsested
;&;&;&;&:&;&;&;&;&;&:&;&

  if keyword_set(do_rebaseline) then begin
     fitmask=cocube*1.
     use_rms=mad(cocube)
     fitpix=where(cocube lt 5.*use_rms)
     use_rms=mad(cocube[fitpix])
     fitpix=where(cocube lt 2.*use_rms)
     fitmask[fitpix]=0.
     mc3_rebaseline,idl_in=cocube,hdr_in=chdr,idl_mask=fitmask,idl_out=cocube_rebase,hdr_out=chdr_rebase,order=1
     cocube=cocube_rebase & chdr=chdr_rebase
  end
  
  
;&;&;&;&:&;&;&;&;&;&:&;&
; weighting -- use LA routines
;&;&;&;&:&;&;&;&;&;&:&;&
  la_icomap=icomap
  infidx=where(finite(icomap) eq 0, infct)
  if infct gt 0 then la_icomap[infidx]=la_undef()

     
; MAKE SAMPLING GRID
  make_axes, icohdr, raxis=raxis,daxis=daxis,rimg=rimg, dimg=dimg
  x_ctr = mean(rimg,/nan) & y_ctr = mean(dimg,/nan)

; SET GRID SPACING AT HALF OF THE TARGET RESOLUTION
     
     hex_grid $
     , ctr_x = x_ctr, ctr_y = y_ctr, /radec $
     , spacing = use_sample_as/(2. * 3600.) $
     , r_limit = 2*use_field_deg $
     , xout = ra_samp, yout = dec_samp $
     , /center

      ; REMOVE POINTS OUTSIDE THE CO MAP
      adxy, icohdr, ra_samp, dec_samp, x_samp, y_samp
      x_samp = round(x_samp)
      y_samp = round(y_samp)
      map_sz = size(icomap)
      keep = where(x_samp ge 0 and $
                   y_samp ge 0 and $
                   x_samp lt map_sz[1] and $
                   y_samp lt map_sz[2] and $
                   finite(icomap[x_samp,y_samp]))
      ra_samp = ra_samp[keep]
      dec_samp = dec_samp[keep]

      ; GENERATE FINAL SAMPLING POINTS
      adxy, icohdr, ra_samp, dec_samp, x_samp, y_samp
      n_pts = n_elements(x_samp)
      message,"Found "+sigfig(n_pts,4)+" sampling points",/info

      if keyword_set(do_show) then begin
         window,use_win & use_win=use_win+1
         loadct,39 & fgcolor=255
         disp,icomap,raxis=raxis,daxis=daxis,min=-2.,max=use_icoplot_max
         window,use_win & use_win=use_win+1
         plot,x_samp,y_samp,psym=3,/xsty,/ysty,xr=[0,map_sz[1]],yr=[0,map_sz[2]]
         window,use_win & use_win=use_win+1
         plot,ra_samp,dec_samp,psym=3,/xsty,/ysty,xr=[max(ra_samp,/nan),min(ra_samp,/nan)],yr=[min(dec_samp,/nan),max(dec_samp,/nan)]
         wait,3
      end
      

      ; SHUFFLE CUBE
      ;    create velocity axis for stacking
     make_axes, chdr, vaxis=vaxis_co, /vonly

     shuffle_cube=shuffle(spec=cocube $
                          ,vaxis=vaxis_co/1.e3 $
                          ,zero=mom1map $
                          ,target_vaxis=shuffle_vaxis/1.e3 $
                          ,missing=!values.f_nan)
     schdr=chdr
     sxaddpar,schdr,'HISTORY','Shuffled cube on moment-1 template'
     writefits,use_scube,shuffle_cube,schdr ; IDL smart enough to update Nchannels in HDR


      sstruct=mc3_empty_sampling_struct()
      this_sstruct=replicate(sstruct,n_pts)

      for j=0,n_pts-1 do begin
         
         ; position information
         this_sstruct[j].scale_as=use_sample_as
         this_sstruct[j].region_size=use_sample_pix
         this_sstruct[j].x=x_samp[j]
         this_sstruct[j].y=y_samp[j]
         this_sstruct[j].ra=ra_samp[j]
         this_sstruct[j].dec=dec_samp[j]

         ; get indices for this region
         dummy=mc3_mean_box(icomap,icohdr,(x_samp[j]),(y_samp[j]),0,use_sample_pix,stdev,npix,indices=idx,/silent)
         this_sstruct[j].npix=npix
         this_sstruct[j].idx=ptr_new(idx)

         if n_elements(idx) lt use_minpix then goto, skip

; weighted 2D high resolution measurements
         this_sstruct[j].ico_wgt=mc3_wavg(icomap[idx],icomap[idx],/weight,/forcepos)
         this_sstruct[j].tpk_wgt=mc3_wavg(pkcomap[idx],icomap[idx],/weight,/forcepos)
         this_sstruct[j].vdisp_wgt=mc3_wavg(mom2map[idx],icomap[idx],/weight,/forcepos)
         this_sstruct[j].ico_avg=mean(icomap[idx],/nan)
         this_sstruct[j].tpk_avg=mean(pkcomap[idx],/nan)
         this_sstruct[j].vdisp_avg=mean(mom2map[idx],/nan)

; combine shuffled spectra for this region
         twod_ij=array_indices(icomap,idx)
         ; extract corresponding region of weights map 
         wgt2d=la_mul(la_icomap,0)
         wgt2d[idx]=la_icomap[idx]
         
         la_shuffle_cube=shuffle_cube
         infidx=where(finite(shuffle_cube) eq 0, infct)
         if infct gt 0 then la_shuffle_cube[infidx]=la_undef()
         cubetmp=la_mul(la_shuffle_cube,0)

         for kk=0L,n_elements(idx)-1 do $
            cubetmp[twod_ij[0,kk],twod_ij[1,kk],*]=la_shuffle_cube[twod_ij[0,kk],twod_ij[1,kk],*]

         dims=size(cubetmp)

         ; strongly downweight spectra with negative I(CO)
         negidx=where(wgt2d lt 0,negct)
         if negct gt 0 then wgt2d[negidx]=use_epsilon
         
         var=cubetmp & wgt=rebin(wgt2d,[dims[1],dims[2],dims[3]])
         wgt_spec=la_div(la_tot(la_mul(var,wgt),dim=1),la_tot(wgt,dim=1)) ; 

         sum_spec=la_tot(var,dim=1)
         avg_spec=la_tot(var,dim=1)/npix
  
         wgt_spec_mom0=1.e-3*deltav*la_tot(wgt_spec)
         wgt_spec_pk=la_max(wgt_spec)
         wgt_spec_eqw=wgt_spec_mom0/(sqrt(2*!pi)*wgt_spec_pk)

         avg_spec_mom0=1.e-3*deltav*la_tot(avg_spec)
         avg_spec_pk=la_max(avg_spec)
         avg_spec_eqw= avg_spec_mom0/(sqrt(2*!pi)*avg_spec_pk)

         if keyword_set(do_show) then begin
            window,use_win & use_win=use_win+1
            plot,shuffle_vaxis/1.e3,avg_spec,/xsty,/yst,yr=[-0.5,use_specplot_max],color=200
            oplot,shuffle_vaxis/1.e3,wgt_spec,color=254
            oplot,shuffle_vaxis/1.e3,wgt_spec*0.0,color=fgcolor,lines=1
            al_legend, /top, /right, box=0, clear=0 $
                       ,'EQW: '+sigfig([wgt_spec_eqw,avg_spec_eqw],3) $
                       ,lines=[0,0] $
                       ,color=[254,200]
         end
         

         this_sstruct[j].avg_spec=ptr_new(avg_spec)
         this_sstruct[j].wgt_spec=ptr_new(wgt_spec)
         this_sstruct[j].wgt_spec_eqw=wgt_spec_eqw
         this_sstruct[j].wgt_spec_mom0=wgt_spec_mom0
         this_sstruct[j].wgt_spec_pk=wgt_spec_pk
         this_sstruct[j].avg_spec_eqw=avg_spec_eqw
         this_sstruct[j].avg_spec_mom0=avg_spec_mom0
         this_sstruct[j].avg_spec_pk=avg_spec_pk
         this_sstruct[j].wgt_spec_b=wgt_spec_mom0/(wgt_spec_eqw*wgt_spec_eqw)
         this_sstruct[j].avg_spec_b=avg_spec_mom0/(avg_spec_eqw*avg_spec_eqw)

         this_sstruct[j].region_size=use_sample_pix
         this_sstruct[j].npix=npix
         this_sstruct[j].idx=ptr_new(idx)

         skip:
         use_win=0L
      endfor
      
;&;&;&;&:&;&;&;&;&;&:&;&
;& Save Output
;&;&;&;&:&;&;&;&;&;&:&;&
;==== ADD STRUCTURE TO END OF EXISTING STRUCTURE
;##############################
      
  if keyword_set(do_append) then begin
     restore, use_save
     samp_str=[samp_str,this_sstruct]
  endif else begin
     samp_str=this_sstruct
  endelse
  save,file=use_save,samp_str
      
  if not keyword_set(nostop) then  stop
  the_end:
end

