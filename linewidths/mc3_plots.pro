pro mc3_plots,save_file=save_file,mom1=mom1,mom2=mom2,namestr=namestr,plotdir=plotdir,png=png,nostop=nostop

  use_win=0L
  define_la_common
  
  use_save_file='../sav_data/GMC1_12CO_3p5as_bosma.idl'
  use_plotdir='../plots/'
  use_mom1='../moment_maps/GMC1_12CO_3p5as_bosma.mom1.fits'
  use_mom2='../moment_maps/GMC1_12CO_3p5as_bosma.mom2.fits'
  use_namestr='GMC1_12CO_3p5as_bosma'
  
  if keyword_set(save_file) then use_save_file=save_file
  if keyword_set(namestr) then use_namestr=namestr
  if keyword_set(plotdir) then use_plotdir=plotdir
  if keyword_set(mom1) then use_mom1=mom1
  if keyword_set(mom2) then use_mom2=mom2

  if keyword_set(png) then begin
     use_m1pngfile=use_plotdir+'/'+use_namestr+'.mom1.png'
     use_m2pngfile=use_plotdir+'/'+use_namestr+'.mom2.png'
     use_vchistpngfile=use_plotdir+'/'+use_namestr+'_m1_vchisto.png'
     use_vdhistpngfile=use_plotdir+'/'+use_namestr+'_m2_vdhisto.png'
     use_lwpngfile=use_plotdir+'/'+use_namestr+'_lwpixprop.png'
  end
  

  multiscale:
  restore,use_save_file

  window,use_win,xsize=600,ysize=600  & use_win=use_win+1
  loadct,39 & reversect

  scales=samp_str.scale_as
  sscales=scales(sort(scales))
  uscales=sscales[rem_dup(sscales)]
  nsca=n_elements(uscales)
  fgcol=255
  cols=range_gen(nsca,[30,255-30])
  xr=[min(samp_str.vdisp_wgt,/nan)>0,max(samp_str.vdisp_wgt,/nan)]
  yr=[min(samp_str.wgt_spec_eqw,/nan)>0,max(samp_str.wgt_spec_eqw,/nan)]
  
  tit=use_namestr+' multi lws'
  for j=0,nsca-1 do begin
     xx=where(samp_str.scale_as eq uscales[j] and finite(samp_str.vdisp_wgt) eq 1 and finite(samp_str.wgt_spec_eqw) eq 1,ct)
     if ct lt 1 then goto,skip
     if j eq 0 then begin
        plot,[samp_str[xx].vdisp_wgt],[samp_str[xx].wgt_spec_eqw],/ysty,/xsty $
               ,/nodata,xtit='2D Wgt Avg Velocity Dispersion [km/s]',ytit='Wgt Avg Shuffle Spectrum EQW [km/s]' $
               ,tit=tit,xr=xr,yr=yr,col=fgcol $
             ,xthick=2,ythick=2,thick=2,charsize=1.5,charthick=1.4
        equality,lines=1,thick=2,col=fgcol
        end
     PLOTSYM, 0 ,2., /FILL     ;Plotting symbol is a filled circle,
      plots,[samp_str[xx].vdisp_wgt],[samp_str[xx].wgt_spec_eqw],psym=8,col=fgcol
     PLOTSYM, 0 ,1.5, /FILL     ;Plotting symbol is a filled circle,
     plots,[samp_str[xx].vdisp_wgt],[samp_str[xx].wgt_spec_eqw],psym=8,col=cols[j]
     skip:
     print,j
  endfor

  PLOTSYM, 0 ,2., /FILL         ;Plotting symbol is a filled circle,
 al_legend, /top, /left, box=0, clear=0 $
             , strtrim(string(fix(uscales)),2)+'as' $
             , lines=0 $
             , thick=4 $
             , psym=8 $
             , textcolor=fgcol $
             , color=fix(cols) $
            , charsize=1.5, charthick=1.4

  if keyword_set(png) then write_png,use_lwpngfile,TVRD(/TRUE)

  histograms:
  m1=readfits(use_mom1,m1hdr)
  m2=readfits(use_mom2,m2hdr)
  la_m1=m1
  badidx=where(finite(la_m1) eq 0 or la_m1 lt -500 or la_m1 gt 500, bct,comp=goodidx)
  la_m1[badidx]=la_undef()
  m1[badidx]=!values.f_nan
  
  la_m2=m2
  badidx=where(finite(la_m2) eq 0 or la_m2 lt 0, bct,comp=goodidx)
  la_m2[badidx]=la_undef()
  m2[badidx]=!values.f_nan

  m1_vcdisp=la_sigma(la_m1,n_sigma=5.,/med,mean=m1_vcavg)
  vmin=la_min(la_m1) 
  vmax=la_max(la_m1)
  vmin=percentile(m1[goodidx],99.5)
  vmax=percentile(m1[goodidx],0.5) 
  
  finidx=where(la_m1 ne la_undef(),finct)
  nbins=round(finct/1000.)
  xr=[0.99*vmin,1.01*vmax]
  !p.position=[0.2,0.2,0.8,0.8]
  hh=histogram(la_m1[finidx],max=vmax,min=vmin,locations=vbins,nbins=nbins)
  tit=use_namestr+' moment-1'

  window,use_win,xsize=600,ysize=600  & use_win=use_win+1
  yr=[0,1.2*(max(hh)/float(total(hh,/nan)))]
  cgplot,vbins,(hh/float(total(hh,/nan))),/ysty,/xsty $
         ,/nodata,xr=xr,yr=yr,xtit='Radial Velocity [km/s]',ytit='% Npixels',tit=tit $
         ,xthick=2,ythick=2,thick=2,charsize=1.4,charthick=1.7
  cgplot,vbins,(hh/float(total(hh,/nan))),psym=10,/overplot,thick=2

  al_legend, /top, /left, box=0, clear=0 $
            , ['Mean M1 Value: '+strtrim(string(m1_vcavg),2), $
              'RMS M1 Values: '+strtrim(string(m1_vcdisp),2)] $
             , lines=-99 $
             , textcolor=fgcol $
            , charsize=1.5, charthick=1.4

  if keyword_set(png) then write_png,use_vchistpngfile,TVRD(/TRUE)

  
  m2_disp=la_sigma(la_m2,n_sigma=5.,/med,mean=m2_avg)
  vmin=la_min(la_m2) 
  vmax=la_max(la_m2)
  finidx=where(la_m2 ne la_undef(),finct)
  nbins=round(finct/1000.)
  xr=[0.99*vmin,1.01*vmax]
  !p.position=[0.2,0.2,0.8,0.8]
  hh=histogram(la_m2[finidx],max=vmax,min=vmin,locations=vbins,nbins=nbins)
  tit=use_namestr+' moment-2'

  window,use_win,xsize=600,ysize=600  & use_win=use_win+1
  yr=[0,1.2*(max(hh)/float(total(hh,/nan)))]
  cgplot,vbins,(hh/float(total(hh,/nan))),/ysty,/xsty $
         ,/nodata,xr=xr,yr=yr,xtit='Velocity Dispersion [km/s]',ytit='% Npixels',tit=tit $
         ,xthick=2,ythick=2,thick=2,charsize=1.4,charthick=1.7
  cgplot,vbins,(hh/float(total(hh,/nan))),psym=10,/overplot,thick=2

    al_legend, /top, /left, box=0, clear=0 $
            , ['Mean M2 Value: '+sigfig(m2_avg,3), $
              'RMS M2 Values: '+sigfig(m2_disp,3)] $
             , lines=-99 $
             , textcolor=fgcol $
            , charsize=1.5, charthick=1.4


  
  if keyword_set(png) then write_png,use_vdhistpngfile,TVRD(/TRUE)

  maps:
    !p.position=[0.15,0.15,0.85,0.85]

  m1=readfits(use_mom1,m1hdr)
  make_axes,m1hdr,raxis=raxis,daxis=daxis

  la_m1=m1
  badidx=where(finite(la_m1) eq 0 or la_m1 lt -500 or la_m1 gt 500, bct,comp=goodidx)
  la_m1[badidx]=la_undef()
  m1[badidx]=!values.f_nan
  
  obp=[1.1,0,1.13,1]
  imrange=[percentile(m1[goodidx],99),percentile(m1[goodidx],1)]

  window,use_win,xsize=600,ysize=600 & use_win=use_win+1
  loadct,39 & fgcolor=255 & reversect
  image_cont20,la_m1,m1hdr,/square,/silent,tit=use_namestr,off_bar_pos=obp,/nologo,imrange=imrange,bar_tit='[km/s]'
  if keyword_set(png) then write_png,use_m1pngfile,tvrd(/true)

  m2=readfits(use_mom2,m2hdr)
  make_axes,m2hdr,raxis=raxis,daxis=daxis

  la_m2=m2
  badidx=where(finite(la_m2) eq 0 or la_m2 lt 0, bct,comp=goodidx)
  la_m2[badidx]=la_undef()
  m2[badidx]=!values.f_nan

  obp=[1.1,0,1.13,1]
  imrange=[percentile(m2[goodidx],99),percentile(m2[goodidx],1)]

  window,use_win,xsize=600,ysize=600 & use_win=use_win+1
  loadct,39 & fgcolor=255 & reversect
  image_cont20,la_m2,m2hdr,/square,/silent,tit=use_namestr,off_bar_pos=obp,/nologo,imrange=imrange,bar_tit='[km/s]'
  if keyword_set(png) then write_png,use_m2pngfile,tvrd(/true)
 
  if not keyword_set(nostop) then stop
  
the_end:
  
end
