pro mc3_rebaseline,datadir=datadir,fits_in=fits_in $
                    ,idl_in=idl_in,hdr_in=hdr_in,idl_mask=idl_mask $
                    ,outdir=outdir,fits_out=fits_out $
                    ,idl_out=idl_out,hdr_out=hdr_out $
                    ,blflags=blflags,blcoeffs=blcoeffs $
                    ,order=order, help=help,verbose=verbose 


;+
; NAME:
;   mc3_rebaseline
;
; PURPOSE:
;	Rebaseline profiles in a spectral line cube
;
; CALLING SEQUENCE:
;   mc3_rebaseline, fits_in=fits_in, idl_in=idl_in, hdr_in=hdr_in, $
;              fits_out=fits_out, fits_mask=fits_mask, $
;              idl_mask=idl_mask, order=order, blflags=blflags, blcoeffs=blcoeffs 
;
; INPUT PARAMETERS:
;    FITS_IN: Name for input data cube in standard fits format (String)
;    IDL_IN:  Name for input data cube in IDL format (String)
;    HDR_IN:     IDL hdr if passing a cube as IDL_IN
;    MASK_FITS: 1/0 Mask in standard fits format (string) 
;               We will fit baselines using pixels where mask=0
;    MASK_IDL: 1/0 Mask in IDL format
;    ORDER: order of polynomial to use for baseline fitting (integer =  0,1,2)
;
; OUTPUTS:
;    FITS_OUT: Rebaselined data cube in standard fits format (string, without .fits) 
;    BLFLAGS: [x,y] 2D array with flags reporting status of baseline
;              fit for that pixel
;    BLCOEFFS: [x,y,ON] array with values of coefficients C0, C1, C2
;              etc of the baseline fit for that pixel
;
; ACCEPTED KEY-WORDS:
;     help = print this help
;     verbose = print extra information to screen
;
; REQUIRES:
;    Goddard Astron, Freudenreich routines
;
; MODIFICATION HISTORY:
;    written 25-11-2016 by AH
; COMMENTS AND TO DO LIST:
;    Minimally tested!
;    Error trapping not implemented.
;    Assumes input cube will have matching (x,y,v) astrometric grid as input mask
;-

  IF keyword_set(help) THEN BEGIN
     doc_library,'mc3_rebaseline'
     goto,the_end
  ENDIF

  nan = !values.f_nan

;==============
; defaults
;==============

  use_outdir='./'
  use_datadir='./'
  use_order=0

; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
; PROCESS USER INPUT
; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  if keyword_set(outdir) then use_outdir=outdir
  if keyword_set(datadir) then use_datadir=datadir
  if keyword_set(order) then use_order=order

;==============
; read data
;==============
  
  if keyword_set(fits_in) then data_in=readfits(use_datadir+fits_in,hdr)
  if keyword_set(fits_mask) then mask_in=readfits(use_datadir+fits_mask,jnk)
  if keyword_set(idl_in) then data_in=idl_in
  if keyword_set(hdr_in) then hdr=hdr_in
  if keyword_set(idl_mask) then mask_in=idl_mask

  if n_elements(data_in) eq 0 then begin
     message, 'Problem reading input data',/info
     goto, the_end
  end

  if n_elements(hdr) eq 0 then begin
     message, 'I need header information',/info
     goto, the_end
  end

; replace NaNs in mask with zeros
  badmask=where(finite(mask_in) eq 0, nbad)
  if nbad gt 0 then mask_in[badmask] = 0

; size of cube/FoV
  sz=size(data_in)
  xmax=sz[1] & ymax=sz[2] & vmax=sz[3]
  nbltot=float(xmax)*ymax

; we will make a 2d map of how the baselines were treated
  blflags=data_in[*,*,0] & blflags[*]=-1.

; we will make maps of the constant + coefficients of x^1, x^2 etc.
  blcoeffs=fltarr(xmax,ymax,use_order+1)


; ---------------------- sanity checking ------------------------------

  smask=size(mask_in)
  if total(smask[1:3]-sz[1:3]) ne 0 then begin
     message,"Mask and cube do not have same dimensions",/info
     goto, the_end
  endif

;--------------------- start the loop -------------------------

  blcount=0.
  data_out=data_in
  data_out[*]=nan

  for j=0, ymax-1 do begin
     for i=0, xmax-1 do begin

      spec=reform(data_in[i,j,*])
      mspec=reform(mask_in[i,j,*])
      
; we use channels without emission (mask=0) for the baselines
      usechans= where(mspec eq 0 and finite(spec), nch)
      emchans= where(mspec eq 1 and finite(spec), nech)
      
      if nch gt 0 then begin
         
         blcount=blcount+1
         
         CASE use_order of
            0: begin
               offset=median(spec(usechans),/even)
               data_out[i,j,*]=spec-offset
               blflags[i,j]=0
               blcoeffs[i,j,0]=offset
            end
            1: begin
               y=spec
               x=indgen(vmax)
               coeff=robust_regress(transpose(x(usechans)),y(usechans),yfit)
               if n_elements(coeff) eq 2 then begin
                  blfit=coeff[0]+coeff[1]*x
                  data_out[i,j,*]=spec-blfit
                  blflags[i,j]=1
                  blcoeffs[i,j,0]=coeff[0]
                  blcoeffs[i,j,1]=coeff[1]
                                ; if regression fails, apply a simple offset
               end else if n_elements(coeff) lt 2 or n_elements(coeff) gt 2 then begin
                  offset=median(spec(usechans),/even)
                  data_out[i,j,*]=spec-offset
                  blflags[i,j]=0
                  blcoeffs[i,j,0]=coeff[0]
               endif
            end
            2: begin
               y=spec
               x=indgen(vmax)
               x2=x*x
               coeff=robust_regress([transpose(x(usechans)),transpose(x2(usechans))],y(usechans),yfit)
               if n_elements(coeff) eq 3 then begin
                  blfit=coeff[0]+coeff[1]*x+coeff[2]*x2
                  data_out[i,j,*]=spec-blfit
                  blflags[i,j]=2
                  blcoeffs[i,j,0]=coeff[0]
                  blcoeffs[i,j,1]=coeff[1]
                  blcoeffs[i,j,2]=coeff[2]
                                ; if regression fails, apply a simple offset
               end else if n_elements(coeff) lt 3 or n_elements(coeff) gt 3 then begin
                  offset=median(spec(usechans),/even)
                  data_out[i,j,*]=spec-offset
                  blflags[i,j]=0
                  blcoeffs[i,j,0]=coeff[0]
               endif
            end
            else: begin
               print,'Only baselines up to order 2 implemented at present'
            end
         endcase
      end else begin
         if keyword_set(verbose) then $
            print,'Skipping this spectrum because there are no emission-free channels'
      end
       
   endfor
  endfor
  
;-------------- WRITE OUT FITS CUBE ---------------
  
sxaddpar,hdr,'DATAMIN',min(data_out,/nan)
sxaddpar,hdr,'DATAMAX',max(data_out,/nan)
sxaddpar,hdr,'HISTORY','Cube rebaselined by mc3_rebaseline.pro'
sxaddpar,hdr,'HISTORY','Baseline order used: '+strtrim(string(use_order),2)

badpix=where(blflags lt 0, nbad)
if nbad gt 0 then blflags[badpix]=!values.f_nan

if keyword_set(verbose) then print,blcount/nbltot," baselines modified"

if keyword_set(fits_out) then writefits,use_outdir+fits_out+'.fits',data_out,hdr

idl_out=data_out
hdr_out=hdr

the_end:
return

end
