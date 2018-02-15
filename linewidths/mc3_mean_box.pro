function mc3_mean_box,im,head,X0,Y0,dx_from,dx_to,st_mean,Nmean, $
                 silent=silent,dist=dist,d_value=d_value,indices=indices

;+
; NAME:
;       mean_box
; CALLING SEQUENCE:
;       mean=mean_box(im,head,X0,Y0,dx_from,dx_to,st_mean,Nmean[,/silent][,dist=dist])
; PURPOSE:
;	Compute the average and stdev in an image into an annulus centered on pixel X0,Y0
; INPUTS:
;	im,head	= image and header.
;	X0,Y0	= center pixel coordinates.
;	dx_from	= 1/2 interior size of the window.
;	dx_to	= 1/2 exterior size of the window
; OPTIONAL INPUT:
;	dist	= If set, boxes are choosed along this distance
; OUTPUTS:
;       Result of the function	= average in the box
;       stmean	= mean standard deviation
;	Nmean   = Number of defined pixels in averaged region
; SUBROUTINE AND PROCEDURE USED
;
;	stdev
; SIDE EFFECTS:
; MODIFICATION HISTORY:
;       WRITTEN, J.P. Bernard, 1-Aug-1997
;-

IF n_params() LT 8 then begin
  print,'truc=mean_box(image,header,X0,Y0,dx_from,dx_to,st_mean,Nmean)'
  print,'Accepted Key-Words: /silent,dist='
  RETURN,1
  GOTO,bad_sortie
ENDIF

si=size(im)
Nx=si(1)
Ny=si(2)

IF not keyword_set(dist) THEN BEGIN
  dist=fltarr(Nx,Ny)
  FOR i=0,Nx-1 DO BEGIN
    FOR j=0,Ny-1 DO BEGIN
      dist(i,j)=max(abs([i-x0,j-y0])) ; box
      dist(i,j)=sqrt((i-x0)^2+(j-y0)^2) ; circle?
    ENDFOR
  ENDFOR
ENDIF

;stop
  
flux_zone=[!values.f_nan]
dist_zone=[!values.f_nan]
ind=where(dist GE dx_from AND dist LT dx_to,count)
IF count NE 0 THEN BEGIN
  flux_zone=im(ind)
  dist_zone=dist(ind)
  indices=ind
ENDIF

nfz=count ; number of pixels in the flux zone
ind=where(finite(flux_zone) eq 1,count)
npix=count ; number of defined pixels

IF (count EQ 0) THEN BEGIN
  if not keyword_set(silent) THEN print,'integration region has only undefined values'
  goto,bad_sortie
ENDIF ELSE BEGIN
  if not keyword_set(silent) THEN BEGIN
    print,'Integration region has',npix,' pixels.'
    print,'Integration region has',nfz-npix,' undefined pixels.'
  ENDIF
ENDELSE
;when Moment can be used
IF npix GT 1 THEN BEGIN
  stdd=moment(flux_zone(ind))
  mmm=stdd(0) & std=sqrt(stdd(1))
;  std=stdev(flux_zone(ind),mmm)
  mean=mmm & st_mean=std 
  med=median(flux_zone(ind))
  iq=percentile(flux_zone(ind),25)-percentile(flux_zone(ind),75)
ENDIF ELSE BEGIN
  st_mean=!values.f_nan
  mean=avg(flux_zone(ind))
ENDELSE
d_value=avg(dist_zone(ind))

;Subtract background
;sum=flux-(2*dx+1)*(2*dy+1)*back
if n_elements(silent) EQ 0 THEN BEGIN
  print,'Mean= ',mean,' +- ',st_mean
  print,'Median= ',med,' +- ',iq
  ENDIF

Nmean=npix

return,[mean,st_mean,med,iq,npix]

bad_sortie:
st_mean=!values.f_nan
Nmean=0
return,[!values.f_nan,!values.f_nan,!values.f_nan,!values.f_nan,0]

end

