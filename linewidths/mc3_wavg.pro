function mc3_wavg,array,worerr,weights=weights,errors=errors,forcepos=forcepos,help=help,verbose=verbose
;+
; quicky weighted average
;
; results=wavg(array,weights,/weights)
; or
; results=wavg(array,errors,/errors)
;
;-

  if (n_params(0) lt 2) or keyword_set(help) then begin
    print,'Call> results=wavg(array,weights,/weights)'
    print,' or'
    print,'Call> results=wavg(array,errors,/errors)'
    return,-1
    endif

  if (n_elements(weights) eq 0) then weights=0
  if (n_elements(errors) eq 0) then errors=0

  if (weights+errors eq 0) then begin
    message,'WAVG: Error.  either /weights or /errors must be set',/info
    return,-1
    endif

  if (n_elements(array) ne n_elements(worerr)) then begin
    message,'WAVG: Error.  array and errors must have same dimensions',/info
    return,-1
    endif

  if keyword_set(forcepos) then begin
     if keyword_set(verbose) then message,'Setting negative values and their weights to NaN',/info
     negidx = where (array lt 0, negct)
     if (negct gt 0) then begin
        array[negidx]=!values.f_nan   
        worerr[negidx]=!values.f_nan
        end
  end

  infidx=where(finite(array) eq 0 or finite(worerr) eq 0,infct)
  if infct gt 0 then begin
     array[infidx] = !values.f_nan & worerr[infidx]=!values.f_nan
  end

  if (errors ne 0) then begin
    npts=n_elements(worerr)
    errorarr=worerr
    weightarr=1/errorarr^2
    wav1=total(array*weightarr,/nan)/total(weightarr,/nan)
    wer1=sqrt(1/total(1/worerr^2),/nan)
    return,[wav1,wer1]

;    wer1=sqrt( (total(array^2*weightarr)/total(weightarr) - wav1^2) * $
;              (npts/(npts-1)) ) / sqrt(npts)
;      1/sqrt(total(1/worerr^2))]
;      sqrt(total(worerr^2*weightarr))/total(weightarr)] ; fixed 980312
    endif

  weightarr=worerr

  return,total(array*weightarr,/nan)/total(weightarr,/nan)

end

  
