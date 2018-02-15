function mc3_empty_sampling_struct

nan=!values.f_nan
s={x:nan,$
   y:nan,$
   ra:nan,$
   dec:nan,$
   scale_as:nan,$
   region_size:nan,$
   npix:nan,$
   ico_wgt:nan,$
   tpk_wgt:nan,$
   vdisp_wgt:nan,$
   ico_avg:nan,$
   tpk_avg:nan,$
   vdisp_avg:nan,$
   rms_val:nan,$
   idx:ptr_new(),$
   wgt_spec:ptr_new(),$
   avg_spec:ptr_new(),$
   wgt_spec_mom0:nan,$
   wgt_spec_pk:nan,$
   wgt_spec_eqw:nan,$
   wgt_spec_b:nan,$
   avg_spec_mom0:nan,$
   avg_spec_pk:nan,$
   avg_spec_eqw:nan, $
   avg_spec_b:nan $
  }

  return,s

end
