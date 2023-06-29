
   var     = "olr"                            ; name of file                   
   diri    = "./"                             ; new input directory
   fili    = "olr.anom.nc"                    ; input OLR anomaly       


;***********************************************************
; Read user specified time and create required yyyyddd                    
;***********************************************************

   f       = addfile (diri+fili , "r")                          
   anom       = f->$var$(:,:,:) 
   printVarSummary( anom ) 

; ***********************************************
; create the filter weights and apply
; ***********************************************

  ihp     = 2                           ; band pass
  sigma   = 1.0                           ; Lanczos sigma

  nWgt    = 141                                                    
  fca     = 1./90.                        ; start freq
  fcb     = 1./25.                        ; last  freq
  wgt     = filwgts_lanczos (nWgt, ihp, fca, fcb, sigma )
  printVarSummary(wgt)
  BPF    = wgt_runave_n_Wrap ( anom, wgt, 0,0 )      ; 25-90 day
  printVarSummary(BPF)

;***********************************************************
; Save the output in netcdf
;***********************************************************

setfileoption("nc","Format","LargeFile")
system( "rm " + "./input.nc" )
output = addfile("./input.nc","c")
output->anom_bf = BPF



