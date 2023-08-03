; No missing value are allowed in the input data
; Attention: we project the intraseasonal time-filtered OLR fields onto each mode 
;            in the entire period
;
begin

print("Reading the input file...")

f = addfile("./input.nc", "r")

dum 	= cd_calendar(f->time,0)
printVarSummary(dum)
print(dum(:3, 0) + "-" + dum(:3, 1) + "-" + dum(:3, 2))
;
i_yrRead   = (/1980, 2020/)
months 	= toint(dum(:,1))
imonths     = ind(months.ge.6 .and. months.le.10 .and. \
              dum(:,0) .ge. i_yrRead(0) .and. dum(:,0) .le. i_yrRead(1)) 
i_dimsSeaDays = 30 + 31 + 31 + 30 + 31

olr     = f->anom_bf
olr_entire = olr(:, {-30:30}, :)
olr    := olr(imonths, {-30:30}, :)
printVarSummary(olr)
;
d_time  = olr&time
printVarSummary(d_time)

print("Reshape...")

nyears	= dimsizes(olr&time)/i_dimsSeaDays   ; ;JJASO, which is 153 days
nlat	= dimsizes(olr&lat)
nlon	= dimsizes(olr&lon)


i_dims := (/nyears, i_dimsSeaDays, nlat, nlon/)
x       = reshape(olr, i_dims)
copy_VarMeta(olr(0, :, :), x(0, 0, :, :))
x!1     = "day"
x&day   = dum(:i_dimsSeaDays-1, 2)
x!0     = "year"
x&year  = ispan(i_yrRead(0), i_yrRead(1), 1)
printVarSummary(x)

d_time2 = reshape(d_time, i_dims(:1))
copy_VarAtts(d_time, d_time2)


print("stacking up in S-mode for different lags...")

lags    = (/-10,-5,0/)

data    = new((/nyears, i_dimsSeaDays+lags(0), 3*nlat, 3*nlon/),float) 
copy_VarMeta(x(:, 10:, 0, 0), data(:, :, 0, 0))
do i = 0, nyears-1
    data(i, :, 0:nlat-1, 0:nlon-1)           = x(i, 0:i_dimsSeaDays-1+lags(0), :, :) 
    data(i, :, nlat:2*nlat-1, nlon:2*nlon-1) = x(i, 5:i_dimsSeaDays-1+lags(1), :, :)
    data(i, :, 2*nlat:3*nlat-1, 2*nlon:3*nlon-1) = x(i, 10:i_dimsSeaDays-1+lags(2), :, :)
end do
r_x_original = x 
delete(x)

printVarSummary(data) ; [year | 41] x [day | 143] x [lat | 75] x [lon | 432]

x  = reshape(data, (/nyears*(i_dimsSeaDays + lags(0)), 3*nlat, 3*nlon/)) ; The size 75 is from 3*nlat
x!0="time"
x!1="lat"
x!2="lon"
copy_VarCoords(data(0, 0, :, :), x(0, :, :))
printMinMax(x, True)
printVarSummary(x) ; [time | 5863] x [lat | 75] x [lon | 432]


print("Calculating EEOFs in S-modes...")

neof        = 2
eof         = eofunc_n_Wrap(x, neof, False, 0)
printVarSummary(eof) ; [evn | 2] x [lat | 75] x [lon | 432]
printMinMax(eof,True)

;; PC series 

i_dimTime   = dimsizes(olr_entire(:, 0, 0))
x_entire    = new((/i_dimTime+lags(0), 3*nlat, 3*nlon/), float)
copy_VarCoords(data(0, 0, :, :), x_entire(0, :, :))
copy_VarMeta(olr_entire(10:, 0, 0), x_entire(:, 0, 0))
x_entire(:, 0:nlat-1, 0:nlon-1)               = olr_entire(0:i_dimTime-1+lags(0), :, :)
x_entire(:, nlat:2*nlat-1, nlon:2*nlon-1)     = olr_entire(5:i_dimTime-1+lags(1), :, :)
x_entire(:, 2*nlat:3*nlat-1, 2*nlon:3*nlon-1) = olr_entire(10:i_dimTime-1+lags(2), :, :)
ev_ts       = eofunc_ts_n_Wrap(x_entire, eof, False, 0)
printVarSummary(ev_ts) ; [evn | 2] x [time | ]
printMinMax(ev_ts, False)
;
; normalize
ev_ts       = dim_standardize_n_Wrap(ev_ts, 1, 1)
dum        := cd_calendar(ev_ts&time, 0)
print(dum(:50, 0) + "-" + dum(:50, 1) + "-" + dum(:50, 2) + "-")
;
;; reshape
r_eof       = new((/neof, 3, nlat, nlon/), float)
copy_VarCoords(r_x_original(0, 0, :, :), r_eof(0, 0, :, :))
copy_VarMeta(eof(:, 0, 0), r_eof(:, 0, 0, 0))
r_eof!1     = "lag_lead"
r_eof&lag_lead      = lags
r_eof(:, 0, :, :)   = eof(:, 0:nlat-1, 0:nlon-1)
r_eof(:, 1, :, :)   = eof(:, nlat:2*nlat-1, nlon:2*nlon-1)
r_eof(:, 2, :, :)   = eof(:, 2*nlat:3*nlat-1, 2*nlon:3*nlon-1)


;; check the sign of pattern 
;;; eof 1 : negative - positive dipole from Indian Ocean to Indian subcontinent
if ( any(r_eof(0, 0, {12:20}, {72:84}) .lt. 0) .or. \
     any(r_eof(0, 0, {-2.5:0.5}, {75:85}) .gt. 0) ) then
    r_eof(0, :, :, :) = r_eof(0, :, :, :) * (-1.)
    ev_ts(0, :)       = ev_ts(0, :) * (-1.)
    print("reverse the eof1!!!")
end if
if ( any(r_eof(1, 0, {-2.5:5}, {80:90}) .lt. 0) .or. \
     any(r_eof(1, 0, {12.5:17.5}, {128:132}).gt. 0) ) then
    r_eof(1, :, :, :) = r_eof(1, :, :, :) * (-1.)
    ev_ts(1, :)       = ev_ts(1, :) * (-1.)
    print("reverse the eof2!!!")
end if

;; output 
s_DirOut = "./"
system("mkdir -p " + s_DirOut)
s_FileOut           = s_DirOut + "cal_eeof_olr_inReanalysis.nc"
system("rm -rf " + s_FileOut)
FileOut             = addfile(s_FileOut, "c")
GAttr               = True
GAttr@history       = systemfunc("date") + "  " + systemfunc("pwd") + \
                    "/cal_eeof_olr.ncl"
fileattdef(FileOut, GAttr)
;
FileOut->eof        = r_eof
FileOut->ev_ts      = ev_ts
;
print(s_FileOut + "")
;
end








