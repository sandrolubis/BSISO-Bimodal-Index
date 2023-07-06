; No missing value are allowed in the input data

begin

print("Reading the input file...")

f = addfile("./input.nc", "r") ; complete year

dum 	= cd_calendar(f->time,0)
months 	= toint(dum(:,1))
imonths = ind(months.ge.6 .and. months.le.10) ;JJASO, which is 153 days


olr     = f->anom_bf
olr    := olr(imonths, {-30:30}, :)
printVarSummary(olr)
;
d_time  = olr&time
printVarSummary(d_time)
; exit()

print("Reshape...")

nyears	= dimsizes(olr&time)/153   ; ;JJASO, which is 153 days
nlat	= dimsizes(olr&lat)
nlon	= dimsizes(olr&lon)

; x	    = onedtond(olr, (/nyears, 153, nlat, nlon/))
;; Ziming on 06/05/2023
i_dims := (/nyears, 153, nlat, nlon/)
x       = reshape(olr, i_dims)
copy_VarMeta(olr(0, :, :), x(0, 0, :, :))
x!1     = "day"
x&day   = dum(:152, 2)
x!0     = "year"
x&year  = ispan(1980, 2020, 1)
printVarSummary(x)
; print(x&day(:30) + "")
d_time2 = reshape(d_time, i_dims(:1))
copy_VarAtts(d_time, d_time2)
; exit()

print("stacking up in S-mode for different lags...")

lags    = (/-10,-5,0/)
; data    = new((/nyears, 143, 3*nlat, nlon/),float) 
;; Ziming Chen on 06/05/2023
data    = new((/nyears, 143, 3*nlat, 3*nlon/),float) 
copy_VarMeta(x(:, 10:, 0, 0), data(:, :, 0, 0))
do i = 0, nyears-1
	; data(i, :, 0:nlat-1, :)     	= x(i, 0:152+lags(0), :, :) 
	; data(i, :, nlat:2*nlat-1, :)	= x(i, 5:152+lags(1), :, :)
	; data(i, :, 2*nlat:3*nlat-1, :)	= x(i, 10:152+lags(2), :, :)
    ;; Ziming Chen on 06/05/2023
    data(i, :, 0:nlat-1, 0:nlon-1)           = x(i, 0:152+lags(0), :, :) 
    data(i, :, nlat:2*nlat-1, nlon:2*nlon-1) = x(i, 5:152+lags(1), :, :)
    data(i, :, 2*nlat:3*nlat-1, 2*nlon:3*nlon-1) = x(i, 10:152+lags(2), :, :)
end do
r_x_original = x 
delete(x)

printVarSummary(data) ; [year | 41] x [day | 143] x [lat | 75] x [lon | 432]

x  = reshape(data, (/nyears*143, 3*nlat, 3*nlon/)) ; The size 75 is from 3*nlat
x!0="time"
x!1="lat"
x!2="lon"
copy_VarCoords(data(0, 0, :, :), x(0, :, :))
printMinMax(x, True)
printVarSummary(x) ; [time | 5863] x [lat | 75] x [lon | 432]
; exit()


print("Calculating EEOFs in S-modes...")

neof        = 2
eof         = eofunc_n_Wrap(x, neof, False, 0)
printVarSummary(eof) ; [evn | 2] x [lat | 75] x [lon | 432]
printMinMax(eof,True)
;
;; PC series 
ev_ts       = eofunc_ts_n_Wrap(x, eof, False, 0)
printVarSummary(ev_ts) ; [evn | 2] x [time | 5863]
printMinMax(ev_ts, False)
;
; normalize
ev_ts       = dim_standardize_n_Wrap(ev_ts, 1, 1)
printVarSummary(ev_ts)
d_time3     = new(dimsizes(ev_ts(0, :)), double)
copy_VarAtts(d_time, d_time3)
do i = 0, nyears-2
    printVarSummary(d_time2(i, 10:))
    printVarSummary(d_time3(i*143:(i+1)*143-1))
    d_time3(i*143:(i+1)*143-1) = d_time2(i, 10:)
end do
ev_ts&time := d_time3
dum        := cd_calendar(ev_ts&time, 0)
print(dum(0:10, 0) + "-" + dum(0:10, 1) + "-" + dum(0:10, 2) + "-")
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
;
; ;; reshape the PC seriess 
; r_ev_ts     = new((/2, nyears, 143/), float)
; copy_VarCoords(data(:, :, 0, 0), r_ev_ts(0, :, :))
; copy_VarMeta(ev_ts(:, 0), r_ev_ts(:, 0, 0))
; do i = 0, nyears-2
;     r_ev_ts(:, i, :)   = ev_ts(:, i*143:(i+1)*143-1 ) + 0
; end do
; r_ev_ts(:, nyears-1, :)= ev_ts(:, (nyears-1)*143: ) + 0
;
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








