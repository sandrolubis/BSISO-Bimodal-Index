undef("my_atan")
function my_atan(y, x)
; In this code, we compute atan using the ratio (y/x)
local r_y, r_x, pi, r_theta, r_Value
begin
    r_y = y
    r_x = x
    pi  = acos(-1.0)

    r_x = where(r_x .eq. 0, r_x@_FillValue, r_x)

    r_Value = r_y / r_x
    r_theta = atan(r_Value)

    ;;; Quadrants

    ; First: y > 0, x > 0
    ; No modification needed

    ; Second: y > 0, x < 0
    r_theta = where(r_y .ge. 0 .and. r_x .lt. 0, r_theta + pi, r_theta)

    ; Third: y < 0, x < 0
    r_theta = where(r_y .le. 0 .and. r_x .lt. 0, r_theta - pi, r_theta)

    ; Fourth: y < 0, x > 0
    ; No modification needed

    return(r_theta)
end


;;; Start from here! ;;;;

begin

s_DirRead = "./"

f1 = addfile(s_DirRead + "cal_eeof_olr_inReanalysis.nc", "r")
system("ls -lh " + s_DirRead + "cal_eeof_olr_inReanalysis.nc")

x = f1->ev_ts
printVarSummary(x)

ntime = dimsizes(x&time)

dates = cd_calendar(x&time, -2)
print(dates(:10))


print("Calculating BSISO index...")


;--------------------------------------------------
; Calculate BSISO amplitude
;--------------------------------------------------

r = sqrt(x(evn|0, time|:) ^ 2.0 + x(evn|1, time|:) ^ 2.0)

r!0 = "time"
r&time = x&time

pi = acos(-1.0)


;--------------------------------------------------
; Extract PC1, PC2, and amplitude
;--------------------------------------------------

PC1 = x(0,:)
PC2 = x(1,:)
amp = r

PC1!0 = "time"
PC2!0 = "time"
amp!0 = "time"

PC1&time = x&time
PC2&time = x&time
amp&time = x&time

printVarSummary(PC1)
printVarSummary(PC2)
printVarSummary(amp)


;--------------------------------------------------
; Calculate BSISO phase angle
;--------------------------------------------------

theta = my_atan(-1. * x(0,:), -1. * x(1,:))
theta = theta * 180.0 / pi

theta!0 = "time"
theta&time = x&time

printVarSummary(r)
printVarSummary(theta)
printMinMax(theta, False)


print("Indexing...")


;--------------------------------------------------
; BSISO phase without amplitude threshold
;--------------------------------------------------

global_phase = new(ntime, integer)

do i = 0, ntime - 1

    ; if (r(i).le.1) then
    ;     global_phase(i) = 0
    ; else

        if (theta(i).ge.-180 .and. theta(i).lt.-135) then
            global_phase(i) = 1
        end if

        if (theta(i).ge.-135 .and. theta(i).lt.-90) then
            global_phase(i) = 2
        end if

        if (theta(i).ge.-90 .and. theta(i).lt.-45) then
            global_phase(i) = 3
        end if

        if (theta(i).ge.-45 .and. theta(i).lt.0) then
            global_phase(i) = 4
        end if

        if (theta(i).ge.0 .and. theta(i).lt.45) then
            global_phase(i) = 5
        end if

        if (theta(i).ge.45 .and. theta(i).lt.90) then
            global_phase(i) = 6
        end if

        if (theta(i).ge.90 .and. theta(i).lt.135) then
            global_phase(i) = 7
        end if

        if (theta(i).ge.135 .and. theta(i).lt.180) then
            global_phase(i) = 8
        end if

    ; end if

end do


;--------------------------------------------------
; BSISO phase with amplitude > 1
;--------------------------------------------------

global_phase_nAmpGe1 = new(ntime, integer)

do i = 0, ntime - 1

    if (r(i).le.1) then

        global_phase_nAmpGe1(i) = 0  ; inactive phase

    else

        if (theta(i).ge.-180 .and. theta(i).lt.-135) then
            global_phase_nAmpGe1(i) = 1
        end if

        if (theta(i).ge.-135 .and. theta(i).lt.-90) then
            global_phase_nAmpGe1(i) = 2
        end if

        if (theta(i).ge.-90 .and. theta(i).lt.-45) then
            global_phase_nAmpGe1(i) = 3
        end if

        if (theta(i).ge.-45 .and. theta(i).lt.0) then
            global_phase_nAmpGe1(i) = 4
        end if

        if (theta(i).ge.0 .and. theta(i).lt.45) then
            global_phase_nAmpGe1(i) = 5
        end if

        if (theta(i).ge.45 .and. theta(i).lt.90) then
            global_phase_nAmpGe1(i) = 6
        end if

        if (theta(i).ge.90 .and. theta(i).lt.135) then
            global_phase_nAmpGe1(i) = 7
        end if

        if (theta(i).ge.135 .and. theta(i).lt.180) then
            global_phase_nAmpGe1(i) = 8
        end if

    end if

end do


;--------------------------------------------------
; Assign time coordinates
;--------------------------------------------------

printVarSummary(global_phase)
print(global_phase(:80) + "")
printMinMax(global_phase, True)

bsiso_phase = global_phase

bsiso_phase!0 = "time"
bsiso_phase&time = x&time

global_phase_nAmpGe1!0 = "time"
global_phase_nAmpGe1&time = x&time

printMinMax(bsiso_phase, True)
printVarSummary(bsiso_phase)


;--------------------------------------------------
; Write output NetCDF
;--------------------------------------------------

print("Opening the output file...")

fout1 = "bsiso_phase.nc"

system("rm -f " + fout1)

g1 = addfile(fout1, "c")

GAttr = True
GAttr@history = systemfunc("date") + "  " + systemfunc("pwd") + \
                "/cal_bsiso_phase.updated.ncl"

fileattdef(g1, GAttr)


;--------------------------------------------------
; Output BSISO phases
;--------------------------------------------------

g1->bsiso_phase = bsiso_phase
g1->bsiso_phase_nAmpGe1 = global_phase_nAmpGe1


;--------------------------------------------------
; Output PC1, PC2, and amplitude
;--------------------------------------------------

g1->PC1 = PC1
g1->PC2 = PC2
g1->amp = amp


print(fout1 + "")
print("Closing the output file...")

end
