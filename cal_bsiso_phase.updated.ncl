;; Ziming CHEN: 2023.07.11
; attention: we cal the bsiso phase in the entire period 
;
undef("my_atan")
function my_atan(y, x)
; in this code, we compute the atan by the ratio of (y/x)
local r_y, r_x, pi, r_theta
begin
;
r_y 				= y 
r_x 				= x 
pi 				= acos(-1.0)
;
r_x 				= where(r_x .eq. 0, r_x@_FillValue, r_x)
;
r_Value 			= r_y / r_x 
r_theta 			= atan(r_Value) ; range of returen values: -pi/2 < theta < pi/2
;
;; the quadrants 
;;; the first : y > 0, x > 0
;;; we do not need to do anything 
;
;;; the second: y > 0, x < 0, then need to be: pi/2 < theta < pi
r_theta 			= where(r_y .ge. 0 .and. r_x .lt. 0, r_theta + pi, r_theta)
;
;;; the third: y < 0, x < 0, then need to be: -pi < theta < -pi/2
r_theta 			= where(r_y .le. 0 .and. r_x .lt. 0, r_theta - pi, r_theta)
;
;;; the fourth: y < 0, and x > 0, then need to be: -pi/2 < theta < 0
;; we do not need to do anything 
;
return(r_theta)
;
end



begin
s_DirRead = "./"
f1 = addfile(s_DirRead + "cal_eeof_olr_inReanalysis.nc", "r")
system("ls -lh " + s_DirRead + "cal_eeof_olr_inReanalysis.nc")

x     = f1->ev_ts
printVarSummary(x)  ; pc

ntime=dimsizes(x&time)

dates=cd_calendar(x&time,-2)
print(dates(:10))


print("Calculating BSISO index...")

r = sqrt(x(evn|0, time|:) ^ 2.0 + x(evn|1, time|:) ^ 2.0)
r!0 = "time"
r&time = ispan (0, ntime-1,1)

pi = acos(-1.0)

;
;; Ziming CHEN: 2023.07.07
; printVarSummary(x) ; [evn | 2] x [time | 5863]
theta 	   = my_atan(-1. * x(0, :), -1. * x(1, :))

theta 	   = theta * 180.0 / pi
theta!0    = "time"
theta&time = ispan (0, ntime-1,1)

printVarSummary(r)
printVarSummary(theta)
printMinMax(theta, False)


print("Indexing...")


global_phase = new(ntime, integer)

;
do i = 0, ntime - 1
	; if (r(i).le.1) then
	; 	global_phase(i) = 0  ; inactive phase
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
;
;; nAmp greater than 1
global_phase_nAmpGe1 = new(ntime, integer)
;
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


printVarSummary(global_phase)
print(global_phase(:80) + "")
printMinMax(global_phase,True)
; exit()

bsiso_phase      			=  global_phase
bsiso_phase!0 				= "time"
bsiso_phase&time 			= x&time
global_phase_nAmpGe1!0  	= "time"
global_phase_nAmpGe1&time   = x&time
printMinMax(bsiso_phase,True)
printVarSummary(bsiso_phase)


print("Opening the output file...")

fout1 = "bsiso_phase.nc"

system("rm " + fout1)
g1 = addfile(fout1, "c")
;; Ziming Chen: 2023.07.09
GAttr 			= True
GAttr@history 	= systemfunc("date") + "  " + systemfunc("pwd") + \
					"/cal_bsiso_phase.updated.ncl"
fileattdef(g1, GAttr)
g1->bsiso_phase = bsiso_phase
g1->bsiso_phase_nAmpGe1 = global_phase_nAmpGe1

print(fout1 + "")
print("Closing the output file...")
end




