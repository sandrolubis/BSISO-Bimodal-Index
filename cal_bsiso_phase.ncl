
s_DirRead = "./"
f1 = addfile(s_DirRead + "cal_eeof_olr_inReanalysis.nc", "r")

system("ls -lh " + s_DirRead + "cal_eeof_olr_inReanalysis.nc")
; exit()

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

theta = atan2(-1 * x(0, :), -1 * x(1, :))

theta = theta * 180.0 / pi
theta!0 = "time"
theta&time = ispan (0, ntime-1,1)

printVarSummary(r)
printVarSummary(theta)


print("Indexing...")


global_phase = new(ntime, integer)

; ;; Ziming CHEN: 2023.07.07
do i = 0, ntime - 1
	if (r(i).le.1) then
		global_phase(i) = 0  ; inactive phase
	else
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
	end if
end do


printVarSummary(global_phase)
print(global_phase(:50) + "")
printMinMax(global_phase,True)
; exit()

bsiso_phase     = global_phase
bsiso_phase!0="time"
bsiso_phase&time=x&time
printMinMax(bsiso_phase,True)
printVarSummary(bsiso_phase)


print("Opening the output file...")

fout1 = "bsiso_phase.nc"

system("rm " + fout1)
g1 = addfile(fout1, "c")
;; Ziming Chen: 2023.07.09
GAttr 			= True
GAttr@history 	= systemfunc("date") + "  " + systemfunc("pwd") + \
					"/cal_bsiso_phase.ncl"
fileattdef(g1, GAttr)
g1->bsiso_phase = bsiso_phase

print("Closing the output file...")





