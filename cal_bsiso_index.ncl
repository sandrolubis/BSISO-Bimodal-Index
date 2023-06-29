
f1 = addfile("EEOF_OLR.nc", "r")

pc     = f1->ev_ts
printVarSummary(pc)


nevn 	= dimsizes(pc&evn)
nyear	= dimsizes(pc&year)
nday 	= dimsizes(pc&day)
ntime   = nyear*nday

x     = onedtond(pc,(/nevn,nyear*nday/))
x!0="evn"
x!1="time"
x&evn=pc&evn
printVarSummary(x)


print("Calculating BSISO index...")

r = sqrt(x(evn|0, time|:) ^ 2.0 + x(evn|1, time|:) ^ 2.0)
r!0 = "time"
r&time = ispan (0, ntime-1,1)

pi = acos(-1.0)

theta = atan2(x(evn|1, time|:), x(evn|0, time|:))

theta = theta * 180.0 / pi
theta!0 = "time"
theta&time = ispan (0, ntime-1,1)

printVarSummary(r)
printVarSummary(theta)


print("Indexing...")


global_phase = new(ntime, integer)

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
printMinMax(global_phase,True)

bsiso_phase     = onedtond(global_phase,(/nyear, nday/))
bsiso_phase!0="year"
bsiso_phase&year=pc&year
bsiso_phase!1="day"
bsiso_phase&day=pc&day
printMinMax(bsiso_phase,True)
printVarSummary(bsiso_phase)


print("Opening the output file...")

fout1 = "bsiso_phase.nc"

system("rm " + fout1)
g1 = addfile(fout1, "c")
g1->bsiso_phase = bsiso_phase

print("Closing the output file...")




