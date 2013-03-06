function util.computeKineticConvRatesForSpace(dom, maxLev, minLev, discType, p, 
													timeScheme, orderOrTheta, startTime, endTime,
													dtmin, dtmax, dtred, exactSol, exactGrad,
													createApproxSpace, createDomainDisc, createSolver,
									 				plotPath, solPath, dataPath, bLinear)

print("\n")
print("---------------------------")
print(" General parameters chosen:")
print("---------------------------")
print("    dim        = " .. dim)
print("    grid       = " .. gridName)
print("    maxLev     = " .. maxLev)
print("    minLev     = " .. minLev)
print("    dtmin      = " .. dtmin)
print("    dtmax      = " .. dtmax)
print("    dtred      = " .. dtred)
print("    timeScheme = " .. timeScheme)
print("    order      = " .. orderOrTheta)
print("    startTime  = " .. startTime)
print("    endTime    = " .. endTime)
print("    type       = " .. discType)
print("    order      = " .. p)
print("\n")

-- create Approximation Space
print(">> Create ApproximationSpace: "..discType..", "..p)
local approxSpace = createApproxSpace(dom, discType, p)

print(">> Create Domain Disc: "..discType..", "..p)
local domainDisc = createDomainDisc(discType, p, approxSpace)

print(">> Create Solver")
local solver = createSolver(approxSpace)

--------------------------------------------------------------------------------
--  Apply Solver
--------------------------------------------------------------------------------

local l2exact = {}
local l2diff = {}
local h1exact = {}
local h1diff = {}
local numDoFs = {}
local dts = {}

local kmax = 0
local dt = dtmax
while dt >= dtmin do
	kmax = kmax + 1
	dt = dt*dtred
end

local uMostAccurate = nil;
for lev = maxLev, minLev, -1 do
	l2exact[lev] = {}
	l2diff[lev] = {}
	h1exact[lev] = {}
	h1diff[lev] = {}

	for k = kmax, 1, -1 do
		local dt = dtmax * (dtred^(k-1))
		dts[k] = dt

		u = GridFunction(approxSpace, lev)
		Interpolate(exactSol, u, "c", startTime)
		write("\n>> Algebra created.\n")
		
		if bLinear == true then
			util.SolveLinearTimeProblem(u, domainDisc, solver, nil, nil,
										timeScheme, orderOrTheta, startTime, endTime, dt);
		else
			util.SolveNonlinearTimeProblem(u, domainDisc, solver, nil, nil,
											timeScheme, orderOrTheta, startTime, endTime, dt);
		end
		
		if lev == maxLev and k == kmax then
			uMostAccurate = u:clone();
		end			
		
		-- compute error
		quadOrder = p+2
		l2exact[lev][k] = L2Error(exactSol, u, "c", endTime, quadOrder)
		l2diff[lev][k] = L2Error(uMostAccurate, "c", u, "c", quadOrder)
		h1exact[lev][k] = H1Error(exactSol, exactGrad, u, "c", endTime, quadOrder)
		h1diff[lev][k] = H1Error(uMostAccurate, "c", u, "c", quadOrder)
		numDoFs[lev] = u:size()
		write(">> L2-Error on Level "..lev..", dt: "..dt.." is "..string.format("%.3e", l2exact[lev][k]) .."\n");
		write(">> L2-Diff  on Level "..lev..", dt: "..dt.." is "..string.format("%.3e", l2diff[lev][k]) .."\n");
		write(">> H1-Error on Level "..lev..", dt: "..dt.." is "..string.format("%.3e", h1exact[lev][k]) .."\n");
		write(">> H1-Diff  on Level "..lev..", dt: "..dt.." is "..string.format("%.3e", h1diff[lev][k]) .."\n");
		write(">> #DoF     on Level "..lev.." is "..numDoFs[lev] .."\n");
	end	
end

return l2exact, l2diff, h1exact, h1diff, numDoFs, dts
end


function util.computeRate(data, min, max, base, style, levOrK)

	local fac = {}
	local rate = {}
	local cpyData = {}

	if style == "L" then
		for lev = min, max do
			cpyData[lev] = data[lev][levOrK]
		end
	elseif style == "dt" then
		for k = min, max do
			cpyData[k] = data[levOrK][k]
		end
	elseif style == "diag" then
		local size = math.min(max-min, #data[min])
		for k = 1, size do
			cpyData[k] = data[min+k-1][k]
		end
		min = 1
		max = size
	else
		print("Style not foung. Use 'L', 'dt' or 'diag'."); exit();
	end

	rate[min] = "--"
	for k = min+1, max do	
		fac[k] = cpyData[k-1]/cpyData[k]
		rate[k] = math.log(fac[k]) / math.log(base)
		if math.abs(rate[k]) < 1e-10 then rate[k] = "--" end
	end
			
	return cpyData, rate
end

function util.computeKineticConvRates(dom, createApproxSpace, createDomainDisc, createSolver, DiscTypes, bLinear)
	
	-- create directories
	plotPath = "plots/"
	solPath  = "sol/"
	dataPath = "data/"
	os.execute("mkdir " .. dataPath)
	os.execute("mkdir " .. plotPath)
	os.execute("mkdir " .. solPath)

	-- compute element size	
	local numRefs = dom:grid():num_levels() - 1;
	local exactSol		= DiscTypes.exactSol
	local exactGrad		= DiscTypes.exactGrad

	-- compute problem
	for type = 1,#DiscTypes do
		local l2_plotdata = {}; local h1_plotdata = {}
		local l2_h_plotdata = {}; local h1_h_plotdata = {}
		local l2_dt_plotdata = {}; local h1_dt_plotdata = {}
		local discType 		= DiscTypes[type].type
		local pmin 			= DiscTypes[type].pmin
		local pmax 			= DiscTypes[type].pmax
		local lmin 			= DiscTypes[type].lmin
		local lmax 			= DiscTypes[type].lmax
		
		local dtmin 		= DiscTypes[type].dtmin		
		local dtmax 		= DiscTypes[type].dtmax		
		local dtred 		= DiscTypes[type].dtred		
		local timeScheme	= DiscTypes[type].timeScheme
		local orderOrTheta 	= DiscTypes[type].orderOrTheta
		local startTime 	= DiscTypes[type].startTime
		local endTime		= DiscTypes[type].endTime
				
		if discType == nil then print("discType required."); exit(); end
		if pmin == nil then pmin = 1 end
		if pmax == nil then pmax = 1 end
		if lmin == nil then lmin = 0 end
		if lmax == nil then lmax = numRefs end
		
		if lmin > lmax then
			print("lmin: "..lmin.." must be less or equal lmax: "..lmax)
			exit()
		end
		if lmax > numRefs then
			print("lmax: "..lmax.." must be less or equal numRefs: "..numRefs)
			exit()
		end

		if dtmin == nil then print("dtmin required."); exit(); end
		if dtmax == nil then print("dtmax required."); exit(); end
		if dtmin > dtmax then
			print("dtmin: "..dtmin.." must be less or equal dtmax: "..dtmax)
			exit()
		end
		if dtred == nil then dtred = 0.5 end
		if dtred >= 1 then
			print("dtred: "..dtred.." must be less than 1")
			exit()
		end
		if timeScheme == nil then print("timeScheme required."); exit(); end
		if orderOrTheta == nil then orderOrTheta = 1; end
		if startTime == nil then print("startTime required."); exit(); end
		if endTime == nil then print("endTime required."); exit(); end
		if exactSol == nil then print("exactSol required."); exit(); end
		
		for p = pmin, pmax do
		
			local maxLev = lmax - math.floor(p/2)
			if pmin == pmax then maxLev = lmax end
			local minLev = lmin

			local h = {}
			for lev = minLev, maxLev do
				h[lev] = MaxElementDiameter(dom, lev)
			end
			
			local l2exact, l2diff, h1exact, h1diff, numDoFs, dts
				 = util.computeKineticConvRatesForSpace(
				 		dom, maxLev, minLev, discType, p, 
				 		timeScheme, orderOrTheta, startTime, endTime,
				 		dtmin, dtmax, dtred, exactSol, exactGrad,
				 		createApproxSpace, createDomainDisc, createSolver,
				 		plotPath, solPath, dataPath, bLinear)
				 		
			-- write: for each time step size show convergence in space, keeping time step constant
			for k = 1, #dts do
				local l2ex, l2exactrate = util.computeRate(l2exact, minLev, maxLev, 2, "L", k)
				local h1ex, h1exactrate = util.computeRate(h1exact, minLev, maxLev, 2, "L", k)
				local l2di, l2diffrate =  util.computeRate(l2diff,  minLev, maxLev, 2, "L", k)
				local h1di, h1diffrate =  util.computeRate(h1diff,  minLev, maxLev, 2, "L", k)
				local level = {}; for lev = minLev, maxLev do level[lev] = lev end
				local dtout = {}; for lev = minLev, maxLev do dtout[lev] = dts[k] end
				
				-- write data to screen
				print("\nConvergence in space: fixed dt = "..dts[k])
				table.print({level, h, dtout, numDoFs, l2di, l2diffrate, l2ex, l2exactrate, h1di, h1diffrate, h1ex, h1exactrate}, 
							{title = {"L", "h", "dt", "#DoFs", "L2 Diff", "rate", "L2 Exact", "rate", "H1 Diff", "rate", "H1 Exact", "rate"}, 
							 format = {"%d", "%.2e", "%.2e", "%d", "%.2e", "%.3f", "%.2e", "%.3f", "%.2e", "%.3f", "%.2e", "%.3f"},
							 hline = true, vline = true})
				
				-- write data to file
				gnuplot.write_data(dataPath.."error_in_h__constdt"..k.."__"..discType..p..".dat", {numDoFs, l2ex, h1ex})
				
				-- schedule l2 plot
				l2_h_plotdata = gnuplot.array_concat(l2_h_plotdata,
				{{label="L_2 error P_"..p.." dt: "..dts[k], file=dataPath.."error_in_h__constdt"..k.."__"..discType..p..".dat", style="linespoints", 1, 2 }})
				
				-- schedule h1 plot
				h1_h_plotdata = gnuplot.array_concat(h1_h_plotdata, 
				{{label="H^1 error P_"..p.." dt: "..dts[k], file=dataPath.."error_in_h__constdt"..k.."__"..discType..p..".dat", style="linespoints", 1, 3 }})	
			end

			-- write: for each grid size show convergence in time, keeping grid size constant
			for lev = minLev, maxLev do
				local l2ex, l2exactrate = util.computeRate(l2exact, 1, #dts, 1/dtred, "dt", lev)
				local h1ex, h1exactrate = util.computeRate(h1exact, 1, #dts, 1/dtred, "dt", lev)
				local l2di, l2diffrate =  util.computeRate(l2diff,  1, #dts, 1/dtred, "dt", lev)
				local h1di, h1diffrate =  util.computeRate(h1diff,  1, #dts, 1/dtred, "dt", lev)
				local numDoFsout = {}; for k = 1, #dts do numDoFsout[k] = numDoFs[lev] end
				local level = {}; for k = 1, #dts do level[k] = lev end
				local hout  = {}; for k = 1, #dts do hout[k] = h[lev] end
				local dtout = {}; for k = 1, #dts do dtout[k] = dts[k] end
				
				-- write data to screen
				print("\nConvergence in time: fixed h = "..string.format("%.3e", h[lev])..", #DoF = "..numDoFs[lev])
				table.print({level, hout, dtout, numDoFsout, l2di, l2diffrate, l2ex, l2exactrate, h1di, h1diffrate, h1ex, h1exactrate}, 
							{title = {"L", "h", "dt", "#DoFs", "L2 Diff", "rate", "L2 Exact", "rate", "H1 Diff", "rate", "H1 Exact", "rate"}, 
							 format = {"%d", "%.2e", "%.2e", "%d", "%.2e", "%.3f", "%.2e", "%.3f", "%.2e", "%.3f", "%.2e", "%.3f"},
							 hline = true, vline = true})
				
				-- write data to file
				gnuplot.write_data(dataPath.."error_in_dt__constdof"..numDoFs[lev].."__"..discType..p..".dat", {dtout, l2ex, h1ex})
				
				-- schedule l2 plot
				l2_dt_plotdata = gnuplot.array_concat(l2_dt_plotdata,
				{{label="L_2 error P_"..p.." #dof: "..numDoFs[lev], file=dataPath.."error_in_dt__constdof"..numDoFs[lev].."__"..discType..p..".dat", style="linespoints", 1, 2 }})
				
				-- schedule h1 plot
				h1_dt_plotdata = gnuplot.array_concat(h1_dt_plotdata, 
				{{label="H^1 error P_"..p.." #dof: "..numDoFs[lev], file=dataPath.."error_in_dt__constdof"..numDoFs[lev].."__"..discType..p..".dat", style="linespoints", 1, 3 }})	
			end

			-- write data to file
			gnuplot.write_data(dataPath.."error__"..discType..p..".dat", {numDoFs, dts, l2exact, h1exact})
			
			-- schedule l2 plot
			l2_plotdata = gnuplot.array_concat(l2_plotdata,
			{{label="L_2 error P_"..p, file=dataPath.."error__"..discType..p..".dat", style="linespoints", 1, 2, 3}})
			
			-- schedule h1 plot
			h1_plotdata = gnuplot.array_concat(h1_plotdata, 
			{{label="H^1 error P_"..p, file=dataPath.."error__"..discType..p..".dat", style="linespoints", 1, 2, 4}})	
		end
		
		-- create plots
		gnuplot.plot(plotPath..discType.."_h_l2.pdf", l2_h_plotdata, {grid = true, logscale = true, xlabel = "#DoF"})
		gnuplot.plot(plotPath..discType.."_h_h1.pdf", h1_h_plotdata, {grid = true, logscale = true, xlabel = "#DoF"})

		gnuplot.plot(plotPath..discType.."_dt_l2.pdf", l2_dt_plotdata, {grid = true, logscale = true, xlabel = "dt"})
		gnuplot.plot(plotPath..discType.."_dt_h1.pdf", h1_dt_plotdata, {grid = true, logscale = true, xlabel = "dt"})
		
		gnuplot.plot(plotPath..discType.."_l2.pdf", l2_plotdata, {grid = true, logscale = true, xlabel = "#DoF", ylabel = "dt"})
		gnuplot.plot(plotPath..discType.."_h1.pdf", h1_plotdata, {grid = true, logscale = true, xlabel = "#DoF", ylabel = "dt"})
		
	end
end

function util.computeLinearKineticConvRates(dom, createApproxSpace, createDomainDisc, createSolver, DiscTypes)
	util.computeKineticConvRates(dom, createApproxSpace, createDomainDisc, createSolver, DiscTypes, true)
end

function util.computeNonlinearKineticConvRates(dom, createApproxSpace, createDomainDisc, createSolver, DiscTypes)
	util.computeKineticConvRates(dom, createApproxSpace, createDomainDisc, createSolver, DiscTypes, false)
end