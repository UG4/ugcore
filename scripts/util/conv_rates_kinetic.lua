util = util or {}
util.rates = util.rates or {}
util.rates.kinetic = util.rates.kinetic or {}

ug_load_script("util/persistence.lua")

--------------------------------------------------------------------------------
-- Std Functions (used as defaults)
--------------------------------------------------------------------------------

function util.rates.kinetic.StdMaxLevelPadding(p)
	return math.floor(p/2)
end

--------------------------------------------------------------------------------
-- Label names
--------------------------------------------------------------------------------
	
util.rates.kinetic.StdLabel = util.rates.kinetic.StdLabel or {}
	
function util.rates.kinetic.StdLabel.MeasLatexP(disc, p)
	return disc.." $\\mathbb{P}_{"..p.."}$"
end

function util.rates.kinetic.StdLabel.MeasLatexQ(disc, p)
	return disc.." $\\mathbb{Q}_{"..p.."}$"
end

function util.rates.kinetic.StdLabel.XLatex(x)
	local gpXLabel ={ DoFs = "Anzahl Unbekannte",	h = "h (Gitterweite)"}
	return gpXLabel[x]
end

function util.rates.kinetic.StdLabel.YLatex(f, t, n)
	local gpType = {	["l-exact"] = 	"{}",		
						["l-lmax"] = 	"h_{\text{min}}",
						["l-prev"] = 	"{h/2}",
					}
	local gpNorm = 	{ l2 = "L_2",	h1 = "H^1"}
	
	if t == "interpol" then
		return "$\\norm{\\mathcal{I}_h("..f..") - "..f.."}_{"..gpNorm[n].."}$"
	else
		return "$\\norm{"..f.."_h - "..f.."_{"..gpType[t].."} }_{ "..gpNorm[n].."}$"
	end
end

function util.rates.kinetic.StdLabel.MeasPdfP(disc, p)
	return disc.." $P_"..p.."$"
end

function util.rates.kinetic.StdLabel.MeasPdfQ(disc, p)
	return disc.." $Q_"..p.."$"
end

function util.rates.kinetic.StdLabel.XPdf(x)
	local gpXLabel ={ DoFs = "Anzahl Unbekannte",	h = "h (Gitterweite)"}
	return gpXLabel[x]
end

function util.rates.kinetic.StdLabel.YPdf(f, t, n)
	local gpType = {	["l-exact"] = 	"{}",		
						["l-lmax"] = 	"h_{\text{min}}",
						["l-prev"] = 	"{h/2}",
					}
	local gpNorm = 	{ l2 = "L_2",	h1 = "H^1"}
	
	if t == "interpol" then
		return "$|| I_h("..f..") - "..f.." ||_{"..gpNorm[n].."}$"
	else
		return "$|| "..f.."_h - "..f.."_{"..gpType[t].."} ||_{"..gpNorm[n].."}$"
	end
end

--------------------------------------------------------------------------------
-- util.rates.kinetic.compute (main-function)
--------------------------------------------------------------------------------

function util.rates.kinetic.compute(ConvRateSetup)

	-- check passed param
	local CRS
	if ConvRateSetup == nil then print("No setup passed."); exit()		
	else CRS = ConvRateSetup end
	
	-- create directories
	local plotPath = CRS.plotPath or "plots/"
	local solPath  = CRS.solPath  or "sol/"
	local dataPath = CRS.dataPath or "data/"

	local gpOptions = CRS.gpOptions or
	{	
		grid = true, 
		logscale = true,
		mtics = true
	 }

	-- check for methods
	local CreateApproxSpace = 	CRS.CreateApproxSpace
	local CreateDomainDisc = 	CRS.CreateDomainDisc
	local CreateSolver = 		CRS.CreateSolver
	local CreateDomain = 		CRS.CreateDomain
	local MaxLevelPadding = 	CRS.MaxLevelPadding 	or util.rates.kinetic.StdMaxLevelPadding
	
	if 	CreateApproxSpace == nil or CreateDomainDisc == nil or 
		CreateSolver == nil or CreateDomain == nil then
		print("You must pass: CreateApproxSpace, CreateDomainDisc, CreateSolver, CreateDomain")
		exit()
	end
	
	local DiscTypes = CRS.DiscTypes

	local maxlevel = CRS.maxlevel; 		if maxlevel == nil then maxlevel = true end
	local prevlevel = CRS.prevlevel; 	if prevlevel == nil then prevlevel = true end
	local exact = CRS.exact; 			if exact == nil then exact = true end
	local interpol = CRS.interpol; 		if interpol == nil then interpol = true end
	local plotSol = CRS.plotSol; 		if plotSol == nil then plotSol = false end

	local ExactSol = CRS.ExactSol 
 	local ExactGrad = CRS.ExactGrad 
	local PlotCmps = CRS.PlotCmps
	
	local MeasLabel = CRS.MeasLabel or util.rates.kinetic.StdLabel.MeasLatexP
	local XLabel = 	  CRS.XLabel 	or util.rates.kinetic.StdLabel.XLatex
	local YLabel = 	  CRS.YLabel 	or util.rates.kinetic.StdLabel.YLatex
	
	local bLinear = CRS.bLinear or true
	
	--------------------------------------------------------------------
	--  Loop Discs
	--------------------------------------------------------------------

	local function ensureDir(name)
		if not(DirectoryExists(name)) then CreateDirectory(name) end
	end
	
	if plotSol then ensureDir(solPath) end
	ensureDir(dataPath)
	ensureDir(plotPath)

	-- compute element size	
	local dom = CreateDomain()
	local numRefs = dom:grid():num_levels() - 1;

	-- to store measurement
	local gpData = {};
	local errors = {};

	-- compute problem
	for _, DiscType in ipairs(DiscTypes) do

		local l2_plotdata = {}; local h1_plotdata = {}
		local l2_h_plotdata = {}; local h1_h_plotdata = {}
		local l2_dt_plotdata = {}; local h1_dt_plotdata = {}

		local disc 			= DiscType.type
		if disc == nil then print("type required."); exit(); end

		local pmin 			= DiscType.pmin or 1
		local pmax 			= DiscType.pmax or 1
		local lmin 			= DiscType.lmin or 0
		local lmax 			= DiscType.lmax or numRefs
		
		if lmin > lmax then print("lmin: "..lmin.." must be less or equal lmax: "..lmax); exit(); end
		if lmax > numRefs then print("lmax: "..lmax.." must be less or equal numRefs: "..numRefs); exit(); end
		
		local dtmin 		= DiscType.dtmin		
		local dtmax 		= DiscType.dtmax		
		local dtred 		= DiscType.dtred or 0.5
		local timeScheme	= DiscType.timeScheme
		local orderOrTheta 	= DiscType.orderOrTheta or 1
		local startTime 	= DiscType.startTime
		local endTime		= DiscType.endTime
				
		if dtmin == nil then print("dtmin required."); exit(); end
		if dtmax == nil then print("dtmax required."); exit(); end
		if dtmin > dtmax then print("dtmin: "..dtmin.." must be less or equal dtmax: "..dtmax); exit(); end
		if dtred >= 1 then print("dtred: "..dtred.." must be less than 1"); exit(); end
		if timeScheme == nil then print("timeScheme required."); exit(); end
		if startTime == nil then print("startTime required."); exit(); end
		if endTime == nil then print("endTime required."); exit(); end
		if ExactSol == nil then print("ExactSol required."); exit(); end
		
		for p = pmin, pmax do
		
			local maxLev = lmax - MaxLevelPadding(p)
			local minLev = lmin
			
			--------------------------------------------------------------------
			--  Print Setup
			--------------------------------------------------------------------
			
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
			print("    type       = " .. disc)
			print("    order      = " .. p)
			print("\n")

			--------------------------------------------------------------------
			--  Create ApproxSpace, Disc and Solver
			--------------------------------------------------------------------
				
			print(">> Create ApproximationSpace: "..disc..", "..p)
			local approxSpace = CreateApproxSpace(dom, disc, p)
			
			print(">> Create Domain Disc: "..disc..", "..p)
			local domainDisc = CreateDomainDisc(approxSpace, disc, p)
			
			print(">> Create Solver")
			local solver = CreateSolver(approxSpace)

			-- get names in approx space
			local SpaceCmp = approxSpace:names()
			
			-- per default compute each Space-cmp
			if PlotCmps == nil then
				PlotCmps = {}
				for f = 1,#SpaceCmp do
					PlotCmps[SpaceCmp[f]] = {SpaceCmp[f]}
				end
			end

			-- check functions to measure
			for _,Cmps in pairs(PlotCmps) do
				for _, cmp in pairs(Cmps) do
					if not table.contains(SpaceCmp, cmp) then
						print("Cmp: '"..cmp.."' not contained in ApproxSpace.");
						exit()
					end
				end
			end

			--------------------------------------------------------------------
			--  Create Solutions on each level
			--------------------------------------------------------------------
			
			write(">> Allocating storage for solution vectors.\n")
			local u = {}
			for lev = minLev, maxLev do
				u[lev] = GridFunction(approxSpace, lev)
				Interpolate(ExactSol, u[lev], "c", startTime)
			end
									
			--------------------------------------------------------------------------------
			--  Apply Solver
			--------------------------------------------------------------------------------
	
			local h = {}
			for lev = minLev, maxLev do
				h[lev] = MaxElementDiameter(dom, lev)
			end	
			
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
			
					if bLinear == true then
						util.SolveLinearTimeProblem(u[lev], domainDisc, solver, nil, nil,
													timeScheme, orderOrTheta, startTime, endTime, dt);
					else
						util.SolveNonlinearTimeProblem(u[lev], domainDisc, solver, nil, nil,
														timeScheme, orderOrTheta, startTime, endTime, dt);
					end
					
					if lev == maxLev and k == kmax then
						uMostAccurate = u[lev]:clone();
					end			
					
					-- compute error
					quadOrder = p+3
					l2exact[lev][k] = L2Error(ExactSol, u[lev], "c", endTime, quadOrder)
					l2diff[lev][k] = L2Error(uMostAccurate, "c", u[lev], "c", quadOrder)
					h1exact[lev][k] = H1Error(ExactSol, ExactGrad, u[lev], "c", endTime, quadOrder)
					h1diff[lev][k] = H1Error(uMostAccurate, "c", u[lev], "c", quadOrder)
					numDoFs[lev] = u[lev]:size()
					write(">> L2-Error on Level "..lev..", dt: "..dt.." is "..string.format("%.3e", l2exact[lev][k]) .."\n");
					write(">> L2-Diff  on Level "..lev..", dt: "..dt.." is "..string.format("%.3e", l2diff[lev][k]) .."\n");
					write(">> H1-Error on Level "..lev..", dt: "..dt.." is "..string.format("%.3e", h1exact[lev][k]) .."\n");
					write(">> H1-Diff  on Level "..lev..", dt: "..dt.." is "..string.format("%.3e", h1diff[lev][k]) .."\n");
					write(">> #DoF     on Level "..lev.." is "..numDoFs[lev] .."\n");
				end	
			end
	
			local function computeRate(data, min, max, base, style, levOrK)
			
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
				else
					print("Style not foung. Use 'L', 'dt'."); exit();
				end
			
				rate[min] = "--"
				for k = min+1, max do	
					fac[k] = cpyData[k-1]/cpyData[k]
					rate[k] = math.log(fac[k]) / math.log(base)
					if math.abs(rate[k]) < 1e-10 then rate[k] = "--" end
				end
						
				return cpyData, rate
			end
							 		
			-- write: for each time step size show convergence in space, keeping time step constant
			for k = 1, #dts do
				local l2ex, l2exactrate = computeRate(l2exact, minLev, maxLev, 2, "L", k)
				local h1ex, h1exactrate = computeRate(h1exact, minLev, maxLev, 2, "L", k)
				local l2di, l2diffrate =  computeRate(l2diff,  minLev, maxLev, 2, "L", k)
				local h1di, h1diffrate =  computeRate(h1diff,  minLev, maxLev, 2, "L", k)
				local level = {}; for lev = minLev, maxLev do level[lev] = lev end
				local dtout = {}; for lev = minLev, maxLev do dtout[lev] = dts[k] end
				
				-- write data to screen
				print("\nConvergence in space: fixed dt = "..dts[k])
				table.print({level, h, dtout, numDoFs, l2di, l2diffrate, l2ex, l2exactrate, h1di, h1diffrate, h1ex, h1exactrate}, 
							{heading = {"L", "h", "dt", "#DoFs", "L2 Diff", "rate", "L2 Exact", "rate", "H1 Diff", "rate", "H1 Exact", "rate"}, 
							 format = {"%d", "%.2e", "%.2e", "%d", "%.2e", "%.3f", "%.2e", "%.3f", "%.2e", "%.3f", "%.2e", "%.3f"},
							 hline = true, vline = true})
				
				-- write data to file
				gnuplot.write_data(dataPath.."error_in_h__constdt"..k.."__"..disc..p..".dat", {numDoFs, l2ex, h1ex})
				
				-- schedule l2 plot
				l2_h_plotdata = gnuplot.array_concat(l2_h_plotdata,
				{{label="L_2 error P_"..p.." dt: "..dts[k], file=dataPath.."error_in_h__constdt"..k.."__"..disc..p..".dat", style="linespoints", 1, 2 }})
				
				-- schedule h1 plot
				h1_h_plotdata = gnuplot.array_concat(h1_h_plotdata, 
				{{label="H^1 error P_"..p.." dt: "..dts[k], file=dataPath.."error_in_h__constdt"..k.."__"..disc..p..".dat", style="linespoints", 1, 3 }})	
			end
			
			-- write: for each grid size show convergence in time, keeping grid size constant
			for lev = minLev, maxLev do
				local l2ex, l2exactrate = computeRate(l2exact, 1, #dts, 1/dtred, "dt", lev)
				local h1ex, h1exactrate = computeRate(h1exact, 1, #dts, 1/dtred, "dt", lev)
				local l2di, l2diffrate =  computeRate(l2diff,  1, #dts, 1/dtred, "dt", lev)
				local h1di, h1diffrate =  computeRate(h1diff,  1, #dts, 1/dtred, "dt", lev)
				local numDoFsout = {}; for k = 1, #dts do numDoFsout[k] = numDoFs[lev] end
				local level = {}; for k = 1, #dts do level[k] = lev end
				local hout  = {}; for k = 1, #dts do hout[k] = h[lev] end
				local dtout = {}; for k = 1, #dts do dtout[k] = dts[k] end
				
				-- write data to screen
				print("\nConvergence in time: fixed h = "..string.format("%.3e", h[lev])..", #DoF = "..numDoFs[lev])
				table.print({level, hout, dtout, numDoFsout, l2di, l2diffrate, l2ex, l2exactrate, h1di, h1diffrate, h1ex, h1exactrate}, 
							{heading = {"L", "h", "dt", "#DoFs", "L2 Diff", "rate", "L2 Exact", "rate", "H1 Diff", "rate", "H1 Exact", "rate"}, 
							 format = {"%d", "%.2e", "%.2e", "%d", "%.2e", "%.3f", "%.2e", "%.3f", "%.2e", "%.3f", "%.2e", "%.3f"},
							 hline = true, vline = true})
				
				-- write data to file
				gnuplot.write_data(dataPath.."error_in_dt__constdof"..numDoFs[lev].."__"..disc..p..".dat", {dtout, l2ex, h1ex})
				
				-- schedule l2 plot
				l2_dt_plotdata = gnuplot.array_concat(l2_dt_plotdata,
				{{label="L_2 error P_"..p.." #dof: "..numDoFs[lev], file=dataPath.."error_in_dt__constdof"..numDoFs[lev].."__"..disc..p..".dat", style="linespoints", 1, 2 }})
				
				-- schedule h1 plot
				h1_dt_plotdata = gnuplot.array_concat(h1_dt_plotdata, 
				{{label="H^1 error P_"..p.." #dof: "..numDoFs[lev], file=dataPath.."error_in_dt__constdof"..numDoFs[lev].."__"..disc..p..".dat", style="linespoints", 1, 3 }})	
			end

			-- write data to file
			gnuplot.write_data(dataPath.."error__"..disc..p..".dat", {numDoFs, dts, l2exact, h1exact})
			
			-- schedule l2 plot
			l2_plotdata = gnuplot.array_concat(l2_plotdata,
			{{label="L_2 error P_"..p, file=dataPath.."error__"..disc..p..".dat", style="linespoints", 1, 2, 3}})
			
			-- schedule h1 plot
			h1_plotdata = gnuplot.array_concat(h1_plotdata, 
			{{label="H^1 error P_"..p, file=dataPath.."error__"..disc..p..".dat", style="linespoints", 1, 2, 4}})	
		end
		
		-- create plots
		gnuplot.plot(plotPath..disc.."_h_l2.pdf", l2_h_plotdata, {grid = true, logscale = true, label = {x = "#DoF"}})
		gnuplot.plot(plotPath..disc.."_h_h1.pdf", h1_h_plotdata, {grid = true, logscale = true, label = {x = "#DoF"}})

		gnuplot.plot(plotPath..disc.."_dt_l2.pdf", l2_dt_plotdata, {grid = true, logscale = true, label = {x = "dt"}})
		gnuplot.plot(plotPath..disc.."_dt_h1.pdf", h1_dt_plotdata, {grid = true, logscale = true, label = {x = "dt"}})
		
		gnuplot.plot(plotPath..disc.."_l2.pdf", l2_plotdata, {grid = true, logscale = true, label = {x = "#DoF", y = "dt"}})
		gnuplot.plot(plotPath..disc.."_h1.pdf", h1_plotdata, {grid = true, logscale = true, label = {x = "#DoF", y = "dt"}})
		
	end
end
