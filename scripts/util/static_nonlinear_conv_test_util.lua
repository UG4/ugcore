function util.computeNonlinearStaticConvRatesForSpace(dom, maxLev, minLev, discType, p,
													exactSol, exactGrad, 
													createApproxSpace, createDomainDisc, createSolver,
													startTimeLoop,
									 				plotPath, solPath, dataPath)

print("\n")
print("---------------------------")
print(" General parameters chosen:")
print("---------------------------")
print("    dim        = " .. dim)
print("    grid       = " .. gridName)
print("    maxLev     = " .. maxLev)
print("    minLev     = " .. minLev)
print("    type       = " .. discType)
print("    order      = " .. p)
print("\n")


-- create Approximation Space
print(">> Create ApproximationSpace: "..discType..", "..p)
local approxSpace = createApproxSpace(dom, discType, p)

print(">> Create Domain Disc: "..discType..", "..p)
local domainDisc = createDomainDisc(discType, p, approxSpace)

print(">> Create Solver")
local newtonSolver = createSolver(approxSpace)

--------------------------------------------------------------------------------
--  Apply Solver
--------------------------------------------------------------------------------

local l2exact = {}
local l2diff = {}
local l2levdiff = {}
local h1exact = {}
local h1diff = {}
local h1levdiff = {}
local numDoFs = {}
local u = {}

for lev = minLev, maxLev do
	u[lev] = GridFunction(approxSpace, lev)
end

-- perform a time loop to find inial guess
if startTimeLoop ~= nil then
	local dtmin 		= startTimeLoop.dtmin		
	local dtmax 		= startTimeLoop.dtmax		
	local dtred 		= startTimeLoop.dtred		
	local timeScheme	= startTimeLoop.timeScheme
	local orderOrTheta 	= startTimeLoop.orderOrTheta
	local startTime 	= startTimeLoop.startTime
	local endTime		= startTimeLoop.endTime
	local interpolStart = startTimeLoop.interpolateStartvalues

	interpolStart(u[minLev], startTime)
	
	util.SolveNonlinearTimeProblem(u[minLev], domainDisc, newtonSolver, nil, nil,
									timeScheme, orderOrTheta, startTime, endTime, 
									dtmax, dtmin, dtred);	
									
	write("\n>> --------------------------------\n")									
	write(">> Preparing inital guess finished.\n")									
	write(">> --------------------------------\n\n")									
end

	write(">> -------------------------------------------\n")									
	write(">> Starting non-linear solution on each level.\n")									
	write(">> -------------------------------------------\n\n")									
	
for lev = minLev, maxLev do

	write("\n>> Solving level "..lev..", #DoFs: "..u[lev]:size().."\n")
	
	newtonSolver:init(AssembledOperator(domainDisc, u[lev]:grid_level()))
	if newtonSolver:apply(u[lev]) == false then
		 print ("Newton solver apply failed."); exit();
	end
	write(">> Newton Solver done.\n")
	--WriteGridFunctionToVTK(u[lev], "Sol"..lev)
	
	if lev < maxLev then	
		Prolongate(u[lev+1], u[lev]);
		write(">> Solution interpolated as start value for finer level.\n")
		--WriteGridFunctionToVTK(u[lev+1], "Sol"..lev+1)
	end
end
	
for lev = maxLev, minLev, -1 do
	quadOrder = p+2
	if exactSol ~= nil then	l2exact[lev] = L2Error(exactSol, u[lev], "c", 0.0, quadOrder)
	else l2exact[lev] = 0.0; end
	l2diff[lev] = L2Error(u[maxLev], "c", u[lev], "c", quadOrder)
	if lev > minLev then l2levdiff[lev] = L2Error(u[lev], "c", u[lev-1], "c", quadOrder)
	else l2levdiff[lev] = 0.0 end
	
	if exactSol ~= nil and exactGrad ~= nil then h1exact[lev] = H1Error(exactSol, exactGrad, u[lev], "c", 0.0, quadOrder)
	else h1exact[lev] = 0.0 end
	h1diff[lev] = H1Error(u[maxLev], "c", u[lev], "c", quadOrder)
	if lev > minLev then h1levdiff[lev] = H1Error(u[lev], "c", u[lev-1], "c", quadOrder)
	else h1levdiff[lev] = 0.0 end
	
	numDoFs[lev] = u[lev]:size()
	write(">> L2 l-exact on Level "..lev.." is "..string.format("%.3e", l2exact[lev]) .."\n");
	write(">> L2 l-lmax  on Level "..lev.." is "..string.format("%.3e", l2diff[lev]) .."\n");
	write(">> L2 l-(l-1) on Level "..lev.." is "..string.format("%.3e", l2levdiff[lev]) .."\n");
	write(">> H1 l-exact on Level "..lev.." is "..string.format("%.3e", h1exact[lev]) .."\n");
	write(">> H1 l-lmax  on Level "..lev.." is "..string.format("%.3e", h1diff[lev]) .."\n");
	write(">> H1 l-(l-1) on Level "..lev.." is "..string.format("%.3e", h1levdiff[lev]) .."\n");
	write(">> #DoF     on Level "..lev.." is "..numDoFs[lev] .."\n");
	
	-- 6. write solution
	WriteGridFunctionToVTK(u[lev], solPath.."sol_"..discType..p.."_l"..lev)
	write(">> Solution written to: "..solPath.."sol_"..discType..p.."_l"..lev.."\n");	
end

return l2exact, l2diff, l2levdiff, h1exact, h1diff, h1levdiff, numDoFs
end



function util.computeNonlinearStaticConvRates(dom, createApproxSpace, createDomainDisc, createSolver, DiscTypes)
	
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
	local startTimeLoop = DiscTypes.startTimeLoop

	-- compute problem
	for type = 1,#DiscTypes do
		local l2plotdata = {}; local h1plotdata = {}
		local discType 	= DiscTypes[type].type
		local pmin 		= DiscTypes[type].pmin
		local pmax 		= DiscTypes[type].pmax
		local lmin 		= DiscTypes[type].lmin
		local lmax 		= DiscTypes[type].lmax		
		
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
		
		for p = pmin, pmax do
		
			local maxLev = lmax - math.floor(p/2)
			local minLev = lmin

			local h = {}
			for lev = minLev, maxLev do
				h[lev] = MaxElementDiameter(dom, lev)
			end
			
			local l2exact, l2diff, l2levdiff, h1exact, h1diff, h1levdiff, numDoFs
				 = util.computeNonlinearStaticConvRatesForSpace(
				 		dom, maxLev, minLev, discType, p, 
				 		exactSol, exactGrad,
				 		createApproxSpace, createDomainDisc, createSolver,
				 		startTimeLoop,
				 		plotPath, solPath, dataPath)
				 		
			local l2exactfac = {}; 
			local h1exactfac = {};
			local l2exactrate = {};
			local h1exactrate = {};
			local l2difffac = {};
			local h1difffac = {};
			local l2diffrate = {};
			local h1diffrate = {};
			local l2levdifffac = {};
			local h1levdifffac = {};
			local l2levdiffrate = {};
			local h1levdiffrate = {};
			local level = {};

			l2exactfac[minLev] = "--"
			h1exactfac[minLev] = "--"
			l2exactrate[minLev] = "--"
			h1exactrate[minLev] = "--"
			l2difffac[minLev] = "--"
			h1difffac[minLev] = "--"
			l2diffrate[minLev] = "--"
			h1diffrate[minLev] = "--"
			l2levdifffac[minLev] = "--"
			h1levdifffac[minLev] = "--"
			l2levdiffrate[minLev] = "--"
			h1levdiffrate[minLev] = "--"
			l2levdifffac[minLev+1] = "--"
			h1levdifffac[minLev+1] = "--"
			l2levdiffrate[minLev+1] = "--"
			h1levdiffrate[minLev+1] = "--"
			level[minLev] = minLev
	
			for lev = minLev+1, maxLev do
				level[lev] = lev
				l2exactfac[lev] = l2exact[lev-1]/l2exact[lev]
				h1exactfac[lev] = h1exact[lev-1]/h1exact[lev]
				l2exactrate[lev] = math.log(l2exactfac[lev]) / math.log(2)
				h1exactrate[lev] = math.log(h1exactfac[lev]) / math.log(2)
			end
			
			for lev = minLev+1 ,maxLev-1 do
				l2difffac[lev] = l2diff[lev-1]/l2diff[lev]
				h1difffac[lev] = h1diff[lev-1]/h1diff[lev]
				l2diffrate[lev] = math.log(l2difffac[lev]) / math.log(2)
				h1diffrate[lev] = math.log(h1difffac[lev]) / math.log(2)
			end
			
			for lev = minLev+2 ,maxLev do
				l2levdifffac[lev] = l2levdiff[lev-1]/l2levdiff[lev]
				h1levdifffac[lev] = h1levdiff[lev-1]/h1levdiff[lev]
				l2levdiffrate[lev] = math.log(l2levdifffac[lev]) / math.log(2)
				h1levdiffrate[lev] = math.log(h1levdifffac[lev]) / math.log(2)
			end
			
				l2levdiff[minLev] = "---"
				h1levdiff[minLev] = "---"
				l2difffac[maxLev] = "--"
				h1difffac[maxLev] = "--"
				l2diffrate[maxLev] = "--"
				h1diffrate[maxLev] = "--"
	
			-- write data to screen
			if exactGrad ~= nil and exactSol ~= nil then
			table.print({level, h, numDoFs, l2diff, l2diffrate, l2levdiff, l2levdiffrate, l2exact, l2exactrate, h1diff, h1diffrate, h1levdiff, h1levdiffrate, h1exact, h1exactrate}, 
						{title = {"L", "h", "#DoFs", "L2 l-lmax", "rate", "L2 l-(l-1)", "rate", "L2 l-exact", "rate", "H1 l-lmax", "rate", "H1 l-(l-1)", "rate","H1 l-exact", "rate"}, 
						 format = {"%d", "%.2e", "%d", "%.2e", "%.3f", "%.2e", "%.3f", "%.2e", "%.3f", "%.2e", "%.3f", "%.2e", "%.3f", "%.2e", "%.3f"},
						 hline = true, vline = true})
			else
			table.print({level, h, numDoFs, l2diff, l2diffrate, l2levdiff, l2levdiffrate, h1diff, h1diffrate, h1levdiff, h1levdiffrate}, 
						{title = {"L", "h", "#DoFs", "L2 l-lmax", "rate", "L2 l-(l-1)", "rate", "H1 l-lmax", "rate", "H1 l-(l-1)", "rate",}, 
						 format = {"%d", "%.2e", "%d", "%.2e", "%.3f", "%.2e", "%.3f", "%.2e", "%.3f", "%.2e", "%.3f"},
						 hline = true, vline = true})
			end
			
			-- write data to file
			gnuplot.write_data(dataPath.."error_"..discType..p..".dat", {numDoFs, l2diff, h1diff})
			
			-- schedule l2 plot
			l2plotdata = gnuplot.array_concat(l2plotdata,
			{{label="L_2 error P_"..p, file=dataPath.."error_"..discType..p..".dat", style="linespoints", 1, 2 }})
			
			-- schedule h1 plot
			h1plotdata = gnuplot.array_concat(h1plotdata, 
			{{label="H^1 error P_"..p, file=dataPath.."error_"..discType..p..".dat", style="linespoints", 1, 3 }})	
		end
		
		-- create plots
		gnuplot.plot(plotPath..discType.."_l2.pdf", l2plotdata, {grid = true, logscale = true})
		gnuplot.plot(plotPath..discType.."_h1.pdf", h1plotdata, {grid = true, logscale = true})
	end
end
