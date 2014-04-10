util = util or {}
util.rates = util.rates or {}
util.rates.static = util.rates.static or {}

ug_load_script("util/persistence.lua")

--------------------------------------------------------------------------------
-- Std Functions (used as defaults)
--------------------------------------------------------------------------------

function util.rates.static.StdPrepareInitialGuess(u, lev, minLev, maxLev,
															domainDisc, solver)

	if lev > minLev then	
		Prolongate(u[lev], u[lev-1]);
		write(">> Solution interpolated as start value from coarser level.\n")
	else
		u[lev]:set(0.0)
		write(">> Solution set to zero on coarsest level.\n")	
	end		
end


function util.rates.static.StdComputeLinearSolution(u, domainDisc, solver)

	-- create operator from discretization
	local A = AssembledLinearOperator(domainDisc)
	local b = u:clone()
	write(">> Algebra created.\n")
	
	-- 1. init operator
	domainDisc:assemble_linear(A, b)
	write(">> Matrix and Rhs assembled.\n")
	
	-- 2. set dirichlet values in start iterate
	domainDisc:adjust_solution(u)
	write(">> Inital guess for solution prepared.\n")
	
	-- 3. init solver for linear Operator
	solver:init(A, u)
	write(">> Linear Solver initialized.\n")
	
	-- 4. apply solver
	if solver:apply_return_defect(u, b) ~= true then
		write(">> Linear solver failed. Aborting."); exit();
	end
end

function util.rates.static.StdComputeNonLinearSolution(u, domainDisc, solver)

	solver:init(AssembledOperator(domainDisc, u:grid_level()))
	if solver:apply(u) == false then
		 print (">> Newton solver apply failed."); exit();
	end
	write(">> Newton Solver done.\n")
end

function util.rates.static.StdMaxLevelPadding(p)
	return math.floor(p/2)
end

function util.rates.static.StdMinLevelPadding(p)
	return 0
end

--------------------------------------------------------------------------------
-- Label names
--------------------------------------------------------------------------------
	
util.rates.static.StdLabel = util.rates.static.StdLabel or {}
	
function util.rates.static.StdLabel.MeasLatexP(disc, p)
	return disc.." $\\mathbb{P}_{"..p.."}$"
end

function util.rates.static.StdLabel.MeasLatexQ(disc, p)
	return disc.." $\\mathbb{Q}_{"..p.."}$"
end

function util.rates.static.StdLabel.XLatex(x)
	local gpXLabel ={ DoFs = "Anzahl Unbekannte",	h = "h (Gitterweite)"}
	return gpXLabel[x]
end

function util.rates.static.StdLabel.YLatex(f, t, n)
	local gpType = {	["l-exact"] = 	"{}",		
						["l-lmax"] = 	"h_{\\text{min}}",
						["l-prev"] = 	"{h/2}",
					}
	local gpNorm = 	{ l2 = "L_2",	h1 = "H^1"}
	
	if t == "interpol" then
		return "$\\norm{\\mathcal{I}_h("..f..") - "..f.."}_{"..gpNorm[n].."}$"
	else
		return "$\\norm{"..f.."_h - "..f.."_{"..gpType[t].."} }_{ "..gpNorm[n].."}$"
	end
end

function util.rates.static.StdLabel.MeasPdfP(disc, p)
	return disc.." $P_"..p.."$"
end

function util.rates.static.StdLabel.MeasPdfQ(disc, p)
	return disc.." $Q_"..p.."$"
end

function util.rates.static.StdLabel.XPdf(x)
	local gpXLabel ={ DoFs = "Anzahl Unbekannte",	h = "h (Gitterweite)"}
	return gpXLabel[x]
end

function util.rates.static.StdLabel.YPdf(f, t, n)
	local gpType = {	["l-exact"] = 	"{}",		
						["l-lmax"] = 	"h_{\\text{min}}",
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
-- util.rates.static.compute (main-function)
--------------------------------------------------------------------------------

--[[!
Computes convergence rates for a static problem

In the convergence rate setup the following parameters can be passed:
- (required) CreateDomain()				
			 	function used to create Domain
- (required) CreateApproxSpace(dom, discType, p)		
			 	function used to create ApproximationSpace
- (required) CreateDomainDisc(approxSpace, discType, p)			
			 	function used to create Domain Discretization
- (required) CreateSolver(approxSpace, discType, p)				
				function used to create Solver
- (required) DiscTypes					
				Array containing types, orders and level to be looped
- (optional) ComputeSolution(u, domainDisc, solver)			
				function used to compute solution
- (optional) PrepareInitialGuess		
				function used to prepare Initial Guess
- (optional) ExactSol					
				Array containing exact solution as a function
- (optional) ExactGrad					
				Array containing exact gradients as a function
 
@param ConvRate Setup setup used 
]]--
function util.rates.static.compute(ConvRateSetup)
	
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
	local PrepareInitialGuess = CRS.PrepareInitialGuess or util.rates.static.StdPrepareInitialGuess
	local ComputeSolution = 	CRS.ComputeSolution     or util.rates.static.StdComputeLinearSolution
	local CreateApproxSpace = 	CRS.CreateApproxSpace
	local CreateDomainDisc = 	CRS.CreateDomainDisc
	local CreateSolver = 		CRS.CreateSolver
	local CreateDomain = 		CRS.CreateDomain
	local MaxLevelPadding = 	CRS.MaxLevelPadding 	or util.rates.static.StdMaxLevelPadding
	local MinLevelPadding = 	CRS.MinLevelPadding 	or util.rates.static.StdMinLevelPadding
	
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
	
	local MeasLabel = CRS.MeasLabel or util.rates.static.StdLabel.MeasLatexP
	local XLabel = 	  CRS.XLabel 	or util.rates.static.StdLabel.XLatex
	local YLabel = 	  CRS.YLabel 	or util.rates.static.StdLabel.YLatex
	
	--------------------------------------------------------------------
	--  Loop Discs
	--------------------------------------------------------------------

	local function ensureDir(name)
		if not(DirectoryExists(name)) then CreateDirectory(name) end
	end
	
	if plotSol then ensureDir(solPath) end

	-- compute element size	
	local dom = CreateDomain()
	local numRefs = dom:grid():num_levels() - 1;
	
	-- to store measurement
	local gpData = {};
	local errors = {};
	
	for _, DiscType in ipairs(DiscTypes) do
	
		local disc 			= DiscType.type
		if disc == nil then print("type required."); exit(); end

		local pmin 			= DiscType.pmin or 1
		local pmax 			= DiscType.pmax or 1
		local lmin 			= DiscType.lmin or 0
		local lmax 			= DiscType.lmax or numRefs
		
		if lmin > lmax then print("lmin: "..lmin.." must be less or equal lmax: "..lmax); exit(); end
		if lmax > numRefs then print("lmax: "..lmax.." must be less or equal numRefs: "..numRefs); exit(); end
	
		errors[disc] = errors[disc] or {}
		
		for p = pmin, pmax do
		
			errors[disc][p] = errors[disc][p] or {}
			
			local maxLev = lmax - MaxLevelPadding(p)
			local minLev = lmin + MinLevelPadding(p)
	
			--------------------------------------------------------------------
			--  Print Setup
			--------------------------------------------------------------------
	
			print("\n")
			print(">> -------------------------------------------------------------------")
			print(">>    Computing solutions and error norms of the following problem")
			print(">> -------------------------------------------------------------------")
			print(">>     dim        = " .. dim)
			print(">>     grid       = " .. gridName)
			print(">>     minLev     = " .. minLev)
			print(">>     maxLev     = " .. maxLev)
			print(">>     type       = " .. disc)
			print(">>     order      = " .. p)
			print("\n")
			
			--------------------------------------------------------------------
			--  Create ApproxSpace, Disc and Solver
			--------------------------------------------------------------------

			print(">> Create ApproximationSpace: "..disc..", "..p)
			local approxSpace = CreateApproxSpace(dom, disc, p)
			
			print(">> Create Domain Disc: "..disc..", "..p)
			local domainDisc = CreateDomainDisc(approxSpace, disc, p)
			
			print(">> Create Solver")
			local solver = CreateSolver(approxSpace, disc, p)
	
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
			end

			--------------------------------------------------------------------
			--  Compute Solution on each level
			--------------------------------------------------------------------
			--local solver = nil
			if exact or maxlevel or prevlevel then
				for lev = minLev, maxLev do
					write("\n>> Computing Level "..lev..", "..disc..", "..p..".\n")
					
					write(">> Preparing inital guess on level "..lev..".\n")
					--solver = CreateSolver(approxSpace, disc, p)
					PrepareInitialGuess(u, lev, minLev, maxLev, domainDisc, solver)
					
					write(">> Start: Computing solution on level "..lev..".\n")
					--solver = CreateSolver(approxSpace, disc, p)
					ComputeSolution(u[lev], domainDisc, solver)
					write(">> End: Solver done.\n")
					
					if plotSol then
						WriteGridFunctionToVTK(u[lev], solPath.."sol_"..disc..p.."_l"..lev)
						write(">> Solution written to: "..solPath.."sol_"..disc..p.."_l"..lev.."\n");	
					end
				end
			end
						
			approxSpace, domainDisc, solver = nil, nil, nil
			collectgarbage()

			--------------------------------------------------------------------
			--  Compute Error Norms on each level
			--------------------------------------------------------------------
			
			-- prepare error measurement			
			local err = errors[disc][p]
			err.h, err.DoFs, err.level = {}, {}, {}	
			for lev = minLev, maxLev do	
				err.h[lev] = MaxElementDiameter(dom, lev) 
				err.level[lev] = lev
				err.DoFs[lev] = u[lev]:num_dofs()
			end	

			-- loop levels and compute error
			for lev = maxLev, minLev, -1 do
				write("\n>> Error Norm values on Level "..lev..".\n")
				
				local quadOrder = p*p+3
				write(">> #DoF       on Level "..lev.." is "..err.DoFs[lev] .."\n");
			
				-- compute for each component
				for f, Cmps in pairs(PlotCmps) do

					-- create component
					err[f] = err[f] or {}
							
					-- help fct to create an measurement
					local function createMeas(f, t, n)
						err[f][t] = err[f][t] or {}
						err[f][t][n] = err[f][t][n] or {}
						err[f][t][n].value = err[f][t][n].value or {}						
						return err[f][t][n].value
					end
										
					-- check for exact solution and grad
					local solAvail, gradAvail = true, true
					for _,cmp in pairs(Cmps) do
						if ExactSol  == nil or ExactSol[cmp]  == nil then solAvail  = false end
						if ExactGrad == nil or ExactGrad[cmp] == nil then gradAvail = false end
					end
															
					-- w.r.t exact solution		
					if exact and solAvail then 					
						local value = createMeas(f, "l-exact", "l2")
						value[lev] = 0.0
						for _,cmp in pairs(Cmps) do
							value[lev] = value[lev] + math.pow(L2Error(ExactSol[cmp], u[lev], cmp, 0.0, quadOrder), 2)
						end
						value[lev] = math.sqrt(value[lev])
						write(">> L2 l-exact for "..f.." on Level "..lev.." is "..string.format("%.3e", value[lev]) .."\n");

						if gradAvail then 					
							local value = createMeas(f, "l-exact", "h1")
							value[lev] = 0.0
							for _,cmp in pairs(Cmps) do
								value[lev] = value[lev] + math.pow(H1Error(ExactSol[cmp], ExactGrad[cmp], u[lev], cmp, 0.0, quadOrder), 2)
							end
							value[lev] = math.sqrt(value[lev])
							write(">> H1 l-exact for "..f.." on Level "..lev.." is "..string.format("%.3e", value[lev]) .."\n");
						end
					end

					if interpol and solAvail then
						local uExact = u[lev]:clone()
						for _,cmp in pairs(Cmps) do
							Interpolate(ExactSol[cmp], uExact, cmp)
						end

						local value = createMeas(f, "interpol", "l2")
						value[lev] = 0.0
						for _,cmp in pairs(Cmps) do
							value[lev] = value[lev] + math.pow(L2Error(ExactSol[cmp], uExact, cmp, 0.0, quadOrder), 2)
						end
						value[lev] = math.sqrt(value[lev])
						write(">> L2 interpol for "..f.." on Level "..lev.." is "..string.format("%.3e", value[lev]) .."\n");
						
						if gradAvail then 					
							local value = createMeas(f, "interpol", "h1")
							value[lev] = 0.0
							for _,cmp in pairs(Cmps) do
								value[lev] = value[lev] + math.pow(H1Error(ExactSol[cmp], ExactGrad[cmp], uExact, cmp, 0.0, quadOrder), 2)
							end
							value[lev] = math.sqrt(value[lev])
							write(">> H1 interpol for "..f.." on Level "..lev.." is "..string.format("%.3e", value[lev]) .."\n");
						end
						
						uExact = nil						
					end
					
					-- w.r.t max level solution
					if maxlevel and lev < maxLev then
						local value = createMeas(f, "l-lmax", "l2")
						value[lev] = 0.0
						for _,cmp in pairs(Cmps) do
							value[lev] = value[lev] + math.pow(L2Error(u[maxLev], cmp, u[lev], cmp, quadOrder, "Inner, Outer"), 2)
						end
						value[lev] = math.sqrt(value[lev])
						write(">> L2 l-lmax  for "..f.." on Level "..lev.." is "..string.format("%.3e", value[lev]) .."\n");

						local value = createMeas(f, "l-lmax", "h1")
						value[lev] = 0.0
						for _,cmp in pairs(Cmps) do
							value[lev] = value[lev] + math.pow(H1Error(u[maxLev], cmp, u[lev], cmp, quadOrder, "Inner, Outer"), 2)
						end
						value[lev] = math.sqrt(value[lev])
						write(">> H1 l-lmax  for "..f.." on Level "..lev.." is "..string.format("%.3e", value[lev]) .."\n");
					end
				
					-- w.r.t previous level solution
					if prevlevel and lev < maxLev then 
						local value = createMeas(f, "l-prev", "l2")
						value[lev] = 0.0
						for _,cmp in pairs(Cmps) do
							value[lev] = value[lev] + math.pow(L2Error(u[lev+1], cmp, u[lev], cmp, quadOrder, "Inner, Outer"), 2)
						end
						value[lev] = math.sqrt(value[lev])
						write(">> L2 l-(l-1) for "..f.." on Level "..lev.." is "..string.format("%.3e", value[lev]) .."\n");
	
						local value = createMeas(f, "l-prev", "h1")
						value[lev] = 0.0
						for _,cmp in pairs(Cmps) do
							value[lev] = value[lev] + math.pow(H1Error(u[lev+1], cmp, u[lev], cmp, quadOrder, "Inner, Outer"), 2)
						end
						value[lev] = math.sqrt(value[lev])
						write(">> H1 l-(l-1) for "..f.." on Level "..lev.." is "..string.format("%.3e", value[lev]) .."\n");
					end
				end -- end fct
								
			end -- end level

			for lev = minLev, maxLev do u[lev] = nil end
			u = nil
			collectgarbage()
				
			--------------------------------------------------------------------
			--  Compute Factors and Rates
			--------------------------------------------------------------------

			for f, _ in pairs(PlotCmps) do
				for t, _ in pairs(err[f]) do
					for n, _ in pairs(err[f][t]) do

				local meas = err[f][t][n]
	
				meas.fac = meas.fac or {}
				meas.rate = meas.rate or {}
				
				local value = meas.value
				local fac = meas.fac
				local rate = meas.rate
				
				for lev, _ in pairs(value) do
					if value[lev] ~= nil and value[lev-1] ~= nil then
						fac[lev] = value[lev-1]/value[lev]
						rate[lev] = math.log(fac[lev]) / math.log(2)
					end
				end
				
					end
				end
			end	

			--------------------------------------------------------------------
			--  Write Data to Screen
			--------------------------------------------------------------------

			for f, Cmps in pairs(PlotCmps) do
			
				write("\n>> Statistic for type: "..disc..", order: "..p..", comp: "..f.." [ ")			
				for _, cmp in pairs(Cmps) do write(cmp.." ") end
				print("]\n")
				
				local values = {err.level, err.h, err.DoFs}
				local heading = {"L", "h", "#DoFs"}
				local format = {"%d", "%.2e", "%d"}

				for t, _ in pairs(err[f]) do
					for n, _ in pairs(err[f][t]) do
						local meas = err[f][t][n]
						table.append(values, {meas.value, meas.rate}) 
						table.append(heading,{n.." "..t, "rate"})
						table.append(format, {"%.2e", "%.3f"})
					end
				end
										
				table.print(values, {heading = heading, format = format, 
									 hline = true, vline = true, forNil = "--"})
			end
			
		end -- end loop over p		
	end -- end loop over type
	
	
	--------------------------------------------------------------------
	--  Write gnuplot data
	--------------------------------------------------------------------

	-- the following is serial and one proc doing it is sufficient
	if ProcRank() ~= 0 then return end
	ensureDir(dataPath)
	ensureDir(plotPath)

	-- help fct to create a plot
	local plots = {}
	local function accessPlot(...)
		local keys = arg
		local plot = plots
		for _, key in ipairs(keys) do
			plot[key] = plot[key] or {}
			plot = plot[key]
		end
		return plot
	end
	
	local function getPlot(...)
		local plot = accessPlot(...)
		local pl = table.ideepcopy( plot )
		pl.label = plot.label
		return pl
	end
	
	for disc, _ in pairs(errors) do
		for p, _ in pairs(errors[disc]) do
			for f, _ in pairs(PlotCmps) do
				for t, _ in pairs(errors[disc][p][f]) do
					for n, _ in pairs(errors[disc][p][f][t]) do
					
		-- write l2 and h1 to data file
		local file = dataPath..table.concat({"error",disc,p,f,t,n},"_")..".dat"
		local cols = {errors[disc][p].DoFs, errors[disc][p].h, errors[disc][p][f][t][n].value}
		gnuplot.write_data(file, cols)
		
		-- store dataset	
		for xCol, x in ipairs({"DoFs", "h"}) do
			local dataset = {label=MeasLabel(disc, p), file=file, style="linespoints", xCol, 3}
			local label = { x = XLabel(x), y = YLabel(f,t,n)}
							
			local function addSet(plot, dataset, label)
				table.insert( plot, dataset)			
				plot.label = label
			end
							
			addSet( accessPlot(disc, p,  f, t, n, x), dataset, label)
			addSet( accessPlot(disc,     f, t, n, x), dataset, label)
			addSet( accessPlot("all",    f, t, n, x), dataset, label)
		end
			
					end
				end
			end	
		end
	end

	--------------------------------------------------------------------
	--  Execute Plot of gnuplot
	--------------------------------------------------------------------
	
	-- one plot
	local validP
	ensureDir(plotPath.."dataset/")
	ensureDir(plotPath.."disc/")
	ensureDir(plotPath.."multi/")
	for disc, _ in pairs(errors) do
		for p, _ in pairs(errors[disc]) do
			for f, _ in pairs(PlotCmps) do
				for t, _ in pairs(errors[disc][p][f]) do
					for n, _ in pairs(errors[disc][p][f][t]) do
						for _, x in ipairs({"DoFs", "h"}) do
		
		-- single dataset			
		local file = plotPath.."dataset/"..table.concat({disc,p,f,t,n,x},"_")
		gpData[file] = getPlot(disc, p, f, t, n, x)

		-- grouping by (disc+p)								
		local file = plotPath.."disc/"..table.concat({f,disc,t,n,x}, "_")	
		gpData[file] = getPlot(disc, f, t, n, x)		

		-- grouping (all discs+p)
		local file = plotPath.."disc/"..table.concat({f,"all",t,n,x}, "_")	
		gpData[file] = getPlot("all", f, t, n, x)	


		validP = p
						end
					end
				end	
			end
		end
	end

	for disc, _ in pairs(errors) do
		local p = validP
			for f, _ in pairs(PlotCmps) do
				for t, _ in pairs(errors[disc][p][f]) do
					for n, _ in pairs(errors[disc][p][f][t]) do
						for _, x in ipairs({"DoFs", "h"}) do

		-- multi-plot: all discs for one norm
		local file = plotPath.."multi/"..table.concat({f,t,n,x}, "_")	
		gpData[file] = gpData[file] or {}
		gpData[file].multiplot = {rows = 1}
		table.insert( gpData[file], getPlot(disc, f, t, n, x) )			

		-- multi-plot: all norms for one disc
		local file = plotPath.."multi/"..table.concat({f,disc,t,x}, "_")	
		gpData[file] = gpData[file] or {}
		gpData[file].multiplot = {cols = 1}
		table.insert( gpData[file], getPlot(disc, f, t, n, x) )			

		-- multi-plot: all types for one disc and one norm
		local file = plotPath.."multi/"..table.concat({f,disc,n,x}, "_")	
		gpData[file] = gpData[file] or {}
		gpData[file].multiplot = {cols = 1}
		table.insert( gpData[file], getPlot(disc, f, t, n, x) )			

		-- multi-plot: all comps for one disc and one norm
		local file = plotPath.."multi/"..table.concat({"all",disc,t,n,x}, "_")	
		gpData[file] = gpData[file] or {}
		gpData[file].multiplot = {cols = 1}
		table.insert( gpData[file], getPlot(disc, f, t, n, x) )			

						end
					end
				end	
			end
	end

	-- multi-plot: all discs for all norm
	for _, n in pairs({"h1", "l2"}) do
	for disc, _ in pairs(errors) do
		local p = validP
			for f, _ in pairs(PlotCmps) do
				for t, _ in pairs(errors[disc][p][f]) do
					if errors[disc][p][f][t][n] then
						for _, x in ipairs({"DoFs", "h"}) do

		local file = plotPath.."multi/"..table.concat({f,t,x}, "_")	
		gpData[file] = gpData[file] or {}
		gpData[file].multiplot = {cols = 2}
		table.insert( gpData[file], getPlot(disc, f, t, n, x) )			

						end
					end
				end	
			end
	end
	end
	
	-- save for reuse
	persistence.store(dataPath.."gp-data-files.lua", gpData);	
	
	-- create scheduled plots
	if CRS.noplot == nil or CRS.noplot == false then
		for plotFile, data in pairs(gpData) do 
			local opt = table.deepcopy(gpOptions)
			if data.multiplot then opt.multiplot = data.multiplot end
			gnuplot.plot(plotFile..".tex", data, opt)
			gnuplot.plot(plotFile..".pdf", data, opt)
		end	
	end
end


function util.rates.static.replot(gpOptions, file)
	if ProcRank() ~= 0 then return end
	
	local dataPath = "data/"
	local plotPath = "plots/"
	
	local function ensureDir(name)
		if not(DirectoryExists(name)) then CreateDirectory(name) end
	end
	
	ensureDir(dataPath)
	ensureDir(plotPath)
	ensureDir(plotPath.."dataset/")
	ensureDir(plotPath.."disc/")
	ensureDir(plotPath.."multi/")
	
	local file = file or dataPath.."gp-data-files.lua"
	local gpData = persistence.load(file);
	
	-- create scheduled plots
	for plotFile, data in pairs(gpData) do 
		local opt = table.deepcopy(gpOptions)
		if data.multiplot then opt.multiplot = data.multiplot end
		gnuplot.plot(plotFile..".tex", data, opt)
		gnuplot.plot(plotFile..".pdf", data, opt)
	end
	
end