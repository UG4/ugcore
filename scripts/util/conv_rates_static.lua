util = util or {}
util.rates = util.rates or {}
util.rates.static = util.rates.static or {}

--------------------------------------------------------------------------------
-- help functions
--------------------------------------------------------------------------------

function util.rates.static.resetStorage(err, minLev, maxLev, FctCmp, defValue)

	if defValue == nil then defValue = "--" end
	
	err.minLev = minLev
	err.maxLev = maxLev
	err.numDoFs = {}
	err.FctCmp = FctCmp
	
	local function createNormStorage(norm)

		local function createNormTypeStorage(normType)

			local function createNormTypeRatesStorage(array)
				for lev = minLev, maxLev do
					array[lev] = defValue
				end
			end

			normType.value = {}
			normType.fac = {}
			normType.rate = {}
			
			for i = 1, #FctCmp do
				local f = FctCmp[i]
				normType.value[f] = {}
				createNormTypeRatesStorage(normType.value[f])
				normType.fac[f] = {}
				createNormTypeRatesStorage(normType.fac[f])
				normType.rate[f] = {}
				createNormTypeRatesStorage(normType.rate[f])
			end
		end
		
		-- error w.r.t to exact solution
		if err.bUse.exact then
		norm.exact = {}
		norm.exact.title = "l-exact"
		createNormTypeStorage(norm.exact)
		end
		
		-- error w.r.t to most refined level  
		if err.bUse.maxlevel then
		norm.maxlevel = {}  		
		norm.maxlevel.title = "l-lmax"
		createNormTypeStorage(norm.maxlevel)
		end
		
		-- error w.r.t to lev-1
		if err.bUse.prevlevel then
		norm.prevlevel = {}
		norm.prevlevel.title = "l-prev"
		createNormTypeStorage(norm.prevlevel)
		end
	end
	
	err.l2 = {}
	err.l2.title = "L2"
	createNormStorage(err.l2)

	err.h1 = {}
	err.h1.title = "H1"
	createNormStorage(err.h1)

	err.level = {}
	for lev = minLev, maxLev do
		err.level[lev] = lev
	end
	
	return err
end

function util.writeAndScheduleGnuplotData(err, gnuplotFiles, discType, p)

	local titles = {l2 = "L2", h1 = "H1",
					exact = "l-exact", maxlevel = "l-lmax", prevlevel = "l-prev"}

	local gpTitle = {	exact = 	"-Error w.r.t. exact Solution",
						maxlevel = 	"-Error w.r.t. finest Solution",
						prevlevel = "-Error w.r.t. previous level Solution"
					}
				
	local gpType = {	exact = 	"exact",		
						maxlevel = 	"L_{max}",
						prevlevel = "L-1"
					}

	local gpNorm = 	{ l2 = "L_2",	h1 = "H^1"}
					
	local gpXLabel ={ DoF = "# DoFs",	h = "h (mesh size)"}

	local function addGnuplotDataType(err, discType, p, f, t)
	
		local function schedule(err, file, data, norm, x)
		
			gnuplotFiles[file] = gnuplotFiles[file] or {} 				
			table.append(gnuplotFiles[file], data)
			gnuplotFiles[file].title = gpNorm[norm]..gpTitle[t].." for Fct "..f
			gnuplotFiles[file].xlabel = gpXLabel[x]
			gnuplotFiles[file].ylabel = "|| "..f.."_L - "..f.."_{"..gpType[t].."} ||_{ "..gpNorm[norm].."}"
		end
			
		-- data values
		local l2value = err["l2"][t].value[f]
		local h1value = err["h1"][t].value[f]
	
		-- write l2 and h1 to data file
		local singleFileName = "error_"..titles[t].."_"..discType.."_"..p.."_"..f
		local file = err.dataPath..singleFileName..".dat"
		local dataCols = {err.numDoFs, err.h, l2value, h1value}
		gnuplot.write_data(file, dataCols)
	
		-- create plot for single run
		local options = {grid = true, logscale = true}
		local style = "linespoints"
		
		for y, yCol in pairs({l2 = 3, h1 = 4}) do
			for x, xCol in pairs({DoF = 1, h = 2}) do
				local Data = {{label=discType.." P_"..p, file=file, style=style, xCol, yCol}}
				gnuplot.plot(err.plotPath.."single/"..singleFileName.."_"..y.."_"..x..".pdf", Data, options)				
				schedule(err, err.plotPath..discType.."_"..titles[t].."_"..f.."_"..y.."_"..x..".pdf", Data, y, x)
				schedule(err, err.plotPath.."all_"..titles[t].."_"..f.."_"..y.."_"..x..".pdf", Data, y, x)
			end	
		end
	end
	
	local FctCmp = err.FctCmp
	
	for _, t in ipairs({"exact", "prevlevel", "maxlevel"}) do

		for i = 1, #FctCmp do
			local f = FctCmp[i]

			-- finest level compared to finest level is not senseful --> remove it
			if t == "maxlevel" then
				err.numDoFs[err.maxLev] = nil
				err.h[err.maxLev] = nil
				err.l2.maxlevel.value[f][err.maxLev] = nil
				err.h1.maxlevel.value[f][err.maxLev] = nil
			
			end
			
			if err.bUse[t] then
				addGnuplotDataType(err, discType, p, f, t)
			end
		end	
	end
end

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
	CRS.plotPath = CRS.plotPath or "plots/"
	CRS.solPath  = CRS.solPath  or "sol/"
	CRS.dataPath = CRS.dataPath or "data/"
	os.execute("mkdir " .. CRS.dataPath)
	os.execute("mkdir " .. CRS.plotPath)
	os.execute("mkdir " .. CRS.plotPath.."single/")
	os.execute("mkdir " .. CRS.solPath)

	-- check for methods
	CRS.PrepareInitialGuess = CRS.PrepareInitialGuess or util.rates.static.StdPrepareInitialGuess
	CRS.ComputeSolution = 		  CRS.ComputeSolution or util.rates.static.StdComputeLinearSolution
	
	-- compute element size	
	local dom = CRS.CreateDomain()
	local numRefs = dom:grid():num_levels() - 1;

	-- create error storage and compute elem diameters
	local err = {}	
	err.dataPath = CRS.dataPath
	err.plotPath = CRS.plotPath
	gnuplotFiles = {};

	-- check for exact solution
	err.bUse = {maxlevel = true, prevlevel = true}
	if CRS.ExactSol ~= nil and CRS.ExactGrad ~= nil then
		err.bUse.exact = true
	else
		err.bUse.exact = false
	end

	--------------------------------------------------------------------
	--  Loop Discs
	--------------------------------------------------------------------
	
	local DiscTypes = CRS.DiscTypes
	for type = 1,#DiscTypes do
	
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
			print(">>     type       = " .. discType)
			print(">>     order      = " .. p)
			print("\n")
			
			--------------------------------------------------------------------
			--  Create ApproxSpace, Disc and Solver
			--------------------------------------------------------------------

			print(">> Create ApproximationSpace: "..discType..", "..p)
			local approxSpace = CRS.CreateApproxSpace(dom, discType, p)
			
			print(">> Create Domain Disc: "..discType..", "..p)
			local domainDisc = CRS.CreateDomainDisc(approxSpace, discType, p)
			
			print(">> Create Solver")
			local solver = CRS.CreateSolver(approxSpace, discType, p)
			
			--------------------------------------------------------------------
			--  Create Solutions on each level
			--------------------------------------------------------------------
			
			local u = {}
			
			write("\n>> Allocating storage for solution vectors.\n")
			for lev = minLev, maxLev do
				u[lev] = GridFunction(approxSpace, lev)
			end

			--------------------------------------------------------------------
			--  Compute Solution on each level
			--------------------------------------------------------------------
			
			for lev = minLev, maxLev do
				write("\n>> Computing Level "..lev..".\n")
			
				write(">> Preparing inital guess on level "..lev..".\n")
				CRS.PrepareInitialGuess(u, lev, minLev, maxLev, domainDisc, solver)
				
				write(">> Computing solution on level "..lev..".\n")
				CRS.ComputeSolution(u[lev], domainDisc, solver)
				write(">> Solver done.\n")
				
				WriteGridFunctionToVTK(u[lev], CRS.solPath.."sol_"..discType..p.."_l"..lev)
				write(">> Solution written to: "..CRS.solPath.."sol_"..discType..p.."_l"..lev.."\n");	
			end

			--------------------------------------------------------------------
			--  Compute Error Norms on each level
			--------------------------------------------------------------------
			
			-- get names in approx space
			local FctCmp = approxSpace:names()

			-- prepare error measurement			
			err.h = {}
			for lev = minLev, maxLev do	err.h[lev] = MaxElementDiameter(dom, lev) end	
			util.rates.static.resetStorage(err, minLev, maxLev, FctCmp)

			-- loop levels and compute error
			for lev = maxLev, minLev, -1 do
				write("\n>> Error Norm values on Level "..lev..".\n")
				
				quadOrder = p+3
				err.numDoFs[lev] 			= u[lev]:size()
				write(">> #DoF       on Level "..lev.." is "..err.numDoFs[lev] .."\n");
			
				-- compute for each component
				for i = 1, #FctCmp do
					local f = FctCmp[i] 
					
					-- w.r.t exact solution		
					if err.bUse.exact then 
					err.l2.exact.value[f][lev] 	= L2Error(CRS.ExactSol[f], u[lev], f, 0.0, quadOrder)
					err.h1.exact.value[f][lev] 	= H1Error(CRS.ExactSol[f], CRS.ExactGrad[f], u[lev], f, 0.0, quadOrder)
					write(">> L2 l-exact for "..f.." on Level "..lev.." is "..string.format("%.3e", err.l2.exact.value[f][lev]) .."\n");
					write(">> H1 l-exact for "..f.." on Level "..lev.." is "..string.format("%.3e", err.h1.exact.value[f][lev]) .."\n");
					end
					
					-- w.r.t max level solution
					if err.bUse.maxlevel and lev < maxLev then 
					err.l2.maxlevel.value[f][lev] 	= L2Error(u[maxLev], f, u[lev], f, quadOrder)
					err.h1.maxlevel.value[f][lev] 	= H1Error(u[maxLev], f, u[lev], f, quadOrder)
					write(">> L2 l-lmax  for "..f.." on Level "..lev.." is "..string.format("%.3e", err.l2.maxlevel.value[f][lev]) .."\n");
					write(">> H1 l-lmax  for "..f.." on Level "..lev.." is "..string.format("%.3e", err.h1.maxlevel.value[f][lev]) .."\n");
					end
				
					-- w.r.t previous level solution
					if err.bUse.prevlevel and lev > minLev then 
					err.l2.prevlevel.value[f][lev] = L2Error(u[lev], f, u[lev-1], f, quadOrder)
					err.h1.prevlevel.value[f][lev] = H1Error(u[lev], f, u[lev-1], f, quadOrder)
					write(">> L2 l-(l-1) for "..f.." on Level "..lev.." is "..string.format("%.3e", err.l2.prevlevel.value[f][lev]) .."\n");
					write(">> H1 l-(l-1) for "..f.." on Level "..lev.." is "..string.format("%.3e", err.h1.prevlevel.value[f][lev]) .."\n");
					end
				end
								
			end
	
			for _, n in ipairs({"l2", "h1"}) do
				for _, t in ipairs({"exact", "maxlevel", "prevlevel"}) do
					for _, f in ipairs(FctCmp) do
			
						local normType = err[n][t]
						local value = normType.value[f]
						
						for lev, _ in ipairs(value) do
							if value[lev] ~= "--" and value[lev-1] ~= "--" then
								normType.fac[f][lev] = value[lev-1]/value[lev]
								normType.rate[f][lev] = math.log(normType.fac[f][lev]) / math.log(2)
							end
						end
					end
				end
			end	

			--------------------------------------------------------------------
			--  Write Erorr Data
			--------------------------------------------------------------------

			local titles = {l2 = "L2", h1 = "H1",
							exact = "l-exact", maxlevel = "l-lmax", prevlevel = "l-prev"}

			-- write data to screen
			for i = 1, #FctCmp do
				local f = FctCmp[i]
			
				print("\n>> Statistic for type: "..discType..", order: "..p..", comp: "..f.."\n")			
				
				local values = {err.level, err.h, err.numDoFs}
				local title = {"L", "h", "#DoFs"}
				local format = {"%d", "%.2e", "%d"}
			
				for _, n in ipairs({"l2", "h1"}) do
					local norm = err[n]
					for _, t in ipairs({"exact", "maxlevel", "prevlevel"}) do
						local type = norm[t]
						if type ~= nil then
							table.append(values, {type.value[f], type.rate[f]}) 
							table.append(title,  {titles[n].." "..titles[t], "rate"})
							table.append(format, {"%.2e", "%.3f"})
						end
					end
				end
										
				table.print(values, {title = title, format = format, 
									 hline = true, vline = true, forNil = "--"})
			end

			-- write data to gnuplot						
			util.writeAndScheduleGnuplotData(err, gnuplotFiles, discType, p)
		end
		
		-- create scheduled plots (maybe overwriting several times)
		for plotFile, data in pairs(gnuplotFiles) do 
			local options = {	grid = true, 
								logscale = true,
								title = data.title,
								xlabel = data.xlabel,
								ylabel = data.ylabel,
								"set key left bottom",
							 }
			gnuplot.plot(plotFile, data, options)
		end
	end
end
