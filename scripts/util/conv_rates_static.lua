util = util or {}
util.rates = util.rates or {}
util.rates.static = util.rates.static or {}

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
	local plotPath = CRS.plotPath or "plots/"
	local solPath  = CRS.solPath  or "sol/"
	local dataPath = CRS.dataPath or "data/"
	os.execute("mkdir " .. dataPath)
	os.execute("mkdir " .. plotPath)
	os.execute("mkdir " .. plotPath.."single/")
	os.execute("mkdir " .. solPath)

	-- check for methods
	CRS.PrepareInitialGuess = CRS.PrepareInitialGuess or util.rates.static.StdPrepareInitialGuess
	CRS.ComputeSolution = 		  CRS.ComputeSolution or util.rates.static.StdComputeLinearSolution
	
	-- compute element size	
	local dom = CRS.CreateDomain()
	local numRefs = dom:grid():num_levels() - 1;

	gnuplotFiles = {};

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
				
				WriteGridFunctionToVTK(u[lev], solPath.."sol_"..discType..p.."_l"..lev)
				write(">> Solution written to: "..solPath.."sol_"..discType..p.."_l"..lev.."\n");	
			end

			--------------------------------------------------------------------
			--  Compute Error Norms on each level
			--------------------------------------------------------------------
			
			-- create error storage and compute elem diameters
			local err = {}	
		
			-- check for exact solution
			err.bUse = {maxlevel = true, prevlevel = true}
			if CRS.ExactSol ~= nil and CRS.ExactGrad ~= nil then err.bUse.exact = true
			else												 err.bUse.exact = false end
			
			-- get names in approx space
			local FctCmp = approxSpace:names()

			-- prepare error measurement			
			err.h = {}
			err.numDoFs = {}
			err.FctCmp = FctCmp
			err.level = {}
			
			for lev = minLev, maxLev do	
				err.h[lev] = MaxElementDiameter(dom, lev) 
				err.level[lev] = lev
			end	

			-- loop levels and compute error
			for lev = maxLev, minLev, -1 do
				write("\n>> Error Norm values on Level "..lev..".\n")
				
				quadOrder = p+3
				err.numDoFs[lev] 			= u[lev]:size()
				write(">> #DoF       on Level "..lev.." is "..err.numDoFs[lev] .."\n");
			
				-- compute for each component
				for _, f in pairs(FctCmp) do

					-- create component
					err[f] = err[f] or {}
										
					-- w.r.t exact solution		
					if err.bUse.exact then 					
					err[f]["exact"] = err[f]["exact"] or {}
					
					err[f]["exact"]["l2"] = err[f]["exact"]["l2"] or {}
					err[f]["exact"]["l2"].value = err[f]["exact"]["l2"].value or {}
					err[f]["exact"]["l2"].value[lev] = L2Error(CRS.ExactSol[f], u[lev], f, 0.0, quadOrder)
					write(">> L2 l-exact for "..f.." on Level "..lev.." is "..string.format("%.3e", err[f]["exact"]["l2"].value[lev]) .."\n");

					err[f]["exact"]["h1"] = err[f]["exact"]["h1"] or {}
					err[f]["exact"]["h1"].value = err[f]["exact"]["h1"].value or {}
					err[f]["exact"]["h1"].value[lev] = H1Error(CRS.ExactSol[f], CRS.ExactGrad[f], u[lev], f, 0.0, quadOrder)
					write(">> H1 l-exact for "..f.." on Level "..lev.." is "..string.format("%.3e", err[f]["exact"]["h1"].value[lev]) .."\n");
					end
					
					-- w.r.t max level solution
					if err.bUse.maxlevel and lev < maxLev then
					err[f]["maxlevel"] = err[f]["maxlevel"] or {}
					
					err[f]["maxlevel"]["l2"] = err[f]["maxlevel"]["l2"] or {}
					err[f]["maxlevel"]["l2"].value = err[f]["maxlevel"]["l2"].value or {}
					err[f]["maxlevel"]["l2"].value[lev]	= L2Error(u[maxLev], f, u[lev], f, quadOrder)
					write(">> L2 l-lmax  for "..f.." on Level "..lev.." is "..string.format("%.3e", err[f]["maxlevel"]["l2"].value[lev]) .."\n");

					err[f]["maxlevel"]["h1"] = err[f]["maxlevel"]["h1"] or {}
					err[f]["maxlevel"]["h1"].value = err[f]["maxlevel"]["h1"].value or {}
					err[f]["maxlevel"]["h1"].value[lev] 	= H1Error(u[maxLev], f, u[lev], f, quadOrder)
					write(">> H1 l-lmax  for "..f.." on Level "..lev.." is "..string.format("%.3e", err[f]["maxlevel"]["h1"].value[lev]) .."\n");
					end
				
					-- w.r.t previous level solution
					if err.bUse.prevlevel and lev > minLev then 
					err[f]["prevlevel"] = err[f]["prevlevel"] or {}
					
					err[f]["prevlevel"]["l2"] = err[f]["prevlevel"]["l2"] or {}
					err[f]["prevlevel"]["l2"].value = err[f]["prevlevel"]["l2"].value or {}
					err[f]["prevlevel"]["l2"].value[lev] = L2Error(u[lev], f, u[lev-1], f, quadOrder)
					write(">> L2 l-(l-1) for "..f.." on Level "..lev.." is "..string.format("%.3e", err[f]["prevlevel"]["l2"].value[lev]) .."\n");

					err[f]["prevlevel"]["h1"] = err[f]["prevlevel"]["h1"] or {}
					err[f]["prevlevel"]["h1"].value = err[f]["prevlevel"]["h1"].value or {}
					err[f]["prevlevel"]["h1"].value[lev] = H1Error(u[lev], f, u[lev-1], f, quadOrder)
					write(">> H1 l-(l-1) for "..f.." on Level "..lev.." is "..string.format("%.3e", err[f]["prevlevel"]["h1"].value[lev]) .."\n");
					end
				end
								
			end

			--------------------------------------------------------------------
			--  Compute Factors and Rates
			--------------------------------------------------------------------
	
			for _, f in ipairs(FctCmp) do
				for _, t in ipairs({"exact", "maxlevel", "prevlevel"}) do
					for _, n in ipairs({"l2", "h1"}) do
			
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
			--  Prepare labels
			--------------------------------------------------------------------

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

			--------------------------------------------------------------------
			--  Write Data to Screen
			--------------------------------------------------------------------

			for _, f in ipairs(FctCmp) do
			
				print("\n>> Statistic for type: "..discType..", order: "..p..", comp: "..f.."\n")			
				
				local values = {err.level, err.h, err.numDoFs}
				local heading = {"L", "h", "#DoFs"}
				local format = {"%d", "%.2e", "%d"}

				for _, t in ipairs({"exact", "maxlevel", "prevlevel"}) do
					for _, n in ipairs({"l2", "h1"}) do
			
						local meas = err[f][t][n]
			
						table.append(values, {meas.value, meas.rate}) 
						table.append(heading,{titles[n].." "..titles[t], "rate"})
						table.append(format, {"%.2e", "%.3f"})
					end
				end
										
				table.print(values, {heading = heading, format = format, 
									 hline = true, vline = true, forNil = "--"})
			end

			--------------------------------------------------------------------
			--  Write Data to GnuPlot
			--------------------------------------------------------------------
			
			for _, f in ipairs(FctCmp) do
				for _, t in ipairs({"exact", "prevlevel", "maxlevel"}) do

					if err[f][t] ~= nil then
						local function schedule(err, file, data, norm, x)
						
							gnuplotFiles[file] = gnuplotFiles[file] or {} 				
							table.append(gnuplotFiles[file], data)
							gnuplotFiles[file].title = gpNorm[norm]..gpTitle[t].." for Fct "..f
							gnuplotFiles[file].xlabel = gpXLabel[x]
							gnuplotFiles[file].ylabel = "|| "..f.."_L - "..f.."_{"..gpType[t].."} ||_{ "..gpNorm[norm].."}"
						end

						for _, n in ipairs({"l2", "h1"}) do
			
							local meas = err[f][t][n]
							if meas ~= nil then
																					
								-- data values
								local l2value = err[f][t]["l2"].value
								local h1value = err[f][t]["h1"].value
							
								-- write l2 and h1 to data file
								local singleFile = table.concat({"error",titles[t],n,discType,p,f},"_")
								local file = dataPath..singleFile..".dat"
								local dataCols = {err.numDoFs, err.h, err[f][t][n].value}
								gnuplot.write_data(file, dataCols)
							
								-- create plot for single run
								local options = {grid = true, logscale = true}
								local style = "linespoints"
								
								for x, xCol in pairs({DoF = 1, h = 2}) do
									local Data = {{label=discType.." P_"..p, file=file, style=style, xCol, 3}}
									
									local file = table.concat({plotPath.."single/"..singleFile,n,x,".pdf"},"_")
									gnuplot.plot(file, Data, options)			
									
									local file = table.concat({plotPath..discType,titles[t],f,n,x,".pdf"}, "_")	
									schedule(err, file, Data, n, x)
									
									local file = table.concat({plotPath.."all",titles[t],f,n,x,".pdf"}, "_")
									schedule(err, file, Data, n, x)
								end	
							end
						end
					end
				end	
			end
			
		end -- end loop over p
		
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
		
	end -- end loop over type
end
