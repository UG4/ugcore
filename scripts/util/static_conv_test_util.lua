function util.computeLinearStaticConvRatesForSpace(dom, maxLev, minLev, discType, p,
													exactSol, exactGrad, 
													createApproxSpace, createDomainDisc, createSolver,
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
	local solver = createSolver(approxSpace)
	
	--------------------------------------------------------------------------------
	--  Apply Solver
	--------------------------------------------------------------------------------
	
	local l2exact = {}
	local l2diff = {}
	local h1exact = {}
	local h1diff = {}
	local numDoFs = {}
	local u = {}
	
	for lev = maxLev, minLev, -1 do
		
		-- create operator from discretization
		local A = AssembledLinearOperator(domainDisc)
		u[lev] = GridFunction(approxSpace, lev)
		local b = GridFunction(approxSpace, lev)
		write("\n>> Algebra created.\n")
		
		-- 1. init operator
		domainDisc:assemble_linear(A, b)
		write(">> Matrix and Rhs assembled.\n")
		
		-- 2. set dirichlet values in start iterate
		u[lev]:set(0.0)
		domainDisc:adjust_solution(u[lev])
		write(">> Inital guess for solution prepared.\n")
		
		-- print matrix for test purpose
		if debug == true then
			SaveMatrixForConnectionViewer(u, A, "Stiffness_"..discType..p.."_l"..lev..".mat")
			SaveVectorForConnectionViewer(b, "Rhs_"..discType..p.."_l"..lev..".vec")
			write(">> Linear Equation system written for debug purpose.\n")
		end
		
		-- 3. init solver for linear Operator
		solver:init(A, u[lev])
		write(">> Linear Solver initialized.\n")
		
		-- 4. apply solver
		if solver:apply_return_defect(u[lev],b) ~= true then
			write(">> Linear solver failed. Aborting."); exit();
		end
		write(">> Linear Solver done.\n")
		
		-- 5. compute error
		quadOrder = p+3
		l2exact[lev] = L2Error(exactSol, u[lev], "c", 0.0, quadOrder)
		l2diff[lev] = L2Error(u[maxLev], "c", u[lev], "c", quadOrder)
		h1exact[lev] = H1Error(exactSol, exactGrad, u[lev], "c", 0.0, quadOrder)
		h1diff[lev] = H1Error(u[maxLev], "c", u[lev], "c", quadOrder)
		numDoFs[lev] = u[lev]:size()
		write(">> L2-Error on Level "..lev.." is "..string.format("%.3e", l2exact[lev]) .."\n");
		write(">> L2-Diff  on Level "..lev.." is "..string.format("%.3e", l2diff[lev]) .."\n");
		write(">> H1-Error on Level "..lev.." is "..string.format("%.3e", h1exact[lev]) .."\n");
		write(">> H1-Diff  on Level "..lev.." is "..string.format("%.3e", h1diff[lev]) .."\n");
		write(">> #DoF     on Level "..lev.." is "..numDoFs[lev] .."\n");
		
		-- 6. write solution
		WriteGridFunctionToVTK(u[lev], solPath.."sol_"..discType..p.."_l"..lev)
		write(">> Solution written to: "..solPath.."sol_"..discType..p.."_l"..lev.."\n");	
	end
	
	return l2exact, l2diff, h1exact, h1diff, numDoFs
end



function util.computeLinearStaticConvRates(dom, createApproxSpace, createDomainDisc, createSolver, DiscTypes)
	
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
			
			local l2exact, l2diff, h1exact, h1diff, numDoFs
				 = util.computeLinearStaticConvRatesForSpace(
				 		dom, maxLev, minLev, discType, p, 
				 		exactSol, exactGrad,
				 		createApproxSpace, createDomainDisc, createSolver,
				 		plotPath, solPath, dataPath)
				 		
			local l2exactfac = {}; 
			local h1exactfac = {};
			local l2exactrate = {};
			local h1exactrate = {};
			local l2difffac = {};
			local h1difffac = {};
			local l2diffrate = {};
			local h1diffrate = {};
			local level = {};

			l2exactfac[minLev] = "--"
			h1exactfac[minLev] = "--"
			l2exactrate[minLev] = "--"
			h1exactrate[minLev] = "--"
			l2difffac[minLev] = "--"
			h1difffac[minLev] = "--"
			l2diffrate[minLev] = "--"
			h1diffrate[minLev] = "--"
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
				l2difffac[maxLev] = "--"
				h1difffac[maxLev] = "--"
				l2diffrate[maxLev] = "--"
				h1diffrate[maxLev] = "--"
	
			-- write data to screen
			table.print({level, h, numDoFs, l2diff, l2diffrate, l2exact, l2exactrate, h1diff, h1diffrate, h1exact, h1exactrate}, 
						{title = {"L", "h", "#DoFs", "L2 Diff", "rate", "L2 Exact", "rate", "H1 Diff", "rate", "H1 Exact", "rate"}, 
						 format = {"%d", "%.2e", "%d", "%.2e", "%.3f", "%.2e", "%.3f", "%.2e", "%.3f", "%.2e", "%.3f"},
						 hline = true, vline = true})
			
			-- write data to file
			gnuplot.write_data(dataPath.."error_"..discType..p..".dat", {numDoFs, l2exact, h1exact})
			
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
