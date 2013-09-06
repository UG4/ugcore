
function util.resetErrorRatesStorage(errRates, minLev, maxLev, FctCmp, defValue)

	if defValue == nil then defValue = "--" end
	
	errRates.minLev = minLev
	errRates.maxLev = maxLev
	errRates.numDoFs = {}
	errRates.FctCmp = FctCmp
	
	function util.createNormStorage(norm)

		function util.createNormTypeStorage(normType)

			function util.createNormTypeRatesStorage(array)
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
				util.createNormTypeRatesStorage(normType.value[f])
				normType.fac[f] = {}
				util.createNormTypeRatesStorage(normType.fac[f])
				normType.rate[f] = {}
				util.createNormTypeRatesStorage(normType.rate[f])
			end
		end
		
		-- error w.r.t to exact solution
		if errRates.bUseExact then
		norm.exact = {}
		norm.exact.title = "l-exact"
		util.createNormTypeStorage(norm.exact)
		end
		
		-- error w.r.t to most refined level  
		if errRates.bMaxLevel then
		norm.maxlevel = {}  		
		norm.maxlevel.title = "l-lmax"
		util.createNormTypeStorage(norm.maxlevel)
		end
		
		-- error w.r.t to lev-1
		if errRates.bPrevLevel then
		norm.prevlevel = {}
		norm.prevlevel.title = "l-prev"
		util.createNormTypeStorage(norm.prevlevel)
		end
	end
	
	errRates.l2 = {}
	errRates.l2.title = "L2"
	util.createNormStorage(errRates.l2)

	errRates.h1 = {}
	errRates.h1.title = "H1"
	util.createNormStorage(errRates.h1)

	errRates.level = {}
	for lev = minLev, maxLev do
		errRates.level[lev] = lev
	end
	
	return errRates
end

function util.computeErrorRates(errRates)
	
	local FctCmp = errRates.FctCmp
	function util.computeNormErrorRates(norm)
	
		function util.computeNormTypeErrorRates(normType, minLev, maxLev)
			for i = 1, #FctCmp do
				local f = FctCmp[i]
				for lev = minLev + 1, maxLev do
					normType.fac[f][lev] = normType.value[f][lev-1]/normType.value[f][lev]
					normType.rate[f][lev] = math.log(normType.fac[f][lev]) / math.log(2)
				end		
			end
		end

		local minLev = errRates.minLev
		local maxLev = errRates.maxLev
		
		if errRates.bUseExact then	util.computeNormTypeErrorRates(norm.exact,     minLev, maxLev) end
		if errRates.bMaxLevel then  util.computeNormTypeErrorRates(norm.maxlevel,  minLev, maxLev - 1) end 
		if errRates.bPrevLevel then util.computeNormTypeErrorRates(norm.prevlevel, minLev + 1, maxLev) end
	end 		
	
	util.computeNormErrorRates(errRates.l2)
	util.computeNormErrorRates(errRates.h1)
end

function util.updateErrorRatesScreenOutput(err, f)

	function util.addConvRateOutputNorm(err, norm)
		
		function util.addConvRateOutputNormType(err, norm, type)
			
			err.screen.values = gnuplot.array_concat(err.screen.values, {type.value[f], type.rate[f]}) 
			err.screen.title  = gnuplot.array_concat(err.screen.title,  {norm.title.." "..type.title, "rate"})
			err.screen.format = gnuplot.array_concat(err.screen.format, {"%.2e", "%.3f"})
		end
									
		if err.bUseExact  then util.addConvRateOutputNormType(err, norm, norm.exact) end
		if err.bMaxLevel  then util.addConvRateOutputNormType(err, norm, norm.maxlevel) end
		if err.bPrevLevel then util.addConvRateOutputNormType(err, norm, norm.prevlevel) end
	end
	
	err.screen = {}
	err.screen.values = {err.level, err.h, err.numDoFs}
	err.screen.title = {"L", "h", "#DoFs"}
	err.screen.format = {"%d", "%.2e", "%d"}
	
	util.addConvRateOutputNorm(err, err.l2)
	util.addConvRateOutputNorm(err, err.h1)
end

function util.writeAndScheduleGnuplotData(err, discType, p)

	local FctCmp = err.FctCmp

	function util.writeAndScheduleGnuplotDataType(err, discType, p, f, errSuffix, l2value, h1value,
												  typeString)
	
		-- write l2 and h1 to data file
		local singleFileName = "error_"..errSuffix.."_"..discType.."_"..p.."_"..f
		local file = err.dataPath..singleFileName..".dat"
		local dataCols = {err.numDoFs, err.h, l2value, h1value}
		gnuplot.write_data(file, dataCols)
	
		-- create plot for single run
		local options = {grid = true, logscale = true}
		local style = "linespoints"
		
		local L2_DoF_Data = {{label=discType.." P_"..p, file=file, style=style, 1, 3 }}
		local H1_DoF_Data = {{label=discType.." P_"..p, file=file, style=style, 1, 4 }}
		local L2_h_Data = {{label=discType.." P_"..p, file=file, style=style, 2, 3 }}
		local H1_h_Data = {{label=discType.." P_"..p, file=file, style=style, 2, 4 }}
		
		gnuplot.plot(err.plotPath.."single/"..singleFileName.."_l2_DoF.pdf", L2_DoF_Data, options)
		gnuplot.plot(err.plotPath.."single/"..singleFileName.."_h1_DoF.pdf", H1_DoF_Data, options)
		gnuplot.plot(err.plotPath.."single/"..singleFileName.."_l2_h.pdf", L2_h_Data, options)
		gnuplot.plot(err.plotPath.."single/"..singleFileName.."_h1_h.pdf", H1_h_Data, options)
	
		function util.scheduleGnuplotFilename(err, plotFileName, data, normString, yLabelString)
		
			if err.gnuplot[plotFileName] == nil then 
				err.gnuplot[plotFileName] = {} 				
			end
			err.gnuplot[plotFileName] = gnuplot.array_concat(err.gnuplot[plotFileName], data)

			local titelWRT = ""
			if typeString == "exact" 		then titelWRT = " w.r.t. exact Solution" end
			if typeString == "maxlevel" 	then titelWRT = " w.r.t. finest Solution" end
			if typeString == "prevlevel" 	then titelWRT = " w.r.t. previous level Solution" end

			local cmpSolSuffix = ""
			if typeString == "exact" 		then cmpSolSuffix = "exact" end
			if typeString == "maxlevel" 	then cmpSolSuffix = "L_{max}" end
			if typeString == "prevlevel" 	then cmpSolSuffix = "L-1" end
			
			err.gnuplot[plotFileName].title = normString.."-Error"..titelWRT.." for Fct "..f
			err.gnuplot[plotFileName].xlabel = yLabelString
			err.gnuplot[plotFileName].ylabel = "|| "..f.."_L - "..f.."_{"..cmpSolSuffix.."} ||_{ "..normString.."}"
		end
		
		-- schedule for plots of same disc type
		util.scheduleGnuplotFilename(err, err.plotPath..discType.."_"..errSuffix.."_"..f.."_l2_DoF.pdf", L2_DoF_Data, "L_2", "# DoFs")
		util.scheduleGnuplotFilename(err, err.plotPath..discType.."_"..errSuffix.."_"..f.."_h1_DoF.pdf", H1_DoF_Data, "H^1", "# DoFs")
		util.scheduleGnuplotFilename(err, err.plotPath..discType.."_"..errSuffix.."_"..f.."_l2_h.pdf", L2_h_Data, "L_2", "h (mesh size)")
		util.scheduleGnuplotFilename(err, err.plotPath..discType.."_"..errSuffix.."_"..f.."_h1_h.pdf", H1_h_Data, "H^1", "h (mesh size)")
	
		-- schedule for plots of all types
		util.scheduleGnuplotFilename(err, err.plotPath.."all_"..errSuffix.."_"..f.."_l2_DoF.pdf", L2_DoF_Data, "L_2", "# DoFs")
		util.scheduleGnuplotFilename(err, err.plotPath.."all_"..errSuffix.."_"..f.."_h1_DoF.pdf", H1_DoF_Data, "H^1", "# DoFs")
		util.scheduleGnuplotFilename(err, err.plotPath.."all_"..errSuffix.."_"..f.."_l2_h.pdf", L2_h_Data, "L_2", "h (mesh size)")
		util.scheduleGnuplotFilename(err, err.plotPath.."all_"..errSuffix.."_"..f.."_h1_h.pdf", H1_h_Data, "H^1", "h (mesh size)")
	end
	
	if err.bUseExact then
		for i = 1, #FctCmp do
		local f = FctCmp[i]
		util.writeAndScheduleGnuplotDataType(
			err, discType, p, f, err.l2.exact.title, err.l2.exact.value[f], err.h1.exact.value[f], "exact")
		end
	end	
	if err.bPrevLevel then
		for i = 1, #FctCmp do
		local f = FctCmp[i]
		util.writeAndScheduleGnuplotDataType(
			err, discType, p, f, err.l2.prevlevel.title, err.l2.prevlevel.value[f], err.h1.prevlevel.value[f], "prevlevel")
		end
	end	
	if err.bMaxLevel then
		for i = 1, #FctCmp do
		local f = FctCmp[i]
		-- finest level compared to finest level is not senseful --> remove it
		err.numDoFs[err.maxLev] = nil
		err.h[err.maxLev] = nil
		err.l2.maxlevel.value[f][err.maxLev] = nil
		err.h1.maxlevel.value[f][err.maxLev] = nil
				
		util.writeAndScheduleGnuplotDataType(
			err, discType, p, f, err.l2.maxlevel.title, err.l2.maxlevel.value[f], err.h1.maxlevel.value[f], "maxlevel")
		end
	end	
end

function util.computeStaticConvRates_StdPrepareInitialGuess(u, lev, minLev, maxLev,
															domainDisc, solver)

	if lev > minLev then	
		Prolongate(u[lev], u[lev-1]);
		write(">> Solution interpolated as start value from coarser level.\n")
	else
		u[lev]:set(0.0)
		write(">> Solution set to zero on coarsest level.\n")	
	end		
end


function util.computeStaticConvRates_StdComputeLinearSolution(u, lev, approxSpace, domainDisc, solver)

	-- create operator from discretization
	local A = AssembledLinearOperator(domainDisc)
	local b = GridFunction(approxSpace, lev)
	write(">> Algebra created.\n")
	
	-- 1. init operator
	domainDisc:assemble_linear(A, b)
	write(">> Matrix and Rhs assembled.\n")
	
	-- 2. set dirichlet values in start iterate
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
end

function util.computeStaticConvRates_StdComputeNonLinearSolution(u, lev, approxSpace, domainDisc, solver)

	solver:init(AssembledOperator(domainDisc, u[lev]:grid_level()))
	if solver:apply(u[lev]) == false then
		 print (">> Newton solver apply failed."); exit();
	end
	write(">> Newton Solver done.\n")
end

function util.computeConvRatesStatic(ConvRateSetup)
	
	-- create directories
	if ConvRateSetup.plotPath == nil then  ConvRateSetup.plotPath = "plots/" end
	if ConvRateSetup.solPath == nil then  ConvRateSetup.solPath = "sol/" end
	if ConvRateSetup.dataPath == nil then  ConvRateSetup.dataPath = "data/" end
	os.execute("mkdir " .. ConvRateSetup.dataPath)
	os.execute("mkdir " .. ConvRateSetup.plotPath)
	os.execute("mkdir " .. ConvRateSetup.plotPath.."single/")
	os.execute("mkdir " .. ConvRateSetup.solPath)

	-- check for methods
	if ConvRateSetup.prepareInitialGuess == nil then
		ConvRateSetup.prepareInitialGuess = 
			util.computeStaticConvRates_StdPrepareInitialGuess
	end
	if ConvRateSetup.computeSolution == nil then
		ConvRateSetup.computeSolution = 
			util.computeStaticConvRates_StdComputeLinearSolution
	end
	
	-- compute element size	
	local dom = ConvRateSetup.createDomain()
	local numRefs = dom:grid():num_levels() - 1;

	-- create error storage and compute elem diameters
	local err = {}	
	err.dataPath = ConvRateSetup.dataPath
	err.plotPath = ConvRateSetup.plotPath
	err.gnuplot = {};

	-- check for exact solution
	if ConvRateSetup.exactSol ~= nil and ConvRateSetup.exactGrad ~= nil then
		err.bUseExact = true
		err.exactSol = ConvRateSetup.exactSol
		err.exactGrad = ConvRateSetup.exactGrad
	else
		err.bUseExact = false
	end
	err.bMaxLevel = true
	err.bPrevLevel = true

	--------------------------------------------------------------------
	--  Loop Discs
	--------------------------------------------------------------------
	
	local DiscTypes = ConvRateSetup.DiscTypes
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
			local approxSpace = ConvRateSetup.createApproxSpace(dom, discType, p)
			
			print(">> Create Domain Disc: "..discType..", "..p)
			local domainDisc = ConvRateSetup.createDomainDisc(discType, p, approxSpace)
			
			print(">> Create Solver")
			local solver = ConvRateSetup.createSolver(approxSpace)
			
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
				ConvRateSetup.prepareInitialGuess(u, lev, minLev, maxLev, domainDisc, solver)
				
				write(">> Computing solution on level "..lev..".\n")
				ConvRateSetup.computeSolution(u, lev, approxSpace, domainDisc, solver)
				write(">> Solver done.\n")
				
				WriteGridFunctionToVTK(u[lev], ConvRateSetup.solPath.."sol_"..discType..p.."_l"..lev)
				write(">> Solution written to: "..ConvRateSetup.solPath.."sol_"..discType..p.."_l"..lev.."\n");	
			end

			--------------------------------------------------------------------
			--  Compute Error Norms on each level
			--------------------------------------------------------------------
			
			-- get names in approx space
			local FctCmp = {}
			for i = 0, approxSpace:num_fct()-1 do FctCmp[i+1] = approxSpace:name(i) end

			-- prepare error measurement			
			err.h = {}
			for lev = minLev, maxLev do	err.h[lev] = MaxElementDiameter(dom, lev) end	
			util.resetErrorRatesStorage(err, minLev, maxLev, FctCmp)

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
					if err.bUseExact then 
					err.l2.exact.value[f][lev] 	= L2Error(err.exactSol, u[lev], f, 0.0, quadOrder)
					err.h1.exact.value[f][lev] 	= H1Error(err.exactSol, err.exactGrad, u[lev], f, 0.0, quadOrder)
					write(">> L2 l-exact for "..f.." on Level "..lev.." is "..string.format("%.3e", err.l2.exact.value[f][lev]) .."\n");
					write(">> H1 l-exact for "..f.." on Level "..lev.." is "..string.format("%.3e", err.h1.exact.value[f][lev]) .."\n");
					end
					
					-- w.r.t max level solution
					if err.bMaxLevel and lev < maxLev then 
					err.l2.maxlevel.value[f][lev] 	= L2Error(u[maxLev], f, u[lev], f, quadOrder)
					err.h1.maxlevel.value[f][lev] 	= H1Error(u[maxLev], f, u[lev], f, quadOrder)
					write(">> L2 l-lmax  for "..f.." on Level "..lev.." is "..string.format("%.3e", err.l2.maxlevel.value[f][lev]) .."\n");
					write(">> H1 l-lmax  for "..f.." on Level "..lev.." is "..string.format("%.3e", err.h1.maxlevel.value[f][lev]) .."\n");
					end
				
					-- w.r.t previous level solution
					if err.bPrevLevel then 
					if lev > minLev then 
					err.l2.prevlevel.value[f][lev] = L2Error(u[lev], f, u[lev-1], f, quadOrder)
					err.h1.prevlevel.value[f][lev] = H1Error(u[lev], f, u[lev-1], f, quadOrder)
					write(">> L2 l-(l-1) for "..f.." on Level "..lev.." is "..string.format("%.3e", err.l2.prevlevel.value[f][lev]) .."\n");
					write(">> H1 l-(l-1) for "..f.." on Level "..lev.." is "..string.format("%.3e", err.h1.prevlevel.value[f][lev]) .."\n");
					else
					write(">> L2 l-(l-1) for "..f.." on Level "..lev.." is "..err.l2.prevlevel.value[f][lev] .."\n");	
					write(">> H1 l-(l-1) for "..f.." on Level "..lev.." is "..err.h1.prevlevel.value[f][lev].."\n");
					end
					end
				end
								
			end
	
			util.computeErrorRates(err)				 		

			--------------------------------------------------------------------
			--  Write Erorr Data
			--------------------------------------------------------------------

			-- write data to screen
			for i = 1, #FctCmp do
				print("\n>> Statistic for type: "..discType..", order: "..p..", comp: "..FctCmp[i].."\n")			
				util.updateErrorRatesScreenOutput(err, FctCmp[i])
				table.print(err.screen.values, {title = err.screen.title, format = err.screen.format, hline = true, vline = true})
			end

			-- write data to gnuplot						
			util.writeAndScheduleGnuplotData(err, discType, p)
		end
		
		-- create scheduled plots (maybe overwriting several times)
		for plotFile, data in pairs(err.gnuplot) do 
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
