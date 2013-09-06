
function util.resetErrorRatesStorage(errRates, minLev, maxLev, defValue)

	if defValue == nil then defValue = "--" end
	
	errRates.minLev = minLev
	errRates.maxLev = maxLev
	errRates.numDoFs = {}
	
	function util.createNormStorage(norm)

		function util.createNormTypeStorage(normType)

			function util.createNormTypeRatesStorage(array)
				for lev = minLev, maxLev do
					array[lev] = defValue
				end
			end

			normType.value = {}
			util.createNormTypeRatesStorage(normType.value)
			normType.fac = {}
			util.createNormTypeRatesStorage(normType.fac)
			normType.rate = {}
			util.createNormTypeRatesStorage(normType.rate)
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
	
	function util.computeNormErrorRates(norm)
	
		function util.computeNormTypeErrorRates(normType, minLev, maxLev)
			for lev = minLev + 1, maxLev do
				normType.fac[lev] = normType.value[lev-1]/normType.value[lev]
				normType.rate[lev] = math.log(normType.fac[lev]) / math.log(2)
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

function util.updateErrorRatesScreenOutput(err)

	function util.addConvRateOutputNorm(err, norm)
		
		function util.addConvRateOutputNormType(err, norm, type)
			
			err.screen.values = gnuplot.array_concat(err.screen.values, {type.value, type.rate}) 
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

	function util.writeAndScheduleGnuplotDataType(err, discType, p, errSuffix, l2value, h1value,
												  typeString)
	
		-- write l2 and h1 to data file
		local singleFileName = "error_"..errSuffix.."_"..discType.."_"..p
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
			
			err.gnuplot[plotFileName].title = normString.."-Error"..titelWRT
			err.gnuplot[plotFileName].xlabel = yLabelString
			err.gnuplot[plotFileName].ylabel = "|| u_L - u_{"..cmpSolSuffix.."} ||_{ "..normString.."}"
		end
		
		-- schedule for plots of same disc type
		util.scheduleGnuplotFilename(err, err.plotPath..discType.."_"..errSuffix.."_l2_DoF.pdf", L2_DoF_Data, "L_2", "# DoFs")
		util.scheduleGnuplotFilename(err, err.plotPath..discType.."_"..errSuffix.."_h1_DoF.pdf", H1_DoF_Data, "H^1", "# DoFs")
		util.scheduleGnuplotFilename(err, err.plotPath..discType.."_"..errSuffix.."_l2_h.pdf", L2_h_Data, "L_2", "h (mesh size)")
		util.scheduleGnuplotFilename(err, err.plotPath..discType.."_"..errSuffix.."_h1_h.pdf", H1_h_Data, "H^1", "h (mesh size)")
	
		-- schedule for plots of all types
		util.scheduleGnuplotFilename(err, err.plotPath.."all_"..errSuffix.."_l2_DoF.pdf", L2_DoF_Data, "L_2", "# DoFs")
		util.scheduleGnuplotFilename(err, err.plotPath.."all_"..errSuffix.."_h1_DoF.pdf", H1_DoF_Data, "H^1", "# DoFs")
		util.scheduleGnuplotFilename(err, err.plotPath.."all_"..errSuffix.."_l2_h.pdf", L2_h_Data, "L_2", "h (mesh size)")
		util.scheduleGnuplotFilename(err, err.plotPath.."all_"..errSuffix.."_h1_h.pdf", H1_h_Data, "H^1", "h (mesh size)")
	end
	
	if err.bUseExact then
		util.writeAndScheduleGnuplotDataType(
			err, discType, p, err.l2.exact.title, err.l2.exact.value, err.h1.exact.value, "exact")
	end	
	if err.bMaxLevel then
		util.writeAndScheduleGnuplotDataType(
			err, discType, p, err.l2.maxlevel.title, err.l2.maxlevel.value, err.h1.maxlevel.value, "maxlevel")
	end	
	if err.bPrevLevel then
		util.writeAndScheduleGnuplotDataType(
			err, discType, p, err.l2.prevlevel.title, err.l2.prevlevel.value, err.h1.prevlevel.value, "prevlevel")
	end	
end


function util.computeLinearStaticConvRatesForSpace(	err, dom, minLev, maxLev, discType, p,
													ConvRateSetup)	
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
	
	-- create Approximation Space
	print(">> Create ApproximationSpace: "..discType..", "..p)
	local approxSpace = ConvRateSetup.createApproxSpace(dom, discType, p)
	
	print(">> Create Domain Disc: "..discType..", "..p)
	local domainDisc = ConvRateSetup.createDomainDisc(discType, p, approxSpace)
	
	print(">> Create Solver")
	local solver = ConvRateSetup.createSolver(approxSpace)
	
	--------------------------------------------------------------------------------
	--  Apply Solver
	--------------------------------------------------------------------------------
	
	local u = {}

	write("\n>> Allocating storage for solution vectors.\n")
	for lev = minLev, maxLev do
		u[lev] = GridFunction(approxSpace, lev)
	end
	
	for lev = minLev, maxLev do
		write("\n>> Computing Level "..lev..".\n")

		-- set start values	
		if lev > minLev then	
			Prolongate(u[lev], u[lev-1]);
			write(">> Solution interpolated as start value for finer level.\n")
		else
			u[lev]:set(0.0)
			write(">> Solution set to zero on coarsest level.\n")	
		end		
		
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
		write(">> Linear Solver done.\n")
	end

	-- 5. compute error
	for lev = maxLev, minLev, -1 do
		write("\n>> Error Norm values on Level "..lev..".\n")
		
		quadOrder = p+3
		err.numDoFs[lev] 			= u[lev]:size()
		write(">> #DoF       on Level "..lev.." is "..err.numDoFs[lev] .."\n");

		-- w.r.t exact solution		
		if err.bUseExact then 
		err.l2.exact.value[lev] 	= L2Error(err.exactSol, u[lev], "c", 0.0, quadOrder)
		err.h1.exact.value[lev] 	= H1Error(err.exactSol, err.exactGrad, u[lev], "c", 0.0, quadOrder)
		write(">> L2 l-exact on Level "..lev.." is "..string.format("%.3e", err.l2.exact.value[lev]) .."\n");
		write(">> H1 l-exact on Level "..lev.." is "..string.format("%.3e", err.h1.exact.value[lev]) .."\n");
		end
		
		-- w.r.t max level solution
		if err.bMaxLevel then 
		err.l2.maxlevel.value[lev] 	= L2Error(u[maxLev], "c", u[lev], "c", quadOrder)
		err.h1.maxlevel.value[lev] 	= H1Error(u[maxLev], "c", u[lev], "c", quadOrder)
		write(">> L2 l-lmax  on Level "..lev.." is "..string.format("%.3e", err.l2.maxlevel.value[lev]) .."\n");
		write(">> H1 l-lmax  on Level "..lev.." is "..string.format("%.3e", err.h1.maxlevel.value[lev]) .."\n");
		end

		-- w.r.t previous level solution
		if err.bPrevLevel then 
		if lev > minLev then 
		err.l2.prevlevel.value[lev] = L2Error(u[lev], "c", u[lev-1], "c", quadOrder)
		err.h1.prevlevel.value[lev] = H1Error(u[lev], "c", u[lev-1], "c", quadOrder)
		write(">> L2 l-(l-1) on Level "..lev.." is "..string.format("%.3e", err.l2.prevlevel.value[lev]) .."\n");
		write(">> H1 l-(l-1) on Level "..lev.." is "..string.format("%.3e", err.h1.prevlevel.value[lev]) .."\n");
		else
		write(">> L2 l-(l-1) on Level "..lev.." is "..err.l2.prevlevel.value[lev] .."\n");	
		write(">> H1 l-(l-1) on Level "..lev.." is "..err.h1.prevlevel.value[lev].."\n");
		end
		end
		
		-- 6. write solution
		WriteGridFunctionToVTK(u[lev], ConvRateSetup.solPath.."sol_"..discType..p.."_l"..lev)
		write(">> Solution written to: "..ConvRateSetup.solPath.."sol_"..discType..p.."_l"..lev.."\n");	
	end
end


function util.computeLinearStaticConvRates(ConvRateSetup)
	
	-- create directories
	if ConvRateSetup.plotPath == nil then  ConvRateSetup.plotPath = "plots/" end
	if ConvRateSetup.solPath == nil then  ConvRateSetup.solPath = "sol/" end
	if ConvRateSetup.dataPath == nil then  ConvRateSetup.dataPath = "data/" end
	os.execute("mkdir " .. ConvRateSetup.dataPath)
	os.execute("mkdir " .. ConvRateSetup.plotPath)
	os.execute("mkdir " .. ConvRateSetup.plotPath.."single/")
	os.execute("mkdir " .. ConvRateSetup.solPath)

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

	-- compute problem
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

			err.h = {}
			for lev = minLev, maxLev do	err.h[lev] = MaxElementDiameter(dom, lev) end
		
			util.resetErrorRatesStorage(err, minLev, maxLev)
		
			util.computeLinearStaticConvRatesForSpace(err,
				 		dom, minLev, maxLev, discType, p, 
				 		ConvRateSetup)
	
			util.computeErrorRates(err)				 		

			-- write data to screen
			print("\n>> Statistic for type: "..discType..", order: "..p.."\n")
			util.updateErrorRatesScreenOutput(err)
			table.print(err.screen.values, {title = err.screen.title, format = err.screen.format, hline = true, vline = true})

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
