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
		
		errors[disc] = errors[disc] or {}
		
		for p = pmin, pmax do
		
			errors[disc][p] = errors[disc][p] or {}
			
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
	
			print(">> Create Time Disc: "..timeScheme..", "..orderOrTheta)
			local timeDisc = util.CreateTimeDisc(domainDisc, timeScheme, orderOrTheta)
			
			print(">> Create Solver")
			local newtonSolver = CreateSolver(approxSpace)

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
			local dts = {}
			local _dt = dtmax
			while _dt >= dtmin do
				dts[#dts+1] = _dt
				_dt = _dt * dtred
			end
			
			local u, uOld, TimeSeries, time, step = {}, {}, {}, {}, {}
			for lev = minLev, maxLev do
				u[lev], uOld[lev], TimeSeries[lev], time[lev], step[lev] = {}, {}, {}, {}, {}
				for k = 1, #dts do
					time[lev][k] = startTime
					step[lev][k] = 0
					u[lev][k] = GridFunction(approxSpace, lev)
					Interpolate(ExactSol["c"], u[lev][k], "c", time[lev][k])				
					TimeSeries[lev][k] = SolutionTimeSeries()
					TimeSeries[lev][k]:push(u[lev][k]:clone(), time[lev][k])
					
					if timeDisc:num_stages() > 1 then uOld[lev][k] = u[lev][k]:clone() end			
				end
			end
									
			--------------------------------------------------------------------------------
			--  Apply Solver
			--------------------------------------------------------------------------------
	
			local h, level, DoFs = {}, {}, {}
			for lev = minLev, maxLev do
				h[lev] = MaxElementDiameter(dom, lev)
				level[lev] = lev
			end	
	
			for lev = maxLev, minLev, -1 do
				for k = 1, #dts do

					DoFs[lev] = u[lev][k]:num_dofs()
	
					-- set order for bdf to 1 (initially)
					if timeScheme:lower() == "bdf" then timeDisc:set_order(1) end			
	
					while time[lev][k] < endTime do
			
						-- update step count
						step[lev][k] = step[lev][k] + 1
															
						-- time step size
						local dt = dts[k]
						local _dt = dt
												
						if time[lev][k] + dt > endTime then _dt = endTime - time[lev][k] end
						if ((endTime - (time[lev][k]+dt))/dt) < 1e-8 then _dt = endTime - time[lev][k] end
				
						-- get old solution if multistage
						if timeDisc:num_stages() > 1 then
							VecScaleAssign(u[lev][k], 1.0, uOld[lev][k])
						end			
	
						for stage = 1, timeDisc:num_stages() do
							timeDisc:set_stage(stage)
							timeDisc:prepare_step(TimeSeries[lev][k], _dt)
	
							-- solve step						
							newtonSolver:init(AssembledOperator(timeDisc, u[lev][k]:grid_level()))
							if newtonSolver:prepare(u[lev][k]) == false then print (">> Newton init failed."); exit(); end 
							if newtonSolver:apply(u[lev][k]) == false then print (">> Newton solver failed."); exit(); end 
											
							-- update new time
							time[lev][k] = timeDisc:future_time()
							if math.abs(endTime - time[lev][k]) < 1e-8*dt then time[lev][k] = endTime end
							
							-- push oldest solutions with new values to front, oldest sol pointer is poped from end	
							if timeScheme:lower() == "bdf" and step[lev][k] < orderOrTheta then
								print("++++++ BDF: Increasing order to "..step[lev][k]+1)
								timeDisc:set_order(step[lev][k]+1)
								solTimeSeries:push(u[lev][k]:clone(), time[lev][k])
							else 
								local oldestSol = TimeSeries[lev][k]:oldest()
								VecScaleAssign(oldestSol, 1.0, u[lev][k])
								TimeSeries[lev][k]:push_discard_oldest(oldestSol, time[lev][k])
							end
						end
			
						-- save this solution if multistage
						if timeDisc:num_stages() > 1 then
							uOld[lev][k] = u[lev][k]
						end			
					
						vtk = VTKOutput()
						vtk:print("Sol"..step[lev][k], u[lev][k])
					
						-- compute for each component
						local quadOrder = p+3
						local tp = time[lev][k]
						errors[disc][p][timeScheme] = errors[disc][p][timeScheme] or {}
						errors[disc][p][timeScheme][tp] = errors[disc][p][timeScheme][tp] or {}
						local err = errors[disc][p][timeScheme][tp]
						for f, Cmps in pairs(PlotCmps) do
		
							-- create component
							err[f] = err[f] or {}
									
							-- help fct to create an measurement
							local function createMeas(f, t, n, lev)
								err[f][t] = err[f][t] or {}
								err[f][t][n] = err[f][t][n] or {}
								err[f][t][n].value = err[f][t][n].value or {}						
								err[f][t][n].value[lev] = err[f][t][n].value[lev] or {}						
								return err[f][t][n].value[lev]
							end
												
							-- check for exact solution and grad
							local solAvail, gradAvail = true, true
							for _,cmp in pairs(Cmps) do
								if ExactSol  == nil or ExactSol[cmp]  == nil then solAvail  = false end
								if ExactGrad == nil or ExactGrad[cmp] == nil then gradAvail = false end
							end
																	
							-- w.r.t exact solution		
							if exact and solAvail then 					
								local value = createMeas(f, "l-exact", "l2", lev)
								value[dt] = 0.0
								for _,cmp in pairs(Cmps) do
									value[dt] = value[dt] + math.pow(L2Error(ExactSol[cmp], u[lev][k], cmp, tp, quadOrder), 2)
								end
								value[dt] = math.sqrt(value[dt])
								write(">> L2 l-exact for "..f.." on Level "..lev..", dt: "
										..dt..": "..string.format("%.3e", value[dt]) ..", at time "..tp.."\n");
		
								if gradAvail then 					
									local value = createMeas(f, "l-exact", "h1", lev)
									value[dt] = 0.0
									for _,cmp in pairs(Cmps) do
										value[dt] = value[dt] + math.pow(H1Error(ExactSol[cmp], ExactGrad[cmp], u[lev][k], cmp, tp, quadOrder), 2)
									end
									value[dt] = math.sqrt(value[dt])
									write(">> H1 l-exact for "..f.." on Level "..lev..", dt: "
											..dt..": "..string.format("%.3e", value[dt]) ..", at time "..tp.."\n");
								end
							end
						end
					end	-- end time interval
					
				end -- end timestep size
			end -- end level

			--------------------------------------------------------------------
			--  Compute Factors and Rates
			--------------------------------------------------------------------
	
			local err = errors[disc][p][timeScheme]
			for tp, _ in pairs(err) do
				for f, _ in pairs(PlotCmps) do
					for t, _ in pairs(err[tp][f]) do
						for n, _ in pairs(err[tp][f][t]) do
	
				local meas = err[tp][f][t][n]
	
				meas.fac = meas.fac or {}
				meas.rate = meas.rate or {}
				
				local value = meas.value
				local fac = meas.fac
				local rate = meas.rate
				value.h, value.dt = {}, {}
				fac.h, fac.dt = {}, {}
				rate.h, rate.dt = {}, {}
				
				-- rate in time
				for lev = minLev, maxLev do
					value.dt[lev], fac.dt[lev], rate.dt[lev] = {}, {}, {}
					for k, dt in ipairs(dts) do
						value.dt[lev][k] = value[lev][dt]
					end

					for k, dt in ipairs(dts) do
						if value.dt[lev][k] ~= nil and value.dt[lev][k-1] ~= nil then
							fac.dt[lev][k] = value.dt[lev][k-1]/value.dt[lev][k]
							rate.dt[lev][k] = math.log(fac.dt[lev][k]) / math.log(1/dtred)
						end
					end
				end

				-- rate in space
				for k, dt in ipairs(dts) do
					value.h[k], fac.h[k], rate.h[k] = {}, {}, {}
					
					for lev = minLev, maxLev do
						value.h[k][lev] = value[lev][dt]
					end
					
					for lev = minLev, maxLev do
						if value.h[k][lev] ~= nil and value.h[k][lev-1] ~= nil then
							fac.h[k][lev] = value.h[k][lev-1]/value.h[k][lev]
							rate.h[k][lev] = math.log(fac.h[k][lev]) / math.log(2)
						end
					end
				end
								
						end						
					end
				end
			end	

			--------------------------------------------------------------------
			--  Write Data to Screen
			--------------------------------------------------------------------
			
			local err = errors[disc][p][timeScheme]
			local tp = endTime
			for f, Cmps in pairs(PlotCmps) do
			
				write("\n>> Statistic for type: "..disc..", order: "..p..", comp: "..f.." [ ")			
				for _, cmp in pairs(Cmps) do write(cmp.." ") end
				print("]")
				
				for k, dt in ipairs(dts) do
					
					-- write data to screen
					print("\n>> Convergence in space: fixed dt = "..dt)

					local values = {level, h, DoFs}
					local heading = {"L", "h", "#DoFs"}
					local format = {"%d", "%.2e", "%d"}

					for t, _ in pairs(err[tp][f]) do
						for n, _ in pairs(err[tp][f][t]) do
							local meas = err[tp][f][t][n]
							table.append(values, {meas.value.h[k], meas.rate.h[k]}) 
							table.append(heading,{n.." "..t, "rate"})
							table.append(format, {"%.2e", "%.3f"})
						end
					end
											
					table.print(values, {heading = heading, format = format, 
										 hline = true, vline = true, forNil = "--"})										 										 
				end
				
				write("\n>> Statistic for type: "..disc..", order: "..p..", comp: "..f.." [ ")			
				for _, cmp in pairs(Cmps) do write(cmp.." ") end
				print("]")
				
				for lev = minLev, maxLev do
					
					-- write data to screen
					print("\nConvergence in time: fixed level "..lev..", h = "
							..string.format("%.3e", h[lev])..", #DoF = "..DoFs[lev])

					local values = {dts}
					local heading = {"dt"}
					local format = {"%.2e"}

					for t, _ in pairs(err[tp][f]) do
						for n, _ in pairs(err[tp][f][t]) do
							local meas = err[tp][f][t][n]
							table.append(values, {meas.value.dt[lev], meas.rate.dt[lev]}) 
							table.append(heading,{n.." "..t, "rate"})
							table.append(format, {"%.2e", "%.3f"})
						end
					end
											
					table.print(values, {heading = heading, format = format, 
										 hline = true, vline = true, forNil = "--"})
				end
				
			end
			
		end -- end p
	end -- end disc type
end
