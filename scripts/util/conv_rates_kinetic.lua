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

function util.rates.kinetic.NoMaxLevelPadding(p)
	return 0
end

function util.rates.kinetic.StdAutoStepSize(lev, h) 
	return h 
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

function util.rates.kinetic.StdLabel.SpaceLatex(x)
	local gpSpaceLabel ={ DoFs = "Anzahl Unbekannte",	h = "h"}
	return gpSpaceLabel[x]
end

function util.rates.kinetic.StdLabel.TimestepLatex()
	return "$\\Delta t$"
end

function util.rates.kinetic.StdLabel.TimeLatex()
	return "Zeit"
end

function util.rates.kinetic.StdLabel.NormLatex(f, t, n)
	local gpType = {	["exact"] = 	"{}",		
						["best"] = 		"h_{\text{min}}",
					}
	local gpNorm = 	{ l2 = "L_2",	h1 = "H^1"}
	
	return "$\\norm{"..f.."_h - "..f.."_{"..gpType[t].."} }_{ "..gpNorm[n].."}$"
end

function util.rates.kinetic.StdLabel.MeasPdfP(disc, p)
	return disc.." $P_"..p.."$"
end

function util.rates.kinetic.StdLabel.MeasPdfQ(disc, p)
	return disc.." $Q_"..p.."$"
end

function util.rates.kinetic.StdLabel.SpacePdf(x)
	local gpSpaceLabel ={ DoFs = "Unbekannte",	h = "h"}
	return gpSpaceLabel[x]
end

function util.rates.kinetic.StdLabel.TimestepPdf()
	return "dt"
end

function util.rates.kinetic.StdLabel.TimePdf()
	return "Zeit"
end

function util.rates.kinetic.StdLabel.NormPdf(f, t, n)
	local gpType = {	["exact"] = 	"{}",		
						["best"] = 		"h_{\text{min}}",
					}
	local gpNorm = 	{ l2 = "L_2",	h1 = "H^1"}
	
	return "$|| "..f.."_h - "..f.."_{"..gpType[t].."} ||_{"..gpNorm[n].."}$"
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

	local gpOptions = CRS.gpOptions or {}

	-- check for methods
	local CreateApproxSpace = 	CRS.CreateApproxSpace
	local CreateDomainDisc = 	CRS.CreateDomainDisc
	local CreateSolver = 		CRS.CreateSolver
	local CreateDomain = 		CRS.CreateDomain
	local SetStartSolution = 	CRS.SetStartSolution
	local MaxLevelPadding = 	CRS.MaxLevelPadding 	or util.rates.kinetic.StdMaxLevelPadding
	local AutoStepSize = 		CRS.AutoStepSize or util.rates.kinetic.StdAutoStepSize
	
	if 	CreateApproxSpace == nil or CreateDomainDisc == nil or 
		CreateSolver == nil or CreateDomain == nil or 
		SetStartSolution == nil then
		print("You must pass: CreateApproxSpace, CreateDomainDisc, CreateSolver, "
				.."CreateDomain, SetStartSolution")
		exit()
	end
	
	local SpaceDiscs = CRS.SpaceDiscs
	local TimeDiscs = CRS.TimeDiscs

	local best = CRS.best; 				if best == nil then best = true end
	local exact = CRS.exact; 			if exact == nil then exact = true end
	local plotSol = CRS.plotSol; 		if plotSol == nil then plotSol = false end
	local onlyLast = CRS.onlyLast; 		if onlyLast == nil then onlyLast = true end
	local plotOnlyTime = CRS.plotOnlyTime; if plotOnlyTime == nil then plotOnlyTime = false end

	local ExactSol = CRS.ExactSol 
 	local ExactGrad = CRS.ExactGrad 
	local PlotCmps = CRS.PlotCmps
	
	local MeasLabel = 	CRS.MeasLabel or util.rates.kinetic.StdLabel.MeasLatexP
	local SpaceLabel = 	CRS.SpaceLabel 	or util.rates.kinetic.StdLabel.SpaceLatex
	local TimeLabel = 	CRS.TimeLabel 	or util.rates.kinetic.StdLabel.TimeLatex
	local TimestepLabel = 	CRS.TimestepLabel 	or util.rates.kinetic.StdLabel.TimestepLatex
	local NormLabel = 	CRS.NormLabel 	or util.rates.kinetic.StdLabel.NormLatex
	
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

	local ts = nil
	
	local StartTime 	= CRS.StartTime or 0.0
	local EndTime		= CRS.EndTime	or 1.0

	--------------------------------------------------------------------
	--  Loop Spacial discs
	--------------------------------------------------------------------
	for _, SpaceDisc in ipairs(SpaceDiscs) do

		local disc 			= SpaceDisc.type
		if disc == nil then print("type required."); exit(); end

		local pmin 			= SpaceDisc.pmin or 1
		local pmax 			= SpaceDisc.pmax or 1
		local lmin 			= SpaceDisc.lmin or 0
		local lmax 			= SpaceDisc.lmax or numRefs
		
		if lmin > lmax then print("lmin: "..lmin.." must be less or equal lmax: "..lmax); exit(); end
		if lmax > numRefs then print("lmax: "..lmax.." must be less or equal numRefs: "..numRefs); exit(); end
				
		errors[disc] = errors[disc] or {}
		
		for p = pmin, pmax do
		
			errors[disc][p] = errors[disc][p] or {}
			
			local maxLev = lmax - MaxLevelPadding(p)
			local minLev = lmin

			--------------------------------------------------------------------
			--  Create ApproxSpace, Disc and Solver
			--------------------------------------------------------------------
				
			print(">> Create ApproximationSpace: "..disc..", "..p)
			local approxSpace = CreateApproxSpace(dom, disc, p)
			
			print(">> Create Domain Disc: "..disc..", "..p)
			local domainDisc = CreateDomainDisc(approxSpace, disc, p)
			
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
			
			local AutoTimeScheme = false
			if not TimeDiscs then AutoTimeScheme = true end
			local UsedTimeDiscs = TimeDiscs or {{type = "alexander", orderOrTheta = p, dt = 1, sub = 2, refs = 0}}
			
			--------------------------------------------------------------------
			--  Loop Time discs
			--------------------------------------------------------------------
			for _, TimeDisc in ipairs(UsedTimeDiscs) do

				local timeScheme	= TimeDisc.type
				local orderOrTheta 	= TimeDisc.orderOrTheta or 1
				if timeScheme == nil then print("type required."); exit(); end
				local ts = timeScheme..orderOrTheta

				local dt	 		= TimeDisc.dt		
				local refs	 		= TimeDisc.refs
				local sub	 		= TimeDisc.sub or 2
				if dt == nil then print("dt required."); exit(); end
				if refs == nil then print("refs required."); exit(); end
				
				----------------------------------------------------------------
				--  Print Setup
				----------------------------------------------------------------
				
				print("\n")
				print("---------------------------")
				print(" General parameters chosen:")
				print("---------------------------")
				print("    dim        = " .. dim)
				print("    grid       = " .. gridName)
				print("    maxLev     = " .. maxLev)
				print("    minLev     = " .. minLev)
				print("    dt         = " .. dt)
				print("    sub        = " .. sub)
				print("    refs       = " .. refs)
				print("    TimeDisc   = " .. timeScheme)
				print("    order      = " .. orderOrTheta)
				print("    StartTime  = " .. StartTime)
				print("    EndTime    = " .. EndTime)
				print("    SpaceDisc  = " .. disc)
				print("    order      = " .. p)
				print("\n")
	
				----------------------------------------------------------------
				--  Create Solutions on each level
				----------------------------------------------------------------
		
				print(">> Create Time Disc: "..timeScheme..", "..orderOrTheta)
				local timeDisc = util.CreateTimeDisc(domainDisc, timeScheme, orderOrTheta)
				
				write(">> Allocating storage for solution vectors.\n")
				local memory = {}
				for lev = minLev, maxLev do
					memory[lev] = {}
					for k = 0, refs do
						memory[lev][k] = {}
						local mem = memory[lev][k]
						if k == 0 then mem.dt = dt;
						else mem.dt = memory[lev][k-1].dt * (1/sub) end
												
						mem.time = StartTime
						mem.step = 0
						mem.u = GridFunction(approxSpace, lev)
						mem.TimeSeries = SolutionTimeSeries()
						SetStartSolution(mem.u, mem.time)				
						mem.TimeSeries:push(mem.u:clone(), mem.time)
					end
				end
										
				----------------------------------------------------------------
				--  Statistics
				----------------------------------------------------------------
		
				-- prepare error measurement		
				errors[disc][p][ts] = errors[disc][p][ts] or {}	
				local err = errors[disc][p][ts]
				
				err.h, err.DoFs, err.level, err.dt, err.time = {}, {}, {}, {}, {}
				for lev = minLev, maxLev do	
					err.h[lev] = MaxElementDiameter(dom, lev) 
					err.level[lev] = lev
					err.DoFs[lev] = memory[lev][0].u:num_dofs()
				end	
					
				if AutoTimeScheme then
					for lev = minLev, maxLev do	
						memory[lev][0].dt = AutoStepSize(lev, err.h[lev])
					end
				end
					
				----------------------------------------------------------------
				--  Compute solution for all level and time step sizes
				----------------------------------------------------------------
				
				-- loop time interval
				local sliceTime = StartTime
				local tp = 0
				while sliceTime < EndTime do
				
					sliceTime = sliceTime + dt
					if sliceTime > EndTime then sliceTime = EndTime end
					if (EndTime-sliceTime)/sliceTime < 1e-8 then sliceTime = EndTime end
					print(">> >>>>> Advancing to time "..sliceTime)
					tp = tp +1
					err.time[tp] = sliceTime

					-- advance all discs to end of global time step
					for lev = maxLev, minLev, -1 do
						for k = refs, 0, -1 do
		
							local mem = memory[lev][k]
		
							err.dt[k] = mem.dt
					
							write(">> dt "..mem.dt.." on lev "..lev..": ")
							while mem.time < sliceTime do
								write(".")
								-- update step count
								mem.step = mem.step + 1
																	
								-- time step size
								local dt = mem.dt
								local dodt = dt

								local usedTimeDisc = timeDisc
								if ts:sub(1,3):lower() == "bdf" then
									timeDisc:set_order(mem.TimeSeries:size())
									if mem.TimeSeries:size() < orderOrTheta then
										usedTimeDisc = util.CreateTimeDisc(domainDisc, "sdirk", orderOrTheta)
										dodt =  dodt
										print(" BDF("..orderOrTheta.."): using  SDIRK("..
												orderOrTheta..") for step "..mem.TimeSeries:size())
									end
								end
														
								if mem.time + dodt > sliceTime then dodt = sliceTime - mem.time end
								if ((sliceTime - (mem.time+dodt))/dodt) < 1e-8 then dodt = sliceTime - mem.time end
						
								for stage = 1, usedTimeDisc:num_stages() do
									usedTimeDisc:set_stage(stage)
									usedTimeDisc:prepare_step(mem.TimeSeries, dodt)
			
									-- solve step						
									newtonSolver:init(AssembledOperator(usedTimeDisc, mem.u:grid_level()))
									if newtonSolver:prepare(mem.u) == false then print (">> Newton init failed."); exit(); end 
									if newtonSolver:apply(mem.u) == false then print (">> Newton solver failed."); exit(); end 
									AdjustMeanValue(mem.u, "p")

									if plotSol then
										vtk = VTKOutput()
										vtk:print(solPath.."Sol_"..ts.."stage_"..stage.."_s"..mem.step, mem.u)
									end
																										
									-- update new time
									mem.time = usedTimeDisc:future_time()
									if math.abs(sliceTime - mem.time) < 1e-8*dt then mem.time = sliceTime end

									if plotSol then
										vtk = VTKOutput()
										vtk:print(solPath.."Sol_"..ts.."_stage"..stage.."_s"..mem.step, mem.u)
									end
									
									-- push oldest solutions with new values to front, oldest sol pointer is poped from end	
									if ts:sub(1,3):lower() == "bdf" and mem.TimeSeries:size() < orderOrTheta and stage == usedTimeDisc:num_stages() then
										--Interpolate(ExactSol["c"], mem.u, "c", mem.time)				
									
										print("++++++ BDF: Increasing order to "..mem.step+1)
										mem.TimeSeries:push(mem.u:clone(), mem.time)
									else 
										if not (ts:sub(1,3):lower() == "bdf" and stage ~= usedTimeDisc:num_stages()) then
											local oldestSol = mem.TimeSeries:oldest()
											VecScaleAssign(oldestSol, 1.0, mem.u)
											mem.TimeSeries:push_discard_oldest(oldestSol, mem.time)
										end
									end
								end
					
							end	-- end slice interval
							write("\n")
							
							if plotSol then
								vtk = VTKOutput()
								vtk:print(solPath.."Sol_"..ts.."_s"..mem.step, mem.u)
							end
						
							---------------------------------------------------- 
							-- compute norms at slice point (start)
							---------------------------------------------------- 
							local quadOrder = p+3
							
							for f, Cmps in pairs(PlotCmps) do
			
								-- create component
								err[f] = err[f] or {}
										
								-- help fct to create an measurement
								local function createMeas(f, t, n, lev)
									err[f][t] = err[f][t] or {}
									err[f][t][n] = err[f][t][n] or {}
									err[f][t][n][tp] = err[f][t][n][tp] or {}
									err[f][t][n][tp].value = err[f][t][n][tp].value or {}						
									err[f][t][n][tp].value[lev] = err[f][t][n][tp].value[lev] or {}						
									return err[f][t][n][tp].value[lev]
								end
													
								-- check for exact solution and grad
								local solAvail, gradAvail = true, true
								for _,cmp in pairs(Cmps) do
									if ExactSol  == nil or ExactSol[cmp]  == nil then solAvail  = false end
									if ExactGrad == nil or ExactGrad[cmp] == nil then gradAvail = false end
								end
																		
								-- w.r.t exact solution		
								if exact and solAvail then 					
									local value = createMeas(f, "exact", "l2", lev)
									value[k] = 0.0
									for _,cmp in pairs(Cmps) do
										value[k] = value[k] + math.pow(L2Error(ExactSol[cmp], mem.u, cmp, mem.time, quadOrder), 2)
									end
									value[k] = math.sqrt(value[k])
									write(">> L2 l-exact for "..f.." on Level "..lev..", dt: "
											..mem.dt..": "..string.format("%.3e", value[k]) ..", at time "..mem.time.."\n");
			
									if gradAvail then 					
										local value = createMeas(f, "exact", "h1", lev)
										value[k] = 0.0
										for _,cmp in pairs(Cmps) do
											value[k] = value[k] + math.pow(H1Error(ExactSol[cmp], ExactGrad[cmp], mem.u, cmp, mem.time, quadOrder), 2)
										end
										value[k] = math.sqrt(value[k])
										write(">> H1 l-exact for "..f.." on Level "..lev..", dt: "
												..mem.dt..": "..string.format("%.3e", value[k]) ..", at time "..mem.time.."\n");
									end
								end
								
								-- w.r.t max level solution
								if best and lev < maxLev and k < refs then
									local value = createMeas(f, "best", "l2", lev)
									value[k] = 0.0
									for _,cmp in pairs(Cmps) do
										value[k] = value[k] + math.pow(L2Error(memory[maxLev][refs].u, cmp, mem.u, cmp, quadOrder), 2)
									end
									value[k] = math.sqrt(value[k])
									write(">> L2 l-lmax  for "..f.." on Level "..lev..", dt: "
											..mem.dt..": "..string.format("%.3e", value[k]) ..", at time "..mem.time.."\n");
			
									local value = createMeas(f, "best", "h1", lev)
									value[k] = 0.0
									for _,cmp in pairs(Cmps) do
										value[k] = value[k] + math.pow(H1Error(memory[maxLev][refs].u, cmp, mem.u, cmp, quadOrder), 2)
									end
									value[k] = math.sqrt(value[k])
									write(">> H1 l-lmax  for "..f.." on Level "..lev..", dt: "
											..mem.dt..": "..string.format("%.3e", value[k]) ..", at time "..mem.time.."\n");
								end
				
								
							end
							---------------------------------------------------- 
							-- compute norms at slice point (end)
							---------------------------------------------------- 
							
						end -- end timestep size
					end -- end level

					--------------------------------------------------------------------
					--  Compute Factors and Rates
					--------------------------------------------------------------------
					for f, _ in pairs(PlotCmps) do
						for t, _ in pairs(err[f]) do
							for n, _ in pairs(err[f][t]) do
			
						local meas = err[f][t][n][tp]
			
						meas.fac = meas.fac or {}
						meas.rate = meas.rate or {}
						
						local value = meas.value
						local fac = meas.fac
						local rate = meas.rate
						value.h, value.dt = {}, {}
						fac.h, fac.dt = {}, {}
						rate.h, rate.dt = {}, {}
						
						-- rate in time
						for lev, _ in iipairs(value) do
							value.dt[lev], fac.dt[lev], rate.dt[lev] = {}, {}, {}
							for k, _ in iipairs(value[lev]) do
								value.h[k] = value.h[k] or {}
								fac.h[k] = fac.h[k] or {} 
								rate.h[k] = rate.h[k] or {}
								value.dt[lev][k] = value[lev][k]
								value.h[k][lev] = value[lev][k]
							end
						end
		
						for lev, _ in iipairs(value.dt) do
							for k, _ in iipairs(value.dt[lev]) do
								if value.dt[lev][k] ~= nil and value.dt[lev][k-1] ~= nil then
									fac.dt[lev][k] = value.dt[lev][k-1]/value.dt[lev][k]
									rate.dt[lev][k] = math.log(fac.dt[lev][k]) / math.log(sub)
								end
							end
						end
		
						-- rate in space
						for k, _ in iipairs(value.h) do
							for lev, _ in iipairs(value.h[k]) do
								if value.h[k][lev] ~= nil and value.h[k][lev-1] ~= nil then
									fac.h[k][lev] = value.h[k][lev-1]/value.h[k][lev]
									rate.h[k][lev] = math.log(fac.h[k][lev]) / math.log(2)
								end
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
						print("] at time "..err.time[tp])
						
						if not plotOnlyTime then
						for k, dt in iipairs(err.dt) do
							
							-- write data to screen
							print("\n>> Convergence in space: fixed dt = "..dt)
		
							local values = {err.level, err.h, err.DoFs}
							local heading = {"L", "h", "#DoFs"}
							local format = {"%d", "%.2e", "%d"}
		
							for t, _ in pairs(err[f]) do
								for n, _ in pairs(err[f][t]) do
									local meas = err[f][t][n][tp]
									table.append(values, {meas.value.h[k], meas.rate.h[k]}) 
									table.append(heading,{n.." "..t, "rate"})
									table.append(format, {"%.2e", "%.3f"})
								end
							end
													
							table.print(values, {heading = heading, format = format, 
												 hline = true, vline = true, forNil = "--"})										 										 
						end
						end
						
						write("\n>> Statistic for type: "..disc..", order: "..p..", comp: "..f.." [ ")			
						for _, cmp in pairs(Cmps) do write(cmp.." ") end
						print("] at time "..err.time[tp])
						
						for lev, _ in iipairs(err.h) do
							
							-- write data to screen
							print("\nConvergence in time: fixed level "..lev..", h = "
									..string.format("%.3e", err.h[lev])..", #DoF = "..err.DoFs[lev])
		
							local values = {err.dt}
							local heading = {"dt"}
							local format = {"%.2e"}
		
							for t, _ in pairs(err[f]) do
								for n, _ in pairs(err[f][t]) do
									local meas = err[f][t][n][tp]
									table.append(values, {meas.value.dt[lev], meas.rate.dt[lev]}) 
									table.append(heading,{n.." "..t, "rate"})
									table.append(format, {"%.2e", "%.3f"})
								end
							end
													
							table.print(values, {heading = heading, format = format, 
												 hline = true, vline = true, forNil = "--"})
						end
					end
					
				end -- end slice time	
			end -- end time disc type 
		end -- end p
	end -- end space disc type
	
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
	
	local function addSet(plot, dataset, label)
		table.insert( plot, dataset)			
		plot.label = label
	end

	local style3d = "lines"

	write(">> Writing measured rates to data files ...")
	for disc, _ in pairs(errors) do
		for p, _ in pairs(errors[disc]) do
			for ts, _ in pairs(errors[disc][p]) do
				for f, _ in pairs(PlotCmps) do
					for t, _ in pairs(errors[disc][p][ts][f]) do
						for n, _ in pairs(errors[disc][p][ts][f][t]) do

		---------------------------------------------
		-- plots for single timestep
		---------------------------------------------
		local err = errors[disc][p][ts]
		for tp, _ in pairs(errors[disc][p][ts][f][t][n]) do
	
			if onlyLast then tp = #errors[disc][p][ts][f][t][n] end
			
			local dir = dataPath..tp.."/"
			ensureDir(dir)		
	
			local plotDir = plotPath..tp.."/"
			ensureDir(plotDir)				
			gpData["dirs"] = gpData["dirs"] or {}
			table.insert(gpData["dirs"], plotDir)
						
			local value = errors[disc][p][ts][f][t][n][tp].value
	
			-- convergence in time
			for lev, _ in iipairs(value.dt) do		
				local levID = "lev"..lev
				local file = dir..table.concat({"error",disc,p,ts,f,t,n,levID},"_")..".dat"
				local cols = {err.dt, value.dt[lev]}
				gnuplot.write_data(file, cols)
				
				local dataset = {label=ts, file=file, style="linespoints", 1, 2}
				local label = { x = TimestepLabel(), y = NormLabel(f,t,n)}
				addSet( accessPlot(disc, p, ts, tp, f, t, n), dataset, label)
				addSet( accessPlot(disc, p, "all", tp, f, t, n), dataset, label)

				-- single dataset			
				local file = plotDir..table.concat({disc,p,ts,f,"__",t,n,"dt"},"_")
				gpData[file] = getPlot(disc, p, ts, tp, f, t, n)
				gpData[file].gpOptions = {logscale = {x = true, y = true}}
				
				-- grouping timediscs			
				local file = plotDir..table.concat({disc,p,"all",f,"__",t,n,"dt"},"_")
				gpData[file] = getPlot(disc, p, "all", tp, f, t, n)
				gpData[file].gpOptions = {logscale = {x = true, y = true}}
			end
				
			-- convergence in space
			if not plotOnlyTime then
			for k, _ in iipairs(value.h) do						
				local dtID = "dt"..k.."["..err.dt[k].."]"
				local file = dir..table.concat({"error",disc,p,ts,f,t,n,dtID},"_")..".dat"
				local cols = {err.DoFs, err.h, value.h[k]}
				gnuplot.write_data(file, cols)		
				
				for xCol, x in ipairs({"DoFs", "h"}) do
					local dataset = {label=MeasLabel(disc, p), file=file, style="linespoints", xCol, 3}
					local label = { x = SpaceLabel(x), y = NormLabel(f,t,n)}
					addSet( accessPlot(disc, p, "all", tp, f, t, n, x), dataset, label)
					addSet( accessPlot(disc,    "all", tp, f, t, n, x), dataset, label)
					addSet( accessPlot("all",   "all", tp, f, t, n, x), dataset, label)

					-- single dataset			
					local file = plotDir..table.concat({disc,p,"all",f,"__",t,n,x},"_")
					gpData[file] = getPlot(disc, p, "all", tp, f, t, n, x)
					gpData[file].gpOptions = {logscale = {x = true, y = true}}
					
					-- grouping space disc			
					local file = plotDir..table.concat({disc,"all",f,"__",t,n,x},"_")
					gpData[file] = getPlot(disc,    "all", tp, f, t, n, x)
					gpData[file].gpOptions = {logscale = {x = true, y = true}}
	
					-- grouping space disc	+ p		
					local file = plotDir..table.concat({"all","all",f,"__",t,n,x},"_")
					gpData[file] = getPlot("all",   "all", tp, f, t, n, x)
					gpData[file].gpOptions = {logscale = {x = true, y = true}}
				end
			end

			-- convergence in space-time
			local file = dir..table.concat({"error",disc,p,ts,f,t,n},"_")..".dat"
			local fio = io.open(file, "w+")
			fio:close()
			for k, _ in iipairs(value.h) do						
				local dt = err.dt[k]
				local dts = {}; for i,_ in pairs(value.h[k]) do dts[i] = dt end
				local cols = {err.DoFs, err.h, dts, value.h[k]}
				gnuplot.write_data(file, cols, false, "a+")
				local fio = io.open(file, "a+")
				fio:write("\n")
				fio:close()
			end
			
			for xCol, x in ipairs({"DoFs", "h"}) do
				local dataset = {label=MeasLabel(disc, p), file=file, style=style3d, xCol, 3, 4}
				local label = { x = SpaceLabel(x), y = TimestepLabel(), z = NormLabel(f,t,n)}
								
				addSet( accessPlot(disc, p, ts, tp, f, t, n, x), dataset, label)
				addSet( accessPlot(disc,    ts, tp, f, t, n, x), dataset, label)
				addSet( accessPlot("all",   ts, tp, f, t, n, x), dataset, label)
							
				-- single dataset			
				local file = plotDir..table.concat({disc,p,ts,f,"__",t,n,"dt",x},"_")
				gpData[file] = getPlot(disc, p, ts, tp, f, t, n, x)
				gpData[file].gpOptions = {logscale = {x = true, y = true, z = true}}
		
				-- grouping by (disc+p)								
				local file = plotDir..table.concat({f,disc,ts,"__",t,n,"dt",x}, "_")	
				gpData[file] = getPlot(disc, ts, tp, f, t, n, x)		
				gpData[file].gpOptions = {logscale = {x = true, y = true, z = true}}
		
				-- grouping (all discs+p)
				local file = plotDir..table.concat({f,"all",ts,"__",t,n,"dt",x}, "_")	
				gpData[file] = getPlot("all", ts, tp, f, t, n, x)	
				gpData[file].gpOptions = {logscale = {x = true, y = true, z = true}}					
			end
			end
			
			if onlyLast then break end
		end
		
		---------------------------------------------
		-- plots for time series
		---------------------------------------------
		if not plotOnlyTime then

		-- 2d: error over time
		local value = errors[disc][p][ts][f][t][n][1].value
		for lev, _ in iipairs(value) do		
			for k, _ in iipairs(value[lev]) do		
				local levID = "lev"..lev
				local dtID = "dt"..k.."["..err.dt[k].."]"
			
				local name = table.concat({disc,p,ts,f,t,n,"time",levID,dtID},"_")
				local file = dataPath.."error_"..name..".dat"

				local fio = io.open(file, "w+")
				for tp, _ in pairs(errors[disc][p][ts][f][t][n]) do
			
					local value = errors[disc][p][ts][f][t][n][tp].value

					if value[lev][k] then
						fio:write(err.time[tp])
						fio:write(" "..value[lev][k])	
						fio:write("\n")					
					end
				end				
				fio:close()
				
				local levID = "lev"..lev
				local dtID = "dt"..k.."["..err.dt[k].."]"

				local dataset = {label = MeasLabel(disc, p), file=file, style="linespoints", 1, 2}
				local label = { x = TimeLabel(), y = NormLabel(f,t,n)}
				addSet( accessPlot(disc, p, ts, f, t, n, "time", levID, dtID), dataset, label)
				addSet( accessPlot(disc,    ts, f, t, n, "time", levID, dtID), dataset, label)
				addSet( accessPlot("all",   ts, f, t, n, "time", levID, dtID), dataset, label)

				local dataset = {label = MeasLabel(disc, p)..", dt: "..err.dt[k], file=file, style="linespoints", 1, 2}
				local label = { x = TimeLabel(), y = NormLabel(f,t,n)}
				addSet( accessPlot(disc, p, ts, f, t, n, "time", levID), dataset, label)
				addSet( accessPlot(disc,    ts, f, t, n, "time", levID), dataset, label)
				addSet( accessPlot("all",   ts, f, t, n, "time", levID), dataset, label)
				
				local dataset = {label = MeasLabel(disc, p)..", lev "..lev, file=file, style="linespoints", 1, 2}
				local label = { x = TimeLabel(), y = NormLabel(f,t,n)}
				addSet( accessPlot(disc, p, ts, f, t, n, "time", dtID), dataset, label)
				addSet( accessPlot(disc,    ts, f, t, n, "time", dtID), dataset, label)
				addSet( accessPlot("all",   ts, f, t, n, "time", dtID), dataset, label)

				local dataset = {label = MeasLabel(disc, p)..", lev "..lev..", dt: "..err.dt[k], file=file, style="linespoints", 1, 2}
				local label = { x = TimeLabel(), y = NormLabel(f,t,n)}
				addSet( accessPlot(disc, p, ts, f, t, n, "time"), dataset, label)
				addSet( accessPlot(disc,    ts, f, t, n, "time"), dataset, label)
				addSet( accessPlot("all",   ts, f, t, n, "time"), dataset, label)

				-- single
				local file = plotPath..table.concat({disc,p,ts,f,"__",t,n,"time","__",levID,dtID},"_")
				gpData[file] = getPlot(disc, p, ts, f, t, n, "time", levID, dtID)
				gpData[file].gpOptions = {logscale = false}
								
				-- single (grouping dts)
				local file = plotPath..table.concat({disc,p,ts,f,"__",t,n,"time","__",levID},"_")
				gpData[file] = getPlot(disc, p, ts, f, t, n, "time", levID)
				gpData[file].gpOptions = {logscale = false}

				-- single (grouping lev)
				local file = plotPath..table.concat({disc,p,ts,f,"__",t,n,"time","__",dtID},"_")
				gpData[file] = getPlot(disc, p, ts, f, t, n, "time", dtID)
				gpData[file].gpOptions = {logscale = false}

				-- single (grouping dts+lev)
				local file = plotPath..table.concat({disc,p,ts,f,"__",t,n,"time"},"_")
				gpData[file] = getPlot(disc, p, ts, f, t, n, "time")
				gpData[file].gpOptions = {logscale = false}


				-- grouping by (disc+p)
				local file = plotPath..table.concat({disc,  ts,f,"__",t,n,"time","__",levID,dtID},"_")
				gpData[file] = getPlot(disc,    ts, f, t, n, "time", levID, dtID)
				gpData[file].gpOptions = {logscale = false}
								
				-- grouping by (disc+p) (grouping dts)
				local file = plotPath..table.concat({disc,  ts,f,"__",t,n,"time","__",levID},"_")
				gpData[file] = getPlot(disc,    ts, f, t, n, "time", levID)
				gpData[file].gpOptions = {logscale = false}

				-- grouping by (disc+p) (grouping lev)
				local file = plotPath..table.concat({disc,  ts,f,"__",t,n,"time","__",dtID},"_")
				gpData[file] = getPlot(disc,    ts, f, t, n, "time", dtID)
				gpData[file].gpOptions = {logscale = false}

				-- grouping by (disc+p) (grouping dts+lev)
				local file = plotPath..table.concat({disc,  ts,f,"__",t,n,"time"},"_")
				gpData[file] = getPlot(disc,    ts, f, t, n, "time")
				gpData[file].gpOptions = {logscale = false}


				-- grouping (all discs+p)
				local file = plotPath..table.concat({"all",  ts,f,"__",t,n,"time","__",levID,dtID},"_")
				gpData[file] = getPlot("all",    ts, f, t, n, "time", levID, dtID)
				gpData[file].gpOptions = {logscale = false}
								
				-- grouping (all discs+p) (grouping dts)
				local file = plotPath..table.concat({"all",  ts,f,"__",t,n,"time","__",levID},"_")
				gpData[file] = getPlot("all",    ts, f, t, n, "time", levID)
				gpData[file].gpOptions = {logscale = false}

				-- grouping (all discs+p) (grouping lev)
				local file = plotPath..table.concat({"all",  ts,f,"__",t,n,"time","__",dtID},"_")
				gpData[file] = getPlot("all",    ts, f, t, n, "time", dtID)
				gpData[file].gpOptions = {logscale = false}

				-- grouping (all discs+p) (grouping dts+lev)
				local file = plotPath..table.concat({"all",  ts,f,"__",t,n,"time"},"_")
				gpData[file] = getPlot("all",    ts, f, t, n, "time")
				gpData[file].gpOptions = {logscale = false}
			end
		end

		-- 3d: error over time and space
		local value = errors[disc][p][ts][f][t][n][1].value
		for k, _ in iipairs(value[table.imin(value)]) do		
			local dtID = "dt"..k.."["..err.dt[k].."]"
			local name = table.concat({disc,p,ts,f,t,n,"time","space",dtID},"_")
			local file = dataPath.."error_"..name..".dat"

			local fio = io.open(file, "w+")
			for lev, _ in iipairs(value) do		
				for tp, _ in pairs(errors[disc][p][ts][f][t][n]) do
			
					local value = errors[disc][p][ts][f][t][n][tp].value
		
					if value[lev][k] then
						fio:write(err.DoFs[lev].." "..err.h[lev])
						fio:write(" "..err.time[tp])
						fio:write(" "..value[lev][k])	
						fio:write("\n")					
					end
				end			
				fio:write("\n")
			end
			fio:close()

			for xCol, x in ipairs({"DoFs", "h"}) do
				local dataset = {label = MeasLabel(disc, p), file=file, style=style3d, xCol, 3, 4}
				local label = { x = SpaceLabel(x), y = TimeLabel(), z = NormLabel(f,t,n)}
				addSet( accessPlot(disc, p, ts, f, t, n, "time", x, dtID), dataset, label)
				addSet( accessPlot(disc,    ts, f, t, n, "time", x, dtID), dataset, label)
				addSet( accessPlot("all",   ts, f, t, n, "time", x, dtID), dataset, label)
				
				local dataset = {label = MeasLabel(disc, p)..", dt: "..err.dt[k], file=file, style="linespoints", xCol, 3, 4}
				local label = { x = SpaceLabel(x), y = TimeLabel(), z = NormLabel(f,t,n)}
				addSet( accessPlot(disc, p, ts, f, t, n, "time", x), dataset, label)
				addSet( accessPlot(disc,    ts, f, t, n, "time", x), dataset, label)
				addSet( accessPlot("all",   ts, f, t, n, "time", x), dataset, label)
				
				-- single
				local file = plotPath..table.concat({disc,p,ts,f,"__",t,n,"time",x,"__",dtID},"_")
				gpData[file] = getPlot(disc, p, ts, f, t, n, "time", x, dtID)
				gpData[file].gpOptions = {logscale = {x = true, y = false, z = true}}
								
				-- single (grouping dts)
				local file = plotPath..table.concat({disc,p,ts,f,"__",t,n,"time",x},"_")
				gpData[file] = getPlot(disc, p, ts, f, t, n, "time", x)
				gpData[file].gpOptions = {logscale = {x = true, y = false, z = true}}
				
				-- grouping by (disc+p)								
				local file = plotPath..table.concat({disc,  ts,f,"__",t,n,"time",x,"__",dtID},"_")	
				gpData[file] = getPlot(disc,    ts, f, t, n, "time", x, dtID)	
				gpData[file].gpOptions = {logscale = {x = true, y = false, z = true}}

				-- grouping by (disc+p)	(grouping dts)	 						
				local file = plotPath..table.concat({disc,  ts,f,"__",t,n,"time",x},"_")	
				gpData[file] = getPlot(disc,    ts, f, t, n, "time", x)	
				gpData[file].gpOptions = {logscale = {x = true, y = false, z = true}}
		
				-- grouping (all discs+p)
				local file = plotPath..table.concat({"all",  ts,f,"__",t,n,"time",x,"__",dtID},"_")		
				gpData[file] = getPlot("all",    ts, f, t, n, "time", x, dtID)	
				gpData[file].gpOptions = {logscale = {x = true, y = false, z = true}}					

				-- grouping (all discs+p) (grouping dts)	 
				local file = plotPath..table.concat({"all",  ts,f,"__",t,n,"time",x},"_")		
				gpData[file] = getPlot("all",    ts, f, t, n, "time", x)	
				gpData[file].gpOptions = {logscale = {x = true, y = false, z = true}}					
				
			end
		end

		-- 3d: error over time and dt
		local value = errors[disc][p][ts][f][t][n][1].value
		for lev, _ in iipairs(value) do		
			local levID = "lev"..lev
			local name = table.concat({disc,p,ts,f,t,n,"time","dt",levID},"_")
			local file = dataPath.."error_"..name..".dat"

			local fio = io.open(file, "w+")
			for k, _ in iipairs(value[table.imin(value)]) do		
				for tp, _ in pairs(errors[disc][p][ts][f][t][n]) do
			
					local value = errors[disc][p][ts][f][t][n][tp].value
		
					if value[lev][k] then
						fio:write(" "..err.dt[k])
						fio:write(" "..err.time[tp])
						fio:write(" "..value[lev][k])	
						fio:write("\n")					
					end
				end			
				fio:write("\n")
			end
			fio:close()

			local dataset = {label = MeasLabel(disc, p), file=file, style=style3d, 1, 2, 3}
			local label = { x = TimestepLabel(), y = TimeLabel(), z = NormLabel(f,t,n)}
			addSet( accessPlot(disc, p, ts, f, t, n, "time", "dt",levID), dataset, label)
			addSet( accessPlot(disc,    ts, f, t, n, "time", "dt",levID), dataset, label)
			addSet( accessPlot("all",   ts, f, t, n, "time", "dt",levID), dataset, label)

			local dataset = {label = MeasLabel(disc, p)..", lev "..lev, file=file, style="linespoints", 1, 2, 3}
			local label = { x = TimestepLabel(), y = TimeLabel(), z = NormLabel(f,t,n)}
			addSet( accessPlot(disc, p, ts, f, t, n, "time", "dt"), dataset, label)
			addSet( accessPlot(disc,    ts, f, t, n, "time", "dt"), dataset, label)
			addSet( accessPlot("all",   ts, f, t, n, "time", "dt"), dataset, label)

			-- single
			local file = plotPath..table.concat({disc,p,ts,f,"__",t,n,"time","dt","__",levID},"_")
			gpData[file] = getPlot(disc, p, ts, f, t, n, "time", "dt", levID)
			gpData[file].gpOptions = {logscale = {x = true, y = false, z = true}}
							
			-- single (grouping lev)
			local file = plotPath..table.concat({disc,p,ts,f,"__",t,n,"time","dt"},"_")
			gpData[file] = getPlot(disc, p, ts, f, t, n, "time", "dt")
			gpData[file].gpOptions = {logscale = {x = true, y = false, z = true}}
			
			-- grouping by (disc+p)								
			local file = plotPath..table.concat({disc,p,ts,f,"__",t,n,"time","dt","__",levID},"_")
			gpData[file] = getPlot(disc,    ts, f, t, n, "time", "dt", levID)
			gpData[file].gpOptions = {logscale = {x = true, y = false, z = true}}

			-- grouping by (disc+p)	(grouping lev)	 						
			local file = plotPath..table.concat({disc,p,ts,f,"__",t,n,"time","dt"},"_")
			gpData[file] = getPlot(disc,    ts, f, t, n, "time", "dt")
			gpData[file].gpOptions = {logscale = {x = true, y = false, z = true}}
	
			-- grouping (all discs+p)
			local file = plotPath..table.concat({disc,p,ts,f,"__",t,n,"time","dt","__",levID},"_")
			gpData[file] = getPlot("all",   ts, f, t, n, "time", "dt", levID)
			gpData[file].gpOptions = {logscale = {x = true, y = false, z = true}}					

			-- grouping (all discs+p) (grouping lev)	 
			local file = plotPath..table.concat({disc,p,ts,f,"__",t,n,"time","dt"},"_")
			gpData[file] = getPlot("all",   ts, f, t, n, "time", "dt")
			gpData[file].gpOptions = {logscale = {x = true, y = false, z = true}}					
		end
		end
						end
					end
				end
			end	
		end
	end
	print("done.")
	
	--------------------------------------------------------------------
	--  Execute Plot of gnuplot
	--------------------------------------------------------------------
	
	-- save for reuse
	persistence.store(dataPath.."gp-data-files.lua", gpData);	
	
	-- create scheduled plots
	if CRS.noplot == nil or CRS.noplot == false then
		write(">> Creating gnuplot plots ... ")
		for plotFile, data in pairs(gpData) do 
			if plotFile ~= "dirs" then
				local opt = table.deepcopy(gpOptions)
				if data.gpOptions then 
					local plotOpt = table.deepcopy(data.gpOptions)
					for k,v in pairs(plotOpt) do
						opt[k] = v
					end
				end
				--gnuplot.plot(plotFile..".tex", data, opt)
				gnuplot.plot(plotFile..".pdf", data, opt)
			end	
		end
		print("done.")
	end
	
end

function util.rates.kinetic.replot(gpOptions, file)
	if ProcRank() ~= 0 then return end
	
	local dataPath = "data/"
	local plotPath = "plots/"
	
	local function ensureDir(name)
		if not(DirectoryExists(name)) then CreateDirectory(name) end
	end
		
	local file = file or dataPath.."gp-data-files.lua"
	local gpData = persistence.load(file);
	
	ensureDir(dataPath)
	ensureDir(plotPath)
	for _, dir in pairs(gpData["dirs"]) do
		ensureDir(dir)
	end
	
	-- create scheduled plots
	for plotFile, data in pairs(gpData) do 
		if plotFile ~= "dirs" then
			local opt = table.deepcopy(gpOptions)
			if data.gpOptions then 
				local plotOpt = table.deepcopy(data.gpOptions)
				for k,v in pairs(plotOpt) do
					opt[k] = v
				end
			end
			gnuplot.plot(plotFile..".tex", data, opt)
			gnuplot.plot(plotFile..".pdf", data, opt)
		end
	end	
end