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

function util.rates.kinetic.StdLabel.SpaceLatex(x)
	local gpSpaceLabel ={ DoFs = "Anzahl Unbekannte",	h = "h (Gitterweite)"}
	return gpSpaceLabel[x]
end

function util.rates.kinetic.StdLabel.TimestepLatex()
	return "$\\Delta t$ (Zeitschrittweite)"
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
	local gpSpaceLabel ={ DoFs = "Anzahl Unbekannte",	h = "h (Gitterweite)"}
	return gpSpaceLabel[x]
end

function util.rates.kinetic.StdLabel.TimestepPdf()
	return "dt (Zeitschrittweite)"
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
	local MaxLevelPadding = 	CRS.MaxLevelPadding 	or util.rates.kinetic.StdMaxLevelPadding
	
	if 	CreateApproxSpace == nil or CreateDomainDisc == nil or 
		CreateSolver == nil or CreateDomain == nil then
		print("You must pass: CreateApproxSpace, CreateDomainDisc, CreateSolver, CreateDomain")
		exit()
	end
	
	local SpaceDiscs = CRS.SpaceDiscs
	local TimeDiscs = CRS.TimeDiscs

	local maxlevel = CRS.maxlevel; 		if maxlevel == nil then maxlevel = true end
	local exact = CRS.exact; 			if exact == nil then exact = true end
	local plotSol = CRS.plotSol; 		if plotSol == nil then plotSol = false end
	local onlyLast = CRS.onlyLast; 		if onlyLast == nil then onlyLast = true end

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
			
			--------------------------------------------------------------------
			--  Loop Time discs
			--------------------------------------------------------------------
			for _, TimeDisc in ipairs(TimeDiscs) do

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
						Interpolate(ExactSol["c"], mem.u, "c", mem.time)				
						mem.TimeSeries:push(mem.u:clone(), mem.time)
						
						if timeDisc:num_stages() > 1 then mem.uOld = mem.u:clone() end			
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
		
							-- set order for bdf to 1 (initially)
							if ts:lower() == "bdf" then timeDisc:set_order(1) end			
			
							while mem.time < sliceTime do
					
								-- update step count
								mem.step = mem.step + 1
																	
								-- time step size
								local dt = mem.dt
								local dodt = dt
														
								if mem.time + dt > sliceTime then dodt = sliceTime - mem.time end
								if ((sliceTime - (mem.time+dt))/dt) < 1e-8 then dodt = sliceTime - mem.time end
						
								-- get old solution if multistage
								if timeDisc:num_stages() > 1 then
									VecScaleAssign(mem.u, 1.0, mem.uOld)
								end			
			
								for stage = 1, timeDisc:num_stages() do
									timeDisc:set_stage(stage)
									timeDisc:prepare_step(mem.TimeSeries, dodt)
			
									-- solve step						
									newtonSolver:init(AssembledOperator(timeDisc, mem.u:grid_level()))
									if newtonSolver:prepare(mem.u) == false then print (">> Newton init failed."); exit(); end 
									if newtonSolver:apply(mem.u) == false then print (">> Newton solver failed."); exit(); end 
													
									-- update new time
									mem.time = timeDisc:future_time()
									if math.abs(sliceTime - mem.time) < 1e-8*dt then mem.time = sliceTime end
									
									-- push oldest solutions with new values to front, oldest sol pointer is poped from end	
									if ts:lower() == "bdf" and mem.step < orderOrTheta then
										print("++++++ BDF: Increasing order to "..mem.step+1)
										timeDisc:set_order(mem.step+1)
										solTimeSeries:push(mem.u:clone(), mem.time)
									else 
										local oldestSol = mem.TimeSeries:oldest()
										VecScaleAssign(oldestSol, 1.0, mem.u)
										mem.TimeSeries:push_discard_oldest(oldestSol, mem.time)
									end
								end
					
								-- save this solution if multistage
								if timeDisc:num_stages() > 1 then
									mem.uOld = mem.u
								end			
							end	-- end slice interval
							
							if plotSol then
								vtk = VTKOutput()
								vtk:print("Sol"..mem.step, mem.u)
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
								if maxlevel and lev < maxLev and k < refs then
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

	write(">> Writing measured rates to data files ...")
	for disc, _ in pairs(errors) do
		for p, _ in pairs(errors[disc]) do
			for ts, _ in pairs(errors[disc][p]) do
				for f, _ in pairs(PlotCmps) do
					for t, _ in pairs(errors[disc][p][ts][f]) do
						for n, _ in pairs(errors[disc][p][ts][f][t]) do

		-- values for each time step
		local err = errors[disc][p][ts]
		for tp, _ in pairs(errors[disc][p][ts][f][t][n]) do
	
			if onlyLast then tp = #errors[disc][p][ts][f][t][n] end
			
			local dir = dataPath..tp.."/"
			ensureDir(dir)				
							
			local value = errors[disc][p][ts][f][t][n][tp].value
	
			-- convergence in time
			for lev, _ in iipairs(value.dt) do		
				local file = dir..table.concat({"error",disc,p,ts,f,t,n,"lev"..lev},"_")..".dat"
				local cols = {err.DoFs, err.h, value.dt[lev]}
				gnuplot.write_data(file, cols)
			end
					
			-- convergence in space
			for k, _ in iipairs(value.h) do						
				local file = dir..table.concat({"error",disc,p,ts,f,t,n,"dt"..k.."["..err.dt[k].."]"},"_")..".dat"
				local cols = {err.DoFs, err.h, value.h[k]}
				gnuplot.write_data(file, cols)
			end
	
			-- convergence in space-time
			local file = dir..table.concat({"error",disc,p,ts,f,t,n},"_")..".dat"
			local fio = io.open(file, "w+")
			fio:close()
			for k, _ in pairs(value.h) do						
				local dt = err.dt[k]
				local dts = {}; for i,_ in pairs(value.h[k]) do dts[i] = dt end
				local cols = {err.DoFs, err.h, dts, value.h[k]}
				gnuplot.write_data(file, cols, false, "a+")
				local fio = io.open(file, "a+")
				fio:write("\n")
				fio:close()
			end
			
			for xCol, x in ipairs({"DoFs", "h"}) do
				local dataset = {label=MeasLabel(disc, p), file=file, style="linespoints", xCol, 3, 4}
				local label = { x = SpaceLabel(x), y = TimestepLabel(), z = NormLabel(f,t,n)}
								
				addSet( accessPlot(disc, p, ts, tp, f, t, n, x), dataset, label)
				addSet( accessPlot(disc,    ts, tp, f, t, n, x), dataset, label)
				addSet( accessPlot("all",   ts, tp, f, t, n, x), dataset, label)
			end
			
			if onlyLast then break end
		end
		
		-- values for the time series
		local value = errors[disc][p][ts][f][t][n][1].value
		for lev, _ in ipairs(value) do		
			for k, _ in pairs(value[lev]) do		
				local name = table.concat({disc,p,ts,f,t,n,"lev"..lev,"dt"..k.."["..err.dt[k].."]"},"_")
				local file = dataPath.."error"..name..".dat"

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
	
				local dataset = {label = MeasLabel(disc, p), file=file, style="linespoints", 1, 2}
				local label = { x = TimeLabel(), y = NormLabel(f,t,n)}
				local gpFile = plotPath..name
				local plot = {}
				table.insert( plot, dataset)			
				plot.label = label
				
				gpData[gpFile] = plot
				gpData[gpFile].gpOptions = {logscale = false}
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
	for disc, _ in pairs(errors) do
		for p, _ in pairs(errors[disc]) do
			for ts, _ in pairs(errors[disc][p]) do
				for f, _ in pairs(PlotCmps) do
					for t, _ in pairs(errors[disc][p][ts][f]) do
						for n, _ in pairs(errors[disc][p][ts][f][t]) do
							for _, x in ipairs({"DoFs", "h"}) do

		for tp, _ in pairs(errors[disc][p][ts][f][t][n]) do
			if onlyLast then tp = #errors[disc][p][ts][f][t][n] end
		
			local dir = plotPath..tp.."/"
			ensureDir(dir)				
						
			-- single dataset			
			local file = dir..table.concat({disc,p,ts,f,t,n,x},"_")
			gpData[file] = getPlot(disc, p, ts, tp, f, t, n, x)
			gpData[file].gpOptions = {logscale = true}
	
			-- grouping by (disc+p)								
			local file = dir..table.concat({f,disc,ts,t,n,x}, "_")	
			--gpData[file] = getPlot(disc, ts, tp, f, t, n, x)		
	
			-- grouping (all discs+p)
			local file = dir..table.concat({f,"all",ts,t,n,x}, "_")	
			--gpData[file] = getPlot("all", ts, tp, f, t, n, x)	

			if onlyLast then break end
		end
							end
						end
					end
				end
			end
		end
	end
	
	-- create scheduled plots
	if CRS.noplot == nil or CRS.noplot == false then
		write(">> Creating gnuplot plots ... ")
		for plotFile, data in pairs(gpData) do 
			local opt = table.deepcopy(gpOptions)
			if data.gpOptions then 
				local plotOpt = table.deepcopy(data.gpOptions)
				for k,v in pairs(plotOpt) do
					opt[k] = v
				end
			end
			--gnuplot.plot(plotFile..".tex", data, opt)
			--print(">> plotting: "..plotFile..".pdf")
			gnuplot.plot(plotFile..".pdf", data, opt)
		end	
		print("done.")
	end
	
end
