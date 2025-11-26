-- Copyright (c) 2011-2013:  G-CSC, Goethe University Frankfurt
-- Author: Sebastian Reiter
-- 
-- This file is part of UG4.
-- 
-- UG4 is free software: you can redistribute it and/or modify it under the
-- terms of the GNU Lesser General Public License version 3 (as published by the
-- Free Software Foundation) with the following additional attribution
-- requirements (according to LGPL/GPL v3 §7):
-- 
-- (1) The following notice must be displayed in the Appropriate Legal Notices
-- of covered and combined works: "Based on UG4 (www.ug4.org/license)".
-- 
-- (2) The following notice must be displayed at a prominent place in the
-- terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
-- 
-- (3) The following bibliography is recommended for citation and must be
-- preserved in all covered files:
-- "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
--   parallel geometric multigrid solver on hierarchically distributed grids.
--   Computing and visualization in science 16, 4 (2013), 151-164"
-- "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
--   flexible software system for simulating pde based models on high performance
--   computers. Computing and visualization in science 16, 4 (2013), 165-179"
-- 
-- This program is distributed in the hope that it will be useful,
-- but WITHOUT ANY WARRANTY; without even the implied warranty of
-- MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
-- GNU Lesser General Public License for more details.

--[[!
\addtogroup scripts_util
\{
\file scaling_analyzer.lua
\author Sebastian Reiter
\brief lua script to compare profiling outputs from different runs of one problem

This lua script can be used to compare profiling outputs from different
runs of one problem. Edit the inFiles-list below, so that it contains your
output-files. Make sure, that all files contain the same profiler output.
Run the script using the lua command-line-interpreter
(call lua scaling_analyzer.lua).
Optionally you can also run it with ug4 (ugshell -ex scaling_analyzer.lua)

The script outputs the timings of associated profiler node in the different
files in one row, together with the speedup factor between neighbored runs.

Note - if the timings do not exactly match and if some nodes in each run
have the same name, then this script may give you unexpected results.

If an input-file contains lines beginning with "#ANALYZER INFO:", the rest
of each such line will be printed during parsing of the file.
]]--

-- switch defining amount of output:
verbose = false

--	the files which will be processed
inFiles = {"log_refs_7.txt", "log_refs_8.txt",
		   "log_refs_9.txt", "log_refs_10.txt"}


--	num digits behind comma
numDigitsBehindComma = 2
minimumTime = 0 --- set this to 0.1 to display only times with > 0.1 s

print()
print("Scaling Analyzer")

--	given a time unit as a string (supported: "s", "ms", "us", "ns")
--	this table returns the factor, which scales the unit to seconds.
timeFactors = {}
timeFactors["s"] = 1
timeFactors["ms"] = 1 / 1000
timeFactors["us"] = 1 / 1000000
timeFactors["ns"] = 1 / 1000000000


--	Adds spaces until the string has size minSize. The result is returned
function FillSpaces(str, minSize)
	local strLen = string.len(str)
	if strLen < minSize then
		return str .. string.rep(" ", minSize - strLen) 
	end
	return str
end

--	Returns an entry in entryList, where entry.name == name.
--	Search starts at guessInd, and then continues in both directions
--	until both search indices point to an invalid entry.
--	If no entry was found, the method then starts a new straight forward
--	search over all entries in entryList
function FindEntry(entryList, name, guessInd)
	guessInd = math.floor(guessInd)
	local doSearch = true
	local offset = 0
	while doSearch == true do
		ind1 = guessInd - offset
		ind2 = guessInd + offset

		ind1Valid = false
		ind2Valid = false
				
		if ind1 >= 0 then
			local entry = entryList[ind1]
			if entry ~= nil then
				ind1Valid = true
				if entry.name == name then
					return entry
				end
			end
		end
		
		if ind2 >= 0 then
			local entry = entryList[ind2]
			if entry ~= nil then
				ind2Valid = true
				if entry.name == name then
					return entry
				end
			end
		end
		
		if ind1Valid == false and ind2Valid == false then
			break
		end
		
		offset = offset + 1
	end
	
--	if we reach this point we didn't find anything. Now perform
--	a straight search in all entries of entryList, in case the
--	list contained some nil entries
	for _, entry in ipairs(entryList) do
		if entry.name == name then
			return entry
		end
	end
	
	return nil
end



--	for each file we'll add a list to timings.
--	This list again contains tuples of {name, time, timeUnit, spaces}
--	for each profiling entry
timings = {}


-- parse all input files and extract profiling information.
-- This information will be written to lists which are then
-- added to the timings array.
fileCounter = 1
print("Files:")
for _, fileName in ipairs(inFiles) do
	print("("..fileCounter .. ")\t" .. fileName .. " ...")

	local f = io.open(fileName, "r")
	if f ~= nil then
		local printedInfo = false
		local readingProfilerOutput = false
		local fileTimings = {}
		
	--	iterate over all lines
		local counter = 1
		for line in f:lines() do
			if readingProfilerOutput == false then
				local i, j = string.find(line, "call tree")
				if i ~= nil then
					readingProfilerOutput = true
				else
					local str = string.match(line, "#ANALYZER INFO:%s*(.+)")
					if str ~= nil then
						if verbose == true then print("", "  - " .. str) end
						printedInfo = true
					end
				end
			else -- (i.e. 'readingProfilerOutput == true')
			--	get entries for each line through pattern matching
				local 	spaces, name, hits, selfTime, selfUnit, selfPerc,
						totalTime, totalUnit, totalPerc
						-- Note: support for large and small numbers containing e+ and e- is only achieved
						--       by adding those characters to the corresp. capture patterns.
						--       The actual conversion of strings to numbers is delegated to 'tonumber()'!
						--       (No '-' in the 'hit' capture pattern since we always expect that #hits is an integer ...)
						= string.match(line, "(%s*)(<*%a[_%.():%w]+>*)%s+([%.%deE+]+)%s+" -- leading spaces, name and hits
								   .."([%.%deE+-]+)%s+(%a+)%s+(%d+)%%%s+"	  -- self:  time, unit, percent
								   .."([%.%deE+-]+)%s+(%a+)%s+(%d+)%%")		  -- total: time, unit, percent
										   
				if name == nil then
				--	we reached the end of the profiler output
					readingProfilerOutput = false
				else
				--	add the timing to the fileTimings list
					local entry = {}
					entry.spaces = spaces
					entry.name = name
					entry.time = tonumber(totalTime)
					entry.timeUnit = totalUnit
					
					fileTimings[counter] = entry
					counter = counter + 1
				end
			end -- (of "else" (i.e. 'readingProfilerOutput == true')
		end -- iterate over all lines
		
	--	add an empty line if info has been printed
		if printedInfo == true then
			if verbose == true then print() end
		end
		
		f:close()
		
	--	add the timings of this file to the global timings
		timings[fileCounter] = fileTimings
		fileCounter = fileCounter + 1
	else
		print("  file " .. fileName .. " not found. Ignoring file.")
	end
end


--	If the profiling contain a "main" entry, we will subtract all timings
--	from direct child entries of "main" and store the result in an additional
--	entry named "unknown"
for _, timing in ipairs(timings) do
	local mainEntry = nil
	local childTimingSum = 0
	local highestIndex = 0
	local spaces = "  "
	for ind, entry in ipairs(timing) do
		if entry.name == "main" then
			mainEntry = entry
			spaces = entry.spaces.."  "
		elseif mainEntry ~=nil and entry.spaces == spaces then
			childTimingSum = childTimingSum + entry.time * timeFactors[entry.timeUnit]		
		end
		highestIndex = ind
	end
	
--	if a "main" entry was found, we will now add a new child entry called "unknown"
	if mainEntry ~= nil then
		local entry = {}
		entry.spaces = spaces
		entry.name = "unknown"
		entry.time = mainEntry.time * timeFactors[mainEntry.timeUnit] - childTimingSum
		entry.timeUnit = "s"
		
	--	make the time 'nice'
		
		
		timing[highestIndex + 1] = entry
	end
end

--------------------------------------------------------------------------------

--	reading input finished - now print all the timings
--	file 1 is the main file - for each entry, we'll check the other timings
--	for comparable entries and calculate the speedup from the timings
if timings[1] == nil or timings[2] == nil then
	print("at least 2 input files are required to compare timings.")
else
	local mainTimings = timings[1]
	
	local title = FillSpaces("NAME", 48)
	for i, entry in ipairs(timings) do
		title = title .. FillSpaces("time "..i, 12)
		if timings[i+1] ~= nil then
			title = title .. FillSpaces("(fac)", 10)
		end
	end
	print()
	print(title)
	
--	iterate over all entries of mainTimings
	for mainEntryInd, mainEntry in ipairs(mainTimings) do
		local bPrint=false
		if minimumTime == 0 then
			bPrint=true
		end
	--	print the line
		local line = FillSpaces(mainEntry.spaces .. mainEntry.name, 48)
		
	--	iterate over all timings and calculate speedups
		local timingInd1 = 1
		local timingInd2 = 2
		
		while timings[timingInd2] ~= nil do
			local timings1 = timings[timingInd1]
			local timings2 = timings[timingInd2]

		-- find entries with the same names as mainEntry in timings1 and timings2 and
		-- store the times in seconds.
			local time1 = nil
			local curEntry1 = FindEntry(timings1, mainEntry.name, mainEntryInd)
			if curEntry1 ~= nil then
				time1 = curEntry1.time * timeFactors[curEntry1.timeUnit]
			end
			
			local time2 = nil
			local curEntry2 = FindEntry(timings2, mainEntry.name, mainEntryInd)
			if curEntry2 ~= nil then
				time2 = curEntry2.time * timeFactors[curEntry2.timeUnit]
			end
			
			local tStr = ""
			if curEntry1 == nil then
				tStr = "----"
			else
				tStr = curEntry1.time .. " " .. curEntry1.timeUnit
			end
			line = line .. FillSpaces(tStr, 12)
			
		--	calculate the speedup factor
			if time1 == nil or time2 == nil then
				tStr = "(???)"
			else			
				if time1 == 0  then
					tStr = "(nan)"
				else
					tStr = "(" ..
						   string.format("%."..numDigitsBehindComma.."f", time2 / time1) ..
						   ")"
				end
			end
			
			if time1 ~= nil and time2 ~= nil and (time2 > minimumTime or time1 > minimumTime) then
			   bPrint = true
			end
			
			line = line .. FillSpaces(tStr, 10)
		
		--	if the next timing is the last, we'll add its times here
			if timings[timingInd2 + 1] == nil then
				if curEntry2 == nil then
					tStr = "----"
				else
					tStr = curEntry2.time .. " " .. curEntry2.timeUnit
				end
				line = line .. FillSpaces(tStr, 12)
			end
				
			
		--	compare the next timings
			timingInd1 = timingInd1 + 1
			timingInd2 = timingInd2 + 1
		end -- "while timings[timingInd2] ~= nil do"
		
		if bPrint then
				print(line)
		end
	end
end

--[[!
\}
]]--
