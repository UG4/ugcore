-- created by Martin Rupp, 06.02.2012, martin.rupp@gcsc.uni-frankfurt.de
util = util or {}

-- This file and the comments are intended to make your life easier when you want to do several runs of a script and
-- you need to evaluate data from the runs, like compare solution times for different number of refinements.

-- First you need to know, what statistics you want to output and compare. 
-- Imagine you want to run laplace.lua for several different numRefs or even num procs, 
-- and results should be saved in a file "stats.txt", so that you can use in a spreadsheet.
-- The stats format is as follows: it is an array, consisting of entries {description, content}
-- So we need to add to our laplace.lua file the following lines:
--     stats = { {"procs", GetNumProcesses()}, {"numRefs", numRefs}, {"steps", convCheck:step()}, {"SVN Revision", GetSVNRevision()} }
--     util.writeFileStats(stats, "stats.txt")
-- after 2 runs, the file stats.txt looks as follows
-- procs 	numRefs 	steps 	SVN Revision
-- 4 	8 	6 	4459M
-- 16 	9 	6 	4459M
-- note that a heading "procs	numRefs..." is added to the file if it is empty
-- note that all items are seperated by default by " \t", but you can choose different seperators
-- you can copy and paste the data in most spreadsheet applications like OpenOffice Calc, Excel and Apple Numbers. 
-- if you wish a brief table result description of the current task, you can use
-- util.printStats(stats)
-- which prints
--   procs: 4
--   numRefs: 16
--   steps: 6
--   SVN Revision: 4459M
-- to the console
-- 
-- for boolean values, use tostring(bExternalCoarsening) -> "true" / "false"
--
--
-- you can seperate your runs by date with the following example bash script:
-- mydate=`date +%Y-%m-%d-%H.%M.%S`
-- mkdir $mydate
-- for i in {3..6}
-- do
--   ./ugshell -ex mycheck.lua -outdir $mydate/ -logtofile $mydate/mycheck${i} -numRefs $i
-- done
-- 
-- and in you lua-file, you could use util.writeFileStats(stats, util.GetParam("-outdir", "").."stats.txt")

-- instead of defining -logtofile from the shell, you could also make something like
-- GetLogAssistant():enable_file_output(true, util.GetParam("-outdir". "").."check_numRefs"..numRefs.."_procs"..procs..".txt")

--- util.getStats
-- internal method
function util.getStats(stats, bHeader, seperator, bStats, seperator2)
	local output=""
	for i,v in ipairs(stats) do
		if bHeader then output = output..v[1]..seperator end
		if bStats then
			if v[2] == nil then
				print("value for item "..v[1].." is nil!")
				output = output.." "..seperator2
			else				
				output = output..tostring(v[2])..seperator2
			end
		end
	end
	output = output.."\n"
	return output
end

--- util.writeFileStats
-- writes stats to a file
-- @param stats
-- @param filename
-- @param seperator (default " \t")
function util.writeFileStats(stats, filename, seperator)
	if seperator == nil then
		seperator = " \t"
	end
	local output = io.open(filename, "a")
	
	if fsize(output) == 0 then
		print("file is empty, writing header...")
		output:write(util.getStats(stats, true, seperator, false, seperator))		
	end
	output:write(util.getStats(stats, false, seperator, true, seperator))
end

--- util.printStats
-- prints stats to the console
-- @param stats
-- @param filename
function util.printStats(stats)
	print(util.getStats(stats, true, ": ", true, "\n"))
end
