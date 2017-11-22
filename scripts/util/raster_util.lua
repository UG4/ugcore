-- Copyright (c) 2017:  G-CSC, Goethe University Frankfurt
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

util = util or {}
util.raster = util.raster or {}


-- Returns a callback function and a raster object.
-- Arguments of the callback depend on the specified dimension:
-- 1d: 'number callback (x)'
-- 2d: 'number callback (x, y)'
-- 3d: 'number callback (x, y, z)'

function util.raster.CreateRasterValueCallback (filename, dim, interpOrder)
--	get the dimension of the underlying raster file
	local fileDim = util.raster.GetRasterFileDimension (filename)

	if fileDim == nil then
		print("ERROR in util.raster.CreateRasterCallback: invalid dimension: " .. fileDim)
		return nil
	end

	local fullFilename = FindFileInStandardPaths(filename)
	if fullFilename == "" then
		print("ERROR in util.raster.GetRasterFileDimension: file not found: '" .. filename .. "'")
	end

	local func = nil
	local raster = nil

	if fileDim == 1 then
		raster = NumberRaster1d()
		raster:load_from_asc(fullFilename)
		if dim == 1 then
			func =	function(x, t, si)
						raster:set_cursor(0, x);
						return raster:interpolate_at_cursor(interpOrder)
					end
		elseif dim == 2 then
			func =	function(x, y, t, si)
						raster:set_cursor(0, x);
						return raster:interpolate_at_cursor(interpOrder)
					end
		elseif dim == 3 then
			func =	function(x, y, z, t, si)
						raster:set_cursor(0, x);
						return raster:interpolate_at_cursor(interpOrder)
					end
		end

	elseif fileDim == 2 then
		raster = NumberRaster2d()
		raster:load_from_asc(fullFilename)
		if dim == 1 then
			func =	function(x, t, si)
						raster:set_cursor(0, x);
						return raster:interpolate_at_cursor(interpOrder)
					end
		elseif dim == 2 then
			func =	function(x, y, t, si)
						raster:set_cursor(0, x);
						raster:set_cursor(1, y);
						return raster:interpolate_at_cursor(interpOrder)
					end
		elseif dim == 3 then
			func =	function(x, y, z, t, si)
						raster:set_cursor(0, x);
						raster:set_cursor(1, y);
						return raster:interpolate_at_cursor(interpOrder)
					end
		end

	elseif fileDim == 3 then
		raster = NumberRaster3d()
		raster:load_from_asc(fullFilename)
		if dim == 1 then
			func =	function(x, t, si)
						raster:set_cursor(0, x);
						return raster:interpolate_at_cursor(interpOrder)
					end
		elseif dim == 2 then
			func =	function(x, y, t, si)
						raster:set_cursor(0, x);
						raster:set_cursor(1, y);
						return raster:interpolate_at_cursor(interpOrder)
					end
		elseif dim == 3 then
			func =	function(x, y, z, t, si)
						raster:set_cursor(0, x);
						raster:set_cursor(1, y);
						raster:set_cursor(2, z);
						return raster:interpolate_at_cursor(interpOrder)
					end
		end
	end

	return func, raster
end

-- Returns a callback function and a raster object.
-- Arguments of the callback depend on the specified dimension:
-- 1d: 'mat1x1 callback (x)'
-- 2d: 'mat2x2 callback (x, y)'
-- 3d: 'mat3x3 callback (x, y, z)'
--
-- where matNxN stands for a sequence (not a table) of N*N number values.

function util.raster.CreateRasterMatrixCallback (filename, dim, interpOrder)
--	get the dimension of the underlying raster file
	local fileDim = util.raster.GetRasterFileDimension (filename)

	if fileDim == nil then
		print("ERROR in util.raster.CreateRasterCallback: invalid dimension: " .. fileDim)
		return nil
	end

	local fullFilename = FindFileInStandardPaths(filename)
	if fullFilename == "" then
		print("ERROR in util.raster.GetRasterFileDimension: file not found: '" .. filename .. "'")
	end

	local func = nil
	local raster = nil

	if fileDim == 1 then
		raster = NumberRaster1d()
		raster:load_from_asc(fullFilename)
		if dim == 1 then
			func =	function(x, t, si)
						raster:set_cursor(0, x);
						return raster:interpolate_at_cursor(interpOrder)
					end
		elseif dim == 2 then
			func =	function(x, y, t, si)
						raster:set_cursor(0, x);
						local d = raster:interpolate_at_cursor(interpOrder)
						return d, 0, 0, d
					end
		elseif dim == 3 then
			func =	function(x, y, z, t, si)
						raster:set_cursor(0, x);
						local d = raster:interpolate_at_cursor(interpOrder)
						return d, 0, 0, 0, d, 0, 0, 0, d
					end
		end

	elseif fileDim == 2 then
		raster = NumberRaster2d()
		raster:load_from_asc(fullFilename)
		if dim == 1 then
			func =	function(x, t, si)
						raster:set_cursor(0, x);
						return raster:interpolate_at_cursor(interpOrder)
					end
		elseif dim == 2 then
			func =	function(x, y, t, si)
						raster:set_cursor(0, x);
						raster:set_cursor(1, y);
						local d = raster:interpolate_at_cursor(interpOrder)
						return d, 0, 0, d
					end
		elseif dim == 3 then
			func =	function(x, y, z, t, si)
						raster:set_cursor(0, x);
						raster:set_cursor(1, y);
						local d = raster:interpolate_at_cursor(interpOrder)
						return d, 0, 0, 0, d, 0, 0, 0, d
					end
		end

	elseif fileDim == 3 then
		raster = NumberRaster3d()
		raster:load_from_asc(fullFilename)
		if dim == 1 then
			func =	function(x, t, si)
						raster:set_cursor(0, x);
						return raster:interpolate_at_cursor(interpOrder)
					end
		elseif dim == 2 then
			func =	function(x, y, t, si)
						raster:set_cursor(0, x);
						raster:set_cursor(1, y);
						local d = raster:interpolate_at_cursor(interpOrder)
						return d, 0, 0, d
					end
		elseif dim == 3 then
			func =	function(x, y, z, t, si)
						raster:set_cursor(0, x);
						raster:set_cursor(1, y);
						raster:set_cursor(2, z);
						local d = raster:interpolate_at_cursor(interpOrder)
						return d, 0, 0, 0, d, 0, 0, 0, d
					end
		end
	end

	return func, raster
end

function util.raster.GetRasterFileDimension (filename)
	local fullFilename = FindFileInStandardPaths(filename)
	if fullFilename == "" then
		print("ERROR in util.raster.GetRasterFileDimension: file not found: '" .. filename .. "'")
	end

--todo: Open the file on one process and scatter the result
	local f = io.open(fullFilename, "r")
	if f == nil then
		print("ERROR in util.raster.GetRasterFileDimension: file not found: '" .. filename .. "'")
		return nil
	end

	local dim = 0
	while true do
		local l = f:read("*line")
		if string.find(l, "ncols") and dim < 1 then dim = 1
		elseif string.find(l, "nrows") and dim < 2 then dim = 2
		elseif string.find(l, "nstacks") and dim < 3 then dim = 3
		else break
		end
	end

	return dim
end
