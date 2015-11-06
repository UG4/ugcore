-- Copyright (c) 2014:  G-CSC, Goethe University Frankfurt
-- Author: Martin Rupp
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

--------------------------------------------------------------------------------
-- User Data utils
--------------------------------------------------------------------------------

--! creates a Const User Matrix 2d 
function util.CreateConstUserMatrix2d(m00, m01, m10, m11)
	local mat = ConstUserMatrix2d()
	mat:set_entry(0, 0, m00)
	mat:set_entry(0, 1, m01)
	mat:set_entry(1, 0, m10)
	mat:set_entry(1, 1, m11)	
	return mat
end

--! creates a Const User Matrix 3d
function util.CreateConstUserMatrix3d(m00, m01, m02, m10, m11, m12, m20, m21, m22)
	local mat = ConstUserMatrix3d()
	mat:set_entry(0, 0, m00)
	mat:set_entry(0, 1, m01)
	mat:set_entry(0, 2, m02)
	mat:set_entry(1, 0, m10)
	mat:set_entry(1, 1, m11)
	mat:set_entry(1, 2, m12)
	mat:set_entry(2, 0, m20)
	mat:set_entry(2, 1, m21)
	mat:set_entry(2, 2, m22)	
	return mat
end

--! creates a Const User Vector 2d
function util.CreateConstUserVector2d(v0, v1)
	local vec = ConstUserVector2d()
	vec:set_entry(0, v0)
	vec:set_entry(1, v1)
	return vec
end

--! creates a Const User Vector 3d
function util.CreateConstUserVector3d(v0, v1, v2)
	local vec = ConstUserVector3d()
	vec:set_entry(0, v0)
	vec:set_entry(1, v1)
	vec:set_entry(2, v2)
	return vec
end


--------------------------------------------------------------------------------
-- some auxiliary functions
--------------------------------------------------------------------------------
--! function returns true if the number is a power of two
function util.IsPowerOfTwo(n)
	local number compare = 1

	while (compare < n) do
		compare = compare*2
	end

	return compare==n
end

--! function returns true if the number is a natural number
function util.IsNaturalNumber(n)
	if n-math.floor(n) == 0 then
		return true
	else
		return false
	end
end

--! function to factorise number which has to be a power of 2 in two factors
--! which differ at most by a factor of 2 and returns both
--! (first the smaller one, then the larger one).
function util.FactorizeInPowersOfTwo(n)
	if not util.IsPowerOfTwo(n) then
		print("Number to factorise must be a power of 2. Aborting.")
		exit()
	end

	local number firstFactor = n
	local number secFactor = 1

	while (firstFactor > 2*secFactor) do
		firstFactor = firstFactor/2
		secFactor = secFactor*2
	end

	return secFactor, firstFactor
end