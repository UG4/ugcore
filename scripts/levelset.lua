----------------------------------------------------------
--
--   Lua - Script to test the Level Set mechanism
--
--   Author: Christian Wehner
--
----------------------------------------------------------

ug_load_script("ug_util.lua")

-- choose algebra
InitAlgebra(CPUAlgebraSelector());

-- constants
dim = 2

if  dim == 2 then 
gridName = util.GetParam("-grid", "unit_square_01/unit_square_01_tri_2x2.ugx")
--gridName = util.GetParam("-grid", "unit_square_01/unit_square_01_quads_2x2.ugx")
elseif dim == 3 then gridName = util.GetParam("-grid", "unit_square/unit_cube_hex.ugx")
end

numPreRefs = util.GetParamNumber("-numPreRefs", 0)
numRefs    = util.GetParamNumber("-numRefs",    3)

print(" Choosen Parater:")
print("    numRefs    = " .. numRefs)
print("    numPreRefs = " .. numPreRefs)
print("    grid       = " .. gridName)

--------------------------------------------------------------------------------
-- Domain Setup
--------------------------------------------------------------------------------

-- Create, Load, Refine and Distribute Domain
neededSubsets = {"Inner", "Boundary"}
dom = util.CreateAndDistributeDomain(gridName, numRefs, numPreRefs, neededSubsets)

-- create Approximation Space
print("Create ApproximationSpace")
approxSpace = util.CreateApproximationSpace(dom)
approxSpace:add_fct("c", "Lagrange", 1)
approxSpace:init()
approxSpace:print_statistic()

--------------------------------------------------------------------------------
-- User Data Setup
--------------------------------------------------------------------------------
	
-- Diffusion Tensor setup
diffusionMatrix = util.CreateConstDiagUserMatrix(1.0, dim)

function ourRhs2d(x, y, t)
	local s = 2*math.pi
	return	s*s*(math.sin(s*x) + math.sin(s*y))
end

-- rhs setup
rhs = util.CreateLuaUserNumber("ourRhs"..dim.."d", dim)

function ExactSolution(x, y, t)
	local s = 2*math.pi
        local s = 10;
        return math.sin(x)+4*math.cos(y);
	-- return math.sin(s*x) + math.sin(s*y)
end

-- The dirichlet condition
function DirichletBnd2d(x, y, t)
	return true, ExactSolution(x, y, t)
end

-- dirichlet setup
dirichlet = util.CreateLuaBoundaryNumber("DirichletBnd"..dim.."d", dim)
	
--------------------------------------------------------------------------------
--  Setup Dirichlet Boundary
--------------------------------------------------------------------------------

dirichletBND = util.CreateDirichletBoundary(approxSpace)
dirichletBND:add(dirichlet, "c", "Boundary")

--------------------------------------------------------------------------------
--  Create a grid function
--------------------------------------------------------------------------------

-- get grid function
u = approxSpace:create_surface_function()

-- set initial value
print("Interpolation start values")-- start value
startValue = util.CreateLuaUserNumber("ExactSolution", dim)
time = 0.0
InterpolateFunction(startValue, u, "c", time)

levDisc = FV1LevelSetDisc2d()
levDisc:set_dt(0.1);
levDisc:add_post_process(dirichletBND);
levDisc:compute_error(u,0);
--------------------------------------------------------------------------------
--  Output
--------------------------------------------------------------------------------

WriteGridFunctionToVTK(u, "Solution")
SaveVectorForConnectionViewer(u, "u.mat")
