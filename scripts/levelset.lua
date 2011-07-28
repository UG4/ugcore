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
gridName = util.GetParam("-grid", "unit_square_01/unit_square_01_tri_2x2.ugx")
gridName = util.GetParam("-grid", "unit_square/unit_square_unstructured_tris_coarse.ugx")
-- gridName = util.GetParam("-grid", "unit_square_01/unit_square_01_quads_2x2.ugx")
elseif dim == 3 then gridName = util.GetParam("-grid", "unit_square/unit_cube_hex.ugx")
end

numPreRefs = util.GetParamNumber("-numPreRefs", 0)
numRefs    = util.GetParamNumber("-numRefs",    6)

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
        return 2*x-y;
--      return math.sqrt(x*x + y*y) - 0.5;
--        return x-t;
        ---
        -- return x - 4*y;
        -- return math.sin(x)+4*math.cos(y);
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
phiNew = approxSpace:create_surface_function();
phiOld = approxSpace:create_surface_function();
phiNew:set(0);

-- set initial value
print("Interpolation start values")-- start value
startValue = util.CreateLuaUserNumber("ExactSolution", dim)
time = 0.0
InterpolateFunction(startValue, phiOld, "c", time)

lsDisc = FV1LevelSetDisc2d()
lsDisc:set_dt(0.2/2/2/2/1.5/2);
lsDisc:add_post_process(dirichletBND);
-- lsDisc:set_info(false);
------- perform all steps in one call
------- lsDisc:set_nr_of_steps(40);
------- lsDisc:advect_lsf(phiNew,phiOld);
-- choose number of time steps
NumTimeSteps =  util.GetParamNumber("-numTimeSteps", 40)
-- lsDisc:compute_error(phiNew);

for step = 1, NumTimeSteps do
    lsDisc:advect_lsf(phiNew,phiOld);
    VecScaleAssign(phiOld,1,phiNew);
end
lsDisc:compute_error(phiNew);
--------------------------------------------------------------------------------
--  Output
--------------------------------------------------------------------------------

WriteGridFunctionToVTK(phiNew, "Solution")
SaveVectorForConnectionViewer(phiNew, "u.mat")



