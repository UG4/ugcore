----------------------------------------------------------
--
--   Lua - Script to test the Level Set mechanism
--
--   Author: Christian Wehner
--
----------------------------------------------------------

ug_load_script("ug_util.lua")

-- constants
dim = 2

movie=false;

neumann = false;

-- setup 0 : use hard-coded functions for velocity/boundary conditions
-- setup 1 : use lua-functions for velocity/boundary conditions
-- setup 2 : use vector field for velocity
-- setup 3 : use constant data
setup = 0;

if  dim == 2 then 
gridName = util.GetParam("-grid", "unit_square_01/unit_square_01_tri_2x2.ugx")
gridName = util.GetParam("-grid", "unit_square/unit_square_quads_2x2.ugx")
if neumann==true then
    gridName = util.GetParam("-grid", "unit_square/dirichlet_neumann_2x2.ugx")
end;
--gridName = util.GetParam("-grid", "unit_square/unit_square_unstructured_tris_coarse.ugx")
--gridName = util.GetParam("-grid", "unit_square/unit_square_unstructured_tris_coarse_neu_dir.ugx")
--gridName = util.GetParam("-grid", "neumann.ugx")
--gridName = util.GetParam("-grid", "unit_square_01/unit_square_01_quads_2x2.ugx")
elseif dim == 3 then gridName = util.GetParam("-grid", "unit_square/unit_cube_hex.ugx")
end

numPreRefs = util.GetParamNumber("-numPreRefs", 0)
numRefs    = util.GetParamNumber("-numRefs",    5)

print(" Choosen Parater:")
print("    numRefs    = " .. numRefs)
print("    numPreRefs = " .. numPreRefs)
print("    grid       = " .. gridName)

InitUG(dim, CPUAlgebraSelector());

--------------------------------------------------------------------------------
-- Domain Setup
--------------------------------------------------------------------------------

-- Create, Load, Refine and Distribute Domain
neededSubsets = {"Inner", "Boundary"}
dom = util.CreateAndDistributeDomain(gridName, numRefs, numPreRefs, neededSubsets)

-- create Approximation Space
print("Create ApproximationSpace")
approxSpace = ApproximationSpace(dom)
approxSpace:add_fct("c", "Lagrange", 1)
approxSpace:init()
approxSpace:print_statistic()


--------------------------------------------------------------------------------
-- User Data Setup
--------------------------------------------------------------------------------

function ourRhs2d(x, y, t)
	local s = 2*math.pi
	return	s*s*(math.sin(s*x) + math.sin(s*y))
end

-- rhs setup
rhs = util.CreateLuaUserNumber("ourRhs"..dim.."d", dim)

function ExactSolution(x, y, t)
--         delta = 0.0;
         xnew = x*math.cos(t)+y*math.sin(t);
	 ynew = -x*math.sin(t)+y*math.cos(t);
--	 xnew = x;
--	 ynew = y;
--	 return math.sqrt( (xnew-0.5)*(xnew-0.5)+ynew*ynew )-(0.3+delta*t);
	  return math.sqrt( (xnew-0.5)*(xnew-0.5)+ynew*ynew )-0.3;
--	local s = 2*math.pi
  --      local s = 10;
        --return 0;
       -- return x*x;
--       return x-t;
        --return math.min(math.min(x+1,1-x),math.min(y+1,1-y));
        --return 2*x-y;
        --return 2*x-y;
--      return math.sqrt(x*x + y*y) - 0.5;
--        return x-t;
        ---
        -- return x - 4*y;
        -- return math.sin(x)+4*math.cos(y);
	-- return math.sin(s*x) + math.sin(s*y)
end

function ConstantVSolution(x,y,t)
 --   delta = 0.0;
    xnew = x - t;
	ynew = y - t;
	return 2*xnew+ynew;
--	return math.sqrt( (xnew-0.5)*(xnew-0.5)+ynew*ynew )-(0.3+delta*t);
end

function vx(x,y,t)
--	return 0;
    return -y;
end

function vy(x,y,t)
--    return 0;
	return x;
end

-- The dirichlet condition
function DirichletBnd2d(x, y, t)
	return true, ExactSolution(x, y, t)
end

-- dirichlet setup
-- dirichlet = util.CreateLuaBoundaryNumber("Boundary"..dim.."d", dim)
	
--------------------------------------------------------------------------------
--  Setup Dirichlet Boundary
--------------------------------------------------------------------------------

--dirichletBND = util.CreateDirichletBoundary(approxSpace)
--dirichletBND:add(dirichlet, "c", "DirichletBnd")

--------------------------------------------------------------------------------
--  Create a grid function
--------------------------------------------------------------------------------

-- get grid function
phiNew = approxSpace:create_surface_function();
phiOld = approxSpace:create_surface_function();
phiNew:set(0);

lsDisc = FV1LevelSetDisc();
time = lsDisc:get_time();

-- add Dirichlet boundary post process
-- lsDisc:add_post_process(dirichletBND);

lsDisc:set_source(0);
lsDisc:set_delta(0.0);

-- lsDisc:set_neumann_boundary("Boundary");
--lsDisc:set_neumann_boundary(phiOld,"Boundary");
lsDisc:set_dirichlet_boundary(phiOld,"Boundary");

if (setup==0) then
	lsDisc:init_function(phiOld);
	lsDisc:set_vel_x();
	lsDisc:set_vel_y();
	lsDisc:set_dirichlet_data();
end;
if (setup==1) then
    solfunctor = util.CreateLuaUserNumber("ExactSolution", dim);
	lsDisc:set_dirichlet_data(solfunctor);
	InterpolateFunction(solfunctor, phiOld, "c", time);
	vxfunctor = util.CreateLuaUserNumber("vx", dim);
	vyfunctor = util.CreateLuaUserNumber("vy", dim);
	lsDisc:set_vel_x(vxfunctor);
	lsDisc:set_vel_y(vyfunctor);
end;
if (setup==2) then
    solfunctor = util.CreateLuaUserNumber("ExactSolution", dim);
	lsDisc:set_dirichlet_data(solfunctor);
--	lsDisc:set_dirichlet_data();
	InterpolateFunction(solfunctor, phiOld, "c", time);
    -- create velocity vectors
    vx = approxSpace:create_surface_function();
    vy = approxSpace:create_surface_function();
    -- fill velocity vectors
    lsDisc:fill_v_vec(vx,0);
    lsDisc:fill_v_vec(vy,1);
    -- set velocity vectors
    lsDisc:set_vel_x(vx);
    lsDisc:set_vel_y(vy);
end;
if (setup==3) then
    solfunctor = util.CreateLuaUserNumber("ConstantVSolution", dim);
	lsDisc:set_dirichlet_data(solfunctor);
	InterpolateFunction(solfunctor, phiOld, "c", time);
	lsDisc:set_vel_x(1);
	lsDisc:set_vel_y(1);
end;

lsDisc:set_limiter(false);
lsDisc:compute_error(phiOld);

VecAssign(phiNew,phiOld);

NumTimeSteps =  util.GetParamNumber("-numTimeSteps", 400)
EndTime = 2*math.pi;
lsDisc:set_dt(EndTime/NumTimeSteps);

if (movie==true) then
    -- filename
    filename = "lssol"

    -- write start solution
    print("Writing start values")
    out = VTKOutput()
	step=0;
	time=lsDisc:get_time();
    out:print(filename, phiOld, step, time)
end;

-- lsDisc:create_ls_subsets(phiOld);
tBefore = os.clock()

for step = 1, NumTimeSteps do
    lsDisc:advect_lsf(phiNew,phiOld);
	VecAssign(phiOld,phiNew);
   --lsDisc:compute_error(phiNew);
	if (movie==true) then
    	time=lsDisc:get_time();
    	out:print(filename, phiNew, step, time)
	end;
end
tAfter = os.clock()
print("ls stepping took " .. tAfter-tBefore .. " seconds");
lsDisc:compute_error(phiNew);
--------------------------------------------------------------------------------
--  Output
--------------------------------------------------------------------------------
if (movie==true) then
    out:write_time_pvd(filename, phiNew);
end;

if (writetovtk==true) then
    WriteGridFunctionToVTK(phiNew, "Solution")
end;

SaveVectorForConnectionViewer(phiNew, "u.vec")
