------------------------------------------------------------------
--
--   Lua - Script to test the Navier-Stokes implementation
--
--	This script sets up a problem for the Navier-Stokes discretization
--	and is intended to test the implementation.
--
--   Author: Josef Dubsky, Andreas Vogel
--
-------------------------------------------------------------------

--------------------------------
--------------------------------
-- Parameter setup
--------------------------------
--------------------------------

-- Right at the beginning we load a lot of util functions, that help us basically
-- to programm a domain independent lua-script and provide some basic helper
-- functions. Most the util functions are in the namespace util, i.e. they
-- are used by 'util.SomeFunction'
ug_load_script("ug_util.lua")

-- Next, we have to choose some algebra. ug4 provides several algebra, e.g.
-- there are block structured matrices or simple double-valued matrices. We
-- decide to use the double-valued CSR Matrix. This is the default case for the
-- Algebra chooser and so we leave the intiallizer of the AlgebraChooser empty.
algChooser = CPUAlgebraSelector()
algChooser:set_fixed_blocksize(3)
InitAlgebra(algChooser);

-- Depending on the dimension we will choose our domain object
-- (either 1d, 2d or 3d) and associated discretization objects. Note that
-- we're using some methods defined in "ug_util.lua" here. The dimesion is
-- read in from the bash command line passing e.g. the option "-dim 2". If no
-- dim-option is specified the default value in the second argument of
-- GetParamNumber is used.
dim = util.GetParamNumber("-dim", 2) -- default dimension is 2.

-- Next, we decide which grid to use. This can again be passed as a command line
-- option or a default value is used.
if 	dim == 2 then gridName = util.GetParam("-grid", "cylinder.ugx")
else print("Choosen Dimension not supported. Exiting."); exit(); end

-- We additionally use parameters which allow to specify the number of
-- pre- and total-refinement steps (wrt domain distribution).
numPreRefs = util.GetParamNumber("-numPreRefs", 0)
numTotalRefs = util.GetParamNumber("-numTotalRefs", 0)

-- Calculate the number of post-refs and make sure that the result makes sense.
numPostRefs = numTotalRefs - numPreRefs
if numPostRefs < 0 then
	print("WARNING:\tnumPreRefs exceeds the number of total refinements.")
	print("\t\t\tNo refinement will be preformed after distribution.")
	numPostRefs = 0
end

-- Lets write some info about the choosen parameter
print(" Choosen Parater:")
print("    dim        	= " .. dim)
print("    numTotalRefs = " .. numTotalRefs)
print("    numPreRefs 	= " .. numPreRefs)
print("    grid       	= " .. gridName)

--------------------------------------------
--------------------------------------------
-- Loading Domain and Domain Refinement
--------------------------------------------
--------------------------------------------

-- Create the domain and load a grid
dom = util.CreateDomain(dim)

if util.LoadDomain(dom, gridName) == false then
	print("Loading of domain " .. gridName .. " failed. Aborting.")
	exit()
end

-- Now that we're here, the domain was successfully loaded
print("Loaded domain from " .. gridName)

-- We will create a refiner now. This is a tool associated with a domain.
-- UG defines factory methods for refiners, which automatically choose
-- the right refiner for the given context, i.e. different refiners are
-- created depending on whether we are in a parallel or in a serial environment.
-- Note that another factory method is HangingNodeDomainRefiner, which is
-- subject to a later tutorial.
refiner = GlobalDomainRefiner(dom)

-- perform pre-refinement
for i = 1, numPreRefs do
	refiner:refine()
end

-- Distribute the refined domain to all involved processes
if DistributeDomain(dom) == false then
	print("Error while Distributing Domain. Aborting.")
	exit()
end

-- perform post-refinement
for i = 1, numPostRefs do
	refiner:refine()
end

--------------------------------
--------------------------------
-- Approximation Space
--------------------------------
--------------------------------

-- We succesfully loaded and refined the domain. Now its time to setup an
-- ApproximationSpace on the domain. First, we check that the domain we use
-- has suitable subsets. This are parts of the domain, that partition the domain.
-- We need them, to handle e.g. boundary conditions on different parts of the
-- domain boundary.

-- Lets define a list of all subsets that we need
neededSubsets = {"Inner", "Inlet", "Outlet", "UpperWall", "LowerWall", "CylinderWall"}

-- Now we loop all subsets an search for it in the SubsetHandler of the domain
if util.CheckSubsets(dom, neededSubsets) == false then print("Wrong subsets detected.") end

-- All subset are ok. So we can create the Approximation Space
approxSpace = util.CreateApproximationSpace(dom)

-- we add the components of the velocity as Lagrange Ansatz functions of first order
if dim >= 1 then approxSpace:add_fct("u", "Lagrange", 1) end
if dim >= 2 then approxSpace:add_fct("v", "Lagrange", 1) end
if dim >= 3 then approxSpace:add_fct("w", "Lagrange", 1) end

-- we add the pressure as Lagrange Ansatz function of first order
approxSpace:add_fct("p", "Lagrange", 1)

-- we init the approximation space. This fixes the dof pattern used and afterwards
-- no function can be added or removed. This is needed to have fast access and
-- less checks during the computation phase.
approxSpace:init()

-- finally we print some statistic on the distributed dofs
approxSpace:print_statistic()

--------------------------------
--------------------------------
-- Discretization
--------------------------------
--------------------------------

-- We set up the discretization. The basic idea is to first create the single
-- components of the discretization (e.g. element discretization, boundary
-- conditions, ...) and then glue everything together in one DomainDiscretization
-- object that does nothing else than grouping the single components. When it
-- comes to assembling this object runs all scheduled discretization component
-- one after the other.

-- Lets begin with the element discretization. This is the core part for the
-- Navier-Stokes equation and provides the local stiffness- and mass matrices,
-- that are then added to the global matrices. First we concatenate a string
-- that holds all function we use. This string must be passed to the navier
-- stokes elem disc to indicate which function is used for velocity and which
-- one for pressure. E.g. for dim==2 this gives: "u, v, p".
fctUsed = "u"
if dim >= 2 then fctUsed = fctUsed .. ", v" end
if dim >= 3 then fctUsed = fctUsed .. ", w" end
fctUsed = fctUsed .. ", p"

-- Now we create the Navier-Stokes disc. Here we use a util function, that
-- creates the disc for the right world dimension according to the loaded domain
-- (this is handled internally, since the approxSpace knows the domain and the
-- domain knows the dimension). We tell the elem Disc to assemble on all elements
-- that are in the subset "Inner". If we would have several domains, where we
-- would like to do the same, this could be done by passing a list of subsets
-- separated by ',', (e.g. "Inner1, Inner2, Inner3").
elemDisc = util.CreateFV1NavierStokes(approxSpace, fctUsed, "Inner")

-- Now, we have to setup the stabilization, that is used for the Continuity Equation.
-- The stabilization is passed to the Navier-Stokes elem disc as an object.
-- Therefore, we created it now and will pass it to the disc. But first, we have
-- to create the Upwind scheme, that is used inside the stabilization. There are
-- several possibilities:

noUpwind = util.CreateNavierStokesNoUpwind(dim);
fullUpwind = util.CreateNavierStokesFullUpwind(dim);
skewedUpwind = util.CreateNavierStokesSkewedUpwind(dim);
LPSUpwind = util.CreateNavierStokesLinearProfileSkewedUpwind(dim);
POSUpwind = util.CreateNavierStokesPositiveUpwind(dim);

-- Now, we create the stabilization ...
fieldsStab = util.CreateNavierStokesFIELDSStabilization(dim)

-- ... and set the upwind
fieldsStab:set_upwind(fullUpwind)

-- We also can choose, how the diffusion length of the stabilization is computed.
-- Under the option we pick on:
fieldsStab:set_diffusion_length("NS_RAW")
--fieldsStab:set_diffusion_length("NS_FIVEPOINT")
--fieldsStab:set_diffusion_length("NS_COR")

-- Next we set the options for the Navier-Stokes elem disc ...
elemDisc:set_stabilization(fieldsStab)
elemDisc:set_conv_upwind(fullUpwind)
elemDisc:set_peclet_blend(false)
elemDisc:set_exact_jacobian(false)

-- ... and finally we choose a value for the kinematic viscosity.
ConstKinViscosity = util.CreateConstUserNumber(1.0e-1, dim)
elemDisc:set_kinematic_viscosity(ConstKinViscosity);


-- Next, lets create the boundary conditions. We will use lua-callback functions
-- to implement the values of the dirichlet-bnd conditions. So lets define
-- some functions in lua: Parameters are the coordinates of the point at which
-- the value should be evaluated and the time of the current time-step.
-- The dirichlet boundary callback has to return two values: a boolean
-- defining whether the given point really is a dirichlet boundary point
-- and the boundary value itself.
-- Note that we use math-functions from luas standard math-library.
-- (here the . operator is used, since math is not an object but a library)
function inletVelX2d(x, y, t)
	local H = 0.41
	local Um = 0.3
	return true, 4 * Um * y * (H-y) / (H*H)
end

function inletVelY2d(x, y, t)
	local H = 0.41
	local Um = 0.2
	return true, Um * math.sin(y/H*2*3.1415)
end

-- Next, we create objects that encapsulate our callback. Those can then
-- be registered at the discretization object. Note that we use the .. operator
-- to concatenate strings and numbers. This saves us from a lot of if dim == 2 ... else ...
-- For the dirichlet callback we use utilCreateLuaBoundaryNumber, where
-- a boolean and a number are returned.
LuaInletVelX2d = util.CreateLuaBoundaryNumber("inletVelX" .. dim .. "d", dim)
LuaInletVelY2d = util.CreateLuaBoundaryNumber("inletVelY" .. dim .. "d", dim)
ConstZeroDirichlet = util.CreateConstBoundaryNumber(0.0, dim)


-- We need neumann boundary for the pressure at Inlet
function inletPressure2d(x, y, t)
	local b, vel = inletVelX2d(x,y,t)
	return true, -1.0 * vel
end

LuaNeumannPressure = util.CreateLuaBoundaryNumber("inletPressure"..dim.."d", dim)
neumannDisc = util.CreateNeumannBoundary(approxSpace, "Inner")
neumannDisc:add(LuaNeumannPressure, "p", "Inlet")

-- Now we create a Dirichlet Boundary object. At this object all boundary conditions
-- are registered.
dirichletBnd = util.CreateDirichletBoundary(approxSpace)
dirichletBnd:add(LuaInletVelX2d, "u", "Inlet")
dirichletBnd:add(LuaInletVelY2d, "v", "Inlet")
dirichletBnd:add(ConstZeroDirichlet, "u", "UpperWall,LowerWall,CylinderWall")
dirichletBnd:add(ConstZeroDirichlet, "v", "UpperWall,LowerWall,CylinderWall")
dirichletBnd:add(ConstZeroDirichlet, "p", "Outlet")

-- Finally we create the discretization object which combines all the
-- separate discretizations into one domain discretization.
domainDisc = DomainDiscretization()
domainDisc:add(elemDisc)
domainDisc:add(neumannDisc)
domainDisc:add(dirichletBnd)

--------------------------------
--------------------------------
-- Operator
--------------------------------
--------------------------------

-- In ug4 we use Operator interfaces. An operator is simply a mapping from in
-- space into the other. A frequently used implementation of such a mapping is
-- the usage of discretizations, that map a solution (grid function) onto some
-- right-hand side. We can create an operator that uses the recently created
-- domainDisc.
op = AssembledOperator()

-- the discretization object from which the operator is assembled
op:set_discretization(domainDisc)
-- since we do not use a multi-grid method here, we want to operate on the
-- unknowns (degrees of freedom - dofs) of the surface grid. Since we
-- didn't refine in this example, this of course is the same as the base grid.
op:set_dof_distribution(approxSpace:get_surface_dof_distribution())

-- the operator is now complete. To perform the discretization call init.
op:init()

--------------------------------
--------------------------------
-- Solution of the Problem
--------------------------------
--------------------------------

-- Now lets solve the problem. Create a vector of unknowns and a vector
-- which contains the right hand side. We will use the approximation space
-- to create the vectors. Make sure to create the vector for the same
-- dofs as set to linOp through linOp:set_dof_distribution.
-- Note that the vectors automatically have the right size.
u = approxSpace:create_surface_function()

-- Init the vector representing the unknowns with 0 to obtain
-- reproducable results.
u:set(0)

-- We could also interpolate some user defined function
-- setup the lua functions ...
function Pressure_StartValue2d(x, y, t)
	return 0.0
end
function VelX_StartValue2d(x, y, t)
	return 0.0
end
function VelY_StartValue2d(x, y, t)
	return 0.0
end

-- ... and wrap the lua-callback
if dim == 2 then
LuaPressureStartValue = util.CreateLuaUserNumber("Pressure_StartValue"..dim.."d", dim)
LuaVelXStartValue = util.CreateLuaUserNumber("VelX_StartValue"..dim.."d", dim)
LuaVelYStartValue = util.CreateLuaUserNumber("VelY_StartValue"..dim.."d", dim)
end

-- Now interpolate the function
time = 0.0
InterpolateFunction(LuaPressureStartValue, u, "p", time);
InterpolateFunction(LuaVelXStartValue, u, "u", time);
InterpolateFunction(LuaVelYStartValue, u, "v", time);

-- we need a linear solver that solves the linearized problem inside of the
-- newton solver iteration. So, we create an exact LU solver here.
linSolver = LU()


-- Next we need a convergence check, that computes the defect within each
-- newton step and stops the iteration when a specified creterion is fullfilled.
-- For our purpose is the StandardConvergenceCheck is sufficient. Please note,
-- that this class derives from a general IConvergenceCheck-Interface and
-- also more specialized or self-coded convergence checks could be used.
newtonConvCheck = StandardConvergenceCheck()
newtonConvCheck:set_maximum_steps(20)
newtonConvCheck:set_minimum_defect(5e-8)
newtonConvCheck:set_reduction(1e-10)
newtonConvCheck:set_verbose_level(true)

-- Within each newton step a line search can be applied. In order to do so an
-- implementation of the ILineSearch-Interface can be passed to the newton
-- solver. Here again we use the standard implementation.
newtonLineSearch = StandardLineSearch()
newtonLineSearch:set_maximum_steps(20)
newtonLineSearch:set_lambda_start(1.0)
newtonLineSearch:set_reduce_factor(0.5)
newtonLineSearch:set_accept_best(true)

-- Sometimes its helpful to write the defect and jacobian of the newton step
-- to debug the implementation. For that, we use the debug writer
dbgWriter = util.CreateGridFunctionDebugWriter(dim)
dbgWriter:set_reference_grid_function(u)
dbgWriter:set_vtk_output(false)

-- Now we can set up the newton solver. We set the linear solver created above
-- as solver for the linearized problem and set the convergence check. If you
-- want to you can also set the line search.
newtonSolver = NewtonSolver()
newtonSolver:set_linear_solver(linSolver)
newtonSolver:set_convergence_check(newtonConvCheck)
newtonSolver:set_line_search(newtonLineSearch)
newtonSolver:set_debug(dbgWriter)

-- Finally we set the non-linear operator created above and initiallize the
-- Newton solver for this operator.
newtonSolver:init(op)


-- In a first step we have to prepare the newton solver for the solution u. This
-- sets e.g. dirichlet values in u.
if newtonSolver:prepare(u) == false then
	print ("Newton solver prepare failed."); exit();
end


SaveVectorForConnectionViewer(u, "StartSolution.mat")

-- Now we can apply the newton solver. A newton itertation is performed to find
-- the solution.
if newtonSolver:apply(u) == false then
	 print ("Newton solver apply failed."); exit();
end

-- Finally we're nearly done. The only thing left to do is to write the
-- solution to a file which can then be examined using e.g. Paraview.
-- (Open "Solution.vtu" in paraview to view the complete domain
vtkWriter = util.CreateVTKWriter(dim)
vtkWriter:select_all(false)
vtkWriter:select_nodal_vector("u,v", "velocity")
vtkWriter:select_nodal_scalar("p", "pressure")
vtkWriter:print("Solution", u)

print("done.")

