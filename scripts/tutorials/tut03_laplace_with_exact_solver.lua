--------------------------------------------------------------------------------
--	tut03_laplace_with_exact_solver.lua
--
--	Based on the preceeding tutorial, we will now solve a laplace equation on
--	a given domain using an exact lu-decomposition.
--
--	Note that the code runs in a parallel evironment too. However the exact
--	LU decomposition can only be used for serial environments. One would
--	thus have to replace the solver before the code can be executed in parallel.
--	See the next tutorials for an overview of alternate solvers.
--
--	The new code can be found below the comment NEW CODE STARTS HERE.
--------------------------------------------------------------------------------

-- We will include a script file which defines some methods often used.
-- Loaded methods are all found in the library util.
ug_load_script("../ug_util.lua")

-- To keep the script flexible we will now define some variables which have
-- a default value that can be overwritten by command line arguments.

-- Depending on the dimension we will choose our domain object
-- (either 1d, 2d or 3d) and associated discretization objects. Note that
-- we're using some methods defined in "ug_util.lua" here.
dim = util.GetParamNumber("-dim", 2) -- default dimension is 2.

-- Since ug supports a bunch of different dimensions and algebra modules 
-- we will choose a combination here. This should always be the first thing 
-- you do in an ug-script. The cpu-algebra is fine for now.
InitUG(dim, AlgebraType("CPU", 1))

-- We also need a filename for the grid that shall be loaded.
gridName = util.GetParam("-grid", "unit_square/unit_square_quads_8x8.ugx")

-- Since we want to save the domains grid to a file, we also need an output file.
-- Note that we only use a prefix here, since we want to attach the process number
-- and file format ourselfs
outFileNamePrefix = util.GetParam("-o", "distributed_domain_")




-- Now its time to create the domain object. We will use an util method here,
-- which automatically creates the right domain for the given dimension.
dom = Domain()

-- Now that we have a domain, we can load a grid into it. We check the return
-- value whether the loading was successful or not.
-- Note that we use the method utilLoadDomain instead of LoadDomain. utilLoadDomain
-- has the benefit that grids are automatically searched in the data/grids folder if
-- they were not found at the default locations (execution-path or a path specified
-- in your environments path-variable).
if LoadDomain(dom, gridName) == false then
	print("Loading of domain " .. gridName .. " failed. Aborting.")
--	call exit to leave the application right away.
	exit() 
end

-- Now that we're here, the domain was successfully loaded
print("Loaded domain from " .. gridName)


-- Distribute the domain to all involved processes
if DistributeDomain(dom) == false then
	print("Error while Distributing Domain. Aborting.")
	exit()
end


-- Lets save the domain on each process
outFileName = outFileNamePrefix .. GetProcessRank() .. ".ugx"
if SaveDomain(dom, outFileName) == false then
	print("Saving of domain to " .. outFileName .. " failed. Aborting.")
	exit()
end

-- Everything seems to went fine.
print("Saved domain to " .. outFileName)



--------------------------------------------------------------------------------
--	NEW CODE STARTS HERE
--------------------------------------------------------------------------------

-- Before we proceed, we want to make sure that the domain is ready for our
-- problem. Our geometry consists of two main parts: an inner part and a
-- boundary part. Those parts are called "subsets". Those subsets are handled
-- by so called "subset handlers". Every domain object contains such a subset
-- handler. It is filled during LoadDomain. Since our algorithms later on
-- require that an "Inner" and a "Boundary" subset exists, we check here
-- for their existance and exit if they are not present.
-- We obtain the subset handler from the domain-object. Note that we use
-- the : operator to access member methods of an object.
sh = dom:subset_handler()
if sh:get_subset_index("Inner") == -1 then
	print("Domain does not contain subset 'Inner'. Aborting.")
	exit()
end

if sh:get_subset_index("Boundary") == -1 then
	print("Domain does not contain subset 'Boundary'. Aborting.")
	exit()
end

-- For the laplace problem we need a right hand side (rhs), a diffusion tensor.
-- and dirichlet boundary values. Several methods exist how to specify those.
-- In this example we will use script callback functions to do this. One could
-- instead use constant values which are faster but less flexible (obviously!)
-- or one could plug exported quantities from other equations into those slots.
-- This is an advanced and powerful feature since it allows the coupling of
-- different equations without touching a line of source-code. We will examine
-- the coupling process in later tutorials.

-- the diffusion tensor for two and three dimensions.
-- Parameters are the coordinates of the point at which the tensor should be
-- evaluated and the time of the current time-step.
-- Note that we simply return multiple values - the coefficients of the diffusion matrix.
function ourDiffTensor2d(x, y, t)
	return	1, 0, 
			0, 1
end

function ourDiffTensor3d(x, y, z, t)
	return	1, 0, 0,
			0, 1, 0,
			0, 0, 1
end

-- Next we need the dirichlet boundary
-- Parameters are the coordinates of the point at which the tensor should be
-- evaluated and the time of the current time-step.
-- The dirichlet boundary callback has to return two values: a boolean
-- defining whether the given point really is a dirichlet boundary point
-- and the boundary value itself.
-- Note that we use math-functions from luas standard math-library.
-- (here the . operator is used, since math is not an object but a library)
function ourDirichletBnd2d(x, y, t)
	local s = 2*math.pi
	return true, math.sin(s*x) + math.sin(s*y)
end

function ourDirichletBnd3d(x, y, z, t)
	local s = 2*math.pi
	return true, math.sin(s*x) + math.sin(s*y) + math.sin(s*z)
end

-- Finally we have to specify the right hand side. Make sure that
-- it is compatible with the specified boundary values.
-- This callback again takes the coordinates of the evaluated point
-- and the current time step as arguments. It has to return the value
-- of the rhs at the given point at the given time.
function ourRhs2d(x, y, t)
	local s = 2*math.pi
	return	s*s*(math.sin(s*x) + math.sin(s*y))
end

function ourRhs3d(x, y, z, t)
	local s = 2*math.pi
	return	s*s*(math.sin(s*x) + math.sin(s*y) + math.sin(s*z))
end


-- Now that all callbacks are set up, we will create an Approximation Space
-- which describes the unknowns in the equation. In our case we will
-- only require one unknown for the concentration ("c").
-- Note that the Approximation Space is build on the domain created above.
print("Create ApproximationSpace")
approxSpace = ApproximationSpace(dom) -- creates new object
approxSpace:add_fct("c", "Lagrange", 1)          -- adds one function


-- Now we will create the discretization objects. This is performed in
-- several steps to guarantee maximal flexibility. First we will create
-- element discretizations. Those perform the discretization on the
-- actual elements. UG features several differen discretization schemes.
-- We use the convection / diffusion scheme with finite volumes.
-- Note that we do not set a convection callback, which will result in
-- a default convection of 0 - the convection / diffusion equation thus
-- simplifies to a pure diffusion equation - the laplace problem.

-- The element discretization
-- Here we create a new instance of a convection diffusion equation.
-- We furthermore tell it that it shall operate on the unknowns associated
-- with function "c" (the concentration defined in the approximation space)
-- and that the discretization shall only operate on elements in the subset
-- "Inner". Note that one could specify multiple subsets here by enumerating
-- them all in the subsets-string separated by , (i.e. "Inner1, Inner2").

-- Select upwind
--upwind = NoUpwind()
--upwind = FullUpwind()
upwind = WeightedUpwind(); upwind:set_weight(0.0)
--upwind = PartialUpwind()
 
elemDisc = ConvectionDiffusion("c", "Inner")
elemDisc:set_disc_scheme("fv1")
elemDisc:set_upwind(upwind)
elemDisc:set_diffusion_tensor("ourDiffTensor"..dim.."d")	-- set the diffusion matrix
elemDisc:set_source("ourRhs"..dim.."d")						-- set the right hand side

-- Note that the dirichlet boundary callback was not registered at the
-- element discretization. This is important since the concept of
-- element discretizations and boundary values are separated in ug4.
-- This means that an element discretization only assembles those
-- parts of the final stiffness matrix, which describe the relation between
-- the elements of the given subset. The boundary conditions are then
-- enforced either by other discretization objects (neumann) or through
-- a post-process (dirichlet).

-- first we create objects that encapsulate our callbacks. Those can then
-- be registered at the discretization object. Note that we use the .. operator
-- to concatenate strings and numbers. This saves us from a lot of if dim == 2 ... else ...
-- For the dirichlet callback we use utilCreateLuaBoundaryNumber, since here
-- a boolean and a number are returned.
dirichletCallback = LuaBoundaryNumber("ourDirichletBnd" .. dim .. "d")

-- Here we set up such a dirichlet boundary condition. We explicitly
-- add subsets on which the boundary callback defined above shall be
-- applied to a given function (as defined in the approximation space).
dirichletBnd = DirichletBoundary()
dirichletBnd:add(dirichletCallback, "c", "Boundary")

-- Finally we create the discretization object which combines all the
-- separate discretizations into one domain discretization.
domainDisc = DomainDiscretization(approxSpace)
domainDisc:add(elemDisc)
domainDisc:add(dirichletBnd)


-- The discretization itself is now completely prepared. It is now time
-- to assemble a matrix from it. To keep things general, we do not
-- create a matrix but a linear operator which represents an assembled system
-- (in most cases of course this will be a matrix).
linOp = AssembledLinearOperator()
-- the discretization object from which the operator is assembled
linOp:set_discretization(domainDisc)

-- Now lets solve the problem. Create a vector of unknowns and a vector
-- which contains the right hand side. We will use the approximation space
-- to create the vectors. Make sure to create the vector for the same
-- dofs as set to linOp through linOp:set_dof_distribution.
-- Note that the vectors automatically have the right size.
u = GridFunction(approxSpace)
b = GridFunction(approxSpace)


-- the linear operator is now complete. To perform the discretization call init.
-- Here we assemble the rhs in one loop to save resources.
linOp:init_op_and_rhs(b)

-- Init the vector representing the unknowns with 0 to obtain
-- reproducable results.
u:set(0)

-- We need a solver to solve the system. Since we want an exact solution,
-- we simply use the LU solver.
-- But be careful... if we are in a parallel environment, the LU decomposition
-- will not work. We use GetNumProcesses() to retrieve the number of active processes.

if GetNumProcesses() > 1 then
	print("Can't apply LU decomposition in parallel environment. Aborting.")
	exit()
end

-- Create the solver and initialize it with the linear operator which
-- shall be inverted.
solver = LU()
solver:init(linOp)

-- The solver is ready. We can now apply it.
solver:apply(u, b)
-- Note that one could also call solver:apply_return_defect(u, b). In this
-- case b would contain the defect after the solution.


-- Finally we're nearly done. The only thing left to do is to write the
-- solution to a file which can then be examined using e.g. Paraview.
-- (Open "Solution_t0000.pvd" in paraview to view the complete domain
-- at timestep 0 (we solved a stationary problem))
WriteGridFunctionToVTK(u, "Solution")

print("")
print("done.")

