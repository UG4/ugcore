--------------------------------------------------------------------------------
--
--   Lua - Script for scalability studies solving the
--   Laplace-Problem.
--
--   For some flexibility several definitions concerning
--   the solver can also be given via command line parameters
--
--   Authors: Ingo Heppner, Sebastian Reiter, Andreas Vogel
--
--------------------------------------------------------------------------------

--------------------------------------------------------------------------------
-- Execute on cekon via e.g.
-- salloc -n 32 mpirun ./ugshell  -ex ../scripts/tests/scalability_test.lua -dim 2 -grid unit_square/unit_square_quads_8x8.ugx -numPreRefs 4 -numRefs 8

-- Execute on JuGene via e.g. ('mpirun' call has to be specified in a LoadLeveler script!):
-- mpirun -np 32 -exe ./ugshell -mode VN -mapfile TXYZ -verbose 2 -env LD_LIBRARY_PATH=/bgsys/drivers/ppcfloor/comm/lib/ -args "-ex ../scripts/tests/scalability_test.lua -dim 2 -grid unit_square/unit_square_quads_8x8.ugx -numPreRefs 4 -numRefs  8"
-- Or the same job interactively:
-- llrun -v -np 32 -exe ./ugshell -mode VN -mapfile TXYZ -verbose 2 -env LD_LIBRARY_PATH=/bgsys/drivers/ppcfloor/comm/lib/ -args "-ex ../scripts/tests/scalability_test.lua -dim 2 -grid  unit_square/unit_square_quads_8x8.ugx -numPreRefs 4 -numRefs 8"

--------------------------------------------------------------------------------
PrintBuildConfiguration()

ug_load_script("ug_util.lua")
ug_load_script("domain_distribution_util.lua")

dim = util.GetParamNumber("-dim", 2)


-- This is an ugly hack. Don't do that at home!!!
ddu.numProcesses = util.GetParamNumber("-numProcs", GetNumProcesses())


-- Here all parameters related to refinement and distribution are parsed.
-- -numPreRefs, -numRefs, -distType, -numPPN, -hRedistFirstLevel,
-- -numRedistNewProcsPerStep, -hRedistStepSize.
-- Check the documentation of domain_distribution_util.lua for more information.
-- distribution related parameters
if(ddu.ParseAndInitializeParameters(dim) == false) then
	print("An error occured during ddu.ParseAndInitializeParameters. Aborting")
	exit()
end

print(" Parallelisation related parameters chosen:")
print("    numProcs   = " .. ddu.numProcesses)
ddu.PrintParameters("    ")

print("")
print("Now performe a 'dry run' of the hierarchical redistribution approach ...")
-- now print the single steps which will be performed...
ddu.PrintSteps()
print("... 'dry run' finished!")
