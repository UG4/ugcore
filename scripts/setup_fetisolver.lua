--[[
   TODO:
   * Resultate von Martin Rupp mit "famg" als Neumann- und Dirichlet-Problem-Solver in 'numrefs?.txt' reproduzieren:

   UGARGS="-ex ../scripts/tests/scalability_test.lua -dim 2 -grid ../data/grids/unit_square_01/unit_square_01_quads_8x8.ugx -lsMaxIter 100 -numPreRefs 3 -lsType feti -numPPSD 4"

   salloc -n 64 mpirun ugshell $UGARGS -numRefs  7 -ds famg -ns famg -logtofile bla -rlf
   salloc -n 64 mpirun ugshell $UGARGS -numRefs  8 -ds famg -ns famg -logtofile bla -rlf - haengt!?

   Noch nicht probiert:
   salloc -n 64 mpirun ugshell $UGARGS -numRefs  9 -ds famg -ns famg -logtofile bla -rlf
   salloc -n 64 mpirun ugshell $UGARGS -numRefs 10 -ds rsamg -ns rsamg -logtofile bla -rlf


   Testweise der hierzu kleinste sinnvolle Job - alles andere stuerzt ab:
   cekon:
   salloc -n 16 mpirun ugshell $UGARGS -numRefs  5 -ds famg -ns famg -logtofile bla -rlf

   local:
   openmpirun -np 16   ugshell $UGARGS -numRefs  5 -ds famg -ns famg -logtofile bla -rlf

   ==> ERROR in 'LagrangeDirichletBoundary:extract_scheduled_data': Approximation Space not set.

   falls 'CreateAMGTestvectorDirichlet0()' zur Erzeugung der Testvektoren verwendet wird?

   Ursache fuer den Fehler ist offenbar, dass die Methode 'LagrangeDirichletBoundary::set_approximation_space()'
   (ugbase/lib_disc/spatial_disc/constraints/dirichlet_boundary/lagrange_dirichlet_boundary.h)
   nirgends aufgerufen wird!

   Dies wird aber durch 'linOp:init_op_and_rhs(b)', das ich in 'scalability_test.lua' mal testweise, statt 'AssembleLinearOperatorRhsAndSolution(linOp, u, b)', aufrufe -
   und siehe da, kein Absturz mehr!
   Der obige Fehler aber nicht, dazu muesste die Reihenfolge umgedreht werden: Erst Assemblieren, dann Aufsetzen FAMG mit Testvektoren!

   ]]
----------------------------------------------------------
--
--   Lua - Script which creates and configures a FETI solver (linear solver).
--
--   Author: Ingo Heppner (using FETI solver specific code from 'fetitest.lua'.
--   Begin: 06072011.
--
--   Provides: Function which set up a FETI solver object.
--
--          Input parameters:
--          'str_problem' -- string describing problem (not yet used)
--          'dim',
--          'linMaxIterations'
--          'numProcs'
--          'dirichletBND', 'approxSpace',
--          'activateDbgWriter', 'verbosity', logfileName.
--          Return value:
--          Reference of the fully configured FETI solver object.
--
--   Description of some parameters / options:
--      -AMGwriteMat:   write testvectors (only if Neumann or Dirichlet problem solver is of type "famg")

--   Usage example:
--	1. In the calling LUA script:
--
--	ug_load_script("setup_fetisolver.lua")
--	solver = SetupFETISolver(str_problem,
--				 dim
--				 linMaxIterations,
--				 numProcs,
--				 dirichletBND, approxSpace, -- for testvector writer for FAMG (created by 'CreateAMGTestvectorDirichlet0()'); 'approxSpace' also for 'GridFunctionDebugWriter()'
--				 activateDbgWriter,
--				 verbosity, logfileName)
--
--      where 'solver' is the solver object which is called by
--      'ApplyLinearSolver(., ., ., <solver>)', i.e.:
--	
--	print("Apply solver.")
--	if ApplyLinearSolver(linOp, u, b, solver) == false then
--		print("Could not apply linear solver.");
--	end
--
--	2. Execution of the calling LUA script:
--[[

# Local:
########
setenv UGARGS "-ex ../scripts/tests/scalability_test.lua -dim 2 -grid ../data/grids/unit_square_01/unit_square_01_quads_8x8.ugx -lsMaxIter 100 -numPreRefs 0 -lsType feti -numPPSD 1"
openmpirun -np 4 ugshell $UGARGS -numRefs 5 # etc.
openmpirun -np 4 ugshell $UGARGS -numRefs 5 -dbgw 1 # with debug writer
	
# cekon:
########
UGARGS="-ex ../scripts/tests/scalability_test.lua -dim 2 -grid ../data/grids/unit_square_01/unit_square_01_quads_8x8.ugx -lsMaxIter 100 -numPreRefs 3 -lsType feti -numPPSD 1"

salloc -n  4 mpirun ./ugshell $UGARGS -numRefs 5 -logtofile feti-sd1_8x8-quad_prerefs3-refs5_pe04.txt - geht!
salloc -n  4 mpirun ./ugshell $UGARGS -numRefs 6 -logtofile feti-sd1_8x8-quad_prerefs3-refs6_pe04.txt - geht, zumindest mit 2000 Schritten Dirichlet-Solver!
salloc -n  4 mpirun ./ugshell $UGARGS -numRefs 7 -logtofile feti-sd1_8x8-quad_prerefs3-refs7_pe04.txt - "Could not solve Dirichlet problem (step 3.b)", auf allen 4 Procs!
salloc -n  4 mpirun ./ugshell $UGARGS -numRefs 7 -ds rsamg - works!
salloc -n  4 mpirun ./ugshell $UGARGS -numRefs 8 -ds rsamg - nun Ausstieg des Neumann-Solvers!

-- Test of solvers:
UGARGS="-ex ../scripts/tests/scalability_test.lua -dim 2 -grid ../data/grids/unit_square_01/unit_square_01_quads_8x8.ugx -lsMaxIter 100 -numPreRefs 3 -lsType feti -numPPSD 1"
-- Default solvers:
salloc -n  4 mpirun ./ugshell $UGARGS -numRefs 7 -ds cg    -ns cg    -cps exact -logtofile laplace_feti-sd1_8x8-quads_prerefs3-refs7_dps-cg_nps-cg_cps-exact_pe04.txt
-- RSAMG for Dirichlet problem:
salloc -n  4 mpirun ./ugshell $UGARGS -numRefs 7 -ds rsamg -ns cg    -cps exact -logtofile laplace_feti-sd1_8x8-quads_prerefs3-refs7_dps-rsamg_nps-cg_cps-exact_pe04.txt
-- RSAMG for Dirichlet+Neumann problem:
salloc -n  4 mpirun ./ugshell $UGARGS -numRefs 7 -ds rsamg -ns rsamg -cps exact -logtofile laplace_feti-sd1_8x8-quads_prerefs3-refs7_dps-rsamg_nps-rsamg_cps-exact_pe04.txt

-- exact for Dirichlet problem - only for one process per subdomain since 'exact' isn't parallelized:
salloc -n  4 mpirun ./ugshell $UGARGS -numRefs 7 -ds exact -ns cg    -cps exact -logtofile laplace_feti-sd1_8x8-quads_prerefs3-refs7_dps-exact_nps-cg_cps-exact_pe04.txt
-- exact for Dirichlet+Neumann problem - only for one process per subdomain since 'exact' isn't parallelized:
salloc -n  4 mpirun ./ugshell $UGARGS -numRefs 7 -ds exact -ns exact -cps exact -logtofile laplace_feti-sd1_8x8-quads_prerefs3-refs7_dps-exact_nps-exact_cps-exact_pe04.txt


salloc -n 16 mpirun ./ugshell $UGARGS -numRefs 6 -logtofile feti-sd1_8x8-quad_prerefs3-refs6_pe16.txt - geht!
salloc -n 16 mpirun ./ugshell $UGARGS -numRefs 7 -logtofile feti-sd1_8x8-quad_prerefs3-refs7_pe16.txt - "Could not solve Dirichlet problem (step 3.b)", nur auf Proc 8!? <-- #DoF's wie bei Klawonn

salloc -n 64 mpirun ./ugshell $UGARGS -numRefs 7 -logtofile feti-sd1_8x8-quad_prerefs3-refs7_pe64.txt - geht!
salloc -n 64 mpirun ./ugshell $UGARGS -numRefs 8 -logtofile feti-sd1_8x8-quad_prerefs3-refs8_pe64.txt - "Could not solve Dirichlet problem (step 3.b)", auf Proc 53, 55, 61, 4, 5, 1 und 16
(Das sind in der Tat alle, bei denen es schief geht - keine weiteren durch Auskommentieren des "return(false)" in der entsprechenden Bedingung!)

#Mit '-distType grid2d':
salloc -n 64 mpirun ./ugshell $UGARGS -numRefs 8 -distType grid2d -logtofile feti-sd1_8x8-quad_prerefs3-refs8_grid2d_pe64.txt - jetzt geht "step 3.b" auf folgenden Procs nicht:
39, 47, 31, 55, 60, 61, 62, 3, 1, 2


################################################################################
# With new naming scheme (27102011):
################################################################################
UGARGS="-ex ../scripts/tests/scalability_test.lua -dim 2 -grid ../data/grids/unit_square_01/unit_square_01_quads_8x8.ugx -lsMaxIter 100 -numPreRefs 3 -lsType feti"
#  -nPPSD 1:
############
salloc          -n   4 mpirun                            ./ugshell $UGARGS -numRefs 5 -dps exact -nps cg -cps exact -logtofile bla -rlf
salloc          -n   4 mpirun                            ./ugshell $UGARGS -numRefs 6 -dps exact -nps cg -cps exact -logtofile bla -rlf - Absturz, wohl wegen Speicherbedarf fuer exact-Solvers
salloc          -n  16 mpirun                            ./ugshell $UGARGS -numRefs 6 -dps exact -nps cg -cps exact -logtofile bla -rlf
salloc          -n  64 mpirun                            ./ugshell $UGARGS -numRefs 7 -dps exact -nps cg -cps exact -logtofile bla -rlf

# Auswertung auf Anzahl FETI-Schritte:
grep " Relative reduction 1" scaltest-laplace-2d_8x8-quads_refs-3-?_feti_nppsd-1_spss-exact-cg-exact_pe-*.txt

salloc          -n   4 mpirun                            ./ugshell $UGARGS -numRefs 6 -dps cg    -nps cg -cps exact -logtofile bla -rlf
salloc          -n  16 mpirun                            ./ugshell $UGARGS -numRefs 6 -dps cg    -nps cg -cps exact -logtofile bla -rlf
salloc          -n  16 mpirun                            ./ugshell $UGARGS -numRefs 7 -dps cg    -nps cg -cps exact -logtofile bla -rlf -- step 3.b failed nur auf Proc 8 (wie bisher)
salloc          -n  64 mpirun                            ./ugshell $UGARGS -numRefs 7 -dps cg    -nps cg -cps exact -logtofile bla -rlf
salloc          -n  64 mpirun                            ./ugshell $UGARGS -numRefs 8 -dps cg    -nps cg -cps exact -logtofile bla -rlf -- TODO! Duerfte wohl ebenfalls an 3.b failen!?
salloc -N 23 -O -n 256 mpirun -mca mpi_yield_when_idle 1 ./ugshell $UGARGS -numRefs 7 -dps cg    -nps cg -cps exact -logtofile bla -rlf
salloc -N 23 -O -n 256 mpirun -mca mpi_yield_when_idle 1 ./ugshell $UGARGS -numRefs 8 -dps cg    -nps cg -cps exact -logtofile bla -rlf
salloc -N 23 -O -n 256 mpirun -mca mpi_yield_when_idle 1 ./ugshell $UGARGS -numRefs 9 -dps cg    -nps cg -cps exact -logtofile bla -rlf -- TODO!

# Auswertung auf Anzahl FETI-Schritte:
grep " Relative reduction 1" scaltest-laplace-2d_8x8-quads_refs-3-?_feti_nppsd-1_spss-cg-cg-exact_pe-*.txt

#  -nPPSD 4:
############
salloc          -n  16 mpirun                            ./ugshell $UGARGS -numRefs 4 -dps cg    -nps cg    -cps exact -nPPSD 4 -logtofile bla -rlf -- Solve Neumann problem failed
salloc          -n  16 mpirun                            ./ugshell $UGARGS -numRefs 5 -dps cg    -nps cg    -cps exact -nPPSD 4 -logtofile bla -rlf -- Solve Neumann problem failed
salloc          -n  16 mpirun                            ./ugshell $UGARGS -numRefs 6 -dps cg    -nps cg    -cps exact -nPPSD 4 -logtofile bla -rlf -- Solve Neumann problem failed

# Da mit diesen Parameter schon 16 PE nicht laufen, das folgende noch gar nicht probiert ...
salloc          -n  64 mpirun                            ./ugshell $UGARGS -numRefs 7 -dps cg    -nps cg    -cps exact -nPPSD 4 -logtofile bla -rlf
salloc -N 23 -O -n 256 mpirun -mca mpi_yield_when_idle 1 ./ugshell $UGARGS -numRefs 7 -dps cg    -nps cg    -cps exact -nPPSD 4 -logtofile bla -rlf
salloc -N 23 -O -n 256 mpirun -mca mpi_yield_when_idle 1 ./ugshell $UGARGS -numRefs 8 -dps cg    -nps cg    -cps exact -nPPSD 4 -logtofile bla -rlf

# ... leider geht das folgende, Neumann-Solver 'exact', nicht, da dieser Loeser nicht parallel ...:
salloc          -n  16 mpirun                            ./ugshell $UGARGS -numRefs 6 -dps cg    -nps exact -cps exact -nPPSD 4 -logtofile bla -rlf -- step 3.b failed
salloc          -n  16 mpirun                            ./ugshell $UGARGS -numRefs 5 -dps cg    -nps exact -cps exact -nPPSD 4 -logtofile bla -rlf -- step 3.b failed, mit 'nans' als last defect! 
salloc          -n  16 mpirun                            ./ugshell $UGARGS -numRefs 4 -dps cg    -nps exact -cps exact -nPPSD 4 -logtofile bla -rlf -- step 3.b failed
salloc          -n  16 mpirun                            ./ugshell $UGARGS -numRefs 3 -dps cg    -nps exact -cps exact -nPPSD 4 -logtofile bla -rlf -- step 3.b failed, mit 'nans' als last defect!
salloc          -n  64 mpirun                            ./ugshell $UGARGS -numRefs 7 -dps cg    -nps exact -cps exact -nPPSD 4 -logtofile bla -rlf -- TODO!
salloc -N 23 -O -n 256 mpirun -mca mpi_yield_when_idle 1 ./ugshell $UGARGS -numRefs 8 -dps cg    -nps exact -cps exact -nPPSD 4 -logtofile bla -rlf -- TODO!
# ... stattdessen sowohl fuer Dirichlet- als auch fuer Neumann-Solver 'exact' gewaehlt:
salloc          -n  16 mpirun                            ./ugshell $UGARGS -numRefs 6 -dps exact -nps exact -cps exact -nPPSD 4 -logtofile bla -rlf -- abgebrochen nach 2 1/2 Schritten (6h)
salloc          -n  64 mpirun                            ./ugshell $UGARGS -numRefs 7 -dps exact -nps exact -cps exact -nPPSD 4 -logtofile bla -rlf -- abgebrochen nach 1 1/2 Schritten ...

################################################################################
# Test, ob S_PiPi identisch fuer konstante Anzahl Subdomains:
################################################################################
UGARGS="-ex ../scripts/tests/scalability_test.lua -dim 2 -grid ../data/grids/unit_square_01/unit_square_01_quads_8x8.ugx -lsMaxIter 100 -numPreRefs 3 -lsType feti"
salloc        -n  16 mpirun ./ugshell $UGARGS -numRefs 3 -ds rsamg -verb 2 -dbgw 1 -numPPSD 1
salloc        -n  64 mpirun ./ugshell $UGARGS -numRefs 3 -ds bicg  -verb 2 -dbgw 1 -numPPSD 4
# mit Overcommit:
salloc -N 4 -O -n 64 mpirun -mca mpi_yield_when_idle 1 ./ugshell $UGARGS -numRefs 3 -ds bicg  -verb 2 -dbgw 1 -numPPSD 4

salloc -N 16 -O -n 256 mpirun -mca mpi_yield_when_idle 1 ./ugshell $UGARGS -numRefs 8 -ds bicg  -verb 2 -dbgw 1 -numPPSD 4
salloc -N 23 -O -n 256 mpirun -mca mpi_yield_when_idle 1 ./ugshell $UGARGS -numRefs 2 -ds bicg  -verb 2 -dbgw 1 -numPPSD 16
==> alle drei Schurkomplement-Matrizen identisch (umbenannt):
RootSchurComplementMatrix_p0000_pe256_sd16_coords.mat RootSchurComplementMatrix_p0000_pe64_sd16_coords.mat RootSchurComplementMatrix_p0000_pe16_sd16_coords.mat


# Vergleich und Test der Ausgaben beim Erzeugen von S_PiPi - identisches Problem, unterschiedliche outprocs:
UGARGS="-ex ../scripts/tests/scalability_test.lua -dim 2 -grid ../data/grids/unit_square_01/unit_square_01_quads_8x8.ugx -lsMaxIter 100 -numPreRefs 3 -lsType feti"
salloc -n 64 mpirun ./ugshell $UGARGS -numRefs 3 -ds bicg  -verb 2 -dbgw 1 -numPPSD 4 -outproc  0 -logtofile feti_pe64_sd16_p00.txt
salloc -n 64 mpirun ./ugshell $UGARGS -numRefs 3 -ds bicg  -verb 2 -dbgw 1 -numPPSD 4 -outproc  1 -logtofile feti_pe64_sd16_p01.txt
salloc -n 64 mpirun ./ugshell $UGARGS -numRefs 3 -ds bicg  -verb 2 -dbgw 1 -numPPSD 4 -outproc  2 -logtofile feti_pe64_sd16_p02.txt
salloc -n 64 mpirun ./ugshell $UGARGS -numRefs 3 -ds bicg  -verb 2 -dbgw 1 -numPPSD 4 -outproc  3 -logtofile feti_pe64_sd16_p03.txt


salloc -n 64 mpirun ./ugshell $UGARGS -numRefs 3 -ds bicg  -verb 2 -dbgw 1 -numPPSD 4 -outproc 12 -logtofile feti_pe64_sd16_p12.txt
salloc -n 64 mpirun ./ugshell $UGARGS -numRefs 3 -ds bicg  -verb 2 -dbgw 1 -numPPSD 4 -outproc 13 -logtofile feti_pe64_sd16_p13.txt
salloc -n 64 mpirun ./ugshell $UGARGS -numRefs 3 -ds bicg  -verb 2 -dbgw 1 -numPPSD 4 -outproc 14 -logtofile feti_pe64_sd16_p14.txt
salloc -n 64 mpirun ./ugshell $UGARGS -numRefs 3 -ds bicg  -verb 2 -dbgw 1 -numPPSD 4 -outproc 15 -logtofile feti_pe64_sd16_p15.txt


# JuGene:
#########
Interaktiv - Job lief jedoch nicht mehr (01092011):
UGARGS="-ex ../scripts/tests/scalability_test.lua -dim 2 -grid ../data/grids/unit_square_01/unit_square_01_quads_8x8.ugx -lsMaxIter 100 -numPreRefs 3 -lsType feti"
llrun -v  -np 16 -exe ./ugshell -mode VN -verbose 2 -env LD_LIBRARY_PATH=/bgsys/drivers/ppcfloor/comm/lib/ -args $UGARGS -numPPSD 1 -outproc 0


UGARGS="-ex ../scripts/tests/scalability_test.lua -dim 2 -grid ../data/grids/unit_square_01/unit_square_01_quads_8x8.ugx -lsMaxIter 100 -numPreRefs 3 -lsType feti -numPPSD 1"
mpirun -np  64 -exe ./ugshell -mode VN -mapfile TXYZ -verbose 2 -env LD_LIBRARY_PATH=/bgsys/drivers/ppcfloor/comm/lib/ \
                       -args "$UGARGS -numRefs 8 -logtofile jugene_feti-sd1_8x8-quad_prerefs3-refs8_pe64.txt"

mpirun -np 256 -exe ./ugshell -mode VN -mapfile TXYZ -verbose 2 -env LD_LIBRARY_PATH=/bgsys/drivers/ppcfloor/comm/lib/ \
                       -args "$UGARGS -numRefs 8 -logtofile jugene_feti-sd1_8x8-quad_prerefs3-refs8_pe256.txt -verb 1"   - geht!

mpirun -np 256 -exe ./ugshell -mode VN -mapfile TXYZ -verbose 2 -env LD_LIBRARY_PATH=/bgsys/drivers/ppcfloor/comm/lib/ \
                       -args "$UGARGS -numRefs 9 -logtofile jugene_feti-sd1_8x8-quad_prerefs3-refs9_pe256.txt" - hier geht "step 3.b" auf folgenden Procs nicht:
95,253,1,117,215,213,127,245,247,223,221,64,65,69,20,21,16,17,4,5
Umgeordnet:
echo "95,253,1,117,215,213,127,245,247,223,221,64,65,69,20,21,16,17,4,5" | tr "," "\012" | sort -n | tr "\012" ","
1,4,5,16,17,20,21,64,65,69,95,117,127,213,215,221,223,245,247,253 - 68, das ebenfalls am linken Rand liegt, taucht nicht auf; 65 liegt am oberen Rand
Automatisch, direkt aus dem Output-File auf JuGene:
grep "Could not solve Dirichlet problem " ug4_laplace_feti.204348.out_feti-sd1_8x8-quad_prerefs3-refs9_bisect_pe256 | awk ' { printf( "%d\n", $2) }' | sort -n | tr "\012" ","
1,4,5,16,17,20,21,64,65,69,95,117,127,213,215,221,223,245,247,253

# Mit '-distType grid2d':
mpirun -np 256 -exe ./ugshell -mode VN -mapfile TXYZ -verbose 2 -env LD_LIBRARY_PATH=/bgsys/drivers/ppcfloor/comm/lib/ \
                       -args "$UGARGS -numRefs 9 -distType grid2d -logtofile jugene_feti-sd1_8x8-quad_prerefs3-refs9_grid2d_pe256.txt -verb 1"
79,176,239,16,111,159,144,96,127,128,191,112,64,143,48,207,32,223,175,80,12,243,11,8,254,253,252,1,3,2,247,244,5,4,7,6,249,250,251,248
Umgeordnet:
echo "79,176,239,16,111,159,144,96,127,128,191,112,64,143,48,207,32,223,175,80,12,243,11,8,254,253,252,1,3,2,247,244,5,4,7,6,249,250,251,248" | tr "," "\012" | sort -n | tr "\012" ","
1,2,3,4,5,6,7,8,11,12,16,32,48,64,79,80,96,111,112,127,128,143,144,159,175,176,191,207,223,239,243,244,247,248,249,250,251,252,253,254
Automatisch, direkt aus dem Output-File auf JuGene:
grep "Could not solve Dirichlet problem " ug4_laplace_feti.204720.out_feti-sd1_8x8-quad_prerefs3-refs9_grid2d_pe256 | awk ' { printf( "%d\n", $2) }' | sort -n | tr "\012" ","
1,2,3,4,5,6,7,8,11,12,16,32,48,64,79,80,96,111,112,127,128,143,144,159,175,176,191,207,223,239,243,244,247,248,249,250,251,252,253,254

  Number of DoF's, starting from coarse grid 8 x 8 quadrilaterals:
  Level  0 ==> (2^3 * 2^{ 0} + 1)^2 nodes = (2^{ 3} + 1)^2 nodes =              81 nodes
  Level  1 ==> (2^3 * 2^{ 1} + 1)^2 nodes = (2^{ 4} + 1)^2 nodes =             289 nodes
  Level  2 ==> (2^3 * 2^{ 2} + 1)^2 nodes = (2^{ 5} + 1)^2 nodes =           1'089 nodes
  Level  3 ==> (2^3 * 2^{ 3} + 1)^2 nodes = (2^{ 6} + 1)^2 nodes =           4'225 nodes
  Level  4 ==> (2^3 * 2^{ 4} + 1)^2 nodes = (2^{ 7} + 1)^2 nodes =          16'641 nodes
  Level  5 ==> (2^3 * 2^{ 5} + 1)^2 nodes = (2^{ 8} + 1)^2 nodes =          66'049 nodes
  Level  6 ==> (2^3 * 2^{ 6} + 1)^2 nodes = (2^{ 9} + 1)^2 nodes =         263'169 nodes
  Level  7 ==> (2^3 * 2^{ 7} + 1)^2 nodes = (2^{10} + 1)^2 nodes =       1'050'625 nodes - das war bei Klawonn & Rheinbach (2010) die Problemgroesse fuer 16 PE ...
  Level  8 ==> (2^3 * 2^{ 8} + 1)^2 nodes = (2^{11} + 1)^2 nodes =       4'198'401 nodes
  Level  9 ==> (2^3 * 2^{ 9} + 1)^2 nodes = (2^{12} + 1)^2 nodes =      16'785'409 nodes
  Level 10 ==> (2^3 * 2^{10} + 1)^2 nodes = (2^{13} + 1)^2 nodes =      67'125'249 nodes
  Level 11 ==> (2^3 * 2^{11} + 1)^2 nodes = (2^{14} + 1)^2 nodes =     268'468'225 nodes
  Level 12 ==> (2^3 * 2^{12} + 1)^2 nodes = (2^{15} + 1)^2 nodes =   1'073'807'361 nodes
  Level 13 ==> (2^3 * 2^{13} + 1)^2 nodes = (2^{16} + 1)^2 nodes =   4'295'098'369 nodes
  Level 14 ==> (2^3 * 2^{14} + 1)^2 nodes = (2^{17} + 1)^2 nodes =  17'180'131'329 nodes
  Level 15 ==> (2^3 * 2^{15} + 1)^2 nodes = (2^{18} + 1)^2 nodes =  68'720'001'025 nodes


# -numPPSD 4 with FAMG (once again - next try 05032012):
UGARGS="-ex ../scripts/tests/modular_scalability_test.lua -dim 2 -grid ../data/grids/unit_square_01/unit_square_01_quads_8x8.ugx -lsMaxIter 100 -numPreRefs 3 -lsType feti -numPPSD 4"

salloc -n 64 ./ugshell $UGARGS -numRefs  7 -ds famg -ns famg
salloc -n 64 ./ugshell $UGARGS -numRefs  8 -ds famg -ns famg
salloc -n 64 ./ugshell $UGARGS -numRefs  9 -ds famg -ns famg
salloc -n 64 ./ugshell $UGARGS -numRefs 10 -ds famg -ns famg

]]


----------------------------------------------------------
-- auxiliary functions for FAMG
-- Testvectors for FAMG ---
----------------------------------------------------------
function CreateAMGTestvector(gridfunction, luaCallbackName, dim)
	local amgTestvector;
	print("          Create writer for testvector via grid function for FAMG ...")
	amgTestvector = GridFunctionVectorWriter()
	amgTestvector:set_reference_grid_function(gridfunction)
	amgTestvector:set_user_data(LuaUserNumber(luaCallbackName))
	return amgTestvector	
end

function CreateAMGTestvectorDirichlet0(dirichletBND, approxSpace)
	local amgDirichlet0
	print("          Create vector writer for 'Dirichlet 0, constant 1 else' testvector for FAMG ...")
	amgDirichlet0 = GridFunctionVectorWriterDirichlet0()
	amgDirichlet0:init(dirichletBND, approxSpace)
	return amgDirichlet0
end

----------------------------------------------------------
-- function 'SetupFETISolver()' (first parameters only for (re)naming of logfile):
----------------------------------------------------------
function SetupFETISolver(str_problem,
			 dim,
			 linMaxIterations,
			 numProcs,
			 dirichletBND, approxSpace, -- for testvector writer for FAMG (created by 'CreateAMGTestvectorDirichlet0()')
			 activateDbgWriter,
			 verbosity, logfileName)

	print("    'setup_fetisolver.lua': Setting up FETI solver ...")

	--------------------------------------------------------------------------------
	-- preconditioners, convergence checks, sub solvers and FETI solver:
	--------------------------------------------------------------------------------

	-- coarse problem solver stuff:
        local cpJac = Jacobi()
        local cpGS  = GaussSeidel()
        local cpILU = ILU()

        local cpConvCheck
        local coarseproblemSolver

	-- Neumann problem solver stuff:
        local npJac = Jacobi()
        local npGS  = GaussSeidel()
        local npILU = ILU()

        local neumannConvCheck
        local neumannSolver

	-- Dirichlet problem solver stuff:
        local dpJac = Jacobi()
        local dpGS  = GaussSeidel()
        local dpILU = ILU()

        local dirichletConvCheck
        local dirichletSolver

	-- Stuff for the FETI-DP solver itself:
        local fetiConvCheck
	local linAbsLimit  = 1e-16 -- before 25102011: 1e-7
	local linReduction = 1e-7  -- acc. to Klawonn&Rheinbach, 2010; before 25102011: 1e-16

        local fetidbgWriter

        local fetiSolver
        
     -- options for AMG
        local bWriteMat = util.HasParamOption("-AMGwriteMat")
        local bAggressiveCoarsening = false


	--------------------------------------------------------------------------------
	-- Checking for FETI solver related parameters (begin)
	-- (preconditioners are fix - edit script!)
	--------------------------------------------------------------------------------
	-- get number of processes per subdomain:
	local numProcsPerSubdomain       = util.GetParamNumber("-numPPSD", 1)

	-- types of sub solvers:
	local coarseProblemSolverType    = util.GetParam("-cps", "exact") -- choose one in ["exact" | "cg" | "hlib" ]
	local neumannProblemSolverType   = util.GetParam("-ns",    "cg") -- choose one in ["exact" | "ls" | "cg" | "bicg" | "rsamg" | "famg" ]
	local dirichletProblemSolverType = util.GetParam("-ds",    "cg") -- choose one in ["exact" | "ls" | "cg" | "bicg" | "rsamg" | "famg" ]


	-- Display parameters (or defaults):
	print("    FETI solver related parameters chosen (or defaults):")
	print("        numPPSD (numProcsPerSubdomain)   = " .. numProcsPerSubdomain)
	
	print("        cps (coarseProblemSolverType)    = " .. coarseProblemSolverType)
	print("        ns  (neumannProblemSolverType)   = " .. neumannProblemSolverType)
	print("        ds  (dirichletProblemSolverType) = " .. dirichletProblemSolverType)
	
	--------------------------------------------------------------------------------
	-- Checking for FETI solver related parameters (end)
	--------------------------------------------------------------------------------
	
	--------------------------------------------------------------------------------
	--------------------------------------------------------------------------------
	-- Gather info for domain decomposition
	--------------------------------------------------------------------------------
	--------------------------------------------------------------------------------

	-- check number of processes
	if numProcs < 2 then
		print("WARNING: number of processes is smaller than 2 - huh??")
	end
	
	if not util.IsPowerOfTwo(numProcsPerSubdomain) then
		print("WARNING: numPPSD = '" .. numProcsPerSubdomain .. "' is not a power of 2!" )
	--	exit()
	end
	
	print("    Check if numPPSD = '" .. numProcsPerSubdomain .. "' process(es) per subdomain makes sense ..." )
	
	-- compute number of subdomains
	numSubdomains = numProcs / numProcsPerSubdomain
	
	-- check that numSubdomains is greater 1 && \in \N && a power of 2.
	if numSubdomains < 2 then
		print("ERROR:   numSubdomains = numProcs / numProcsPerSubdomain = '" .. numSubdomains .. "' is smaller than 2!? Aborting!" )
		exit()
	end

	if not util.IsNaturalNumber(numSubdomains) then
		print("ERROR:   numSubdomains = numProcs / numProcsPerSubdomain = '" .. numSubdomains .. "' is NOT a natural number!? Aborting!" )
		exit()
	end
	
	if not util.IsPowerOfTwo(numSubdomains) then
		print("WARNING: numSubdomains = numProcs / numProcsPerSubdomain = '" .. numSubdomains .. "' is not a power of 2! Continuing ..." )
	-- TODO: Maybe switch to a default value then
	--	exit() -- in this case the partition can be quite erratic (at least on small (triangular) grids)..
	end
	
	print("    NumProcs is " .. numProcs .. ", NumSubDomains is " .. numSubdomains )
	--------------------------------------------------------------------------------
	
	-- create subdomain info
	print("    Create domainDecompInfo")
	domainDecompInfo = StandardDomainDecompositionInfo()
	domainDecompInfo:set_num_subdomains(numSubdomains)

	-- test one to many interface creation
	if verbosity >= 1 then
		for i=0,numProcs-1 do
			print("    subdom of proc " .. i .. ": " .. domainDecompInfo:map_proc_id_to_subdomain_id(i))
		end
	end
	--------------------------------------------------------------------------------
	--------------------------------------------------------------------------------
	-- Setup of FETI solver:
	--------------------------------------------------------------------------------
	--------------------------------------------------------------------------------
	
	--if verbosity >= 1 then
		print("")
		print("    Setup of FETI coarse and sub problem solvers ...")
	--end
	----------------------------------------------------------
	-- preconditioners used by FETI sub solvers
	----------------------------------------------------------
	cpJac = Jacobi()
	cpJac:set_damp(0.8)
	cpGS  = GaussSeidel()
	cpILU = ILU()
	
	npJac = Jacobi()
	npJac:set_damp(0.8)
	npGS  = GaussSeidel()
	npILU = ILU()
	
	dpJac = Jacobi()
	dpJac:set_damp(0.8)
	dpGS  = GaussSeidel()
	dpILU = ILU()
	
	----------------------------------------------------------
	-- create and configure coarse problem solver
	----------------------------------------------------------
	--if verbosity >= 1 then
		print("       Create and configure coarse problem solver of type '" .. coarseProblemSolverType .. "' ...")
	--end
	-- choose solver for coarse problem
	if coarseProblemSolverType == "exact" then
	
		coarseproblemSolver = LU()
	
	elseif coarseProblemSolverType == "cg" then
	
		coarseproblemSolver = CG()
		coarseproblemSolver:set_preconditioner(cpILU) -- Absturz: "Signal: Bus error (10)", "Signal code:  (5583)", etc ...
		--coarseproblemSolver:set_preconditioner(cpJac) -- "ERROR in 'JacobiPreconditioner::apply': Cannot change parallel status of correction to consistent."
		--coarseproblemSolver:set_preconditioner(cpGS) -- "ERROR in 'Gauss-Seidel::apply': Cannot change parallel storage type of correction to consistent."
		-- Ganz ohne preconditioner: "Cannot convert z to consistent vector."
	
	elseif coarseProblemSolverType == "hlib" then
	
		coarseproblemSolver = HLIB() -- create HLIB Solver
		coarseproblemSolver:set_hlib_accuracy_H(1.e-4)
	
	else
		print ("ERROR: coarse problem solver not specified ==> exit")
		exit()
	end
	
	-- define convergence criteria for the coarse problem solver
	cpConvCheck = StandardConvergenceCheck()
	cpConvCheck:set_maximum_steps(20)
	cpConvCheck:set_minimum_defect(1e-10)
	cpConvCheck:set_reduction(1e-16)
	cpConvCheck:set_verbose_level(false)
	
	coarseproblemSolver:set_convergence_check(cpConvCheck)
	
	----------------------------------------------------------
	-- create and configure sub problem solvers
	----------------------------------------------------------

	----------------------------------------------------------
	-- create and configure Neumann problem solver
	----------------------------------------------------------
	--if verbosity >= 1 then
		print("       Create and configure Neumann problem solver of type '" .. neumannProblemSolverType .. "' ...")
	--end
	-- choose solver for Neumann problems
	if neumannProblemSolverType == "exact" then
	
		neumannSolver = LU()
	
	elseif neumannProblemSolverType == "ls" then
	
		neumannSolver = LinearSolver()
		neumannSolver:set_preconditioner(npILU) -- npJac
	
	elseif neumannProblemSolverType == "cg" then
	
		neumannSolver = CG()
		neumannSolver:set_preconditioner(npILU) -- npJac
	
	elseif neumannProblemSolverType == "bicg" then
	
		neumannSolver = BiCGStab() -- ERROR in 'PrimalSubassembledMatrixInverse::init': Could not solve local problem to compute Schur complement w.r.t. primal unknowns. - war das bei 1x1?
		-- BiCGStab statt CG hatte z.T. eine Verdopplung bis Verdreifachung der Iterationszahlen
		-- (apply_F), z.T. eine Verringerung von einem Viertel bis zur Haelfte (backsolve; compute_d) bewirkt!?
		neumannSolver:set_preconditioner(npILU) -- npJac
	
	elseif neumannProblemSolverType == "rsamg" or neumannProblemSolverType == "famg" then
--		print ("TMP: RSAMG or FAMG as Neumann problem solver!")

		if neumannProblemSolverType == "famg" then
			--print ("    Create FAMG as Neumann problem solver ... ")

			npAMG = FAMGPreconditioner()	
			npAMG:set_delta(0.5)
			npAMG:set_theta(0.95)
			npAMG:set_aggressive_coarsening(bAggressiveCoarsening)

--[[
---- {
			-- add testvector which is 1 everywhere and only 0 on the dirichlet Boundary.
			testvectorwriter = CreateAMGTestvectorDirichlet0(dirichletBND, approxSpace)
print("          (1)")
			testvector = GridFunction(approxSpace)
print("          (2)")
			testvectorwriter:update(testvector)	
print("          'Dirichlet 0, constant 1 else' testvector for FAMG created (TMP)!")
			npAMG:add_testvector(testvectorwriter, 1.0)
---- }
   ]]				
			local FAMGtestvectorSmoother = Jacobi()
			FAMGtestvectorSmoother:set_damp(0.66)
				
			npAMG:set_testvector_smooths(1) -- was: 'set_testvector_damps(1)' (17022012ih)
			npAMG:set_damping_for_smoother_in_interpolation_calculation(0.66)
			npAMG:set_testvector_smoother(FAMGtestvectorSmoother) -- was: 'set_testvectorsmoother_damps(1)' (17022012ih)
			npAMG:set_testvector_from_matrix_rows(true)
				
			if bExternalCoarsening then
				npAMG:set_external_coarsening(true)
				npAMG:set_parallel_coarsening(GetColorCoarsening())
				-- npAMG:set_parallel_coarsening(GetFullSubdomainBlockingCoarsening())
				-- npAMG:set_parallel_coarsening(GetRS3Coarsening())
				if dim == 2 then
					npAMG:set_strong_connection_external(0.2)
				else
					npAMG:set_strong_connection_external(0.1)
				end				
			end
			
			if bWriteMat then
				npAMG:write_testvectors(true)
			end
			npAMG:set_galerkin_truncation(1e-9)
			npAMG:set_H_reduce_interpolation_nodes_parameter(0.1)
			npAMG:set_prereduce_A_parameter(0.01)
		else
			print ("    Create RSMG as Neumann problem solver ... ")
			npAMG = RSAMGPreconditioner()
			-- amg:set_parallel_coarsening(GetFullSubdomainBlockingCoarsening())
			-- amg:set_parallel_coarsening(GetColorCoarsening()) --
			npAMG:set_parallel_coarsening(GetRS3Coarsening()) --
			-- amg:set_parallel_coarsening(GetSimpleParallelCoarsening())
			if bAggressiveCoarsening then
				npAMG:enable_aggressive_coarsening_A(2)
			end
		end
		
		maxBase = util.GetParamNumber("-maxBase", 1000) -- TODO (maybe): make 'maxBase' different for "np solver" and "dp solver"!
		
		if bWriteMat then
			vectorWriter = GridFunctionPositionProvider()
			vectorWriter:set_reference_grid_function(u)
			npAMG:set_position_provider(vectorWriter)
			npAMG:set_matrix_write_path("/Users/mrupp/matrices/")
		end
		
		
		npAMGGS = GaussSeidel()
		npAMGBase = LU()
		
		npAMG:set_num_presmooth(3)
		npAMG:set_num_postsmooth(3)
		npAMG:set_cycle_type(1)
		npAMG:set_presmoother(npAMGGS)
		npAMG:set_postsmoother(npAMGGS)
		npAMG:set_base_solver(npAMGBase)
		npAMG:set_max_levels(20)
		npAMG:set_max_nodes_for_base(maxBase)
		npAMG:set_max_fill_before_base(0.7)
		npAMG:set_fsmoothing(true)
		
		print("\nThe AMG solver:")
		npAMG:tostring()
		print("")
		
		neumannSolver = LinearSolver()
		neumannSolver:set_preconditioner(npAMG)
	
	else
		print ("ERROR: Neumann problem solver not specified ==> exit")
		exit()
	end
	
	-- define convergence criteria for the Neumann problem solver
	neumannConvCheck = StandardConvergenceCheck()
	neumannConvCheck:set_maximum_steps(2000)
	neumannConvCheck:set_minimum_defect(1e-10)
	neumannConvCheck:set_reduction(1e-16)
	neumannConvCheck:set_verbose_level(false)
	
	neumannSolver:set_convergence_check(neumannConvCheck)
	
	----------------------------------------------------------
	-- create and configure Dirichlet problem solver
	----------------------------------------------------------
	--if verbosity >= 1 then
		print("       Create and configure Dirichlet problem solver of type '" .. dirichletProblemSolverType .. "' ...")
	--end
	-- choose solver for Dirichlet problems
	if dirichletProblemSolverType == "exact" then
	
		dirichletSolver = LU()
	
	elseif dirichletProblemSolverType == "ls" then
	
		dirichletSolver = LinearSolver()
		dirichletSolver:set_preconditioner(dpILU) -- dpJac
	
	elseif dirichletProblemSolverType == "cg" then
	
		dirichletSolver = CG()
		dirichletSolver:set_preconditioner(dpILU) -- dpJac
	
	elseif dirichletProblemSolverType == "bicg" then
	
		dirichletSolver = BiCGStab() -- not yet tested for Dirichlet problem!
		dirichletSolver:set_preconditioner(dpILU) -- dpJac
	
	elseif dirichletProblemSolverType == "rsamg" or dirichletProblemSolverType == "famg" then
--		print ("TMP: RSAMG or FAMG as Dirichlet problem solver!")

		if dirichletProblemSolverType == "famg" then
			--print ("create FAMG as Dirichlet problem solver ... ")

			dpAMG = FAMGPreconditioner()	
			dpAMG:set_delta(0.5)
			dpAMG:set_theta(0.95)
			dpAMG:set_aggressive_coarsening(bAggressiveCoarsening)
				
--[[
---- {
			-- add testvector which is 1 everywhere and only 0 on the dirichlet Boundary.
			testvectorwriter = CreateAMGTestvectorDirichlet0(dirichletBND, approxSpace)
print("          (1)")
			testvector = GridFunction(approxSpace)
print("          (2)")
			testvectorwriter:update(testvector)	
print("          'Dirichlet 0, constant 1 else' testvector for FAMG created (TMP)!")
			dpAMG:add_testvector(testvectorwriter, 1.0)
---- }
]]
			if bWriteMat then
				dpAMG:write_testvectors(true)
			end
			local FAMGtestvectorSmoother = Jacobi()
			FAMGtestvectorSmoother:set_damp(0.66)
				
			dpAMG:set_testvector_smooths(1) -- was: 'set_testvector_damps(1)' (17022012ih)
			dpAMG:set_damping_for_smoother_in_interpolation_calculation(0.66)
			dpAMG:set_testvector_smoother(FAMGtestvectorSmoother) -- was: 'set_testvectorsmoother_damps(1)' (17022012ih)
			dpAMG:set_testvector_from_matrix_rows(true)
				
			if bExternalCoarsening then
				dpAMG:set_external_coarsening(true)
				dpAMG:set_parallel_coarsening(GetColorCoarsening())
				-- dpAMG:set_parallel_coarsening(GetFullSubdomainBlockingCoarsening())
				-- dpAMG:set_parallel_coarsening(GetRS3Coarsening())
				if dim == 2 then
					dpAMG:set_strong_connection_external(0.2)
				else
					dpAMG:set_strong_connection_external(0.1)
				end				
			end
			
			dpAMG:set_galerkin_truncation(1e-9)
			dpAMG:set_H_reduce_interpolation_nodes_parameter(0.1)
			dpAMG:set_prereduce_A_parameter(0.01)
		else
			print ("    Create RSMG as Dirichlet problem solver ... ")
			dpAMG = RSAMGPreconditioner()
			-- amg:set_parallel_coarsening(GetFullSubdomainBlockingCoarsening())
			-- amg:set_parallel_coarsening(GetColorCoarsening()) --
			dpAMG:set_parallel_coarsening(GetRS3Coarsening()) --
			-- amg:set_parallel_coarsening(GetSimpleParallelCoarsening())
			if bAggressiveCoarsening then
				dpAMG:enable_aggressive_coarsening_A(2)
			end
		end
		
		maxBase = util.GetParamNumber("-maxBase", 1000) -- TODO (maybe): make 'maxBase' different for "np solver" and "dp solver"!
		
		if bWriteMat then
			vectorWriter = GridFunctionPositionProvider()
			vectorWriter:set_reference_grid_function(u)
			dpAMG:set_position_provider(vectorWriter)
			dpAMG:set_matrix_write_path("/Users/mrupp/matricesNew/")
		end
		
		dpAMGGS = GaussSeidel()
		dpAMGBase = LU()
		
		dpAMG:set_num_presmooth(3)
		dpAMG:set_num_postsmooth(3)
		dpAMG:set_cycle_type(1)
		dpAMG:set_presmoother(dpAMGGS)
		dpAMG:set_postsmoother(dpAMGGS)
		dpAMG:set_base_solver(dpAMGBase)
		dpAMG:set_max_levels(20)
		dpAMG:set_max_nodes_for_base(maxBase)
		dpAMG:set_max_fill_before_base(0.7)
		dpAMG:set_fsmoothing(true)
		
		print("\nThe AMG solver:")
		dpAMG:tostring()
		print("")
		
		dirichletSolver = LinearSolver()
		dirichletSolver:set_preconditioner(dpAMG)
	
	else
		print ("ERROR: Dirichlet problem solver not specified ==> exit")
		exit()
	end
	
	-- define convergence criteria for the Dirichlet problem solver
	dirichletConvCheck = StandardConvergenceCheck()
	dirichletConvCheck:set_maximum_steps(2000)
	dirichletConvCheck:set_minimum_defect(1e-10)
	dirichletConvCheck:set_reduction(1e-16)
	dirichletConvCheck:set_verbose_level(false)
	
	dirichletSolver:set_convergence_check(dirichletConvCheck)
	
	----------------------------------------------------------
	-- create and configure FETI Solver
	----------------------------------------------------------
	fetiSolver = FETI()

	fetiSolver:set_domain_decomp_info(domainDecompInfo)
	
	fetiSolver:set_neumann_solver(neumannSolver)
	fetiSolver:set_dirichlet_solver(dirichletSolver)
	fetiSolver:set_coarse_problem_solver(coarseproblemSolver)
	
	-- define convergence criteria for the FETI solver
	fetiConvCheck = StandardConvergenceCheck()

	print("    'setup_fetisolver.lua': linMaxIterations = " .. linMaxIterations .. ", linAbsLimit = " .. linAbsLimit .. ", linReduction = " .. linReduction)

	fetiConvCheck:set_maximum_steps(linMaxIterations)
	fetiConvCheck:set_minimum_defect(linAbsLimit)
	fetiConvCheck:set_reduction(linReduction)
	
	fetiSolver:set_convergence_check(fetiConvCheck)

	if activateDbgWriter >= 1 then
		print("    Setting debug writer for 'fetiSolver' (no make consistent (former: 'raw data')")
		-- debug writer
		fetidbgWriter = GridFunctionDebugWriter(approxSpace)
		fetidbgWriter:set_vtk_output(true)
		-- alt:fetidbgWriter:set_print_raw_data(true) -- if 'true' print "raw" data (no "make consistent" before printing): for checking temporary results in FETI
		fetidbgWriter:set_print_consistent(false) -- if 'false' print "raw" data (no "make consistent" before printing): for checking temporary results in FETI

		fetiSolver:set_debug(fetidbgWriter)
	end
	
	fetiSolver:set_test_one_to_many_layouts(true)

--[[ Layout tests sind irgendwann rausgeflogen ...
	if verbosity >= 2 then
		print("    Performing layout tests:")
		--BuildDomainDecompositionLayoutsTest2d(u, domainDecompInfo);
		--OneToManyTests2d(u)
	end
]]

	----------------------------------------------------------
-- add feti specific name parts to new logfile name
	str_spsolvers = "spss" -- "sub problem solvers"
	str_nppsd     = "nppsd" -- "number of processes per Subdomain"

	str_nppsd     = str_nppsd     .. "-" .. numProcsPerSubdomain
	logfileName   = logfileName   .. "_" .. str_nppsd

	str_spsolvers = str_spsolvers .. "-" .. dirichletProblemSolverType
                                      .. "-" .. neumannProblemSolverType
                                      .. "-" .. coarseProblemSolverType
	logfileName   = logfileName   .. "_" .. str_spsolvers

	print("    'setup_fetisolver.lua': returning FETI solver 'fetiSolver', ready for application!")
	return fetiSolver, logfileName
end
