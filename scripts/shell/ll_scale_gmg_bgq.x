# TO CHECK: Job-Name and comment:
# @ job_name = ug4rev7229_laplace_2d_gmg_numprerefs_3
# @ comment = "BGQ GMG Scalability Laplace 2D"
# @ error = $(job_name).$(jobid).out
# @ output = $(job_name).$(jobid).out
# @ environment = COPY_ALL;

########################################
# TO CHECK: Minimal "wall clock time"
# @ wall_clock_limit = 00:15:00
########################################

# @ notification = error
# @ notify_user = ingo.heppner@gcsc.uni-frankfurt.de
# @ job_type = bluegene

########################################
# Connection type in [TORUS| MESH | PREFER_TORUS] - NOTE: New keyword 'bg_connectivity' on JuQueen instead of 'bg_connection'!
# Possible values according to 'http://www.fz-juelich.de/SharedDocs/Downloads/IAS/JSC/EN/JUQUEEN/BGQJSCAccessDocter.pdf':
# bg_connectivity = MESH | TORUS | Either | <Xa Xb Xc Xd>
#                                            Xa,Xb,Xc,Xd - Torus or Mesh for each dimension
# @ bg_connectivity = TORUS

# TO CHECK: - Achtung: Highmessage vom 18.10.2012:
# * 18.10.2012         16 Blue Gene/Q racks are available for production now
# *              The smallest job size is 128 nodes until the IO rich rack will be
# *	              reinstalled in the 24 rack system.
# @ bg_size = 128
# 'http://www.fz-juelich.de/ias/jsc/EN/Expertise/Supercomputers/JUQUEEN/UserInfo/LoadLeveler.html':
# The 'bg_size' keyword specifies the number of compute nodes the job should use.
# Blue Gene/Q only allows blocks including 32, 64, 128, 256 and multiples of 512 compute nodes.
# Thus e.g. bg_size of 1 specifies a block of size 32 and bg_size of 129 requests a partition of size 256.
#
# A bg_size of 32 comprises  512 PE,
# a bg_size of 64 comprises 1024 PE, ...!

# @ queue

################################################################################
# Output of runtime environment variables beginning with 'LOADL_' (for example
# 'LOADL_BG_CONNECTION', to check if connection type is really "TORUS"; see e.g.
# http://bluegene.epfl.ch/LoadLeveler_UsingAdministrating_a2278818.pdf, p. 420):
env | grep LOADL_

################################################################################
################################################################################
### My jobs so far ...
### 'runjob' calls commented out just after submitting via 'll_submit <ll-file>',
### to get a history of all our jobs; job id numbers added to the job comments
### are given by the 'll_submit' command.
###
### Raw timing analysis of these jobs:
### egrep 'NumProcs|main   |unknown|MPI_Init| ASS_AssembleLinearOperatorRhsAndSolution| ALS_InitLinearSolver| ALS_ApplyLinearSolver' ug4_fullystatic_laplace*
###
### Check start and end norms of residuum, if solver converged and number of iterations:
### grep -A3 " Iterative Linear Solver "  ug4_fullystatic_laplace*
### grep -B1 " Iteration not successful " ug4_fullystatic_laplace*
### grep -B1 " reached after "            ug4_fullystatic_laplace*
### egrep 'Iterative Linear Solver | reached after ' ug4_fullystatic_laplace*
### or, shorter
### grep "ANALYZER INFO: linear solver converged in" ug4_fullystatic_laplace*
################################################################################
# Syntax for 'runjob' (from  "Message-of-today"):
# runjob [--np <n>] [--ranks-per-node <r>] --exe <excutable> [--args <args>]
#
# See also:
# runjob --help, man runjob.
# Most important parameters (if not self-explaining):
# --np <n>            : positive number of ranks in the entire job
# --ranks-per-node <r>: number of ranks per node: 1, 2, 4, 8, 16, 32, or 64 (default: 1)
# --mapping <arg>     : permutation of 'ABCDET' or path to mapping file (default: ABCDET) - T increments first, then E (which can get maximal 2 on BG/Q)
################################################################################

################################################################################
# Serie 1:   4 bis  16k PE with 'numPreRefs 3'           ==> TO CHECK: Job Name!
################################################################################
################################################################################

##########################################
# Set ug executable, arguments (to shorten the lengthy commands):
##########################################
#UGSHELL=ugshell

# 2-D arguments:
UGARGS2D="-ex ../apps/scaling_tests/modular_scalability_test.lua -dim 2 -grid ../data/grids/unit_square_01/unit_square_01_quads_8x8.ugx -lsIterator gmg -lsMaxIter 100 -verb 0 -numPreRefs 3"
# 3-D arguments:
UGARGS3D="-ex ../apps/scaling_tests/modular_scalability_test.lua -dim 3 -grid ../data/grids/unit_square_01/unit_cube_01_hex_2x2x2.ugx   -lsIterator gmg -lsMaxIter 100 -verb 0 -numPreRefs 2"

echo "Content of UGARGS2D is: '$UGARGS2D'!"
echo "Content of UGARGS3D is: '$UGARGS3D'!"

################################################################################
# 17102012 (alle Schritte bis hierher, insbesondere die Ergaenzung von
# '-fno-strict-aliasing' in 'ug_cmake_includes_included.txt' zwecks Behebung
# des Abstuerzens in 'll_scale_gmg_bgq_first-try.x'):
################################################################################

##########################################
#### Until  256 PE 'bg_size =  32'!   ####
##########################################
# job "juqueen2c1.zam.kfa-juelich.de.28617":
#runjob --np      4 --ranks-per-node  1 --mapping ABCDET --verbose 3 : ./ugshell $UGARGS2D -numRefs  6

# job "juqueen1c1.zam.kfa-juelich.de.28618":
#runjob --np     16 --ranks-per-node 16 --mapping ABCDET --verbose 3 : ./ugshell $UGARGS2D -numRefs  7

# job "juqueen1c1.zam.kfa-juelich.de.28619":
#runjob --np     64 --ranks-per-node 16 --mapping ABCDET --verbose 3 : ./ugshell $UGARGS2D -numRefs  8

# job "juqueen1c1.zam.kfa-juelich.de.28620":
#runjob --np    256 --ranks-per-node 16 --mapping ABCDET --verbose 3 : ./ugshell $UGARGS2D -numRefs  9

##########################################
####  1024 PE: 'bg_size =   64'!      ####
##########################################
# Hierarchische Nachverteilung: Letzter Schritt von 64 auf 1k, erster von 1 auf 64:
# job "juqueen1c1.zam.kfa-juelich.de.28621":
#runjob --np   1024 --ranks-per-node 16 --mapping ABCDET --verbose 3 : ./ugshell $UGARGS2D -numRefs 10 -hRedistFirstLevel 5 -hRedistStepSize 100 -hRedistNewProcsPerStep  16

##########################################
####  4096 PE: 'bg_size =  256'!      ####
##########################################
# Hierarchische Nachverteilung: Letzter Schritt von 64 auf 4k, erster von 1 auf 64:
# job "juqueen1c1.zam.kfa-juelich.de.28622":
#runjob --np   4096 --ranks-per-node 16 --mapping ABCDET --verbose 3 : ./ugshell $UGARGS2D -numRefs 11 -hRedistFirstLevel 5 -hRedistStepSize 100 -hRedistNewProcsPerStep  64

##########################################
#### 16384 PE: 'bg_size = 1024'!      ####
##########################################
# Hierarchische Nachverteilung: Letzter Schritt von 64 auf 16k, erster von 1 auf 64:
# job "juqueen1c1.zam.kfa-juelich.de.28623":
#runjob --np  16384 --ranks-per-node 16 --mapping ABCDET --verbose 3 : ./ugshell $UGARGS2D -numRefs 12 -hRedistFirstLevel 5 -hRedistStepSize 100 -hRedistNewProcsPerStep 256

##########################################
####  65536 PE: 'bg_size = 4096'!      ####
##########################################
# job "juqueen1c1.zam.kfa-juelich.de.28624":
# Hierarchische Nachverteilung: Letzter Schritt von 256 auf 64k, erster von 1 auf 256:
#runjob --np  65536 --ranks-per-node 16 --mapping ABCDET --verbose 3 : ./ugshell $UGARGS2D -numRefs 13 -hRedistFirstLevel 6 -hRedistStepSize 100 -hRedistNewProcsPerStep 256
