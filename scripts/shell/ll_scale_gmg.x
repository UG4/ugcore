#!/bin/bash
################################################################################
# Begin:   06022012ih
# Purpose: Loadleveler command file for scalability test of GMG iterator.
# Usage:   Submit it (after careful editing to specify your job ...) via
#          'llsubmit ll_scale_gmg.x' in your 'bin/' directory.
#
# Please do not modify this file for submitting your jobs - work on your own copy!
#
################################################################################

################################################################################
# I.  For further infos about IBMs LoadLeveler syntax see 'll_template.x'!
################################################################################

################################################################################
# II. What we are interested in:
#
# Weak scaling measurements, starting from 4 PE, coarse grid 8 x 8 quadrilaterals in 2-D
#     4 PE: Level  7 ==> (2^3 * 2^{ 7} + 1)^2 nodes = (2^{10} + 1)^2 nodes =       1'050'625 nodes, ~262'656 nodes/PE
#    16 PE: Level  8 ==> (2^3 * 2^{ 8} + 1)^2 nodes = (2^{11} + 1)^2 nodes =       4'198'401 nodes, ~262'400 nodes/PE
#    64 PE: Level  9 ==> (2^3 * 2^{ 9} + 1)^2 nodes = (2^{12} + 1)^2 nodes =      16'785'409 nodes, ~262'272 nodes/PE
#   256 PE: Level 10 ==> (2^3 * 2^{10} + 1)^2 nodes = (2^{13} + 1)^2 nodes =      67'125'249 nodes, ~262'208 nodes/PE
#  1024 PE: Level 11 ==> (2^3 * 2^{11} + 1)^2 nodes = (2^{14} + 1)^2 nodes =     268'468'225 nodes, ~262'176 nodes/PE
#  4096 PE: Level 12 ==> (2^3 * 2^{12} + 1)^2 nodes = (2^{15} + 1)^2 nodes =   1'073'807'361 nodes, ~262'160 nodes/PE
#  16 k PE: Level 13 ==> (2^3 * 2^{13} + 1)^2 nodes = (2^{16} + 1)^2 nodes =   4'295'098'369 nodes, ~262'152 nodes/PE
#  64 k PE: Level 14 ==> (2^3 * 2^{14} + 1)^2 nodes = (2^{17} + 1)^2 nodes =  17'180'131'329 nodes, ~262'148 nodes/PE
# 256 k PE: Level 15 ==> (2^3 * 2^{15} + 1)^2 nodes = (2^{18} + 1)^2 nodes =  68'720'001'025 nodes, ~262'146 nodes/PE
# etc. (nur nicht mehr auf JuGene ...)
#
# It turned out that we have to refine one level less to avoid memory problems:
#     4 PE: Level  6 ==> (2^3 * 2^{ 6} + 1)^2 nodes = (2^{ 9} + 1)^2 nodes =         263'169 nodes, ~ 65'792 nodes/PE
#    16 PE: Level  7 etc.
#
# Maybe we are able to extend the numbers of refinements with the improvements in
# memory consumption introduced in January 2012 ...
#
################################################################################
################################################################################

################################################################################
# III. Job definitions (note that keywords must not be followed by comments ==> put them on separate lines):
################################################################################
# TO CHECK: Job-Name and comment:
# @ job_name = ug4_fullystatic_laplace_3d
# @ comment = "BGP GMG test 3d"
# @ error = $(job_name).$(jobid).out
# @ output = $(job_name).$(jobid).out
# @ environment = COPY_ALL;
########################################
# TO CHECK: Minimal 'wall_clock_limit' depends obviously on job size:
#           Up to 'bg_size =  16384' (nec. for  64k PE) 30 min are enough,
#           for   'bg_size =  65536' (nec. for 256k PE) minimum is 45 min (acc. to error message of LoadLeveler which refuses to submit the job with only 30 min wall_clock_limit)
# @ wall_clock_limit = 00:30:00
########################################
# @ notification = error
# @ notify_user = ingo.heppner@gcsc.uni-frankfurt.de
# @ job_type = bluegene

########################################
# Connection type in [TORUS| MESH | PREFER_TORUS] - default: 'bg_connection=MESH'.
# For 'bg_connection=TORUS' 'bg_size' must be >= 512 (i.e. one midplane)
# http://www.prace-ri.eu/IMG/pdf/Best-practise-guide-JUGENE-v0-3.pdf states on p. 8:
# A midplane is the smallest building block for the TORUS network.
# And on p. 9:
# "A three- dimensional torus topology is available to the user application only if
# the assigned partition of the machine is a multiple of the [a] midplane (16 node cards, half of a rack)."
# @ bg_connection = TORUS
########################################
# IMPORTANT: Specification of 'bg_size'! Statement given by JuGene-Support-Team:
#    "Durch das Batch-Script wird die in bg_size angegebene Partitionsgroesse
#    angefordert. Diese muss dann daher auch vom Scheduler reserviert werden,
#    bestimmt die Klasse des Jobs und wird auch anschliessend abgerechnet".
#
#    Best-practise-guide-JUGENE-v0-3.pdf, p. 32:
#    'bg_size' is the size of the Blue Gene job in *number of compute nodes* to be used.
#    'bg_size = 32' (= one node board) is the minimum on JUGENE [sc.: if 'connection_type' != 'TORUS'].
#
# Examples for 'bg_size' (for 2-D scalability test):
#    'bg_size =    512': enough for 0.5 Ki*4 =   2'048 PE, i.e. for 4, 16, 64, 256, 1024 PE
#    'bg_size =   1024': enough for   1 Ki*4 =   4'096 PE, i.e. only for   4096 PE
#    'bg_size =   4096': enough for   4 Ki*4 =  16'384 PE, i.e. only for  16384 PE
#    'bg_size =  16384': enough for  16 Ki*4 =  65'536 PE, i.e. only for  65536 PE
#    'bg_size =  65536': enough for  64 Ki*4 = 262'144 PE, i.e. only for 262144 PE.
# Other (i.e. 3-D scalability test):
#    'bg_size =   8192': enough for   8 Ki*4 =  32'768 PE,
#    'bg_size =  32768': enough for  32 Ki*4 = 131'072 PE.
#
# TO CHECK:
# @ bg_size = 8192
########################################
# @ queue

# Note: Letztlich bestimmt der 'mpirun'-Parameter '-np' (wie gehabt) die Anzahl der PE's!
# Die oben spezifizierte Groesse 'bg_size' bestimmt lediglich eine bestimmte (Mindest-)
# Partitionsgroesse (da deucht mir die Dokumentation etwas missverstaendlich/widerspruechlich
# zu sein ...)!

################################################################################
# Output of runtime environment variables beginning with 'LOADL_' (for example
# 'LOADL_BG_CONNECTION', to check if connection type is really "TORUS"; see e.g.
# http://bluegene.epfl.ch/LoadLeveler_UsingAdministrating_a2278818.pdf, p. 420):
env | grep LOADL_



################################################################################
################################################################################
### My jobs so far ...
### 'mpirun' calls commented out just after submitting via 'll_submit <ll-file>',
### to get a history of all our jobs; job id numbers added to the job comments
### are given by the 'll_submit' command (maybe this turns out to be useful one day)!
###
### Please note that this file is a "cleaned up" and slightly modified/improved
### version!
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
################################################################################

##########################################
# Set ug executable, arguments (to shorten the lengthy commands):
##########################################
UGSHELL=ugshell_rev4354_fully_static

# 2-D arguments:
UGARGS2D="-ex ../scripts/tests/modular_scalability_test.lua -dim 2 -grid ../data/grids/unit_square_01/unit_square_01_quads_8x8.ugx -lsIterator gmg -lsMaxIter 100 -verb 0 -numPreRefs 3"
# 3-D arguments:
UGARGS3D="-ex ../scripts/tests/modular_scalability_test.lua -dim 3 -grid ../data/grids/unit_square_01/unit_cube_01_hex_2x2x2.ugx -lsIterator gmg -lsMaxIter 100 -verb 0 -numPreRefs 2"

##########################################
# load path - only necessary (in my experience) for dynamically linked executables:
# 'mpirun' parameter '-env LD_LIBRARY_PATH=/bgsys/drivers/ppcfloor/comm/lib/'
##########################################

##########################################
##########################################
###             2-D                    ###
###        4 PE  -  256 Ki PE          ###
##########################################
##########################################

##########################################
#### Until 1024 PE 'bg_size = 512'!   ####
# job "jugene4b.292614":
#mpirun -np      4 -exe ./$UGSHELL -mode VN -mapfile TXYZ -verbose 2 -args "$UGARGS2D -numRefs  6"
# job "jugene4b.292615":
#mpirun -np     16 -exe ./$UGSHELL -mode VN -mapfile TXYZ -verbose 2 -args "$UGARGS2D -numRefs  7"
# job "jugene4b.292616":
#mpirun -np     64 -exe ./$UGSHELL -mode VN -mapfile TXYZ -verbose 2 -args "$UGARGS2D -numRefs  8"
# job "jugene4b.292617":
#mpirun -np    256 -exe ./$UGSHELL -mode VN -mapfile TXYZ -verbose 2 -args "$UGARGS2D -numRefs  9"
# job "jugene4b.292618":
# Hierarchische Nachverteilung: Letzter Schritt von 64 auf 1k, erster von 1 auf 64:
#mpirun -np   1024 -exe ./$UGSHELL -mode VN -mapfile TXYZ -verbose 2 -args "$UGARGS2D -numRefs 10 -hRedistFirstLevel 5 -hRedistStepSize 100 -hRedistNewProcsPerStep  16"
#### Up to this with 'bg_size = 512'! ####
##########################################

##########################################
#### For this job 'bg_size =   1024'! ####
# job "jugene4b.292619":
# Hierarchische Nachverteilung: Letzter Schritt von 64 auf 4k, erster von 1 auf 64:
#mpirun -np   4096 -exe ./$UGSHELL -mode VN -mapfile TXYZ -verbose 2 -args "$UGARGS2D -numRefs 11 -hRedistFirstLevel 5 -hRedistStepSize 100 -hRedistNewProcsPerStep  64"

#### For this job 'bg_size =   4096'! ####
# job "jugene4b.292621":
# Hierarchische Nachverteilung: Letzter Schritt von 64 auf 16k, erster von 1 auf 64:
#mpirun -np  16384 -exe ./$UGSHELL -mode VN -mapfile TXYZ -verbose 2 -args "$UGARGS2D -numRefs 12 -hRedistFirstLevel 5 -hRedistStepSize 100 -hRedistNewProcsPerStep 256"

#### For this job 'bg_size =  16384'! ####
# job "jugene4b.292622":
# Hierarchische Nachverteilung: Letzter Schritt von 256 auf 64k, erster von 1 auf 256:
#mpirun -np  65536 -exe ./$UGSHELL -mode VN -mapfile TXYZ -verbose 2 -args "$UGARGS2D -numRefs 13 -hRedistFirstLevel 6 -hRedistStepSize 100 -hRedistNewProcsPerStep 256"

##########################################
# And now ... the next big step (sadly the last step on JuGene also ...), with 256 Ki PE:
##########################################

#### For this job 'bg_size = 65536' - and wall_time 45 min!  ####
# job "jugene4b.289331.0":
#mpirun -np 262144 -exe ./$UGSHELL -mode VN -mapfile TXYZ -verbose 2 -args "$UGARGS2D -numRefs 14 -hRedistFirstLevel 6 -hRedistStepSize 3 -hRedistNewProcsPerStep 64 -hRedistMaxSteps 2"


################################################################################
################################################################################

##########################################
##########################################
###             3-D                    ###
###        1 PE  -  256 Ki PE          ###
##########################################
##########################################

#### Until  512 PE 'bg_size = 512'!   ####
# job "jugene3b.289034":
#mpirun -np      1 -exe ./$UGSHELL -mode VN -mapfile TXYZ -verbose 2 -args "$UGARGS3D -numRefs  4"

# job "jugene3b.289035":
#mpirun -np      8 -exe ./$UGSHELL -mode VN -mapfile TXYZ -verbose 2 -args "$UGARGS3D -numRefs  5"

# job "jugene3b.289036":
#mpirun -np     64 -exe ./$UGSHELL -mode VN -mapfile TXYZ -verbose 2 -args "$UGARGS3D -numRefs  6 -hRedistFirstLevel 4 -hRedistStepSize 100 -hRedistNewProcsPerStep   8"

# job "jugene3b.289037":
#mpirun -np    512 -exe ./$UGSHELL -mode VN -mapfile TXYZ -verbose 2 -args "$UGARGS3D -numRefs  7 -hRedistFirstLevel 4 -hRedistStepSize 100 -hRedistNewProcsPerStep  32"
#### Up to this with 'bg_size = 512'! ####
##########################################

##########################################
#### For this job 'bg_size =   1024'! ####
# job "jugene3b.289038":
#mpirun -np   4096 -exe ./$UGSHELL -mode VN -mapfile TXYZ -verbose 2 -args "$UGARGS3D -numRefs  8 -hRedistFirstLevel 5 -hRedistStepSize 100 -hRedistNewProcsPerStep  64"


#### For this job 'bg_size =   8192'! ####
# job "jugene4b.289733":
#mpirun -np  32768 -exe ./$UGSHELL -mode VN -mapfile TXYZ -verbose 2 -args "$UGARGS3D -numRefs 9 -hRedistFirstLevel 4 -hRedistStepSize 2 -hRedistNewProcsPerStep  64 -hRedistMaxSteps 2"


##########################################
# And now ... the next (and last) big step also in 3-D, with 256 Ki PE:
##########################################

#### For this job 'bg_size =   65536' - and wall_time 45 min! ####
# job "jugene4b.289332.0":
#mpirun -np 262144 -exe ./$UGSHELL -mode VN -mapfile TXYZ -verbose 2 -args "$UGARGS3D -numRefs 10 -hRedistFirstLevel 4 -hRedistStepSize   2 -hRedistNewProcsPerStep  64 -hRedistMaxSteps 2" 




