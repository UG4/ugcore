#!/bin/bash
################################################################################
# Begin: 14062011ih
# Purpose: Example loadleveler command file with some documentation of loadleveler
#          command file and 'mpirun' syntax.
# Usage:   Submit it (after careful editing to specify your job ...) via
#          'llsubmit ll_template' in ug4's 'bin/' directory.
#
# Please do not modify this file for submitting your jobs - work on your own copy!
#
################################################################################

################################################################################
# I. Some infos about the syntax of a load leveler command file (load leveler script)
#    for the "IBM Tivoli Workload Scheduler LoadLeveler" (TWS LoadLeveler)
################################################################################
# 0. The job command file can include the following (cf. p. 200,
#    http://bluegene.epfl.ch/LoadLeveler_UsingAdministrating_a2278818.pdf):
################################################################################
# * A LoadLeveler script is an ordinary unix shell (bash) script.
#
# * Comment statements: You can use comments to document your job command files.
#   You can add comment lines to the file as you would in a shell script. 
#
# * LoadLeveler keyword statements: A keyword is a word that can appear in job command  files.
#   A keyword statement is a statement that begins with a LoadLeveler keyword.
#
# * Shell command statements: If you use a shell script as the executable,
#   the job command file can include shell commands. 
#   According to Best-practise-guide-JUGENE-v0-3.pdf, sec. 3.5.2.3:
#   The application block *is* a regular shell script containing any shell commands!
#
# * LoadLeveler variables (gemeint ist wohl: Zugriff auf deren Inhalt, wie durch die "keyword statements" zuvor gesetzt).
#
################################################################################
# 1. In some detail:
# * Loadleveler keyword statements begin with '# @':
#   '# @ <keyword> = value'
#
#   There can be any number of blanks between the # and the @.
#
# * comments in this script: Everthing else after a hashmark is a comment.
#
#   NOTE: Ein Kommentar muss immer in einer eigenen Zeile stehen
#   (so http://bluegene.epfl.ch/LoadLeveler_UsingAdministrating_a2278818.pdf)!
#
# * 'comment': Comments associated with the job. These comments help to distinguish one job from another job.
#              This parameter is optional.

# * last LL line: '# @ queue'
#
# * remaining part: shell script to be executed
#
# Main jugene specific loadleveler variables (see also 
# http://www2.fz-juelich.de/jsc/jugene/usage/loadl/keywords;
# http://www.prace-project.eu/hpc-training/training_pdfs/talk.pdf):
#
# * 'job_name': Specifies the name of the job. This keyword must be specified in the first job step.
#               If it is specified in other job steps in the job command file, it is ignored. 
#               The job_name only appears in the long reports of the llq, llstatus, and llsummary commands, and in mail related to the job.
#
# * 'job_type': specifies the type of job step to process. 
#
#   Syntax: job_type = serial | parallel | bluegene | MPICH
#   Must be set to 'bluegene', otherwise a serial job on the front end is started without allocating a bluegene partition.
#   Default value: serial
#
# * 'wall_clock_limit': Specifies the "wall clock time" (i.e., the amount of
#                       real time that elapses from start to end of a job).
#                       The minimum 'wall_clock_time' is 30 min (values less than this
#                       produce errors, e.g.: "Wall_Clock_Limit (00:03:00) too small! Minimum is 30 min").
#
# TO CHECK: Is it possible to avoid wasting of compute time with 'job_cpu_limit'?
#       But this limit might be also quite high, as the limit for 'wall_clock_limit' is ...
#
# * 'bg_connection': sets the topology:
#
#   + TORUS        Specifies that the admissible partitions must be wireable as a torus.
#   + MESH         Specifies that the admissible partitions must be wireable as a mesh.
#   + PREFER_TORUS Specifies that the admissible partitions should be wireable as a torus,
#                  but if there are no such partitions then the selected partition must be wireable as a mesh. 
#
#   This keyword is only valid for job type bluegene.
#   This keyword cannot be used if the bg_partition keyword is specified.
#   This keyword is not inherited by other job steps. 
#   Default value: MESH is the default value. 
#
#   http://www.prace-ri.eu/IMG/pdf/Best-practise-guide-JUGENE-v0-3.pdf:
#   + The TORUS network interconnects all compute nodes and has a topology of a three-dimensional torus,
#     i.e. each node has six nearest neighbours 
#     sect. 1.5.1: "... torus topology is available to the user application only if the assigned partition of the 
#                machine is a multiple of the midplane (16 node cards, half of a rack)"
#
#   + The MESH network is a global collective tree network. It interconnects all compute and I/O nodes.
#
# * 'bg_size': requested number of compute nodes.
#
# * 'bg_shape': requested "torus shape" specified in number of midplanes in x, y and z direction.
#
#               Syntax: bg_shape = XxYxZ 
#
#               where X, Y, and Z are positive integers indicating the number of base partitions (midplanes) in the X-direction,
#               Y-direction, and Z-direction, respectively, of the requested job shape.
#               The values of X, Y, and Z or their rotations, if bg_rotate is true, must not be greater than the corresponding
#               X, Y, and Z sizes of the Blue Gene system, otherwise the job will never be able to start.
#
#               This keyword is only valid for job type bluegene.
#               This keyword can not be used if the bg_partition or bg_size keyword is specified.
#               This keyword is not inherited by other job steps.
#
#               NOTE: The X, Y, and Z dimensions of the allocated partition will be exactly as defined by the 'bg_shape' job
#                     command file keyword unless the job command file keyword 'bg_rotate' is specified as 'true',
#                     in which case all possible rotations of the dimensions are possible. 
#
#               Default value: No default is set.
#
#               http://www.prace-ri.eu/IMG/pdf/Best-practise-guide-JUGENE-v0-3.pdf: The max shape on JUGENE is 9x4x4.
#               In der Tat ergibt das 144 midplanes (72 racks mal 2 midplanes/rack)!
#
# NOTE 1: The size of a job has to be specified by '# @ bg_size'   OR by  '# @ bg_shape'!!!
#         (Cf. http://www2.fz-juelich.de/jsc/jugene/usage/loadleveler/).
#
#         Falls man aber mit 'mpirun'-Parameter '-mapfile' ein *explizites* Mapfile (siehe unten) spezifizieren moechte,
#         MUSS man die Job-Groesse per 'bg_shape' angeben!
#
# NOTE 2: Im ersten Teil des Outputs, produziert durch 'env | grep LOADL_' (s.u.), wird fuer beides eine 'LOADL_*'-Variable angezeigt,
#         auch wenn der Job durch z. B. '@ bg_size = 2048' spezifiziert wurde:
#            ...
#            LOADL_BG_SIZE=2048
#            ...
#            LOADL_BG_SHAPE=2x1x2
#            ...
#         Und spaeter im Output des Backend-MPI heisst es fuer einen solchen Job dann, gemaess der im 'mpirun'-Command gewaehlten Anzahl Prozesse:
#         "... FE_MPI (Debug): jobenvs[]      = <BG_SIZE=512 ..."
#
# * 'bg_rotate': allow the requested shape to be rotated (damit der LoadLeveler leichter/schneller eine freie "Luecke" finden kann).
#
# * 'res_id': specifies the reservation to run on (not yet used til now).
#
# * Evaluation of LoadLeveler variables (defined by one of these commands): '$(<variable>)' (thus, same syntax as for make-variables).
#
# ONLY Multistep loadleveler scripts:
# dependency - Specifies  the dependencies between job steps. A job dependency, if used in a given  job step, must be explicitly specified for that step.
# Syntax: 
#    dependency = step_name operator value
#
#    'step_name' is the name of a previously defined job step (as described in the step_name keyword). 
#    'operator' is one of the usual C boolean operators, and
#    'value' is usually a number that specifies the job return code to which the 'step_name' is set -- etwas unpraezise Beschreibung ... Siehe Beispiele!
################################################################################
# 2. mpirun main parameters:
# -np <#PE> - specifies the number of MPI processes, as usual.
#
# -exe <executable> - specifies the name of the executable to be processed in parallel.
#
# -mode {VN | DUAL | SMP} - specifies the "execution mode":
#        VN:   "Virtual Mode" (official name is "Quad mode"): 4 MPI processes per CN,
#              1/4 RAM of the CN per core/process, i.e. ~512 Mbyte.
#        DUAL: "Dual Mode": 2 MPI processes per CN; the other two cores could be idle
#              or execute OpenMP threads (=> hybride parallelisation).
#              1/2 RAM of the CN per core/process, i.e. ~1 Gbyte.
#        SMP:  "Symmetrical Multi Processing Mode": 1 MPI process per CN; the other
#              three cores could be idle or execute OpenMP threads.
#              The complete RAM (~2 Gbyte) of the CN is available for the one MPI process.
#
#        Obviously "VN" is the preferred execution mode if large numbers of processes should be achieved -
#        and ug4 works with VN mode (at least up to ~64 Ki DoFs per process)! :-)
#
# -verbose {0 | 1 | 2 | 3 | 4} - specifies a level of verbosity (larger number => more output).
#                                http://www.scc.acad.bg/documentation/Application_User_Worshop.pdf:
#                                The default is 0 which means that mpirun will not output any
#                                status or diagnostic messages unless a severe error occurs.
#                                If you are curious as to what is happening try levels 1 or 2.
#                                All mpirun generated status and error messages appear on STDERR.
#
# -env NAME=VALUE              - for passing environment variables
#
# -h                           - for extended help, run at front-end
#
# -mapfile <mapfile|mapping>   - mapfile contains a user specified MPI topology,
#                              - mapping is a permutation of XYZT
#
# http://www.prace-project.eu/hpc-training/training_pdfs/talk.pdf:
# explicit map file 
#   usage: mpirun -mapfile filename
#   each line contains the torus (4D;x y z t) coordinates of one rank;
#   line 1 sets the coordinates of rank 0, 2 of 1, ... 
#   must specify BG_SHAPE=XxYxZ and BG_ROTATE=FALSE.
# Jedoch: Laut Support muss man nicht notwendigerweise 'bg_rotate = false' setzen, sobald man mit 'bg_shape' arbeitet,
# wie man eine entsprechende Aussage in 'http://www.prace-project.eu/hpc-training/training_pdfs/talk.pdf' interpretieren koennte!
# Ich vermute, das ist erst noetig, wenn man mit *expliziten* Mapfiles arbeitet!
#
# http://www2.fz-juelich.de/jsc/jugene/usage/mpirun/:
# The user can specify with this option the order in which the MPI tasks are distributed across the nodes and cores of the partition reserved.
# The standard mapping on JUGENE is to place the tasks in XYZT order,
# where X,Y, and Z are the torus coordinates of the nodes in the partition and T is the number of the cores within each node (T=0,1,2,3).
#
# When the tasks are distributed across the nodes the first dimension is increased first,
# i.e. for XYZT the first three tasks would be executed by the nodes with the torus coordinates <0,0,0,0>, <1,0,0,0> and <2,0,0,0>
# 'mapfile' can either be a permutation of X,Y,Z and T or the name of a mapfile in which the distribution of the tasks is specified.
#
# NOTE: Nach 'http://www.ichec.ie/support/tutorials/bg_port.pdf' gilt dieses Beispiel allgemein, nicht nur fuer SMP-Mode!
# Zitat: Once each node of the partition will have a running process allocated on, then the T coordinate will be used
# (T probably standing for thread).
#
# The drawback of this approach is that processes of nearby ranks will be spread across nodes instead of being packed into the same nodes.
################################################################################

################################################################################
################################################################################
# II. What we are interested in:
#
# Weak scaling measurements, A. Vogel, Nehalem cluster: 32 - 2048 PE's, 131'841 dof's per PE, coarse grid 8 x 8 quadrilaterals
#    8 PE: Level  7 ==> (2^3 * 2^{ 7} + 1)^2 nodes =   1'050'625 nodes, ~131'328 nodes/PE (output of proc 0: 131841 dofs in total)
#   32 PE: Level  8 ==> (2^3 * 2^{ 8} + 1)^2 nodes =   4'198'401 nodes, ~131'200 nodes/PE (output of proc 0: 131841 dofs in total)
#  128 PE: Level  9 ==> (2^3 * 2^{ 9} + 1)^2 nodes =  16'785'409 nodes, ~131'136 nodes/PE (output of proc 0: 131841 dofs in total)
#  512 PE: Level 10 ==> (2^3 * 2^{10} + 1)^2 nodes =  67'125'249 nodes, ~131'104 nodes/PE (output of proc 0: 131841 dofs in total)
# 2048 PE: Level 11 ==> (2^3 * 2^{11} + 1)^2 nodes = 268'468'225 nodes, ~131'088 nodes/PE (output of proc 0: 131841 dofs in total)
# additional (at least ...):
# 8192 PE: Level 12 ==> (2^3 * 2^{12} + 1)^2 nodes = 1'073'807'361 nodes, ~131'080 nodes/PE (output of proc 0: 131841 dofs in total)
#
################################################################################
################################################################################

################################################################################
# III. Example job (note again, that comments must be on a separate line):
################################################################################
# TO CHECK: Job-Name and comment:
# @ job_name = ug4_laplace
# @ comment = "BGP First test"
# @ error = $(job_name).$(jobid).out
# @ output = $(job_name).$(jobid).out
# @ environment = COPY_ALL;
########################################
# TO CHECK: Minimal 'wall_clock_limit' depends obviously on job size:
#           Up to 'bg_size =  16384' (nec. for  64k PE) 30 min are enough,
#           for   'bg_size =  65536' (nec. for 256k PE) minimum is 45 min (acc. to error message of LoadLeveler which refuses to submit the job)
# @ wall_clock_limit = 00:30:00
########################################
# @ notification = error
# @ notify_user = ingo.heppner@gcsc.uni-frankfurt.de
# @ job_type = bluegene
# Connection type in [TORUS| MESH | PREFER_TORUS] - default: MESH; for TORUS 'bg_size' must be >= 512
# @ bg_connection = TORUS
########################################
# IMPORTANT: Specification of 'bg_size':
#    "Durch das Batch-Script wird die in bg_size angegebene Partitionsgroesse
#    angefordert. Diese muss dann daher auch vom Scheduler reserviert werden,
#    bestimmt die Klasse des Jobs und wird auch anschliessend abgerechnet"
#    (statement by JuGene-Support-Team).
#
#    Best-practise-guide-JUGENE-v0-3.pdf, p. 32:
#    'bg_size' is the size of the Blue Gene job in *number of compute nodes* to be used.
#    'bg_size = 32' is the minimum on JUGENE [sc.: if 'connection_type' != 'TORUS'].
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
# @ bg_size = 512
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
# Set ug executable, arguments and load path (to shorten the lengthy commands):
##########################################
UGSHELL=ugshell
UGARGS="-ex ../scripts/tests/modular_scalability_test.lua -dim 2 -grid ../data/grids/unit_square_01/unit_square_01_quads_8x8.ugx -lsIterator gmg -lsMaxIter 100 -verb 0 -numPreRefs 3"

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
#mpirun -np      4 -exe ./$UGSHELL -mode VN -mapfile TXYZ -verbose 2 -args "$UGARGS -numRefs  6"
# job "jugene4b.292615":
#mpirun -np     16 -exe ./$UGSHELL -mode VN -mapfile TXYZ -verbose 2 -args "$UGARGS -numRefs  7"
# job "jugene4b.292616":
#mpirun -np     64 -exe ./$UGSHELL -mode VN -mapfile TXYZ -verbose 2 -args "$UGARGS -numRefs  8"
# job "jugene4b.292617":
#mpirun -np    256 -exe ./$UGSHELL -mode VN -mapfile TXYZ -verbose 2 -args "$UGARGS -numRefs  9"
# job "jugene4b.292618":
# Hierarchical redistribution: Last stage from 64 to 1 Ki, first from 1 to 64:
#mpirun -np   1024 -exe ./$UGSHELL -mode VN -mapfile TXYZ -verbose 2 -args "$UGARGS -numRefs 10 -hRedistFirstLevel 5 -hRedistStepSize 100 -hRedistNewProcsPerStep  16"
#### Up to this with 'bg_size = 512'! ####
##########################################

# ...