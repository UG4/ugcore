#!/bin/sh
# execute UG test-suite
# ------------------------
# testsuite arguments
args='--output_format=XML --log_level=all --report_level=no --run_test=*NumProc$np --log_sink=utf_log_np$np.xml'

# has to be set to find libug4.so
LD_LIBRARY_PATH=$WORKSPACE/trunk/lib/
# path to testsuite
ts=$WORKSPACE/trunk/bin/testsuite
# should be "serial" or "parallel"
mode=$1
# run testsuite in temporary directory in jobs workspace
dir=$(mktemp -d --tmpdir=$WORKSPACE)
cd $dir


case "$mode" in

# run testsuite in serial mode
serial)
	$ts $args
;;

# run testsuite in parallel mode for 1 till 16 processors
parallel)
	# if run_test filter matches no tests, continue
	for np in {1..16}; do 
  		mpirun -n $np $ts $args || continue
	done
;;
esac

# delete empty log files (necessary due to empty testsuites and we do not want to archive them)
find -empty -name "ug_test_numprocs*.log" -exec rm {} +
find -empty -name "*.xml" -exec rm {} +

# copy utf logs and ug logs to workspace for archivation
mv *.xml *.log $WORKSPACE

# delete temporary run directory
rm -rf $dir

# return success, no matter what was the last return code...
#exit 0
