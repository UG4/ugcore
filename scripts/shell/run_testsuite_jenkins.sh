#!/bin/sh -ex
# execute UG test-suite
# ------------------------
# xml input files
ugroot=$WORKSPACE
unit_test_data=$ugroot/apps/unit_tests/data/script_tests
core_tests=$unit_test_data/script_test_param.xml
experimental_test=$unit_test_data/experimental_plugins.xml
validate_schema=$unit_test_data/ScriptParamMappingSchema.xsd
 
# testsuite arguments
defargs='--output_format=XML --log_level=all --report_level=no --log_sink=utf_log_${mode}_np$np.xml'
# note script params defaults to $core_tests
testcore_args='-script_params $core_tests --run_test=*NumProc$np'
# run testsuite LuaScripts for plugins
testplugins_args='-script_params $experimental_test --run_test=/LUAScriptsNumProc$np'

# has to be set to find libug4.so
export LD_LIBRARY_PATH=$ugroot/lib/
# path to testsuite
ts=$ugroot/bin/testsuite
# should be "serial" or "parallel"
mode=$1
# should experimental plugins tested?
plugins=$2

if [ "$mode" != "parallel" -a "$mode" != "serial" ]; then
	echo "mode not set: should be parallel or serial"
	exit 1
fi

# if plugins is set to something, use plugin argument for testsuite
# else use core arguments
if [ $plugins ]; then
	echo 'use plugin args'
	additional_args=$testplugins_args
else
	additional_args=$testcore_args
fi

# first of all validate input xml files againt schema
which xmllint
if [ $? -eq 0 ]; then
	if [ $plugins ]; then
		file=$experimental_test
	else
		file=$core_tests
	fi
	xmllint --schema $validate_schema --noout $file  
	if [ $? -gt 0 ]; then
		echo "$file does not validate against schema"
		exit 1
	fi
else
	echo "warning: input files did not get validated, because xmllint is missing"
fi

# run testsuite in temporary directory in jobs workspace
dir=$(mktemp -d --tmpdir=$WORKSPACE)
cd $dir

case "$mode" in
# run testsuite in serial mode
serial)
	np=1
	eval $ts $defargs $additional_args || true 
;;

# run testsuite in parallel mode for 1 till 16 processors
parallel)
	# if run_test filter matches no tests, continue
	for np in {1..16}; do 
		eval mpirun -n $np $ts $defargs $additional_args || continue
	done
;;
esac

# delete empty log files (necessary due to empty testsuites and we do not want to archive them)
find -empty -name "ug_test_numprocs*.log" -exec rm {} +
find -empty -name "*.xml" -exec rm {} +

# copy utf logs and ug logs to workspace for archivation with jenkins
mv *.xml *.log $WORKSPACE

# delete temporary run directory
rm -rf $dir
