touch --date="next minute" CMakeCache.txt

if [[ ! $(svnversion -n ..) ]];
then
	echo unknown > svnversion.txt
else
	echo `svnversion -n ..` > svnversion.txt
fi

if [[ ! $(hostname) ]];
then
	echo `hostname` > build_host.txt
else
	echo unknown > build_host.txt
fi

if [[ ! $(date) ]];
then
	echo `date` > compile_date.txt
else
	echo unknown > compile_date.txt
fi


cat > compiledstrings.cpp << EOF
namespace ug{
const char *SVN_REVISION = "$MY_SVN_VERSION";
const char *BUILD_HOST = "$MY_HOSTNAME";
const char *COMPILE_DATE = "$MY_COMPILE_DATE";
}
EOF

exit

if [[ -e compiledstrings.cpp ]]; then
	if [[ $(diff compiledstrings.compiledstrings svnversion.cpp_test) ]];
	then
		mv compiledstrings.cpp_test compiledstrings.cpp
	else
		rm compiledstrings.cpp_test
	fi
else
	mv compiledstrings.cpp_test compiledstrings.cpp
fi
