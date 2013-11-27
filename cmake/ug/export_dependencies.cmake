################################################################################
# Declare a method that allows all sub-cmake-files to add their dependencies
# to a common library.
#######################
# Export dependencies to global variable.
# PURPOSE: use this function to add local package lib dependencies to global
#          ugLibDependencies property which is used to build libug4
# @param sources sources list to export
function(ExportDependencies deps)
    # iterate over all arguments and insert given prefix
    foreach(l ${ARGV})
            # retrieve the global property ugDependencies and store it
            # in tmp variable
            # NOTE: properties must be assigned to variables before being used
            get_property(tmp GLOBAL PROPERTY ugDependencies)
            # append tmp to the global ugDependencies property using the correct prefix
            set_property(GLOBAL PROPERTY ugDependencies ${tmp} "${l}")
    endforeach(l)
endfunction(ExportDependencies)