################################################################################
# Declare a method that allows all sub-cmake-files to add their include paths
# to a common library.
#######################
# Export include paths to global variable.
# PURPOSE: use this function to add local package lib dependencies to global
#          ugLibDependencies property which is used to build libug4
# @param sources sources list to export
function(ExportIncludes includes)
    # iterate over all arguments and insert given prefix
    foreach(includeDir ${ARGV})
            # retrieve the global property ugDependencies and store it
            # in temp variable
            # NOTE: properties must be assigned to variables before being used
            get_property(temp GLOBAL PROPERTY ugIncludes)
            # append tmp to the global ugDependencies property using the correct prefix
            set_property(GLOBAL PROPERTY ugIncludes ${temp} "${includeDir}")
    endforeach(includeDir)
endfunction(ExportIncludes)
