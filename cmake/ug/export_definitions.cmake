################################################################################
# Declare a method that allows all sub-cmake-files to add their own definitions
# to the main project (P_UG4) defines
#######################
# Export defines to global variable.
# PURPOSE: use this function to add local project defines to global
#          ug4LibDefinitions property which is used to build libug4 project
# @param definitions list of definitions to be exported
function(ExportDefinitions definitions)
    # iterate over all arguments and insert given prefix
    foreach(definition ${ARGV})
            # retrieve the global property ugdefinitions and store its values
            # in a temp variable
            # NOTE: properties must be assigned to variables before being used
            get_property(temp GLOBAL PROPERTY ugDefinitions)
            # append tmp to the global ugDependencies property using the correct prefix
            set_property(GLOBAL PROPERTY ugDefinitions ${temp} "${definition}")
    endforeach(definition)
endfunction(ExportDefinitions)
