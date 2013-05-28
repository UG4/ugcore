################################################################################
# Declare a method that allows all sub-cmake-files to add their sources
# to a common library.

#######################
# Export sources to global variable.
# PURPOSE: use this function to add local package sources to global
#          ugSources property which is used to build libug4
# @param prefix current directory prefix
# @param sources sources list to export
function(ExportSources prefix sources)
    # iterate over all arguments and insert given prefix
    foreach(l ${ARGV})
        # FIXME: this is a hack to omit the first argument
        #        which is the prefix. Shall we use boolean or index variable?
        if(NOT "${l}" STREQUAL "${prefix}")
            # retrieve the global property ugSources and store it
            # in tmp variable
            # NOTE: properties must be assigned to variables before being used
            get_property(tmp GLOBAL PROPERTY ugSources)
            # append tmp to the global ugSources property using the correct prefix
            if("${prefix}" STREQUAL "")
				set_property(GLOBAL PROPERTY ugSources ${tmp} "${l}")
			else("${prefix}" STREQUAL "")
				set_property(GLOBAL PROPERTY ugSources ${tmp} "${prefix}/${l}")
			endif("${prefix}" STREQUAL "")
        endif(NOT "${l}" STREQUAL "${prefix}")
    endforeach(l)
endfunction(ExportSources)