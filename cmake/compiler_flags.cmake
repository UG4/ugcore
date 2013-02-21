# This function sets the CMAKE_C_FLAGS and CMAKE_CXX_FLAGS variabels.
# 
# TODO: maybe add sanity checks for flags
#
function(add_cxx_flag flag)
	foreach(lang C CXX)
		set(CMAKE_${lang}_FLAGS
			"${CMAKE_${lang}_FLAGS} ${flag}"
			CACHE STRING "overriden flags!" FORCE)
	endforeach()
endfunction(add_cxx_flag)


# used to clear c/cxx flags
function(reset_cxx_flags)
	foreach(lang C CXX)
			set(CMAKE_${lang}_FLAGS "" CACHE STRING "clear flags" FORCE)
	endforeach()
endfunction(reset_cxx_flags)
