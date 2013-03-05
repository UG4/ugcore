# This function sets the CMAKE_C_FLAGS and CMAKE_CXX_FLAGS variabels.
# 
# TODO: maybe add sanity checks for flags
#
function (add_c_flag flag)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${flag}" CACHE STRING "overriden flags!" FORCE)
endfunction(add_c_flag)

function (add_cpp_flag flag)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${flag}" CACHE STRING "overriden flags!" FORCE)
endfunction(add_cpp_flag)

function(add_cxx_flag flag)
	add_c_flag(${flag})
	add_cpp_flag(${flag})
endfunction(add_cxx_flag)


# used to clear c/cxx flags
function(reset_cxx_flags)
	foreach(lang C CXX)
			set(CMAKE_${lang}_FLAGS "" CACHE STRING "clear flags" FORCE)
	endforeach()
endfunction(reset_cxx_flags)
