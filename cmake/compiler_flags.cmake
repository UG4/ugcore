function(add_cxx_flag flag)
    if(NOT ${CMAKE_BUILD_TYPE})
       message(fatalerror, "No build type set. Your cflags would have been ignored.")
       message(fatalerror, "do not call add_cxx_flag before calling add_cxx_flag")
    endif()
    
    # reuse initial cflag of build type, if set, once.
    # eg. debug cflag -g etc.
    if(NOT FIRST_CALL_ADD_CXX_FLAG)
        set(FIRST_CALL_ADD_CXX_FLAG 42 PARENT_SCOPE)
            foreach(lang C CXX)
                 set(CMAKE_${lang}_FLAGS_${CMAKE_BUILD_TYPE} 
                    "${CMAKE_${lang}_FLAGS_${CMAKE_BUILD_TYPE}_INIT}"
                    CACHE STRING "overriden flags!" FORCE)
                 #message("set init flag to: ${CMAKE_${lang}_FLAGS_${CMAKE_BUILD_TYPE}_INIT}")
            endforeach()
    endif()

    # now append given flag to build type flags
    foreach(lang C CXX)
      set(CMAKE_${lang}_FLAGS_${CMAKE_BUILD_TYPE} 
        "${CMAKE_${lang}_FLAGS_${CMAKE_BUILD_TYPE}} ${flag}" 
        CACHE STRING "overriden flags!" FORCE)
    endforeach()
endfunction(add_cxx_flag)