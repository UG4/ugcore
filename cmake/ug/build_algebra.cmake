# included from ug_includes.cmake
########################################
# buildAlgebra
if(buildAlgebra)
	add_definitions(-DUG_ALGEBRA)
	
	# we'll check for lapack and blas here
	if(LAPACK OR BLAS)
		# On OSX, we know where lapack and blas are located. To avoid errors
		# with Fortran compilers on OSX, we'll add APPLE as a special case here.
		if(APPLE)
			if(LAPACK)
				set(LAPACK_LIBRARIES "-framework vecLib" CACHE STRING "LAPACK library" FORCE)
				FIND_PATH(LAPACK_INCLUDE_PATH clapack.h /usr/local/include/ /usr/include /include)
				if(LAPACK_INCLUDE_PATH)
		  			set(LAPACK_FOUND YES)
		  		endif(LAPACK_INCLUDE_PATH)
			endif(LAPACK)
			
			if(BLAS)
				set(BLAS_LIBRARIES "-framework vecLib" CACHE STRING "CBLAS library" FORCE)
				find_path(BLAS_INCLUDE_PATH cblas.h /usr/local/include/ /usr/include /include)
				set(BLAS_FOUND YES)
				if(BLAS_INCLUDE_PATH)
					set(BLAS_FOUND YES)
				endif(BLAS_INCLUDE_PATH)
			endif(BLAS)
			
		elseif(NOT BUILTIN_BLAS AND NOT BUILTIN_LAPACK)			
			# a fortran compiler is required to find the packages
			# Sadly there seems to be a cmake-bug, hence the following workaround.
			# ENABLE_LANGUAGE workaround begin (issue 0009220)
			message(STATUS "Info: If problems with the Fortran compiler occur, consider deactivating LAPACK and BLAS")
			if(DEFINED CMAKE_Fortran_COMPILER AND CMAKE_Fortran_COMPILER MATCHES "^$")
			  set(CMAKE_Fortran_COMPILER CMAKE_Fortran_COMPILER-NOTFOUND)
			endif(DEFINED CMAKE_Fortran_COMPILER AND CMAKE_Fortran_COMPILER MATCHES "^$")
			# ENABLE_LANGUAGE workaround end (issue 0009220)
			
			ENABLE_LANGUAGE(Fortran OPTIONAL)
			if(CMAKE_Fortran_COMPILER_WORKS)
				if(LAPACK)
					find_package(LAPACK)
				endif(LAPACK)
				
				if(BLAS)
					find_package(BLAS)
				endif(BLAS)
				
			endif(CMAKE_Fortran_COMPILER_WORKS)
		endif(APPLE)
		
	# We'll output whether lapack and blas are used, to avoid misconceptions
		if(LAPACK_FOUND)
			message(STATUS "Info: Using Lapack")
			include_directories (${LAPACK_INCLUDE_PATH})
			set(linkLibraries ${linkLibraries} ${LAPACK_LIBRARIES})
			add_definitions(-DLAPACK_AVAILABLE)
		elseif(BUILTIN_LAPACK)
			message(STATUS "Info: Using Builtin Lapack")
			add_definitions(-DLAPACK_AVAILABLE)
		else(LAPACK_FOUND)	
			message(STATUS "Info: Not using Lapack. No package found.")
		endif(LAPACK_FOUND)
		
		if(BLAS_FOUND)
			message(STATUS "Info: Using Blas")
			include_directories (${BLAS_INCLUDE_PATH})
			set(linkLibraries ${linkLibraries} ${BLAS_LIBRARIES})
			add_definitions(-DBLAS_AVAILABLE)
		elseif(BUILTIN_BLAS)
			message(STATUS "Info: Using Builtin Blas")
			add_definitions(-DBLAS_AVAILABLE)
		else(BLAS_FOUND)	
			message(STATUS "Info: Not using Blas. No package found.")
		endif(BLAS_FOUND)
		
		if(BLAS_goto2_LIBRARY)
			message(STATUS "Info: GotoBLAS2 found: (${BLAS_goto2_LIBRARY}). Adding -lgfortran.")
			SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lgfortran")
		endif()
	endif(LAPACK OR BLAS)
endif(buildAlgebra)