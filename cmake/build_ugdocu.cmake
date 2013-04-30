if(BUILD_UGDOCU)
	add_custom_target(updateUGDocu ALL )
	add_custom_command(TARGET updateUGDocu
						 PRE_BUILD
						 COMMAND "${UG_ROOT_PATH}/bin/ugdocu" -silent -list
						 WORKING_DIRECTORY ${UG_ROOT_PATH}/bin)
	add_dependencies(updateUGDocu ${TARGET} ugdocu)
endif(BUILD_UGDOCU)