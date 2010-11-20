########################################################################
How to compile ugbase:
- a recent version of cmake is required (min v 2.6), see www.cmake.org

- a recent version of boost should be located somewhere on your machine.
  At least boost 1.34 is required.
  UG only uses boost headers. The boost library thus doesn't have to be build.
  Consider setting an environment-variable BOOST_ROOT to the path that contains boost or
  to add the boost path to your PATH environment variable.
  (this could look like this: export PATH=$PATH:$HOME/libs/boost_1_41_0)

- start in the root folder of UG4 (.../ug4/)

- create a new folder and call it build (mkdir build)

- change path to the build folder (cd build)

- execute cmake on the ug4 path (cmake ../)

- execute make (make)

- you now should find all libs of ug4 in ug4/lib

- if you want to recompile the libs later on, you wont have to
  execute cmake again. make will automatically recognise changes
  to the source files and to CMakeLists.txt.


########################################################################
How to compile applications:
- you should find a readme that describes compilation in
  the source-path of each application.

- If you want to write a new application please take a look at ug4/ug_cmake_includes.txt


########################################################################
How to compile algebra and discretization module:
- note: lib_algebra uses two external libraries: boost and hyprelib. If you want to 
        compile lib_algebra correctly please add the directories of those libraries 
        to your standard search path, e.g.

	export PATH=$PATH:/Users/andreasvogel/Software/Boost
	export PATH=$PATH:/Users/andreasvogel/Software/hypre (optional)


########################################################################
You can enable and disable some parts of ug.

You can enable and disable some parts of ug, in order to reduce build-times or
unwanted dependencies.
It is good practice to disable everything you don't need.
The following options can be enabled / disabled:

- UG_DEBUG: Produces slow but debuggable applications.
- NO_MPI: Disables MPI support.
- NO_ALGEBRA: Disables support for algebra and discretization.
- NO_BOOST: Disables support for boost (should only be disabled if required).
- NO_PCL: Disables creation of the pcl-library.

To enable or disable options during the cmake call, call cmake like this:
cmake -DUG_DEBUG=ON -DNO_MPI=OFF SomePath

If you want to disable some parts of ug in the CMakeLists.txt file of your
application (a good idea if your application doesn't need some part of ug)
you can set those options before including the ug_cmake_includes file using
set(NO_MPI ON)
set(NO_ALGEBRA OFF)
include(ug_cmake_includes.txt)

Please note that UG_DEBUG should not be set in any CMakeLists.txt file. It should
only be used as a command-line parameter to cmake.

