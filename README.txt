How to compile ugbase:
- default build setting includes: common, lib_grid.
  For parallel support see: 'How to enable parallel support'.
  For lib_algebra, lib_discretization see: 'How to compile algebra and discretization module'
- a recent version of cmake is required (min v 2.6), see www.cmake.org
- start in the root folder of UG4 (.../ug4/)
- create a new folder and call it build (mkdir build)
- change path to the build folder (cd build)
- execute cmake on the ug4 path (cmake ../)
  (for parallel support see 'How to enable parallel support')
- execute make (make)
- you now should find all libs of ug4 in ug4/lib
- if you want to recompile the libs later on, you wont have to
  execute cmake again. make will automatically recognise changes
  to the source files and to CMakeLists.txt.

How to compile applications:
- you should find a readme that describes compilation in
  the source-path of each application.

How to enable parallel support:
- to compile ugbase with support for parallelization, you have to
  specify a ParallelToolchain file. To do this proceed as in
  'How to compile ugbase' but instead of calling 'cmake ../' call

  cmake -DCMAKE_TOOLCHAIN_FILE=../ParallelToolchain_###.txt ../

  where ### is replaced by a supported MPI implementaion.
  Choose the ParallelToolchain file that matches your platform:
  
  * ParallelToolchain_OMPI.txt: MacOSX default MPI installation.
  * ParallelToolchain_MPICH.txt: support for the MPICH MPI implementation.
  * ParallelToolchain_MPT.txt: This MPI implementation can be found on some
							linux-clusters (HLRB2-Munich).

How to compile algebra and discretization module:
- per default only the modules common, lib_grid and pcl (if choosen) 
  are build using cmake. In order to build lib_algebra and lib_discretization use:
  cmake - DWITH_ALGEBRA=ON
