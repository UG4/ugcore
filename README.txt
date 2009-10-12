How to compile ugbase:
- a recent version of cmake is required (min v 2.6), see www.cmake.org
- start in the root folder of UG4 (.../ug4/)
- create a new folder and call it build (mkdir build)
- change path to the build folder (cd build)
- execute cmake on the ug4 path (cmake ../)
- execute make (make)
- you now should find all libs of ug4 in ug4/lib
- if you want to recompile the libs later on, you wont have to
  execute cmake again. make will automatically recognise changes
  to the source files and to CMakeLists.txt.

how to compile applications:
- you should find a readme that describes compilation in
  the source-path of each application.
