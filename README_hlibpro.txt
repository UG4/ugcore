################################################################################
Coupling of HLibPro (http://www.hlibpro.com/) with ug4 (begin march 2011)
################################################################################

I.   General:
#############
HLibPro is a parallel H-matrix library for shared and distributed memory machines.
Hierarchical matrices, or H-matrices for short, are a powerful tool for numerical
applications allowing the usage of the complete matrix algebra, e.g. matrix
multiplication and inversion, with almost linear complexity (short description
from http://www.hlibpro.com/).

NOTE: HLibPro is an external library, for which we do not provide the sources!
Nevertheless our group has access to the HLibPro sources within "ASIL" project.
Otherwise you may get a license contacting Dr. Ronald Kriemann, rok@mis.mpg.de.


II.  Short description of HLibPro functionality provided in ug4 so far:
#######################################################################
Some functionality of HLibPro is used to define a sparse matrix solver (type
'IMatrixOperatorInverse'), based on the functions of the HLibPro C interface.
See 'ugbase/lib_algebra/operator/linear_solver/hlibpro.h'.

For importing the system matrix assembled by UG to HLibPro a sparse matrix format
is used. This format is implemented by a class 'CRSMatrix', which is also defined
in 'hlibpro.h'.

The solver is registered as 'HLIBSolver' in 'lib_algebra_bridge.cpp'.


TODO:
Import of coordinates and construction of cluster trees based on coordinates
aren't implemented yet.



III. Installation/Usage (Unix only):
####################################
To use this HLibPro functionality from within ug4,

1. get a copy of HLibPro (sources are not provided with ug4!) and place it
   somewhere on your filesystem, e.g. ~/bin/hlibpro-0.13.6/, then configure
   and build HLibPro. See below for some hints to do this!

2. Check 'cmake/ug_includes.cmake': If necessary adapt path after 'PATH_SUFFIXES'
   (last argument of 'find_path()' call). In the moment this path is given
   relative to the location of 'cmake'.

   E.g., if your HLibPro main directory resides in '~/bin/', then use

   find_path (HLIBPROLIB_DIR libhpro.a
              PATHS ENV PATH
              PATH_SUFFIXES ../bin/hlibpro-0.13.6/lib/ )

   This works e.g. on 'cekon.gcsc.uni-frankfurt.de'. Sometimes it seems to be
   necessary to give 'PATH_SUFFIXES' relative to the ug4 main directory (observed
   e.g. on 'quadruped.gcsc.uni-frankfurt.de'):

   find_path (HLIBPROLIB_DIR libhpro.a
              PATHS ENV PATH
              PATH_SUFFIXES ../../bin/hlibpro-0.13.6/lib/ )

   Everything else (e.g. include paths, including the path to the C interface
   header file 'hlib-c.hh') depends on 'PATH_SUFFIXES' and is automatically set.


3. Configure ug4 for using HLibPro, e.g.:

   % cmake <other options> -DHLIBPRO=ON  ../ # Enable HLibPro

   To disable usage of HLibPro:
   % cmake <other options> -DHLIBPRO=OFF ../ # Disable HLibPro

4. Compile ug4.

5. Script code for testing HLibPro is provided in

      scripts/hlibtest.lua

   Execute this script by running e.g.

   % mpirun -n 4 ./ugshell -ex scripts/hlibtest.lua

   For (preliminary) tests for systems of pdes see
      scripts/systemlaplace_hlib.lua
      scripts/navierstokes_hlib.lua

IV. Installation of HLibPro:
############################
To use the HLibPro configuration system, a Python interpreter is needed.
Furthermore, the build system for HLibpro is 'SCons' (http://www.scons.org),
which also is based on Python. 

For further information to the build process of HLibpro we refer to the
documentation included in the HLibpro sources, but for convenience we sketch the
typical steps here (version 0.13.6; we assume, that the tar ball is placed in
'~/bin/' (which probably isn't the most elegant place ...); consider also the
"remarks" below):

(a) extract tar ball:

    % cd ~/bin/
    % tar xzf hlibpro-0.13.6.tgz

(b) change into HLibPro main directory and configure HLibPro:

    % cd hlibpro-0.13.6
    % ./configure

(c) compile HLibPro via 'scons' by simply typing

    % scons 

(d) Maybe execute the stand-alone test executables now available in the sub
    directory 'examples/':

    % cd examples
    % bem1d 100

    etc.


REMARKS:

1. Python:
HLibPro's 'configure' script needs a "not to old python version", especially
one whichs knows "subprocesses"! E.g. Python 2.4.3 is ok, Python 2.3.4 is not!
So, if you get a message like the following

   quadruped@/home3/ingo/bin/hlibpro-0.13.6> ./configure
   Traceback (most recent call last):
     File "./configure", line 11, in ?
       import subprocess
   ImportError: No module named subprocess

You have to update your Python installation first!

To do this (or to install it for the first time), get a copy of its current
sources (see http://www.python.org) and proceed as follows (this method
does not require root access and is also feasible if there is no package
available for your machine; in this example we place it in '~/bin/'):

   % cd ~/bin/
   % wget http://www.python.org/ftp/python/2.7.1/Python-2.7.1.tar.bz2
   % bunzip Python-2.7.1.tar.bz2
   % tar xvf Python-2.7.1.tar
   % cd Python-2.7.1

The build process consists of the usual

   % ./configure
   % make
   % make install # (if you have appropriate rights)

in the "python main directory".


After all, to finally configure HLibPro you have to check the first line in
HLibPro's 'configure' script,

   #!/usr/bin/env /usr/bin/python

and adapt it if needed to point to your new Python interpreter, e.g. change
it to something like:

   #!/usr/bin/env /home3/ingo/bin/Python-2.7.1/python

Then you can type (as mentioned above):

    % cd ~/bin/hlibpro-0.13.6
    % ./configure

to configure HLibPro.



2. SCons:
If you have to build and install 'SCons', get a copy of its current sources
unpack them, change into the "scons main directory" and execute its Python-
based setup script, e.g. (again we place everything in '~/bin/'):


   % cd ~/bin/
   % wget http://prdownloads.sourceforge.net/scons/scons-2.0.1.tar.gz
   % gunzip scons-2.0.1.tar.gz 
   % tar xvf scons-2.0.1.tar

   % cd scons-2.0.1
   % python setup.py install --prefix=<install dir> # e.g. '--prefix=~/bin/scons'!

To finally compile HLibPro execute this version of 'scons' in the "HLibPro
main directory" e.g. by typing:

   % cd hlibpro-0.13.6
   % ~/bin/scons/bin/scons



3. LAPACK
Internally, HLibpro uses LAPACK for most arithmetic operations. For the case no
such implementation is available, HLibpro can use a modified version of CLAPACK
as a substitute for LAPACK. CLAPACK is contained in the sources of HLibPro.

By default HLibpro is linked against the substitute library, 'libclapack.a'
("by default" according to the provided documentation, 'hlibpro-user.pdf',
but I'm in doubt about that, at least on OS X  ...), which might result in a
reduced performance of HLibpro.

In case an (probably) optimised version of LAPACK is available for your system
it is therefore highly recommended to use this instead of CLAPACK. To do so,
configure HLibpro this way:

   % ./configure --lapack=-Llapack


ATTENTION: There are cases where the substitute library 'libclapack.a' is
needed (see "Known problems below"). To use CLAPACK, type

   % ./configure --lapack=-CLAPACK'



4. Known problems:
4.1
While building the executables in the 'examples/' directory the linker may
drop the following error message:

    /usr/bin/ld: Undefined symbols:
    __Unwind_Resume
    collect2: ld returned 1 exit status
    scons: *** [examples/bem1d] Error 1
    scons: building terminated because of errors.

This seems to be a quite common problem with older versions of OS X / XCode when
linking C++ software with gcc instead of g++.


Remedy:
Use 'g++' for linking the HLibPro examples! Reason: 'g++' automatically links
other system libraries than 'gcc' does, so that the otherwise undefined symbol
can be resolved!

To achieve this, do:

   % ./configure --cc=g++


4.2
(a) When attempting to start one of HLibPro's example executables maybe the
following error message shows up, e.g.

   % bem1d 100
   % bem1d: symbol lookup error: bem1d: undefined symbol: slamch_

This behavior seems to be due to a too old (?) version of LAPACK on your system!
It was observed with 3.0.3.(currently installed on 'quadruped'). In contrast,
e.g. with Lapack 3.1.1 (currently installed on 'cekon') this problem does not
occur.

Remedy: Use the substitute library 'libclapack.a' instead, i.e. configure
HLibPro accordingly (cf. above):

   % ./configure --lapack=CLAPACK'

to enforce building of 'libclapack.a'.

(b) Independent of (but obviously directly related to) the problem above you may
run in a similar problem when attempting to start 'ugshell' compiled with HLibPro:

   % ugshell -ex ../scripts/laplace.lua
   % ugshell: symbol lookup error: ugshell: undefined symbol: slamch_

This may happen even if the stand-alone examples of HLibPro work properly - which
is the case if HLibPro is linked against CLAPACK, but ug4's configuration via
'cmake' has automatically chosen the (too old) version of LAPACK installed on
your system!


Remedy:
Link ug4 also against the "LAPACK substitute", 'libclapack.a'!

To achieve this, add also 'libclapack.a' to the 'UG_LIBRARIES', i.e., edit the
appropriate part in 'ug_cmake_includes.txt':

   if(HLIBPRO AND HLIBPROLIB_DIR)
      set(UG_LIBRARIES ${UG_LIBRARIES} ${HLIBPROLIB_DIR}/libhpro.a ${HLIBPROLIB_DIR}/libclapack.a) # 'libclapack.a' added


(Perhaps/probably there exists a more elegant way, but it works ...)

