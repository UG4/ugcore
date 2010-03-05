INCLUDES = -IAMG -IAlgebra -Inumprocs -Iblocks -I. 

DEBUG_FLAGS = -O0 -mfix-and-continue -DDEBUG= -DUSE_COMMAND_LINE
RELEASE_FLAGS = -O3 -DUSE_COMMAND_LINE -funroll-loops -ftree-vectorize -msse3 -mssse3

ifdef VARIABLE_STORAGE
FLAGS := $(FLAGS) -DVARIABLE_STORAGE
endif

ifdef UNKNOWN_NR
FLAGS := $(FLAGS) -DUNKNOWN_NR=$(UNKNOWN_NR)
else
FLAGS := $(FLAGS) -DUNKNOWN_NR=1
endif


DEPENDENCIES = main.cpp disc.cpp algebra.h blocks/arrayStorage.h blocks/blockDenseMatrix.h blocks/blockSparseMatrix.h blocks/blockMatrix.h blocks/blockVector.h \
Algebra/matrixRow.h Algebra/matrixRow.hpp Algebra/sparseMatrix.h Algebra/sparseMatrix.hpp Algebra/Vector.h Algebra/Vector.hpp Algebra/TemplateExpressions.h \
numprocs/linearsolver.h numprocs/preconditioner.h \
AMG/graph.h AMG/CoarseSolver.h AMG/CoarseSolver.hpp AMG/amg.hpp AMG/amg.h AMG/maxheap.h

release: $(DEPENDENCIES)
	gcc $(RELEASE_FLAGS) $(INCLUDES) $(FLAGS) \
	-o flexamg main.cpp misc.cpp disc.cpp -lstdc++ -framework vecLib 

debug: $(DEPENDENCIES)
	gcc $(DEBUG_FLAGS) $(INCLUDES) $(FLAGS) \
	-o flexamg main.cpp misc.cpp disc.cpp -lstdc++ -framework vecLib

clean:
	rm -f *.o
	rm flexamg