
SET (Boost_SRCPATH /Users/anaegel/Software/boost-trunk/)
SET (Boost_LIBPATH /Users/anaegel/Software/boost-trunk/)

INCLUDE_DIRECTORIES (../solver/ ${Boost_SRCPATH})
ADD_EXECUTABLE(dsolver demo_solver.cc)
